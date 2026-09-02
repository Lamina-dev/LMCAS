#include "series_engine.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/series_support.hpp"
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <variant>

static bool series_is_number(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !lamina::detail::node(expr)) return false;
    return expr->is_number();
}

static double series_get_double(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !lamina::detail::node(expr)) return 0.0;
    auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr));
    if (!num) return 0.0;
    if (std::holds_alternative<BigInt>(num->value()))
        return std::get<BigInt>(num->value()).to_double();
    if (std::holds_alternative<Rational>(num->value()))
        return std::get<Rational>(num->value()).to_double();
    return static_cast<double>(std::get<lmmc_real_t>(num->value()));
}

static bool series_is_infinity(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !lamina::detail::node(expr)) return false;
    auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr));
    if (func && func->type() == FunctionNode::FuncType::Infinity) return true;
    if (series_is_number(expr)) return std::isinf(series_get_double(expr));
    return false;
}

[[maybe_unused]] static bool series_depends_on(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    if (!expr || !lamina::detail::node(expr)) return false;
    return expr->to_string().find(var) != std::string::npos;
}

static bool series_extract_alternating(const std::shared_ptr<const SymbolicNode>& node,
                                       const std::string& n,
                                       std::shared_ptr<SymbolicExpr>& remainder) {
    auto pow = std::dynamic_pointer_cast<const PowerNode>(node);
    if (pow) {
        auto base_num = std::dynamic_pointer_cast<const NumberNode>(pow->base());
        auto exp_var = std::dynamic_pointer_cast<const VariableNode>(pow->exponent());
        if (base_num && exp_var && exp_var->name() == n) {
            double base_val = 0.0;
            if (std::holds_alternative<BigInt>(base_num->value()))
                base_val = std::get<BigInt>(base_num->value()).to_double();
            else if (std::holds_alternative<Rational>(base_num->value()))
                base_val = std::get<Rational>(base_num->value()).to_double();
            else
                base_val = static_cast<double>(std::get<lmmc_real_t>(base_num->value()));
            if (std::abs(base_val + 1.0) < 1e-12) {
                remainder = SymbolicExpr::number(1);
                return true;
            }
        }
    }
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (mul) {
        for (size_t i = 0; i < mul->operands().size(); ++i) {
            auto pw = std::dynamic_pointer_cast<const PowerNode>(mul->operands()[i]);
            if (!pw) continue;
            auto base_num = std::dynamic_pointer_cast<const NumberNode>(pw->base());
            auto exp_var = std::dynamic_pointer_cast<const VariableNode>(pw->exponent());
            if (!base_num || !exp_var || exp_var->name() != n) continue;
            double base_val = 0.0;
            if (std::holds_alternative<BigInt>(base_num->value()))
                base_val = std::get<BigInt>(base_num->value()).to_double();
            else if (std::holds_alternative<Rational>(base_num->value()))
                base_val = std::get<Rational>(base_num->value()).to_double();
            else
                base_val = static_cast<double>(std::get<lmmc_real_t>(base_num->value()));
            if (std::abs(base_val + 1.0) < 1e-12) {
                std::vector<std::shared_ptr<const SymbolicNode>> rest;
                for (size_t j = 0; j < mul->operands().size(); ++j)
                    if (j != i) rest.push_back(mul->operands()[j]);
                if (rest.empty()) remainder = SymbolicExpr::number(1);
                else if (rest.size() == 1) remainder = lamina::detail::make_expression_ptr(rest[0]);
                else remainder = lamina::detail::make_expression_ptr(lamina::detail::make_node<MultiplyNode>(rest));
                return true;
            }
        }
    }
    return false;
}

static bool series_detect_trig_oscillation(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& n,
    std::shared_ptr<SymbolicExpr>& amplitude,
    lamina::ComputationContext& context) {
    auto func = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (func && (func->type() == FunctionNode::FuncType::Sin || func->type() == FunctionNode::FuncType::Cos)) {
        if (!func->arguments().empty()) {
            auto arg_expr = lamina::detail::make_expression_ptr(func->arguments()[0]);
            auto arg_limit = lamina::detail::propagate_result(
                lamina::limit_expression_checked(
                    arg_expr, n, SymbolicExpr::infinity(),
                    LimitDirection::Both, context));
            if (series_is_infinity(arg_limit)) { amplitude = SymbolicExpr::number(1); return true; }
        }
    }
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!mul) return false;
    for (size_t i = 0; i < mul->operands().size(); ++i) {
        auto f = std::dynamic_pointer_cast<const FunctionNode>(mul->operands()[i]);
        if (!f || (f->type() != FunctionNode::FuncType::Sin && f->type() != FunctionNode::FuncType::Cos)) continue;
        if (f->arguments().empty()) continue;
        auto arg_expr = lamina::detail::make_expression_ptr(f->arguments()[0]);
        auto arg_limit = lamina::detail::propagate_result(
            lamina::limit_expression_checked(
                arg_expr, n, SymbolicExpr::infinity(),
                LimitDirection::Both, context));
        if (!series_is_infinity(arg_limit)) continue;
        std::vector<std::shared_ptr<const SymbolicNode>> rest;
        for (size_t j = 0; j < mul->operands().size(); ++j)
            if (j != i) rest.push_back(mul->operands()[j]);
        if (rest.empty()) amplitude = SymbolicExpr::number(1);
        else if (rest.size() == 1) amplitude = lamina::detail::make_expression_ptr(rest[0]);
        else amplitude = lamina::detail::make_expression_ptr(lamina::detail::make_node<MultiplyNode>(rest));
        return true;
    }
    return false;
}

static std::shared_ptr<SymbolicExpr> series_abs(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr) return nullptr;
    if (series_is_number(expr)) return SymbolicExpr::number(std::abs(series_get_double(expr)));
    return lamina::detail::make_expression_ptr(lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(expr)}));
}

static std::shared_ptr<SymbolicExpr> series_negate(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr) return nullptr;
    return SymbolicExpr::multiply(SymbolicExpr::number(-1), expr)->simplify();
}
namespace lamina {

using detail::series_support::validate_power_series_coefficients;
using detail::series_support::supported_laurent_integer_power;
using detail::series_support::validate_series_variable;

static ConvergenceInfo convergence_test_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&,
    ComputationContext&);
ExpressionResult convergence_radius_checked(
    const std::shared_ptr<SymbolicExpr>& general_coefficient,
    const std::string& index_var,
    ComputationContext& context) {
    constexpr const char* operation = "convergence_radius";
    if (!general_coefficient || !lamina::detail::node(general_coefficient) ||
        index_var.empty()) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument,
            "convergence radius requires a coefficient term and index variable",
            operation);
    }
    auto step = context.consume_steps(2, operation);
    if (!step) return ExpressionResult::failure(step.error());
    auto node = lamina::detail::node(general_coefficient);
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return ExpressionResult::success(
            number->is_zero() ? SymbolicExpr::infinity()
                              : SymbolicExpr::number(1));
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent_variable = std::dynamic_pointer_cast<const VariableNode>(
            power->exponent());
        if (exponent_variable && exponent_variable->name() == index_var &&
            !expression_depends_on_variable(power->base(), index_var)) {
            auto absolute = lamina::detail::make_expression_ptr(
                lamina::detail::make_node<FunctionNode>(
                    FunctionNode::FuncType::Abs,
                    std::vector<std::shared_ptr<const SymbolicNode>>{
                        power->base()}));
            return ExpressionResult::success(
                SymbolicExpr::divide(SymbolicExpr::number(1), absolute)->simplify());
        }
        auto base_variable = std::dynamic_pointer_cast<const VariableNode>(
            power->base());
        if (base_variable && base_variable->name() == index_var &&
            !expression_depends_on_variable(power->exponent(), index_var)) {
            return ExpressionResult::success(SymbolicExpr::number(1));
        }
    }
    return ExpressionResult::failure(
        CasErrc::Inconclusive,
        "coefficient-term limit is outside the supported ratio/root-test domain",
        operation);
}

ExpressionResult convergence_radius_checked(
    const std::shared_ptr<SymbolicExpr>& general_coefficient,
    const std::string& index_var) {
    ComputationContext context;
    return convergence_radius_checked(
        general_coefficient, index_var, context);
}

ExpressionResult convergence_radius_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coefficients,
    const std::string& var,
    ComputationContext& context)
{
    const std::string operation = "convergence_radius";
    auto var_check = validate_series_variable(var, context, operation);
    if (!var_check) return ExpressionResult::failure(var_check.error());
    if (coefficients.empty()) {
        return ExpressionResult::failure(CasErrc::InvalidArgument,
                                         "coefficient list cannot be empty",
                                         operation);
    }
    auto coeff_check = validate_power_series_coefficients(
        coefficients, operation, "polynomial");
    if (!coeff_check) return ExpressionResult::failure(coeff_check.error());
    auto budget = context.consume_steps(coefficients.size() + 1, operation);
    if (!budget) return ExpressionResult::failure(budget.error());

    for (size_t index = 0; index < coefficients.size(); ++index) {
        if (expression_depends_on_variable(
                lamina::detail::node(coefficients[index]), var)) {
            return ExpressionResult::failure(
                CasErrc::InvalidArgument,
                "coefficient at index " + std::to_string(index) +
                    " depends on the series variable",
                operation);
        }
    }

    try {
        /// 有限系数列表完整地表示一个多项式,其收敛半径恒为无穷.
        /// 无限级数必须使用通项重载,有限前缀不能证明其尾项行为.
        return ExpressionResult::success(SymbolicExpr::infinity());
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while constructing convergence radius",
            operation);
    } catch (const std::exception& ex) {
        return ExpressionResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

ExpressionResult convergence_radius_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coefficients,
    const std::string& var)
{
    ComputationContext context;
    return convergence_radius_checked(coefficients, var, context);
}


ConvergenceInfoResult convergence_test_checked(
    const std::shared_ptr<SymbolicExpr>& general_term,
    const std::string& index_var,
    ComputationContext& context)
{
    const std::string operation = "convergence_test";
    auto var_check = validate_series_variable(index_var, context, operation);
    if (!var_check) return ConvergenceInfoResult::failure(var_check.error());
    if (!general_term || !lamina::detail::node(general_term)) {
        return ConvergenceInfoResult::failure(CasErrc::InvalidArgument,
                                              "general term cannot be null",
                                              operation);
    }
    auto budget = context.consume_steps(32, operation);
    if (!budget) return ConvergenceInfoResult::failure(budget.error());
    try {
        auto info = convergence_test_impl(general_term, index_var, context);
        if (info.result == ConvergenceResult::Inconclusive) {
            return ConvergenceInfoResult::failure(
                CasErrc::Inconclusive,
                "convergence test is outside the current supported domain",
                operation);
        }
        return ConvergenceInfoResult::success(std::move(info));
    } catch (const detail::ResultPropagation& propagation) {
        return ConvergenceInfoResult::failure(propagation.error());
    } catch (const std::bad_alloc&) {
        return ConvergenceInfoResult::failure(CasErrc::ResourceLimit,
                                              "allocation failed while testing convergence",
                                              operation);
    } catch (const std::exception& ex) {
        return ConvergenceInfoResult::failure(CasErrc::InternalInvariant,
                                              ex.what(),
                                              operation);
    }
}

ConvergenceInfoResult convergence_test_checked(
    const std::shared_ptr<SymbolicExpr>& general_term,
    const std::string& index_var)
{
    ComputationContext context;
    return convergence_test_checked(general_term, index_var, context);
}

static ConvergenceInfo convergence_test_impl(
    const std::shared_ptr<SymbolicExpr>& general_term,
    const std::string& index_var,
    ComputationContext& context) {
    if (!general_term || !lamina::detail::node(general_term)) return {ConvergenceResult::Inconclusive, ""};
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(general_term))) {
        auto base_var = std::dynamic_pointer_cast<const VariableNode>(power->base());
        auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
        if (base_var && base_var->name() == index_var && exponent) {
            double p = 0.0;
            if (std::holds_alternative<BigInt>(exponent->value())) {
                p = std::get<BigInt>(exponent->value()).to_double();
            } else if (std::holds_alternative<Rational>(exponent->value())) {
                p = std::get<Rational>(exponent->value()).to_double();
            } else {
                p = static_cast<double>(std::get<lmmc_real_t>(exponent->value()));
            }
            if (p < -1.0) return {ConvergenceResult::Convergent, "p-series"};
            if (p >= -1.0 && p < 0.0) return {ConvergenceResult::Divergent, "p-series"};
        }
    }
    /// Laurent 单项式 c*n^e 直接应用 p 级数判据.
    /// 该路径使倒数幂保持在线性规则内,并跳过 Abs 比值极限的递归化简.
    if (const auto laurent_power = supported_laurent_integer_power(
            lamina::detail::node(general_term), index_var)) {
        if (*laurent_power < -1) {
            return {ConvergenceResult::Convergent, "p-series"};
        }
        if (*laurent_power >= -1 && *laurent_power < 0) {
            return {ConvergenceResult::Divergent, "p-series"};
        }
    }
    auto n = SymbolicExpr::variable(index_var);
    auto n1 = SymbolicExpr::add(n, SymbolicExpr::number(1));
    auto inf = SymbolicExpr::infinity();
    auto a_n1 = general_term->substitute(index_var, n1);
    if (a_n1) {
        a_n1 = a_n1->simplify();
        auto ratio = SymbolicExpr::divide(a_n1, general_term);
        if (ratio) {
            ratio = ratio->simplify();
            auto abs_r = lamina::detail::make_expression_ptr(lamina::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Abs, std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(ratio)}));
            abs_r = abs_r->simplify();
            auto lim = detail::propagate_result(limit_expression_checked(
                abs_r, index_var, inf, LimitDirection::Both, context));
            if (lim) {
                auto ls = lim->simplify();
                if (ls && series_is_number(ls)) {
                    double v = series_get_double(ls);
                    if (v < 1.0 - 1e-12) return {ConvergenceResult::Convergent, "ratio"};
                    if (v > 1.0 + 1e-12) return {ConvergenceResult::Divergent, "ratio"};
                }
                if (ls && ls->is_zero()) return {ConvergenceResult::Convergent, "ratio"};
                if (series_is_infinity(ls)) return {ConvergenceResult::Divergent, "ratio"};
            }
        }
    }
    return {ConvergenceResult::Inconclusive, ""};
}
ExpressionResult lim_sup_checked(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n,
    ComputationContext& context) {
    constexpr const char* operation = "lim_sup";
    if (!a_n || !detail::node(a_n) || n.empty()) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument,
            "upper limit requires a sequence term and index variable",
            operation);
    }
    auto budget = context.consume_steps(1, operation);
    if (!budget) return ExpressionResult::failure(budget.error());
    try {
        std::shared_ptr<SymbolicExpr> remainder;
        if (series_extract_alternating(detail::node(a_n), n, remainder)) {
            auto limit = detail::propagate_result(limit_expression_checked(
                remainder, n, SymbolicExpr::infinity(),
                LimitDirection::Both, context));
            auto simplified = limit ? limit->simplify() : nullptr;
            if (simplified && !series_is_infinity(simplified)) {
                return ExpressionResult::success(series_abs(simplified));
            }
            if (series_is_infinity(simplified)) {
                return ExpressionResult::success(SymbolicExpr::infinity());
            }
            return ExpressionResult::success(series_abs(limit));
        }

        std::shared_ptr<SymbolicExpr> amplitude;
        if (series_detect_trig_oscillation(
                detail::node(a_n), n, amplitude, context)) {
            auto limit = detail::propagate_result(limit_expression_checked(
                amplitude, n, SymbolicExpr::infinity(),
                LimitDirection::Both, context));
            auto simplified = limit ? limit->simplify() : nullptr;
            if (simplified && series_is_number(simplified)) {
                return ExpressionResult::success(SymbolicExpr::number(
                    std::abs(series_get_double(simplified))));
            }
            if (simplified) {
                return ExpressionResult::success(series_abs(simplified));
            }
        }

        auto limit = detail::propagate_result(limit_expression_checked(
            a_n, n, SymbolicExpr::infinity(),
            LimitDirection::Both, context));
        auto simplified = limit ? limit->simplify() : nullptr;
        return ExpressionResult::success(
            simplified ? std::move(simplified) : std::move(limit));
    } catch (const detail::ResultPropagation& propagation) {
        return ExpressionResult::failure(propagation.error());
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while calculating upper limit", operation);
    } catch (const std::exception& ex) {
        return ExpressionResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

ExpressionResult lim_sup_checked(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n) {
    ComputationContext context;
    return lim_sup_checked(a_n, n, context);
}

ExpressionResult lim_inf_checked(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n,
    ComputationContext& context) {
    constexpr const char* operation = "lim_inf";
    if (!a_n || !detail::node(a_n) || n.empty()) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument,
            "lower limit requires a sequence term and index variable",
            operation);
    }
    auto budget = context.consume_steps(1, operation);
    if (!budget) return ExpressionResult::failure(budget.error());
    try {
        std::shared_ptr<SymbolicExpr> remainder;
        if (series_extract_alternating(detail::node(a_n), n, remainder)) {
            auto limit = detail::propagate_result(limit_expression_checked(
                remainder, n, SymbolicExpr::infinity(),
                LimitDirection::Both, context));
            auto simplified = limit ? limit->simplify() : nullptr;
            if (simplified && !series_is_infinity(simplified)) {
                return ExpressionResult::success(
                    series_negate(series_abs(simplified)));
            }
            if (series_is_infinity(simplified)) {
                return ExpressionResult::success(
                    SymbolicExpr::infinity(-1));
            }
            return ExpressionResult::success(
                series_negate(series_abs(limit)));
        }

        std::shared_ptr<SymbolicExpr> amplitude;
        if (series_detect_trig_oscillation(
                detail::node(a_n), n, amplitude, context)) {
            auto limit = detail::propagate_result(limit_expression_checked(
                amplitude, n, SymbolicExpr::infinity(),
                LimitDirection::Both, context));
            auto simplified = limit ? limit->simplify() : nullptr;
            if (simplified && series_is_number(simplified)) {
                return ExpressionResult::success(SymbolicExpr::number(
                    -std::abs(series_get_double(simplified))));
            }
            if (simplified) {
                return ExpressionResult::success(
                    series_negate(series_abs(simplified)));
            }
        }

        auto limit = detail::propagate_result(limit_expression_checked(
            a_n, n, SymbolicExpr::infinity(),
            LimitDirection::Both, context));
        auto simplified = limit ? limit->simplify() : nullptr;
        return ExpressionResult::success(
            simplified ? std::move(simplified) : std::move(limit));
    } catch (const detail::ResultPropagation& propagation) {
        return ExpressionResult::failure(propagation.error());
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while calculating lower limit", operation);
    } catch (const std::exception& ex) {
        return ExpressionResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

ExpressionResult lim_inf_checked(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n) {
    ComputationContext context;
    return lim_inf_checked(a_n, n, context);
}
} // namespace lamina
