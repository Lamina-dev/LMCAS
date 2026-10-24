#include <stdexcept>
#include "internal/expression_transform.hpp"
#include "limit_result.hpp"
#include "series_engine.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/series_support.hpp"
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <variant>

namespace LMCAS {

namespace detail::series_support {
double series_number_value(const NumberNode& number) {
    if (std::holds_alternative<BigInt>(number.value())) {
        return std::get<BigInt>(number.value()).to_double();
    }
    return std::get<Rational>(number.value()).to_double();
}

bool series_is_number(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return false;
    }
    return expr->is_number();
}

double series_get_double(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return 0.0;
    }
    auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr));
    if (!num) {
        return 0.0;
    }
    return series_number_value(*num);
}

bool series_is_infinity(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return false;
    }
    auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(expr));
    if (func && func->type() == FunctionNode::FuncType::Infinity) {
        return true;
    }
    if (series_is_number(expr)) {
        return std::isinf(series_get_double(expr));
    }
    return false;
}
}

using detail::series_support::series_number_value;
using detail::series_support::series_is_number;
using detail::series_support::series_get_double;
using detail::series_support::series_is_infinity;

using detail::series_support::validate_power_series_coefficients;
using detail::series_support::supported_laurent_integer_power;
using detail::series_support::validate_series_variable;

static ConvergenceInfoResult convergence_test_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&,
    ComputationContext&);
static bool series_index_variable(const std::shared_ptr<const VariableNode>& variable,
                                  const std::string& index_var) {
    return variable && !variable->is_constant() && variable->name() == index_var;
}

static std::optional<std::shared_ptr<SymbolicExpr>> coefficient_power_radius(
    const PowerNode& power, const std::string& index_var, ComputationContext& context) {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    auto exponent_variable = std::dynamic_pointer_cast<const VariableNode>(
        power.exponent());
    if (series_index_variable(exponent_variable, index_var) &&
        !expression_depends_on_variable(power.base(), index_var)) {
        auto absolute = LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Abs,
                std::vector<std::shared_ptr<const SymbolicNode>>{
                    power.base()}));
        return detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(SymbolicExpr::number(1), absolute)), transform_budget));
    }
    auto base_variable = std::dynamic_pointer_cast<const VariableNode>(
        power.base());
    if (series_index_variable(base_variable, index_var) &&
        !expression_depends_on_variable(power.exponent(), index_var)) {
        return SymbolicExpr::number(1);
    }
    return std::nullopt;
}

ExpressionResult convergence_radius_checked(
    const std::shared_ptr<SymbolicExpr>& general_coefficient,
    const std::string& index_var,
    ComputationContext& context) {
    constexpr const char* operation = "convergence_radius";
    if (!general_coefficient || !LMCAS::detail::node(general_coefficient) ||
        index_var.empty()) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument,
            "convergence radius requires a coefficient term and index variable",
            operation);
    }
    auto step = context.consume_steps(2, operation);
    if (!step) {
        return ExpressionResult::failure(step.error());
    }
    auto node = LMCAS::detail::node(general_coefficient);
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return ExpressionResult::success(
            number->is_zero() ? SymbolicExpr::infinity()
                              : SymbolicExpr::number(1));
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto radius = coefficient_power_radius(*power, index_var, context);
        if (radius) {
            return ExpressionResult::success(std::move(*radius));
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
    if (!var_check) {
        return ExpressionResult::failure(var_check.error());
    }
    if (coefficients.empty()) {
        return ExpressionResult::failure(CasErrc::InvalidArgument,
                                         "coefficient list cannot be empty",
                                         operation);
    }
    auto coeff_check = validate_power_series_coefficients(
        coefficients, operation, "polynomial");
    if (!coeff_check) {
        return ExpressionResult::failure(coeff_check.error());
    }
    auto budget = context.consume_steps(coefficients.size() + 1, operation);
    if (!budget) {
        return ExpressionResult::failure(budget.error());
    }

    for (size_t index = 0; index < coefficients.size(); ++index) {
        if (expression_depends_on_variable(
                LMCAS::detail::node(coefficients[index]), var)) {
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
    if (!var_check) {
        return ConvergenceInfoResult::failure(var_check.error());
    }
    if (!general_term || !LMCAS::detail::node(general_term)) {
        return ConvergenceInfoResult::failure(CasErrc::InvalidArgument,
                                              "general term cannot be null",
                                              operation);
    }
    auto budget = context.consume_steps(32, operation);
    if (!budget) {
        return ConvergenceInfoResult::failure(budget.error());
    }
    try {
        auto analyzed =
            convergence_test_impl(general_term, index_var, context);
        if (!analyzed) {
            return analyzed;
        }
        auto info = std::move(analyzed.value());
        if (info.result == ConvergenceResult::Inconclusive) {
            return ConvergenceInfoResult::failure(
                CasErrc::Inconclusive,
                "convergence test is outside the current supported domain",
                operation);
        }
        return info;
    } catch (const CasError& error) {
        return ConvergenceInfoResult::failure(error);
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

static std::optional<ConvergenceInfo> geometric_term_convergence(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& index_var) {
    auto power = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!power) {
        return std::nullopt;
    }
    auto exponent = std::dynamic_pointer_cast<const VariableNode>(power->exponent());
    auto base = std::dynamic_pointer_cast<const NumberNode>(power->base());
    if (!exponent || exponent->is_constant() || exponent->name() != index_var || !base) {
        return std::nullopt;
    }

    bool inside_unit_circle = false;
    if (const auto* integer = std::get_if<BigInt>(&base->value())) {
        inside_unit_circle = *integer > BigInt(-1) && *integer < BigInt(1);
    } else if (const auto* rational = std::get_if<Rational>(&base->value())) {
        inside_unit_circle = *rational > Rational(-1) && *rational < Rational(1);
    }
    return ConvergenceInfo{
        inside_unit_circle ? ConvergenceResult::Convergent
                           : ConvergenceResult::Divergent,
        "geometric"};
}

static std::optional<ConvergenceInfo> power_term_convergence(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& index_var) {
    auto power = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!power) {
        return std::nullopt;
    }
    auto base_var = std::dynamic_pointer_cast<const VariableNode>(power->base());
    auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
    if (!base_var || base_var->is_constant() || base_var->name() != index_var || !exponent) {
        return std::nullopt;
    }

    bool below_negative_one = false;
    bool at_or_above_negative_one = false;
    if (const auto* integer = std::get_if<BigInt>(&exponent->value())) {
        below_negative_one = *integer < BigInt(-1);
        at_or_above_negative_one = !below_negative_one;
    } else if (const auto* rational = std::get_if<Rational>(&exponent->value())) {
        below_negative_one = *rational < Rational(-1);
        at_or_above_negative_one = !below_negative_one;
    }
    if (below_negative_one) {
        return ConvergenceInfo{
            ConvergenceResult::Convergent, "p-series"};
    }
    if (at_or_above_negative_one) {
        return ConvergenceInfo{
            ConvergenceResult::Divergent, "p-series"};
    }
    return std::nullopt;
}

static ConvergenceInfo classify_ratio_limit(
    const std::shared_ptr<SymbolicExpr>& lim, ComputationContext& context) {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    if (!lim) {
        return {ConvergenceResult::Inconclusive, ""};
    }
    auto ls = detail::make_expression_ptr(detail::normalize_node(detail::node(lim), transform_budget));
    auto infinity = ls
        ? std::dynamic_pointer_cast<const FunctionNode>(detail::node(ls))
        : nullptr;
    if (infinity && infinity->type() == FunctionNode::FuncType::Infinity) {
        return {ConvergenceResult::Divergent, "ratio"};
    }
    auto number = ls
        ? std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(ls))
        : nullptr;
    if (!number) {
        return {ConvergenceResult::Inconclusive, ""};
    }
    if (const auto* integer = std::get_if<BigInt>(&number->value())) {
        return {*integer < BigInt(1) ? ConvergenceResult::Convergent
                                     : *integer > BigInt(1) ? ConvergenceResult::Divergent
                                                            : ConvergenceResult::Inconclusive, "ratio"};
    }
    if (const auto* rational = std::get_if<Rational>(&number->value())) {
        return {*rational < Rational(1) ? ConvergenceResult::Convergent
                                        : *rational > Rational(1) ? ConvergenceResult::Divergent
                                                                  : ConvergenceResult::Inconclusive, "ratio"};
    }
    return {ConvergenceResult::Inconclusive, ""};
}

static ConvergenceInfoResult ratio_convergence_test(
    const std::shared_ptr<SymbolicExpr>& general_term,
    const std::string& index_var, ComputationContext& context) try {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    auto n = SymbolicExpr::variable(index_var);
    auto n1 = SymbolicExpr::add(n, SymbolicExpr::number(1));
    auto inf = SymbolicExpr::infinity();
    auto a_n1 = detail::make_expression_ptr(detail::substitute_node(detail::node(general_term), index_var, detail::node(n1), transform_budget));
    if (a_n1) {
        a_n1 = detail::make_expression_ptr(detail::normalize_node(detail::node(a_n1), transform_budget));
        auto ratio = SymbolicExpr::divide(a_n1, general_term);
        if (ratio) {
            ratio = detail::make_expression_ptr(detail::normalize_node(detail::node(ratio), transform_budget));
            auto abs_r = LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Abs, std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(ratio)}));
            abs_r = detail::make_expression_ptr(detail::normalize_node(detail::node(abs_r), transform_budget));
            auto limited = limit_expression_checked(
                abs_r, index_var, inf, LimitDirection::Both, context);
            if (!limited) {
                return ConvergenceInfoResult::failure(limited.error());
            }
            auto lim = std::move(limited.value());
            return classify_ratio_limit(lim, context);
        }
    }
    return ConvergenceInfo{ConvergenceResult::Inconclusive, ""};
} catch (const CasError& error) {
    return ConvergenceInfoResult::failure(error);
} catch (const std::bad_alloc&) {
    return ConvergenceInfoResult::failure(CasErrc::ResourceLimit, "allocation failed", "ratio_convergence_test");
} catch (const std::length_error& error) {
    return ConvergenceInfoResult::failure(CasErrc::ResourceLimit, error.what(), "ratio_convergence_test");
} catch (const std::exception& error) {
    return ConvergenceInfoResult::failure(CasErrc::InternalInvariant, error.what(), "ratio_convergence_test");
}

static ConvergenceInfoResult convergence_test_impl(
    const std::shared_ptr<SymbolicExpr>& general_term,
    const std::string& index_var,
    ComputationContext& context) {
    if (!general_term || !LMCAS::detail::node(general_term)) {
        return ConvergenceInfo{
            ConvergenceResult::Inconclusive, ""};
    }
    if (auto geometric = geometric_term_convergence(
            detail::node(general_term), index_var)) {
        return std::move(*geometric);
    }
    auto power_info = power_term_convergence(detail::node(general_term), index_var);
    if (power_info) {
        return std::move(*power_info);
    }
    /**
     * @brief 对 Laurent 单项式 c*n^e 直接应用 p 级数判据。
     * 倒数幂沿用线性规则，省去 Abs 比值极限的递归化简。
     */
    if (const auto laurent_power = supported_laurent_integer_power(
            LMCAS::detail::node(general_term), index_var)) {
        if (*laurent_power < -1) {
            return ConvergenceInfo{
                ConvergenceResult::Convergent, "p-series"};
        }
        if (*laurent_power >= -1) {
            return ConvergenceInfo{
                ConvergenceResult::Divergent, "p-series"};
        }
    }
    return ratio_convergence_test(general_term, index_var, context);
}

}
