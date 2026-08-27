/**
 * @file series_engine.cpp
 */

#include "series_engine.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "integration.hpp"
#include <cmath>
#include <memory>
#include <string>
#include <cstdio>
#include <typeinfo>
#include <algorithm>
#include <variant>
#include <limits>
#include <optional>

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

static bool series_detect_trig_oscillation(const std::shared_ptr<const SymbolicNode>& node,
                                           const std::string& n,
                                           std::shared_ptr<SymbolicExpr>& amplitude) {
    auto func = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (func && (func->type() == FunctionNode::FuncType::Sin || func->type() == FunctionNode::FuncType::Cos)) {
        if (!func->arguments().empty()) {
            auto arg_expr = lamina::detail::make_expression_ptr(func->arguments()[0]);
            auto arg_limit = arg_expr->limit(n, SymbolicExpr::infinity());
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
        auto arg_limit = arg_expr->limit(n, SymbolicExpr::infinity());
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

static bool try_get_int(const std::shared_ptr<SymbolicExpr>& expr, long long& out) {
    if (!expr || !lamina::detail::node(expr)) return false;
    auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr));
    if (!num) return false;
    if (std::holds_alternative<BigInt>(num->value())) {
        out = static_cast<long long>(std::get<BigInt>(num->value()).to_int());
        return true;
    }
    if (std::holds_alternative<Rational>(num->value())) {
        auto& r = std::get<Rational>(num->value());
        if (r.get_denominator() == BigInt(1)) {
            out = static_cast<long long>(r.get_numerator().to_int());
            return true;
        }
        return false;
    }
    double d = static_cast<double>(std::get<lmmc_real_t>(num->value()));
    if (d == std::floor(d) && std::abs(d) < 1e15) { out = static_cast<long long>(d); return true; }
    return false;
}

[[maybe_unused]] static std::shared_ptr<SymbolicExpr> make_pi() { return SymbolicExpr::number(LMMC_CONST_PI); }

static int detect_parity(const std::shared_ptr<SymbolicExpr>& f, const std::string& var) {
    if (!f) return 0;
    auto neg_x = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::variable(var));
    auto f_neg = f->substitute(var, neg_x);
    if (!f_neg) return 0;
    f_neg = f_neg->simplify();
    auto f_s = f->simplify();
    auto diff_even = SymbolicExpr::add(f_neg, SymbolicExpr::multiply(SymbolicExpr::number(-1), f_s));
    if (diff_even) { diff_even = diff_even->simplify(); if (diff_even->is_zero()) return 1; }
    auto diff_odd = SymbolicExpr::add(f_neg, f_s);
    if (diff_odd) { diff_odd = diff_odd->simplify(); if (diff_odd->is_zero()) return -1; }
    return 0;
}

namespace lamina {

namespace {

Result<void> validate_power_series_coefficients(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coeffs,
    const std::string& operation,
    const std::string& name)
{
    for (size_t i = 0; i < coeffs.size(); ++i) {
        if (!coeffs[i] || !lamina::detail::node(coeffs[i])) {
            return Result<void>::failure(
                CasErrc::InvalidArgument,
                name + " contains a null coefficient at index " + std::to_string(i),
                operation);
        }
    }
    return Result<void>::success();
}

Result<void> validate_power_series_order(
    int order,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (order <= 0) {
        return Result<void>::failure(
            CasErrc::InvalidArgument,
            "power series truncation order must be positive",
            operation);
    }
    return Result<void>::success();
}

Result<void> validate_series_variable(const std::string& var,
                                      ComputationContext& context,
                                      const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (var.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "series variable name cannot be empty",
                                     operation);
    }
    return Result<void>::success();
}

Result<void> validate_series_expr_point(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& center,
    const std::string& var,
    ComputationContext& context,
    const std::string& operation)
{
    auto var_check = validate_series_variable(var, context, operation);
    if (!var_check) return var_check;
    if (!expr || !lamina::detail::node(expr) || !center || !lamina::detail::node(center)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "series expression and center cannot be null",
                                     operation);
    }
    return Result<void>::success();
}

bool has_symbolic_nonzero_coefficient(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coefficients)
{
    for (const auto& coefficient : coefficients) {
        if (!coefficient || !lamina::detail::node(coefficient) || coefficient->is_zero()) continue;
        if (!series_is_number(coefficient)) return true;
    }
    return false;
}

Result<void> validate_laurent_orders(int order_neg,
                                     int order_pos,
                                     const std::string& operation)
{
    if (order_neg < 0 || order_pos < 0 || order_neg > 64 || order_pos > 64) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "Laurent truncation orders must be between 0 and 64",
                                     operation);
    }
    return Result<void>::success();
}

std::optional<int> integer_value_from_node(const std::shared_ptr<const SymbolicNode>& node)
{
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return std::nullopt;
    if (std::holds_alternative<BigInt>(number->value())) {
        return std::get<BigInt>(number->value()).to_int();
    }
    if (std::holds_alternative<Rational>(number->value())) {
        const auto& rational = std::get<Rational>(number->value());
        if (rational.get_denominator() != BigInt(1)) return std::nullopt;
        return rational.get_numerator().to_int();
    }
    double value = static_cast<double>(std::get<lmmc_real_t>(number->value()));
    if (!std::isfinite(value) || value != std::floor(value) ||
        std::abs(value) > static_cast<double>(std::numeric_limits<int>::max())) {
        return std::nullopt;
    }
    return static_cast<int>(value);
}

std::optional<int> supported_laurent_integer_power(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var)
{
    if (!node) return std::nullopt;
    if (std::dynamic_pointer_cast<const NumberNode>(node)) return 0;
    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        return variable->name() == var ? std::optional<int>(1) : std::nullopt;
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent = integer_value_from_node(power->exponent());
        if (!exponent) return std::nullopt;
        auto base = std::dynamic_pointer_cast<const VariableNode>(power->base());
        if (base && base->name() == var) return *exponent;
        auto nested = supported_laurent_integer_power(power->base(), var);
        if (!nested) return std::nullopt;
        return *nested * *exponent;
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        int total_power = 0;
        for (const auto& operand : multiply->operands()) {
            if (std::dynamic_pointer_cast<const NumberNode>(operand)) continue;
            auto power = supported_laurent_integer_power(operand, var);
            if (!power) return std::nullopt;
            total_power += *power;
        }
        return total_power;
    }
    return std::nullopt;
}

} // namespace

/// Stubs for functions not yet fully implemented in this task
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
    auto coeff_check = validate_power_series_coefficients(coefficients, operation, "series");
    if (!coeff_check) return ExpressionResult::failure(coeff_check.error());
    if (has_symbolic_nonzero_coefficient(coefficients)) {
        return ExpressionResult::failure(
            CasErrc::Inconclusive,
            "convergence radius for symbolic coefficients is outside the current support domain",
            operation);
    }
    auto budget = context.consume_steps(coefficients.size() * 4 + 4, operation);
    if (!budget) return ExpressionResult::failure(budget.error());
    try {
        auto radius = convergence_radius(coefficients, var);
        if (!radius || !lamina::detail::node(radius)) {
            return ExpressionResult::failure(
                CasErrc::Inconclusive,
                "convergence radius could not be constructed in the supported domain",
                operation);
        }
        return ExpressionResult::success(radius);
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                         "allocation failed while calculating convergence radius",
                                         operation);
    } catch (const std::exception& ex) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                         ex.what(),
                                         operation);
    }
}

ExpressionResult convergence_radius_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coefficients,
    const std::string& var)
{
    ComputationContext context;
    return convergence_radius_checked(coefficients, var, context);
}

std::shared_ptr<SymbolicExpr> convergence_radius(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coefficients, const std::string& var) {
    (void)var;
    if (coefficients.size() < 2) return SymbolicExpr::infinity();
    std::vector<size_t> nz;
    for (size_t i = 0; i < coefficients.size(); ++i)
        if (coefficients[i] && !coefficients[i]->is_zero()) nz.push_back(i);
    if (nz.size() < 2) return SymbolicExpr::infinity();
    std::vector<double> ratios;
    for (size_t i = 0; i + 1 < nz.size(); ++i) {
        auto a = coefficients[nz[i]]; auto b = coefficients[nz[i+1]];
        if (!a || !b || b->is_zero()) continue;
        if (series_is_number(a) && series_is_number(b)) {
            double va = std::abs(series_get_double(a)), vb = std::abs(series_get_double(b));
            if (vb > 1e-300) ratios.push_back(va / vb);
        }
    }
    if (!ratios.empty()) {
        double r = ratios.back();
        if (std::isinf(r) || r > 1e15) return SymbolicExpr::infinity();
        return SymbolicExpr::number(r);
    }
    return SymbolicExpr::infinity();
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
        auto info = convergence_test(general_term, index_var);
        if (info.result == ConvergenceResult::Inconclusive) {
            return ConvergenceInfoResult::failure(
                CasErrc::Inconclusive,
                "convergence test is outside the current supported domain",
                operation);
        }
        return ConvergenceInfoResult::success(std::move(info));
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

ConvergenceInfo convergence_test(
    const std::shared_ptr<SymbolicExpr>& general_term, const std::string& index_var) {
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
    // Fast path: a Laurent monomial c*n^e resolves by the p-series rule
    // without the ratio limit, whose Abs simplification can fail to
    // terminate for reciprocal powers such as 1/(n^2).
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
            auto lim = abs_r->limit(index_var, inf);
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

std::vector<std::shared_ptr<SymbolicExpr>> power_series_add(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b) {
    size_t len = std::max(a.size(), b.size());
    std::vector<std::shared_ptr<SymbolicExpr>> result(len);
    for (size_t i = 0; i < len; ++i) {
        auto ai = (i < a.size() && a[i]) ? a[i] : SymbolicExpr::number(0);
        auto bi = (i < b.size() && b[i]) ? b[i] : SymbolicExpr::number(0);
        result[i] = SymbolicExpr::add(ai, bi)->simplify();
    }
    return result;
}

PowerSeriesResult power_series_multiply_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b, int order,
    ComputationContext& context) {
    const std::string operation = "power_series_multiply";
    auto order_check = validate_power_series_order(order, context, operation);
    if (!order_check) return PowerSeriesResult::failure(order_check.error());
    auto a_check = validate_power_series_coefficients(a, operation, "left series");
    if (!a_check) return PowerSeriesResult::failure(a_check.error());
    auto b_check = validate_power_series_coefficients(b, operation, "right series");
    if (!b_check) return PowerSeriesResult::failure(b_check.error());

    size_t n = static_cast<size_t>(order);
    std::vector<std::shared_ptr<SymbolicExpr>> result(n, SymbolicExpr::number(0));
    for (size_t k = 0; k < n; ++k) {
        auto step = context.consume_steps(1, operation);
        if (!step) return PowerSeriesResult::failure(step.error());
        std::vector<std::shared_ptr<const SymbolicNode>> terms;
        for (size_t j = 0; j <= k; ++j) {
            if (j >= a.size() || (k-j) >= b.size()) continue;
            auto aj = a[j] ? a[j] : SymbolicExpr::number(0);
            auto bkj = b[k-j] ? b[k-j] : SymbolicExpr::number(0);
            terms.push_back(lamina::detail::node(SymbolicExpr::multiply(aj, bkj)));
        }
        if (terms.empty()) result[k] = SymbolicExpr::number(0);
        else if (terms.size() == 1) result[k] = lamina::detail::make_expression_ptr(terms[0])->simplify();
        else result[k] = lamina::detail::make_expression_ptr(lamina::detail::make_node<AddNode>(terms))->simplify();
    }
    return PowerSeriesResult::success(std::move(result));
}

PowerSeriesResult power_series_multiply_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b, int order) {
    ComputationContext context;
    return power_series_multiply_checked(a, b, order, context);
}

std::vector<std::shared_ptr<SymbolicExpr>> power_series_multiply(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b, int order) {
    if (order <= 0) return {};
    size_t n = static_cast<size_t>(order);
    std::vector<std::shared_ptr<SymbolicExpr>> result(n, SymbolicExpr::number(0));
    for (size_t k = 0; k < n; ++k) {
        std::vector<std::shared_ptr<const SymbolicNode>> terms;
        for (size_t j = 0; j <= k; ++j) {
            if (j >= a.size() || (k-j) >= b.size()) continue;
            auto aj = a[j] ? a[j] : SymbolicExpr::number(0);
            auto bkj = b[k-j] ? b[k-j] : SymbolicExpr::number(0);
            terms.push_back(lamina::detail::node(SymbolicExpr::multiply(aj, bkj)));
        }
        if (terms.empty()) result[k] = SymbolicExpr::number(0);
        else if (terms.size() == 1) result[k] = lamina::detail::make_expression_ptr(terms[0])->simplify();
        else result[k] = lamina::detail::make_expression_ptr(lamina::detail::make_node<AddNode>(terms))->simplify();
    }
    return result;
}

PowerSeriesResult power_series_compose_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& g, int order,
    ComputationContext& context) {
    const std::string operation = "power_series_compose";
    auto order_check = validate_power_series_order(order, context, operation);
    if (!order_check) return PowerSeriesResult::failure(order_check.error());
    auto f_check = validate_power_series_coefficients(f, operation, "outer series");
    if (!f_check) return PowerSeriesResult::failure(f_check.error());
    auto g_check = validate_power_series_coefficients(g, operation, "inner series");
    if (!g_check) return PowerSeriesResult::failure(g_check.error());

    size_t n = static_cast<size_t>(order);
    if (!g.empty() && !g[0]->is_zero()) {
        return PowerSeriesResult::failure(
            CasErrc::DomainError,
            "power series composition around zero requires g(0) = 0",
            operation);
    }
    std::vector<std::shared_ptr<SymbolicExpr>> result(n, SymbolicExpr::number(0));
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> gp(n);
    gp[0].resize(n, SymbolicExpr::number(0)); gp[0][0] = SymbolicExpr::number(1);
    if (n > 1) { gp[1].resize(n, SymbolicExpr::number(0)); for (size_t i = 0; i < std::min(g.size(), n); ++i) gp[1][i] = g[i] ? g[i] : SymbolicExpr::number(0); }
    for (size_t k = 2; k < n; ++k) {
        auto step = context.consume_steps(1, operation);
        if (!step) return PowerSeriesResult::failure(step.error());
        auto multiplied = power_series_multiply_checked(gp[k-1], gp[1], order, context);
        if (!multiplied) return multiplied;
        gp[k] = std::move(multiplied.value());
    }
    for (size_t k = 0; k < std::min(f.size(), n); ++k) {
        auto step = context.consume_steps(1, operation);
        if (!step) return PowerSeriesResult::failure(step.error());
        auto fk = f[k] ? f[k] : SymbolicExpr::number(0);
        if (fk->is_zero()) continue;
        for (size_t i = 0; i < n && i < gp[k].size(); ++i)
            result[i] = SymbolicExpr::add(result[i], SymbolicExpr::multiply(fk, gp[k][i]))->simplify();
    }
    return PowerSeriesResult::success(std::move(result));
}

PowerSeriesResult power_series_compose_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& g, int order) {
    ComputationContext context;
    return power_series_compose_checked(f, g, order, context);
}

std::vector<std::shared_ptr<SymbolicExpr>> power_series_compose(
    const std::vector<std::shared_ptr<SymbolicExpr>>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& g, int order) {
    if (order <= 0) return {};
    size_t n = static_cast<size_t>(order);
    if (!g.empty() && g[0] && !g[0]->is_zero()) return {};
    std::vector<std::shared_ptr<SymbolicExpr>> result(n, SymbolicExpr::number(0));
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> gp(n);
    gp[0].resize(n, SymbolicExpr::number(0)); gp[0][0] = SymbolicExpr::number(1);
    if (n > 1) { gp[1].resize(n, SymbolicExpr::number(0)); for (size_t i = 0; i < std::min(g.size(), n); ++i) gp[1][i] = g[i] ? g[i] : SymbolicExpr::number(0); }
    for (size_t k = 2; k < n; ++k) gp[k] = power_series_multiply(gp[k-1], gp[1], order);
    for (size_t k = 0; k < std::min(f.size(), n); ++k) {
        auto fk = f[k] ? f[k] : SymbolicExpr::number(0);
        if (fk->is_zero()) continue;
        for (size_t i = 0; i < n && i < gp[k].size(); ++i)
            result[i] = SymbolicExpr::add(result[i], SymbolicExpr::multiply(fk, gp[k][i]))->simplify();
    }
    return result;
}

std::shared_ptr<SymbolicExpr> fourier_series(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& period, int n_terms) {
    if (!f || !period || n_terms < 0) return nullptr;

    /// 周期 T，半周期 L = T/2，基频 w = 2π/T
    auto T = period;
    auto pi = lamina::detail::make_expression_ptr(lamina::detail::make_node<VariableNode>("pi"));
    auto two = SymbolicExpr::number(2);
    auto L = SymbolicExpr::divide(T, two);          // 半周期
    auto half_lo = SymbolicExpr::multiply(SymbolicExpr::number(-1), L);
    auto half_hi = L;
    auto w = SymbolicExpr::divide(SymbolicExpr::multiply(two, pi), T); // 2π/T

    int parity = detect_parity(f, var); // 1=even, -1=odd, 0=unknown

    Integrator integrator;
    auto x = SymbolicExpr::variable(var);

    /// a0 = (1/L) ∫_{-L}^{L} f dx ; 常数项 a0/2
    std::shared_ptr<SymbolicExpr> a0 = SymbolicExpr::number(0);
    {
        auto integ = integrator.integrate_def(*f, var, *half_lo, *half_hi);
        a0 = SymbolicExpr::divide(lamina::detail::make_expression_ptr(lamina::detail::node(integ)), L)->simplify();
    }
    auto result = SymbolicExpr::divide(a0, two); // a0/2

    for (int k = 1; k <= n_terms; ++k) {
        auto kw = SymbolicExpr::multiply(SymbolicExpr::number(k), w); // k·2π/T
        auto arg = SymbolicExpr::multiply(kw, x);

        /// a_k：奇函数时为 0
        if (parity != -1) {
            auto integrand = SymbolicExpr::multiply(f, SymbolicExpr::cos(arg));
            auto integ = integrator.integrate_def(*integrand, var, *half_lo, *half_hi);
            auto ak = SymbolicExpr::divide(lamina::detail::make_expression_ptr(lamina::detail::node(integ)), L)->simplify();
            if (!(lamina::detail::node(ak) && lamina::detail::node(ak)->is_zero())) {
                result = SymbolicExpr::add(result, SymbolicExpr::multiply(ak, SymbolicExpr::cos(arg)));
            }
        }
        /// b_k：偶函数时为 0
        if (parity != 1) {
            auto integrand = SymbolicExpr::multiply(f, SymbolicExpr::sin(arg));
            auto integ = integrator.integrate_def(*integrand, var, *half_lo, *half_hi);
            auto bk = SymbolicExpr::divide(lamina::detail::make_expression_ptr(lamina::detail::node(integ)), L)->simplify();
            if (!(lamina::detail::node(bk) && lamina::detail::node(bk)->is_zero())) {
                result = SymbolicExpr::add(result, SymbolicExpr::multiply(bk, SymbolicExpr::sin(arg)));
            }
        }
    }
    return result->simplify();
}

std::shared_ptr<SymbolicExpr> laurent_series(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center, int order_neg, int order_pos) {
    if (f && lamina::detail::node(f) && center && center->is_zero()) {
        if (auto power = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(f))) {
            auto base_var = std::dynamic_pointer_cast<const VariableNode>(power->base());
            auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
            if (base_var && base_var->name() == var && exponent) {
                double exponent_value = 0.0;
                if (std::holds_alternative<BigInt>(exponent->value())) {
                    exponent_value = std::get<BigInt>(exponent->value()).to_double();
                } else if (std::holds_alternative<Rational>(exponent->value())) {
                    exponent_value = std::get<Rational>(exponent->value()).to_double();
                } else {
                    exponent_value = static_cast<double>(std::get<lmmc_real_t>(exponent->value()));
                }
                if (exponent_value < 0.0 && exponent_value == std::floor(exponent_value)) {
                    return f->simplify();
                }
            }
        }
    }
    auto full = laurent_series_full(f, var, center, order_neg, order_pos);
    return full.series;
}

ExpressionResult laurent_series_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos,
    ComputationContext& context)
{
    auto full = laurent_series_full_checked(f, var, center, order_neg, order_pos, context);
    if (!full) return ExpressionResult::failure(full.error());
    return ExpressionResult::success(full.value().series);
}

ExpressionResult laurent_series_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos)
{
    ComputationContext context;
    return laurent_series_checked(f, var, center, order_neg, order_pos, context);
}

LaurentSeriesResult laurent_series_full_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos,
    ComputationContext& context)
{
    const std::string operation = "laurent_series_full";
    auto input = validate_series_expr_point(f, center, var, context, operation);
    if (!input) return LaurentSeriesResult::failure(input.error());
    auto order_check = validate_laurent_orders(order_neg, order_pos, operation);
    if (!order_check) return LaurentSeriesResult::failure(order_check.error());
    auto budget = context.consume_steps(static_cast<std::size_t>(order_neg + order_pos + 1) * 16 + 16,
                                        operation);
    if (!budget) return LaurentSeriesResult::failure(budget.error());
    if (!center->is_zero()) {
        return LaurentSeriesResult::failure(
            CasErrc::Inconclusive,
            "checked Laurent support is currently limited to expansions around zero",
            operation);
    }
    auto simplified_input = f->simplify();
    if (!simplified_input || !lamina::detail::node(simplified_input)) {
        return LaurentSeriesResult::failure(
            CasErrc::Inconclusive,
            "checked Laurent input could not be simplified in the supported domain",
            operation);
    }
    auto supported_power = supported_laurent_integer_power(lamina::detail::node(simplified_input), var);
    if (!supported_power) {
        return LaurentSeriesResult::failure(
            CasErrc::Inconclusive,
            "checked Laurent support is currently limited to integer-power monomials",
            operation);
    }
    try {
        auto result = laurent_series_full(f, var, center, order_neg, order_pos);
        if (!result.series || !lamina::detail::node(result.series)) {
            return LaurentSeriesResult::failure(
                CasErrc::Inconclusive,
                "Laurent series could not be constructed in the supported domain",
                operation);
        }
        if (!result.residue || !lamina::detail::node(result.residue)) {
            return LaurentSeriesResult::failure(
                CasErrc::InternalInvariant,
                "Laurent series produced a null residue expression",
                operation);
        }
        int expected_pole_order = *supported_power < 0 ? -*supported_power : 0;
        if (result.pole_order != expected_pole_order) {
            return LaurentSeriesResult::failure(
                CasErrc::Inconclusive,
                "Laurent pole order could not be verified in the supported domain",
                operation);
        }
        return LaurentSeriesResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return LaurentSeriesResult::failure(CasErrc::ResourceLimit,
                                            "allocation failed while calculating Laurent series",
                                            operation);
    } catch (const std::exception& ex) {
        return LaurentSeriesResult::failure(CasErrc::InternalInvariant,
                                            ex.what(),
                                            operation);
    }
}

LaurentSeriesResult laurent_series_full_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos)
{
    ComputationContext context;
    return laurent_series_full_checked(f, var, center, order_neg, order_pos, context);
}

LaurentResult laurent_series_full(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center, int order_neg, int order_pos) {
    LaurentResult res{nullptr, SingularityType::Removable, 0, SymbolicExpr::number(0)};
    if (!f || !center) return res;

    auto x = SymbolicExpr::variable(var);
    /// (x - center)
    auto shift = SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), center));

    /// 寻找极点阶数 m：使 (x-c)^m · f 在 c 处解析（有限）。
    int pole_order = 0;
    std::shared_ptr<SymbolicExpr> regular = f;
    {
        auto test = f;
        const int MAX_M = (order_neg > 0 ? order_neg : 5);
        /// 先判断 f 在 center 是否已解析
        auto f_at = f->substitute(var, center);
        if (f_at) f_at = f_at->simplify();
        bool finite = f_at && f_at->is_number();
        if (!finite) {
            for (int m = 1; m <= MAX_M; ++m) {
                auto powm = SymbolicExpr::power(shift, SymbolicExpr::number(m));
                auto g = SymbolicExpr::multiply(powm, f)->simplify();
                auto g_at = g->substitute(var, center);
                if (g_at) g_at = g_at->simplify();
                if (g_at && g_at->is_number()) {
                    pole_order = m;
                    regular = g;       // (x-c)^m f(x)，在 c 处解析
                    break;
                }
            }
        }
    }

    /// 对 regular（解析部分）做 Taylor 展开
    int taylor_order = order_pos + pole_order + 1;
    if (taylor_order < 1) taylor_order = 1;
    auto taylor = regular->series(var, center, taylor_order);
    if (!taylor) return res;

    /// Laurent = taylor / (x-c)^pole_order
    std::shared_ptr<SymbolicExpr> laurent = taylor;
    if (pole_order > 0) {
        auto powm = SymbolicExpr::power(shift, SymbolicExpr::number(pole_order));
        laurent = SymbolicExpr::divide(taylor, powm)->simplify();
    }

    res.series = laurent->simplify();
    res.pole_order = pole_order;
    if (pole_order == 0) res.singularity = SingularityType::Removable;
    else res.singularity = SingularityType::Pole;

    /// 留数 = Taylor 展开中 (x-c)^(pole_order-1) 的系数
    if (pole_order >= 1) {
        /// a_{m-1} = (1/(m-1)!) lim_{x->c} d^{m-1}/dx^{m-1} [regular]
        auto deriv = regular;
        for (int i = 0; i < pole_order - 1; ++i) deriv = deriv->differentiate(var);
        /// 用极限求值，稳健处理 0/0 形式（如 z/(z(z+1)) 在 z=0）。
        auto val = deriv->limit(var, center);
        if (!val) val = deriv->substitute(var, center);
        long long fact = 1;
        for (int i = 2; i <= pole_order - 1; ++i) fact *= i;
        res.residue = SymbolicExpr::divide(val, SymbolicExpr::number((int)fact))->simplify();
    }

    return res;
}

std::shared_ptr<SymbolicExpr> asymptotic_expand(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var, int order) {
    if (!f || order < 0) return nullptr;
    /// 通过 x = 1/t 替换，在 t=0 处做 Taylor 展开，再回代 t = 1/x，
    /// 得到按 x 递减幂次的渐近展开。
    auto x = SymbolicExpr::variable(var);
    std::string tname = var + "__asym_t";
    auto t = SymbolicExpr::variable(tname);
    auto inv_t = SymbolicExpr::divide(SymbolicExpr::number(1), t);
    auto g = f->substitute(var, inv_t);
    if (!g) return nullptr;
    g = g->simplify();
    auto taylor_t = g->series(tname, SymbolicExpr::number(0), order + 1);
    if (!taylor_t) return nullptr;
    /// 回代 t = 1/x
    auto inv_x = SymbolicExpr::divide(SymbolicExpr::number(1), x);
    auto back = taylor_t->substitute(tname, inv_x);
    if (!back) return nullptr;
    return back->simplify();
}

std::shared_ptr<SymbolicExpr> symbolic_sum(
    const std::shared_ptr<SymbolicExpr>& body, const std::string& index,
    const std::shared_ptr<SymbolicExpr>& lower, const std::shared_ptr<SymbolicExpr>& upper) {
    if (!body || !lower || !upper) return nullptr;
    long long lv = 0, uv = 0;
    if (try_get_int(lower, lv) && try_get_int(upper, uv)) {
        if (uv < lv) return SymbolicExpr::number(0);
        if (uv - lv < 100) {
            std::vector<std::shared_ptr<const SymbolicNode>> terms;
            for (long long k = lv; k <= uv; ++k) {
                auto val = SymbolicExpr::number(static_cast<int>(k));
                auto term = body->substitute(index, val);
                if (term) { term = term->simplify(); terms.push_back(lamina::detail::node(term)); }
            }
            if (terms.empty()) return SymbolicExpr::number(0);
            if (terms.size() == 1) return lamina::detail::make_expression_ptr(terms[0])->simplify();
            return lamina::detail::make_expression_ptr(lamina::detail::make_node<AddNode>(terms))->simplify();
        }
    }
    auto k_var = SymbolicExpr::variable(index);
    if (body->to_string() == k_var->to_string()) {
        auto u = upper; auto l = lower;
        auto u1 = SymbolicExpr::add(u, SymbolicExpr::number(1));
        auto su = SymbolicExpr::divide(SymbolicExpr::multiply(u, u1), SymbolicExpr::number(2));
        auto lm1 = SymbolicExpr::add(l, SymbolicExpr::number(-1));
        auto sl = SymbolicExpr::divide(SymbolicExpr::multiply(lm1, l), SymbolicExpr::number(2));
        return SymbolicExpr::add(su, SymbolicExpr::multiply(SymbolicExpr::number(-1), sl))->simplify();
    }
    return lamina::detail::make_expression_ptr(lamina::detail::make_node<SummationNode>(lamina::detail::node(body), index, lamina::detail::node(lower), lamina::detail::node(upper)));
}

std::shared_ptr<SymbolicExpr> symbolic_product(
    const std::shared_ptr<SymbolicExpr>& body, const std::string& index,
    const std::shared_ptr<SymbolicExpr>& lower, const std::shared_ptr<SymbolicExpr>& upper) {
    if (!body || !lower || !upper) return nullptr;
    long long lv = 0, uv = 0;
    if (try_get_int(lower, lv) && try_get_int(upper, uv)) {
        if (uv < lv) return SymbolicExpr::number(1);
        if (uv - lv < 50) {
            auto result = SymbolicExpr::number(1);
            for (long long k = lv; k <= uv; ++k) {
                auto val = SymbolicExpr::number(static_cast<int>(k));
                auto term = body->substitute(index, val);
                if (term) { term = term->simplify(); result = SymbolicExpr::multiply(result, term); }
            }
            return result->simplify();
        }
    }
    return lamina::detail::make_expression_ptr(lamina::detail::make_node<ProductNode>(lamina::detail::node(body), index, lamina::detail::node(lower), lamina::detail::node(upper)));
}


std::shared_ptr<SymbolicExpr> lim_sup(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n) {
    if (!a_n || !lamina::detail::node(a_n)) return nullptr;

    /// Detect (-1)^n oscillation pattern
    std::shared_ptr<SymbolicExpr> remainder;
    if (series_extract_alternating(lamina::detail::node(a_n), n, remainder)) {
        auto inf = SymbolicExpr::infinity();
        auto lim_f = remainder->limit(n, inf);
        if (lim_f) {
            auto ls = lim_f->simplify();
            if (ls && !series_is_infinity(ls)) return series_abs(ls);
            if (series_is_infinity(ls)) return SymbolicExpr::infinity();
        }
        return series_abs(lim_f);
    }

    /// Detect sin/cos oscillation
    std::shared_ptr<SymbolicExpr> amplitude;
    if (series_detect_trig_oscillation(lamina::detail::node(a_n), n, amplitude)) {
        auto inf = SymbolicExpr::infinity();
        auto amp_lim = amplitude->limit(n, inf);
        if (amp_lim) {
            auto als = amp_lim->simplify();
            if (als && series_is_number(als)) return SymbolicExpr::number(std::abs(series_get_double(als)));
            if (als) return series_abs(als);
        }
    }

    /// For convergent/monotone sequences: lim sup = lim
    auto inf = SymbolicExpr::infinity();
    auto lim = a_n->limit(n, inf);
    if (lim) { auto s = lim->simplify(); if (s) return s; return lim; }
    return nullptr;
}

std::shared_ptr<SymbolicExpr> lim_inf(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n) {
    if (!a_n || !lamina::detail::node(a_n)) return nullptr;

    /// Detect (-1)^n oscillation pattern
    std::shared_ptr<SymbolicExpr> remainder;
    if (series_extract_alternating(lamina::detail::node(a_n), n, remainder)) {
        auto inf = SymbolicExpr::infinity();
        auto lim_f = remainder->limit(n, inf);
        if (lim_f) {
            auto ls = lim_f->simplify();
            if (ls && !series_is_infinity(ls)) return series_negate(series_abs(ls));
            if (series_is_infinity(ls)) return SymbolicExpr::infinity(-1);
        }
        return series_negate(series_abs(lim_f));
    }

    /// Detect sin/cos oscillation
    std::shared_ptr<SymbolicExpr> amplitude;
    if (series_detect_trig_oscillation(lamina::detail::node(a_n), n, amplitude)) {
        auto inf = SymbolicExpr::infinity();
        auto amp_lim = amplitude->limit(n, inf);
        if (amp_lim) {
            auto als = amp_lim->simplify();
            if (als && series_is_number(als)) return SymbolicExpr::number(-std::abs(series_get_double(als)));
            if (als) return series_negate(series_abs(als));
        }
    }

    /// For convergent/monotone sequences: lim inf = lim
    auto inf = SymbolicExpr::infinity();
    auto lim = a_n->limit(n, inf);
    if (lim) { auto s = lim->simplify(); if (s) return s; return lim; }
    return nullptr;
}

} // namespace lamina
