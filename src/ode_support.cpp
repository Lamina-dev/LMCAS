/**
 * @file ode_support.cpp
 * @brief ODE 输入检查与结果处理的共享约定。
 */
#include "symbolic_ode_engine.hpp"
#include "internal/symbolic_ast.hpp"
#include "symbolic.hpp"
#include "lmmc/config.h"
#include "internal/ode_support.hpp"
#include <cmath>
#include <memory>
#include <string>
#include <limits>
#include <vector>

namespace LMCAS {


/**
 * @internal
 * @brief 尝试将表达式求值为 double。
 * @return 若表达式为纯数值返回其值，否则返回 NaN。
 */
double try_eval_double(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (expr->is_number()) {
        auto node = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr));
        if (!node) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (std::holds_alternative<BigInt>(node->value())) {
            return std::get<BigInt>(node->value()).to_double();
        }
        if (std::holds_alternative<Rational>(node->value())) {
            return std::get<Rational>(node->value()).to_double();
        }
        return static_cast<double>(std::get<lmmc_real_t>(node->value()));
    }
    return std::numeric_limits<double>::quiet_NaN();
}

Result<void> validate_ode_expr_var_pair(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& x,
    const std::string& y,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) {
        return step;
    }
    if (!expr || !LMCAS::detail::node(expr)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "ODE expression cannot be null",
                                     operation);
    }
    if (x.empty() || y.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "ODE variable names cannot be empty",
                                     operation);
    }
    if (x == y) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "ODE independent and dependent variables must be distinct",
                                     operation);
    }
    return Result<void>::success();
}

Result<void> validate_ode_pair_var_pair(
    const std::shared_ptr<SymbolicExpr>& first,
    const std::shared_ptr<SymbolicExpr>& second,
    const std::string& x,
    const std::string& y,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) {
        return step;
    }
    if (!first || !LMCAS::detail::node(first) || !second || !LMCAS::detail::node(second)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "ODE expressions cannot be null",
                                     operation);
    }
    if (x.empty() || y.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "ODE variable names cannot be empty",
                                     operation);
    }
    if (x == y) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "ODE independent and dependent variables must be distinct",
                                     operation);
    }
    return Result<void>::success();
}

ODESolutionResult wrap_ode_solution(ODESolution solution,
                                           ODEType expected_method,
                                           const std::string& operation)
{
    if (!solution.general_solution || !LMCAS::detail::node(solution.general_solution)) {
        return ODESolutionResult::failure(
            CasErrc::Inconclusive,
            "ODE solver produced no solution in the supported domain",
            operation);
    }
    if (solution.method_used != expected_method) {
        return ODESolutionResult::failure(
            CasErrc::InternalInvariant,
            "ODE solver reported an unexpected method",
            operation);
    }
    return ODESolutionResult::success(std::move(solution));
}

Result<void> validate_ode_variables(
    const std::string& x,
    const std::string& y,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) {
        return step;
    }
    if (x.empty() || y.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "ODE variable names cannot be empty",
                                     operation);
    }
    if (x == y) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "ODE independent and dependent variables must be distinct",
                                     operation);
    }
    return Result<void>::success();
}

Result<void> validate_numeric_ode_coefficients(
    const std::vector<double>& coeffs,
    std::size_t min_size,
    std::size_t max_size,
    const std::string& operation)
{
    if (coeffs.size() < min_size || coeffs.size() > max_size) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "ODE coefficient list has unsupported size",
                                     operation);
    }
    if (!std::isfinite(coeffs.front()) || coeffs.front() == 0.0) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "ODE leading coefficient must be finite and nonzero",
                                     operation);
    }
    for (double coeff : coeffs) {
        if (!std::isfinite(coeff)) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "ODE coefficients must be finite",
                                         operation);
        }
    }
    return Result<void>::success();
}

Result<void> validate_ode_three_expr_one_var(
    const std::shared_ptr<SymbolicExpr>& first,
    const std::shared_ptr<SymbolicExpr>& second,
    const std::shared_ptr<SymbolicExpr>& third,
    const std::string& x,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) {
        return step;
    }
    if (!first || !LMCAS::detail::node(first) || !second || !LMCAS::detail::node(second) ||
        !third || !LMCAS::detail::node(third)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "ODE expressions cannot be null",
                                     operation);
    }
    if (x.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "ODE variable name cannot be empty",
                                     operation);
    }
    return Result<void>::success();
}

Result<void> validate_ode_two_expr_point(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) {
        return step;
    }
    if (!p || !LMCAS::detail::node(p) || !q || !LMCAS::detail::node(q) || !x0 || !LMCAS::detail::node(x0)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "Frobenius inputs cannot be null",
                                     operation);
    }
    if (x.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "ODE variable name cannot be empty",
                                     operation);
    }
    return Result<void>::success();
}

}
