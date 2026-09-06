/**
 * @file symbolic_ode_engine.cpp
 * @brief 统一 ODE 求解引擎实现：类型检测与分类。
 */
#include "../include/symbolic_ode_engine.hpp"
#include "symbolic_ast.hpp"
#include "../include/symbolic.hpp"
#include "../include/poly_utils.hpp"
#include "../include/residual_verification.hpp"
#include "internal/expression_analysis.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "internal/ode_support.hpp"
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <limits>
#include <map>
#include <optional>
#include <vector>

namespace LMCAS {


/**
 * @internal
 * @brief 尝试将表达式求值为 double。
 * @return 若表达式为纯数值返回其值，否则返回 NaN。
 */
double try_eval_double(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) return std::numeric_limits<double>::quiet_NaN();
    if (expr->is_number()) {
        auto node = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr));
        if (!node) return std::numeric_limits<double>::quiet_NaN();
        if (std::holds_alternative<BigInt>(node->value()))
            return std::get<BigInt>(node->value()).to_double();
        if (std::holds_alternative<Rational>(node->value()))
            return std::get<Rational>(node->value()).to_double();
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
    if (!step) return step;
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
    if (!step) return step;
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
    if (!step) return step;
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
    if (!step) return step;
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
    if (!step) return step;
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

/**
 * @internal
 * @brief 检查表达式是否为纯数值常量（不依赖任何变量）。
 */
[[maybe_unused]] static bool is_constant_expr(const std::shared_ptr<SymbolicExpr>& expr,
                                              const std::string& x,
                                              const std::string& y) {
    if (!expr || !LMCAS::detail::node(expr)) return true;
    return !expression_depends_on_variable(LMCAS::detail::node(expr), x) && !expression_depends_on_variable(LMCAS::detail::node(expr), y);
}

/**
 * @internal
 * @brief 检查表达式是否仅依赖指定变量（不依赖另一个变量）。
 */
[[maybe_unused]] static bool depends_only_on(const std::shared_ptr<SymbolicExpr>& expr,
                                             const std::string&,
                                             const std::string& other_var) {
    if (!expr || !LMCAS::detail::node(expr)) return true;
    return !expression_depends_on_variable(LMCAS::detail::node(expr), other_var);
}

static bool valid_classifier_variables(
    const std::string& x,
    const std::string& y)
{
    return !x.empty() && !y.empty() && x != y;
}


bool is_separable(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y)
{
    if (!rhs || !LMCAS::detail::node(rhs)) return false;
    if (!valid_classifier_variables(x, y)) return false;

    bool has_x = expression_depends_on_variable(LMCAS::detail::node(rhs), x);
    bool has_y = expression_depends_on_variable(LMCAS::detail::node(rhs), y);

    /// 若只依赖一个变量或都不依赖，则可分离
    if (!has_x || !has_y) return true;

    /// 检查乘法形式 f(x)*g(y)
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(rhs));
    if (mul) {
        /// 将因子分为仅含 x 的和仅含 y 的
        bool all_separable = true;
        for (const auto& factor : mul->operands()) {
            bool fx = expression_depends_on_variable(factor, x);
            bool fy = expression_depends_on_variable(factor, y);
            if (fx && fy) {
                all_separable = false;
                break;
            }
        }
        if (all_separable) return true;
    }

    /// 检查除法形式 f(x)/g(y) 或 g(y)/f(x)
    auto pow = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(rhs));
    if (pow) {
        auto exp_node = std::dynamic_pointer_cast<const NumberNode>(pow->exponent());
        if (exp_node) {
            double exp_val = try_eval_double(LMCAS::detail::make_expression_ptr(pow->exponent()));
            int eq;
            lmmc_double_nearly_equal_tol(exp_val, -1.0, 1e-12, 1e-12, &eq);
            if (eq) {
                /// rhs = base^(-1) = 1/base; 若 base 仅含一个变量则可分离
                bool base_x = expression_depends_on_variable(pow->base(), x);
                bool base_y = expression_depends_on_variable(pow->base(), y);
                if ((base_x && !base_y) || (!base_x && base_y)) return true;
            }
        }
    }

    return false;
}


bool is_linear_first_order(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y,
    std::shared_ptr<SymbolicExpr>& P,
    std::shared_ptr<SymbolicExpr>& Q)
{
    P.reset();
    Q.reset();
    if (!rhs || !LMCAS::detail::node(rhs) ||
        !valid_classifier_variables(x, y)) {
        return false;
    }

    /// 方程形式: y' = rhs(x, y)
    /// 线性形式: y' + P(x)*y = Q(x)  →  y' = Q(x) - P(x)*y  →  rhs = Q(x) - P(x)*y
    /// 即 rhs 关于 y 是线性的: rhs = A(x) + B(x)*y，其中 Q = A, P = -B

    /// 若 rhs 不依赖 y，则 P=0, Q=rhs
    if (!expression_depends_on_variable(LMCAS::detail::node(rhs), y)) {
        P = SymbolicExpr::number(0);
        Q = rhs;
        return true;
    }

    /// 对 rhs 关于 y 求导，若结果不依赖 y，则 rhs 关于 y 是线性的
    auto drhs_dy = rhs->differentiate(y);
    if (!drhs_dy || expression_depends_on_variable(LMCAS::detail::node(drhs_dy), y)) {
        return false;
    }

    /// rhs = A(x) + B(x)*y，其中 B(x) = ∂rhs/∂y
    /// 计算 A(x) = rhs|_{y=0}
    auto A = rhs->substitute(y, SymbolicExpr::number(0));
    if (!A) return false;

    /// 验证 A 不依赖 y
    if (expression_depends_on_variable(LMCAS::detail::node(A), y)) return false;

    /// B(x) = drhs_dy（已验证不依赖 y）
    /// 线性形式: y' = A(x) + B(x)*y  →  y' - B(x)*y = A(x)  →  y' + (-B(x))*y = A(x)
    /// 所以 P = -B(x), Q = A(x)
    P = SymbolicExpr::multiply(SymbolicExpr::number(-1), drhs_dy)->simplify();
    Q = A->simplify();
    return true;
}


bool is_homogeneous_ode(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y)
{
    if (!rhs || !LMCAS::detail::node(rhs)) return false;
    if (!valid_classifier_variables(x, y)) return false;

    /// 齐次方程: f(tx, ty) = f(x, y) 对所有 t 成立
    /// 用 t=2 进行数值测试：f(2x, 2y) 应等于 f(x, y)
    auto t_val = SymbolicExpr::number(2);
    auto tx = SymbolicExpr::multiply(t_val, SymbolicExpr::variable(x));
    auto ty = SymbolicExpr::multiply(t_val, SymbolicExpr::variable(y));

    auto f_scaled = rhs->substitute(x, tx);
    f_scaled = f_scaled->substitute(y, ty);
    f_scaled = f_scaled->simplify();

    /// 计算 f_scaled - rhs，若化简为零则齐次
    auto diff = SymbolicExpr::add(f_scaled,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), rhs));
    diff = diff->simplify();

    if (diff->is_zero()) return true;

    ComputationContext context;
    auto verified = check_zero_residual(diff, context);
    return verified &&
        std::holds_alternative<ProvedZeroResidual>(verified.value());
}


static std::optional<int> bernoulli_integer_exponent(const NumberNode& number)
{
    std::optional<std::int64_t> value;
    if (std::holds_alternative<BigInt>(number.value())) {
        value = std::get<BigInt>(number.value()).try_to_int64();
    } else if (std::holds_alternative<Rational>(number.value())) {
        const auto& rational = std::get<Rational>(number.value());
        if (!rational.is_integer()) return std::nullopt;
        value = rational.to_BigInt().try_to_int64();
    } else {
        const lmmc_real_t approximate =
            std::get<lmmc_real_t>(number.value());
        if (!std::isfinite(approximate) ||
            approximate != std::floor(approximate) ||
            approximate < std::numeric_limits<int>::min() ||
            approximate > std::numeric_limits<int>::max()) {
            return std::nullopt;
        }
        return static_cast<int>(approximate);
    }
    if (!value || *value < std::numeric_limits<int>::min() ||
        *value > std::numeric_limits<int>::max()) {
        return std::nullopt;
    }
    return static_cast<int>(*value);
}

static bool extract_bernoulli_monomial(
    const std::shared_ptr<const SymbolicNode>& term,
    const std::string& y,
    int& exponent,
    std::shared_ptr<SymbolicExpr>& coefficient)
{
    exponent = 0;
    coefficient = SymbolicExpr::number(1);
    std::vector<std::shared_ptr<const SymbolicNode>> factors;
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(term)) {
        factors = product->operands();
    } else {
        factors.push_back(term);
    }

    for (const auto& factor : factors) {
        if (auto variable =
                std::dynamic_pointer_cast<const VariableNode>(factor);
            variable && variable->name() == y) {
            if (exponent == std::numeric_limits<int>::max()) return false;
            ++exponent;
            continue;
        }

        if (auto power = std::dynamic_pointer_cast<const PowerNode>(factor)) {
            auto base =
                std::dynamic_pointer_cast<const VariableNode>(power->base());
            auto power_value =
                std::dynamic_pointer_cast<const NumberNode>(power->exponent());
            if (base && base->name() == y && power_value) {
                auto value = bernoulli_integer_exponent(*power_value);
                if (!value ||
                    (*value > 0 &&
                     exponent > std::numeric_limits<int>::max() - *value) ||
                    (*value < 0 &&
                     exponent < std::numeric_limits<int>::min() - *value)) {
                    return false;
                }
                exponent += *value;
                continue;
            }
        }

        if (expression_depends_on_variable(factor, y)) return false;
        coefficient = SymbolicExpr::multiply(
            coefficient, LMCAS::detail::make_expression_ptr(factor))->simplify();
    }
    return true;
}

bool is_bernoulli_ode(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y,
    std::shared_ptr<SymbolicExpr>& P,
    std::shared_ptr<SymbolicExpr>& Q,
    int& n)
{
    P.reset();
    Q.reset();
    n = 0;
    if (!rhs || !LMCAS::detail::node(rhs) ||
        !valid_classifier_variables(x, y)) return false;
    if (!expression_depends_on_variable(LMCAS::detail::node(rhs), y)) {
        return false;
    }

    std::vector<std::shared_ptr<const SymbolicNode>> terms;
    if (auto sum =
            std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(rhs))) {
        terms = sum->operands();
    } else {
        terms.push_back(LMCAS::detail::node(rhs));
    }

    std::map<int, std::shared_ptr<SymbolicExpr>> coefficients;
    for (const auto& term : terms) {
        int exponent = 0;
        std::shared_ptr<SymbolicExpr> coefficient;
        if (!extract_bernoulli_monomial(
                term, y, exponent, coefficient)) {
            return false;
        }
        auto& combined = coefficients[exponent];
        combined = combined
            ? SymbolicExpr::add(combined, coefficient)->simplify()
            : coefficient->simplify();
    }

    auto linear_coefficient = SymbolicExpr::number(0);
    std::shared_ptr<SymbolicExpr> nonlinear_coefficient;
    int nonlinear_exponent = 0;
    for (auto& [exponent, coefficient] : coefficients) {
        coefficient = coefficient->simplify();
        if (coefficient->is_zero()) continue;
        if (exponent == 1) {
            linear_coefficient = coefficient;
            continue;
        }
        if (exponent == 0 || nonlinear_coefficient) return false;
        nonlinear_exponent = exponent;
        nonlinear_coefficient = coefficient;
    }

    if (!nonlinear_coefficient ||
        nonlinear_exponent == 0 || nonlinear_exponent == 1) {
        return false;
    }

    P = SymbolicExpr::multiply(
        SymbolicExpr::number(-1), linear_coefficient)->simplify();
    Q = nonlinear_coefficient->simplify();
    n = nonlinear_exponent;
    return true;
}


bool is_exact_ode(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x,
    const std::string& y)
{
    if (!M || !LMCAS::detail::node(M) ||
        !N || !LMCAS::detail::node(N) ||
        !valid_classifier_variables(x, y)) return false;

    /// 恰当条件: ∂M/∂y = ∂N/∂x
    auto dM_dy = M->differentiate(y);
    auto dN_dx = N->differentiate(x);

    if (!dM_dy || !dN_dx) return false;

    dM_dy = dM_dy->simplify();
    dN_dx = dN_dx->simplify();
    if (dM_dy->compare(dN_dx) == 0) return true;

    /// 计算差值并化简
    auto diff = SymbolicExpr::add(dM_dy,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dN_dx));
    diff = diff->simplify();

    if (diff->is_zero()) return true;

    ComputationContext context;
    auto verified = check_zero_residual(diff, context);
    return verified &&
        std::holds_alternative<ProvedZeroResidual>(verified.value());

}


bool is_constant_coefficient(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coeffs,
    const std::string& x)
{
    if (coeffs.empty() || x.empty()) return false;
    for (const auto& c : coeffs) {
        if (!c || !LMCAS::detail::node(c)) return false;
        if (expression_depends_on_variable(LMCAS::detail::node(c), x)) return false;
    }
    return true;
}


bool is_euler_equation(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coeffs,
    const std::string& x,
    std::vector<double>& euler_consts)
{
    euler_consts.clear();
    if (coeffs.empty() || x.empty()) return false;
    if (!coeffs.front() || !LMCAS::detail::node(coeffs.front()) ||
        coeffs.front()->is_zero()) {
        return false;
    }

    const int highest_order = static_cast<int>(coeffs.size()) - 1;
    std::vector<double> extracted(coeffs.size(), 0.0);
    for (int index = 0; index <= highest_order; ++index) {
        const int derivative_order = highest_order - index;
        if (!coeffs[index] || !LMCAS::detail::node(coeffs[index])) {
            return false;
        }
        if (coeffs[index]->is_zero()) continue;

        std::shared_ptr<SymbolicExpr> ratio;
        if (derivative_order == 0) {
            ratio = coeffs[index]->simplify();
        } else {
            auto x_power = SymbolicExpr::power(
                SymbolicExpr::variable(x),
                SymbolicExpr::number(derivative_order));
            ratio = SymbolicExpr::divide(
                coeffs[index], x_power)->simplify();
        }
        if (!ratio || !LMCAS::detail::node(ratio)) return false;

        if (expression_depends_on_variable(
                LMCAS::detail::node(ratio), x)) {
            return false;
        }

        const double value = try_eval_double(ratio);
        if (!std::isfinite(value)) return false;
        extracted[index] = value;
    }

    euler_consts = std::move(extracted);
    return true;
}


ODEClassification classify_first_order_ode(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y)
{
    ODEClassification result;
    result.order = 1;

    if (!rhs || !LMCAS::detail::node(rhs)) {
        result.type = ODEType::Unknown;
        return result;
    }

    /// 1. 检测可分离变量
    if (is_separable(rhs, x, y)) {
        result.type = ODEType::Separable;
        return result;
    }

    /// 2. 检测一阶线性
    std::shared_ptr<SymbolicExpr> P, Q;
    if (is_linear_first_order(rhs, x, y, P, Q)) {
        result.type = ODEType::Linear1;
        result.P_coeff = P;
        result.Q_coeff = Q;
        return result;
    }

    /// 3. 检测齐次方程
    if (is_homogeneous_ode(rhs, x, y)) {
        result.type = ODEType::Homogeneous;
        return result;
    }

    /// 4. 检测 Bernoulli
    std::shared_ptr<SymbolicExpr> bP, bQ;
    int bn;
    if (is_bernoulli_ode(rhs, x, y, bP, bQ, bn)) {
        result.type = ODEType::Bernoulli;
        result.bernoulli_P = bP;
        result.bernoulli_Q = bQ;
        result.bernoulli_n = bn;
        return result;
    }

    /// 5. 检测恰当方程
    /// 将 y' = rhs 改写为 M + N*y' = 0 形式: -rhs + y' = 0
    /// 即 M = -rhs, N = 1（标准形式 M*dx + N*dy = 0 → M + N*dy/dx = 0）
    /// 更一般地: 若 rhs = f(x,y)，则 -f(x,y)*dx + dy = 0
    /// M = -f(x,y), N = 1
    auto neg_rhs = SymbolicExpr::multiply(SymbolicExpr::number(-1), rhs)->simplify();
    auto one = SymbolicExpr::number(1);
    if (is_exact_ode(neg_rhs, one, x, y)) {
        result.type = ODEType::Exact;
        result.exact_M = neg_rhs;
        result.exact_N = one;
        return result;
    }

    result.type = ODEType::Unknown;
    return result;
}


ODEClassification classify_higher_order_ode(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string&)
{
    ODEClassification result;
    result.order = static_cast<int>(coeffs.size()) - 1;

    if (coeffs.empty()) {
        result.type = ODEType::Unknown;
        return result;
    }

    /// 1. 检测当前数值分类结果能够表示的常系数
    if (is_constant_coefficient(coeffs, x)) {
        std::vector<double> numeric_coeffs;
        numeric_coeffs.reserve(coeffs.size());
        for (const auto& c : coeffs) {
            double val = try_eval_double(c);
            if (!std::isfinite(val)) {
                result.type = ODEType::Unknown;
                return result;
            }
            numeric_coeffs.push_back(val);
        }

        result.type = result.order == 2
            ? ODEType::Linear2_ConstCoeff
            : ODEType::HigherOrder_ConstCoeff;
        result.const_coeffs = std::move(numeric_coeffs);
        result.forcing_func = forcing;
        return result;
    }

    /// 2. 检测 Euler 方程
    std::vector<double> euler_consts;
    if (is_euler_equation(coeffs, x, euler_consts)) {
        result.type = ODEType::Euler;
        result.euler_coeffs = euler_consts;
        result.euler_forcing = forcing;
        return result;
    }

    result.type = ODEType::Unknown;
    return result;
}

} // namespace LMCAS
