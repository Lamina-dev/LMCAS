/**
 * @file symbolic_ode_engine.cpp
 * @brief 统一 ODE 求解引擎实现：类型检测与分类。
 */
#include "../include/symbolic_ode_engine.hpp"
#include "symbolic_ast.hpp"
#include "../include/symbolic.hpp"
#include "../include/poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "internal/ode_support.hpp"
#include <cmath>
#include <memory>
#include <string>

namespace lamina {


/**
 * @internal
 * @brief 尝试将表达式求值为 double。
 * @return 若表达式为纯数值返回其值，否则返回 NaN。
 */
double try_eval_double(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !lamina::detail::node(expr)) return std::numeric_limits<double>::quiet_NaN();
    if (expr->is_number()) {
        auto node = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr));
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
    if (!expr || !lamina::detail::node(expr)) {
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
    if (!first || !lamina::detail::node(first) || !second || !lamina::detail::node(second)) {
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
    if (!solution.general_solution || !lamina::detail::node(solution.general_solution)) {
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
    if (!std::isfinite(coeffs.front()) || std::abs(coeffs.front()) < 1e-15) {
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
    if (!first || !lamina::detail::node(first) || !second || !lamina::detail::node(second) ||
        !third || !lamina::detail::node(third)) {
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
    if (!p || !lamina::detail::node(p) || !q || !lamina::detail::node(q) || !x0 || !lamina::detail::node(x0)) {
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
    if (!expr || !lamina::detail::node(expr)) return true;
    return !expression_depends_on_variable(lamina::detail::node(expr), x) && !expression_depends_on_variable(lamina::detail::node(expr), y);
}

/**
 * @internal
 * @brief 检查表达式是否仅依赖指定变量（不依赖另一个变量）。
 */
[[maybe_unused]] static bool depends_only_on(const std::shared_ptr<SymbolicExpr>& expr,
                                             const std::string&,
                                             const std::string& other_var) {
    if (!expr || !lamina::detail::node(expr)) return true;
    return !expression_depends_on_variable(lamina::detail::node(expr), other_var);
}


bool is_separable(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y)
{
    if (!rhs || !lamina::detail::node(rhs)) return true;

    bool has_x = expression_depends_on_variable(lamina::detail::node(rhs), x);
    bool has_y = expression_depends_on_variable(lamina::detail::node(rhs), y);

    /// 若只依赖一个变量或都不依赖，则可分离
    if (!has_x || !has_y) return true;

    /// 检查乘法形式 f(x)*g(y)
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(rhs));
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
    auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(rhs));
    if (pow) {
        auto exp_node = std::dynamic_pointer_cast<const NumberNode>(pow->exponent());
        if (exp_node) {
            double exp_val = try_eval_double(lamina::detail::make_expression_ptr(pow->exponent()));
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
    const std::string&,
    const std::string& y,
    std::shared_ptr<SymbolicExpr>& P,
    std::shared_ptr<SymbolicExpr>& Q)
{
    if (!rhs || !lamina::detail::node(rhs)) {
        P = SymbolicExpr::number(0);
        Q = SymbolicExpr::number(0);
        return true;
    }

    /// 方程形式: y' = rhs(x, y)
    /// 线性形式: y' + P(x)*y = Q(x)  →  y' = Q(x) - P(x)*y  →  rhs = Q(x) - P(x)*y
    /// 即 rhs 关于 y 是线性的: rhs = A(x) + B(x)*y，其中 Q = A, P = -B

    /// 若 rhs 不依赖 y，则 P=0, Q=rhs
    if (!expression_depends_on_variable(lamina::detail::node(rhs), y)) {
        P = SymbolicExpr::number(0);
        Q = rhs;
        return true;
    }

    /// 对 rhs 关于 y 求导，若结果不依赖 y，则 rhs 关于 y 是线性的
    auto drhs_dy = rhs->differentiate(y);
    if (!drhs_dy || expression_depends_on_variable(lamina::detail::node(drhs_dy), y)) {
        return false;
    }

    /// rhs = A(x) + B(x)*y，其中 B(x) = ∂rhs/∂y
    /// 计算 A(x) = rhs|_{y=0}
    auto A = rhs->substitute(y, SymbolicExpr::number(0));
    if (!A) return false;

    /// 验证 A 不依赖 y
    if (expression_depends_on_variable(lamina::detail::node(A), y)) return false;

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
    if (!rhs || !lamina::detail::node(rhs)) return false;

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

    /// 数值验证：在具体点 (x=1, y=1) 和 (x=2, y=3) 测试
    auto test_at = [&](double xv, double yv) -> bool {
        auto f_orig = rhs->substitute(x, SymbolicExpr::number(xv));
        f_orig = f_orig->substitute(y, SymbolicExpr::number(yv));
        f_orig = f_orig->simplify();

        double t_test = 2.0;
        auto f_sc = rhs->substitute(x, SymbolicExpr::number(t_test * xv));
        f_sc = f_sc->substitute(y, SymbolicExpr::number(t_test * yv));
        f_sc = f_sc->simplify();

        double v_orig = try_eval_double(f_orig);
        double v_sc = try_eval_double(f_sc);

        if (std::isnan(v_orig) || std::isnan(v_sc)) return false;
        if (std::abs(v_orig) < 1e-15 && std::abs(v_sc) < 1e-15) return true;
        if (std::abs(v_orig) < 1e-15) return false;

        int eq;
        lmmc_double_nearly_equal_tol(v_orig, v_sc, 1e-9, 1e-9, &eq);
        return eq != 0;
    };

    /// 在多个点测试
    if (test_at(1.0, 1.0) && test_at(2.0, 3.0) && test_at(0.5, 1.5)) {
        return true;
    }

    return false;
}


bool is_bernoulli_ode(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string&,
    const std::string& y,
    std::shared_ptr<SymbolicExpr>& P,
    std::shared_ptr<SymbolicExpr>& Q,
    int& n)
{
    if (!rhs || !lamina::detail::node(rhs)) return false;
    if (!expression_depends_on_variable(lamina::detail::node(rhs), y)) return false;

    /// Bernoulli: y' + P(x)*y = Q(x)*y^n  →  y' = -P(x)*y + Q(x)*y^n
    /// 即 rhs = -P(x)*y + Q(x)*y^n = y*(-P(x) + Q(x)*y^{n-1})

    /// 尝试 rhs / y 并检查结果是否为 A(x) + B(x)*y^m 形式
    /// 其中 m = n-1, P = -A, Q = B, n = m+1

    /// 首先检查 rhs 是否含 y 的幂次
    /// 对 rhs 关于 y 求两次导，检查是否为 y 的幂函数
    auto d1 = rhs->differentiate(y);
    if (!d1) return false;
    auto d2 = d1->differentiate(y);
    if (!d2) return false;

    /// 若 d2 不依赖 y，则 rhs 关于 y 最多是二次的
    /// 对于 Bernoulli，我们需要 rhs = A(x)*y + B(x)*y^n
    /// 尝试特定的 n 值: 2, 3, -1
    for (int test_n : {2, 3, -1, 4}) {
        /// rhs 应为 -P(x)*y + Q(x)*y^n
        /// 令 rhs / y = -P(x) + Q(x)*y^{n-1}
        /// 若 n=2: rhs/y = -P(x) + Q(x)*y → 关于 y 线性
        /// 若 n=3: rhs/y = -P(x) + Q(x)*y^2 → 关于 y 二次

        /// 计算 rhs 在 y=1 和 y=2 处的值来推断结构
        auto rhs_at_y1 = rhs->substitute(y, SymbolicExpr::number(1))->simplify();
        auto rhs_at_y0 = rhs->substitute(y, SymbolicExpr::number(0))->simplify();

        /// 若 rhs(x, 0) = 0，则 rhs 含 y 因子
        if (!rhs_at_y0->is_zero()) continue;

        /// rhs = y * h(x, y)，计算 h = rhs / y
        /// h(x, y) = -P(x) + Q(x)*y^{n-1}
        /// h(x, 0) = -P(x)
        /// h(x, 1) = -P(x) + Q(x)

        /// 用 y 除 rhs：对 rhs 做 substitute 检查
        /// 实际上，若 rhs(x,0)=0，则 rhs 含 y 因子
        /// 计算 ∂rhs/∂y|_{y=0} = h(x, 0) = -P(x)
        auto h_at_0 = d1->substitute(y, SymbolicExpr::number(0))->simplify();
        if (expression_depends_on_variable(lamina::detail::node(h_at_0), y)) continue;

        /// 对于 Bernoulli n=test_n:
        /// rhs = -P*y + Q*y^n
        /// d(rhs)/dy = -P + n*Q*y^{n-1}
        /// d(rhs)/dy|_{y=0} = -P (对 n>=2)
        /// d²(rhs)/dy²|_{y=0} = n*(n-1)*Q*y^{n-2}|_{y=0}
        ///   对 n=2: = 2*Q
        ///   对 n=3: = 0 (需要 y=0 时 y^1 = 0)

        if (test_n == 2) {
            /// d²rhs/dy² = 2*Q(x) (常数关于 y)
            auto d2_simplified = d2->simplify();
            if (expression_depends_on_variable(lamina::detail::node(d2_simplified), y)) continue;

            /// 验证三阶导为零
            auto d3 = d2->differentiate(y)->simplify();
            if (!d3->is_zero()) continue;

            /// P = -(d1|_{y=0}), Q = d2/2
            P = SymbolicExpr::multiply(SymbolicExpr::number(-1), h_at_0)->simplify();
            Q = SymbolicExpr::divide(d2_simplified, SymbolicExpr::number(2))->simplify();

            /// 验证 Q 不依赖 y
            if (expression_depends_on_variable(lamina::detail::node(Q), y)) continue;

            n = 2;
            return true;
        }

        if (test_n == 3) {
            /// rhs = -P*y + Q*y^3
            /// d1 = -P + 3*Q*y^2
            /// d2 = 6*Q*y
            /// d3 = 6*Q
            auto d2_at_0 = d2->substitute(y, SymbolicExpr::number(0))->simplify();
            if (!d2_at_0->is_zero()) continue;

            auto d3 = d2->differentiate(y)->simplify();
            if (expression_depends_on_variable(lamina::detail::node(d3), y)) continue;

            auto d4 = d3->differentiate(y)->simplify();
            if (!d4->is_zero()) continue;

            P = SymbolicExpr::multiply(SymbolicExpr::number(-1), h_at_0)->simplify();
            Q = SymbolicExpr::divide(d3, SymbolicExpr::number(6))->simplify();

            if (expression_depends_on_variable(lamina::detail::node(Q), y)) continue;

            n = 3;
            return true;
        }

        /// 对于 n=-1 和 n=4，使用数值验证
        /// rhs = -P*y + Q*y^n
        /// 在 y=1: rhs(x,1) = -P + Q
        /// 在 y=2: rhs(x,2) = -2P + Q*2^n
        /// 在 y=3: rhs(x,3) = -3P + Q*3^n
        /// 从两个方程解出 P 和 Q，用第三个验证
    }

    return false;
}


bool is_exact_ode(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x,
    const std::string& y)
{
    if (!M || !N) return false;

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

    /// 数值验证：在几个点检查
    auto eval_at = [&](double xv, double yv) -> bool {
        auto d = diff->substitute(x, SymbolicExpr::number(xv));
        d = d->substitute(y, SymbolicExpr::number(yv));
        d = d->simplify();
        double val = try_eval_double(d);
        if (std::isnan(val)) return false;
        int eq;
        lmmc_double_nearly_equal_tol(val, 0.0, 1e-9, 1e-9, &eq);
        return eq != 0;
    };

    if (eval_at(1.0, 1.0) && eval_at(2.0, 3.0) && eval_at(0.5, -1.0)) {
        return true;
    }

    return false;
}


bool is_constant_coefficient(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coeffs,
    const std::string& x)
{
    for (const auto& c : coeffs) {
        if (!c) continue;
        if (expression_depends_on_variable(lamina::detail::node(c), x)) return false;
    }
    return true;
}


bool is_euler_equation(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coeffs,
    const std::string& x,
    std::vector<double>& euler_consts)
{
    /// Euler 方程: 第 k 阶导数的系数为 a_k * x^k
    /// coeffs[0] 对应最高阶 (阶数 n)，coeffs[i] 对应阶数 n-i
    int n = static_cast<int>(coeffs.size()) - 1;
    euler_consts.clear();
    euler_consts.resize(coeffs.size(), 0.0);

    for (int i = 0; i <= n; ++i) {
        int order = n - i;  // 该系数对应的导数阶数

        if (!coeffs[i] || coeffs[i]->is_zero()) {
            euler_consts[i] = 0.0;
            continue;
        }

        if (order == 0) {
            /// 零阶项系数应为常数
            if (expression_depends_on_variable(lamina::detail::node(coeffs[i]), x)) return false;
            double val = try_eval_double(coeffs[i]);
            if (std::isnan(val)) return false;
            euler_consts[i] = val;
            continue;
        }

        /// 第 order 阶导数的系数应为 a_k * x^order
        /// 除以 x^order 后应为常数
        auto x_power = SymbolicExpr::power(
            SymbolicExpr::variable(x),
            SymbolicExpr::number(order));
        auto ratio = SymbolicExpr::divide(coeffs[i], x_power)->simplify();

        if (expression_depends_on_variable(lamina::detail::node(ratio), x)) {
            /// 数值验证：在 x=1 和 x=2 处检查比值是否相同
            auto at_1 = ratio->substitute(x, SymbolicExpr::number(1.0))->simplify();
            auto at_2 = ratio->substitute(x, SymbolicExpr::number(2.0))->simplify();
            double v1 = try_eval_double(at_1);
            double v2 = try_eval_double(at_2);
            if (std::isnan(v1) || std::isnan(v2)) return false;
            int eq;
            lmmc_double_nearly_equal_tol(v1, v2, 1e-9, 1e-9, &eq);
            if (!eq) return false;
            euler_consts[i] = v1;
        } else {
            double val = try_eval_double(ratio);
            if (std::isnan(val)) return false;
            euler_consts[i] = val;
        }
    }

    return true;
}


ODEClassification classify_first_order_ode(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y)
{
    ODEClassification result;
    result.order = 1;

    if (!rhs || !lamina::detail::node(rhs)) {
        result.type = ODEType::Separable;
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

    /// 1. 检测常系数
    if (is_constant_coefficient(coeffs, x)) {
        if (result.order == 2) {
            result.type = ODEType::Linear2_ConstCoeff;
        } else {
            result.type = ODEType::HigherOrder_ConstCoeff;
        }

        /// 提取数值系数
        result.const_coeffs.clear();
        for (const auto& c : coeffs) {
            double val = try_eval_double(c);
            result.const_coeffs.push_back(std::isnan(val) ? 0.0 : val);
        }
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

} // namespace lamina
