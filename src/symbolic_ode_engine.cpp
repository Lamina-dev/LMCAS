/**
 * @file symbolic_ode_engine.cpp
 * @brief 统一 ODE 求解引擎实现：类型检测与分类。
 */
#include "../include/symbolic_ode_engine.hpp"
#include "symbolic_ast.hpp"
#include "../include/symbolic.hpp"
#include "../include/poly_utils.hpp"
#include "poly_utils_internal.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include <cmath>
#include <memory>
#include <string>

namespace lamina {

// ============================================================================
/// 辅助函数
// ============================================================================

/**
 * @internal
 * @brief 尝试将表达式求值为 double。
 * @return 若表达式为纯数值返回其值，否则返回 NaN。
 */
static double try_eval_double(const std::shared_ptr<SymbolicExpr>& expr) {
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

static Result<void> validate_ode_expr_var_pair(
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

static Result<void> validate_ode_pair_var_pair(
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

static ODESolutionResult wrap_ode_solution(ODESolution solution,
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

static Result<void> validate_ode_variables(
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

static Result<void> validate_numeric_ode_coefficients(
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

static Result<void> validate_ode_three_expr_one_var(
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

static Result<void> validate_ode_two_expr_point(
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
    return !depends_on_var(lamina::detail::node(expr), x) && !depends_on_var(lamina::detail::node(expr), y);
}

/**
 * @internal
 * @brief 检查表达式是否仅依赖指定变量（不依赖另一个变量）。
 */
[[maybe_unused]] static bool depends_only_on(const std::shared_ptr<SymbolicExpr>& expr,
                                             const std::string&,
                                             const std::string& other_var) {
    if (!expr || !lamina::detail::node(expr)) return true;
    return !depends_on_var(lamina::detail::node(expr), other_var);
}

// ============================================================================
/// is_separable
// ============================================================================

bool is_separable(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y)
{
    if (!rhs || !lamina::detail::node(rhs)) return true;

    bool has_x = depends_on_var(lamina::detail::node(rhs), x);
    bool has_y = depends_on_var(lamina::detail::node(rhs), y);

    /// 若只依赖一个变量或都不依赖，则可分离
    if (!has_x || !has_y) return true;

    /// 检查乘法形式 f(x)*g(y)
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(rhs));
    if (mul) {
        /// 将因子分为仅含 x 的和仅含 y 的
        bool all_separable = true;
        for (const auto& factor : mul->operands()) {
            bool fx = depends_on_var(factor, x);
            bool fy = depends_on_var(factor, y);
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
                bool base_x = depends_on_var(pow->base(), x);
                bool base_y = depends_on_var(pow->base(), y);
                if ((base_x && !base_y) || (!base_x && base_y)) return true;
            }
        }
    }

    return false;
}

// ============================================================================
/// is_linear_first_order
// ============================================================================

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
    if (!depends_on_var(lamina::detail::node(rhs), y)) {
        P = SymbolicExpr::number(0);
        Q = rhs;
        return true;
    }

    /// 对 rhs 关于 y 求导，若结果不依赖 y，则 rhs 关于 y 是线性的
    auto drhs_dy = rhs->differentiate(y);
    if (!drhs_dy || depends_on_var(lamina::detail::node(drhs_dy), y)) {
        return false;
    }

    /// rhs = A(x) + B(x)*y，其中 B(x) = ∂rhs/∂y
    /// 计算 A(x) = rhs|_{y=0}
    auto A = rhs->substitute(y, SymbolicExpr::number(0));
    if (!A) return false;

    /// 验证 A 不依赖 y
    if (depends_on_var(lamina::detail::node(A), y)) return false;

    /// B(x) = drhs_dy（已验证不依赖 y）
    /// 线性形式: y' = A(x) + B(x)*y  →  y' - B(x)*y = A(x)  →  y' + (-B(x))*y = A(x)
    /// 所以 P = -B(x), Q = A(x)
    P = SymbolicExpr::multiply(SymbolicExpr::number(-1), drhs_dy)->simplify();
    Q = A->simplify();
    return true;
}

// ============================================================================
/// is_homogeneous_ode
// ============================================================================

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

// ============================================================================
/// is_bernoulli_ode
// ============================================================================

bool is_bernoulli_ode(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string&,
    const std::string& y,
    std::shared_ptr<SymbolicExpr>& P,
    std::shared_ptr<SymbolicExpr>& Q,
    int& n)
{
    if (!rhs || !lamina::detail::node(rhs)) return false;
    if (!depends_on_var(lamina::detail::node(rhs), y)) return false;

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
        if (depends_on_var(lamina::detail::node(h_at_0), y)) continue;

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
            if (depends_on_var(lamina::detail::node(d2_simplified), y)) continue;

            /// 验证三阶导为零
            auto d3 = d2->differentiate(y)->simplify();
            if (!d3->is_zero()) continue;

            /// P = -(d1|_{y=0}), Q = d2/2
            P = SymbolicExpr::multiply(SymbolicExpr::number(-1), h_at_0)->simplify();
            Q = SymbolicExpr::divide(d2_simplified, SymbolicExpr::number(2))->simplify();

            /// 验证 Q 不依赖 y
            if (depends_on_var(lamina::detail::node(Q), y)) continue;

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
            if (depends_on_var(lamina::detail::node(d3), y)) continue;

            auto d4 = d3->differentiate(y)->simplify();
            if (!d4->is_zero()) continue;

            P = SymbolicExpr::multiply(SymbolicExpr::number(-1), h_at_0)->simplify();
            Q = SymbolicExpr::divide(d3, SymbolicExpr::number(6))->simplify();

            if (depends_on_var(lamina::detail::node(Q), y)) continue;

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

// ============================================================================
/// is_exact_ode
// ============================================================================

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

// ============================================================================
/// is_constant_coefficient
// ============================================================================

bool is_constant_coefficient(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coeffs,
    const std::string& x)
{
    for (const auto& c : coeffs) {
        if (!c) continue;
        if (depends_on_var(lamina::detail::node(c), x)) return false;
    }
    return true;
}

// ============================================================================
/// is_euler_equation
// ============================================================================

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
            if (depends_on_var(lamina::detail::node(coeffs[i]), x)) return false;
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

        if (depends_on_var(lamina::detail::node(ratio), x)) {
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

// ============================================================================
/// classify_first_order_ode
// ============================================================================

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

// ============================================================================
/// classify_higher_order_ode
// ============================================================================

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

// ============================================================================
/// 一阶 ODE 求解方法实现
// ============================================================================

namespace lamina {

// ============================================================================
/// solve_homogeneous_ode
// ============================================================================

ODESolutionResult solve_homogeneous_ode_checked(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y,
    ComputationContext& context)
{
    const std::string operation = "solve_homogeneous_ode";
    auto valid = validate_ode_expr_var_pair(rhs, x, y, context, operation);
    if (!valid) return ODESolutionResult::failure(valid.error());
    auto step = context.consume_steps(12, operation);
    if (!step) return ODESolutionResult::failure(step.error());

    try {
        return wrap_ode_solution(solve_homogeneous_ode(rhs, x, y),
                                 ODEType::Homogeneous,
                                 operation);
    } catch (const std::bad_alloc&) {
        return ODESolutionResult::failure(CasErrc::ResourceLimit,
                                          "homogeneous ODE allocation failed",
                                          operation);
    } catch (const std::exception& e) {
        return ODESolutionResult::failure(CasErrc::InternalInvariant,
                                          e.what(),
                                          operation);
    }
}

ODESolutionResult solve_homogeneous_ode_checked(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y)
{
    ComputationContext context;
    return solve_homogeneous_ode_checked(rhs, x, y, context);
}

ODESolution solve_homogeneous_ode(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y)
{
    ODESolution result;
    result.method_used = ODEType::Homogeneous;
    result.constants = {"C"};

    /// 齐次方程 y' = f(y/x)
    /// 令 v = y/x，则 y = v*x，y' = v + x*v'
    /// 代入: v + x*v' = f(v)
    /// 即: x*v' = f(v) - v
    /// 分离变量: dv/(f(v) - v) = dx/x

    /// 将 rhs 中的 y 替换为 v*x，得到 f 关于 v 和 x 的表达式
    /// 由于齐次性，f(x, vx) = f(1, v)（令 t=1/x 缩放）
    /// 所以直接用 x=1, y=v 代入得到 f(v)
    std::string v_name = "v";
    auto v_var = SymbolicExpr::variable(v_name);
    auto x_var = SymbolicExpr::variable(x);

    /// f(v) = rhs(1, v)：将 x=1, y=v 代入
    auto f_v = rhs->substitute(x, SymbolicExpr::number(1));
    f_v = f_v->substitute(y, v_var);
    f_v = f_v->simplify();

    /// 分离变量方程: dv/(f(v) - v) = dx/x
    /// 积分: ∫ dv/(f(v) - v) = ∫ dx/x = ln|x| + C
    auto f_minus_v = SymbolicExpr::add(f_v,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), v_var));
    f_minus_v = f_minus_v->simplify();

    /// 计算 ∫ 1/(f(v) - v) dv
    auto integrand = SymbolicExpr::divide(SymbolicExpr::number(1), f_minus_v);
    auto lhs_integral = integrand->integrate(v_name);

    /// 右端: ln(x)
    auto rhs_integral = SymbolicExpr::ln(x_var);

    /// 解为: lhs_integral = ln(x) + C
    /// 即: lhs_integral - ln(x) = C
    /// 回代 v = y/x
    auto v_replacement = SymbolicExpr::divide(
        SymbolicExpr::variable(y), x_var);

    auto solution = lhs_integral->substitute(v_name, v_replacement);
    solution = SymbolicExpr::add(solution,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), rhs_integral));
    solution = solution->simplify();

    /// 返回隐式解形式: F(x,y) = C
    result.general_solution = solution;
    return result;
}

// ============================================================================
/// solve_bernoulli_ode
// ============================================================================

ODESolutionResult solve_bernoulli_ode_checked(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    int n,
    const std::string& x,
    const std::string& y,
    ComputationContext& context)
{
    const std::string operation = "solve_bernoulli_ode";
    auto valid = validate_ode_pair_var_pair(P, Q, x, y, context, operation);
    if (!valid) return ODESolutionResult::failure(valid.error());
    if (n == 0 || n == 1) {
        return ODESolutionResult::failure(CasErrc::InvalidArgument,
                                          "Bernoulli exponent must not be 0 or 1",
                                          operation);
    }
    auto step = context.consume_steps(14, operation);
    if (!step) return ODESolutionResult::failure(step.error());

    try {
        return wrap_ode_solution(solve_bernoulli_ode(P, Q, n, x, y),
                                 ODEType::Bernoulli,
                                 operation);
    } catch (const std::bad_alloc&) {
        return ODESolutionResult::failure(CasErrc::ResourceLimit,
                                          "Bernoulli ODE allocation failed",
                                          operation);
    } catch (const std::exception& e) {
        return ODESolutionResult::failure(CasErrc::InternalInvariant,
                                          e.what(),
                                          operation);
    }
}

ODESolutionResult solve_bernoulli_ode_checked(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    int n,
    const std::string& x,
    const std::string& y)
{
    ComputationContext context;
    return solve_bernoulli_ode_checked(P, Q, n, x, y, context);
}

ODESolution solve_bernoulli_ode(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    int n,
    const std::string& x,
    const std::string&)
{
    ODESolution result;
    result.method_used = ODEType::Bernoulli;
    result.constants = {"C"};

    /// Bernoulli 方程: y' + P(x)*y = Q(x)*y^n  (n ≠ 0, 1)
    /// 令 v = y^(1-n)，则 v' = (1-n)*y^(-n)*y'
    /// 从原方程: y' = -P(x)*y + Q(x)*y^n
    /// 两边乘以 (1-n)*y^(-n):
    ///   (1-n)*y^(-n)*y' = (1-n)*(-P(x)*y^(1-n) + Q(x))
    ///   v' = (1-n)*(-P(x)*v + Q(x))
    ///   v' + (1-n)*P(x)*v = (1-n)*Q(x)
    /// 这是关于 v 的一阶线性 ODE

    int one_minus_n = 1 - n;
    auto coeff = SymbolicExpr::number(one_minus_n);

    /// 线性 ODE 的系数: P_linear = (1-n)*P(x), Q_linear = (1-n)*Q(x)
    auto P_linear = SymbolicExpr::multiply(coeff, P)->simplify();
    auto Q_linear = SymbolicExpr::multiply(coeff, Q)->simplify();

    /// 用积分因子法求解线性 ODE: v' + P_linear*v = Q_linear
    /// 积分因子 μ = exp(∫P_linear dx)
    auto intP = P_linear->integrate(x);
    auto mu = SymbolicExpr::exp(intP);

    /// v = (1/μ) * (∫ Q_linear * μ dx + C)
    auto Q_mu = SymbolicExpr::multiply(Q_linear, mu)->simplify();
    auto int_Q_mu = Q_mu->integrate(x);

    auto C_const = SymbolicExpr::variable("C");
    auto numerator = SymbolicExpr::add(int_Q_mu, C_const);
    auto v_solution = SymbolicExpr::divide(numerator, mu)->simplify();

    /// 回代: v = y^(1-n)，所以 y = v^(1/(1-n))
    /// y^(1-n) = v_solution
    /// y = v_solution^(1/(1-n))
    auto exponent = SymbolicExpr::divide(
        SymbolicExpr::number(1),
        SymbolicExpr::number(one_minus_n));
    auto y_solution = SymbolicExpr::power(v_solution, exponent)->simplify();

    result.general_solution = y_solution;
    return result;
}

// ============================================================================
/// find_integrating_factor
// ============================================================================

std::shared_ptr<SymbolicExpr> find_integrating_factor(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x,
    const std::string& y)
{
    if (!M || !N) return nullptr;

    /// 计算 ∂M/∂y - ∂N/∂x
    auto dM_dy = M->differentiate(y);
    auto dN_dx = N->differentiate(x);
    if (!dM_dy || !dN_dx) return nullptr;

    auto diff = SymbolicExpr::add(dM_dy,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dN_dx));
    diff = diff->simplify();

    if (diff->is_zero()) {
        /// 已经恰当，积分因子为 1
        return SymbolicExpr::number(1);
    }

    /// 尝试 μ = μ(x): (∂M/∂y - ∂N/∂x) / N 仅依赖 x
    if (N && !N->is_zero()) {
        auto ratio_x = SymbolicExpr::divide(diff, N)->simplify();
        if (!depends_on_var(lamina::detail::node(ratio_x), y)) {
            /// μ(x) = exp(∫ ratio_x dx)
            auto int_ratio = ratio_x->integrate(x);
            return SymbolicExpr::exp(int_ratio);
        }
    }

    /// 尝试 μ = μ(y): (∂N/∂x - ∂M/∂y) / M 仅依赖 y
    if (M && !M->is_zero()) {
        auto neg_diff = SymbolicExpr::multiply(SymbolicExpr::number(-1), diff)->simplify();
        auto ratio_y = SymbolicExpr::divide(neg_diff, M)->simplify();
        if (!depends_on_var(lamina::detail::node(ratio_y), x)) {
            /// μ(y) = exp(∫ ratio_y dy)
            auto int_ratio = ratio_y->integrate(y);
            return SymbolicExpr::exp(int_ratio);
        }
    }

    return nullptr;
}

// ============================================================================
/// solve_exact_ode
// ============================================================================

ODESolutionResult solve_exact_ode_checked(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x,
    const std::string& y,
    ComputationContext& context)
{
    const std::string operation = "solve_exact_ode";
    auto valid = validate_ode_pair_var_pair(M, N, x, y, context, operation);
    if (!valid) return ODESolutionResult::failure(valid.error());
    auto step = context.consume_steps(18, operation);
    if (!step) return ODESolutionResult::failure(step.error());

    try {
        return wrap_ode_solution(solve_exact_ode(M, N, x, y),
                                 ODEType::Exact,
                                 operation);
    } catch (const std::bad_alloc&) {
        return ODESolutionResult::failure(CasErrc::ResourceLimit,
                                          "exact ODE allocation failed",
                                          operation);
    } catch (const std::exception& e) {
        return ODESolutionResult::failure(CasErrc::InternalInvariant,
                                          e.what(),
                                          operation);
    }
}

ODESolutionResult solve_exact_ode_checked(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x,
    const std::string& y)
{
    ComputationContext context;
    return solve_exact_ode_checked(M, N, x, y, context);
}

ODESolution solve_exact_ode(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x,
    const std::string& y)
{
    ODESolution result;
    result.method_used = ODEType::Exact;
    result.constants = {"C"};

    auto M_eff = M;
    auto N_eff = N;

    /// 检查是否恰当，若不恰当则尝试寻找积分因子
    if (!is_exact_ode(M, N, x, y)) {
        auto mu = find_integrating_factor(M, N, x, y);
        if (!mu) {
            /// 无法找到积分因子
            result.general_solution = nullptr;
            return result;
        }
        /// 乘以积分因子
        M_eff = SymbolicExpr::multiply(mu, M)->simplify();
        N_eff = SymbolicExpr::multiply(mu, N)->simplify();
    }

    /// 恰当方程: ∂F/∂x = M, ∂F/∂y = N
    /// F(x,y) = ∫ M dx + g(y)
    /// 对 F 关于 y 求导: ∂F/∂y = ∂/∂y(∫ M dx) + g'(y) = N
    /// 所以 g'(y) = N - ∂/∂y(∫ M dx)

    /// 步骤 1: 计算 ∫ M dx（将 y 视为常数）
    auto F_partial = M_eff->integrate(x);

    /// 步骤 2: 对 F_partial 关于 y 求导
    auto dF_dy = F_partial->differentiate(y);
    dF_dy = dF_dy ? dF_dy->simplify() : SymbolicExpr::number(0);

    /// 步骤 3: g'(y) = N - dF_dy
    auto g_prime = SymbolicExpr::add(N_eff,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dF_dy));
    g_prime = g_prime->simplify();

    /// 步骤 4: g(y) = ∫ g'(y) dy
    auto g_y = g_prime->integrate(y);

    /// 步骤 5: F(x,y) = F_partial + g(y)
    auto F_total = SymbolicExpr::add(F_partial, g_y)->simplify();

    /// 解为 F(x,y) = C
    result.general_solution = F_total;
    return result;
}

} // namespace lamina

// ============================================================================
/// 高阶常系数 ODE 求解
// ============================================================================

namespace lamina {

/**
 * @internal
 * @brief 表示特征多项式的一个根及其重数。
 */
struct CharRoot {
    double real_part;   ///< 实部
    double imag_part;   ///< 虚部（为零表示实根）
    int multiplicity;   ///< 重数
    bool is_complex;    ///< 是否为复根
};

/**
 * @internal
 * @brief 使用数值方法求解特征多项式的所有根。
 *
 * 对于度数 ≤ 4 的多项式使用解析公式，
 * 对于度数 5-6 使用 Durand-Kerner 迭代法。
 */
[[maybe_unused]] static std::vector<CharRoot> find_characteristic_roots(
    const std::vector<double>& coeffs)
{
    int n = static_cast<int>(coeffs.size()) - 1;
    if (n <= 0) return {};

    std::vector<CharRoot> roots;

    /// 归一化系数（使最高次项系数为 1）
    double leading = coeffs[0];
    if (std::abs(leading) < 1e-15) return {};

    std::vector<double> norm_coeffs(coeffs.size());
    for (size_t i = 0; i < coeffs.size(); ++i) {
        norm_coeffs[i] = coeffs[i] / leading;
    }

    /// 对于低阶多项式，使用解析公式
    if (n == 1) {
        /// r + norm_coeffs[1] = 0
        double r = -norm_coeffs[1];
        roots.push_back({r, 0.0, 1, false});
        return roots;
    }

    if (n == 2) {
        double b = norm_coeffs[1];
        double c = norm_coeffs[2];
        double D = b * b - 4.0 * c;
        int eq;
        lmmc_double_nearly_equal_tol(D, 0.0, 1e-10, 1e-10, &eq);
        if (eq) {
            roots.push_back({-b / 2.0, 0.0, 2, false});
        } else if (D > 0) {
            double r1 = (-b + std::sqrt(D)) / 2.0;
            double r2 = (-b - std::sqrt(D)) / 2.0;
            roots.push_back({r1, 0.0, 1, false});
            roots.push_back({r2, 0.0, 1, false});
        } else {
            double re = -b / 2.0;
            double im = std::sqrt(-D) / 2.0;
            roots.push_back({re, im, 1, true});
        }
        return roots;
    }

    /// 对于 n >= 3，使用 Durand-Kerner 方法求所有根
    /// 初始化：在单位圆上均匀分布初始猜测
    struct Complex {
        double re, im;
        Complex(double r = 0, double i = 0) : re(r), im(i) {}
        Complex operator*(const Complex& o) const {
            return {re * o.re - im * o.im, re * o.im + im * o.re};
        }
        Complex operator+(const Complex& o) const {
            return {re + o.re, im + o.im};
        }
        Complex operator-(const Complex& o) const {
            return {re - o.re, im - o.im};
        }
        Complex operator/(const Complex& o) const {
            double d = o.re * o.re + o.im * o.im;
            if (d < 1e-30) return {0, 0};
            return {(re * o.re + im * o.im) / d,
                    (im * o.re - re * o.im) / d};
        }
        double mag() const { return std::sqrt(re * re + im * im); }
    };

    std::vector<Complex> z(n);
    /// 初始猜测：使用不同半径的点避免对称性问题
    for (int i = 0; i < n; ++i) {
        double angle = 2.0 * 3.14159265358979323846 * i / n + 0.1;
        double radius = 1.0 + 0.3 * i;
        z[i] = {radius * std::cos(angle), radius * std::sin(angle)};
    }

    /// 求值多项式 p(z)
    auto eval_poly = [&](const Complex& val) -> Complex {
        Complex result = {1.0, 0.0};
        for (int i = 1; i <= n; ++i) {
            result = result * val + Complex{norm_coeffs[i], 0.0};
        }
        return result;
    };

    /// Durand-Kerner 迭代
    for (int iter = 0; iter < 1000; ++iter) {
        double max_change = 0.0;
        for (int i = 0; i < n; ++i) {
            Complex num = eval_poly(z[i]);
            Complex denom = {1.0, 0.0};
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    denom = denom * (z[i] - z[j]);
                }
            }
            Complex delta = num / denom;
            z[i] = z[i] - delta;
            double change = delta.mag();
            if (change > max_change) max_change = change;
        }
        if (max_change < 1e-12) break;
    }

    /// 将数值根分类为实根和复共轭对，并检测重根
    std::vector<bool> used(n, false);
    for (int i = 0; i < n; ++i) {
        if (used[i]) continue;

        /// 检查是否为实根（虚部接近零）
        int is_real;
        lmmc_double_nearly_equal_tol(z[i].im, 0.0, 1e-8, 1e-8, &is_real);

        if (is_real) {
            double r = z[i].re;
            int mult = 1;
            used[i] = true;
            /// 检查重根
            for (int j = i + 1; j < n; ++j) {
                if (used[j]) continue;
                int j_real;
                lmmc_double_nearly_equal_tol(z[j].im, 0.0, 1e-8, 1e-8, &j_real);
                if (!j_real) continue;
                int same;
                lmmc_double_nearly_equal_tol(z[j].re, r, 1e-8, 1e-8, &same);
                if (same) {
                    mult++;
                    used[j] = true;
                }
            }
            roots.push_back({r, 0.0, mult, false});
        } else {
            /// 复根：找共轭对
            double re = z[i].re;
            double im = std::abs(z[i].im);
            used[i] = true;
            int mult = 1;

            /// 找到共轭根并标记
            for (int j = i + 1; j < n; ++j) {
                if (used[j]) continue;
                int same_re, conj_im;
                lmmc_double_nearly_equal_tol(z[j].re, re, 1e-8, 1e-8, &same_re);
                lmmc_double_nearly_equal_tol(z[j].im, -z[i].im, 1e-8, 1e-8, &conj_im);
                if (same_re && conj_im) {
                    used[j] = true;
                    break;
                }
            }
            /// 检查重复的复共轭对
            for (int j = i + 1; j < n; ++j) {
                if (used[j]) continue;
                int same_re, same_im;
                lmmc_double_nearly_equal_tol(z[j].re, re, 1e-8, 1e-8, &same_re);
                lmmc_double_nearly_equal_tol(std::abs(z[j].im), im, 1e-8, 1e-8, &same_im);
                if (same_re && same_im) {
                    mult++;
                    used[j] = true;
                    /// 也标记其共轭
                    for (int k = j + 1; k < n; ++k) {
                        if (used[k]) continue;
                        int k_re, k_im;
                        lmmc_double_nearly_equal_tol(z[k].re, re, 1e-8, 1e-8, &k_re);
                        lmmc_double_nearly_equal_tol(z[k].im, -z[j].im, 1e-8, 1e-8, &k_im);
                        if (k_re && k_im) {
                            used[k] = true;
                            break;
                        }
                    }
                }
            }
            roots.push_back({re, im, mult, true});
        }
    }

    return roots;
}

/**
 * @internal
 * @brief 将 double 值转为"干净"的数值表达式。
 *
 * 若值接近整数或简单分数，使用精确表示。
 */
static std::shared_ptr<SymbolicExpr> clean_number(double val) {
    /// 检查是否接近整数
    double rounded = std::round(val);
    int eq;
    lmmc_double_nearly_equal_tol(val, rounded, 1e-10, 1e-10, &eq);
    if (eq) {
        return SymbolicExpr::number(static_cast<int>(rounded));
    }

    /// 检查是否接近简单分数 p/q (q <= 12)
    for (int q = 2; q <= 12; ++q) {
        double p = val * q;
        double p_rounded = std::round(p);
        lmmc_double_nearly_equal_tol(p, p_rounded, 1e-10, 1e-10, &eq);
        if (eq) {
            return SymbolicExpr::divide(
                SymbolicExpr::number(static_cast<int>(p_rounded)),
                SymbolicExpr::number(q));
        }
    }

    return SymbolicExpr::number(val);
}

static bool has_nonzero_forcing(const std::shared_ptr<SymbolicExpr>& forcing) {
    return forcing && lamina::detail::node(forcing) && !forcing->is_zero();
}

/**
 * @internal
 * @brief 根据特征根构造齐次通解。
 *
 * - 实根 r（重数 m）：C_k * x^k * e^(rx)，k = 0, ..., m-1
 * - 复根 α±βi（重数 m）：x^k * e^(αx) * (C_a*cos(βx) + C_b*sin(βx))
 */
[[maybe_unused]] static std::shared_ptr<SymbolicExpr> build_homogeneous_solution(
    const std::vector<CharRoot>& roots,
    const std::string& x,
    std::vector<std::string>& constants)
{
    auto x_var = SymbolicExpr::variable(x);
    std::shared_ptr<SymbolicExpr> solution = nullptr;
    int const_idx = 1;

    for (const auto& root : roots) {
        for (int k = 0; k < root.multiplicity; ++k) {
            if (!root.is_complex) {
                /// 实根: C_i * x^k * e^(r*x)
                std::string c_name = "C" + std::to_string(const_idx++);
                constants.push_back(c_name);
                auto C = SymbolicExpr::variable(c_name);

                std::shared_ptr<SymbolicExpr> term = C;

                /// x^k 因子
                if (k > 0) {
                    auto x_pow = SymbolicExpr::power(x_var, SymbolicExpr::number(k));
                    term = SymbolicExpr::multiply(term, x_pow);
                }

                /// e^(r*x) 因子
                int r_zero;
                lmmc_double_nearly_equal_tol(root.real_part, 0.0, 1e-10, 1e-10, &r_zero);
                if (!r_zero) {
                    auto r_expr = clean_number(root.real_part);
                    auto exp_arg = SymbolicExpr::multiply(r_expr, x_var);
                    auto exp_term = SymbolicExpr::exp(exp_arg);
                    term = SymbolicExpr::multiply(term, exp_term);
                }

                solution = solution ? SymbolicExpr::add(solution, term) : term;
            } else {
                /// 复根 α±βi: 产生两个基本解
                /// C_a * x^k * e^(αx) * cos(βx)
                /// C_b * x^k * e^(αx) * sin(βx)
                std::string ca_name = "C" + std::to_string(const_idx++);
                std::string cb_name = "C" + std::to_string(const_idx++);
                constants.push_back(ca_name);
                constants.push_back(cb_name);
                auto Ca = SymbolicExpr::variable(ca_name);
                auto Cb = SymbolicExpr::variable(cb_name);

                auto beta_expr = clean_number(root.imag_part);
                auto beta_x = SymbolicExpr::multiply(beta_expr, x_var);
                auto cos_term = SymbolicExpr::cos(beta_x);
                auto sin_term = SymbolicExpr::sin(beta_x);

                /// e^(αx) 因子
                std::shared_ptr<SymbolicExpr> exp_factor = nullptr;
                int alpha_zero;
                lmmc_double_nearly_equal_tol(root.real_part, 0.0, 1e-10, 1e-10, &alpha_zero);
                if (!alpha_zero) {
                    auto alpha_expr = clean_number(root.real_part);
                    auto exp_arg = SymbolicExpr::multiply(alpha_expr, x_var);
                    exp_factor = SymbolicExpr::exp(exp_arg);
                }

                /// x^k 因子
                std::shared_ptr<SymbolicExpr> x_pow_factor = nullptr;
                if (k > 0) {
                    x_pow_factor = SymbolicExpr::power(x_var, SymbolicExpr::number(k));
                }

                /// 构造 cos 项: Ca * x^k * e^(αx) * cos(βx)
                auto term_cos = Ca;
                if (x_pow_factor) term_cos = SymbolicExpr::multiply(term_cos, x_pow_factor);
                if (exp_factor) term_cos = SymbolicExpr::multiply(term_cos, exp_factor);
                term_cos = SymbolicExpr::multiply(term_cos, cos_term);

                /// 构造 sin 项: Cb * x^k * e^(αx) * sin(βx)
                auto term_sin = Cb;
                if (x_pow_factor) term_sin = SymbolicExpr::multiply(term_sin, x_pow_factor);
                if (exp_factor) term_sin = SymbolicExpr::multiply(term_sin, exp_factor);
                term_sin = SymbolicExpr::multiply(term_sin, sin_term);

                solution = solution ? SymbolicExpr::add(solution, term_cos) : term_cos;
                solution = SymbolicExpr::add(solution, term_sin);
            }
        }
    }

    return solution;
}

static std::shared_ptr<SymbolicExpr> build_euler_solution(
    const std::vector<CharRoot>& roots,
    const std::string& x,
    std::vector<std::string>& constants)
{
    auto x_var = SymbolicExpr::variable(x);
    auto ln_x = SymbolicExpr::ln(x_var);
    std::shared_ptr<SymbolicExpr> solution = nullptr;
    int const_idx = 1;

    for (const auto& root : roots) {
        for (int k = 0; k < root.multiplicity; ++k) {
            std::shared_ptr<SymbolicExpr> log_factor = nullptr;
            if (k > 0) {
                log_factor = SymbolicExpr::power(ln_x, SymbolicExpr::number(k));
            }

            if (!root.is_complex) {
                std::string c_name = "C" + std::to_string(const_idx++);
                constants.push_back(c_name);
                auto term = SymbolicExpr::variable(c_name);

                auto exponent = clean_number(root.real_part);
                term = SymbolicExpr::multiply(
                    term,
                    SymbolicExpr::power(x_var, exponent));
                if (log_factor) {
                    term = SymbolicExpr::multiply(term, log_factor);
                }
                solution = solution ? SymbolicExpr::add(solution, term) : term;
            } else {
                std::string ca_name = "C" + std::to_string(const_idx++);
                std::string cb_name = "C" + std::to_string(const_idx++);
                constants.push_back(ca_name);
                constants.push_back(cb_name);

                auto alpha_expr = clean_number(root.real_part);
                auto beta_expr = clean_number(root.imag_part);
                auto beta_ln_x = SymbolicExpr::multiply(beta_expr, ln_x);
                auto trig = SymbolicExpr::add(
                    SymbolicExpr::multiply(SymbolicExpr::variable(ca_name),
                                           SymbolicExpr::cos(beta_ln_x)),
                    SymbolicExpr::multiply(SymbolicExpr::variable(cb_name),
                                           SymbolicExpr::sin(beta_ln_x)));
                auto term = SymbolicExpr::multiply(
                    SymbolicExpr::power(x_var, alpha_expr),
                    trig);
                if (log_factor) {
                    term = SymbolicExpr::multiply(term, log_factor);
                }
                solution = solution ? SymbolicExpr::add(solution, term) : term;
            }
        }
    }

    return solution;
}

ODESolutionResult solve_higher_order_ode_checked(
    const std::vector<double>& coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string& y,
    ComputationContext& context)
{
    const std::string operation = "solve_higher_order_ode";
    auto variables = validate_ode_variables(x, y, context, operation);
    if (!variables) return ODESolutionResult::failure(variables.error());

    auto coeff_check = validate_numeric_ode_coefficients(coeffs, 2, 7, operation);
    if (!coeff_check) return ODESolutionResult::failure(coeff_check.error());

    auto budget = context.consume_steps(coeffs.size() * 20 + 20, operation);
    if (!budget) return ODESolutionResult::failure(budget.error());

    if (has_nonzero_forcing(forcing)) {
        return ODESolutionResult::failure(
            CasErrc::Inconclusive,
            "checked higher-order ODE currently supports homogeneous constant-coefficient equations only",
            operation);
    }

    try {
        return wrap_ode_solution(
            solve_higher_order_ode(coeffs, forcing, x, y),
            ODEType::HigherOrder_ConstCoeff,
            operation);
    } catch (const std::bad_alloc&) {
        return ODESolutionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while solving higher-order ODE",
            operation);
    } catch (const std::exception& ex) {
        return ODESolutionResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            operation);
    }
}

ODESolutionResult solve_higher_order_ode_checked(
    const std::vector<double>& coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string& y)
{
    ComputationContext context;
    return solve_higher_order_ode_checked(coeffs, forcing, x, y, context);
}

ODESolution solve_higher_order_ode(
    const std::vector<double>& coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string&)
{
    ODESolution result;
    result.method_used = ODEType::HigherOrder_ConstCoeff;

    if (coeffs.size() < 2 || has_nonzero_forcing(forcing)) {
        result.general_solution = nullptr;
        return result;
    }

    auto roots = find_characteristic_roots(coeffs);
    if (roots.empty()) {
        result.general_solution = nullptr;
        return result;
    }

    result.general_solution = build_homogeneous_solution(roots, x, result.constants);
    if (result.general_solution) {
        result.general_solution = result.general_solution->simplify();
    }
    return result;
}

ODESolutionResult solve_euler_ode_checked(
    const std::vector<double>& euler_coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string& y,
    ComputationContext& context)
{
    const std::string operation = "solve_euler_ode";
    auto variables = validate_ode_variables(x, y, context, operation);
    if (!variables) return ODESolutionResult::failure(variables.error());

    auto coeff_check = validate_numeric_ode_coefficients(euler_coeffs, 3, 4, operation);
    if (!coeff_check) return ODESolutionResult::failure(coeff_check.error());

    auto budget = context.consume_steps(euler_coeffs.size() * 20 + 20, operation);
    if (!budget) return ODESolutionResult::failure(budget.error());

    if (has_nonzero_forcing(forcing)) {
        return ODESolutionResult::failure(
            CasErrc::Inconclusive,
            "checked Euler ODE currently supports homogeneous equations only",
            operation);
    }

    try {
        return wrap_ode_solution(
            solve_euler_ode(euler_coeffs, forcing, x, y),
            ODEType::Euler,
            operation);
    } catch (const std::bad_alloc&) {
        return ODESolutionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while solving Euler ODE",
            operation);
    } catch (const std::exception& ex) {
        return ODESolutionResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            operation);
    }
}

ODESolutionResult solve_euler_ode_checked(
    const std::vector<double>& euler_coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string& y)
{
    ComputationContext context;
    return solve_euler_ode_checked(euler_coeffs, forcing, x, y, context);
}

ODESolution solve_euler_ode(
    const std::vector<double>& euler_coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string&)
{
    ODESolution result;
    result.method_used = ODEType::Euler;

    if ((euler_coeffs.size() != 3 && euler_coeffs.size() != 4) ||
        has_nonzero_forcing(forcing)) {
        result.general_solution = nullptr;
        return result;
    }

    std::vector<double> characteristic;
    if (euler_coeffs.size() == 3) {
        const double a = euler_coeffs[0];
        const double b = euler_coeffs[1];
        const double c = euler_coeffs[2];
        characteristic = {a, b - a, c};
    } else {
        const double a = euler_coeffs[0];
        const double b = euler_coeffs[1];
        const double c = euler_coeffs[2];
        const double d = euler_coeffs[3];
        characteristic = {a, b - 3.0 * a, 2.0 * a - b + c, d};
    }

    auto roots = find_characteristic_roots(characteristic);
    if (roots.empty()) {
        result.general_solution = nullptr;
        return result;
    }

    result.general_solution = build_euler_solution(roots, x, result.constants);
    if (result.general_solution) {
        result.general_solution = result.general_solution->simplify();
    }
    return result;
}

/**
 * @internal
 * @brief 检测非齐次项的类型，用于待定系数法。
 *
 * 支持的类型：
 * - 多项式: x^n
 * - 指数: e^(ax)
 * - 三角: cos(bx), sin(bx)
 */

} // namespace lamina

// ============================================================================
/// 参数变分法与 Frobenius 级数解实现
// ============================================================================

namespace lamina {

// ============================================================================
/// solve_variation_of_parameters
// ============================================================================

ODESolutionResult solve_variation_of_parameters_checked(
    const std::shared_ptr<SymbolicExpr>& y1,
    const std::shared_ptr<SymbolicExpr>& y2,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& x,
    ComputationContext& context)
{
    const std::string operation = "solve_variation_of_parameters";
    auto valid = validate_ode_three_expr_one_var(y1, y2, g, x, context, operation);
    if (!valid) return ODESolutionResult::failure(valid.error());

    auto budget = context.consume_steps(24, operation);
    if (!budget) return ODESolutionResult::failure(budget.error());

    try {
        return wrap_ode_solution(
            solve_variation_of_parameters(y1, y2, g, x),
            ODEType::HigherOrder_ConstCoeff,
            operation);
    } catch (const std::bad_alloc&) {
        return ODESolutionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while applying variation of parameters",
            operation);
    } catch (const std::exception& ex) {
        return ODESolutionResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            operation);
    }
}

ODESolutionResult solve_variation_of_parameters_checked(
    const std::shared_ptr<SymbolicExpr>& y1,
    const std::shared_ptr<SymbolicExpr>& y2,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& x)
{
    ComputationContext context;
    return solve_variation_of_parameters_checked(y1, y2, g, x, context);
}

ODESolution solve_variation_of_parameters(
    const std::shared_ptr<SymbolicExpr>& y1,
    const std::shared_ptr<SymbolicExpr>& y2,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& x)
{
    ODESolution result;
    result.method_used = ODEType::HigherOrder_ConstCoeff;
    result.constants = {};

    if (!y1 || !y2 || !g) {
        result.general_solution = nullptr;
        return result;
    }

    /// 计算 y₁' 和 y₂'
    auto y1_prime = y1->differentiate(x);
    auto y2_prime = y2->differentiate(x);

    if (!y1_prime || !y2_prime) {
        result.general_solution = nullptr;
        return result;
    }

    /// 计算 Wronskian: W = y₁·y₂' - y₂·y₁'
    auto term1 = SymbolicExpr::multiply(y1, y2_prime);
    auto term2 = SymbolicExpr::multiply(y2, y1_prime);
    auto W = SymbolicExpr::add(term1,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), term2));
    W = W->simplify();

    if (!W || W->is_zero()) {
        result.general_solution = nullptr;
        return result;
    }

    /// 计算 u₁' = -y₂·g(x)/W
    auto u1_prime = SymbolicExpr::divide(
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(y2, g)),
        W);
    u1_prime = u1_prime->simplify();

    /// 计算 u₂' = y₁·g(x)/W
    auto u2_prime = SymbolicExpr::divide(
        SymbolicExpr::multiply(y1, g),
        W);
    u2_prime = u2_prime->simplify();

    /// 积分得 u₁ 和 u₂
    auto u1 = u1_prime->integrate(x);
    auto u2 = u2_prime->integrate(x);

    if (!u1 || !u2) {
        result.general_solution = nullptr;
        return result;
    }

    /// 特解: y_p = u₁·y₁ + u₂·y₂
    auto y_p = SymbolicExpr::add(
        SymbolicExpr::multiply(u1, y1),
        SymbolicExpr::multiply(u2, y2));
    y_p = y_p->simplify();

    result.general_solution = y_p;
    return result;
}

// ============================================================================
/// classify_singular_point
// ============================================================================

ODESingularityType classify_singular_point(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x)
{
    if (!p || !q) return ODESingularityType::Ordinary;

    /// 检查 p(x) 和 q(x) 在 x₀ 处是否解析（有限）
    auto p_at_x0 = p->substitute(x, x0)->simplify();
    auto q_at_x0 = q->substitute(x, x0)->simplify();

    double p_val = try_eval_double(p_at_x0);
    double q_val = try_eval_double(q_at_x0);

    if (!std::isnan(p_val) && !std::isinf(p_val) &&
        !std::isnan(q_val) && !std::isinf(q_val)) {
        return ODESingularityType::Ordinary;
    }

    /// 检查正则奇点条件
    auto x_var = SymbolicExpr::variable(x);
    auto x_minus_x0 = SymbolicExpr::add(x_var,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x0));

    auto xp = SymbolicExpr::multiply(x_minus_x0, p)->simplify();
    auto x2q = SymbolicExpr::multiply(
        SymbolicExpr::power(x_minus_x0, SymbolicExpr::number(2)), q)->simplify();

    /// 尝试用极限计算 (x-x₀)·p(x) 和 (x-x₀)²·q(x) 在 x₀ 处的值
    auto xp_limit = xp->limit(x, x0);
    auto x2q_limit = x2q->limit(x, x0);

    double xp_val = xp_limit ? try_eval_double(xp_limit) : std::numeric_limits<double>::quiet_NaN();
    double x2q_val = x2q_limit ? try_eval_double(x2q_limit) : std::numeric_limits<double>::quiet_NaN();

    if (!std::isnan(xp_val) && !std::isinf(xp_val) &&
        !std::isnan(x2q_val) && !std::isinf(x2q_val)) {
        return ODESingularityType::RegularSingular;
    }

    return ODESingularityType::IrregularSingular;
}

// ============================================================================
/// solve_frobenius
// ============================================================================

static Result<void> validate_frobenius_regular_singular_domain(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x,
    const std::string& operation)
{
    auto x_var = SymbolicExpr::variable(x);
    auto x_minus_x0 = SymbolicExpr::add(x_var,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x0));
    auto xp_expr = SymbolicExpr::multiply(x_minus_x0, p)->simplify();
    auto x2q_expr = SymbolicExpr::multiply(
        SymbolicExpr::power(x_minus_x0, SymbolicExpr::number(2)), q)->simplify();

    auto P0_expr = xp_expr->limit(x, x0);
    auto Q0_expr = x2q_expr->limit(x, x0);
    double P0 = P0_expr ? try_eval_double(P0_expr) : std::numeric_limits<double>::quiet_NaN();
    double Q0 = Q0_expr ? try_eval_double(Q0_expr) : std::numeric_limits<double>::quiet_NaN();
    if (!std::isfinite(P0) || !std::isfinite(Q0)) {
        return Result<void>::failure(
            CasErrc::Inconclusive,
            "Frobenius regular-singular coefficients are outside the checked numeric support domain",
            operation);
    }

    double discriminant = (P0 - 1.0) * (P0 - 1.0) - 4.0 * Q0;
    int zero_discriminant;
    lmmc_double_nearly_equal_tol(discriminant, 0.0, 1e-12, 1e-12,
                                 &zero_discriminant);
    if (discriminant < 0.0 && !zero_discriminant) {
        return Result<void>::failure(
            CasErrc::Inconclusive,
            "Frobenius checked API currently supports real indicial roots only",
            operation);
    }

    return Result<void>::success();
}

FrobeniusSolutionResult solve_frobenius_checked(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x,
    int order,
    ComputationContext& context)
{
    const std::string operation = "solve_frobenius";
    auto valid = validate_ode_two_expr_point(p, q, x0, x, context, operation);
    if (!valid) return FrobeniusSolutionResult::failure(valid.error());
    if (order < 0 || order > 64) {
        return FrobeniusSolutionResult::failure(
            CasErrc::InvalidArgument,
            "Frobenius truncation order must be between 0 and 64",
            operation);
    }
    auto budget = context.consume_steps(static_cast<std::size_t>(order + 1) * 24 + 24,
                                        operation);
    if (!budget) return FrobeniusSolutionResult::failure(budget.error());

    try {
        auto point_type = classify_singular_point(p, q, x0, x);
        if (point_type == ODESingularityType::IrregularSingular) {
            return FrobeniusSolutionResult::failure(
                CasErrc::Inconclusive,
                "Frobenius checked API does not support irregular singular points",
                operation);
        }
        if (point_type == ODESingularityType::RegularSingular) {
            auto regular_domain =
                validate_frobenius_regular_singular_domain(p, q, x0, x, operation);
            if (!regular_domain) {
                return FrobeniusSolutionResult::failure(regular_domain.error());
            }
        }

        auto solution = solve_frobenius(p, q, x0, x, order);
        if (!solution.series_solution || !lamina::detail::node(solution.series_solution)) {
            return FrobeniusSolutionResult::failure(
                CasErrc::Inconclusive,
                "Frobenius solver produced no series in the supported domain",
                operation);
        }
        if (solution.point_type != point_type) {
            return FrobeniusSolutionResult::failure(
                CasErrc::InternalInvariant,
                "Frobenius solver reported an unexpected singularity type",
                operation);
        }
        return FrobeniusSolutionResult::success(std::move(solution));
    } catch (const std::bad_alloc&) {
        return FrobeniusSolutionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while solving Frobenius series",
            operation);
    } catch (const std::exception& ex) {
        return FrobeniusSolutionResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            operation);
    }
}

FrobeniusSolutionResult solve_frobenius_checked(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x,
    int order)
{
    ComputationContext context;
    return solve_frobenius_checked(p, q, x0, x, order, context);
}

FrobeniusSolution solve_frobenius(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x,
    int order)
{
    FrobeniusSolution result;
    result.truncation_order = order;
    result.point_type = classify_singular_point(p, q, x0, x);

    auto x_var = SymbolicExpr::variable(x);
    auto x_minus_x0 = SymbolicExpr::add(x_var,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x0));

    if (result.point_type == ODESingularityType::IrregularSingular) {
        result.series_solution = nullptr;
        return result;
    }

    if (result.point_type == ODESingularityType::Ordinary) {
        /// 常点：y = ∑aₙ(x-x₀)ⁿ
        /// 通过在 x₀ 处求导来提取 p, q 的 Taylor 系数
        std::vector<double> p_coeffs(order + 2, 0.0);
        std::vector<double> q_coeffs(order + 2, 0.0);

        auto p_current = p;
        auto q_current = q;
        double factorial = 1.0;
        for (int k = 0; k <= order + 1; ++k) {
            if (k > 0) factorial *= k;
            auto p_val_expr = p_current->substitute(x, x0)->simplify();
            auto q_val_expr = q_current->substitute(x, x0)->simplify();
            double pv = try_eval_double(p_val_expr);
            double qv = try_eval_double(q_val_expr);
            p_coeffs[k] = std::isnan(pv) ? 0.0 : pv / factorial;
            q_coeffs[k] = std::isnan(qv) ? 0.0 : qv / factorial;
            p_current = p_current->differentiate(x);
            q_current = q_current->differentiate(x);
            if (!p_current) p_current = SymbolicExpr::number(0);
            if (!q_current) q_current = SymbolicExpr::number(0);
        }

        /// 递推确定系数: a₀=1, a₁=0
        std::vector<double> a(order + 1, 0.0);
        a[0] = 1.0;

        /// a_{n+2} = -1/((n+2)(n+1)) * ∑_{k=0}^{n}[(k+1)·a_{k+1}·p_{n-k} + a_k·q_{n-k}]
        for (int n_idx = 0; n_idx + 2 <= order; ++n_idx) {
            double sum = 0.0;
            for (int k = 0; k <= n_idx; ++k) {
                sum += (k + 1) * a[k + 1] * p_coeffs[n_idx - k];
                sum += a[k] * q_coeffs[n_idx - k];
            }
            double denom = static_cast<double>((n_idx + 2) * (n_idx + 1));
            a[n_idx + 2] = -sum / denom;
        }

        /// 构造级数解
        auto series_sol = SymbolicExpr::number(0);
        for (int k = 0; k <= order; ++k) {
            if (std::abs(a[k]) < 1e-15) continue;
            auto coeff = SymbolicExpr::number(a[k]);
            auto power_term = (k == 0) ? SymbolicExpr::number(1)
                : SymbolicExpr::power(x_minus_x0, SymbolicExpr::number(k));
            series_sol = SymbolicExpr::add(series_sol,
                SymbolicExpr::multiply(coeff, power_term));
        }

        result.series_solution = series_sol->simplify();
        result.indicial_roots = {0.0};
        return result;
    }

    /// 正则奇点：Frobenius 方法
    /// 计算 P₀ = lim (x-x₀)·p(x), Q₀ = lim (x-x₀)²·q(x)
    auto xp_expr = SymbolicExpr::multiply(x_minus_x0, p)->simplify();
    auto x2q_expr = SymbolicExpr::multiply(
        SymbolicExpr::power(x_minus_x0, SymbolicExpr::number(2)), q)->simplify();

    auto P0_expr = xp_expr->limit(x, x0);
    auto Q0_expr = x2q_expr->limit(x, x0);

    double P0 = P0_expr ? try_eval_double(P0_expr) : 0.0;
    double Q0 = Q0_expr ? try_eval_double(Q0_expr) : 0.0;
    if (std::isnan(P0)) P0 = 0.0;
    if (std::isnan(Q0)) Q0 = 0.0;

    /// 指标方程: r(r-1) + P₀·r + Q₀ = 0  →  r² + (P₀-1)·r + Q₀ = 0
    double ind_b = P0 - 1.0;
    double ind_c = Q0;
    double ind_D = ind_b * ind_b - 4.0 * ind_c;

    double r1, r2;
    int eq;
    lmmc_double_nearly_equal_tol(ind_D, 0.0, 1e-12, 1e-12, &eq);
    if (!eq && ind_D > 0) {
        r1 = (-ind_b + std::sqrt(ind_D)) / 2.0;
        r2 = (-ind_b - std::sqrt(ind_D)) / 2.0;
    } else if (eq) {
        r1 = -ind_b / 2.0;
        r2 = r1;
    } else {
        r1 = -ind_b / 2.0;
        r2 = r1;
    }
    if (r1 < r2) std::swap(r1, r2);
    result.indicial_roots = {r1, r2};

    /// 展开 (x-x₀)·p(x) 和 (x-x₀)²·q(x) 的 Taylor 系数
    std::vector<double> pn_coeffs(order + 1, 0.0);
    std::vector<double> qn_coeffs(order + 1, 0.0);

    auto xp_current = xp_expr;
    auto x2q_current = x2q_expr;
    double fact = 1.0;
    for (int k = 0; k <= order; ++k) {
        if (k > 0) fact *= k;
        double pv, qv;
        if (k == 0) {
            pv = P0;
            qv = Q0;
        } else {
            auto pv_expr = xp_current->substitute(x, x0)->simplify();
            auto qv_expr = x2q_current->substitute(x, x0)->simplify();
            pv = try_eval_double(pv_expr);
            qv = try_eval_double(qv_expr);
            if (std::isnan(pv)) {
                auto lim = xp_current->limit(x, x0);
                pv = lim ? try_eval_double(lim) : 0.0;
            }
            if (std::isnan(qv)) {
                auto lim = x2q_current->limit(x, x0);
                qv = lim ? try_eval_double(lim) : 0.0;
            }
        }
        pn_coeffs[k] = std::isnan(pv) ? 0.0 : pv / fact;
        qn_coeffs[k] = std::isnan(qv) ? 0.0 : qv / fact;

        xp_current = xp_current->differentiate(x);
        x2q_current = x2q_current->differentiate(x);
        if (!xp_current) xp_current = SymbolicExpr::number(0);
        if (!x2q_current) x2q_current = SymbolicExpr::number(0);
    }

    /// 递推: a_n = -1/F(r1+n) * ∑_{k=0}^{n-1} [(r1+k)·p_{n-k} + q_{n-k}]·a_k
    /// 其中 F(s) = s(s-1) + P₀·s + Q₀
    auto indicial_poly = [&](double s) -> double {
        return s * (s - 1.0) + P0 * s + Q0;
    };

    std::vector<double> a(order + 1, 0.0);
    a[0] = 1.0;

    for (int n_idx = 1; n_idx <= order; ++n_idx) {
        double F_val = indicial_poly(r1 + n_idx);
        if (std::abs(F_val) < 1e-15) {
            a[n_idx] = 0.0;
            continue;
        }
        double sum = 0.0;
        for (int k = 0; k < n_idx; ++k) {
            int idx = n_idx - k;
            if (idx > order) continue;
            double p_term = (idx < static_cast<int>(pn_coeffs.size())) ? pn_coeffs[idx] : 0.0;
            double q_term = (idx < static_cast<int>(qn_coeffs.size())) ? qn_coeffs[idx] : 0.0;
            sum += ((r1 + k) * p_term + q_term) * a[k];
        }
        a[n_idx] = -sum / F_val;
    }

    /// 构造级数解: y = (x-x₀)^r₁ · ∑aₙ·(x-x₀)ⁿ
    auto power_prefix = SymbolicExpr::power(x_minus_x0, SymbolicExpr::number(r1));
    auto series_part = SymbolicExpr::number(0);

    for (int k = 0; k <= order; ++k) {
        if (std::abs(a[k]) < 1e-15) continue;
        auto coeff = SymbolicExpr::number(a[k]);
        auto power_term = (k == 0) ? SymbolicExpr::number(1)
            : SymbolicExpr::power(x_minus_x0, SymbolicExpr::number(k));
        series_part = SymbolicExpr::add(series_part,
            SymbolicExpr::multiply(coeff, power_term));
    }

    result.series_solution = SymbolicExpr::multiply(power_prefix, series_part)->simplify();
    return result;
}

} // namespace lamina

