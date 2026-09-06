/**
 * @file symbolic_ode_engine.cpp
 * @brief 统一 ODE 求解引擎实现:类型检测与分类.
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

namespace LMCAS {

static ODESolution solve_homogeneous_ode_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&, const std::string&);
static ODESolution solve_bernoulli_ode_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::shared_ptr<SymbolicExpr>&,
    int, const std::string&, const std::string&);
static ODESolution solve_exact_ode_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&, const std::string&);


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
        return wrap_ode_solution(solve_homogeneous_ode_impl(rhs, x, y),
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

static ODESolution solve_homogeneous_ode_impl(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y)
{
    ODESolution result;
    result.method_used = ODEType::Homogeneous;
    result.constants = {"C"};

    /// 齐次方程 y' = f(y/x)
    /// 令 v = y/x,则 y = v*x,y' = v + x*v'
    /// 代入: v + x*v' = f(v)
    /// 即: x*v' = f(v) - v
    /// 分离变量: dv/(f(v) - v) = dx/x

    /// 将 rhs 中的 y 替换为 v*x,得到 f 关于 v 和 x 的表达式
    /// 由于齐次性,f(x, vx) = f(1, v)(令 t=1/x 缩放)
    /// 所以直接用 x=1, y=v 代入得到 f(v)
    std::string v_name = "v";
    auto v_var = SymbolicExpr::variable(v_name);
    auto x_var = SymbolicExpr::variable(x);

    /// f(v) = rhs(1, v):将 x=1, y=v 代入
    auto f_v = rhs->substitute(x, SymbolicExpr::number(1));
    f_v = f_v->substitute(y, v_var);
    f_v = f_v->simplify();

    /// 分离变量方程: dv/(f(v) - v) = dx/x
    /// 积分: integral dv/(f(v) - v) = integral dx/x = ln|x| + C
    auto f_minus_v = SymbolicExpr::add(f_v,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), v_var));
    f_minus_v = f_minus_v->simplify();

    /// 计算 integral 1/(f(v) - v) dv
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
        return wrap_ode_solution(solve_bernoulli_ode_impl(P, Q, n, x, y),
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

static ODESolution solve_bernoulli_ode_impl(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    int n,
    const std::string& x,
    const std::string&)
{
    ODESolution result;
    result.method_used = ODEType::Bernoulli;
    result.constants = {"C"};

    /// Bernoulli 方程: y' + P(x)*y = Q(x)*y^n  (n != 0, 1)
    /// 令 v = y^(1-n),则 v' = (1-n)*y^(-n)*y'
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
    /// 积分因子 mu = exp(integralP_linear dx)
    auto intP = P_linear->integrate(x);
    auto mu = SymbolicExpr::exp(intP);

    /// v = (1/mu) * (integral Q_linear * mu dx + C)
    auto Q_mu = SymbolicExpr::multiply(Q_linear, mu)->simplify();
    auto int_Q_mu = Q_mu->integrate(x);

    auto C_const = SymbolicExpr::variable("C");
    auto numerator = SymbolicExpr::add(int_Q_mu, C_const);
    auto v_solution = SymbolicExpr::divide(numerator, mu)->simplify();

    /// 回代: v = y^(1-n),所以 y = v^(1/(1-n))
    /// y^(1-n) = v_solution
    /// y = v_solution^(1/(1-n))
    auto exponent = SymbolicExpr::divide(
        SymbolicExpr::number(1),
        SymbolicExpr::number(one_minus_n));
    auto y_solution = SymbolicExpr::power(v_solution, exponent)->simplify();

    result.general_solution = y_solution;
    return result;
}


std::shared_ptr<SymbolicExpr> find_integrating_factor(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x,
    const std::string& y)
{
    if (!M || !N) return nullptr;

    /// 计算 partialM/partialy - partialN/partialx
    auto dM_dy = M->differentiate(y);
    auto dN_dx = N->differentiate(x);
    if (!dM_dy || !dN_dx) return nullptr;

    auto diff = SymbolicExpr::add(dM_dy,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dN_dx));
    diff = diff->simplify();

    if (diff->is_zero()) {
        /// 已经恰当,积分因子为 1
        return SymbolicExpr::number(1);
    }

    /// 尝试 mu = mu(x): (partialM/partialy - partialN/partialx) / N 仅依赖 x
    if (N && !N->is_zero()) {
        auto ratio_x = SymbolicExpr::divide(diff, N)->simplify();
        if (!expression_depends_on_variable(LMCAS::detail::node(ratio_x), y)) {
            /// mu(x) = exp(integral ratio_x dx)
            auto int_ratio = ratio_x->integrate(x);
            return SymbolicExpr::exp(int_ratio);
        }
    }

    /// 尝试 mu = mu(y): (partialN/partialx - partialM/partialy) / M 仅依赖 y
    if (M && !M->is_zero()) {
        auto neg_diff = SymbolicExpr::multiply(SymbolicExpr::number(-1), diff)->simplify();
        auto ratio_y = SymbolicExpr::divide(neg_diff, M)->simplify();
        if (!expression_depends_on_variable(LMCAS::detail::node(ratio_y), x)) {
            /// mu(y) = exp(integral ratio_y dy)
            auto int_ratio = ratio_y->integrate(y);
            return SymbolicExpr::exp(int_ratio);
        }
    }

    return nullptr;
}


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
        return wrap_ode_solution(solve_exact_ode_impl(M, N, x, y),
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

static ODESolution solve_exact_ode_impl(
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

    /// 检查是否恰当,若不恰当则尝试寻找积分因子
    if (!is_exact_ode(M, N, x, y)) {
        auto mu = find_integrating_factor(M, N, x, y);
        if (!mu) {
            /// 积分因子搜索未决时以空 general_solution 表示.
            result.general_solution = nullptr;
            return result;
        }
        /// 乘以积分因子
        M_eff = SymbolicExpr::multiply(mu, M)->simplify();
        N_eff = SymbolicExpr::multiply(mu, N)->simplify();
    }

    /// 恰当方程: partialF/partialx = M, partialF/partialy = N
    /// F(x,y) = integral M dx + g(y)
    /// 对 F 关于 y 求导: partialF/partialy = partial/partialy(integral M dx) + g'(y) = N
    /// 所以 g'(y) = N - partial/partialy(integral M dx)

    /// 步骤 1: 计算 integral M dx(将 y 视为常数)
    auto F_partial = M_eff->integrate(x);

    /// 步骤 2: 对 F_partial 关于 y 求导
    auto dF_dy = F_partial->differentiate(y);
    dF_dy = dF_dy ? dF_dy->simplify() : SymbolicExpr::number(0);

    /// 步骤 3: g'(y) = N - dF_dy
    auto g_prime = SymbolicExpr::add(N_eff,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dF_dy));
    g_prime = g_prime->simplify();

    /// 步骤 4: g(y) = integral g'(y) dy
    auto g_y = g_prime->integrate(y);

    /// 步骤 5: F(x,y) = F_partial + g(y)
    auto F_total = SymbolicExpr::add(F_partial, g_y)->simplify();

    /// 解为 F(x,y) = C
    result.general_solution = F_total;
    return result;
}

} // namespace LMCAS
