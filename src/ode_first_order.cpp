#include "internal/expression_transform.hpp"
/**
 * @file symbolic_ode_engine.cpp
 * @brief 统一 ODE 求解引擎实现:类型检测与分类.
 */
#include "symbolic_ode_engine.hpp"
#include "internal/symbolic_ast.hpp"
#include "symbolic.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "internal/ode_support.hpp"
#include "integrator.hpp"
#include "residual_verification.hpp"
#include <cmath>
#include <memory>
#include <string>
#include <stdexcept>

namespace LMCAS {

static ODESolutionResult solve_homogeneous_ode_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&, const std::string&, ComputationContext& context);
static ODESolutionResult solve_bernoulli_ode_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::shared_ptr<SymbolicExpr>&,
    int, const std::string&, const std::string&, ComputationContext& context);
static ODESolutionResult solve_exact_ode_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&, const std::string&, ComputationContext& context);

static std::shared_ptr<SymbolicExpr> integrating_factor_from_antiderivative(
    const std::shared_ptr<SymbolicExpr>& antiderivative,
    const std::string& variable, ComputationContext& context)
{

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    const auto logarithm_of_variable =
        [&variable](const std::shared_ptr<const SymbolicNode>& node) {
            const auto function =
                std::dynamic_pointer_cast<const FunctionNode>(node);
            if (!function ||
                function->type() != FunctionNode::FuncType::Ln ||
                function->arguments().size() != 1) {
                return false;
            }
            const auto argument = std::dynamic_pointer_cast<const VariableNode>(
                function->arguments().front());
            return argument && !argument->is_constant() && argument->name() == variable;
        };

    const auto& root = detail::node(antiderivative);
    if (logarithm_of_variable(root)) {
        return SymbolicExpr::variable(variable);
    }
    const auto product = std::dynamic_pointer_cast<const MultiplyNode>(root);
    if (product && product->operands().size() == 2) {
        for (std::size_t index = 0; index < 2; ++index) {
            if (logarithm_of_variable(product->operands()[index])) {
                auto exponent = detail::make_expression_ptr(
                    product->operands()[1 - index]);
                return detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::power(
                    SymbolicExpr::variable(variable), exponent)), transform_budget));
            }
        }
    }
    return detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::exp(antiderivative)), transform_budget));
}


ODESolutionResult solve_homogeneous_ode_checked(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y,
    ComputationContext& context)
{
    const std::string operation = "solve_homogeneous_ode";
    auto valid = validate_ode_expr_var_pair(rhs, x, y, context, operation);
    if (!valid) {
        return ODESolutionResult::failure(valid.error());
    }
    auto step = context.consume_steps(12, operation);
    if (!step) {
        return ODESolutionResult::failure(step.error());
    }

    try {
        auto solved = solve_homogeneous_ode_impl(rhs, x, y, context);
        if (!solved) { return solved; }
        return wrap_ode_solution(std::move(solved.value()), ODEType::Homogeneous, operation);
    } catch (const CasError& error) {
        return ODESolutionResult::failure(error);
    } catch (const std::bad_alloc&) {
        return ODESolutionResult::failure(CasErrc::ResourceLimit,
                                          "homogeneous ODE allocation failed",
                                          operation);
    } catch (const std::length_error& error) {
        return ODESolutionResult::failure(CasErrc::ResourceLimit, error.what(), operation);
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

static ODESolutionResult solve_homogeneous_ode_impl(
    const std::shared_ptr<SymbolicExpr>& rhs,
    const std::string& x,
    const std::string& y, ComputationContext& context)
{

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
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
    auto f_v = detail::make_expression_ptr(detail::substitute_node(detail::node(rhs), x, detail::node(SymbolicExpr::number(1)), transform_budget));
    f_v = detail::make_expression_ptr(detail::substitute_node(detail::node(f_v), y, detail::node(v_var), transform_budget));
    f_v = detail::make_expression_ptr(detail::normalize_node(detail::node(f_v), transform_budget));

    /// 分离变量方程: dv/(f(v) - v) = dx/x
    /// 积分: integral dv/(f(v) - v) = integral dx/x = ln|x| + C
    auto f_minus_v = detail::make_expression_ptr(detail::normalize_node(detail::node(detail::make_expression_ptr(detail::expand_node(detail::node(SymbolicExpr::add(f_v,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), v_var))), transform_budget))
        ->cancel()), transform_budget));

    auto v_replacement = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(
        SymbolicExpr::variable(y), x_var)), transform_budget));
    if (f_minus_v->is_zero()) {
        result.general_solution = v_replacement;
        return ODESolutionResult::success(std::move(result));
    }

    std::shared_ptr<SymbolicExpr> integrand;
    const auto reciprocal =
        std::dynamic_pointer_cast<const PowerNode>(detail::node(f_minus_v));
    const auto minus_one = SymbolicExpr::number(-1);
    if (reciprocal &&
        reciprocal->exponent()->compare(*detail::node(minus_one)) == 0) {
        integrand = detail::make_expression_ptr(reciprocal->base());
    } else {
        integrand = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(
            SymbolicExpr::number(1), f_minus_v)), transform_budget));
    }
    Integrator integrator;
    auto integrated = integrator.integrate_checked(*integrand, v_name, context);
    if (!integrated) { return ODESolutionResult::failure(integrated.error()); }
    auto lhs_integral = std::make_shared<SymbolicExpr>(std::move(integrated.value()));

    /// 右端: ln(x)
    auto rhs_integral = SymbolicExpr::ln(x_var);

    /// 解为: lhs_integral = ln(x) + C
    /// 即: lhs_integral - ln(x) = C
    /// 回代 v = y/x

    auto solution = detail::make_expression_ptr(detail::substitute_node(detail::node(lhs_integral), v_name, detail::node(v_replacement), transform_budget));
    solution = SymbolicExpr::add(solution,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), rhs_integral));
    solution = detail::make_expression_ptr(detail::normalize_node(detail::node(solution), transform_budget));

    /// 返回隐式解形式: F(x,y) = C
    result.general_solution = solution;
    return ODESolutionResult::success(std::move(result));
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
    if (!valid) {
        return ODESolutionResult::failure(valid.error());
    }
    if (n == 0 || n == 1) {
        return ODESolutionResult::failure(CasErrc::InvalidArgument,
                                          "Bernoulli exponent must not be 0 or 1",
                                          operation);
    }
    auto step = context.consume_steps(14, operation);
    if (!step) {
        return ODESolutionResult::failure(step.error());
    }

    try {
        auto solved = solve_bernoulli_ode_impl(P, Q, n, x, y, context);
        if (!solved) { return solved; }
        return wrap_ode_solution(std::move(solved.value()), ODEType::Bernoulli, operation);
    } catch (const CasError& error) {
        return ODESolutionResult::failure(error);
    } catch (const std::bad_alloc&) {
        return ODESolutionResult::failure(CasErrc::ResourceLimit,
                                          "Bernoulli ODE allocation failed",
                                          operation);
    } catch (const std::length_error& error) {
        return ODESolutionResult::failure(CasErrc::ResourceLimit, error.what(), operation);
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

static ODESolutionResult solve_bernoulli_ode_impl(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    int n,
    const std::string& x,
    const std::string&, ComputationContext& context)
{

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
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
    auto P_linear = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(coeff, P)->cancel()), transform_budget));
    auto Q_linear = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(coeff, Q)->cancel()), transform_budget));

    /// 用积分因子法求解线性 ODE: v' + P_linear*v = Q_linear
    /// 积分因子 mu = exp(integralP_linear dx)
    Integrator integrator;
    auto integrated_p = integrator.integrate_checked(*P_linear, x, context);
    if (!integrated_p) { return ODESolutionResult::failure(integrated_p.error()); }
    auto intP = std::make_shared<SymbolicExpr>(std::move(integrated_p.value()));
    auto mu = integrating_factor_from_antiderivative(intP, x, context);
    auto Q_mu = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(Q_linear, mu)->cancel()), transform_budget));
    auto integrated_q = integrator.integrate_checked(*Q_mu, x, context);
    if (!integrated_q) { return ODESolutionResult::failure(integrated_q.error()); }
    auto int_Q_mu = std::make_shared<SymbolicExpr>(std::move(integrated_q.value()));

    auto C_const = SymbolicExpr::variable("C");
    auto numerator = SymbolicExpr::add(int_Q_mu, C_const);
    auto v_solution = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(numerator, mu)), transform_budget));

    /// 回代: v = y^(1-n),所以 y = v^(1/(1-n))
    /// y^(1-n) = v_solution
    /// y = v_solution^(1/(1-n))
    auto exponent = SymbolicExpr::divide(
        SymbolicExpr::number(1),
        SymbolicExpr::number(one_minus_n));
    auto y_solution = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::power(v_solution, exponent)), transform_budget));

    result.general_solution = y_solution;
    return ODESolutionResult::success(std::move(result));
}


static ExpressionResult integrating_factor_from_ratio(
    const std::shared_ptr<SymbolicExpr>& ratio,
    const std::string& variable, const std::string& excluded_variable,
    Integrator& integrator, ComputationContext& context) {
    if (expression_depends_on_variable(
            LMCAS::detail::node(ratio), excluded_variable)) {
        return std::shared_ptr<SymbolicExpr>{};
    }
    auto integrated = integrator.integrate_checked(*ratio, variable, context);
    if (!integrated) { return ExpressionResult::failure(integrated.error()); }
    return integrating_factor_from_antiderivative(
        std::make_shared<SymbolicExpr>(std::move(integrated.value())),
        variable, context);
}

static ExpressionResult find_integrating_factor_impl(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x, const std::string& y, ComputationContext& context)
{
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    if (!M || !N) {
        return std::shared_ptr<SymbolicExpr>{};
    }

    /// 计算 partialM/partialy - partialN/partialx
    auto dM_dy = detail::make_expression_ptr(detail::differentiate_node(detail::node(M), y, transform_budget));
    auto dN_dx = detail::make_expression_ptr(detail::differentiate_node(detail::node(N), x, transform_budget));

    auto diff = SymbolicExpr::add(dM_dy,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dN_dx));
    diff = detail::make_expression_ptr(detail::normalize_node(detail::node(diff), transform_budget));

    auto exact = check_zero_residual(diff, context);
    if (!exact) { return ExpressionResult::failure(exact.error()); }
    if (std::holds_alternative<ProvedZeroResidual>(exact.value())) {
        /// 已经恰当,积分因子为 1
        return SymbolicExpr::number(1);
    }

    Integrator integrator;
    /// 尝试 mu = mu(x): (partialM/partialy - partialN/partialx) / N 仅依赖 x
    if (N && !N->is_zero()) {
        auto ratio_x = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(diff, N)), transform_budget));
        auto factor = integrating_factor_from_ratio(ratio_x, x, y, integrator, context);
        if (!factor || factor.value()) { return factor; }
    }

    /// 尝试 mu = mu(y): (partialN/partialx - partialM/partialy) / M 仅依赖 y
    if (M && !M->is_zero()) {
        auto neg_diff = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(SymbolicExpr::number(-1), diff)), transform_budget));
        auto ratio_y = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(neg_diff, M)), transform_budget));
        auto factor = integrating_factor_from_ratio(ratio_y, y, x, integrator, context);
        if (!factor || factor.value()) { return factor; }
    }

    return std::shared_ptr<SymbolicExpr>{};
}

std::shared_ptr<SymbolicExpr> find_integrating_factor(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x, const std::string& y)
{
    ComputationContext context;
    auto result = find_integrating_factor_impl(M, N, x, y, context);
    if (!result) { throw result.error(); }
    return std::move(result.value());
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
    if (!valid) {
        return ODESolutionResult::failure(valid.error());
    }
    auto step = context.consume_steps(18, operation);
    if (!step) {
        return ODESolutionResult::failure(step.error());
    }

    try {
        auto solved = solve_exact_ode_impl(M, N, x, y, context);
        if (!solved) { return solved; }
        return wrap_ode_solution(std::move(solved.value()), ODEType::Exact, operation);
    } catch (const CasError& error) {
        return ODESolutionResult::failure(error);
    } catch (const std::bad_alloc&) {
        return ODESolutionResult::failure(CasErrc::ResourceLimit,
                                          "exact ODE allocation failed",
                                          operation);
    } catch (const std::length_error& error) {
        return ODESolutionResult::failure(CasErrc::ResourceLimit, error.what(), operation);
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

static ODESolutionResult solve_exact_ode_impl(
    const std::shared_ptr<SymbolicExpr>& M,
    const std::shared_ptr<SymbolicExpr>& N,
    const std::string& x,
    const std::string& y, ComputationContext& context)
{

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    ODESolution result;
    result.method_used = ODEType::Exact;
    result.constants = {"C"};

    auto M_eff = M;
    auto N_eff = N;

    auto factor = find_integrating_factor_impl(M, N, x, y, context);
    if (!factor) { return ODESolutionResult::failure(factor.error()); }
    if (!factor.value()) {
        return ODESolutionResult::success(std::move(result));
    }
    auto mu = std::move(factor.value());
    if (!mu->is_one()) {
        M_eff = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(mu, M)), transform_budget));
        N_eff = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(mu, N)), transform_budget));
    }

    /// 恰当方程: partialF/partialx = M, partialF/partialy = N
    /// F(x,y) = integral M dx + g(y)
    /// 对 F 关于 y 求导: partialF/partialy = partial/partialy(integral M dx) + g'(y) = N
    /// 所以 g'(y) = N - partial/partialy(integral M dx)

    /// 步骤 1: 计算 integral M dx(将 y 视为常数)
    Integrator integrator;
    auto integrated_m = integrator.integrate_checked(*M_eff, x, context);
    if (!integrated_m) { return ODESolutionResult::failure(integrated_m.error()); }
    auto F_partial = std::make_shared<SymbolicExpr>(std::move(integrated_m.value()));

    /// 步骤 2: 对 F_partial 关于 y 求导
    auto dF_dy = detail::make_expression_ptr(detail::differentiate_node(detail::node(F_partial), y, transform_budget));
    dF_dy = detail::make_expression_ptr(detail::normalize_node(detail::node(dF_dy), transform_budget));

    /// 步骤 3: g'(y) = N - dF_dy
    auto g_prime = SymbolicExpr::add(N_eff,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dF_dy));
    g_prime = detail::make_expression_ptr(detail::normalize_node(detail::node(g_prime), transform_budget));

    /// 步骤 4: g(y) = integral g'(y) dy
    auto integrated_g = integrator.integrate_checked(*g_prime, y, context);
    if (!integrated_g) { return ODESolutionResult::failure(integrated_g.error()); }
    auto g_y = std::make_shared<SymbolicExpr>(std::move(integrated_g.value()));

    /// 步骤 5: F(x,y) = F_partial + g(y)
    auto F_total = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::add(F_partial, g_y)), transform_budget));

    /// 解为 F(x,y) = C
    result.general_solution = F_total;
    return ODESolutionResult::success(std::move(result));
}

} // namespace LMCAS
