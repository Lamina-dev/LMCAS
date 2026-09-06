/**
 * @file ode_linear.cpp
 * @brief 常系数高阶与 Euler ODE 的解表达式构造和受检 API。
 */
#include "../include/symbolic_ode_engine.hpp"
#include "symbolic_ast.hpp"
#include "../include/symbolic.hpp"
#include "../include/poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "internal/ode_support.hpp"
#include "internal/ode_characteristic_roots.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace LMCAS {

static Result<ODESolution> solve_higher_order_ode_impl(
    const std::vector<double>&,
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&, const std::string&);
static Result<ODESolution> solve_euler_ode_impl(
    const std::vector<double>&,
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&, const std::string&);


using namespace ode_root_detail;

/**
 * @internal
 * @brief 将 double 值转为"干净"的数值表达式.
 *
 * 若值接近整数或简单分数,使用精确表示.
 */
static std::shared_ptr<SymbolicExpr> clean_number(double val) {
    /// 只在转换可表示时提升为精确小整数，避免浮点到整数的越界转换。
    double rounded = std::round(val);
    int eq = 0;
    if (std::isfinite(rounded) &&
        rounded >= static_cast<double>(std::numeric_limits<int>::min()) &&
        rounded <= static_cast<double>(std::numeric_limits<int>::max())) {
        lmmc_double_nearly_equal_tol(
            val, rounded, 1e-10, 1e-10, &eq);
        if (eq && (rounded != 0.0 || val == 0.0)) {
            return SymbolicExpr::number(static_cast<int>(rounded));
        }
    }

    /// 检查是否接近简单分数 p/q (q <= 12)。
    for (int q = 2; q <= 12; ++q) {
        double p = val * q;
        double p_rounded = std::round(p);
        if (!std::isfinite(p_rounded) ||
            p_rounded < static_cast<double>(
                std::numeric_limits<int>::min()) ||
            p_rounded > static_cast<double>(
                std::numeric_limits<int>::max())) {
            continue;
        }
        lmmc_double_nearly_equal_tol(
            p, p_rounded, 1e-10, 1e-10, &eq);
        if (eq && (p_rounded != 0.0 || p == 0.0)) {
            return SymbolicExpr::divide(
                SymbolicExpr::number(static_cast<int>(p_rounded)),
                SymbolicExpr::number(q));
        }
    }

    return SymbolicExpr::number(val);
}

static bool has_nonzero_forcing(const std::shared_ptr<SymbolicExpr>& forcing) {
    return forcing && LMCAS::detail::node(forcing) && !forcing->is_zero();
}

/**
 * @internal
 * @brief 根据特征根构造齐次通解.
 *
 * - 实根 r(重数 m):C_k * x^k * e^(rx),k = 0, ..., m-1
 * - 复根 alpha+/-betai(重数 m):x^k * e^(alphax) * (C_a*cos(betax) + C_b*sin(betax))
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

                /// e^(r*x) 因子。精确零才可省略；小的非零根仍改变解。
                if (root.real_part != 0.0) {
                    auto r_expr = clean_number(root.real_part);
                    auto exp_arg = SymbolicExpr::multiply(r_expr, x_var);
                    auto exp_term = SymbolicExpr::exp(exp_arg);
                    term = SymbolicExpr::multiply(term, exp_term);
                }

                solution = solution ? SymbolicExpr::add(solution, term) : term;
            } else {
                /// 复根 alpha+/-betai: 产生两个基本解
                /// C_a * x^k * e^(alphax) * cos(betax)
                /// C_b * x^k * e^(alphax) * sin(betax)
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

                /// e^(alphax) 因子。精确零才可省略。
                std::shared_ptr<SymbolicExpr> exp_factor = nullptr;
                if (root.real_part != 0.0) {
                    auto alpha_expr = clean_number(root.real_part);
                    auto exp_arg = SymbolicExpr::multiply(alpha_expr, x_var);
                    exp_factor = SymbolicExpr::exp(exp_arg);
                }

                /// x^k 因子
                std::shared_ptr<SymbolicExpr> x_pow_factor = nullptr;
                if (k > 0) {
                    x_pow_factor = SymbolicExpr::power(x_var, SymbolicExpr::number(k));
                }

                /// 构造 cos 项: Ca * x^k * e^(alphax) * cos(betax)
                auto term_cos = Ca;
                if (x_pow_factor) term_cos = SymbolicExpr::multiply(term_cos, x_pow_factor);
                if (exp_factor) term_cos = SymbolicExpr::multiply(term_cos, exp_factor);
                term_cos = SymbolicExpr::multiply(term_cos, cos_term);

                /// 构造 sin 项: Cb * x^k * e^(alphax) * sin(betax)
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
    for (double coefficient : coeffs) {
        if (!std::isfinite(coefficient / coeffs.front())) {
            return ODESolutionResult::failure(
                CasErrc::NumericFailure,
                "normalizing the characteristic polynomial exceeds the finite numeric domain",
                operation);
        }
    }

    auto budget = context.consume_steps(coeffs.size() * 20 + 20, operation);
    if (!budget) return ODESolutionResult::failure(budget.error());

    const bool nonzero_forcing = has_nonzero_forcing(forcing);
    if (nonzero_forcing &&
        (!forcing->is_number() || coeffs.back() == 0.0)) {
        return ODESolutionResult::failure(
            CasErrc::Inconclusive,
            "non-homogeneous constant-coefficient ODE currently requires a constant forcing and nonzero y coefficient",
            operation);
    }

    try {
        auto solution_result = solve_higher_order_ode_impl(
            coeffs, nullptr, x, y);
        if (!solution_result)
            return ODESolutionResult::failure(solution_result.error());
        auto solution = std::move(solution_result.value());
        if (nonzero_forcing && solution.general_solution) {
            auto particular = SymbolicExpr::divide(
                forcing, SymbolicExpr::number(coeffs.back()))->simplify();
            solution.general_solution = SymbolicExpr::add(
                solution.general_solution, particular)->simplify();
        }
        return wrap_ode_solution(
            std::move(solution),
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

static Result<ODESolution> solve_higher_order_ode_impl(
    const std::vector<double>& coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string&)
{
    const std::string operation = "solve_higher_order_ode";
    ODESolution result;
    result.method_used = ODEType::HigherOrder_ConstCoeff;

    if (coeffs.size() < 2 || has_nonzero_forcing(forcing)) {
        return Result<ODESolution>::failure(
            CasErrc::Inconclusive,
            "higher-order ODE is outside the implemented support domain",
            operation);
    }

    auto roots = find_characteristic_roots(coeffs, operation);
    if (!roots) return Result<ODESolution>::failure(roots.error());
    result.general_solution =
        build_homogeneous_solution(roots.value(), x, result.constants);
    if (result.general_solution) {
        result.general_solution = result.general_solution->simplify();
    }
    return Result<ODESolution>::success(std::move(result));
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

    const bool nonzero_forcing = has_nonzero_forcing(forcing);
    if (nonzero_forcing &&
        (!forcing->is_number() || euler_coeffs.back() == 0.0)) {
        return ODESolutionResult::failure(
            CasErrc::Inconclusive,
            "non-homogeneous Euler ODE currently requires a constant forcing and nonzero y coefficient",
            operation);
    }

    try {
        auto solution_result = solve_euler_ode_impl(
            euler_coeffs, nullptr, x, y);
        if (!solution_result)
            return ODESolutionResult::failure(solution_result.error());
        auto solution = std::move(solution_result.value());
        if (nonzero_forcing && solution.general_solution) {
            auto particular = SymbolicExpr::divide(
                forcing,
                SymbolicExpr::number(euler_coeffs.back()))->simplify();
            solution.general_solution = SymbolicExpr::add(
                solution.general_solution, particular)->simplify();
        }
        return wrap_ode_solution(
            std::move(solution),
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

static Result<ODESolution> solve_euler_ode_impl(
    const std::vector<double>& euler_coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string&)
{
    const std::string operation = "solve_euler_ode";
    ODESolution result;
    result.method_used = ODEType::Euler;

    if ((euler_coeffs.size() != 3 && euler_coeffs.size() != 4) ||
        has_nonzero_forcing(forcing)) {
        return Result<ODESolution>::failure(
            CasErrc::Inconclusive,
            "Euler ODE is outside the implemented support domain",
            operation);
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

    auto roots = find_characteristic_roots(characteristic, operation);
    if (!roots) return Result<ODESolution>::failure(roots.error());
    result.general_solution =
        build_euler_solution(roots.value(), x, result.constants);
    if (result.general_solution) {
        result.general_solution = result.general_solution->simplify();
    }
    return Result<ODESolution>::success(std::move(result));
}

/**
 * @internal
 * @brief 检测非齐次项的类型,用于待定系数法.
 *
 * 支持的类型:
 * - 多项式: x^n
 * - 指数: e^(ax)
 * - 三角: cos(bx), sin(bx)
 */

} // namespace LMCAS
