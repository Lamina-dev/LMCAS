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
#include "../include/integrator.hpp"
#include <cmath>
#include <memory>
#include <string>

namespace LMCAS {

static FrobeniusSolutionResult solve_frobenius_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::shared_ptr<SymbolicExpr>&,
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&, int, ODESingularityType, ComputationContext&);


static ODESolutionResult solve_variation_of_parameters_core(
    const std::shared_ptr<SymbolicExpr>& y1,
    const std::shared_ptr<SymbolicExpr>& y2,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& x,
    ComputationContext& context);

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
        auto solved = solve_variation_of_parameters_core(
            y1, y2, g, x, context);
        if (!solved) return solved;
        return wrap_ode_solution(
            std::move(solved.value()),
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

static ODESolutionResult solve_variation_of_parameters_core(
    const std::shared_ptr<SymbolicExpr>& y1,
    const std::shared_ptr<SymbolicExpr>& y2,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& x,
    ComputationContext& context)
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

    Integrator integrator;
    auto u1_value =
        integrator.integrate_checked(*u1_prime, x, context);
    if (!u1_value) {
        return ODESolutionResult::failure(u1_value.error());
    }
    auto u2_value =
        integrator.integrate_checked(*u2_prime, x, context);
    if (!u2_value) {
        return ODESolutionResult::failure(u2_value.error());
    }
    auto u1 = std::make_shared<SymbolicExpr>(
        std::move(u1_value.value()));
    auto u2 = std::make_shared<SymbolicExpr>(
        std::move(u2_value.value()));

    /// 特解: y_p = u₁·y₁ + u₂·y₂
    auto y_p = SymbolicExpr::add(
        SymbolicExpr::multiply(u1, y1),
        SymbolicExpr::multiply(u2, y2));
    y_p = y_p->simplify();

    result.general_solution = y_p;
    return result;
}



ODESingularityResult classify_singular_point_checked(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x,
    ComputationContext& context)
{
    const std::string operation = "classify_singular_point";
    auto valid = validate_ode_two_expr_point(p, q, x0, x, context, operation);
    if (!valid) return ODESingularityResult::failure(valid.error());
    try {
        auto p_at_x0 = p->substitute(x, x0)->simplify();
        auto q_at_x0 = q->substitute(x, x0)->simplify();
        const double p_val = try_eval_double(p_at_x0);
        const double q_val = try_eval_double(q_at_x0);
        if (std::isfinite(p_val) && std::isfinite(q_val)) {
            return ODESingularityResult::success(
                ODESingularityType::Ordinary);
        }
        if ((!std::isfinite(p_val) &&
             !free_variables(detail::node(p_at_x0)).empty()) ||
            (!std::isfinite(q_val) &&
             !free_variables(detail::node(q_at_x0)).empty())) {
            return ODESingularityResult::failure(
                CasErrc::Inconclusive,
                "singularity classification requires numeric coefficient values",
                operation);
        }

        auto x_var = SymbolicExpr::variable(x);
        auto x_minus_x0 = SymbolicExpr::add(
            x_var, SymbolicExpr::multiply(SymbolicExpr::number(-1), x0));
        auto xp = SymbolicExpr::multiply(x_minus_x0, p)->simplify();
        auto x2q = SymbolicExpr::multiply(
            SymbolicExpr::power(x_minus_x0, SymbolicExpr::number(2)),
            q)->simplify();
        auto xp_limit = limit_expression_checked(
            xp, x, x0, LimitDirection::Both, context);
        if (!xp_limit) {
            return ODESingularityResult::failure(xp_limit.error());
        }
        auto x2q_limit = limit_expression_checked(
            x2q, x, x0, LimitDirection::Both, context);
        if (!x2q_limit) {
            return ODESingularityResult::failure(x2q_limit.error());
        }
        const double xp_val = try_eval_double(xp_limit.value());
        const double x2q_val = try_eval_double(x2q_limit.value());
        return ODESingularityResult::success(
            std::isfinite(xp_val) && std::isfinite(x2q_val)
                ? ODESingularityType::RegularSingular
                : ODESingularityType::IrregularSingular);
    } catch (const std::bad_alloc&) {
        return ODESingularityResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while classifying ODE singularity",
            operation);
    } catch (const std::exception& ex) {
        return ODESingularityResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

ODESingularityResult classify_singular_point_checked(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x)
{
    ComputationContext context;
    return classify_singular_point_checked(p, q, x0, x, context);
}


static Result<void> validate_frobenius_regular_singular_domain(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x,
    ComputationContext& context,
    const std::string& operation)
{
    auto x_var = SymbolicExpr::variable(x);
    auto x_minus_x0 = SymbolicExpr::add(x_var,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x0));
    auto xp_expr = SymbolicExpr::multiply(x_minus_x0, p)->simplify();
    auto x2q_expr = SymbolicExpr::multiply(
        SymbolicExpr::power(x_minus_x0, SymbolicExpr::number(2)), q)->simplify();

    auto P0_result = limit_expression_checked(
        xp_expr, x, x0, LimitDirection::Both, context);
    if (!P0_result) return Result<void>::failure(P0_result.error());
    auto Q0_result = limit_expression_checked(
        x2q_expr, x, x0, LimitDirection::Both, context);
    if (!Q0_result) return Result<void>::failure(Q0_result.error());
    const auto& P0_expr = P0_result.value();
    const auto& Q0_expr = Q0_result.value();
    double P0 = P0_expr
        ? try_eval_double(P0_expr)
        : std::numeric_limits<double>::quiet_NaN();
    double Q0 = Q0_expr
        ? try_eval_double(Q0_expr)
        : std::numeric_limits<double>::quiet_NaN();
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
        auto classified =
            classify_singular_point_checked(p, q, x0, x, context);
        if (!classified) {
            return FrobeniusSolutionResult::failure(classified.error());
        }
        const auto point_type = classified.value();
        if (point_type == ODESingularityType::IrregularSingular) {
            return FrobeniusSolutionResult::failure(
                CasErrc::Inconclusive,
                "Frobenius checked API does not support irregular singular points",
                operation);
        }
        if (point_type == ODESingularityType::RegularSingular) {
            auto regular_domain = validate_frobenius_regular_singular_domain(
                p, q, x0, x, context, operation);
            if (!regular_domain) {
                return FrobeniusSolutionResult::failure(regular_domain.error());
            }
        }

        auto solved =
            solve_frobenius_impl(
                p, q, x0, x, order, point_type, context);
        if (!solved) return solved;
        auto solution = std::move(solved.value());
        if (!solution.series_solution ||
            !LMCAS::detail::node(solution.series_solution)) {
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

static FrobeniusSolutionResult solve_frobenius_impl(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x,
    int order,
    ODESingularityType point_type,
    ComputationContext& context)
{
    FrobeniusSolution result;
    result.truncation_order = order;
    result.point_type = point_type;

    auto unsupported_numeric = [] {
        return FrobeniusSolutionResult::failure(
            CasErrc::Inconclusive,
            "Frobenius series coefficients are outside the checked numeric support domain",
            "solve_frobenius");
    };
    auto arithmetic_failure = [] {
        return FrobeniusSolutionResult::failure(
            CasErrc::NumericFailure,
            "Frobenius coefficient recurrence produced a non-finite value",
            "solve_frobenius");
    };

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
            if (!std::isfinite(pv) || !std::isfinite(qv)) {
                return unsupported_numeric();
            }
            p_coeffs[k] = pv / factorial;
            q_coeffs[k] = qv / factorial;
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
            if (!std::isfinite(a[n_idx + 2])) {
                return arithmetic_failure();
            }
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

    auto P0_result = limit_expression_checked(
        xp_expr, x, x0, LimitDirection::Both, context);
    if (!P0_result) {
        return FrobeniusSolutionResult::failure(P0_result.error());
    }
    auto Q0_result = limit_expression_checked(
        x2q_expr, x, x0, LimitDirection::Both, context);
    if (!Q0_result) {
        return FrobeniusSolutionResult::failure(Q0_result.error());
    }
    auto P0_expr = std::move(P0_result.value());
    auto Q0_expr = std::move(Q0_result.value());

    double P0 = P0_expr
        ? try_eval_double(P0_expr)
        : std::numeric_limits<double>::quiet_NaN();
    double Q0 = Q0_expr
        ? try_eval_double(Q0_expr)
        : std::numeric_limits<double>::quiet_NaN();
    if (!std::isfinite(P0) || !std::isfinite(Q0)) {
        return unsupported_numeric();
    }

    /// 指标方程: r(r-1) + P₀·r + Q₀ = 0  →  r² + (P₀-1)·r + Q₀ = 0
    double ind_b = P0 - 1.0;
    double ind_c = Q0;
    double ind_D = ind_b * ind_b - 4.0 * ind_c;
    if (!std::isfinite(ind_D)) {
        return arithmetic_failure();
    }

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
            if (!std::isfinite(pv)) {
                auto limited = limit_expression_checked(
                    xp_current, x, x0, LimitDirection::Both, context);
                if (!limited) {
                    return FrobeniusSolutionResult::failure(
                        limited.error());
                }
                const auto& lim = limited.value();
                pv = lim
                    ? try_eval_double(lim)
                    : std::numeric_limits<double>::quiet_NaN();
            }
            if (!std::isfinite(qv)) {
                auto limited = limit_expression_checked(
                    x2q_current, x, x0, LimitDirection::Both, context);
                if (!limited) {
                    return FrobeniusSolutionResult::failure(
                        limited.error());
                }
                const auto& lim = limited.value();
                qv = lim
                    ? try_eval_double(lim)
                    : std::numeric_limits<double>::quiet_NaN();
            }
        }
        if (!std::isfinite(pv) || !std::isfinite(qv)) {
            return unsupported_numeric();
        }
        pn_coeffs[k] = pv / fact;
        qn_coeffs[k] = qv / fact;

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
        if (!std::isfinite(F_val)) {
            return arithmetic_failure();
        }
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
        if (!std::isfinite(a[n_idx])) {
            return arithmetic_failure();
        }
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

} // namespace LMCAS
