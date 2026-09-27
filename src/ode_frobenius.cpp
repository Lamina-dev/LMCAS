#include "limit_result.hpp"
/**
 * @file ode_frobenius.cpp
 * @brief 奇点分类与 Frobenius 系数递推。
 */
#include "symbolic_ode_engine.hpp"
#include "internal/symbolic_ast.hpp"
#include "symbolic.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "internal/ode_support.hpp"
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace LMCAS {

static FrobeniusSolutionResult solve_frobenius_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::shared_ptr<SymbolicExpr>&,
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&, int, ODESingularityType, ComputationContext&);

static CasError frobenius_unsupported_numeric() {
    return {CasErrc::Inconclusive,
            "Frobenius series coefficients are outside the checked numeric support domain",
            "solve_frobenius"};
}

static CasError frobenius_arithmetic_failure() {
    return {CasErrc::NumericFailure,
            "Frobenius coefficient recurrence produced a non-finite value",
            "solve_frobenius"};
}

struct OdeTaylorCoefficients {
    std::vector<double> p;
    std::vector<double> q;
};

static Result<OdeTaylorCoefficients> ordinary_taylor_coefficients(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x, int order) {
        std::vector<double> p_coeffs(order + 2, 0.0);
        std::vector<double> q_coeffs(order + 2, 0.0);

        auto p_current = p;
        auto q_current = q;
        double factorial = 1.0;
        for (int k = 0; k <= order + 1; ++k) {
            if (k > 0) {
                factorial *= k;
            }
            auto p_val_expr = p_current->substitute(x, x0)->simplify();
            auto q_val_expr = q_current->substitute(x, x0)->simplify();
            double pv = try_eval_double(p_val_expr);
            double qv = try_eval_double(q_val_expr);
            if (!std::isfinite(pv) || !std::isfinite(qv)) {
                return Result<OdeTaylorCoefficients>::failure(frobenius_unsupported_numeric());
            }
            p_coeffs[k] = pv / factorial;
            q_coeffs[k] = qv / factorial;
            p_current = p_current->differentiate(x);
            q_current = q_current->differentiate(x);
            if (!p_current) {
                p_current = SymbolicExpr::number(0);
            }
            if (!q_current) {
                q_current = SymbolicExpr::number(0);
            }
        }

    return OdeTaylorCoefficients{std::move(p_coeffs), std::move(q_coeffs)};
}

static Result<std::vector<double>> ordinary_series_recurrence(
    const OdeTaylorCoefficients& coefficients, int order) {
    const auto& p_coeffs = coefficients.p;
    const auto& q_coeffs = coefficients.q;
        std::vector<double> a(order + 1, 0.0);
        a[0] = 1.0;
        for (int n_idx = 0; n_idx + 2 <= order; ++n_idx) {
            double sum = 0.0;
            for (int k = 0; k <= n_idx; ++k) {
                sum += (k + 1) * a[k + 1] * p_coeffs[n_idx - k];
                sum += a[k] * q_coeffs[n_idx - k];
            }
            double denom = static_cast<double>((n_idx + 2) * (n_idx + 1));
            a[n_idx + 2] = -sum / denom;
            if (!std::isfinite(a[n_idx + 2])) {
                return Result<std::vector<double>>::failure(frobenius_arithmetic_failure());
            }
        }

    return a;
}

static std::shared_ptr<SymbolicExpr> assemble_frobenius_polynomial(
    const std::vector<double>& a,
    const std::shared_ptr<SymbolicExpr>& x_minus_x0) {
    const int order = static_cast<int>(a.size()) - 1;
        auto series_sol = SymbolicExpr::number(0);
        for (int k = 0; k <= order; ++k) {
            if (std::abs(a[k]) < 1e-15) {
                continue;
            }
            auto coeff = SymbolicExpr::number(a[k]);
            auto power_term = (k == 0) ? SymbolicExpr::number(1)
                : SymbolicExpr::power(x_minus_x0, SymbolicExpr::number(k));
            series_sol = SymbolicExpr::add(series_sol,
                SymbolicExpr::multiply(coeff, power_term));
        }

    return series_sol;
}

static Result<std::pair<double, double>> frobenius_indicial_roots(
    double P0, double Q0) {
    double ind_b = P0 - 1.0;
    double ind_c = Q0;
    double ind_D = ind_b * ind_b - 4.0 * ind_c;
    if (!std::isfinite(ind_D)) {
        return Result<std::pair<double, double>>::failure(frobenius_arithmetic_failure());
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
    if (r1 < r2) {
        std::swap(r1, r2);
    }

    return std::make_pair(r1, r2);
}

static Result<double> frobenius_coefficient_limit(
    const std::shared_ptr<SymbolicExpr>& coefficient,
    const std::string& x, const std::shared_ptr<SymbolicExpr>& x0,
    double value, ComputationContext& context) {
    if (std::isfinite(value)) {
        return value;
    }
    auto limited = limit_expression_checked(
        coefficient, x, x0, LimitDirection::Both, context);
    if (!limited) {
        return Result<double>::failure(limited.error());
    }
    const auto& lim = limited.value();
    return lim ? try_eval_double(lim)
               : std::numeric_limits<double>::quiet_NaN();
}

static Result<OdeTaylorCoefficients> singular_taylor_coefficients(
    const std::shared_ptr<SymbolicExpr>& xp_expr,
    const std::shared_ptr<SymbolicExpr>& x2q_expr,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x, int order, double P0, double Q0,
    ComputationContext& context) {
    std::vector<double> pn_coeffs(order + 1, 0.0);
    std::vector<double> qn_coeffs(order + 1, 0.0);

    auto xp_current = xp_expr;
    auto x2q_current = x2q_expr;
    double fact = 1.0;
    for (int k = 0; k <= order; ++k) {
        if (k > 0) {
            fact *= k;
        }
        double pv, qv;
        if (k == 0) {
            pv = P0;
            qv = Q0;
        } else {
            auto pv_expr = xp_current->substitute(x, x0)->simplify();
            auto qv_expr = x2q_current->substitute(x, x0)->simplify();
            pv = try_eval_double(pv_expr);
            qv = try_eval_double(qv_expr);

            auto p_limit = frobenius_coefficient_limit(xp_current, x, x0, pv, context);
            if (!p_limit) {
                return Result<OdeTaylorCoefficients>::failure(p_limit.error());
            }
            pv = p_limit.value();
            auto q_limit = frobenius_coefficient_limit(x2q_current, x, x0, qv, context);
            if (!q_limit) {
                return Result<OdeTaylorCoefficients>::failure(q_limit.error());
            }
            qv = q_limit.value();
        }
        if (!std::isfinite(pv) || !std::isfinite(qv)) {
            return Result<OdeTaylorCoefficients>::failure(frobenius_unsupported_numeric());
        }
        pn_coeffs[k] = pv / fact;
        qn_coeffs[k] = qv / fact;

        xp_current = xp_current->differentiate(x);
        x2q_current = x2q_current->differentiate(x);
        if (!xp_current) {
            xp_current = SymbolicExpr::number(0);
        }
        if (!x2q_current) {
            x2q_current = SymbolicExpr::number(0);
        }
    }

    return OdeTaylorCoefficients{std::move(pn_coeffs), std::move(qn_coeffs)};
}

static Result<std::vector<double>> singular_series_recurrence(
    const OdeTaylorCoefficients& coefficients,
    double P0, double Q0, double r1, int order) {
    const auto& pn_coeffs = coefficients.p;
    const auto& qn_coeffs = coefficients.q;
    /**
     * @brief 递推 a_n = -1/F(r1+n) * ∑_{k=0}^{n-1} [(r1+k)·p_{n-k} + q_{n-k}]·a_k。
     * @note F(s) = s(s-1) + P₀·s + Q₀。
     */
    auto indicial_poly = [&](double s) -> double {
        return s * (s - 1.0) + P0 * s + Q0;
    };

    std::vector<double> a(order + 1, 0.0);
    a[0] = 1.0;

    for (int n_idx = 1; n_idx <= order; ++n_idx) {
        double F_val = indicial_poly(r1 + n_idx);
        if (!std::isfinite(F_val)) {
            return Result<std::vector<double>>::failure(frobenius_arithmetic_failure());
        }
        if (std::abs(F_val) < 1e-15) {
            a[n_idx] = 0.0;
            continue;
        }
        double sum = 0.0;
        for (int k = 0; k < n_idx; ++k) {
            int idx = n_idx - k;
            if (idx > order) {
                continue;
            }
            double p_term = (idx < static_cast<int>(pn_coeffs.size())) ? pn_coeffs[idx] : 0.0;
            double q_term = (idx < static_cast<int>(qn_coeffs.size())) ? qn_coeffs[idx] : 0.0;
            sum += ((r1 + k) * p_term + q_term) * a[k];
        }
        a[n_idx] = -sum / F_val;
        if (!std::isfinite(a[n_idx])) {
            return Result<std::vector<double>>::failure(frobenius_arithmetic_failure());
        }
    }

    return a;
}

static FrobeniusSolutionResult ordinary_frobenius_solution(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x, const std::shared_ptr<SymbolicExpr>& x_minus_x0,
    FrobeniusSolution result) {
    auto coefficients = ordinary_taylor_coefficients(p, q, x0, x, result.truncation_order);
    if (!coefficients) {
        return FrobeniusSolutionResult::failure(coefficients.error());
    }
    auto recurrence = ordinary_series_recurrence(coefficients.value(), result.truncation_order);
    if (!recurrence) {
        return FrobeniusSolutionResult::failure(recurrence.error());
    }
    result.series_solution = assemble_frobenius_polynomial(
        recurrence.value(), x_minus_x0)->simplify();
    result.indicial_roots = {0.0};
    return result;
}

static FrobeniusSolutionResult singular_frobenius_solution(
    const std::shared_ptr<SymbolicExpr>& p,
    const std::shared_ptr<SymbolicExpr>& q,
    const std::shared_ptr<SymbolicExpr>& x0,
    const std::string& x, const std::shared_ptr<SymbolicExpr>& x_minus_x0,
    FrobeniusSolution result, ComputationContext& context) {
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
        return FrobeniusSolutionResult::failure(frobenius_unsupported_numeric());
    }
    /**
     * @brief 在原定义域内确认去心极限有限后约分，提取可去奇点解析延拓的 Taylor 系数。
     * @note 原乘积未经约分的导数仍含缺失点，不能直接用于系数计算。
     */
    xp_expr = xp_expr->cancel()->simplify();
    x2q_expr = x2q_expr->cancel()->simplify();

    auto roots = frobenius_indicial_roots(P0, Q0);
    if (!roots) {
        return FrobeniusSolutionResult::failure(roots.error());
    }
    const double r1 = roots.value().first;
    result.indicial_roots = {r1, roots.value().second};
    auto coefficients = singular_taylor_coefficients(
        xp_expr, x2q_expr, x0, x, result.truncation_order, P0, Q0, context);
    if (!coefficients) {
        return FrobeniusSolutionResult::failure(coefficients.error());
    }
    auto recurrence = singular_series_recurrence(
        coefficients.value(), P0, Q0, r1, result.truncation_order);
    if (!recurrence) {
        return FrobeniusSolutionResult::failure(recurrence.error());
    }
    auto power_prefix = SymbolicExpr::power(x_minus_x0, SymbolicExpr::number(r1));
    auto series_part = assemble_frobenius_polynomial(recurrence.value(), x_minus_x0);
    result.series_solution = SymbolicExpr::multiply(power_prefix, series_part)->simplify();
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
    if (!valid) {
        return ODESingularityResult::failure(valid.error());
    }
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
        auto xp_limit = limit_checked(xp, x, x0, LimitDirection::Both, context);
        if (!xp_limit) return ODESingularityResult::failure(xp_limit.error());
        auto x2q_limit = limit_checked(x2q, x, x0, LimitDirection::Both, context);
        if (!x2q_limit) return ODESingularityResult::failure(x2q_limit.error());
        return ODESingularityResult::success(
            std::holds_alternative<FiniteLimit>(xp_limit.value().value) &&
            std::holds_alternative<FiniteLimit>(x2q_limit.value().value)
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
    if (!P0_result) {
        return Result<void>::failure(P0_result.error());
    }
    auto Q0_result = limit_expression_checked(
        x2q_expr, x, x0, LimitDirection::Both, context);
    if (!Q0_result) {
        return Result<void>::failure(Q0_result.error());
    }
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
    if (!valid) {
        return FrobeniusSolutionResult::failure(valid.error());
    }
    if (order < 0 || order > 64) {
        return FrobeniusSolutionResult::failure(
            CasErrc::InvalidArgument,
            "Frobenius truncation order must be between 0 and 64",
            operation);
    }
    auto budget = context.consume_steps(static_cast<std::size_t>(order + 1) * 24 + 24,
                                        operation);
    if (!budget) {
        return FrobeniusSolutionResult::failure(budget.error());
    }

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
        if (!solved) {
            return solved;
        }
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
    auto x_var = SymbolicExpr::variable(x);
    auto x_minus_x0 = SymbolicExpr::add(x_var,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x0));

    if (result.point_type == ODESingularityType::IrregularSingular) {
        result.series_solution = nullptr;
        return result;
    }

    if (result.point_type == ODESingularityType::Ordinary) {
        return ordinary_frobenius_solution(p, q, x0, x, x_minus_x0, std::move(result));
    }
    return singular_frobenius_solution(
        p, q, x0, x, x_minus_x0, std::move(result), context);
}

}
