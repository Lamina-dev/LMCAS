#include "internal/polynomial_solver_support.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace LMCAS::polynomial_solver_detail {

static std::optional<NumericQuarticDepression> depress_numeric_quartic(
    double av, double bv, double cv, double dv, double ev) {

    double shift_val = bv / (4.0 * av);
    double p_val = (8.0*av*cv - 3.0*bv*bv) / (8.0*av*av);
    double q_val = (bv*bv*bv - 4.0*av*bv*cv + 8.0*av*av*dv) / (8.0*av*av*av);
    double r_val = (-3.0*bv*bv*bv*bv + 256.0*av*av*av*ev - 64.0*av*av*bv*dv + 16.0*av*bv*bv*cv) / (256.0*av*av*av*av);
    if (!std::isfinite(shift_val) || !std::isfinite(p_val) ||
        !std::isfinite(q_val) || !std::isfinite(r_val)) {
        return std::nullopt;
    }
    return NumericQuarticDepression{p_val, q_val, r_val, shift_val};
}

static bool quartic_coefficients_are_numeric(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e) {
    return is_purely_numeric(a) && is_purely_numeric(b) &&
        is_purely_numeric(c) && is_purely_numeric(d) && is_purely_numeric(e);
}

std::optional<NumericQuarticDepression> numeric_quartic_depression(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e) {
    const bool all_numeric = quartic_coefficients_are_numeric(a, b, c, d, e);
    auto av = finite_numeric_value(a);
    auto bv = finite_numeric_value(b);
    auto cv = finite_numeric_value(c);
    auto dv = finite_numeric_value(d);
    auto ev = finite_numeric_value(e);
    if (!all_numeric || !av || !bv || !cv || !dv || !ev || *av == 0.0) {
        return std::nullopt;
    }
    return depress_numeric_quartic(*av, *bv, *cv, *dv, *ev);
}

static std::optional<double> select_numeric_resolvent(
    const std::vector<std::shared_ptr<SymbolicExpr>>& cubic_roots) {
    double m_val = 0.0;
    bool found_m = false;
    for (const auto& root : cubic_roots) {
        auto maybe_val = finite_numeric_value(root);
        if (maybe_val && *maybe_val > 0.0 &&
            (!found_m || *maybe_val > m_val)) {
            m_val = *maybe_val;
            found_m = true;
        }
    }
    if (!found_m) {
        for (const auto& root : cubic_roots) {
            auto maybe_val = finite_numeric_value(root);
            if (maybe_val && *maybe_val != 0.0 &&
                (!found_m ||
                 std::abs(*maybe_val) > std::abs(m_val))) {
                m_val = *maybe_val;
                found_m = true;
            }
        }
    }
    if (!found_m) {
        return std::nullopt;
    }
    return m_val;
}

static std::vector<std::shared_ptr<SymbolicExpr>> shifted_biquadratic_roots(
    double p_val, double r_val, double shift_val) {
    auto u_roots = solve_quadratic_internal(
        num(1),
        SymbolicExpr::number(p_val),
        SymbolicExpr::number(r_val));
    if (u_roots.size() == 1) {
        u_roots.push_back(u_roots.front());
    }

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    auto shift_expr = SymbolicExpr::number(shift_val);
    for (const auto& u_root : u_roots) {
        auto square_root = SymbolicExpr::sqrt(u_root);
        results.push_back(sub(square_root, shift_expr)->simplify());
        results.push_back(
            sub(negate(square_root), shift_expr)->simplify());
    }
    return results;
}

static std::vector<std::shared_ptr<SymbolicExpr>> negative_resolvent_roots(
    double m_val, double p_normalized, double q_normalized,
    double root_scale, double shift_val) {
    auto shift_expr = SymbolicExpr::number(shift_val);
    auto scale_expr = SymbolicExpr::number(root_scale);
    auto p_expr = SymbolicExpr::number(p_normalized);
    auto q_expr = SymbolicExpr::number(q_normalized);
    auto s_expr = SymbolicExpr::sqrt(
        SymbolicExpr::number(2.0 * m_val));
    auto m_plus_p_half = SymbolicExpr::number(
        m_val + p_normalized / 2.0);
    auto q_over_2s = SymbolicExpr::divide(
        q_expr, SymbolicExpr::multiply(num(2), s_expr));
    auto quad1_c_expr =
        sub(m_plus_p_half, q_over_2s)->simplify();
    auto quad2_c_expr =
        SymbolicExpr::add(
            m_plus_p_half, q_over_2s)->simplify();
    auto y_roots1 = solve_quadratic_internal(
        num(1), s_expr, quad1_c_expr);
    auto y_roots2 = solve_quadratic_internal(
        num(1), negate(s_expr), quad2_c_expr);

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (const auto& y : y_roots1) {
        auto scaled_y =
            SymbolicExpr::multiply(scale_expr, y)->simplify();
        results.push_back(
            sub(scaled_y, shift_expr)->simplify());
    }
    for (const auto& y : y_roots2) {
        auto scaled_y =
            SymbolicExpr::multiply(scale_expr, y)->simplify();
        results.push_back(
            sub(scaled_y, shift_expr)->simplify());
    }
    return results;
}

static void append_resolvent_quadratic(
    std::vector<std::shared_ptr<SymbolicExpr>>& results,
    double linear, double constant, double root_scale, double shift_val) {
    const double discriminant = linear * linear - 4.0 * constant;
    if (discriminant >= 0.0) {
        const double square_root = std::sqrt(discriminant);
        const double first = (-linear + square_root) / 2.0;
        const double second = (-linear - square_root) / 2.0;
        results.push_back(SymbolicExpr::number(first * root_scale - shift_val));
        results.push_back(SymbolicExpr::number(second * root_scale - shift_val));
    } else {
        const double real = -linear * root_scale / 2.0 - shift_val;
        const double imaginary = std::sqrt(-discriminant) * root_scale / 2.0;
        auto unit = SymbolicExpr::sqrt(num(-1));
        results.push_back(SymbolicExpr::add(SymbolicExpr::number(real),
            SymbolicExpr::multiply(SymbolicExpr::number(imaginary), unit))->simplify());
        results.push_back(sub(SymbolicExpr::number(real),
            SymbolicExpr::multiply(SymbolicExpr::number(imaginary), unit))->simplify());
    }
}

static std::vector<std::shared_ptr<SymbolicExpr>> positive_resolvent_roots(
    double m_val, double p_normalized, double q_normalized,
    double root_scale, double shift_val) {
    const double s_val = std::sqrt(2.0 * m_val);
    const double quad1_c_val =
        m_val + p_normalized / 2.0 -
        q_normalized / (2.0 * s_val);
    const double quad2_c_val =
        m_val + p_normalized / 2.0 +
        q_normalized / (2.0 * s_val);
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    append_resolvent_quadratic(results, s_val, quad1_c_val, root_scale, shift_val);
    append_resolvent_quadratic(results, -s_val, quad2_c_val, root_scale, shift_val);
    return results;
}

std::vector<std::shared_ptr<SymbolicExpr>> numeric_quartic_roots(
    const NumericQuarticDepression& depression, const std::string& var) {
    const auto& [p_val, q_val, r_val, shift_val] = depression;
    if (q_val == 0.0) {
        return shifted_biquadratic_roots(p_val, r_val, shift_val);
    }
    const double root_scale = std::max({
        std::sqrt(std::abs(p_val)),
        std::cbrt(std::abs(q_val)),
        std::sqrt(std::sqrt(std::abs(r_val)))
    });
    const double p_normalized =
        (p_val / root_scale) / root_scale;
    const double q_normalized =
        ((q_val / root_scale) / root_scale) / root_scale;
    const double r_normalized =
        (((r_val / root_scale) / root_scale) / root_scale) / root_scale;
    if (q_normalized == 0.0) {
        return {};
    }
    auto cubic_roots = solve_cubic(
        SymbolicExpr::number(8.0),
        SymbolicExpr::number(8.0 * p_normalized),
        SymbolicExpr::number(
            2.0 * p_normalized * p_normalized -
            8.0 * r_normalized),
        SymbolicExpr::number(-(q_normalized * q_normalized)),
        var);
    auto selected = select_numeric_resolvent(cubic_roots);
    if (!selected || *selected == 0.0) {
        return {};
    }
    if (*selected < 0.0) {
        return negative_resolvent_roots(
            *selected, p_normalized, q_normalized, root_scale, shift_val);
    }
    return positive_resolvent_roots(
        *selected, p_normalized, q_normalized, root_scale, shift_val);
}

}
