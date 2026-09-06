#include "solve_polynomial.hpp"
#include "symbolic_ast.hpp"
#include "root_of_utils.hpp"
#include "poly_utils.hpp"
#include "numeric_evaluation.hpp"
#include "internal/numeric_probe.hpp"
#include "internal/expression_analysis.hpp"
#include <cmath>
#include <algorithm>
#include <optional>
#include <stdexcept>
#include <limits>

namespace LMCAS {

static bool is_purely_numeric(const std::shared_ptr<SymbolicExpr>& expr) {
    return !expr || free_variables(LMCAS::detail::node(expr)).empty();
}

static std::shared_ptr<SymbolicExpr> cbrt_expr(const std::shared_ptr<SymbolicExpr>& x) {
    return SymbolicExpr::power(x, SymbolicExpr::number(Rational(1, 3)));
}

[[maybe_unused]] static std::shared_ptr<SymbolicExpr> acos_expr(const std::shared_ptr<SymbolicExpr>& x) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::ArcCos,
            std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(x)}));
}

static std::shared_ptr<SymbolicExpr> negate(const std::shared_ptr<SymbolicExpr>& x) {
    return SymbolicExpr::multiply(SymbolicExpr::number(-1), x);
}

static std::shared_ptr<SymbolicExpr> sub(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    return SymbolicExpr::add(a, negate(b));
}

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }
[[maybe_unused]] static std::shared_ptr<SymbolicExpr> num_r(int p, int q) { return SymbolicExpr::number(Rational(p, q)); }

static std::optional<double> finite_numeric_value(
    const std::shared_ptr<SymbolicExpr>& expr) {
    auto evaluated = detail::try_finite_numeric(expr);
    if (!evaluated) return std::nullopt;
    auto simplified = expr->simplify();
    if (simplified && LMCAS::detail::node(simplified)) {
        if (auto num = std::dynamic_pointer_cast<const NumberNode>(
                LMCAS::detail::node(simplified))) {
            if (std::holds_alternative<Rational>(num->value())) {
                const auto& value = std::get<Rational>(num->value());
                if (!value.get_numerator().is_zero() && *evaluated == 0.0) {
                    return std::nullopt;
                }
            } else if (std::holds_alternative<BigInt>(num->value())) {
                const auto& value = std::get<BigInt>(num->value());
                if (!value.is_zero() && *evaluated == 0.0) {
                    return std::nullopt;
                }
            }
        }
    }
    return evaluated;
}

static std::vector<std::shared_ptr<SymbolicExpr>> solve_quadratic_internal(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c);

static bool convert_to_rational_poly(
    const Polynomial<SymbolicPolyCoeff>& sym_poly,
    Polynomial<Rational>& out_poly);

static std::vector<std::shared_ptr<SymbolicExpr>> solve_linear_internal(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b);

static std::vector<std::shared_ptr<SymbolicExpr>> solve_linear_internal(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b) {

    auto neg_b = negate(b);
    return { SymbolicExpr::divide(neg_b, a)->simplify() };
}

static std::vector<std::shared_ptr<SymbolicExpr>> solve_quadratic_internal(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c) {
    auto a_checked = finite_numeric_value(a);
    auto b_checked = finite_numeric_value(b);
    auto c_checked = finite_numeric_value(c);
    if (is_purely_numeric(a) && is_purely_numeric(b) &&
        is_purely_numeric(c) && a_checked && b_checked && c_checked &&
        *a_checked != 0.0) {
        const double av = *a_checked;
        const double bv = *b_checked;
        const double cv = *c_checked;
        const double scale =
            std::max({std::abs(av), std::abs(bv), std::abs(cv)});
        const double an = av / scale;
        const double bn = bv / scale;
        const double cn = cv / scale;
        const bool normalization_lost_coefficient =
            (av != 0.0 && an == 0.0) ||
            (bv != 0.0 && bn == 0.0) ||
            (cv != 0.0 && cn == 0.0);

        if (!normalization_lost_coefficient) {
            const double discriminant =
                std::fma(-4.0 * an, cn, bn * bn);
            if (std::isfinite(discriminant) && discriminant >= 0.0) {
                const double sqrt_discriminant = std::sqrt(discriminant);
                if (sqrt_discriminant == 0.0) {
                    return {
                        SymbolicExpr::number(-bn / (2.0 * an))
                    };
                }

                const double q =
                    -0.5 * (bn + std::copysign(sqrt_discriminant, bn));
                const double x1 = q / an;
                const double x2 = cn / q;
                const bool representable =
                    std::isfinite(x1) && std::isfinite(x2) &&
                    (cv == 0.0 || (x1 != 0.0 && x2 != 0.0));
                if (representable) {
                    return {
                        SymbolicExpr::number(x1),
                        SymbolicExpr::number(x2)
                    };
                }
            }
        }
    }


    auto b2 = SymbolicExpr::power(b, num(2));
    auto four_ac = SymbolicExpr::multiply(num(4), SymbolicExpr::multiply(a, c));
    auto delta = sub(b2, four_ac)->simplify();

    auto neg_b = negate(b);
    auto two_a = SymbolicExpr::multiply(num(2), a);
    if (delta && delta->is_zero()) {
        return { SymbolicExpr::divide(neg_b, two_a)->simplify() };
    }

    auto sqrt_delta = SymbolicExpr::sqrt(delta);
    auto x1 = SymbolicExpr::divide(SymbolicExpr::add(neg_b, sqrt_delta), two_a)->simplify();
    auto x2 = SymbolicExpr::divide(sub(neg_b, sqrt_delta), two_a)->simplify();

    return { x1, x2 };
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_cubic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::string&) {

    auto a_simp = a->simplify();
    if (a_simp->get_number_value_is_zero()) {
        auto b_simp = b->simplify();
        if (b_simp->get_number_value_is_zero()) {
            auto c_simp = c->simplify();
            if (c_simp->get_number_value_is_zero()) {
                auto d_simp = d->simplify();
                if (d_simp->get_number_value_is_zero()) {
                    throw std::invalid_argument(
                        "solve_cubic: identically zero equation has no finite root list");
                }
                if (is_purely_numeric(d_simp)) {
                    return {};
                }
                throw std::invalid_argument(
                    "solve_cubic: equation is independent of the requested variable");
            }
            return solve_linear_internal(c, d);
        }
        return solve_quadratic_internal(b, c, d);
    }


    bool all_numeric = is_purely_numeric(a) && is_purely_numeric(b) &&
                       is_purely_numeric(c) && is_purely_numeric(d);

    std::vector<std::shared_ptr<SymbolicExpr>> roots;
    double p_val = 0.0;
    double q_val = 0.0;
    double shift_val = 0.0;
    bool numeric_depression = false;

    auto a_checked = finite_numeric_value(a);
    auto b_checked = finite_numeric_value(b);
    auto c_checked = finite_numeric_value(c);
    auto d_checked = finite_numeric_value(d);
    if (all_numeric && a_checked && b_checked && c_checked && d_checked &&
        *a_checked != 0.0) {
        const double max_coefficient = std::max({
            std::abs(*a_checked), std::abs(*b_checked),
            std::abs(*c_checked), std::abs(*d_checked)});
        /**
         * Common binary scaling places the largest coefficient in [1, 2).
         * Normal-range scaled values retain their binary significands.
         * @see David Goldberg, "What Every Computer Scientist Should Know About
         * Floating-Point Arithmetic" (1991), Theorem 7 proof.
         * https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
         * @see S. Ghaderpanah and S. Klasa, "Polynomial Scaling" (1990),
         * SIAM Journal on Numerical Analysis 27(1), 117-135 (scaling background).
         * https://doi.org/10.1137/0727007
         */
        const int coefficient_exponent = std::ilogb(max_coefficient);
        const double an = std::scalbn(*a_checked, -coefficient_exponent);
        const double bn = std::scalbn(*b_checked, -coefficient_exponent);
        const double cn = std::scalbn(*c_checked, -coefficient_exponent);
        const double dn = std::scalbn(*d_checked, -coefficient_exponent);
        if (an != 0.0) {
            const double an_squared = an * an;
            const double an_cubed = an_squared * an;
            const double bn_squared = bn * bn;
            const double p_left = 3.0 * an * cn;
            double p_numerator = std::fma(-bn, bn, p_left);
            const double p_roundoff =
                16.0 * std::numeric_limits<double>::epsilon() *
                (std::abs(p_left) + std::abs(bn_squared));
            if (std::abs(p_numerator) <= p_roundoff) {
                p_numerator = 0.0;
            }

            const double q_first = 2.0 * bn_squared * bn;
            const double q_second = -9.0 * an * bn * cn;
            const double q_third = 27.0 * an_squared * dn;
            double q_numerator = (q_first + q_second) + q_third;
            const double q_roundoff =
                32.0 * std::numeric_limits<double>::epsilon() *
                (std::abs(q_first) + std::abs(q_second) +
                 std::abs(q_third));
            if (std::abs(q_numerator) <= q_roundoff) {
                q_numerator = 0.0;
            }

            const double p_denominator = 3.0 * an_squared;
            const double q_denominator = 27.0 * an_cubed;
            if (p_denominator != 0.0 && q_denominator != 0.0) {
                p_val = p_numerator / p_denominator;
                q_val = q_numerator / q_denominator;
                shift_val = bn / (3.0 * an);
                numeric_depression =
                    std::isfinite(p_val) && std::isfinite(q_val) &&
                    std::isfinite(shift_val);
            }
        }
    }

    if (numeric_depression) {
        const double root_scale = std::max(
            std::sqrt(std::abs(p_val)),
            std::cbrt(std::abs(q_val)));

        if (root_scale == 0.0) {
            auto root = SymbolicExpr::number(-shift_val);
            roots = { root, root, root };
        } else {
            const double p_normalized = (p_val / root_scale) / root_scale;
            const double q_normalized =
                ((q_val / root_scale) / root_scale) / root_scale;
            const double q_half = q_normalized / 2.0;
            const double p_third = p_normalized / 3.0;
            const double discriminant_q = q_half * q_half;
            const double discriminant_p = p_third * p_third * p_third;
            const double discriminant = discriminant_q + discriminant_p;
            const double discriminant_tolerance =
                32.0 * std::numeric_limits<double>::epsilon() *
                (std::abs(discriminant_q) + std::abs(discriminant_p));

            if (std::abs(discriminant) <= discriminant_tolerance) {
                const double cbrt_q_half = std::cbrt(q_half);
                const double t1 = -2.0 * cbrt_q_half * root_scale;
                const double t2 = cbrt_q_half * root_scale;

                auto x1 = SymbolicExpr::number(t1 - shift_val);
                auto x2 = SymbolicExpr::number(t2 - shift_val);
                roots = { x1, x2, x2 };
            } else if (discriminant > 0.0) {
                const double sqrt_discriminant = std::sqrt(discriminant);
                const double u = std::cbrt(-q_half + sqrt_discriminant);
                const double v = std::cbrt(-q_half - sqrt_discriminant);
                const double real_sum = (u + v) * root_scale;
                const double x1_value = real_sum - shift_val;
                const double real_part = -real_sum / 2.0 - shift_val;
                const double imaginary_part =
                    std::sqrt(3.0) * (u - v) * root_scale / 2.0;

                auto x1 = SymbolicExpr::number(x1_value);
                auto i_unit = SymbolicExpr::sqrt(num(-1));
                auto x2 = SymbolicExpr::add(
                    SymbolicExpr::number(real_part),
                    SymbolicExpr::multiply(
                        SymbolicExpr::number(imaginary_part), i_unit))->simplify();
                auto x3 = sub(
                    SymbolicExpr::number(real_part),
                    SymbolicExpr::multiply(
                        SymbolicExpr::number(imaginary_part), i_unit))->simplify();
                roots = { x1, x2, x3 };
            } else {
                const double radius =
                    std::sqrt(-(p_normalized * p_normalized * p_normalized) / 27.0);
                double cosine_argument = -q_normalized / (2.0 * radius);
                cosine_argument = std::clamp(cosine_argument, -1.0, 1.0);
                const double theta = std::acos(cosine_argument);
                const double amplitude = 2.0 * std::cbrt(radius) * root_scale;

                for (int k = 0; k < 3; ++k) {
                    const double angle =
                        (theta + 2.0 * k * LMMC_CONST_PI) / 3.0;
                    roots.push_back(SymbolicExpr::number(
                        amplitude * std::cos(angle) - shift_val));
                }
            }
        }
    } else {
        auto a2 = SymbolicExpr::power(a, num(2));
        auto a3 = SymbolicExpr::power(a, num(3));
        auto b2 = SymbolicExpr::power(b, num(2));
        auto b3 = SymbolicExpr::power(b, num(3));
        auto three_ac = SymbolicExpr::multiply(
            num(3), SymbolicExpr::multiply(a, c));
        auto p_num = sub(three_ac, b2);
        auto p_den = SymbolicExpr::multiply(num(3), a2);
        auto p = SymbolicExpr::divide(p_num, p_den)->simplify();
        auto two_b3 = SymbolicExpr::multiply(num(2), b3);
        auto nine_abc = SymbolicExpr::multiply(
            num(9), SymbolicExpr::multiply(
                a, SymbolicExpr::multiply(b, c)));
        auto twentyseven_a2d = SymbolicExpr::multiply(
            num(27), SymbolicExpr::multiply(a2, d));
        auto q_num = SymbolicExpr::add(
            sub(two_b3, nine_abc), twentyseven_a2d);
        auto q_den = SymbolicExpr::multiply(num(27), a3);
        auto q = SymbolicExpr::divide(q_num, q_den)->simplify();
        auto shift = SymbolicExpr::divide(
            b, SymbolicExpr::multiply(num(3), a))->simplify();


        auto q_half = SymbolicExpr::divide(q, num(2));
        auto p_third = SymbolicExpr::divide(p, num(3));
        auto D_expr = SymbolicExpr::add(
            SymbolicExpr::power(q_half, num(2)),
            SymbolicExpr::power(p_third, num(3)))->simplify();

        auto neg_q_half = negate(q_half)->simplify();
        auto sqrt_D = SymbolicExpr::sqrt(D_expr);

        auto u = cbrt_expr(SymbolicExpr::add(neg_q_half, sqrt_D))->simplify();
        auto v = cbrt_expr(sub(neg_q_half, sqrt_D))->simplify();

        auto t1 = SymbolicExpr::add(u, v);
        auto x1 = sub(t1, shift)->simplify();

        auto neg_half_sum = SymbolicExpr::divide(negate(SymbolicExpr::add(u, v)), num(2));
        auto diff_uv = sub(u, v);
        auto sqrt3_half = SymbolicExpr::divide(SymbolicExpr::sqrt(num(3)), num(2));
        auto imag_coeff = SymbolicExpr::multiply(sqrt3_half, diff_uv);

        auto i_unit = SymbolicExpr::sqrt(num(-1));
        auto x2 = sub(SymbolicExpr::add(neg_half_sum, SymbolicExpr::multiply(i_unit, imag_coeff)), shift)->simplify();
        auto x3 = sub(sub(neg_half_sum, SymbolicExpr::multiply(i_unit, imag_coeff)), shift)->simplify();

        roots = { x1, x2, x3 };
    }

    return roots;
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_biquadratic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::string&) {
    auto a_simplified = a->simplify();
    if (a_simplified->get_number_value_is_zero()) {
        auto b_simplified = b->simplify();
        if (!b_simplified->get_number_value_is_zero()) {
            return solve_quadratic_internal(b, num(0), c);
        }

        auto c_simplified = c->simplify();
        if (c_simplified->get_number_value_is_zero()) {
            throw std::invalid_argument(
                "solve_biquadratic: identically zero equation has no finite root list");
        }
        if (is_purely_numeric(c_simplified)) {
            return {};
        }
        throw std::invalid_argument(
            "solve_biquadratic: equation is independent of the requested variable");
    }


    bool all_numeric = is_purely_numeric(a) && is_purely_numeric(b) && is_purely_numeric(c);

    auto a_checked = finite_numeric_value(a);
    auto b_checked = finite_numeric_value(b);
    auto c_checked = finite_numeric_value(c);

    if (all_numeric && a_checked && b_checked && c_checked && *a_checked != 0.0) {
        double av = *a_checked;
        double bv = *b_checked;
        double cv = *c_checked;
        double scale = std::max({std::abs(av), std::abs(bv), std::abs(cv)});
        double an = av / scale;
        double bn = bv / scale;
        double cn = cv / scale;
        bool normalization_lost_coefficient =
            (av != 0.0 && an == 0.0) ||
            (bv != 0.0 && bn == 0.0) ||
            (cv != 0.0 && cn == 0.0);

        if (!normalization_lost_coefficient) {
            double disc = std::fma(-4.0 * an, cn, bn * bn);
            if (std::isfinite(disc) && disc >= 0.0) {
                double sqrt_disc = std::sqrt(disc);
                double u1;
                double u2;
                if (sqrt_disc == 0.0) {
                    u1 = -bn / (2.0 * an);
                    u2 = u1;
                } else {
                    double q = -0.5 * (bn + std::copysign(sqrt_disc, bn));
                    u1 = q / an;
                    u2 = cn / q;
                }

                bool representable =
                    std::isfinite(u1) && std::isfinite(u2) &&
                    (cv == 0.0 || (u1 != 0.0 && u2 != 0.0));
                if (representable) {
                    std::vector<std::shared_ptr<SymbolicExpr>> results;
                    double u_vals[2] = {u1, u2};
                    for (int i = 0; i < 2; ++i) {
                        if (u_vals[i] >= 0.0) {
                            double x_pos = std::sqrt(u_vals[i]);
                            results.push_back(SymbolicExpr::number(x_pos));
                            results.push_back(SymbolicExpr::number(-x_pos));
                        } else {
                            auto u_expr = SymbolicExpr::number(u_vals[i]);
                            results.push_back(SymbolicExpr::sqrt(u_expr)->simplify());
                            results.push_back(negate(SymbolicExpr::sqrt(u_expr))->simplify());
                        }
                    }
                    return results;
                }
            }
        }
    }

    auto u_roots = solve_quadratic_internal(a, b, c);

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (const auto& u : u_roots) {
        auto pos_root = SymbolicExpr::sqrt(u)->simplify();
        auto neg_root = negate(SymbolicExpr::sqrt(u))->simplify();
        results.push_back(pos_root);
        results.push_back(neg_root);
    }

    return results;
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_quartic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e,
    const std::string& var) {

    auto a_simp = a->simplify();
    if (a_simp->get_number_value_is_zero()) {

        return solve_cubic(b, c, d, e, var);
    }

    {
        Polynomial<SymbolicPolyCoeff> symbolic_poly(
            {SymbolicPolyCoeff(e), SymbolicPolyCoeff(d), SymbolicPolyCoeff(c),
             SymbolicPolyCoeff(b), SymbolicPolyCoeff(a)},
            var);
        Polynomial<Rational> rational_poly;
        if (convert_to_rational_poly(symbolic_poly, rational_poly) &&
            rational_poly.degree() == 4) {
            auto rational_roots = find_rational_roots(rational_poly);
            if (rational_roots.size() == 4) {
                std::vector<std::shared_ptr<SymbolicExpr>> proven_roots;
                proven_roots.reserve(rational_roots.size());
                for (const auto& root : rational_roots) {
                    proven_roots.push_back(SymbolicExpr::number(root));
                }
                return proven_roots;
            }
        }
    }

    auto b_simp = b->simplify();
    auto d_simp = d->simplify();
    if (b_simp->get_number_value_is_zero() && d_simp->get_number_value_is_zero()) {
        return solve_biquadratic(a, c, e, var);
    }

    bool all_numeric = is_purely_numeric(a) && is_purely_numeric(b) &&
                       is_purely_numeric(c) && is_purely_numeric(d) && is_purely_numeric(e);

    auto a_checked = finite_numeric_value(a);
    auto b_checked = finite_numeric_value(b);
    auto c_checked = finite_numeric_value(c);
    auto d_checked = finite_numeric_value(d);
    auto e_checked = finite_numeric_value(e);

    if (all_numeric && a_checked && b_checked && c_checked && d_checked && e_checked &&
        *a_checked != 0.0) {

        double av = *a_checked;
        double bv = *b_checked;
        double cv = *c_checked;
        double dv = *d_checked;
        double ev = *e_checked;

        double shift_val = bv / (4.0 * av);
        double p_val = (8.0*av*cv - 3.0*bv*bv) / (8.0*av*av);
        double q_val = (bv*bv*bv - 4.0*av*bv*cv + 8.0*av*av*dv) / (8.0*av*av*av);
        double r_val = (-3.0*bv*bv*bv*bv + 256.0*av*av*av*ev - 64.0*av*av*bv*dv + 16.0*av*bv*bv*cv) / (256.0*av*av*av*av);

        if (!std::isfinite(shift_val) || !std::isfinite(p_val) ||
            !std::isfinite(q_val) || !std::isfinite(r_val)) {
            all_numeric = false;
        }

        if (all_numeric) {

        if (q_val == 0.0) {
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

        if (q_normalized != 0.0) {
            auto cubic_roots = solve_cubic(
                SymbolicExpr::number(8.0),
                SymbolicExpr::number(8.0 * p_normalized),
                SymbolicExpr::number(
                    2.0 * p_normalized * p_normalized -
                    8.0 * r_normalized),
                SymbolicExpr::number(-(q_normalized * q_normalized)),
                var);

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

            if (found_m && m_val != 0.0) {
                if (m_val < 0.0) {
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

                const double s_val = std::sqrt(2.0 * m_val);
                const double quad1_c_val =
                    m_val + p_normalized / 2.0 -
                    q_normalized / (2.0 * s_val);
                const double quad2_c_val =
                    m_val + p_normalized / 2.0 +
                    q_normalized / (2.0 * s_val);
                std::vector<std::shared_ptr<SymbolicExpr>> results;

                const double disc1 =
                    s_val * s_val - 4.0 * quad1_c_val;
                if (disc1 >= 0.0) {
                    const double sqrt_disc1 = std::sqrt(disc1);
                    const double y1 = (-s_val + sqrt_disc1) / 2.0;
                    const double y2 = (-s_val - sqrt_disc1) / 2.0;
                    results.push_back(SymbolicExpr::number(
                        y1 * root_scale - shift_val));
                    results.push_back(SymbolicExpr::number(
                        y2 * root_scale - shift_val));
                } else {
                    const double real_part =
                        -s_val * root_scale / 2.0 - shift_val;
                    const double imag_part =
                        std::sqrt(-disc1) * root_scale / 2.0;
                    auto i_unit = SymbolicExpr::sqrt(num(-1));
                    results.push_back(SymbolicExpr::add(
                        SymbolicExpr::number(real_part),
                        SymbolicExpr::multiply(
                            SymbolicExpr::number(imag_part),
                            i_unit))->simplify());
                    results.push_back(sub(
                        SymbolicExpr::number(real_part),
                        SymbolicExpr::multiply(
                            SymbolicExpr::number(imag_part),
                            i_unit))->simplify());
                }

                const double disc2 =
                    s_val * s_val - 4.0 * quad2_c_val;
                if (disc2 >= 0.0) {
                    const double sqrt_disc2 = std::sqrt(disc2);
                    const double y3 = (s_val + sqrt_disc2) / 2.0;
                    const double y4 = (s_val - sqrt_disc2) / 2.0;
                    results.push_back(SymbolicExpr::number(
                        y3 * root_scale - shift_val));
                    results.push_back(SymbolicExpr::number(
                        y4 * root_scale - shift_val));
                } else {
                    const double real_part =
                        s_val * root_scale / 2.0 - shift_val;
                    const double imag_part =
                        std::sqrt(-disc2) * root_scale / 2.0;
                    auto i_unit = SymbolicExpr::sqrt(num(-1));
                    results.push_back(SymbolicExpr::add(
                        SymbolicExpr::number(real_part),
                        SymbolicExpr::multiply(
                            SymbolicExpr::number(imag_part),
                            i_unit))->simplify());
                    results.push_back(sub(
                        SymbolicExpr::number(real_part),
                        SymbolicExpr::multiply(
                            SymbolicExpr::number(imag_part),
                            i_unit))->simplify());
                }
                return results;
            }
        }
    }
    }

    auto a2 = SymbolicExpr::power(a, num(2));
    auto a3 = SymbolicExpr::power(a, num(3));
    auto a4 = SymbolicExpr::power(a, num(4));
    auto b2 = SymbolicExpr::power(b, num(2));
    auto b3 = SymbolicExpr::power(b, num(3));
    auto b4 = SymbolicExpr::power(b, num(4));

    auto eight_ac = SymbolicExpr::multiply(num(8), SymbolicExpr::multiply(a, c));
    auto three_b2 = SymbolicExpr::multiply(num(3), b2);
    auto p = SymbolicExpr::divide(sub(eight_ac, three_b2), SymbolicExpr::multiply(num(8), a2))->simplify();

    auto four_abc = SymbolicExpr::multiply(num(4), SymbolicExpr::multiply(a, SymbolicExpr::multiply(b, c)));
    auto eight_a2d = SymbolicExpr::multiply(num(8), SymbolicExpr::multiply(a2, d));
    auto q_num = SymbolicExpr::add(sub(b3, four_abc), eight_a2d);
    auto q = SymbolicExpr::divide(q_num, SymbolicExpr::multiply(num(8), a3))->simplify();

    auto neg3_b4 = SymbolicExpr::multiply(num(-3), b4);
    auto t256_a3e = SymbolicExpr::multiply(num(256), SymbolicExpr::multiply(a3, e));
    auto neg64_a2bd = SymbolicExpr::multiply(num(-64), SymbolicExpr::multiply(a2, SymbolicExpr::multiply(b, d)));
    auto t16_ab2c = SymbolicExpr::multiply(num(16), SymbolicExpr::multiply(a, SymbolicExpr::multiply(b2, c)));
    auto r_num = SymbolicExpr::add(SymbolicExpr::add(SymbolicExpr::add(neg3_b4, t256_a3e), neg64_a2bd), t16_ab2c);
    auto r = SymbolicExpr::divide(r_num, SymbolicExpr::multiply(num(256), a4))->simplify();

    auto shift = SymbolicExpr::divide(b, SymbolicExpr::multiply(num(4), a))->simplify();

    auto q_simp = q->simplify();
    if (q_simp->get_number_value_is_zero()) {

        auto u_roots = solve_quadratic_internal(num(1), p, r);
        std::vector<std::shared_ptr<SymbolicExpr>> results;
        for (const auto& u : u_roots) {

            auto pos_y = SymbolicExpr::sqrt(u)->simplify();
            auto neg_y = negate(SymbolicExpr::sqrt(u))->simplify();
            results.push_back(sub(pos_y, shift)->simplify());
            results.push_back(sub(neg_y, shift)->simplify());
        }
        return results;
    }

    auto p2 = SymbolicExpr::power(p, num(2));
    auto q2 = SymbolicExpr::power(q, num(2));
    auto eight_p = SymbolicExpr::multiply(num(8), p);
    auto two_p2_minus_8r = sub(SymbolicExpr::multiply(num(2), p2), SymbolicExpr::multiply(num(8), r))->simplify();
    auto neg_q2 = negate(q2)->simplify();

    auto cubic_roots = solve_cubic(num(8), eight_p, two_p2_minus_8r, neg_q2, var);

    std::shared_ptr<SymbolicExpr> m = nullptr;
    if (!cubic_roots.empty()) {
        m = cubic_roots[0];
    }

    if (!m) return {};

    auto two_m = SymbolicExpr::multiply(num(2), m);
    auto s = SymbolicExpr::sqrt(two_m)->simplify();

    auto p_half = SymbolicExpr::divide(p, num(2))->simplify();

    auto m_plus_p_half = SymbolicExpr::add(m, p_half)->simplify();

    auto q_over_2s = SymbolicExpr::divide(q, SymbolicExpr::multiply(num(2), s))->simplify();

    auto quad1_c = sub(m_plus_p_half, q_over_2s)->simplify();

    auto neg_s = negate(s)->simplify();
    auto quad2_c = SymbolicExpr::add(m_plus_p_half, q_over_2s)->simplify();

    auto y_roots1 = solve_quadratic_internal(num(1), s, quad1_c);
    auto y_roots2 = solve_quadratic_internal(num(1), neg_s, quad2_c);

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (const auto& y : y_roots1) {
        results.push_back(sub(y, shift)->simplify());
    }
    for (const auto& y : y_roots2) {
        results.push_back(sub(y, shift)->simplify());
    }

    return results;
}

static std::vector<BigInt> positive_divisors(const BigInt& n) {
    std::vector<BigInt> divs;
    if (n.is_zero()) {

        return divs;
    }
    BigInt abs_n = n.Abs();

    BigInt i(1);

    while (i * i <= abs_n) {
        BigInt rem = abs_n % i;
        if (rem.is_zero()) {
            divs.push_back(i);
            BigInt counterpart = abs_n / i;
            if (counterpart != i) {
                divs.push_back(counterpart);
            }
        }
        i = i + BigInt(1);
    }

    std::sort(divs.begin(), divs.end());
    return divs;
}

std::vector<std::pair<Polynomial<Rational>, int>> square_free_factorization(
    const Polynomial<Rational>& poly) {

    std::vector<std::pair<Polynomial<Rational>, int>> factors;

    if (poly.is_zero() || poly.degree() <= 0) {
        return factors;
    }

    Polynomial<Rational> f = poly.make_monic();

    if (f.degree() == 1) {
        factors.push_back({f, 1});
        return factors;
    }

    Polynomial<Rational> f_prime = f.differentiate();
    Polynomial<Rational> g = Polynomial<Rational>::gcd(f, f_prime);

    if (g.degree() <= 0) {
        factors.push_back({f, 1});
        return factors;
    }

    auto [w, w_rem] = f.div_mod(g);

    int multiplicity = 1;

    while (w.degree() > 0) {

        Polynomial<Rational> y = Polynomial<Rational>::gcd(g, w);

        auto [z, z_rem] = w.div_mod(y);

        if (z.degree() > 0) {

            factors.push_back({z.make_monic(), multiplicity});
        }

        auto [g_next, g_rem] = g.div_mod(y);
        g = g_next;
        w = y;
        multiplicity++;
    }

    if (g.degree() > 0) {
        factors.push_back({g.make_monic(), multiplicity});
    }

    return factors;
}

static std::vector<std::shared_ptr<SymbolicExpr>> solve_closed_form_poly(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var)
{
    int deg = poly.degree();
    if (deg <= 0) return {};

    auto get_coeff = [&](int d) -> std::shared_ptr<SymbolicExpr> {
        if (d < 0 || d > deg) return SymbolicExpr::number(0);
        return poly.coeffs[d].val ? poly.coeffs[d].val : SymbolicExpr::number(0);
    };

    if (deg == 1) {
        auto a = get_coeff(1);
        auto b = get_coeff(0);
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        return { SymbolicExpr::divide(neg_b, a)->simplify() };
    } else if (deg == 2) {
        return solve_quadratic_internal(get_coeff(2), get_coeff(1), get_coeff(0));
    } else if (deg == 3) {
        return solve_cubic(get_coeff(3), get_coeff(2), get_coeff(1), get_coeff(0), var);
    } else if (deg == 4) {
        return solve_quartic(get_coeff(4), get_coeff(3), get_coeff(2), get_coeff(1), get_coeff(0), var);
    }

    return {};
}

static bool convert_to_rational_poly(
    const Polynomial<SymbolicPolyCoeff>& sym_poly,
    Polynomial<Rational>& out_poly)
{
    out_poly = Polynomial<Rational>(sym_poly.variable_name);
    out_poly.coeffs.resize(sym_poly.coeffs.size(), Rational(0));

    for (size_t i = 0; i < sym_poly.coeffs.size(); ++i) {
        auto coeff_expr = sym_poly.coeffs[i].val;
        if (!coeff_expr) {
            out_poly.coeffs[i] = Rational(0);
            continue;
        }

        auto simplified = coeff_expr->simplify();
        auto num_node = std::dynamic_pointer_cast<const NumberNode>(
            LMCAS::detail::node(simplified));
        if (!num_node) {
            return false;
        }
        if (std::holds_alternative<Rational>(num_node->value())) {
            out_poly.coeffs[i] = std::get<Rational>(num_node->value());
        } else if (std::holds_alternative<BigInt>(num_node->value())) {
            out_poly.coeffs[i] =
                Rational(std::get<BigInt>(num_node->value()));
        } else {
            return false;
        }
    }

    out_poly.trim();
    return true;
}

static Polynomial<SymbolicPolyCoeff> rational_to_symbolic_poly(
    const Polynomial<Rational>& rat_poly)
{
    std::vector<SymbolicPolyCoeff> sym_coeffs;
    sym_coeffs.reserve(rat_poly.coeffs.size());
    for (const auto& c : rat_poly.coeffs) {
        sym_coeffs.push_back(SymbolicPolyCoeff(SymbolicExpr::number(c)));
    }
    return Polynomial<SymbolicPolyCoeff>(sym_coeffs, rat_poly.variable_name);
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_by_factoring(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var) {

    if (poly.is_zero() || poly.degree() < 1) return {};

    if (poly.degree() <= 2) {
        return solve_closed_form_poly(poly, var);
    }

    std::vector<std::shared_ptr<SymbolicExpr>> results;

    Polynomial<Rational> rat_poly;
    if (!convert_to_rational_poly(poly, rat_poly)) {
        if (poly.degree() <= 4) {
            return solve_closed_form_poly(poly, var);
        }
        return make_rootof_solutions(poly, var);
    }

    // Keep rational boundaries exact before cubic/quartic floating formulas.
    auto rational_roots = find_rational_roots(rat_poly);

    for (const auto& r : rational_roots) {
        results.push_back(SymbolicExpr::number(r));
    }

    Polynomial<Rational> quotient = rat_poly;
    for (const auto& r : rational_roots) {
        Polynomial<Rational> linear_factor({-r, Rational(1)}, var);
        auto [q, rem] = quotient.div_mod(linear_factor);
        quotient = q;
    }

    if (quotient.degree() <= 0) {
        return results;
    }

    if (quotient.degree() <= 4) {
        auto sym_quotient = rational_to_symbolic_poly(quotient);
        auto factor_roots = solve_closed_form_poly(sym_quotient, var);
        results.insert(results.end(), factor_roots.begin(), factor_roots.end());
        return results;
    }

    auto sqfree_factors = square_free_factorization(quotient);

    for (const auto& [factor, multiplicity] : sqfree_factors) {
        if (factor.degree() <= 0) continue;

        std::vector<std::shared_ptr<SymbolicExpr>> factor_solutions;

        if (factor.degree() <= 4) {

            auto sym_factor = rational_to_symbolic_poly(factor);
            factor_solutions = solve_closed_form_poly(sym_factor, var);
        } else {

            auto sym_factor = rational_to_symbolic_poly(factor);
            factor_solutions = make_rootof_solutions(sym_factor, var);
        }

        for (int m = 0; m < multiplicity; ++m) {
            results.insert(results.end(), factor_solutions.begin(), factor_solutions.end());
        }
    }

    return results;
}

std::vector<Rational> find_rational_roots(const Polynomial<Rational>& poly) {
    std::vector<Rational> roots;

    if (poly.is_zero() || poly.degree() < 1) {
        return roots;
    }

    // Monic normalization gives a scale-independent rational polynomial.
    // Clearing its denominators gives a primitive integer candidate polynomial.
    Polynomial<Rational> current = poly.make_monic();

    while (current.degree() >= 1 && current.coeffs[0] == Rational(0)) {
        roots.push_back(Rational(0));

        std::vector<Rational> new_coeffs(current.coeffs.begin() + 1, current.coeffs.end());
        current = Polynomial<Rational>(new_coeffs, current.variable_name);
    }

    while (current.degree() >= 1) {

        Rational a0 = current.coeffs[0];
        Rational an = current.lead_coeff();

        // A linear remainder is solved directly by exact rational division.
        if (current.degree() == 1) {
            roots.push_back(-a0 / an);
            break;
        }

        if (a0 == Rational(0)) {

            roots.push_back(Rational(0));
            std::vector<Rational> new_coeffs(current.coeffs.begin() + 1, current.coeffs.end());
            current = Polynomial<Rational>(new_coeffs, current.variable_name);
            continue;
        }

        BigInt lcm_den(1);
        for (const auto& c : current.coeffs) {
            BigInt d = c.get_denominator();
            lcm_den = BigInt::lcm(lcm_den, d);
        }

        BigInt int_a0 = (a0 * Rational(lcm_den)).get_numerator();
        BigInt int_an = (an * Rational(lcm_den)).get_numerator();

        std::vector<BigInt> p_divs = positive_divisors(int_a0);
        std::vector<BigInt> q_divs = positive_divisors(int_an);

        if (p_divs.empty() || q_divs.empty()) {
            break;
        }

        bool found_root = false;
        Rational found_r;

        for (const auto& p : p_divs) {
            for (const auto& q : q_divs) {

                Rational candidate_pos(p, q);
                Rational candidate_neg(-p, q);

                if (current.eval(candidate_pos) == Rational(0)) {
                    found_root = true;
                    found_r = candidate_pos;
                    break;
                }
                if (current.eval(candidate_neg) == Rational(0)) {
                    found_root = true;
                    found_r = candidate_neg;
                    break;
                }
            }
            if (found_root) break;
        }

        if (!found_root) {
            break;
        }

        Polynomial<Rational> linear_factor({-found_r, Rational(1)}, current.variable_name);

        while (current.degree() >= 1) {
            auto [quotient, remainder] = current.div_mod(linear_factor);

            if (!remainder.is_zero()) {
                break;
            }
            roots.push_back(found_r);
            current = quotient;

            if (current.degree() < 1 || current.eval(found_r) != Rational(0)) {
                break;
            }
        }
    }

    return roots;
}

}
