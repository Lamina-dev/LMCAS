#include "internal/polynomial_solver_support.hpp"
#include "root_of_identity.hpp"
#include "rational_polynomial.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace LMCAS {
using namespace polynomial_solver_detail;

static std::shared_ptr<SymbolicExpr> cbrt_expr(const std::shared_ptr<SymbolicExpr>& x) {
    return SymbolicExpr::power(x, SymbolicExpr::number(Rational(1, 3)));
}

static std::vector<std::shared_ptr<SymbolicExpr>> solve_degenerate_cubic(
    const std::shared_ptr<SymbolicExpr>& b, const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d) {
    if (!b->simplify()->get_number_value_is_zero()) {
        return solve_quadratic_internal(b, c, d);
    }
    if (!c->simplify()->get_number_value_is_zero()) {
        return {SymbolicExpr::divide(negate(d), c)->simplify()};
    }
    auto constant = d->simplify();
    if (constant->get_number_value_is_zero()) {
        throw std::invalid_argument(
            "solve_cubic: identically zero equation has no finite root list");
    }
    if (is_purely_numeric(constant)) return {};
    throw std::invalid_argument(
        "solve_cubic: equation is independent of the requested variable");
}

static std::vector<std::shared_ptr<SymbolicExpr>> symbolic_cubic_roots(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d) {
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

    return {x1, x2, x3};
}

bool polynomial_solver_detail::bounded_rational_root_search(
    const Polynomial<Rational>& polynomial) {
    if (polynomial.degree() <= 1) return true;
    if (polynomial.lead_coeff() != Rational(1)) {
        return bounded_rational_root_search(polynomial.make_monic());
    }
    BigInt denominator_lcm(1);
    for (const auto& coefficient : polynomial.coeffs) {
        denominator_lcm = BigInt::lcm(
            denominator_lcm, coefficient.get_denominator());
        if (denominator_lcm.bit_length() > 32) return false;
    }
    /**
     * @brief find_rational_roots 在枚举因子前移除零根。
     * 大小限制针对本原整数多项式的实际首尾系数；
     * 各分母虽小，其最小公倍数仍可能很大。
     */
    for (const auto& coefficient : polynomial.coeffs) {
        if (coefficient != Rational(0)) {
            return (coefficient * Rational(denominator_lcm))
                       .get_numerator().bit_length() <= 32;
        }
    }
    return true;
}

static std::vector<std::shared_ptr<SymbolicExpr>> exact_factor_roots(
    Polynomial<Rational> polynomial) {
    std::vector<std::shared_ptr<SymbolicExpr>> roots;
    if (polynomial.degree() == 3 && polynomial.coeffs[0] == Rational(0)) {
        roots.push_back(SymbolicExpr::number(0));
        polynomial.coeffs.erase(polynomial.coeffs.begin());
    }
    if (polynomial.degree() == 3) {
        if (bounded_rational_root_search(polynomial)) {
            for (const auto& root : find_rational_roots(polynomial)) {
                Polynomial<Rational> linear(
                    {-root, Rational(1)}, polynomial.variable_name);
                auto [quotient, remainder] = polynomial.div_mod(linear);
                if (!remainder.is_zero()) {
                    throw std::logic_error("solve_cubic: inexact rational root division");
                }
                roots.push_back(SymbolicExpr::number(root));
                polynomial = std::move(quotient);
            }
        }
    }
    if (polynomial.degree() == 1) {
        roots.push_back(SymbolicExpr::number(
            -polynomial.coeffs[0] / polynomial.coeffs[1]));
    } else if (polynomial.degree() == 2) {
        const auto& a = polynomial.coeffs[2];
        const auto& b = polynomial.coeffs[1];
        const auto& c = polynomial.coeffs[0];
        auto radical = SymbolicExpr::sqrt(
            SymbolicExpr::number(b * b - Rational(4) * a * c));
        auto negative_b = SymbolicExpr::number(-b);
        auto denominator = SymbolicExpr::number(Rational(2) * a);
        roots.push_back(SymbolicExpr::divide(
            SymbolicExpr::add(negative_b, radical), denominator)->simplify());
        roots.push_back(SymbolicExpr::divide(
            sub(negative_b, radical), denominator)->simplify());
    } else if (polynomial.degree() == 3) {
        auto expression = poly_to_symbolic(polynomial);
        for (int index = 0; index < polynomial.degree(); ++index) {
            roots.push_back(SymbolicExpr::root_of(
                expression, polynomial.variable_name, index));
        }
    }
    return roots;
}

static std::vector<std::shared_ptr<SymbolicExpr>> exact_cubic_roots(
    const Polynomial<Rational>& polynomial) {
    if (polynomial.is_zero()) {
        throw std::invalid_argument(
            "solve_cubic: identically zero equation has no finite root list");
    }
    std::vector<std::shared_ptr<SymbolicExpr>> roots;
    for (const auto& [factor, multiplicity] :
         square_free_factorization(polynomial)) {
        auto factor_roots = exact_factor_roots(factor);
        for (int repeat = 0; repeat < multiplicity; ++repeat) {
            roots.insert(roots.end(), factor_roots.begin(), factor_roots.end());
        }
    }
    return roots;
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_cubic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::string& variable) {
    Polynomial<SymbolicPolyCoeff> symbolic(
        {SymbolicPolyCoeff(d), SymbolicPolyCoeff(c),
         SymbolicPolyCoeff(b), SymbolicPolyCoeff(a)}, variable);
    Polynomial<Rational> exact;
    if (convert_to_rational_poly(symbolic, exact)) {
        return exact_cubic_roots(exact);
    }
    if (a->simplify()->get_number_value_is_zero()) {
        return solve_degenerate_cubic(b, c, d);
    }
    auto depression = numeric_cubic_depression(a, b, c, d);
    if (depression) return numeric_cubic_roots(*depression);
    return symbolic_cubic_roots(a, b, c, d);
}

}
