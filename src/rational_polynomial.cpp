#include "rational_polynomial.hpp"
#include <algorithm>
#include <optional>

namespace LMCAS {

static std::vector<BigInt> positive_divisors(const BigInt& n) {
    std::vector<BigInt> divs;
    if (n.is_zero()) {

        return divs;
    }
    BigInt abs_n = n.abs();

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

static std::optional<Rational> rational_root_candidate(
    const Polynomial<Rational>& current) {
    BigInt denominator_lcm(1);
    for (const auto& coefficient : current.coeffs) {
        denominator_lcm = BigInt::lcm(
            denominator_lcm, coefficient.get_denominator());
    }
    const Rational scale(denominator_lcm);
    const auto numerators = positive_divisors(
        (current.coeffs[0] * scale).get_numerator());
    const auto denominators = positive_divisors(
        (current.lead_coeff() * scale).get_numerator());
    for (const auto& numerator : numerators) {
        for (const auto& denominator : denominators) {
            Rational positive(numerator, denominator);
            if (current.eval(positive) == Rational(0)) return positive;
            Rational negative(-numerator, denominator);
            if (current.eval(negative) == Rational(0)) return negative;
        }
    }
    return std::nullopt;
}

static void deflate_rational_root(Polynomial<Rational>& current,
                                  const Rational& root,
                                  std::vector<Rational>& roots) {
    Polynomial<Rational> linear({-root, Rational(1)}, current.variable_name);
    while (current.degree() >= 1) {
        auto [quotient, remainder] = current.div_mod(linear);
        if (!remainder.is_zero()) break;
        roots.push_back(root);
        current = std::move(quotient);
        if (current.degree() < 1 || current.eval(root) != Rational(0)) break;
    }
}

std::vector<Rational> find_rational_roots(const Polynomial<Rational>& poly) {
    std::vector<Rational> roots;
    if (poly.is_zero() || poly.degree() < 1) return roots;
    Polynomial<Rational> current = poly.make_monic();
    while (current.degree() >= 1) {
        if (current.coeffs[0] == Rational(0)) {
            roots.push_back(Rational(0));
            std::vector<Rational> coefficients(
                current.coeffs.begin() + 1, current.coeffs.end());
            current = Polynomial<Rational>(coefficients, current.variable_name);
            continue;
        }
        if (current.degree() == 1) {
            roots.push_back(-current.coeffs[0] / current.lead_coeff());
            break;
        }
        auto candidate = rational_root_candidate(current);
        if (!candidate) break;
        deflate_rational_root(current, *candidate, roots);
    }
    return roots;
}

}
