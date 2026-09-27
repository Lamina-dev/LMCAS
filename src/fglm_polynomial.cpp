#include "fglm.hpp"
#include <algorithm>
#include <stdexcept>

namespace LMCAS {

namespace {

Monomial subtract_monomials(
    const Monomial& left, const Monomial& right, std::size_t num_vars) {
    Monomial result(num_vars, 0);
    for (std::size_t i = 0; i < num_vars; ++i) {
        const int left_value = i < left.size() ? left[i] : 0;
        const int right_value = i < right.size() ? right[i] : 0;
        result[i] = left_value - right_value;
    }
    return result;
}

Monomial add_monomials(
    const Monomial& left, const Monomial& right, std::size_t num_vars) {
    Monomial result(num_vars, 0);
    for (std::size_t i = 0; i < num_vars; ++i) {
        const int left_value = i < left.size() ? left[i] : 0;
        const int right_value = i < right.size() ? right[i] : 0;
        result[i] = left_value + right_value;
    }
    return result;
}

void subtract_leading_multiple(
    FGLMPoly& polynomial, const FGLMPoly& divisor,
    const Monomial& divisor_leading, const Monomial& leading,
    const MonomialOrder& order) {
    const Monomial quotient = subtract_monomials(
        leading, divisor_leading, polynomial.num_vars);
    const Rational coefficient =
        polynomial.lead_coeff() / divisor.lead_coeff();
    for (const auto& [monomial, value] : divisor.terms) {
        polynomial.add_term(
            add_monomials(quotient, monomial, polynomial.num_vars),
            -(coefficient * value));
    }
    polynomial.sort_terms(order);
    polynomial.normalize();
}

bool reduce_leading_term(
    FGLMPoly& polynomial, const std::vector<FGLMPoly>& basis,
    const Monomial& leading, const MonomialOrder& order) {
    for (const auto& divisor : basis) {
        if (divisor.is_zero()) {
            continue;
        }
        const Monomial divisor_leading = divisor.lead_monomial();
        if (!divides_monomial(divisor_leading, leading)) {
            continue;
        }
        subtract_leading_multiple(
            polynomial, divisor, divisor_leading, leading, order);
        return true;
    }
    return false;
}

}


    void FGLMPoly::sort_terms(const MonomialOrder& order) {
        std::sort(terms.begin(), terms.end(),
            [&order](const std::pair<Monomial, Rational>& a,
                     const std::pair<Monomial, Rational>& b) {
                return order(a.first, b.first);
            });
    }

    void FGLMPoly::normalize() {

        std::vector<std::pair<Monomial, Rational>> cleaned;
        for (auto& [m, c] : terms) {
            if (c.is_zero()) continue;
            if (!cleaned.empty() && cleaned.back().first == m) {
                cleaned.back().second = cleaned.back().second + c;
                if (cleaned.back().second.is_zero()) {
                    cleaned.pop_back();
                }
            } else {
                cleaned.emplace_back(m, c);
            }
        }
        terms = std::move(cleaned);
    }

    FGLMPoly FGLMPoly::reduce(const std::vector<FGLMPoly>& basis,
                    const MonomialOrder& order) const {
        FGLMPoly remainder(num_vars);
        FGLMPoly polynomial = *this;
        polynomial.sort_terms(order);
        polynomial.normalize();

        while (!polynomial.is_zero()) {
            const Monomial previous_leading = polynomial.lead_monomial();
            if (!reduce_leading_term(
                    polynomial, basis, previous_leading, order)) {
                remainder.add_term(
                    previous_leading, polynomial.lead_coeff());
                polynomial.terms.erase(polynomial.terms.begin());
                polynomial.normalize();
            }
            if (!polynomial.is_zero() &&
                !order(previous_leading, polynomial.lead_monomial())) {
                throw std::runtime_error(
                    "FGLM polynomial reduction did not decrease the leading monomial");
            }
        }
        remainder.sort_terms(order);
        remainder.normalize();
        return remainder;
    }

    FGLMPoly FGLMPoly::from_monomial(const Monomial& m, size_t num_vars) {
        FGLMPoly p(num_vars);
        p.add_term(m, Rational(1));
        return p;
    }

FGLMPoly normal_form(const FGLMPoly& f,
                            const std::vector<FGLMPoly>& basis,
                            const MonomialOrder& order) {
    return f.reduce(basis, order);
}

}
