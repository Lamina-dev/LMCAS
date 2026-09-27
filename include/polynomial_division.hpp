#ifndef LMCAS_POLYNOMIAL_DIVISION_HPP
#define LMCAS_POLYNOMIAL_DIVISION_HPP

#include "polynomial.hpp"

namespace LMCAS {

template <typename CoeffType>
inline CoeffType Polynomial<CoeffType>::content() const {
    if (is_zero()) return CoeffType(0);
    CoeffType g = coeffs[0];

    if constexpr (std::is_same_v<CoeffType, Rational>) {
         return CoeffType(1);
    } else {

        for (size_t i = 1; i < coeffs.size(); ++i) {
            g = gcd_coeff(g, coeffs[i]);
            if (g == CoeffType(1)) break;
        }
        return g;
    }
}

template <typename CoeffType>
inline Polynomial<CoeffType> Polynomial<CoeffType>::primitive_part() const {
    if (is_zero()) return *this;
    CoeffType c = content();
    if (c == CoeffType(0)) return *this;

    Polynomial res = *this;
    if (c != CoeffType(1)) {
        for (auto& val : res.coeffs) {
            val = val / c;
        }
    }

    if (res.lead_coeff() < CoeffType(0)) {
        for (auto& val : res.coeffs) {
            val = CoeffType(0) - val;
        }
    }
    return res;
}

template <typename CoeffType>
inline std::pair<Polynomial<CoeffType>, Polynomial<CoeffType>>
Polynomial<CoeffType>::pseudo_div_mod(const Polynomial<CoeffType>& other) const {
    if (other.is_zero()) throw std::runtime_error("Division by zero polynomial");
    require_compatible_variables(other, "Polynomial::pseudo_div_mod");
    if (is_zero()) {
        return {Polynomial(other.variable_name), Polynomial(other.variable_name)};
    }

    int deg_rem = degree();
    int deg_div = other.degree();

    if (deg_rem < deg_div) {
         return {Polynomial(variable_name), *this};
    }

    Polynomial remainder = *this;
    Polynomial quotient(variable_name);
    CoeffType lc_div = other.lead_coeff();

    int delta = deg_rem - deg_div;
    quotient.coeffs.resize(delta + 1, CoeffType(0));

    while (deg_rem >= deg_div && !remainder.is_zero()) {
         int current_diff = deg_rem - deg_div;
         CoeffType lc_rem = remainder.lead_coeff();

         for (auto& c : quotient.coeffs) c = c * lc_div;

         quotient.coeffs[current_diff] = quotient.coeffs[current_diff] + lc_rem;

         Polynomial term(variable_name);
         term.coeffs.resize(current_diff + 1, CoeffType(0));
         term.coeffs[current_diff] = lc_rem;

         for(auto& c : remainder.coeffs) c = c * lc_div;

         for (size_t i = 0; i < other.coeffs.size(); ++i) {
             size_t target_idx = i + current_diff;

             if (target_idx < remainder.coeffs.size()) {
                remainder.coeffs[target_idx] = remainder.coeffs[target_idx] - other.coeffs[i] * lc_rem;
             }

         }

         remainder.trim();
         deg_rem = remainder.degree();
    }

    return {quotient, remainder};
}

template <typename CoeffType>
inline Polynomial<CoeffType>
Polynomial<CoeffType>::pseudo_div_mod_rem(const Polynomial<CoeffType>& other) const {

    return pseudo_div_mod(other).second;
}

template <typename CoeffType>
inline Polynomial<CoeffType> Polynomial<CoeffType>::gcd(
    Polynomial<CoeffType> a, Polynomial<CoeffType> b) {
    a.require_compatible_variables(b, "Polynomial::gcd");
    if (b.is_zero()) return a;
    if (a.is_zero()) return b;

    if constexpr (std::is_same_v<CoeffType, Rational>) {

        while (!b.is_zero()) {
            auto [q, r] = a.div_mod(b);
            a = b;
            b = r;
        }

        return a.make_monic();
    } else {

        CoeffType cA = a.content();
        CoeffType cB = b.content();
        CoeffType c  = gcd_coeff_impl(cA, cB);

        a = a.primitive_part();
        b = b.primitive_part();

        while (!b.is_zero()) {

            Polynomial r = a.pseudo_div_mod_rem(b);

            if (r.is_zero()) {
                a = b;
                b = r;
            } else {

                a = b;
                b = r.primitive_part();
            }
        }

        if (c == CoeffType(1)) return a;

        for (auto& val : a.coeffs) val = val * c;
        return a;
    }
}

template <typename CoeffType>
inline Polynomial<CoeffType> Polynomial<CoeffType>::square_free_part() const {

    if (degree() <= 0) return *this;

    Polynomial deriv = differentiate();
    Polynomial g = gcd(*this, deriv);

    auto [q, r] = div_mod(g);

    if constexpr (std::is_same_v<CoeffType, Rational>) {
        return q.make_monic();
    }
    return q;
}

} // namespace LMCAS

#endif
