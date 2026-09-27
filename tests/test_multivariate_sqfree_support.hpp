#pragma once

#include "test_multivariate_support.hpp"

using namespace LMCAS;

inline MultiPoly test_formal_derivative(const MultiPoly& poly, const std::string& main_var)
{
    if (poly.is_zero()) { return poly; }

    const auto& vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    if (var_idx < 0) {
        return MultiPoly(Rational(0), vars);
    }

    std::vector<MultiPoly::Term> result_terms;
    for (const auto& term : poly.terms()) {
        const Monomial& mono = term.first;
        int exp = (static_cast<size_t>(var_idx) < mono.size()) ? mono[var_idx] : 0;
        if (exp == 0) { continue; }

        Rational new_coeff = term.second * Rational(exp);
        Monomial new_mono = mono;
        new_mono[var_idx] = exp - 1;
        result_terms.emplace_back(std::move(new_mono), std::move(new_coeff));
    }

    if (result_terms.empty()) { return MultiPoly(Rational(0), vars); }
    return MultiPoly(std::move(result_terms), vars);
}

inline bool is_square_free(const MultiPoly& poly, const std::string& main_var)
{
    if (poly.is_zero() || poly.is_constant()) { return true; }

    MultiPoly deriv = test_formal_derivative(poly, main_var);
    if (deriv.is_zero()) { return true; }

    MultiPoly g = multivariate_gcd(poly, deriv);
    return g.is_constant();
}

inline MultiPoly poly_pow(const MultiPoly& base, int exp)
{
    if (exp == 0) {
        return MultiPoly(Rational(1), base.variables());
    }
    MultiPoly result = base;
    for (int i = 1; i < exp; ++i) {
        result = result * base;
    }
    return result;
}

/**
 * @brief 验证无平方分解的重构关系与分量性质。
 * @return f₁*f₂²*f₃³*... 与 f 至多差常数倍，且各 fᵢ 无平方时返回 true。
 */
inline bool verify_sqfree_decomp(const MultiPoly& original,
                                 const SquareFreeDecomp& decomp,
                                 const std::string& main_var)
{
    const auto& comps = decomp.components;
    if (comps.empty()) { return original.is_zero(); }

    MultiPoly product = poly_pow(comps[0], 1);
    for (size_t i = 1; i < comps.size(); ++i) {
        product = product * poly_pow(comps[i], static_cast<int>(i + 1));
    }

    if (original.is_zero() && product.is_zero()) { return true; }
    if (original.is_zero() || product.is_zero()) { return false; }

    MultiPoly orig_prim = original.make_primitive();
    MultiPoly prod_prim = product.make_primitive();

    bool product_matches = (orig_prim == prod_prim) ||
                           (orig_prim == (prod_prim * Rational(-1)));

    if (!product_matches) { return false; }

    for (size_t i = 0; i < comps.size(); ++i) {
        if (!comps[i].is_constant() && !is_square_free(comps[i], main_var)) {
            return false;
        }
    }

    return true;
}
