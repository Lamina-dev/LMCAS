#include "internal/polynomial_solver_support.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "root_of_utils.hpp"

namespace LMCAS {
using namespace polynomial_solver_detail;

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

static void append_factor_solutions(const Polynomial<Rational>& factor, int multiplicity,
    const std::string& var, std::vector<std::shared_ptr<SymbolicExpr>>& results) {
    if (factor.degree() <= 0) { return; }
    auto symbolic = rational_to_symbolic_poly(factor);
    auto solutions = factor.degree() <= 4 ?
        solve_closed_form_poly(symbolic, var) : make_rootof_solutions(symbolic, var);
    for (int copy = 0; copy < multiplicity; ++copy) {
        results.insert(results.end(), solutions.begin(), solutions.end());
    }
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_by_factoring(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var) {

    if (poly.is_zero() || poly.degree() < 1) { return {}; }

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
    auto rational_roots = bounded_rational_root_search(rat_poly) ?
        find_rational_roots(rat_poly) : std::vector<Rational>{};

    for (const auto& r : rational_roots) {
        results.push_back(SymbolicExpr::number(r));
    }

    Polynomial<Rational> quotient = rat_poly;
    for (const auto& r : rational_roots) {
        Polynomial<Rational> linear_factor({-r, Rational(1)}, var);
        auto [q, rem] = quotient.div_mod(linear_factor);
        if (!rem.is_zero()) { throw std::logic_error("solve_by_factoring: inexact root division"); }
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
        append_factor_solutions(factor, multiplicity, var, results);
    }

    return results;
}

}
