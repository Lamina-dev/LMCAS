#include "internal/expression_transform.hpp"
#include "internal/expression_analysis.hpp"
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
    const std::string& var, ComputationContext& context)
{
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
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
        return { detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(neg_b, a)), transform_budget)) };
    } else if (deg == 2) {
        return solve_quadratic_internal(get_coeff(2), get_coeff(1), get_coeff(0), context);
    } else if (deg == 3) {
        return solve_cubic_impl(get_coeff(3), get_coeff(2), get_coeff(1), get_coeff(0), var, context);
    } else if (deg == 4) {
        return solve_quartic_impl(get_coeff(4), get_coeff(3), get_coeff(2), get_coeff(1), get_coeff(0), var, context);
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
    const std::string& var, std::vector<std::shared_ptr<SymbolicExpr>>& results,
    ComputationContext& context) {
    if (factor.degree() <= 0) { return; }
    auto symbolic = rational_to_symbolic_poly(factor);
    auto solutions = factor.degree() <= 4 ?
        solve_closed_form_poly(symbolic, var, context) :
        detail::make_rootof_solutions_impl(symbolic, var, context);
    for (int copy = 0; copy < multiplicity; ++copy) {
        results.insert(results.end(), solutions.begin(), solutions.end());
    }
}

std::vector<std::shared_ptr<SymbolicExpr>> polynomial_solver_detail::solve_by_factoring_impl(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var, ComputationContext& context) {

    if (poly.is_zero() || poly.degree() < 1) { return {}; }

    if (poly.degree() <= 2) {
        return solve_closed_form_poly(poly, var, context);
    }

    std::vector<std::shared_ptr<SymbolicExpr>> results;

    Polynomial<Rational> rat_poly;
    if (!convert_to_rational_poly(poly, rat_poly, context)) {
        if (poly.degree() <= 4) {
            return solve_closed_form_poly(poly, var, context);
        }
        return detail::make_rootof_solutions_impl(poly, var, context);
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
        auto factor_roots = solve_closed_form_poly(sym_quotient, var, context);
        results.insert(results.end(), factor_roots.begin(), factor_roots.end());
        return results;
    }

    auto sqfree_factors = square_free_factorization(quotient);

    for (const auto& [factor, multiplicity] : sqfree_factors) {
        append_factor_solutions(factor, multiplicity, var, results, context);
    }

    return results;
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_by_factoring(
    const Polynomial<SymbolicPolyCoeff>& poly, const std::string& variable) {
    ComputationContext context;
    return solve_by_factoring_impl(poly, variable, context);
}

}
