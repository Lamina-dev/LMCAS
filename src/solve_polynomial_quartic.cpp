#include "internal/expression_transform.hpp"
#include "internal/polynomial_solver_support.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace LMCAS {
using namespace polynomial_solver_detail;

std::vector<std::shared_ptr<SymbolicExpr>> polynomial_solver_detail::solve_quartic_impl(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e,
    const std::string& var, ComputationContext& context) {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");

    auto a_simp = detail::make_expression_ptr(detail::normalize_node(detail::node(a), transform_budget));
    if (a_simp->get_number_value_is_zero()) {

        return solve_cubic_impl(b, c, d, e, var, context);
    }

    {
        Polynomial<SymbolicPolyCoeff> symbolic_poly(
            {SymbolicPolyCoeff(e), SymbolicPolyCoeff(d), SymbolicPolyCoeff(c),
             SymbolicPolyCoeff(b), SymbolicPolyCoeff(a)},
            var);
        Polynomial<Rational> rational_poly;
        if (convert_to_rational_poly(symbolic_poly, rational_poly, context) &&
            rational_poly.degree() == 4 && bounded_rational_root_search(rational_poly)) {
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

    auto b_simp = detail::make_expression_ptr(detail::normalize_node(detail::node(b), transform_budget));
    auto d_simp = detail::make_expression_ptr(detail::normalize_node(detail::node(d), transform_budget));
    if (b_simp->get_number_value_is_zero() && d_simp->get_number_value_is_zero()) {
        return solve_biquadratic_impl(a, c, e, var, context);
    }

    auto depression = numeric_quartic_depression(a, b, c, d, e, context);
    if (depression) {
        auto numeric_roots = numeric_quartic_roots(*depression, var, context);
        if (!numeric_roots.empty()) return numeric_roots;
    }
    return symbolic_quartic_roots(
        symbolic_quartic_depression(a, b, c, d, e, context), var, context);
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_quartic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e,
    const std::string& variable) {
    ComputationContext context;
    return solve_quartic_impl(a, b, c, d, e, variable, context);
}

}
