#include "root_of_utils.hpp"
#include "solve_polynomial.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/exact_root.hpp"

namespace LMCAS {






static std::vector<std::shared_ptr<SymbolicExpr>> solve_closed_form_from_poly(
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
        auto a = get_coeff(2);
        auto b = get_coeff(1);
        auto c = get_coeff(0);

        auto b2 = SymbolicExpr::power(b, SymbolicExpr::number(2));
        auto four_ac = SymbolicExpr::multiply(SymbolicExpr::number(4), SymbolicExpr::multiply(a, c));
        auto delta = SymbolicExpr::add(b2, SymbolicExpr::multiply(four_ac, SymbolicExpr::number(-1)));
        auto sqrt_delta = SymbolicExpr::sqrt(delta);
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        auto two_a = SymbolicExpr::multiply(SymbolicExpr::number(2), a);

        auto x1 = SymbolicExpr::divide(SymbolicExpr::add(neg_b, sqrt_delta), two_a)->simplify();
        auto x2 = SymbolicExpr::divide(SymbolicExpr::add(neg_b, SymbolicExpr::multiply(sqrt_delta, SymbolicExpr::number(-1))), two_a)->simplify();
        return { x1, x2 };
    } else if (deg == 3) {
        auto a = get_coeff(3);
        auto b = get_coeff(2);
        auto c = get_coeff(1);
        auto d = get_coeff(0);
        return solve_cubic(a, b, c, d, var);
    } else if (deg == 4) {
        auto a = get_coeff(4);
        auto b = get_coeff(3);
        auto c = get_coeff(2);
        auto d = get_coeff(1);
        auto e = get_coeff(0);
        return solve_quartic(a, b, c, d, e, var);
    }

    return {};
}




std::shared_ptr<SymbolicExpr> rootof_simplify(
    const std::shared_ptr<SymbolicExpr>& rootof_expr)
{
    if (!rootof_expr || !LMCAS::detail::node(rootof_expr)) return rootof_expr;

    auto root = std::dynamic_pointer_cast<const RootOfNode>(
        LMCAS::detail::node(rootof_expr));
    if (!root) return rootof_expr;

    const std::string& var = root->variable();
    const int k = static_cast<int>(root->index());
    std::vector<SymbolicPolyCoeff> coefficients;
    coefficients.reserve(root->exact_id().polynomial.coeffs.size());
    for (const auto& coefficient : root->exact_id().polynomial.coeffs) {
        coefficients.emplace_back(SymbolicExpr::number(coefficient));
    }
    Polynomial<SymbolicPolyCoeff> sym_poly(coefficients, var);

    int degree = sym_poly.degree();
    if (degree <= 0) return rootof_expr;

    if (k < 0 || k >= degree) {
        return rootof_expr;
    }

    if (degree == 1) {
        auto roots = solve_closed_form_from_poly(sym_poly, var);
        return roots.size() == 1 ? roots.front() : rootof_expr;
    }
    if (degree == 2) {
        auto roots = solve_closed_form_from_poly(sym_poly, var);
        if (roots.size() != 2) return rootof_expr;
        return roots[static_cast<std::size_t>(1 - k)];
    }

    return rootof_expr;
}

}
