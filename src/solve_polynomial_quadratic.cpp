#include "internal/expression_transform.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/polynomial_solver_support.hpp"
#include <stdexcept>

namespace LMCAS {
using namespace polynomial_solver_detail;

namespace polynomial_solver_detail {
std::vector<std::shared_ptr<SymbolicExpr>> solve_quadratic_internal(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, ComputationContext& context) {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    auto b2 = SymbolicExpr::power(b, num(2));
    auto four_ac = SymbolicExpr::multiply(num(4), SymbolicExpr::multiply(a, c));
    auto delta = detail::make_expression_ptr(detail::normalize_node(detail::node(sub(b2, four_ac)), transform_budget));

    auto neg_b = negate(b);
    auto two_a = SymbolicExpr::multiply(num(2), a);
    if (delta && delta->is_zero()) {
        return { detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(neg_b, two_a)), transform_budget)) };
    }

    auto sqrt_delta = SymbolicExpr::sqrt(delta);
    auto exact_b = detail::exact_rational_value(detail::node(b));
    auto exact_delta = detail::exact_rational_value(detail::node(delta));
    if (exact_b && *exact_b != Rational(0) &&
        exact_delta && *exact_delta > Rational(0)) {
        auto numerator = *exact_b < Rational(0)
            ? SymbolicExpr::add(neg_b, sqrt_delta)
            : sub(neg_b, sqrt_delta);
        auto x1 = detail::make_expression_ptr(detail::normalize_node(
            detail::node(SymbolicExpr::divide(numerator, two_a)), transform_budget));
        auto x2 = detail::make_expression_ptr(detail::normalize_node(
            detail::node(SymbolicExpr::divide(
                SymbolicExpr::multiply(num(2), c), numerator)), transform_budget));
        return {x1, x2};
    }
    auto x1 = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(SymbolicExpr::add(neg_b, sqrt_delta), two_a)), transform_budget));
    auto x2 = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(sub(neg_b, sqrt_delta), two_a)), transform_budget));

    return { x1, x2 };
}
}

std::vector<std::shared_ptr<SymbolicExpr>> polynomial_solver_detail::solve_biquadratic_impl(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::string&, ComputationContext& context) {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    auto a_simplified = detail::make_expression_ptr(detail::normalize_node(detail::node(a), transform_budget));
    if (a_simplified->get_number_value_is_zero()) {
        auto b_simplified = detail::make_expression_ptr(detail::normalize_node(detail::node(b), transform_budget));
        if (!b_simplified->get_number_value_is_zero()) {
            return solve_quadratic_internal(b, num(0), c, context);
        }

        auto c_simplified = detail::make_expression_ptr(detail::normalize_node(detail::node(c), transform_budget));
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

    auto u_roots = solve_quadratic_internal(a, b, c, context);

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (const auto& u : u_roots) {
        auto pos_root = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::sqrt(u)), transform_budget));
        auto neg_root = detail::make_expression_ptr(detail::normalize_node(detail::node(negate(SymbolicExpr::sqrt(u))), transform_budget));
        results.push_back(pos_root);
        results.push_back(neg_root);
    }

    return results;
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_biquadratic(
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c,
    const std::string& variable) {
    ComputationContext context;
    return solve_biquadratic_impl(a, b, c, variable, context);
}

}
