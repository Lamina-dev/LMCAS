#include "internal/expression_transform.hpp"
#include "internal/polynomial_solver_support.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace LMCAS::polynomial_solver_detail {

SymbolicQuarticDepression symbolic_quartic_depression(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e, ComputationContext& context) {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    auto a2 = SymbolicExpr::power(a, num(2));
    auto a3 = SymbolicExpr::power(a, num(3));
    auto a4 = SymbolicExpr::power(a, num(4));
    auto b2 = SymbolicExpr::power(b, num(2));
    auto b3 = SymbolicExpr::power(b, num(3));
    auto b4 = SymbolicExpr::power(b, num(4));

    auto eight_ac = SymbolicExpr::multiply(num(8), SymbolicExpr::multiply(a, c));
    auto three_b2 = SymbolicExpr::multiply(num(3), b2);
    auto p = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(sub(eight_ac, three_b2), SymbolicExpr::multiply(num(8), a2))), transform_budget));

    auto four_abc = SymbolicExpr::multiply(num(4), SymbolicExpr::multiply(a, SymbolicExpr::multiply(b, c)));
    auto eight_a2d = SymbolicExpr::multiply(num(8), SymbolicExpr::multiply(a2, d));
    auto q_num = SymbolicExpr::add(sub(b3, four_abc), eight_a2d);
    auto q = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(q_num, SymbolicExpr::multiply(num(8), a3))), transform_budget));

    auto neg3_b4 = SymbolicExpr::multiply(num(-3), b4);
    auto t256_a3e = SymbolicExpr::multiply(num(256), SymbolicExpr::multiply(a3, e));
    auto neg64_a2bd = SymbolicExpr::multiply(num(-64), SymbolicExpr::multiply(a2, SymbolicExpr::multiply(b, d)));
    auto t16_ab2c = SymbolicExpr::multiply(num(16), SymbolicExpr::multiply(a, SymbolicExpr::multiply(b2, c)));
    auto r_num = SymbolicExpr::add(SymbolicExpr::add(SymbolicExpr::add(neg3_b4, t256_a3e), neg64_a2bd), t16_ab2c);
    auto r = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(r_num, SymbolicExpr::multiply(num(256), a4))), transform_budget));

    auto shift = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(b, SymbolicExpr::multiply(num(4), a))), transform_budget));
    return {p, q, r, shift};
}

std::vector<std::shared_ptr<SymbolicExpr>> symbolic_quartic_roots(
    const SymbolicQuarticDepression& depression, const std::string& var,
    ComputationContext& context) {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    const auto& [p, q, r, shift] = depression;
    auto q_simp = detail::make_expression_ptr(detail::normalize_node(detail::node(q), transform_budget));
    if (q_simp->get_number_value_is_zero()) {

        auto u_roots = solve_quadratic_internal(num(1), p, r, context);
        std::vector<std::shared_ptr<SymbolicExpr>> results;
        for (const auto& u : u_roots) {

            auto pos_y = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::sqrt(u)), transform_budget));
            auto neg_y = detail::make_expression_ptr(detail::normalize_node(detail::node(negate(SymbolicExpr::sqrt(u))), transform_budget));
            results.push_back(detail::make_expression_ptr(detail::normalize_node(detail::node(sub(pos_y, shift)), transform_budget)));
            results.push_back(detail::make_expression_ptr(detail::normalize_node(detail::node(sub(neg_y, shift)), transform_budget)));
        }
        return results;
    }

    auto p2 = SymbolicExpr::power(p, num(2));
    auto q2 = SymbolicExpr::power(q, num(2));
    auto eight_p = SymbolicExpr::multiply(num(8), p);
    auto two_p2_minus_8r = detail::make_expression_ptr(detail::normalize_node(detail::node(sub(SymbolicExpr::multiply(num(2), p2), SymbolicExpr::multiply(num(8), r))), transform_budget));
    auto neg_q2 = detail::make_expression_ptr(detail::normalize_node(detail::node(negate(q2)), transform_budget));

    auto cubic_roots = solve_cubic_impl(num(8), eight_p, two_p2_minus_8r, neg_q2, var, context);

    std::shared_ptr<SymbolicExpr> m = nullptr;
    if (!cubic_roots.empty()) {
        m = cubic_roots[0];
    }

    if (!m) return {};

    auto two_m = SymbolicExpr::multiply(num(2), m);
    auto s = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::sqrt(two_m)), transform_budget));

    auto p_half = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(p, num(2))), transform_budget));

    auto m_plus_p_half = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::add(m, p_half)), transform_budget));

    auto q_over_2s = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(q, SymbolicExpr::multiply(num(2), s))), transform_budget));

    auto quad1_c = detail::make_expression_ptr(detail::normalize_node(detail::node(sub(m_plus_p_half, q_over_2s)), transform_budget));

    auto neg_s = detail::make_expression_ptr(detail::normalize_node(detail::node(negate(s)), transform_budget));
    auto quad2_c = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::add(m_plus_p_half, q_over_2s)), transform_budget));

    auto y_roots1 = solve_quadratic_internal(num(1), s, quad1_c, context);
    auto y_roots2 = solve_quadratic_internal(num(1), neg_s, quad2_c, context);

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (const auto& y : y_roots1) {
        results.push_back(detail::make_expression_ptr(detail::normalize_node(detail::node(sub(y, shift)), transform_budget)));
    }
    for (const auto& y : y_roots2) {
        results.push_back(detail::make_expression_ptr(detail::normalize_node(detail::node(sub(y, shift)), transform_budget)));
    }

    return results;
}

}
