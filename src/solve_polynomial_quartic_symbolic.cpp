#include "internal/polynomial_solver_support.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace LMCAS::polynomial_solver_detail {

SymbolicQuarticDepression symbolic_quartic_depression(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d,
    const std::shared_ptr<SymbolicExpr>& e) {
    auto a2 = SymbolicExpr::power(a, num(2));
    auto a3 = SymbolicExpr::power(a, num(3));
    auto a4 = SymbolicExpr::power(a, num(4));
    auto b2 = SymbolicExpr::power(b, num(2));
    auto b3 = SymbolicExpr::power(b, num(3));
    auto b4 = SymbolicExpr::power(b, num(4));

    auto eight_ac = SymbolicExpr::multiply(num(8), SymbolicExpr::multiply(a, c));
    auto three_b2 = SymbolicExpr::multiply(num(3), b2);
    auto p = SymbolicExpr::divide(sub(eight_ac, three_b2), SymbolicExpr::multiply(num(8), a2))->simplify();

    auto four_abc = SymbolicExpr::multiply(num(4), SymbolicExpr::multiply(a, SymbolicExpr::multiply(b, c)));
    auto eight_a2d = SymbolicExpr::multiply(num(8), SymbolicExpr::multiply(a2, d));
    auto q_num = SymbolicExpr::add(sub(b3, four_abc), eight_a2d);
    auto q = SymbolicExpr::divide(q_num, SymbolicExpr::multiply(num(8), a3))->simplify();

    auto neg3_b4 = SymbolicExpr::multiply(num(-3), b4);
    auto t256_a3e = SymbolicExpr::multiply(num(256), SymbolicExpr::multiply(a3, e));
    auto neg64_a2bd = SymbolicExpr::multiply(num(-64), SymbolicExpr::multiply(a2, SymbolicExpr::multiply(b, d)));
    auto t16_ab2c = SymbolicExpr::multiply(num(16), SymbolicExpr::multiply(a, SymbolicExpr::multiply(b2, c)));
    auto r_num = SymbolicExpr::add(SymbolicExpr::add(SymbolicExpr::add(neg3_b4, t256_a3e), neg64_a2bd), t16_ab2c);
    auto r = SymbolicExpr::divide(r_num, SymbolicExpr::multiply(num(256), a4))->simplify();

    auto shift = SymbolicExpr::divide(b, SymbolicExpr::multiply(num(4), a))->simplify();
    return {p, q, r, shift};
}

std::vector<std::shared_ptr<SymbolicExpr>> symbolic_quartic_roots(
    const SymbolicQuarticDepression& depression, const std::string& var) {
    const auto& [p, q, r, shift] = depression;
    auto q_simp = q->simplify();
    if (q_simp->get_number_value_is_zero()) {

        auto u_roots = solve_quadratic_internal(num(1), p, r);
        std::vector<std::shared_ptr<SymbolicExpr>> results;
        for (const auto& u : u_roots) {

            auto pos_y = SymbolicExpr::sqrt(u)->simplify();
            auto neg_y = negate(SymbolicExpr::sqrt(u))->simplify();
            results.push_back(sub(pos_y, shift)->simplify());
            results.push_back(sub(neg_y, shift)->simplify());
        }
        return results;
    }

    auto p2 = SymbolicExpr::power(p, num(2));
    auto q2 = SymbolicExpr::power(q, num(2));
    auto eight_p = SymbolicExpr::multiply(num(8), p);
    auto two_p2_minus_8r = sub(SymbolicExpr::multiply(num(2), p2), SymbolicExpr::multiply(num(8), r))->simplify();
    auto neg_q2 = negate(q2)->simplify();

    auto cubic_roots = solve_cubic(num(8), eight_p, two_p2_minus_8r, neg_q2, var);

    std::shared_ptr<SymbolicExpr> m = nullptr;
    if (!cubic_roots.empty()) {
        m = cubic_roots[0];
    }

    if (!m) return {};

    auto two_m = SymbolicExpr::multiply(num(2), m);
    auto s = SymbolicExpr::sqrt(two_m)->simplify();

    auto p_half = SymbolicExpr::divide(p, num(2))->simplify();

    auto m_plus_p_half = SymbolicExpr::add(m, p_half)->simplify();

    auto q_over_2s = SymbolicExpr::divide(q, SymbolicExpr::multiply(num(2), s))->simplify();

    auto quad1_c = sub(m_plus_p_half, q_over_2s)->simplify();

    auto neg_s = negate(s)->simplify();
    auto quad2_c = SymbolicExpr::add(m_plus_p_half, q_over_2s)->simplify();

    auto y_roots1 = solve_quadratic_internal(num(1), s, quad1_c);
    auto y_roots2 = solve_quadratic_internal(num(1), neg_s, quad2_c);

    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (const auto& y : y_roots1) {
        results.push_back(sub(y, shift)->simplify());
    }
    for (const auto& y : y_roots2) {
        results.push_back(sub(y, shift)->simplify());
    }

    return results;
}

}
