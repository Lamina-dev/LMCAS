#include "internal/integration_rational_support.hpp"

namespace LMCAS {

namespace {
bool rd_factor_quadratic(Polynomial<Rational>& current, int mult,
    std::vector<std::pair<Polynomial<Rational>, int>>& factors_out) {
                auto roots = find_rational_roots(current);
                if (!roots.empty()) {
                    for (const auto& r : roots) {
                        Polynomial<Rational> linear({Rational(0) - r, Rational(1)}, current.variable_name);
                        auto qr = current.div_mod(linear);
                        if (!qr.second.is_zero()) {
                            return false;
                        }
                        factors_out.push_back({linear, mult});
                        current = qr.first;
                        if (current.degree() == 0) { break; }
                    }
                }
                if (current.degree() == 2) {
                    factors_out.push_back({current.make_monic(), mult});
                    current = Polynomial<Rational>({Rational(1)}, current.variable_name);
                } else if (current.degree() == 1) {
                    factors_out.push_back({current.make_monic(), mult});
                    current = Polynomial<Rational>({Rational(1)}, current.variable_name);
                }
    return true;
}

bool rd_factor_piece(const Polynomial<Rational>& piece, int mult,
    std::vector<std::pair<Polynomial<Rational>, int>>& factors_out) {
    auto current = piece.make_monic();
    if (current.degree() == 0) { return true; }
    if (current.degree() == 1) {
        factors_out.push_back({current, mult});
        return true;
    }
    if (current.degree() == 2) { return rd_factor_quadratic(current, mult, factors_out); }
        auto roots = find_rational_roots(current);
        for (const auto& r : roots) {
            Polynomial<Rational> linear({Rational(0) - r, Rational(1)}, current.variable_name);
            auto qr = current.div_mod(linear);
            if (qr.second.is_zero()) {
                factors_out.push_back({linear, mult});
                current = qr.first;
                if (current.degree() <= 0) { break; }
            }
        }

        if (current.degree() <= 0) { return true; }
        if (current.degree() == 1) {
            factors_out.push_back({current.make_monic(), mult});
            return true;
        }
        if (current.degree() == 2) {
            factors_out.push_back({current.make_monic(), mult});
            return true;
        }
        return false;
}
}

bool RationalDecompositionStrategy::factor_denominator(
    const Polynomial<Rational>& Q,
    std::vector<std::pair<Polynomial<Rational>, int>>& factors_out) {
    factors_out.clear();
    if (Q.degree() <= 0) { return false; }
    auto sqfree = square_free_factorization(Q);
    if (sqfree.empty()) { return false; }

    /**
     * @brief 逐次提取各无平方因子的有理根。
     * 剩余因子次数至多为 2 时匹配成功，否则返回未匹配。
     */
    for (const auto& [piece, mult] : sqfree) {
        if (!rd_factor_piece(piece, mult, factors_out)) { return false; }
    }
    return true;
}

}
