#include "internal/integration_rational_support.hpp"

namespace LMCAS {
namespace {
Result<bool> rd_collect_rational(const std::shared_ptr<const SymbolicNode>& node,
                         const std::string& var,
                         std::vector<Polynomial<Rational>>& num,
                         std::vector<Polynomial<Rational>>& den);

Result<bool> rd_collect_sum(const AddNode& add, const std::string& var,
    std::vector<Polynomial<Rational>>& num, std::vector<Polynomial<Rational>>& den) {
        std::vector<std::pair<std::vector<Polynomial<Rational>>,
                              std::vector<Polynomial<Rational>>>> parts;
        parts.reserve(add.operands().size());
        for (const auto& op : add.operands()) {
            std::vector<Polynomial<Rational>> sub_num, sub_den;
            auto collected = rd_collect_rational(op, var, sub_num, sub_den);
            if (!collected || !collected.value()) { return collected; }
            parts.push_back({std::move(sub_num), std::move(sub_den)});
        }

        Polynomial<Rational> total_num(var);
        Polynomial<Rational> total_den({Rational(1)}, var);

        std::vector<Polynomial<Rational>> Q_per_term;
        Q_per_term.reserve(parts.size());
        for (const auto& [pn, pd] : parts) {
            Polynomial<Rational> q({Rational(1)}, var);
            for (const auto& d : pd) q = q * d;
            Q_per_term.push_back(q);
            total_den = total_den * q;
        }

        for (size_t i = 0; i < parts.size(); ++i) {
            Polynomial<Rational> p({Rational(1)}, var);
            for (const auto& nn : parts[i].first) p = p * nn;
            for (size_t j = 0; j < parts.size(); ++j) {
                if (j == i) { continue; }
                p = p * Q_per_term[j];
            }
            total_num = total_num + p;
        }

        num.push_back(total_num);
        den.push_back(total_den);
        return true;
}

Result<bool> rd_collect_power(const PowerNode& pw, const std::string& var,
    std::vector<Polynomial<Rational>>& num, std::vector<Polynomial<Rational>>& den) {
        auto en = std::dynamic_pointer_cast<const NumberNode>(pw.exponent());
        if (!en) { return false; }

        const auto bounded_exponent = integration_bounded_exponent(*en, -64, 64);
        if (!bounded_exponent) { return false; }

        std::vector<Polynomial<Rational>> bn, bd;
        auto collected = rd_collect_rational(pw.base(), var, bn, bd);
        if (!collected || !collected.value()) { return collected; }

        Polynomial<Rational> base_num({Rational(1)}, var);
        for (const auto& p : bn) base_num = base_num * p;
        Polynomial<Rational> base_den({Rational(1)}, var);
        for (const auto& p : bd) base_den = base_den * p;

        if (*bounded_exponent == 0) {
            num.push_back(Polynomial<Rational>({Rational(1)}, var));
            return true;
        }
        const int power_count = *bounded_exponent < 0 ? -*bounded_exponent : *bounded_exponent;
        Polynomial<Rational> n_pow({Rational(1)}, var);
        Polynomial<Rational> d_pow({Rational(1)}, var);
        for (int i = 0; i < power_count; ++i) {
            n_pow = n_pow * base_num;
            d_pow = d_pow * base_den;
        }
        if (*bounded_exponent > 0) {
            num.push_back(n_pow);
            den.push_back(d_pow);
        } else {
            num.push_back(d_pow);
            den.push_back(n_pow);
        }
        return true;
}

Result<bool> rd_collect_rational(const std::shared_ptr<const SymbolicNode>& node,
                         const std::string& var,
                         std::vector<Polynomial<Rational>>& num,
                         std::vector<Polynomial<Rational>>& den) {
    if (!node) { return false; }

    if (auto n = std::dynamic_pointer_cast<const NumberNode>(node)) {
        auto p = symbolic_to_poly<Rational>(LMCAS::detail::make_expression_ptr(n), var);
        if (!p) {
            if (p.error().code == CasErrc::UnsupportedExpression) { return false; }
            return Result<bool>::failure(p.error());
        }
        num.push_back(std::move(p.value()));
        return true;
    }
    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (v->name() == var) {
            num.push_back(Polynomial<Rational>({Rational(0), Rational(1)}, var));
            return true;
        }
        return false;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return rd_collect_sum(*add, var, num, den);
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            auto collected = rd_collect_rational(op, var, num, den);
            if (!collected || !collected.value()) { return collected; }
        }
        return true;
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return rd_collect_power(*power, var, num, den);
    }
    return false;
}

void rd_normalize_fraction(Polynomial<Rational>& P, Polynomial<Rational>& Q,
    const std::string& var) {
    if (!P.is_zero()) {
        Polynomial<Rational> g = Polynomial<Rational>::gcd(P, Q);
        if (!g.is_zero() && g.degree() >= 1) {
            auto pq = P.div_mod(g);
            auto qq = Q.div_mod(g);
            if (pq.second.is_zero() && qq.second.is_zero()) {
                P = pq.first;
                Q = qq.first;
            }
        }
    }

    if (Q.degree() >= 0) {
        Rational lc = Q.lead_coeff();
        if (!(lc == Rational(1)) && !(lc == Rational(0))) {
            Polynomial<Rational> P_scaled(var);
            P_scaled.coeffs.reserve(P.coeffs.size());
            for (const auto& c : P.coeffs) { P_scaled.coeffs.push_back(c / lc); }
            Polynomial<Rational> Q_scaled(var);
            Q_scaled.coeffs.reserve(Q.coeffs.size());
            for (const auto& c : Q.coeffs) { Q_scaled.coeffs.push_back(c / lc); }
            P_scaled.trim();
            Q_scaled.trim();
            P = P_scaled;
            Q = Q_scaled;
        }
    }
}

}

Result<bool> RationalDecompositionStrategy::extract_rational(
    const SymbolicExpr& expr, const std::string& var,
    Polynomial<Rational>& P_out, Polynomial<Rational>& Q_out) {

    std::vector<Polynomial<Rational>> nums, dens;
    auto collected = rd_collect_rational(LMCAS::detail::node(expr), var, nums, dens);
    if (!collected || !collected.value()) { return collected; }

    Polynomial<Rational> P({Rational(1)}, var);
    for (const auto& p : nums) { P = P * p; }
    Polynomial<Rational> Q({Rational(1)}, var);
    for (const auto& p : dens) { Q = Q * p; }

    if (Q.is_zero()) { return false; } /** @brief 分母为零时返回 false。 */

    rd_normalize_fraction(P, Q, var);

    P_out = P;
    Q_out = Q;
    return true;
}

void RationalDecompositionStrategy::poly_divide(
    const Polynomial<Rational>& P, const Polynomial<Rational>& Q,
    Polynomial<Rational>& quotient_out,
    Polynomial<Rational>& remainder_out) {
    auto qr = P.div_mod(Q);
    quotient_out = qr.first;
    remainder_out = qr.second;
}

}
