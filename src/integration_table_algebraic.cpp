#include "internal/integration_table_rules.hpp"

namespace LMCAS::integration_table_detail {

void load_algebraic_elementary_rules(IntegrationTable& table) {
    {
        auto u = wildcard("_u");
        auto u_sq = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto neg_u_sq = SymbolicExpr::multiply(SymbolicExpr::number(-1), u_sq);
        auto one_minus_u_sq = SymbolicExpr::add(SymbolicExpr::number(1), neg_u_sq);
        auto pat = *SymbolicExpr::power(SymbolicExpr::sqrt(one_minus_u_sq), SymbolicExpr::number(-1));
        auto res = *make_unary_function(FT::ArcSin, detail::make_expression_ptr(u));
        table.add_entry(Category::Algebraic, IntegrationEntry(
            "1/sqrt(1-x^2)", pat, res, {"_u"}, u_is_var("_u"), 40));
    }
    {
        auto u = wildcard("_u");
        auto pat = *make_unary_function(FT::Abs, detail::make_expression_ptr(u));
        auto res = *SymbolicExpr::multiply(
            SymbolicExpr::number(Rational(1, 2)),
            SymbolicExpr::multiply(detail::make_expression_ptr(u), make_unary_function(FT::Abs, detail::make_expression_ptr(u))));
        table.add_entry(Category::Algebraic, IntegrationEntry(
            "|x|", pat, res, {"_u"}, u_is_var("_u"), 40));
    }
    {
        auto u = wildcard("_u");
        auto u_sq = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto one_plus_u_sq = SymbolicExpr::add(SymbolicExpr::number(1), u_sq);
        auto pat = *SymbolicExpr::power(one_plus_u_sq, SymbolicExpr::number(-1));
        auto res = *make_unary_function(FT::ArcTan, detail::make_expression_ptr(u));
        table.add_entry(Category::Algebraic, IntegrationEntry(
            "1/(1+x^2)", pat, res, {"_u"}, u_is_var("_u"), 40));
    }
}

void load_algebraic_radicals_rules(IntegrationTable& table) {
    {
        auto u = wildcard("_u");
        auto a = wildcard("_a");
        auto u_sq = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto a_sq = SymbolicExpr::power(detail::make_expression_ptr(a), SymbolicExpr::number(2));
        auto sum = SymbolicExpr::add(u_sq, a_sq);
        auto pat = *SymbolicExpr::power(SymbolicExpr::sqrt(sum), SymbolicExpr::number(-1));

        auto u_sq_r = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto a_sq_r = SymbolicExpr::power(detail::make_expression_ptr(a), SymbolicExpr::number(2));
        auto sum_r = SymbolicExpr::add(u_sq_r, a_sq_r);
        auto inner = SymbolicExpr::add(detail::make_expression_ptr(u), SymbolicExpr::sqrt(sum_r));
        auto res = *SymbolicExpr::ln(inner);
        table.add_entry(Category::Algebraic, IntegrationEntry(
            "1/sqrt(x^2+a^2)", pat, res, {"_u", "_a"},
            u_is_var_a_indep("_u", "_a"), 50));
    }
    {
        auto u = wildcard("_u");
        auto a = wildcard("_a");
        auto u_sq = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto a_sq = SymbolicExpr::power(detail::make_expression_ptr(a), SymbolicExpr::number(2));
        auto neg_a_sq = SymbolicExpr::multiply(SymbolicExpr::number(-1), a_sq);
        auto diff = SymbolicExpr::add(u_sq, neg_a_sq);
        auto pat = *SymbolicExpr::power(SymbolicExpr::sqrt(diff), SymbolicExpr::number(-1));

        auto u_sq_r = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto a_sq_r = SymbolicExpr::power(detail::make_expression_ptr(a), SymbolicExpr::number(2));
        auto neg_a_sq_r = SymbolicExpr::multiply(SymbolicExpr::number(-1), a_sq_r);
        auto diff_r = SymbolicExpr::add(u_sq_r, neg_a_sq_r);
        auto inner = SymbolicExpr::add(detail::make_expression_ptr(u), SymbolicExpr::sqrt(diff_r));
        auto res = *SymbolicExpr::ln(inner);
        table.add_entry(Category::Algebraic, IntegrationEntry(
            "1/sqrt(x^2-a^2)", pat, res, {"_u", "_a"},
            u_is_var_a_indep("_u", "_a"), 50));
    }
}

void load_algebraic_reciprocals_rules(IntegrationTable& table) {
    {
        auto u = wildcard("_u");
        auto a = wildcard("_a");
        auto u_sq = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto a_sq = SymbolicExpr::power(detail::make_expression_ptr(a), SymbolicExpr::number(2));
        auto sum = SymbolicExpr::add(u_sq, a_sq);
        auto pat = *SymbolicExpr::power(sum, SymbolicExpr::number(-1));

        auto inv_a = SymbolicExpr::power(detail::make_expression_ptr(a), SymbolicExpr::number(-1));
        auto u_over_a = SymbolicExpr::multiply(detail::make_expression_ptr(u), SymbolicExpr::power(detail::make_expression_ptr(a), SymbolicExpr::number(-1)));
        auto res = *SymbolicExpr::multiply(inv_a, make_unary_function(FT::ArcTan, u_over_a));
        table.add_entry(Category::Algebraic, IntegrationEntry(
            "1/(x^2+a^2)", pat, res, {"_u", "_a"},
            u_is_var_a_indep("_u", "_a"), 50));
    }
    {
        auto u = wildcard("_u");
        auto a = wildcard("_a");
        auto u_sq = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto a_sq = SymbolicExpr::power(detail::make_expression_ptr(a), SymbolicExpr::number(2));
        auto neg_a_sq = SymbolicExpr::multiply(SymbolicExpr::number(-1), a_sq);
        auto diff = SymbolicExpr::add(u_sq, neg_a_sq);
        auto pat = *SymbolicExpr::power(diff, SymbolicExpr::number(-1));

        auto two_a = SymbolicExpr::multiply(SymbolicExpr::number(2), detail::make_expression_ptr(a));
        auto inv_2a = SymbolicExpr::power(two_a, SymbolicExpr::number(-1));
        auto neg_a = SymbolicExpr::multiply(SymbolicExpr::number(-1), detail::make_expression_ptr(a));
        auto u_minus_a = SymbolicExpr::add(detail::make_expression_ptr(u), neg_a);
        auto u_plus_a = SymbolicExpr::add(detail::make_expression_ptr(u), detail::make_expression_ptr(a));
        auto frac = SymbolicExpr::divide(u_minus_a, u_plus_a);
        auto res = *SymbolicExpr::multiply(inv_2a, SymbolicExpr::ln(frac));
        table.add_entry(Category::Algebraic, IntegrationEntry(
            "1/(x^2-a^2)", pat, res, {"_u", "_a"},
            u_is_var_a_indep("_u", "_a"), 50));
    }
}

void load_polynomial_rules(IntegrationTable& table) {
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(-1));
        auto res = *SymbolicExpr::ln(detail::make_expression_ptr(u));
        table.add_entry(Category::Polynomial, IntegrationEntry(
            "1/x", pat, res, {"_u"}, u_is_var("_u"), 50));
    }
    {
        auto u = wildcard("_u");
        auto n = wildcard("_n");
        auto pat = *SymbolicExpr::power(detail::make_expression_ptr(u), detail::make_expression_ptr(n));
        auto n_plus_1 = SymbolicExpr::add(detail::make_expression_ptr(n), SymbolicExpr::number(1));
        auto res = *SymbolicExpr::divide(
            SymbolicExpr::power(detail::make_expression_ptr(u), n_plus_1),
            n_plus_1);
        table.add_entry(Category::Polynomial, IntegrationEntry(
            "x^n", pat, res, {"_u", "_n"},
            [](const MatchMap& m, const std::string& var) {
                auto it_u = m.find("_u");
                if (it_u == m.end()) { return false; }
                auto v = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(it_u->second));
                if (!v || v->is_constant() || v->name() != var) { return false; }
                auto it_n = m.find("_n");
                if (it_n == m.end()) { return false; }
                if (expression_depends_on_variable(LMCAS::detail::node(it_n->second), var)) { return false; }
                auto n_simp = it_n->second.simplify();
                if (!n_simp) { return false; }
                auto neg_one = SymbolicExpr::number(-1);
                auto diff = SymbolicExpr::add(n_simp, neg_one)->simplify();
                if (diff && diff->is_zero()) { return false; }
                return true;
            }, 80));
    }
    {
        auto u = wildcard("_u");
        auto pat = LMCAS::detail::expression_from_node(LMCAS::detail::node(u));
        auto res = *SymbolicExpr::multiply(
            SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2)),
            sym_rational(1, 2));
        table.add_entry(Category::Polynomial, IntegrationEntry(
            "x", pat, res, {"_u"}, u_is_var("_u"), 90));
    }
}

}
