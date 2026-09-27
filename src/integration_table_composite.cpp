#include "internal/integration_table_rules.hpp"

namespace LMCAS::integration_table_detail {

void load_logarithmic_rules(IntegrationTable& table) {
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::ln(detail::make_expression_ptr(u));
        auto x_ln_x = SymbolicExpr::multiply(detail::make_expression_ptr(u), SymbolicExpr::ln(detail::make_expression_ptr(u)));
        auto res = *sym_sub(*x_ln_x, u);
        table.add_entry(Category::Logarithmic, IntegrationEntry(
            "ln(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
    {
        auto u = wildcard("_u");
        auto inv_u = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(-1));
        auto pat = *SymbolicExpr::multiply(SymbolicExpr::ln(detail::make_expression_ptr(u)), inv_u);
        auto ln_sq = SymbolicExpr::power(SymbolicExpr::ln(detail::make_expression_ptr(u)), SymbolicExpr::number(2));
        auto res = *SymbolicExpr::multiply(ln_sq, sym_rational(1, 2));
        table.add_entry(Category::Logarithmic, IntegrationEntry(
            "ln(x)/x", pat, res, {"_u"}, u_is_var("_u"), 50));
    }
}

void load_exponential_base_rules(IntegrationTable& table) {
    {
        auto a = wildcard("_a");
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(detail::make_expression_ptr(a), detail::make_expression_ptr(u));
        auto inv_ln_a = SymbolicExpr::power(SymbolicExpr::ln(detail::make_expression_ptr(a)), SymbolicExpr::number(-1));
        auto a_pow_u = SymbolicExpr::power(detail::make_expression_ptr(a), detail::make_expression_ptr(u));
        auto res = *SymbolicExpr::multiply(a_pow_u, inv_ln_a);
        table.add_entry(Category::Exponential, IntegrationEntry(
            "a^x", pat, res, {"_a", "_u"}, u_is_var_a_indep("_u", "_a"), 70));
    }
}

void load_scaled_radical_rules(IntegrationTable& table) {
    {
        auto u = wildcard("_u");
        auto a = wildcard("_a");
        auto a_sq = SymbolicExpr::power(detail::make_expression_ptr(a), SymbolicExpr::number(2));
        auto u_sq = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto neg_u_sq = SymbolicExpr::multiply(SymbolicExpr::number(-1), u_sq);
        auto a_sq_minus_u_sq = SymbolicExpr::add(a_sq, neg_u_sq);
        auto pat = *SymbolicExpr::power(SymbolicExpr::sqrt(a_sq_minus_u_sq), SymbolicExpr::number(-1));

        auto u_over_a = SymbolicExpr::multiply(detail::make_expression_ptr(u),
            SymbolicExpr::power(detail::make_expression_ptr(a), SymbolicExpr::number(-1)));
        auto res = *make_unary_function(FT::ArcSin, u_over_a);
        table.add_entry(Category::Algebraic, IntegrationEntry(
            "1/sqrt(a^2-x^2)", pat, res, {"_u", "_a"},
            u_is_var_a_indep("_u", "_a"), 50));
    }
}

void load_scaled_trig_rules(IntegrationTable& table) {
    {
        auto a = wildcard("_a");
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::sin(SymbolicExpr::multiply(detail::make_expression_ptr(a), detail::make_expression_ptr(u)));
        auto inv_a = SymbolicExpr::power(detail::make_expression_ptr(a), SymbolicExpr::number(-1));
        auto neg_inv_a = SymbolicExpr::multiply(SymbolicExpr::number(-1), inv_a);
        auto cos_ax = SymbolicExpr::cos(SymbolicExpr::multiply(detail::make_expression_ptr(a), detail::make_expression_ptr(u)));
        auto res = *SymbolicExpr::multiply(neg_inv_a, cos_ax);
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "sin(a*x)", pat, res, {"_a", "_u"},
            u_is_var_a_indep("_u", "_a"), 55));
    }
    {
        auto a = wildcard("_a");
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::cos(SymbolicExpr::multiply(detail::make_expression_ptr(a), detail::make_expression_ptr(u)));
        auto inv_a = SymbolicExpr::power(detail::make_expression_ptr(a), SymbolicExpr::number(-1));
        auto sin_ax = SymbolicExpr::sin(SymbolicExpr::multiply(detail::make_expression_ptr(a), detail::make_expression_ptr(u)));
        auto res = *SymbolicExpr::multiply(inv_a, sin_ax);
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "cos(a*x)", pat, res, {"_a", "_u"},
            u_is_var_a_indep("_u", "_a"), 55));
    }
    {
        auto u = wildcard("_u");
        auto u_sq = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto pat = *SymbolicExpr::multiply(detail::make_expression_ptr(u), SymbolicExpr::cos(u_sq));

        auto u_sq_r = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto res = *SymbolicExpr::multiply(sym_rational(1, 2), SymbolicExpr::sin(u_sq_r));
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "x*cos(x^2)", pat, res, {"_u"}, u_is_var("_u"), 35));
    }
}

void load_exponential_trig_rules(IntegrationTable& table) {
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::multiply(
            SymbolicExpr::exp(detail::make_expression_ptr(u)),
            SymbolicExpr::sin(detail::make_expression_ptr(u)));

        auto sin_u = SymbolicExpr::sin(detail::make_expression_ptr(u));
        auto neg_cos_u = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::cos(detail::make_expression_ptr(u)));
        auto sin_minus_cos = SymbolicExpr::add(sin_u, neg_cos_u);
        auto exp_u = SymbolicExpr::exp(detail::make_expression_ptr(u));
        auto res = *SymbolicExpr::multiply(
            sym_rational(1, 2),
            SymbolicExpr::multiply(exp_u, sin_minus_cos));
        table.add_entry(Category::Exponential, IntegrationEntry(
            "exp(x)*sin(x)", pat, res, {"_u"}, u_is_var("_u"), 35));
    }
}

}
