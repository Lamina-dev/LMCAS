#include "internal/integration_table_rules.hpp"

namespace LMCAS::integration_table_detail {

void load_exponential_rules(IntegrationTable& table) {
    {
        auto a = wildcard("_a");
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::exp(SymbolicExpr::multiply(detail::make_expression_ptr(a), detail::make_expression_ptr(u)));
        auto res = *SymbolicExpr::multiply(
            SymbolicExpr::power(detail::make_expression_ptr(a), SymbolicExpr::number(-1)),
            SymbolicExpr::exp(SymbolicExpr::multiply(detail::make_expression_ptr(a), detail::make_expression_ptr(u))));
        table.add_entry(Category::Exponential, IntegrationEntry(
            "exp(a*x)", pat, res, {"_a", "_u"},
            u_is_var_a_indep("_u", "_a"), 50));
    }
    {
        auto u = wildcard("_u");
        auto u_sq = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto pat = *SymbolicExpr::multiply(u_sq, SymbolicExpr::exp(detail::make_expression_ptr(u)));

        auto u_sq_r = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto two_u = SymbolicExpr::multiply(SymbolicExpr::number(-2), detail::make_expression_ptr(u));
        auto poly = SymbolicExpr::add(SymbolicExpr::add(u_sq_r, two_u), SymbolicExpr::number(2));
        auto res = *SymbolicExpr::multiply(poly, SymbolicExpr::exp(detail::make_expression_ptr(u)));
        table.add_entry(Category::Exponential, IntegrationEntry(
            "x^2*exp(x)", pat, res, {"_u"}, u_is_var("_u"), 40));
    }
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::multiply(detail::make_expression_ptr(u), SymbolicExpr::exp(detail::make_expression_ptr(u)));

        auto u_minus_1 = SymbolicExpr::add(detail::make_expression_ptr(u), SymbolicExpr::number(-1));
        auto res = *SymbolicExpr::multiply(u_minus_1, SymbolicExpr::exp(detail::make_expression_ptr(u)));
        table.add_entry(Category::Exponential, IntegrationEntry(
            "x*exp(x)", pat, res, {"_u"}, u_is_var("_u"), 50));
    }
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::exp(detail::make_expression_ptr(u));
        auto res = *SymbolicExpr::exp(detail::make_expression_ptr(u));
        table.add_entry(Category::Exponential, IntegrationEntry(
            "exp(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
}

void load_trig_derivatives_rules(IntegrationTable& table) {
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(make_unary_function(FT::Sec, detail::make_expression_ptr(u)), SymbolicExpr::number(2));
        auto res = *make_unary_function(FT::Tan, detail::make_expression_ptr(u));
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "sec(x)^2", pat, res, {"_u"}, u_is_var("_u"), 30));
    }
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(make_unary_function(FT::Csc, detail::make_expression_ptr(u)), SymbolicExpr::number(2));
        auto res = *SymbolicExpr::multiply(SymbolicExpr::number(-1), make_unary_function(FT::Cot, detail::make_expression_ptr(u)));
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "csc(x)^2", pat, res, {"_u"}, u_is_var("_u"), 30));
    }
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::multiply(make_unary_function(FT::Sec, detail::make_expression_ptr(u)), make_unary_function(FT::Tan, detail::make_expression_ptr(u)));
        auto res = *make_unary_function(FT::Sec, detail::make_expression_ptr(u));
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "sec(x)*tan(x)", pat, res, {"_u"}, u_is_var("_u"), 30));
    }
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::multiply(make_unary_function(FT::Csc, detail::make_expression_ptr(u)), make_unary_function(FT::Cot, detail::make_expression_ptr(u)));
        auto res = *SymbolicExpr::multiply(SymbolicExpr::number(-1), make_unary_function(FT::Csc, detail::make_expression_ptr(u)));
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "csc(x)*cot(x)", pat, res, {"_u"}, u_is_var("_u"), 30));
    }
}

void load_trig_squares_rules(IntegrationTable& table) {
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(SymbolicExpr::sin(detail::make_expression_ptr(u)), SymbolicExpr::number(2));

        auto u_half = SymbolicExpr::multiply(detail::make_expression_ptr(u), sym_rational(1, 2));
        auto two_u = SymbolicExpr::multiply(SymbolicExpr::number(2), detail::make_expression_ptr(u));
        auto sin_2u_over_4 = SymbolicExpr::multiply(SymbolicExpr::sin(two_u), sym_rational(-1, 4));
        auto res = *SymbolicExpr::add(u_half, sin_2u_over_4);
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "sin(x)^2", pat, res, {"_u"}, u_is_var("_u"), 35));
    }
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(SymbolicExpr::cos(detail::make_expression_ptr(u)), SymbolicExpr::number(2));

        auto u_half = SymbolicExpr::multiply(detail::make_expression_ptr(u), sym_rational(1, 2));
        auto two_u = SymbolicExpr::multiply(SymbolicExpr::number(2), detail::make_expression_ptr(u));
        auto sin_2u_over_4 = SymbolicExpr::multiply(SymbolicExpr::sin(two_u), sym_rational(1, 4));
        auto res = *SymbolicExpr::add(u_half, sin_2u_over_4);
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "cos(x)^2", pat, res, {"_u"}, u_is_var("_u"), 35));
    }
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(SymbolicExpr::tan(detail::make_expression_ptr(u)), SymbolicExpr::number(2));
        auto neg_u = SymbolicExpr::multiply(SymbolicExpr::number(-1), detail::make_expression_ptr(u));
        auto res = *SymbolicExpr::add(SymbolicExpr::tan(detail::make_expression_ptr(u)), neg_u);
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "tan(x)^2", pat, res, {"_u"}, u_is_var("_u"), 35));
    }
}

void load_trig_primitives_rules(IntegrationTable& table) {
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::sin(detail::make_expression_ptr(u));
        auto res = *SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::cos(detail::make_expression_ptr(u)));
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "sin(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::cos(detail::make_expression_ptr(u));
        auto res = *SymbolicExpr::sin(detail::make_expression_ptr(u));
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "cos(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::tan(detail::make_expression_ptr(u));
        auto res = *SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::ln(SymbolicExpr::cos(detail::make_expression_ptr(u))));
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "tan(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
    {
        auto u = wildcard("_u");
        auto sin_u = SymbolicExpr::sin(detail::make_expression_ptr(u));
        auto cos_u = SymbolicExpr::cos(detail::make_expression_ptr(u));
        auto pat = *SymbolicExpr::multiply(
            sin_u,
            SymbolicExpr::power(cos_u, SymbolicExpr::number(-1)));
        auto res = *SymbolicExpr::multiply(
            SymbolicExpr::number(-1), SymbolicExpr::ln(cos_u));
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "sin(x)/cos(x)", pat, res, {"_u"}, u_is_var("_u"), 61));
    }
    {
        auto u = wildcard("_u");
        auto pat = *make_unary_function(FT::Cot, detail::make_expression_ptr(u));
        auto res = *SymbolicExpr::ln(SymbolicExpr::sin(detail::make_expression_ptr(u)));
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "cot(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
    {
        auto u = wildcard("_u");
        auto sin_u = SymbolicExpr::sin(detail::make_expression_ptr(u));
        auto cos_u = SymbolicExpr::cos(detail::make_expression_ptr(u));
        auto pat = *SymbolicExpr::multiply(
            cos_u,
            SymbolicExpr::power(sin_u, SymbolicExpr::number(-1)));
        auto res = *SymbolicExpr::ln(sin_u);
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "cos(x)/sin(x)", pat, res, {"_u"}, u_is_var("_u"), 61));
    }
    {
        auto u = wildcard("_u");
        auto pat = *make_unary_function(FT::Sec, detail::make_expression_ptr(u));
        auto inner = SymbolicExpr::add(make_unary_function(FT::Sec, detail::make_expression_ptr(u)), make_unary_function(FT::Tan, detail::make_expression_ptr(u)));
        auto res = *SymbolicExpr::ln(inner);
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "sec(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
    {
        auto u = wildcard("_u");
        auto pat = *make_unary_function(FT::Csc, detail::make_expression_ptr(u));
        auto inner = SymbolicExpr::add(make_unary_function(FT::Csc, detail::make_expression_ptr(u)), make_unary_function(FT::Cot, detail::make_expression_ptr(u)));
        auto res = *SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::ln(inner));
        table.add_entry(Category::Trigonometric, IntegrationEntry(
            "csc(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
}

void load_inverse_trig_rules(IntegrationTable& table) {
    {
        auto u = wildcard("_u");
        auto pat = *make_unary_function(FT::ArcSin, detail::make_expression_ptr(u));
        auto x_asin = SymbolicExpr::multiply(detail::make_expression_ptr(u), make_unary_function(FT::ArcSin, detail::make_expression_ptr(u)));
        auto u_sq = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto neg_u_sq = SymbolicExpr::multiply(SymbolicExpr::number(-1), u_sq);
        auto one_minus_u_sq = SymbolicExpr::add(SymbolicExpr::number(1), neg_u_sq);
        auto root = SymbolicExpr::sqrt(one_minus_u_sq);
        auto res = *SymbolicExpr::add(x_asin, root);
        table.add_entry(Category::InverseTrig, IntegrationEntry(
            "arcsin(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
    {
        auto u = wildcard("_u");
        auto pat = *make_unary_function(FT::ArcCos, detail::make_expression_ptr(u));
        auto x_acos = SymbolicExpr::multiply(detail::make_expression_ptr(u), make_unary_function(FT::ArcCos, detail::make_expression_ptr(u)));
        auto u_sq = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto neg_u_sq = SymbolicExpr::multiply(SymbolicExpr::number(-1), u_sq);
        auto one_minus_u_sq = SymbolicExpr::add(SymbolicExpr::number(1), neg_u_sq);
        auto neg_root = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::sqrt(one_minus_u_sq));
        auto res = *SymbolicExpr::add(x_acos, neg_root);
        table.add_entry(Category::InverseTrig, IntegrationEntry(
            "arccos(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
    {
        auto u = wildcard("_u");
        auto pat = *make_unary_function(FT::ArcTan, detail::make_expression_ptr(u));
        auto x_atan = SymbolicExpr::multiply(detail::make_expression_ptr(u), make_unary_function(FT::ArcTan, detail::make_expression_ptr(u)));
        auto u_sq = SymbolicExpr::power(detail::make_expression_ptr(u), SymbolicExpr::number(2));
        auto one_plus_u_sq = SymbolicExpr::add(SymbolicExpr::number(1), u_sq);
        auto half_ln = SymbolicExpr::multiply(SymbolicExpr::ln(one_plus_u_sq), sym_rational(-1, 2));
        auto res = *SymbolicExpr::add(x_atan, half_ln);
        table.add_entry(Category::InverseTrig, IntegrationEntry(
            "arctan(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
}

void load_hyperbolic_rules(IntegrationTable& table) {
    {
        auto u = wildcard("_u");
        auto pat = *make_unary_function(FT::Sinh, detail::make_expression_ptr(u));
        auto res = *make_unary_function(FT::Cosh, detail::make_expression_ptr(u));
        table.add_entry(Category::Hyperbolic, IntegrationEntry(
            "sinh(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
    {
        auto u = wildcard("_u");
        auto pat = *make_unary_function(FT::Cosh, detail::make_expression_ptr(u));
        auto res = *make_unary_function(FT::Sinh, detail::make_expression_ptr(u));
        table.add_entry(Category::Hyperbolic, IntegrationEntry(
            "cosh(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
    {
        auto u = wildcard("_u");
        auto pat = *make_unary_function(FT::Tanh, detail::make_expression_ptr(u));
        auto res = *SymbolicExpr::ln(make_unary_function(FT::Cosh, detail::make_expression_ptr(u)));
        table.add_entry(Category::Hyperbolic, IntegrationEntry(
            "tanh(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
    {
        auto u = wildcard("_u");
        auto sinh_inv = SymbolicExpr::power(make_unary_function(FT::Sinh, detail::make_expression_ptr(u)), SymbolicExpr::number(-1));
        auto pat = *SymbolicExpr::multiply(make_unary_function(FT::Cosh, detail::make_expression_ptr(u)), sinh_inv);
        auto res = *SymbolicExpr::ln(make_unary_function(FT::Sinh, detail::make_expression_ptr(u)));
        table.add_entry(Category::Hyperbolic, IntegrationEntry(
            "coth(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(make_unary_function(FT::Cosh, detail::make_expression_ptr(u)), SymbolicExpr::number(-1));
        auto res = *make_unary_function(FT::ArcTan, make_unary_function(FT::Sinh, detail::make_expression_ptr(u)));
        table.add_entry(Category::Hyperbolic, IntegrationEntry(
            "sech(x)", pat, res, {"_u"}, u_is_var("_u"), 50));
    }
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(make_unary_function(FT::Sinh, detail::make_expression_ptr(u)), SymbolicExpr::number(-1));
        auto u_half = SymbolicExpr::multiply(detail::make_expression_ptr(u), sym_rational(1, 2));
        auto res = *SymbolicExpr::ln(make_unary_function(FT::Tanh, u_half));
        table.add_entry(Category::Hyperbolic, IntegrationEntry(
            "csch(x)", pat, res, {"_u"}, u_is_var("_u"), 50));
    }
}

}
