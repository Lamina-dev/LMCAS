#include "internal/integration_support.hpp"

namespace lamina {

const std::vector<IntegrationEntry> IntegrationTable::empty_entries_;

IntegrationTable::IntegrationTable() {
    load_defaults();
}

void IntegrationTable::add_entry(Category cat, const IntegrationEntry& entry) {
    entries_[static_cast<int>(cat)].push_back(entry);

    auto& vec = entries_[static_cast<int>(cat)];
    std::sort(vec.begin(), vec.end(), [](const IntegrationEntry& a, const IntegrationEntry& b) {
        return a.priority < b.priority;
    });
}

void IntegrationTable::clear_category(Category cat) {
    entries_[static_cast<int>(cat)].clear();
}

const std::vector<IntegrationEntry>& IntegrationTable::get_entries(Category cat) const {
    auto it = entries_.find(static_cast<int>(cat));
    if (it == entries_.end()) return empty_entries_;
    return it->second;
}

std::vector<const IntegrationEntry*> IntegrationTable::get_all_sorted() const {
    std::vector<const IntegrationEntry*> all;
    for (const auto& [cat, vec] : entries_) {
        for (const auto& entry : vec) {
            all.push_back(&entry);
        }
    }
    std::sort(all.begin(), all.end(), [](const IntegrationEntry* a, const IntegrationEntry* b) {
        return a->priority < b->priority;
    });
    return all;
}

void IntegrationTable::load_defaults() {

    using FT = FunctionNode::FuncType;

    auto make_fn = [](FT t, const std::shared_ptr<SymbolicExpr>& arg) {
        return lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(
                t, std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(arg)}));
    };

    // Condition: wildcard `_u` is bound to the integration variable itself.
    auto u_is_var = [](const std::string& wc) {
        return [wc](const MatchMap& m, const std::string& var) -> bool {
            auto it = m.find(wc);
            if (it == m.end()) return false;
            auto v = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(it->second));
            return v && v->name() == var;
        };
    };

    // Condition: `_u` is the integration variable AND `_a` does not depend on it.
    auto u_is_var_a_indep = [](const std::string& u_wc, const std::string& a_wc) {
        return [u_wc, a_wc](const MatchMap& m, const std::string& var) -> bool {
            auto it_u = m.find(u_wc);
            if (it_u == m.end()) return false;
            auto v = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(it_u->second));
            if (!v || v->name() != var) return false;
            auto it_a = m.find(a_wc);
            if (it_a == m.end()) return false;
            return !expression_depends_on_variable(lamina::detail::node(it_a->second), var);
        };
    };


    // exp(a*x) -> exp(a*x)/a
    {
        auto a = wildcard("_a");
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::exp(SymbolicExpr::multiply(make_expr_ptr(a), make_expr_ptr(u)));
        auto res = *SymbolicExpr::multiply(
            SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(-1)),
            SymbolicExpr::exp(SymbolicExpr::multiply(make_expr_ptr(a), make_expr_ptr(u))));
        add_entry(Category::Exponential, IntegrationEntry(
            "exp(a*x)", pat, res, {"_a", "_u"},
            u_is_var_a_indep("_u", "_a"), 50));
    }

    // x^2 * exp(x) -> (x^2 - 2x + 2) * exp(x)
    {
        auto u = wildcard("_u");
        auto u_sq = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto pat = *SymbolicExpr::multiply(u_sq, SymbolicExpr::exp(make_expr_ptr(u)));

        auto u_sq_r = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto two_u = SymbolicExpr::multiply(SymbolicExpr::number(-2), make_expr_ptr(u));
        auto poly = SymbolicExpr::add(SymbolicExpr::add(u_sq_r, two_u), SymbolicExpr::number(2));
        auto res = *SymbolicExpr::multiply(poly, SymbolicExpr::exp(make_expr_ptr(u)));
        add_entry(Category::Exponential, IntegrationEntry(
            "x^2*exp(x)", pat, res, {"_u"}, u_is_var("_u"), 40));
    }

    // x * exp(x) -> (x - 1) * exp(x)
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::multiply(make_expr_ptr(u), SymbolicExpr::exp(make_expr_ptr(u)));

        auto u_minus_1 = SymbolicExpr::add(make_expr_ptr(u), SymbolicExpr::number(-1));
        auto res = *SymbolicExpr::multiply(u_minus_1, SymbolicExpr::exp(make_expr_ptr(u)));
        add_entry(Category::Exponential, IntegrationEntry(
            "x*exp(x)", pat, res, {"_u"}, u_is_var("_u"), 50));
    }

    // exp(x) -> exp(x)
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::exp(make_expr_ptr(u));
        auto res = *SymbolicExpr::exp(make_expr_ptr(u));
        add_entry(Category::Exponential, IntegrationEntry(
            "exp(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }


    // sec(x)^2 -> tan(x)
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(make_fn(FT::Sec, make_expr_ptr(u)), SymbolicExpr::number(2));
        auto res = *make_fn(FT::Tan, make_expr_ptr(u));
        add_entry(Category::Trigonometric, IntegrationEntry(
            "sec(x)^2", pat, res, {"_u"}, u_is_var("_u"), 30));
    }

    // csc(x)^2 -> -cot(x)
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(make_fn(FT::Csc, make_expr_ptr(u)), SymbolicExpr::number(2));
        auto res = *SymbolicExpr::multiply(SymbolicExpr::number(-1), make_fn(FT::Cot, make_expr_ptr(u)));
        add_entry(Category::Trigonometric, IntegrationEntry(
            "csc(x)^2", pat, res, {"_u"}, u_is_var("_u"), 30));
    }

    // sec(x)*tan(x) -> sec(x)
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::multiply(make_fn(FT::Sec, make_expr_ptr(u)), make_fn(FT::Tan, make_expr_ptr(u)));
        auto res = *make_fn(FT::Sec, make_expr_ptr(u));
        add_entry(Category::Trigonometric, IntegrationEntry(
            "sec(x)*tan(x)", pat, res, {"_u"}, u_is_var("_u"), 30));
    }

    // csc(x)*cot(x) -> -csc(x)
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::multiply(make_fn(FT::Csc, make_expr_ptr(u)), make_fn(FT::Cot, make_expr_ptr(u)));
        auto res = *SymbolicExpr::multiply(SymbolicExpr::number(-1), make_fn(FT::Csc, make_expr_ptr(u)));
        add_entry(Category::Trigonometric, IntegrationEntry(
            "csc(x)*cot(x)", pat, res, {"_u"}, u_is_var("_u"), 30));
    }

    // sin(x)^2 -> x/2 - sin(2x)/4
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(SymbolicExpr::sin(make_expr_ptr(u)), SymbolicExpr::number(2));

        auto u_half = SymbolicExpr::multiply(make_expr_ptr(u), sym_rational(1, 2));
        auto two_u = SymbolicExpr::multiply(SymbolicExpr::number(2), make_expr_ptr(u));
        auto sin_2u_over_4 = SymbolicExpr::multiply(SymbolicExpr::sin(two_u), sym_rational(-1, 4));
        auto res = *SymbolicExpr::add(u_half, sin_2u_over_4);
        add_entry(Category::Trigonometric, IntegrationEntry(
            "sin(x)^2", pat, res, {"_u"}, u_is_var("_u"), 35));
    }

    // cos(x)^2 -> x/2 + sin(2x)/4
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(SymbolicExpr::cos(make_expr_ptr(u)), SymbolicExpr::number(2));

        auto u_half = SymbolicExpr::multiply(make_expr_ptr(u), sym_rational(1, 2));
        auto two_u = SymbolicExpr::multiply(SymbolicExpr::number(2), make_expr_ptr(u));
        auto sin_2u_over_4 = SymbolicExpr::multiply(SymbolicExpr::sin(two_u), sym_rational(1, 4));
        auto res = *SymbolicExpr::add(u_half, sin_2u_over_4);
        add_entry(Category::Trigonometric, IntegrationEntry(
            "cos(x)^2", pat, res, {"_u"}, u_is_var("_u"), 35));
    }

    // tan(x)^2 -> tan(x) - x
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(SymbolicExpr::tan(make_expr_ptr(u)), SymbolicExpr::number(2));
        auto neg_u = SymbolicExpr::multiply(SymbolicExpr::number(-1), make_expr_ptr(u));
        auto res = *SymbolicExpr::add(SymbolicExpr::tan(make_expr_ptr(u)), neg_u);
        add_entry(Category::Trigonometric, IntegrationEntry(
            "tan(x)^2", pat, res, {"_u"}, u_is_var("_u"), 35));
    }

    // sin(x) -> -cos(x)
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::sin(make_expr_ptr(u));
        auto res = *SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::cos(make_expr_ptr(u)));
        add_entry(Category::Trigonometric, IntegrationEntry(
            "sin(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }

    // cos(x) -> sin(x)
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::cos(make_expr_ptr(u));
        auto res = *SymbolicExpr::sin(make_expr_ptr(u));
        add_entry(Category::Trigonometric, IntegrationEntry(
            "cos(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }

    // tan(x) -> -ln(cos(x))
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::tan(make_expr_ptr(u));
        auto res = *SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::ln(SymbolicExpr::cos(make_expr_ptr(u))));
        add_entry(Category::Trigonometric, IntegrationEntry(
            "tan(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }

    // cot(x) -> ln(sin(x))
    {
        auto u = wildcard("_u");
        auto pat = *make_fn(FT::Cot, make_expr_ptr(u));
        auto res = *SymbolicExpr::ln(SymbolicExpr::sin(make_expr_ptr(u)));
        add_entry(Category::Trigonometric, IntegrationEntry(
            "cot(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }

    // sec(x) -> ln(sec(x) + tan(x))
    {
        auto u = wildcard("_u");
        auto pat = *make_fn(FT::Sec, make_expr_ptr(u));
        auto inner = SymbolicExpr::add(make_fn(FT::Sec, make_expr_ptr(u)), make_fn(FT::Tan, make_expr_ptr(u)));
        auto res = *SymbolicExpr::ln(inner);
        add_entry(Category::Trigonometric, IntegrationEntry(
            "sec(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }

    // csc(x) -> -ln(csc(x) + cot(x))
    {
        auto u = wildcard("_u");
        auto pat = *make_fn(FT::Csc, make_expr_ptr(u));
        auto inner = SymbolicExpr::add(make_fn(FT::Csc, make_expr_ptr(u)), make_fn(FT::Cot, make_expr_ptr(u)));
        auto res = *SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::ln(inner));
        add_entry(Category::Trigonometric, IntegrationEntry(
            "csc(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }


    // arcsin(x) -> x*arcsin(x) + sqrt(1 - x^2)
    {
        auto u = wildcard("_u");
        auto pat = *make_fn(FT::ArcSin, make_expr_ptr(u));
        auto x_asin = SymbolicExpr::multiply(make_expr_ptr(u), make_fn(FT::ArcSin, make_expr_ptr(u)));
        auto u_sq = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto neg_u_sq = SymbolicExpr::multiply(SymbolicExpr::number(-1), u_sq);
        auto one_minus_u_sq = SymbolicExpr::add(SymbolicExpr::number(1), neg_u_sq);
        auto root = SymbolicExpr::sqrt(one_minus_u_sq);
        auto res = *SymbolicExpr::add(x_asin, root);
        add_entry(Category::InverseTrig, IntegrationEntry(
            "arcsin(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }

    // arccos(x) -> x*arccos(x) - sqrt(1 - x^2)
    {
        auto u = wildcard("_u");
        auto pat = *make_fn(FT::ArcCos, make_expr_ptr(u));
        auto x_acos = SymbolicExpr::multiply(make_expr_ptr(u), make_fn(FT::ArcCos, make_expr_ptr(u)));
        auto u_sq = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto neg_u_sq = SymbolicExpr::multiply(SymbolicExpr::number(-1), u_sq);
        auto one_minus_u_sq = SymbolicExpr::add(SymbolicExpr::number(1), neg_u_sq);
        auto neg_root = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::sqrt(one_minus_u_sq));
        auto res = *SymbolicExpr::add(x_acos, neg_root);
        add_entry(Category::InverseTrig, IntegrationEntry(
            "arccos(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }

    // arctan(x) -> x*arctan(x) - ln(1 + x^2) / 2
    {
        auto u = wildcard("_u");
        auto pat = *make_fn(FT::ArcTan, make_expr_ptr(u));
        auto x_atan = SymbolicExpr::multiply(make_expr_ptr(u), make_fn(FT::ArcTan, make_expr_ptr(u)));
        auto u_sq = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto one_plus_u_sq = SymbolicExpr::add(SymbolicExpr::number(1), u_sq);
        auto half_ln = SymbolicExpr::multiply(SymbolicExpr::ln(one_plus_u_sq), sym_rational(-1, 2));
        auto res = *SymbolicExpr::add(x_atan, half_ln);
        add_entry(Category::InverseTrig, IntegrationEntry(
            "arctan(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }


    // sinh(x) -> cosh(x)
    {
        auto u = wildcard("_u");
        auto pat = *make_fn(FT::Sinh, make_expr_ptr(u));
        auto res = *make_fn(FT::Cosh, make_expr_ptr(u));
        add_entry(Category::Hyperbolic, IntegrationEntry(
            "sinh(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }

    // cosh(x) -> sinh(x)
    {
        auto u = wildcard("_u");
        auto pat = *make_fn(FT::Cosh, make_expr_ptr(u));
        auto res = *make_fn(FT::Sinh, make_expr_ptr(u));
        add_entry(Category::Hyperbolic, IntegrationEntry(
            "cosh(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }

    // tanh(x) -> ln(cosh(x))
    {
        auto u = wildcard("_u");
        auto pat = *make_fn(FT::Tanh, make_expr_ptr(u));
        auto res = *SymbolicExpr::ln(make_fn(FT::Cosh, make_expr_ptr(u)));
        add_entry(Category::Hyperbolic, IntegrationEntry(
            "tanh(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }

    // coth(x) (= cosh(x)/sinh(x)) -> ln(sinh(x))
    {
        auto u = wildcard("_u");
        auto sinh_inv = SymbolicExpr::power(make_fn(FT::Sinh, make_expr_ptr(u)), SymbolicExpr::number(-1));
        auto pat = *SymbolicExpr::multiply(make_fn(FT::Cosh, make_expr_ptr(u)), sinh_inv);
        auto res = *SymbolicExpr::ln(make_fn(FT::Sinh, make_expr_ptr(u)));
        add_entry(Category::Hyperbolic, IntegrationEntry(
            "coth(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }

    // sech(x) (= 1/cosh(x)) -> arctan(sinh(x))
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(make_fn(FT::Cosh, make_expr_ptr(u)), SymbolicExpr::number(-1));
        auto res = *make_fn(FT::ArcTan, make_fn(FT::Sinh, make_expr_ptr(u)));
        add_entry(Category::Hyperbolic, IntegrationEntry(
            "sech(x)", pat, res, {"_u"}, u_is_var("_u"), 50));
    }

    // csch(x) (= 1/sinh(x)) -> ln(tanh(x/2))
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(make_fn(FT::Sinh, make_expr_ptr(u)), SymbolicExpr::number(-1));
        auto u_half = SymbolicExpr::multiply(make_expr_ptr(u), sym_rational(1, 2));
        auto res = *SymbolicExpr::ln(make_fn(FT::Tanh, u_half));
        add_entry(Category::Hyperbolic, IntegrationEntry(
            "csch(x)", pat, res, {"_u"}, u_is_var("_u"), 50));
    }


    // 1/sqrt(1 - x^2) -> arcsin(x)
    {
        auto u = wildcard("_u");
        auto u_sq = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto neg_u_sq = SymbolicExpr::multiply(SymbolicExpr::number(-1), u_sq);
        auto one_minus_u_sq = SymbolicExpr::add(SymbolicExpr::number(1), neg_u_sq);
        auto pat = *SymbolicExpr::power(SymbolicExpr::sqrt(one_minus_u_sq), SymbolicExpr::number(-1));
        auto res = *make_fn(FT::ArcSin, make_expr_ptr(u));
        add_entry(Category::Algebraic, IntegrationEntry(
            "1/sqrt(1-x^2)", pat, res, {"_u"}, u_is_var("_u"), 40));
    }

    /// |x| -> x*|x|/2   (since d/dx[x|x|/2] = |x|)
    {
        auto u = wildcard("_u");
        auto pat = *make_fn(FT::Abs, make_expr_ptr(u));
        auto res = *SymbolicExpr::multiply(
            SymbolicExpr::number(Rational(1, 2)),
            SymbolicExpr::multiply(make_expr_ptr(u), make_fn(FT::Abs, make_expr_ptr(u))));
        add_entry(Category::Algebraic, IntegrationEntry(
            "|x|", pat, res, {"_u"}, u_is_var("_u"), 40));
    }

    // 1/(1 + x^2) -> arctan(x)
    {
        auto u = wildcard("_u");
        auto u_sq = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto one_plus_u_sq = SymbolicExpr::add(SymbolicExpr::number(1), u_sq);
        auto pat = *SymbolicExpr::power(one_plus_u_sq, SymbolicExpr::number(-1));
        auto res = *make_fn(FT::ArcTan, make_expr_ptr(u));
        add_entry(Category::Algebraic, IntegrationEntry(
            "1/(1+x^2)", pat, res, {"_u"}, u_is_var("_u"), 40));
    }

    // 1/sqrt(x^2 + a^2) -> ln(x + sqrt(x^2 + a^2))
    {
        auto u = wildcard("_u");
        auto a = wildcard("_a");
        auto u_sq = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto a_sq = SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(2));
        auto sum = SymbolicExpr::add(u_sq, a_sq);
        auto pat = *SymbolicExpr::power(SymbolicExpr::sqrt(sum), SymbolicExpr::number(-1));

        auto u_sq_r = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto a_sq_r = SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(2));
        auto sum_r = SymbolicExpr::add(u_sq_r, a_sq_r);
        auto inner = SymbolicExpr::add(make_expr_ptr(u), SymbolicExpr::sqrt(sum_r));
        auto res = *SymbolicExpr::ln(inner);
        add_entry(Category::Algebraic, IntegrationEntry(
            "1/sqrt(x^2+a^2)", pat, res, {"_u", "_a"},
            u_is_var_a_indep("_u", "_a"), 50));
    }

    // 1/sqrt(x^2 - a^2) -> ln(x + sqrt(x^2 - a^2))
    {
        auto u = wildcard("_u");
        auto a = wildcard("_a");
        auto u_sq = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto a_sq = SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(2));
        auto neg_a_sq = SymbolicExpr::multiply(SymbolicExpr::number(-1), a_sq);
        auto diff = SymbolicExpr::add(u_sq, neg_a_sq);
        auto pat = *SymbolicExpr::power(SymbolicExpr::sqrt(diff), SymbolicExpr::number(-1));

        auto u_sq_r = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto a_sq_r = SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(2));
        auto neg_a_sq_r = SymbolicExpr::multiply(SymbolicExpr::number(-1), a_sq_r);
        auto diff_r = SymbolicExpr::add(u_sq_r, neg_a_sq_r);
        auto inner = SymbolicExpr::add(make_expr_ptr(u), SymbolicExpr::sqrt(diff_r));
        auto res = *SymbolicExpr::ln(inner);
        add_entry(Category::Algebraic, IntegrationEntry(
            "1/sqrt(x^2-a^2)", pat, res, {"_u", "_a"},
            u_is_var_a_indep("_u", "_a"), 50));
    }

    // 1/(x^2 + a^2) -> (1/a) * arctan(x/a)
    {
        auto u = wildcard("_u");
        auto a = wildcard("_a");
        auto u_sq = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto a_sq = SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(2));
        auto sum = SymbolicExpr::add(u_sq, a_sq);
        auto pat = *SymbolicExpr::power(sum, SymbolicExpr::number(-1));

        auto inv_a = SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(-1));
        auto u_over_a = SymbolicExpr::multiply(make_expr_ptr(u), SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(-1)));
        auto res = *SymbolicExpr::multiply(inv_a, make_fn(FT::ArcTan, u_over_a));
        add_entry(Category::Algebraic, IntegrationEntry(
            "1/(x^2+a^2)", pat, res, {"_u", "_a"},
            u_is_var_a_indep("_u", "_a"), 50));
    }

    // 1/(x^2 - a^2) -> (1/(2a)) * ln((x - a)/(x + a))
    {
        auto u = wildcard("_u");
        auto a = wildcard("_a");
        auto u_sq = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto a_sq = SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(2));
        auto neg_a_sq = SymbolicExpr::multiply(SymbolicExpr::number(-1), a_sq);
        auto diff = SymbolicExpr::add(u_sq, neg_a_sq);
        auto pat = *SymbolicExpr::power(diff, SymbolicExpr::number(-1));

        auto two_a = SymbolicExpr::multiply(SymbolicExpr::number(2), make_expr_ptr(a));
        auto inv_2a = SymbolicExpr::power(two_a, SymbolicExpr::number(-1));
        auto neg_a = SymbolicExpr::multiply(SymbolicExpr::number(-1), make_expr_ptr(a));
        auto u_minus_a = SymbolicExpr::add(make_expr_ptr(u), neg_a);
        auto u_plus_a = SymbolicExpr::add(make_expr_ptr(u), make_expr_ptr(a));
        auto frac = SymbolicExpr::divide(u_minus_a, u_plus_a);
        auto res = *SymbolicExpr::multiply(inv_2a, SymbolicExpr::ln(frac));
        add_entry(Category::Algebraic, IntegrationEntry(
            "1/(x^2-a^2)", pat, res, {"_u", "_a"},
            u_is_var_a_indep("_u", "_a"), 50));
    }


    // 1/x  (= x^(-1))  -> ln(x)
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(-1));
        auto res = *SymbolicExpr::ln(make_expr_ptr(u));
        add_entry(Category::Polynomial, IntegrationEntry(
            "1/x", pat, res, {"_u"}, u_is_var("_u"), 50));
    }

    // x^n  (n a constant != -1)  -> x^(n+1) / (n+1)
    {
        auto u = wildcard("_u");
        auto n = wildcard("_n");
        auto pat = *SymbolicExpr::power(make_expr_ptr(u), make_expr_ptr(n));
        auto n_plus_1 = SymbolicExpr::add(make_expr_ptr(n), SymbolicExpr::number(1));
        auto res = *SymbolicExpr::divide(
            SymbolicExpr::power(make_expr_ptr(u), n_plus_1),
            n_plus_1);
        add_entry(Category::Polynomial, IntegrationEntry(
            "x^n", pat, res, {"_u", "_n"},
            [](const MatchMap& m, const std::string& var) {
                auto it_u = m.find("_u");
                if (it_u == m.end()) return false;
                auto v = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(it_u->second));
                if (!v || v->name() != var) return false;
                auto it_n = m.find("_n");
                if (it_n == m.end()) return false;
                if (expression_depends_on_variable(lamina::detail::node(it_n->second), var)) return false;

                // Reject n = -1 (handled by the dedicated 1/x rule).
                auto n_simp = it_n->second.simplify();
                if (!n_simp) return false;
                auto neg_one = SymbolicExpr::number(-1);
                auto diff = SymbolicExpr::add(n_simp, neg_one)->simplify();
                if (diff && diff->is_zero()) return false;
                return true;
            }, 80));
    }

    // x  -> x^2 / 2  (anchor entry for the bare integration variable)
    {
        auto u = wildcard("_u");
        auto pat = lamina::detail::expression_from_node(lamina::detail::node(u));
        auto res = *SymbolicExpr::multiply(
            SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2)),
            sym_rational(1, 2));
        add_entry(Category::Polynomial, IntegrationEntry(
            "x", pat, res, {"_u"}, u_is_var("_u"), 90));
    }


    // ln(x) -> x*ln(x) - x
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::ln(make_expr_ptr(u));
        auto x_ln_x = SymbolicExpr::multiply(make_expr_ptr(u), SymbolicExpr::ln(make_expr_ptr(u)));
        auto res = *sym_sub(*x_ln_x, u);
        add_entry(Category::Logarithmic, IntegrationEntry(
            "ln(x)", pat, res, {"_u"}, u_is_var("_u"), 60));
    }

    // ln(x)/x -> (ln(x))^2 / 2
    {
        auto u = wildcard("_u");
        auto inv_u = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(-1));
        auto pat = *SymbolicExpr::multiply(SymbolicExpr::ln(make_expr_ptr(u)), inv_u);
        auto ln_sq = SymbolicExpr::power(SymbolicExpr::ln(make_expr_ptr(u)), SymbolicExpr::number(2));
        auto res = *SymbolicExpr::multiply(ln_sq, sym_rational(1, 2));
        add_entry(Category::Logarithmic, IntegrationEntry(
            "ln(x)/x", pat, res, {"_u"}, u_is_var("_u"), 50));
    }


    // a^x -> a^x / ln(a)
    {
        auto a = wildcard("_a");
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::power(make_expr_ptr(a), make_expr_ptr(u));
        auto inv_ln_a = SymbolicExpr::power(SymbolicExpr::ln(make_expr_ptr(a)), SymbolicExpr::number(-1));
        auto a_pow_u = SymbolicExpr::power(make_expr_ptr(a), make_expr_ptr(u));
        auto res = *SymbolicExpr::multiply(a_pow_u, inv_ln_a);
        add_entry(Category::Exponential, IntegrationEntry(
            "a^x", pat, res, {"_a", "_u"}, u_is_var_a_indep("_u", "_a"), 70));
    }

    // 1/sqrt(a^2 - x^2) -> arcsin(x/a)
    {
        auto u = wildcard("_u");
        auto a = wildcard("_a");
        auto a_sq = SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(2));
        auto u_sq = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto neg_u_sq = SymbolicExpr::multiply(SymbolicExpr::number(-1), u_sq);
        auto a_sq_minus_u_sq = SymbolicExpr::add(a_sq, neg_u_sq);
        auto pat = *SymbolicExpr::power(SymbolicExpr::sqrt(a_sq_minus_u_sq), SymbolicExpr::number(-1));

        auto u_over_a = SymbolicExpr::multiply(make_expr_ptr(u),
            SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(-1)));
        auto res = *make_fn(FT::ArcSin, u_over_a);
        add_entry(Category::Algebraic, IntegrationEntry(
            "1/sqrt(a^2-x^2)", pat, res, {"_u", "_a"},
            u_is_var_a_indep("_u", "_a"), 50));
    }

    // sin(a*x) -> -cos(a*x)/a
    {
        auto a = wildcard("_a");
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::sin(SymbolicExpr::multiply(make_expr_ptr(a), make_expr_ptr(u)));
        auto inv_a = SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(-1));
        auto neg_inv_a = SymbolicExpr::multiply(SymbolicExpr::number(-1), inv_a);
        auto cos_ax = SymbolicExpr::cos(SymbolicExpr::multiply(make_expr_ptr(a), make_expr_ptr(u)));
        auto res = *SymbolicExpr::multiply(neg_inv_a, cos_ax);
        add_entry(Category::Trigonometric, IntegrationEntry(
            "sin(a*x)", pat, res, {"_a", "_u"},
            u_is_var_a_indep("_u", "_a"), 55));
    }

    // cos(a*x) -> sin(a*x)/a
    {
        auto a = wildcard("_a");
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::cos(SymbolicExpr::multiply(make_expr_ptr(a), make_expr_ptr(u)));
        auto inv_a = SymbolicExpr::power(make_expr_ptr(a), SymbolicExpr::number(-1));
        auto sin_ax = SymbolicExpr::sin(SymbolicExpr::multiply(make_expr_ptr(a), make_expr_ptr(u)));
        auto res = *SymbolicExpr::multiply(inv_a, sin_ax);
        add_entry(Category::Trigonometric, IntegrationEntry(
            "cos(a*x)", pat, res, {"_a", "_u"},
            u_is_var_a_indep("_u", "_a"), 55));
    }

    // x*cos(x^2) -> sin(x^2)/2
    {
        auto u = wildcard("_u");
        auto u_sq = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto pat = *SymbolicExpr::multiply(make_expr_ptr(u), SymbolicExpr::cos(u_sq));

        auto u_sq_r = SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2));
        auto res = *SymbolicExpr::multiply(sym_rational(1, 2), SymbolicExpr::sin(u_sq_r));
        add_entry(Category::Trigonometric, IntegrationEntry(
            "x*cos(x^2)", pat, res, {"_u"}, u_is_var("_u"), 35));
    }

    // exp(x)*sin(x) -> exp(x)*(sin(x)-cos(x))/2
    {
        auto u = wildcard("_u");
        auto pat = *SymbolicExpr::multiply(
            SymbolicExpr::exp(make_expr_ptr(u)),
            SymbolicExpr::sin(make_expr_ptr(u)));

        auto sin_u = SymbolicExpr::sin(make_expr_ptr(u));
        auto neg_cos_u = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::cos(make_expr_ptr(u)));
        auto sin_minus_cos = SymbolicExpr::add(sin_u, neg_cos_u);
        auto exp_u = SymbolicExpr::exp(make_expr_ptr(u));
        auto res = *SymbolicExpr::multiply(
            sym_rational(1, 2),
            SymbolicExpr::multiply(exp_u, sin_minus_cos));
        add_entry(Category::Exponential, IntegrationEntry(
            "exp(x)*sin(x)", pat, res, {"_u"}, u_is_var("_u"), 35));
    }
}

} // namespace lamina
