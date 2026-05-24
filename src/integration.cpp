#include "integration.hpp"
#include "symbolic_ast.hpp"
#include "polynomial.hpp"
#include "poly_utils.hpp"
#include "solve_polynomial.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <variant>

namespace lamina {

static std::shared_ptr<SymbolicExpr> make_expr_ptr(const SymbolicExpr& e) {
    return std::make_shared<SymbolicExpr>(e);
}

static bool valid_dependency(const SymbolicExpr& expr, const std::string& var) {
    auto diff = expr.differentiate(var);
    if (!diff) return false;
    auto simp_diff = diff->simplify();
    return !simp_diff->is_zero();
}

static std::shared_ptr<SymbolicExpr> sym_sub(const SymbolicExpr& a, const SymbolicExpr& b) {
    auto neg_b = SymbolicExpr::multiply(SymbolicExpr::number(-1), std::make_shared<SymbolicExpr>(b));
    return SymbolicExpr::add(std::make_shared<SymbolicExpr>(a), neg_b);
}

static std::shared_ptr<SymbolicExpr> sym_rational(long long num, long long den) {
    return SymbolicExpr::number(Rational(BigInt(num), BigInt(den)));
}

static std::shared_ptr<SymbolicExpr> make_arctan(const std::shared_ptr<SymbolicExpr>& op) {
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(
            FunctionNode::FuncType::ArcTan,
            std::vector<std::shared_ptr<SymbolicNode>>{op->root}));
}

static bool has_integral_node_check(const std::shared_ptr<SymbolicNode>& node) {
    if (!node) return false;
    if (auto fn = std::dynamic_pointer_cast<FunctionNode>(node)) {
        if (fn->type == FunctionNode::FuncType::Calculus_Integral) return true;
        for (auto& arg : fn->arguments)
            if (has_integral_node_check(arg)) return true;
    } else if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        for (auto& op : add->operands) if (has_integral_node_check(op)) return true;
    } else if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        for (auto& op : mul->operands) if (has_integral_node_check(op)) return true;
    } else if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (has_integral_node_check(pow->base)) return true;
        if (has_integral_node_check(pow->exponent)) return true;
    }
    return false;
}

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

    // Helper: build a single-argument FunctionNode wrapped in a SymbolicExpr ptr.
    auto make_fn = [](FT t, const std::shared_ptr<SymbolicExpr>& arg) {
        return std::make_shared<SymbolicExpr>(
            std::make_shared<FunctionNode>(
                t, std::vector<std::shared_ptr<SymbolicNode>>{arg->root}));
    };

    // Condition: wildcard `_u` is bound to the integration variable itself.
    auto u_is_var = [](const std::string& wc) {
        return [wc](const MatchMap& m, const std::string& var) -> bool {
            auto it = m.find(wc);
            if (it == m.end()) return false;
            auto v = std::dynamic_pointer_cast<VariableNode>(it->second.root);
            return v && v->name == var;
        };
    };

    // Condition: `_u` is the integration variable AND `_a` does not depend on it.
    auto u_is_var_a_indep = [](const std::string& u_wc, const std::string& a_wc) {
        return [u_wc, a_wc](const MatchMap& m, const std::string& var) -> bool {
            auto it_u = m.find(u_wc);
            if (it_u == m.end()) return false;
            auto v = std::dynamic_pointer_cast<VariableNode>(it_u->second.root);
            if (!v || v->name != var) return false;
            auto it_a = m.find(a_wc);
            if (it_a == m.end()) return false;
            return !depends_on_var(it_a->second.root, var);
        };
    };

    // ---------------------------------------------------------------
    // Exponential
    // ---------------------------------------------------------------

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

    // ---------------------------------------------------------------
    // Trigonometric (specific composite patterns first - lower priority)
    // ---------------------------------------------------------------

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

    // ---------------------------------------------------------------
    // Inverse Trigonometric
    // ---------------------------------------------------------------

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

    // ---------------------------------------------------------------
    // Hyperbolic
    // ---------------------------------------------------------------

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

    // ---------------------------------------------------------------
    // Algebraic forms
    // ---------------------------------------------------------------

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

    // ---------------------------------------------------------------
    // Polynomial
    // ---------------------------------------------------------------

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
                auto v = std::dynamic_pointer_cast<VariableNode>(it_u->second.root);
                if (!v || v->name != var) return false;
                auto it_n = m.find("_n");
                if (it_n == m.end()) return false;
                if (depends_on_var(it_n->second.root, var)) return false;

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
        auto pat = SymbolicExpr(u.root);
        auto res = *SymbolicExpr::multiply(
            SymbolicExpr::power(make_expr_ptr(u), SymbolicExpr::number(2)),
            sym_rational(1, 2));
        add_entry(Category::Polynomial, IntegrationEntry(
            "x", pat, res, {"_u"}, u_is_var("_u"), 90));
    }

    // ---------------------------------------------------------------
    // Logarithmic
    // ---------------------------------------------------------------

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

    // ---------------------------------------------------------------
    // Additional standard-table forms
    // ---------------------------------------------------------------

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
}

std::shared_ptr<SymbolicExpr> TableLookupStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth) {

    auto all_entries = ctx.table().get_all_sorted();

    for (const auto* entry : all_entries) {
        MatchMap bindings;
        if (Matcher::match(entry->pattern, expr, entry->wildcards, bindings)) {

            // The matcher allows partial matches on commutative operations
            // (Add/Multiply): unmatched operands are bound to __Add_REST__ /
            // __Mul_REST__ for use by rewrite rules. For integration table
            // lookup that behavior is wrong - a pattern like (1 + _u^2)^-1
            // would otherwise match (1 + x^2 + x + x^3)^-1 with _u=x and the
            // rest x + x^3 silently dropped, producing arctan(x) for an
            // integrand that has nothing to do with arctan. Require an exact
            // match (no leftover) so the result genuinely equals the
            // integrand under this binding.
            if (bindings.find("__Add_REST__") != bindings.end() ||
                bindings.find("__Mul_REST__") != bindings.end()) {
                continue;
            }

            if (entry->condition && !entry->condition(bindings, var)) {
                continue;
            }

            SymbolicExpr result = Matcher::replace(entry->result, bindings, false);
            auto simplified = result.simplify();
            if (simplified && !has_integral_node_check(simplified->root)) {
                return simplified;
            }
        }
    }
    return nullptr;
}

std::shared_ptr<SymbolicExpr> PowerRuleStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth) {

    if (auto v_node = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
        if (v_node->name == var) {
            return SymbolicExpr::multiply(
                SymbolicExpr::power(make_expr_ptr(expr), SymbolicExpr::number(2)),
                sym_rational(1, 2));
        }
    }

    if (auto p_node = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        SymbolicExpr base(p_node->base);
        SymbolicExpr exp_expr(p_node->exponent);

        if (auto b_var = std::dynamic_pointer_cast<VariableNode>(base.root)) {
            if (b_var->name == var && !valid_dependency(exp_expr, var)) {
                auto n_plus_1 = SymbolicExpr::add(make_expr_ptr(exp_expr), SymbolicExpr::number(1))->simplify();
                if (n_plus_1->is_zero()) {

                    return SymbolicExpr::ln(make_expr_ptr(base));
                }
                return SymbolicExpr::divide(
                    SymbolicExpr::power(make_expr_ptr(base), n_plus_1),
                    n_plus_1);
            }
        }
    }
    return nullptr;
}

std::shared_ptr<SymbolicExpr> SubstitutionStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth) {

    std::vector<std::shared_ptr<SymbolicNode>> ops;
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        ops = mul->operands;
    } else {
        ops.push_back(expr.root);
    }

    for (size_t i = 0; i < ops.size(); ++i) {
        SymbolicExpr candidate_term(ops[i]);
        SymbolicExpr u;
        bool possible = false;

        if (auto pow = std::dynamic_pointer_cast<PowerNode>(candidate_term.root)) {
            u = SymbolicExpr(pow->base);
            possible = true;
        } else if (auto func = std::dynamic_pointer_cast<FunctionNode>(candidate_term.root)) {
            if (!func->arguments.empty()) {
                u = SymbolicExpr(func->arguments[0]);
                possible = true;
            }
        }

        if (possible && valid_dependency(u, var)) {
            auto d_ptr = u.differentiate(var);
            if (!d_ptr) continue;
            auto du = d_ptr->simplify();
            if (du->is_zero()) continue;

            auto f_u = candidate_term;
            auto term_times_du = SymbolicExpr::multiply(make_expr_ptr(f_u), make_expr_ptr(*du));
            auto ratio = SymbolicExpr::divide(make_expr_ptr(expr), term_times_du)->simplify();

            if (!valid_dependency(*ratio, var)) {
                std::shared_ptr<SymbolicExpr> prim = nullptr;

                if (auto pow = std::dynamic_pointer_cast<PowerNode>(candidate_term.root)) {
                    SymbolicExpr n(pow->exponent);
                    auto np1 = SymbolicExpr::add(make_expr_ptr(n), SymbolicExpr::number(1))->simplify();
                    if (np1->is_zero()) {
                        prim = SymbolicExpr::ln(make_expr_ptr(u));
                    } else {
                        prim = SymbolicExpr::divide(
                            SymbolicExpr::power(make_expr_ptr(u), np1), np1);
                    }
                } else if (auto func = std::dynamic_pointer_cast<FunctionNode>(candidate_term.root)) {
                    if (func->type == FunctionNode::FuncType::Cos) {
                        prim = SymbolicExpr::sin(make_expr_ptr(u));
                    } else if (func->type == FunctionNode::FuncType::Sin) {
                        prim = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::cos(make_expr_ptr(u)));
                    } else if (func->type == FunctionNode::FuncType::Exp) {
                        prim = SymbolicExpr::exp(make_expr_ptr(u));
                    }
                }

                if (prim) {
                    return SymbolicExpr::multiply(ratio, prim);
                }
            }
        }
    }
    return nullptr;
}

// ---------------------------------------------------------------
// LinearSubstitutionStrategy
// ---------------------------------------------------------------
//
// Recognises an integrand whose outermost wrapper is a single-argument
// FunctionNode (e.g. sin/cos/exp/sec(...)) or a PowerNode (e.g. (2x+1)^5,
// sec(2x+1)^2), and whose inner argument decomposes as a*var + b with a, b
// constant w.r.t. var and a != 0. It then rewrites the integrand into the
// "base form" (substituting the inner argument by a fresh dummy variable),
// looks the base form up via the existing table, and reconstructs the answer
// as (1/a) * F(a*var + b) where F is the antiderivative returned by the table.

bool LinearSubstitutionStrategy::extract_linear_arg(
    const SymbolicExpr& arg,
    const std::string& var,
    std::shared_ptr<SymbolicExpr>& a_out,
    std::shared_ptr<SymbolicExpr>& b_out) {

    // Argument must actually depend on the integration variable.
    if (!depends_on_var(arg.root, var)) return false;

    Polynomial<SymbolicPolyCoeff> poly;
    try {
        poly = symbolic_to_poly<SymbolicPolyCoeff>(make_expr_ptr(arg), var);
    } catch (...) {
        return false;
    }

    // Need to be exactly degree 1 in var: arg = a * var + b with a, b constant.
    if (poly.degree() != 1) return false;
    if (poly.coeffs.size() < 2) return false;

    auto a_expr = poly.coeffs[1].val;
    auto b_expr = poly.coeffs[0].val;
    if (!a_expr || !b_expr) return false;

    // Coefficients must not depend on the integration variable.
    if (depends_on_var(a_expr->root, var)) return false;
    if (depends_on_var(b_expr->root, var)) return false;

    auto a_simp = a_expr->simplify();
    auto b_simp = b_expr->simplify();
    if (!a_simp) a_simp = a_expr;
    if (!b_simp) b_simp = b_expr;

    // a must be demonstrably non-zero (refuse if simplification yields zero).
    if (a_simp->is_zero()) return false;

    a_out = a_simp;
    b_out = b_simp;
    return true;
}

std::shared_ptr<SymbolicExpr> LinearSubstitutionStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth) {

    using FT = FunctionNode::FuncType;

    // ---------------------------------------------------------------
    // Step 1: identify the outer wrapper and the inner argument we'll
    //         try to express as a*var + b.
    // ---------------------------------------------------------------
    enum class Wrapper { None, Function, PowerOfFunction, PowerOfLinear };
    Wrapper kind = Wrapper::None;

    std::shared_ptr<FunctionNode> fn;       // outer FunctionNode (Function case)
    std::shared_ptr<PowerNode>    pn;       // outer PowerNode (Power* cases)
    std::shared_ptr<FunctionNode> base_fn;  // inner FunctionNode of a PowerNode
    SymbolicExpr arg_expr;                  // the linear-candidate sub-expression

    if ((fn = std::dynamic_pointer_cast<FunctionNode>(expr.root))) {
        // Skip wrappers that are not "real" single-argument functions.
        switch (fn->type) {
            case FT::Calculus_Integral:
            case FT::Infinity:
            case FT::Limit:
            case FT::RootOf:
            case FT::Atan2:
            case FT::Log:
                return nullptr;
            default:
                break;
        }
        if (fn->arguments.size() != 1) return nullptr;
        arg_expr = SymbolicExpr(fn->arguments[0]);
        kind = Wrapper::Function;
    } else if ((pn = std::dynamic_pointer_cast<PowerNode>(expr.root))) {
        // Exponent must be independent of the integration variable so that
        // substituting the linear argument back is sound.
        SymbolicExpr exp_expr(pn->exponent);
        if (depends_on_var(exp_expr.root, var)) return nullptr;

        if ((base_fn = std::dynamic_pointer_cast<FunctionNode>(pn->base))) {
            switch (base_fn->type) {
                case FT::Calculus_Integral:
                case FT::Infinity:
                case FT::Limit:
                case FT::RootOf:
                case FT::Atan2:
                case FT::Log:
                    return nullptr;
                default:
                    break;
            }
            if (base_fn->arguments.size() != 1) return nullptr;
            arg_expr = SymbolicExpr(base_fn->arguments[0]);
            kind = Wrapper::PowerOfFunction;
        } else {
            // Treat the whole base as the linear-candidate (e.g. (2x+1)^n).
            arg_expr = SymbolicExpr(pn->base);
            kind = Wrapper::PowerOfLinear;
        }
    } else {
        return nullptr;
    }

    // ---------------------------------------------------------------
    // Step 2: extract a and b from arg = a*var + b.
    // ---------------------------------------------------------------
    std::shared_ptr<SymbolicExpr> a_coeff, b_coeff;
    if (!extract_linear_arg(arg_expr, var, a_coeff, b_coeff)) return nullptr;

    // Degenerate case: arg == var (a = 1, b = 0). The plain TableLookup /
    // PowerRule strategies would have handled this earlier in the chain;
    // returning nullptr lets the next strategy try.
    bool a_is_one  = a_coeff && a_coeff->is_one();
    bool b_is_zero = b_coeff && b_coeff->is_zero();
    if (a_is_one && b_is_zero) return nullptr;

    // ---------------------------------------------------------------
    // Step 3: build the "base-form" test expression with a fresh dummy
    //         variable in place of the linear argument.
    // ---------------------------------------------------------------
    // Use a name unlikely to clash with anything in user input or table entries.
    const std::string dummy_name = "__lin_sub_u__";
    auto dummy_var = SymbolicExpr::variable(dummy_name);

    std::shared_ptr<SymbolicExpr> test_expr;
    if (kind == Wrapper::Function) {
        std::vector<std::shared_ptr<SymbolicNode>> new_args{dummy_var->root};
        test_expr = std::make_shared<SymbolicExpr>(
            std::make_shared<FunctionNode>(fn->type, new_args));
    } else if (kind == Wrapper::PowerOfFunction) {
        std::vector<std::shared_ptr<SymbolicNode>> new_args{dummy_var->root};
        auto new_fn = std::make_shared<FunctionNode>(base_fn->type, new_args);
        test_expr = std::make_shared<SymbolicExpr>(
            std::make_shared<PowerNode>(new_fn, pn->exponent));
    } else { // Wrapper::PowerOfLinear
        test_expr = std::make_shared<SymbolicExpr>(
            std::make_shared<PowerNode>(dummy_var->root, pn->exponent));
    }
    if (!test_expr) return nullptr;

    // ---------------------------------------------------------------
    // Step 4: look up the base form in the table, treating dummy_name as
    //         the integration variable. We deliberately use only the
    //         table-lookup strategy here so that we don't recurse back
    //         through the full strategy chain.
    // ---------------------------------------------------------------
    TableLookupStrategy table_only;
    auto F_dummy = table_only.try_integrate(*test_expr, dummy_name, ctx, depth);
    if (!F_dummy) return nullptr;

    // ---------------------------------------------------------------
    // Step 5: substitute dummy_name -> (a*var + b) in the antiderivative
    //         and multiply by 1/a.
    // ---------------------------------------------------------------
    auto F_substituted = F_dummy->substitute(dummy_name, make_expr_ptr(arg_expr));
    if (!F_substituted) return nullptr;

    auto inv_a = SymbolicExpr::power(a_coeff, SymbolicExpr::number(-1));
    auto result = SymbolicExpr::multiply(inv_a, F_substituted);

    auto simplified = result->simplify();
    if (simplified && !has_integral_node_check(simplified->root)) {
        return simplified;
    }
    return result;
}

std::shared_ptr<SymbolicExpr> PartialFractionStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth) {

    std::shared_ptr<SymbolicExpr> den = nullptr;

    if (auto p = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        double exp_val = 0;
        bool is_inv = false;
        if (auto num_node = std::dynamic_pointer_cast<NumberNode>(p->exponent)) {
            if (std::holds_alternative<lmmc_real_t>(num_node->value))
                exp_val = std::get<lmmc_real_t>(num_node->value);
            else if (std::holds_alternative<Rational>(num_node->value))
                exp_val = std::get<Rational>(num_node->value).to_double();
            else if (std::holds_alternative<BigInt>(num_node->value))
                exp_val = std::get<BigInt>(num_node->value).to_double();
            int eq;
            lmmc_double_nearly_equal_tol(exp_val, -1.0, 1e-9, 1e-9, &eq);
            if (eq) is_inv = true;
        }
        if (is_inv) {
            den = make_expr_ptr(SymbolicExpr(p->base));
        }
    }

    if (!den) return nullptr;

    try {
        Polynomial<SymbolicPolyCoeff> Q = symbolic_to_poly<SymbolicPolyCoeff>(den, var);

        if (Q.degree() == 2) {
            SymbolicExpr c_expr = *(Q.coeffs[0].val);
            SymbolicExpr b_expr = *(Q.coeffs[1].val);
            SymbolicExpr a_expr = *(Q.coeffs[2].val);

            if (!a_expr.is_number() || !b_expr.is_number() || !c_expr.is_number()) {
                return nullptr;
            }

            double a = a_expr.to_numeric();
            double b = b_expr.to_numeric();
            double c = c_expr.to_numeric();

            int eq_a;
            lmmc_double_nearly_equal_tol(a, 0.0, 1e-9, 1e-9, &eq_a);
            if (eq_a) return nullptr;

            double delta = b * b - 4 * a * c;
            int eq_delta;
            lmmc_double_nearly_equal_tol(delta, 0.0, 1e-9, 1e-9, &eq_delta);

            if (!eq_delta && delta > 0) {
                double sqrt_delta = std::sqrt(delta);
                auto scalar = SymbolicExpr::number(1.0 / sqrt_delta);
                auto two_a = SymbolicExpr::number(2.0 * a);
                auto b_num = SymbolicExpr::number(b);
                auto two_a_x = SymbolicExpr::multiply(make_expr_ptr(*two_a), SymbolicExpr::variable(var));
                auto two_a_x_plus_b = SymbolicExpr::add(make_expr_ptr(*two_a_x), make_expr_ptr(*b_num));
                auto term1_arg = sym_sub(*two_a_x_plus_b, *SymbolicExpr::number(sqrt_delta));
                auto term2_arg = SymbolicExpr::add(make_expr_ptr(*two_a_x_plus_b), SymbolicExpr::number(sqrt_delta));
                auto term1 = SymbolicExpr::ln(make_expr_ptr(*term1_arg));
                auto term2 = SymbolicExpr::ln(make_expr_ptr(*term2_arg));
                return SymbolicExpr::multiply(scalar, sym_sub(*term1, *term2));
            } else if (!eq_delta && delta < 0) {
                double sqrt_neg_delta = std::sqrt(-delta);
                auto scalar = SymbolicExpr::number(2.0 / sqrt_neg_delta);
                auto two_a = SymbolicExpr::number(2.0 * a);
                auto b_num = SymbolicExpr::number(b);
                auto num = SymbolicExpr::add(
                    SymbolicExpr::multiply(make_expr_ptr(*two_a), SymbolicExpr::variable(var)),
                    make_expr_ptr(*b_num));
                auto inner = SymbolicExpr::divide(make_expr_ptr(*num), SymbolicExpr::number(sqrt_neg_delta));
                return SymbolicExpr::multiply(scalar, make_arctan(inner));
            }
        }
    } catch (...) {}

    return nullptr;
}

std::shared_ptr<SymbolicExpr> IBPStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth) {

    std::vector<std::shared_ptr<SymbolicNode>> ops;
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        ops = mul->operands;
    } else {
        ops.push_back(expr.root);
    }

    int best_u_idx = -1;
    int best_score = 100;

    auto get_score = [&](const std::shared_ptr<SymbolicNode>& node) -> int {
        SymbolicExpr e(node);
        if (auto fn = std::dynamic_pointer_cast<FunctionNode>(node)) {
            if (fn->type == FunctionNode::FuncType::Ln || fn->type == FunctionNode::FuncType::Log) return 1;
            if (fn->type == FunctionNode::FuncType::ArcSin || fn->type == FunctionNode::FuncType::ArcTan) return 2;
            if (fn->type == FunctionNode::FuncType::Sin || fn->type == FunctionNode::FuncType::Cos) return 4;
            if (fn->type == FunctionNode::FuncType::Exp) return 5;
        }
        if (!valid_dependency(e, var)) return 10;
        if (std::dynamic_pointer_cast<VariableNode>(node)) return 3;
        if (std::dynamic_pointer_cast<PowerNode>(node)) return 3;
        return 10;
    };

    if (ops.size() == 1) {
        int s = get_score(ops[0]);
        if (s <= 2) best_u_idx = 0;
    } else {
        for (size_t i = 0; i < ops.size(); ++i) {
            if (!valid_dependency(SymbolicExpr(ops[i]), var)) continue;
            int s = get_score(ops[i]);
            if (s < best_score) {
                best_score = s;
                best_u_idx = (int)i;
            }
        }
    }

    if (best_u_idx == -1) return nullptr;

    SymbolicExpr u(ops[best_u_idx]);

    std::vector<std::shared_ptr<SymbolicNode>> dv_ops;
    for (size_t i = 0; i < ops.size(); ++i) {
        if ((int)i != best_u_idx) dv_ops.push_back(ops[i]);
    }

    std::shared_ptr<SymbolicExpr> dv;
    if (dv_ops.empty()) dv = SymbolicExpr::number(1);
    else if (dv_ops.size() == 1) dv = make_expr_ptr(SymbolicExpr(dv_ops[0]));
    else dv = std::make_shared<SymbolicExpr>(std::make_shared<MultiplyNode>(dv_ops));

    auto v = ctx.integrate_recursive(*dv, var, depth + 1);
    if (has_integral_node_check(v->root)) return nullptr;

    auto du_ptr = u.differentiate(var);
    if (!du_ptr) return nullptr;
    auto du = make_expr_ptr(*du_ptr->simplify());

    auto uv = SymbolicExpr::multiply(make_expr_ptr(u), v);
    auto vdu = SymbolicExpr::multiply(v, du);
    auto int_vdu = ctx.integrate_recursive(*vdu, var, depth + 1);

    return sym_sub(*uv, *int_vdu);
}

// ---------------------------------------------------------------
// TrigCombinationStrategy
// ---------------------------------------------------------------
//
// Handles three families of trigonometric integrands whose argument is
// exactly the integration variable (linear arguments are handled earlier
// by LinearSubstitutionStrategy):
//
//   1. sin^m(var) * cos^n(var)  (m,n integers >= 0, m+n <= 8)
//   2. tan^n(var)                (n integer in [2,8])
//   3. sec^n(var)                (n even integer in [2,8])
//
// The sin/cos case dispatches on parity:
//   - At least one odd power -> peel off one factor, apply 1 = sin^2 + cos^2
//     to convert the remaining (even) part to a polynomial in u = cos/sin,
//     and integrate term-by-term.
//   - Both even -> apply the half-angle identities
//       sin^2(c x) = (1 - cos(2c x))/2,  cos^2(c x) = (1 + cos(2c x))/2,
//     expand the resulting binomial product into a polynomial in cos(2c x),
//     and recurse on each cos^k(2c x) term (the total degree halves each
//     time, so recursion terminates quickly).
// tan^n uses tan^n = tan^(n-2)*(sec^2 - 1) recursively, with bases
// tan^0 = 1 -> x and tan^1 -> -ln(cos(x)).
// Even-power sec^n uses the reduction
//     int sec^n = sec^(n-2)*tan/(n-1) + (n-2)/(n-1) * int sec^(n-2)
// with base sec^2 -> tan(x).

namespace {

// Build sin(c*var) / cos(c*var) where c is a non-zero integer scale.
std::shared_ptr<SymbolicExpr> make_sin_scaled(long long c, const std::string& var) {
    auto v = SymbolicExpr::variable(var);
    if (c == 1) return SymbolicExpr::sin(v);
    auto cx = SymbolicExpr::multiply(SymbolicExpr::number(static_cast<long long>(c)), v);
    return SymbolicExpr::sin(cx);
}

std::shared_ptr<SymbolicExpr> make_cos_scaled(long long c, const std::string& var) {
    auto v = SymbolicExpr::variable(var);
    if (c == 1) return SymbolicExpr::cos(v);
    auto cx = SymbolicExpr::multiply(SymbolicExpr::number(static_cast<long long>(c)), v);
    return SymbolicExpr::cos(cx);
}

// Convenience: a SymbolicExpr representing the rational num/den.
std::shared_ptr<SymbolicExpr> trig_rational(long long num, long long den) {
    return SymbolicExpr::number(Rational(BigInt(num), BigInt(den)));
}

// Binomial coefficient C(n,k) for small n; safe for n up to ~62.
long long binomial_ll(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    long long c = 1;
    for (int i = 0; i < k; ++i) {
        c = c * static_cast<long long>(n - i) / static_cast<long long>(i + 1);
    }
    return c;
}

// Match a single-argument FunctionNode whose argument is the integration
// variable itself. Returns 0 (sin), 1 (cos), 2 (tan), 3 (sec), or -1 on
// no match for this strategy.
int trig_match_of_var(const std::shared_ptr<SymbolicNode>& node, const std::string& var) {
    auto fn = std::dynamic_pointer_cast<FunctionNode>(node);
    if (!fn) return -1;
    if (fn->arguments.size() != 1) return -1;
    auto v = std::dynamic_pointer_cast<VariableNode>(fn->arguments[0]);
    if (!v || v->name != var) return -1;
    using FT = FunctionNode::FuncType;
    switch (fn->type) {
        case FT::Sin: return 0;
        case FT::Cos: return 1;
        case FT::Tan: return 2;
        case FT::Sec: return 3;
        default: return -1;
    }
}

// Try to extract a non-negative integer exponent.
bool trig_extract_nonneg_int(const std::shared_ptr<SymbolicNode>& exp_node, int& n_out) {
    SymbolicExpr e(exp_node);
    auto simp = e.simplify();
    if (!simp || !simp->is_int()) return false;
    int n = simp->get_int();
    if (n < 0) return false;
    n_out = n;
    return true;
}

// Match a single factor of the form trig(var) or trig(var)^k where trig is
// sin/cos/tan/sec and k is a non-negative integer. Returns true on success
// with the kind (0..3) and the integer power.
bool trig_extract_factor(const std::shared_ptr<SymbolicNode>& node,
                         const std::string& var,
                         int& kind_out, int& power_out) {
    int kind = trig_match_of_var(node, var);
    if (kind >= 0) {
        kind_out = kind;
        power_out = 1;
        return true;
    }
    auto pn = std::dynamic_pointer_cast<PowerNode>(node);
    if (!pn) return false;
    int base_kind = trig_match_of_var(pn->base, var);
    if (base_kind < 0) return false;
    int p = 0;
    if (!trig_extract_nonneg_int(pn->exponent, p)) return false;
    kind_out = base_kind;
    power_out = p;
    return true;
}

} // anonymous namespace

bool TrigCombinationStrategy::extract_sin_cos_powers(
    const SymbolicExpr& expr, const std::string& var, int& m_out, int& n_out) {

    int m = 0, n = 0;
    std::vector<std::shared_ptr<SymbolicNode>> factors;
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        factors = mul->operands;
    } else {
        factors.push_back(expr.root);
    }

    for (const auto& f : factors) {
        int kind = -1, p = 0;
        if (!trig_extract_factor(f, var, kind, p)) return false;
        // Only sin/cos contribute to (m,n); tan/sec disqualify this form.
        if (kind == 0) m += p;
        else if (kind == 1) n += p;
        else return false;
    }

    if (m == 0 && n == 0) return false; // not a non-trivial sin/cos product
    m_out = m;
    n_out = n;
    return true;
}

bool TrigCombinationStrategy::extract_tan_power(
    const SymbolicExpr& expr, const std::string& var, int& n_out) {
    int kind = -1, p = 0;
    if (!trig_extract_factor(expr.root, var, kind, p)) return false;
    if (kind != 2) return false;
    n_out = p;
    return true;
}

bool TrigCombinationStrategy::extract_sec_power(
    const SymbolicExpr& expr, const std::string& var, int& n_out) {
    int kind = -1, p = 0;
    if (!trig_extract_factor(expr.root, var, kind, p)) return false;
    if (kind != 3) return false;
    n_out = p;
    return true;
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth) {

    // sin^m(x) * cos^n(x)
    int m = 0, n = 0;
    if (extract_sin_cos_powers(expr, var, m, n)) {
        if (m < 0 || n < 0) return nullptr;
        if (m + n > 8) return nullptr;
        return integrate_sin_m_cos_n(m, n, 1, var, ctx, depth);
    }

    // tan^n(x)
    int tn = 0;
    if (extract_tan_power(expr, var, tn)) {
        if (tn < 2 || tn > 8) return nullptr;
        return integrate_tan_power(tn, var, ctx, depth);
    }

    // sec^n(x) -- only even powers per design.
    int sn = 0;
    if (extract_sec_power(expr, var, sn)) {
        if (sn < 2 || sn > 8) return nullptr;
        if (sn % 2 != 0) return nullptr;
        return integrate_sec_power(sn, var, ctx, depth);
    }

    return nullptr;
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_sin_m_cos_n(
    int m, int n, long long scale,
    const std::string& var, Integrator& ctx, int depth) {

    if (depth > ctx.max_depth()) return nullptr;
    if (m < 0 || n < 0) return nullptr;

    if (m == 0 && n == 0) {
        // Constant 1 with respect to var: integrand is essentially 1.
        return SymbolicExpr::variable(var);
    }
    if ((m % 2 == 1) || (n % 2 == 1)) {
        return integrate_odd_case(m, n, scale, var, ctx, depth);
    }
    return integrate_even_case(m, n, scale, var, ctx, depth);
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_odd_case(
    int m, int n, long long scale,
    const std::string& var, Integrator& ctx, int depth) {

    if (depth > ctx.max_depth()) return nullptr;
    if (scale == 0) return nullptr;

    // Pick which factor to peel:
    //   m odd  -> peel one sin, set u = cos(scale*x), du = -scale*sin(scale*x)dx
    //             so sin(scale*x)dx = -du/scale.
    //             remaining: sin^(m-1) cos^n = (1-u^2)^k * u^n, k=(m-1)/2.
    //   n odd  -> peel one cos, set u = sin(scale*x), du = +scale*cos(scale*x)dx
    //             so cos(scale*x)dx = du/scale.
    //             remaining: sin^m cos^(n-1) = u^m * (1-u^2)^k, k=(n-1)/2.
    bool peel_sin = (m % 2 == 1);
    int k = 0;
    int other_pow = 0;
    long long sign_factor = 1; // includes the sign from du sign.
    bool u_is_cos = false;

    if (peel_sin) {
        k = (m - 1) / 2;
        other_pow = n;
        sign_factor = -1;
        u_is_cos = true;
    } else {
        k = (n - 1) / 2;
        other_pow = m;
        sign_factor = 1;
        u_is_cos = false;
    }

    auto u_expr = u_is_cos ? make_cos_scaled(scale, var)
                            : make_sin_scaled(scale, var);

    // Result = sum_{i=0..k} C(k,i)*(-1)^i * u^(2i+other_pow+1) / (scale*(2i+other_pow+1)) * sign_factor
    std::vector<std::shared_ptr<SymbolicNode>> add_terms;
    for (int i = 0; i <= k; ++i) {
        long long bin = binomial_ll(k, i);
        long long alt_sign = ((i % 2) == 0) ? 1 : -1;
        long long num = sign_factor * alt_sign * bin;
        int u_pow = 2 * i + other_pow + 1; // always >= 1
        long long den = scale * static_cast<long long>(u_pow);
        if (den == 0) return nullptr;

        std::shared_ptr<SymbolicExpr> u_to_pow;
        if (u_pow == 1) {
            u_to_pow = u_expr;
        } else {
            u_to_pow = SymbolicExpr::power(u_expr, SymbolicExpr::number(u_pow));
        }
        // Normalize sign so the rational denominator is positive.
        if (den < 0) { num = -num; den = -den; }
        auto coeff = SymbolicExpr::number(Rational(BigInt(num), BigInt(den)));
        auto term = SymbolicExpr::multiply(coeff, u_to_pow);
        add_terms.push_back(term->root);
    }

    if (add_terms.empty()) return SymbolicExpr::number(0);
    if (add_terms.size() == 1) {
        return std::make_shared<SymbolicExpr>(add_terms[0]);
    }
    return std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(add_terms));
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_even_case(
    int m, int n, long long scale,
    const std::string& var, Integrator& ctx, int depth) {

    if (depth > ctx.max_depth()) return nullptr;
    if ((m % 2 != 0) || (n % 2 != 0)) return nullptr;
    int p = m / 2;
    int q = n / 2;
    int K = p + q;

    // sin^(2p)(c*x)*cos^(2q)(c*x)
    //   = (1/2^(p+q)) * sum_{i,j} C(p,i)(-1)^i C(q,j) cos^(i+j)(2c*x)
    // Compute coefficients a[k] = sum_{i+j=k} C(p,i)*C(q,j)*(-1)^i for k=0..K.
    std::vector<long long> a(K + 1, 0);
    for (int i = 0; i <= p; ++i) {
        long long bp = binomial_ll(p, i);
        long long si = ((i % 2) == 0) ? 1 : -1;
        for (int j = 0; j <= q; ++j) {
            long long bq = binomial_ll(q, j);
            a[i + j] += si * bp * bq;
        }
    }

    long long pow2 = 1;
    for (int t = 0; t < K; ++t) pow2 *= 2;
    if (pow2 == 0) pow2 = 1;

    std::vector<std::shared_ptr<SymbolicNode>> add_terms;
    for (int k = 0; k <= K; ++k) {
        if (a[k] == 0) continue;
        // Recurse: integrate cos^k(2c*x) at the new scale 2c.
        auto inner = integrate_sin_m_cos_n(0, k, scale * 2, var, ctx, depth + 1);
        if (!inner) return nullptr;
        long long num = a[k];
        long long den = pow2;
        if (den < 0) { num = -num; den = -den; }
        auto coeff = SymbolicExpr::number(Rational(BigInt(num), BigInt(den)));
        auto term = SymbolicExpr::multiply(coeff, inner);
        add_terms.push_back(term->root);
    }

    if (add_terms.empty()) return SymbolicExpr::number(0);
    if (add_terms.size() == 1) {
        return std::make_shared<SymbolicExpr>(add_terms[0]);
    }
    return std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(add_terms));
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_tan_power(
    int n, const std::string& var, Integrator& ctx, int depth) {

    if (depth > ctx.max_depth()) return nullptr;
    if (n < 0) return nullptr;
    auto v = SymbolicExpr::variable(var);

    if (n == 0) {
        // int 1 dx = x
        return v;
    }
    if (n == 1) {
        // int tan(x) dx = -ln(cos(x))
        auto cos_x = SymbolicExpr::cos(v);
        auto ln_cos = SymbolicExpr::ln(cos_x);
        return SymbolicExpr::multiply(SymbolicExpr::number(-1), ln_cos);
    }

    // n >= 2: tan^n = tan^(n-2)*(sec^2 - 1)
    //   int tan^n = tan^(n-1)/(n-1) - int tan^(n-2)
    auto tan_x = SymbolicExpr::tan(v);
    std::shared_ptr<SymbolicExpr> tan_pow_term;
    if (n - 1 == 1) {
        tan_pow_term = tan_x;
    } else {
        tan_pow_term = SymbolicExpr::power(tan_x, SymbolicExpr::number(n - 1));
    }
    auto first = SymbolicExpr::multiply(trig_rational(1, n - 1), tan_pow_term);

    auto rest = integrate_tan_power(n - 2, var, ctx, depth + 1);
    if (!rest) return nullptr;
    auto neg_rest = SymbolicExpr::multiply(SymbolicExpr::number(-1), rest);
    return SymbolicExpr::add(first, neg_rest);
}

std::shared_ptr<SymbolicExpr> TrigCombinationStrategy::integrate_sec_power(
    int n, const std::string& var, Integrator& ctx, int depth) {

    if (depth > ctx.max_depth()) return nullptr;
    if (n < 2 || (n % 2) != 0) return nullptr;

    using FT = FunctionNode::FuncType;
    auto v = SymbolicExpr::variable(var);
    auto tan_x = SymbolicExpr::tan(v);
    auto make_sec_x = [&]() -> std::shared_ptr<SymbolicExpr> {
        return std::make_shared<SymbolicExpr>(
            std::make_shared<FunctionNode>(FT::Sec,
                std::vector<std::shared_ptr<SymbolicNode>>{v->root}));
    };

    if (n == 2) {
        return tan_x;
    }

    auto sec_x = make_sec_x();
    std::shared_ptr<SymbolicExpr> sec_pow;
    if (n - 2 == 1) {
        sec_pow = sec_x;
    } else {
        sec_pow = SymbolicExpr::power(sec_x, SymbolicExpr::number(n - 2));
    }

    // First term: sec^(n-2)(x) * tan(x) / (n-1)
    auto inner = SymbolicExpr::multiply(sec_pow, tan_x);
    auto first = SymbolicExpr::multiply(trig_rational(1, n - 1), inner);

    // Second term: ((n-2)/(n-1)) * int sec^(n-2)
    auto rest = integrate_sec_power(n - 2, var, ctx, depth + 1);
    if (!rest) return nullptr;
    auto coeff = SymbolicExpr::number(Rational(BigInt(static_cast<long long>(n - 2)),
                                               BigInt(static_cast<long long>(n - 1))));
    auto second = SymbolicExpr::multiply(coeff, rest);

    return SymbolicExpr::add(first, second);
}

Integrator::Integrator() {

    strategies_.push_back(std::make_unique<TableLookupStrategy>());
    strategies_.push_back(std::make_unique<PowerRuleStrategy>());
    strategies_.push_back(std::make_unique<LinearSubstitutionStrategy>());
    strategies_.push_back(std::make_unique<SubstitutionStrategy>());
    strategies_.push_back(std::make_unique<TrigCombinationStrategy>());
    strategies_.push_back(std::make_unique<RationalDecompositionStrategy>());
    strategies_.push_back(std::make_unique<SpecialFunctionStrategy>());
    strategies_.push_back(std::make_unique<PartialFractionStrategy>());
    strategies_.push_back(std::make_unique<IBPStrategy>());
}

void Integrator::add_strategy(std::unique_ptr<IntegrationStrategy> strategy, int position) {
    if (!strategy) {
        throw std::invalid_argument("Integrator::add_strategy: strategy must not be null");
    }
    if (position < 0 || position >= (int)strategies_.size()) {
        strategies_.push_back(std::move(strategy));
    } else {
        strategies_.insert(strategies_.begin() + position, std::move(strategy));
    }
}

bool Integrator::depends_on(const SymbolicExpr& expr, const std::string& var) {
    return depends_on_var(expr.root, var);
}

std::shared_ptr<SymbolicExpr> Integrator::make_integral_node(
    const SymbolicExpr& expr, const std::string& var) {
    std::vector<std::shared_ptr<SymbolicNode>> args;
    args.push_back(expr.root);
    args.push_back(SymbolicExpr::variable(var)->root);
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args));
}

std::shared_ptr<SymbolicExpr> Integrator::check_cycle(
    const SymbolicExpr& expr, const std::string& var) {
    for (size_t i = 0; i < cycle_state_.history.size(); ++i) {
        auto ratio = SymbolicExpr::divide(make_expr_ptr(expr), make_expr_ptr(cycle_state_.history[i]))->simplify();
        if (!valid_dependency(*ratio, var)) {
            return SymbolicExpr::multiply(ratio, SymbolicExpr::variable("INT_CYCLE_" + std::to_string(i)));
        }
    }
    return nullptr;
}

void Integrator::resolve_cycle(std::shared_ptr<SymbolicExpr>& result, size_t cycle_idx) {
    std::string cycle_var = "INT_CYCLE_" + std::to_string(cycle_idx);
    if (valid_dependency(*result, cycle_var)) {
        auto B = result->differentiate(cycle_var)->simplify();
        auto A = result->substitute(cycle_var, SymbolicExpr::number(0))->simplify();
        auto one_minus_B = sym_sub(*SymbolicExpr::number(1), *B)->simplify();
        if (!one_minus_B->is_zero()) {
            result = SymbolicExpr::divide(A, one_minus_B);
        }
    }
}

std::shared_ptr<SymbolicExpr> Integrator::apply_linearity(
    const SymbolicExpr& expr, const std::string& var) {

    auto simp_expr = expr.simplify();

    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(simp_expr->root)) {
        std::vector<std::shared_ptr<SymbolicNode>> constants;
        std::vector<std::shared_ptr<SymbolicNode>> dependents;

        for (auto& op : mul->operands) {
            SymbolicExpr term(op);
            if (!valid_dependency(term, var)) {
                constants.push_back(op);
            } else {
                dependents.push_back(op);
            }
        }

        if (!constants.empty() && dependents.size() < mul->operands.size()) {
            SymbolicExpr const_part = (constants.size() == 1) ?
                SymbolicExpr(constants[0]) : SymbolicExpr(std::make_shared<MultiplyNode>(constants));
            SymbolicExpr dep_part = (dependents.empty()) ?
                *SymbolicExpr::number(1) :
                ((dependents.size() == 1) ? SymbolicExpr(dependents[0]) : SymbolicExpr(std::make_shared<MultiplyNode>(dependents)));

            auto int_part = integrate(dep_part, var);
            return SymbolicExpr::multiply(std::make_shared<SymbolicExpr>(const_part), std::make_shared<SymbolicExpr>(int_part));
        }
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(simp_expr->root)) {
        std::vector<std::shared_ptr<SymbolicNode>> results;
        for (auto& op : add->operands) {
            SymbolicExpr term(op);
            auto int_term = integrate(term, var);
            results.push_back(int_term.root);
        }
        return std::make_shared<SymbolicExpr>(std::make_shared<AddNode>(results));
    }

    return nullptr;
}

std::shared_ptr<SymbolicExpr> Integrator::dispatch_strategies(
    const SymbolicExpr& expr, const std::string& var, int depth) {

    for (auto& strategy : strategies_) {
        auto result = strategy->try_integrate(expr, var, *this, depth);
        if (result) {
            err_stream << "[Integration] Strategy '" << strategy->name() << "' succeeded\n";
            return result;
        }
    }
    return nullptr;
}

std::shared_ptr<SymbolicExpr> Integrator::integrate_recursive(
    const SymbolicExpr& expr, const std::string& var, int depth) {

    if (depth > max_depth_) {
        return make_integral_node(expr, var);
    }

    auto cycle_result = check_cycle(expr, var);
    if (cycle_result) return cycle_result;

    cycle_state_.history.push_back(expr);
    size_t my_idx = cycle_state_.history.size() - 1;

    if (expr.is_number() || !valid_dependency(expr, var)) {
        cycle_state_.history.pop_back();
        return SymbolicExpr::multiply(make_expr_ptr(expr), SymbolicExpr::variable(var));
    }

    auto result = dispatch_strategies(expr, var, depth);

    cycle_state_.history.pop_back();

    if (!result) {
        return make_integral_node(expr, var);
    }

    resolve_cycle(result, my_idx);

    return result;
}

SymbolicExpr Integrator::integrate(const SymbolicExpr& expr, const std::string& var_name) {

    cycle_state_.history.clear();

    auto linear_result = apply_linearity(expr, var_name);
    if (linear_result) return *linear_result;

    return *integrate_recursive(expr, var_name, 0);
}

SymbolicExpr Integrator::integrate_def(const SymbolicExpr& expr, const std::string& var_name,
                                        const SymbolicExpr& lower, const SymbolicExpr& upper) {

    SymbolicExpr simp_expr_val = *expr.simplify();
    bool is_inv_x = false;

    if (auto pow = std::dynamic_pointer_cast<PowerNode>(simp_expr_val.root)) {
        if (auto v = std::dynamic_pointer_cast<VariableNode>(pow->base)) {
            if (v->name == var_name) {
                if (auto en = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                    int eq_minus_one = 0;
                    if (std::holds_alternative<lmmc_real_t>(en->value)) {
                        lmmc_double_nearly_equal_tol(std::get<lmmc_real_t>(en->value), -1.0, 1e-9, 1e-9, &eq_minus_one);
                    }
                    if ((std::holds_alternative<lmmc_real_t>(en->value) && eq_minus_one != 0) ||
                        (std::holds_alternative<BigInt>(en->value) && std::get<BigInt>(en->value).to_int() == -1) ||
                        (std::holds_alternative<Rational>(en->value) && std::get<Rational>(en->value).to_double() == -1.0)) {
                        is_inv_x = true;
                    }
                }
            }
        }
    }

    bool numeric_bounds = (lower.root && std::dynamic_pointer_cast<NumberNode>(lower.root)) &&
                          (upper.root && std::dynamic_pointer_cast<NumberNode>(upper.root));

    if (is_inv_x && numeric_bounds) {
        double l_val = lower.to_numeric();
        double u_val = upper.to_numeric();
        int eq_l, eq_u;
        lmmc_double_nearly_equal_tol(l_val, 0.0, 1e-9, 1e-9, &eq_l);
        lmmc_double_nearly_equal_tol(u_val, 0.0, 1e-9, 1e-9, &eq_u);
        if (!eq_l && l_val < 0 && !eq_u && u_val > 0) {
            auto t = std::make_shared<SymbolicExpr>(*SymbolicExpr::variable("t"));
            auto zero = std::make_shared<SymbolicExpr>(*SymbolicExpr::number(0));
            auto int_left = integrate_def(expr, var_name, lower, *t);
            auto lim_left = int_left.limit("t", zero, "-");
            auto int_right = integrate_def(expr, var_name, *t, upper);
            auto lim_right = int_right.limit("t", zero, "+");
            if (lim_left && lim_right) {
                return *SymbolicExpr::add(lim_left, lim_right);
            }
        }
    }

    SymbolicExpr indefinite = integrate(expr, var_name);

    if (auto func = std::dynamic_pointer_cast<FunctionNode>(indefinite.root)) {
        if (func->type == FunctionNode::FuncType::Calculus_Integral) {
            std::vector<std::shared_ptr<SymbolicNode>> args;
            args.push_back(expr.root);
            args.push_back(SymbolicExpr::variable(var_name)->root);
            args.push_back(lower.root);
            args.push_back(upper.root);
            return SymbolicExpr(std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args));
        }
    }

    auto F_b = indefinite.substitute(var_name, make_expr_ptr(upper));
    auto F_a = indefinite.substitute(var_name, make_expr_ptr(lower));
    auto result = sym_sub(*F_b, *F_a);
    return *result->simplify();
}

// ---------------------------------------------------------------
// RationalDecompositionStrategy
// ---------------------------------------------------------------
//
// Handles rational integrands P(x)/Q(x) where deg(Q) >= 3 (lower-degree
// quadratic cases are caught earlier by PartialFractionStrategy). The
// algorithm follows the textbook recipe:
//
//   1. extract_rational : factor the integrand into (numerator polynomial,
//      denominator polynomial) over Q. Returns false if the input is not a
//      rational function (e.g. contains sin, exp, ln, irrational powers).
//   2. poly_divide      : long-divide P by Q when deg(P) >= deg(Q). The
//      polynomial quotient is integrated term-by-term via the power rule.
//   3. factor_denominator : square-free factor Q, peel off rational linear
//      roots, leaving each leftover factor either linear or an irreducible
//      quadratic. Refuses to handle higher-degree irreducible factors.
//   4. solve_coefficients : set up the partial-fraction ansatz with one
//      unknown per linear power and two unknowns per irreducible quadratic
//      power, expand into a linear system whose unknowns are the partial
//      fraction coefficients, and solve via gaussian_eliminate.
//   5. integrate_term  : integrate each component analytically.
//
// On any failure that is not "expression is not rational" we fall back to
// returning an unevaluated integral node so that the rest of the pipeline
// is preserved.

namespace {

// Canonical builder utilities used throughout the rational strategy.
inline std::shared_ptr<SymbolicExpr> rd_num_rat(const Rational& r) {
    return SymbolicExpr::number(r);
}

inline std::shared_ptr<SymbolicExpr> rd_num_int(long long v) {
    return SymbolicExpr::number(BigInt(v));
}

// Build (var - r) as a SymbolicExpr.
inline std::shared_ptr<SymbolicExpr> rd_var_minus(const std::string& var, const Rational& r) {
    auto v = SymbolicExpr::variable(var);
    if (r == Rational(0)) return v;
    auto neg_r = SymbolicExpr::number(Rational(0) - r);
    return SymbolicExpr::add(v, neg_r);
}

// Test whether a rational-coefficient poly is the zero polynomial.
inline bool rd_is_zero_poly(const Polynomial<Rational>& p) {
    if (p.coeffs.empty()) return true;
    for (const auto& c : p.coeffs) if (!(c == Rational(0))) return false;
    return true;
}

// Convert a Polynomial<Rational> into a SymbolicExpr in the given variable.
inline std::shared_ptr<SymbolicExpr> rd_poly_to_sym(
    const Polynomial<Rational>& p, const std::string& var) {
    if (rd_is_zero_poly(p)) return SymbolicExpr::number(0);
    std::vector<std::shared_ptr<SymbolicExpr>> terms;
    auto v = SymbolicExpr::variable(var);
    for (size_t i = 0; i < p.coeffs.size(); ++i) {
        if (p.coeffs[i] == Rational(0)) continue;
        std::shared_ptr<SymbolicExpr> t;
        if (i == 0) {
            t = rd_num_rat(p.coeffs[i]);
        } else if (i == 1) {
            if (p.coeffs[i] == Rational(1)) {
                t = v;
            } else {
                t = SymbolicExpr::multiply(rd_num_rat(p.coeffs[i]), v);
            }
        } else {
            auto pw = SymbolicExpr::power(v, rd_num_int(static_cast<long long>(i)));
            if (p.coeffs[i] == Rational(1)) {
                t = pw;
            } else {
                t = SymbolicExpr::multiply(rd_num_rat(p.coeffs[i]), pw);
            }
        }
        terms.push_back(t);
    }
    if (terms.empty()) return SymbolicExpr::number(0);
    if (terms.size() == 1) return terms[0];
    auto res = terms[0];
    for (size_t i = 1; i < terms.size(); ++i) res = SymbolicExpr::add(res, terms[i]);
    return res;
}

// Recursively decompose a SymbolicNode into (numerator-polys, denominator-polys)
// over the rationals. Each polynomial is assumed to be a polynomial in `var`
// with rational coefficients. Returns false on any non-rational sub-expression
// (functions, irrational/symbolic exponents, variables other than the integration
// variable). Multiplication aggregates by appending; division comes from
// PowerNode with negative integer exponent. Numbers and the integration
// variable are converted directly via symbolic_to_poly.
bool rd_collect_rational(const std::shared_ptr<SymbolicNode>& node,
                         const std::string& var,
                         std::vector<Polynomial<Rational>>& num,
                         std::vector<Polynomial<Rational>>& den) {
    if (!node) return false;

    // A NumberNode is a constant; convert directly.
    if (auto n = std::dynamic_pointer_cast<NumberNode>(node)) {
        Polynomial<Rational> p =
            symbolic_to_poly<Rational>(std::make_shared<SymbolicExpr>(n), var);
        num.push_back(p);
        return true;
    }

    // The integration variable, or any other variable that's actually a constant.
    if (auto v = std::dynamic_pointer_cast<VariableNode>(node)) {
        if (v->name == var) {
            num.push_back(Polynomial<Rational>({Rational(0), Rational(1)}, var));
            return true;
        }
        // Unknown symbolic constant: cannot lift safely into rational poly.
        return false;
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        // Addition: combine all summands over a common denominator.
        // We collect each operand as P_i/Q_i and assemble
        //   sum P_i * (prod_{j != i} Q_j)  /  prod_j Q_j
        std::vector<std::pair<std::vector<Polynomial<Rational>>,
                              std::vector<Polynomial<Rational>>>> parts;
        parts.reserve(add->operands.size());
        for (const auto& op : add->operands) {
            std::vector<Polynomial<Rational>> sub_num, sub_den;
            if (!rd_collect_rational(op, var, sub_num, sub_den)) return false;
            parts.push_back({std::move(sub_num), std::move(sub_den)});
        }

        // Combine: numerator = sum_i (prod_k num_k_i) * prod_{j != i} (prod_k den_k_j)
        //          denominator = prod_i (prod_k den_k_i)
        Polynomial<Rational> total_num(var);    // start as 0
        Polynomial<Rational> total_den({Rational(1)}, var);

        // Pre-compute Q_i (denominator product per term) and total denominator.
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
            // multiply by prod_{j != i} Q_j
            for (size_t j = 0; j < parts.size(); ++j) {
                if (j == i) continue;
                p = p * Q_per_term[j];
            }
            total_num = total_num + p;
        }

        num.push_back(total_num);
        den.push_back(total_den);
        return true;
    }

    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        for (const auto& op : mul->operands) {
            if (!rd_collect_rational(op, var, num, den)) return false;
        }
        return true;
    }

    if (auto pw = std::dynamic_pointer_cast<PowerNode>(node)) {
        // Only integer exponents (negative or non-negative) are allowed.
        auto en = std::dynamic_pointer_cast<NumberNode>(pw->exponent);
        if (!en) return false;

        long long exp_v = 0;
        bool ok = false;
        if (std::holds_alternative<BigInt>(en->value)) {
            const auto& bi = std::get<BigInt>(en->value);
            // Bound to keep matrix sizes sane.
            int v = bi.to_int();
            if (v >= -64 && v <= 64) { ok = true; exp_v = v; }
        } else if (std::holds_alternative<Rational>(en->value)) {
            const auto& r = std::get<Rational>(en->value);
            if (r.is_integer()) {
                int v = r.to_BigInt().to_int();
                if (v >= -64 && v <= 64) { ok = true; exp_v = v; }
            }
        } else if (std::holds_alternative<lmmc_real_t>(en->value)) {
            lmmc_real_t d = std::get<lmmc_real_t>(en->value);
            if (std::isfinite(d) && d == std::floor(d) && d >= -64.0 && d <= 64.0) {
                ok = true; exp_v = static_cast<long long>(d);
            }
        }
        if (!ok) return false;

        // Build the base polynomial.
        std::vector<Polynomial<Rational>> bn, bd;
        if (!rd_collect_rational(pw->base, var, bn, bd)) return false;

        // base_num = product of bn ; base_den = product of bd
        Polynomial<Rational> base_num({Rational(1)}, var);
        for (const auto& p : bn) base_num = base_num * p;
        Polynomial<Rational> base_den({Rational(1)}, var);
        for (const auto& p : bd) base_den = base_den * p;

        if (exp_v == 0) {
            num.push_back(Polynomial<Rational>({Rational(1)}, var));
            return true;
        }
        long long e_abs = std::llabs(exp_v);
        Polynomial<Rational> n_pow({Rational(1)}, var);
        Polynomial<Rational> d_pow({Rational(1)}, var);
        for (long long i = 0; i < e_abs; ++i) {
            n_pow = n_pow * base_num;
            d_pow = d_pow * base_den;
        }
        if (exp_v > 0) {
            num.push_back(n_pow);
            den.push_back(d_pow);
        } else {
            // (P/Q)^(-k) = (Q/P)^k
            num.push_back(d_pow);
            den.push_back(n_pow);
        }
        return true;
    }

    // FunctionNode etc. -> not a rational function in `var`.
    return false;
}

// Returns true and sets out_n if `c` is a non-negative integer NumberNode-like
// constant (kept simple; we only need this for safety checks).
inline bool rd_is_const_in(const std::shared_ptr<SymbolicExpr>& expr,
                           const std::string& var) {
    return !depends_on_var(expr->root, var);
}

} // anonymous namespace

bool RationalDecompositionStrategy::extract_rational(
    const SymbolicExpr& expr, const std::string& var,
    Polynomial<Rational>& P_out, Polynomial<Rational>& Q_out) {

    std::vector<Polynomial<Rational>> nums, dens;
    if (!rd_collect_rational(expr.root, var, nums, dens)) return false;

    Polynomial<Rational> P({Rational(1)}, var);
    for (const auto& p : nums) P = P * p;
    Polynomial<Rational> Q({Rational(1)}, var);
    for (const auto& p : dens) Q = Q * p;

    if (rd_is_zero_poly(Q)) return false; // 0 in denominator, refuse

    // Reduce by GCD to keep things small (and keep variable name consistent).
    if (!rd_is_zero_poly(P)) {
        Polynomial<Rational> g = Polynomial<Rational>::gcd(P, Q);
        if (!rd_is_zero_poly(g) && g.degree() >= 1) {
            auto pq = P.div_mod(g);
            auto qq = Q.div_mod(g);
            if (rd_is_zero_poly(pq.second) && rd_is_zero_poly(qq.second)) {
                P = pq.first;
                Q = qq.first;
            }
        }
    }

    // Make Q monic (rescale numerator accordingly).
    if (Q.degree() >= 0) {
        Rational lc = Q.lead_coeff();
        if (!(lc == Rational(1)) && !(lc == Rational(0))) {
            Polynomial<Rational> P_scaled(var);
            P_scaled.coeffs.reserve(P.coeffs.size());
            for (const auto& c : P.coeffs) P_scaled.coeffs.push_back(c / lc);
            Polynomial<Rational> Q_scaled(var);
            Q_scaled.coeffs.reserve(Q.coeffs.size());
            for (const auto& c : Q.coeffs) Q_scaled.coeffs.push_back(c / lc);
            P_scaled.trim();
            Q_scaled.trim();
            P = P_scaled;
            Q = Q_scaled;
        }
    }

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

bool RationalDecompositionStrategy::factor_denominator(
    const Polynomial<Rational>& Q,
    std::vector<std::pair<Polynomial<Rational>, int>>& factors_out) {
    factors_out.clear();
    if (Q.degree() <= 0) return false;

    // Step 1: square-free factorization to peel off multiplicities.
    auto sqfree = square_free_factorization(Q);
    if (sqfree.empty()) return false;

    // Step 2: for each square-free factor, peel rational roots, leaving an
    //         irreducible-quadratic cofactor at most. Anything that cannot be
    //         reduced to factors of degree <= 2 forces this step to fail.
    for (auto& [piece, mult] : sqfree) {
        Polynomial<Rational> current = piece.make_monic();

        if (current.degree() == 0) {
            // pure constant factor; nothing to integrate, ignore
            continue;
        }
        if (current.degree() == 1 || current.degree() == 2) {
            // Irreducible quadratics can stay as-is; quadratic discriminant
            // determines splitting later in solve_coefficients/integrate_term.
            // For deg 2 we still try rational roots to split into two linears
            // if possible.
            if (current.degree() == 2) {
                auto roots = find_rational_roots(current);
                if (!roots.empty()) {
                    for (const auto& r : roots) {
                        Polynomial<Rational> linear({Rational(0) - r, Rational(1)}, current.variable_name);
                        auto qr = current.div_mod(linear);
                        if (!rd_is_zero_poly(qr.second)) {
                            return false;
                        }
                        factors_out.push_back({linear, mult});
                        current = qr.first;
                        if (current.degree() == 0) break;
                    }
                }
                if (current.degree() == 2) {
                    // Confirmed irreducible quadratic over Q.
                    factors_out.push_back({current.make_monic(), mult});
                    current = Polynomial<Rational>({Rational(1)}, current.variable_name);
                } else if (current.degree() == 1) {
                    factors_out.push_back({current.make_monic(), mult});
                    current = Polynomial<Rational>({Rational(1)}, current.variable_name);
                }
                continue;
            }
            // degree 1
            factors_out.push_back({current, mult});
            continue;
        }

        // degree >= 3: try to split off rational linear factors.
        auto roots = find_rational_roots(current);
        for (const auto& r : roots) {
            Polynomial<Rational> linear({Rational(0) - r, Rational(1)}, current.variable_name);
            auto qr = current.div_mod(linear);
            if (rd_is_zero_poly(qr.second)) {
                factors_out.push_back({linear, mult});
                current = qr.first;
                if (current.degree() <= 0) break;
            }
        }

        if (current.degree() <= 0) continue;
        if (current.degree() == 1) {
            factors_out.push_back({current.make_monic(), mult});
            continue;
        }
        if (current.degree() == 2) {
            factors_out.push_back({current.make_monic(), mult});
            continue;
        }
        // Higher-degree irreducible factor over Q: cannot handle in this
        // strategy. Caller will return an unevaluated integral.
        return false;
    }
    return true;
}

bool RationalDecompositionStrategy::solve_coefficients(
    const Polynomial<Rational>& P,
    const Polynomial<Rational>& Q,
    const std::vector<std::pair<Polynomial<Rational>, int>>& factors,
    std::vector<Polynomial<Rational>>& numerators_out) {

    numerators_out.clear();

    // Map each (factor index, power) to a list of unknowns.
    // Linear factor (x - r)^k contributes k unknowns A_1..A_k (constant numerators).
    // Quadratic factor q(x)^k contributes 2k unknowns: (B_l, C_l) per power l.
    //
    // We also remember, for each (factor index, power l), the partial-fraction
    // term's *full* numerator polynomial (constant or linear) as a function of
    // the unknowns. Then we multiply each term by the missing portion of Q to
    // get a polynomial whose coefficients (in x) are linear combinations of the
    // unknowns; equating to P gives the linear system.

    // Total number of unknowns.
    size_t num_unknowns = 0;
    std::vector<std::pair<int, int>> ifac_ipow; // (factor_index, power index l in [1..mult])
    std::vector<int> kind_per_unknown_block;    // 1 -> linear, 2 -> quadratic
    std::vector<size_t> unknown_offset_per_term;
    unknown_offset_per_term.reserve(factors.size());

    for (size_t i = 0; i < factors.size(); ++i) {
        unknown_offset_per_term.push_back(num_unknowns);
        const auto& [fpoly, mult] = factors[i];
        if (fpoly.degree() == 1) {
            num_unknowns += static_cast<size_t>(mult);
        } else if (fpoly.degree() == 2) {
            num_unknowns += static_cast<size_t>(2 * mult);
        } else {
            return false;
        }
    }
    if (num_unknowns == 0) return false;

    const std::string& var = Q.variable_name;
    int N = Q.degree();
    if (N < 0) return false;

    // Build the matrix: each column is the contribution (coefficient vector
    // of x^0..x^{N-1}) of one unknown to the LHS polynomial; each row is one
    // equation indexed by the power of x.
    // Number of equations = N (size of polynomial of degree < N for the RHS P
    // when properly partial-fractioned: deg(P) < N).
    size_t rows = static_cast<size_t>(N);
    size_t cols = num_unknowns + 1; // augmented column for P

    // Initialize augmented matrix as Rational.
    std::vector<std::vector<Rational>> M(rows, std::vector<Rational>(cols, Rational(0)));

    // RHS column: coefficients of P (degree < N expected; if deg(P) >= N
    // then we have a problem — caller must run poly_divide first).
    if (P.degree() >= N) return false;
    for (size_t k = 0; k < rows; ++k) {
        Rational c = (k < P.coeffs.size()) ? P.coeffs[k] : Rational(0);
        M[k][num_unknowns] = c;
    }

    // For each factor, for each multiplicity power l, build the contribution
    // poly = numerator_term * (Q / fpoly^l). Each unknown in this term
    // contributes one column.
    for (size_t i = 0; i < factors.size(); ++i) {
        const auto& [fpoly, mult] = factors[i];
        // Q / fpoly^l for l = 1..mult; build incrementally: start with Q/fpoly^mult
        // by dividing Q by fpoly mult times.
        Polynomial<Rational> remaining = Q;
        for (int t = 0; t < mult; ++t) {
            auto qr = remaining.div_mod(fpoly);
            if (!rd_is_zero_poly(qr.second)) {
                // Q does not divide cleanly by this power of factor -> bug.
                return false;
            }
            remaining = qr.first;
        }
        // remaining now equals Q / fpoly^mult.

        // Now build for each power l = 1..mult:
        //   term coefficient pattern is multiplied by Q / fpoly^l
        //   = remaining * fpoly^(mult - l)
        Polynomial<Rational> fp_pow({Rational(1)}, var);
        // fp_pow starts as fpoly^0 = 1; we will iterate l from mult down to 1
        // so we increment fp_pow by * fpoly each time after using it.
        // Pre-build the array of coefficient polynomials for each l in 1..mult.
        std::vector<Polynomial<Rational>> q_over_fpow_l(mult + 1, Polynomial<Rational>(var));
        // index by l from 1..mult: q_over_fpow_l[l] = remaining * fpoly^(mult - l)
        Polynomial<Rational> cur({Rational(1)}, var);
        for (int l = mult; l >= 1; --l) {
            // when l = mult: cur = fpoly^0 = 1, factor = remaining * 1 = remaining
            // when decreasing l: factor = remaining * fpoly^(mult - l)
            Polynomial<Rational> factor_poly = remaining * cur;
            q_over_fpow_l[l] = factor_poly;
            cur = cur * fpoly;
        }

        size_t col_off = unknown_offset_per_term[i];

        if (fpoly.degree() == 1) {
            // Linear factor. Each l contributes one unknown A_l (constant
            // numerator). Coefficient column = q_over_fpow_l[l].
            for (int l = 1; l <= mult; ++l) {
                const auto& qpl = q_over_fpow_l[l];
                size_t this_col = col_off + static_cast<size_t>(l - 1);
                for (size_t k = 0; k < rows; ++k) {
                    Rational c = (k < qpl.coeffs.size()) ? qpl.coeffs[k] : Rational(0);
                    M[k][this_col] = c;
                }
            }
        } else if (fpoly.degree() == 2) {
            // Quadratic factor. Each l contributes two unknowns (B_l, C_l).
            // Numerator at power l = B_l * x + C_l ; the coefficient column for
            // B_l is q_over_fpow_l[l] shifted by 1 (multiplied by x), for C_l
            // is q_over_fpow_l[l] (multiplied by 1).
            for (int l = 1; l <= mult; ++l) {
                const auto& qpl = q_over_fpow_l[l];
                size_t base_col = col_off + static_cast<size_t>(2 * (l - 1));
                size_t b_col = base_col;     // for B_l
                size_t c_col = base_col + 1; // for C_l
                // C_l contributes qpl directly.
                for (size_t k = 0; k < rows; ++k) {
                    Rational c = (k < qpl.coeffs.size()) ? qpl.coeffs[k] : Rational(0);
                    M[k][c_col] = c;
                }
                // B_l contributes qpl * x = qpl with index shifted by 1.
                for (size_t k = 0; k < rows; ++k) {
                    if (k == 0) {
                        M[k][b_col] = Rational(0);
                    } else {
                        size_t src = k - 1;
                        Rational c = (src < qpl.coeffs.size()) ? qpl.coeffs[src] : Rational(0);
                        M[k][b_col] = c;
                    }
                }
            }
        }
    }

    // Solve via Gaussian elimination on Rational matrix.
    size_t m = rows;
    size_t n = cols;
    int sign = 1;
    (void)sign;
    size_t pivot_row = 0;
    std::vector<size_t> pivot_col_for_row(m, static_cast<size_t>(-1));
    for (size_t col = 0; col + 1 < n && pivot_row < m; ++col) {
        size_t max_row = pivot_row;
        bool found = false;
        while (max_row < m) {
            if (!(M[max_row][col] == Rational(0))) { found = true; break; }
            ++max_row;
        }
        if (!found) continue;
        if (max_row != pivot_row) std::swap(M[pivot_row], M[max_row]);
        pivot_col_for_row[pivot_row] = col;
        Rational pivot = M[pivot_row][col];
        // Normalize pivot row.
        for (size_t k = col; k < n; ++k) {
            M[pivot_row][k] = M[pivot_row][k] / pivot;
        }
        // Eliminate in other rows.
        for (size_t i = 0; i < m; ++i) {
            if (i == pivot_row) continue;
            Rational f = M[i][col];
            if (f == Rational(0)) continue;
            for (size_t k = col; k < n; ++k) {
                M[i][k] = M[i][k] - f * M[pivot_row][k];
            }
        }
        ++pivot_row;
    }

    // Check for inconsistency: any row with all-zero LHS but nonzero RHS.
    for (size_t r = 0; r < m; ++r) {
        bool all_zero_lhs = true;
        for (size_t c = 0; c + 1 < n; ++c) {
            if (!(M[r][c] == Rational(0))) { all_zero_lhs = false; break; }
        }
        if (all_zero_lhs && !(M[r][n - 1] == Rational(0))) return false;
    }

    // Extract solution: requires that each unknown has a pivot column.
    std::vector<Rational> sol(num_unknowns, Rational(0));
    std::vector<bool> pivoted(num_unknowns, false);
    for (size_t r = 0; r < m; ++r) {
        size_t pc = pivot_col_for_row[r];
        if (pc == static_cast<size_t>(-1)) continue;
        if (pc >= num_unknowns) continue;
        sol[pc] = M[r][n - 1];
        pivoted[pc] = true;
    }
    for (size_t i = 0; i < num_unknowns; ++i) {
        if (!pivoted[i]) return false; // underdetermined -> can't form unique partial fractions
    }

    // Re-package solution into per-(factor, power) numerator polynomials.
    numerators_out.reserve(num_unknowns);
    for (size_t i = 0; i < factors.size(); ++i) {
        const auto& [fpoly, mult] = factors[i];
        size_t col_off = unknown_offset_per_term[i];
        if (fpoly.degree() == 1) {
            for (int l = 1; l <= mult; ++l) {
                Rational A = sol[col_off + static_cast<size_t>(l - 1)];
                numerators_out.push_back(Polynomial<Rational>({A}, Q.variable_name));
            }
        } else { // degree 2
            for (int l = 1; l <= mult; ++l) {
                size_t base_col = col_off + static_cast<size_t>(2 * (l - 1));
                Rational B = sol[base_col];
                Rational C = sol[base_col + 1];
                numerators_out.push_back(
                    Polynomial<Rational>({C, B}, Q.variable_name));
            }
        }
    }
    return true;
}

std::shared_ptr<SymbolicExpr> RationalDecompositionStrategy::integrate_term(
    const Polynomial<Rational>& numerator,
    const Polynomial<Rational>& factor,
    int power, const std::string& var) {

    if (rd_is_zero_poly(numerator)) return SymbolicExpr::number(0);

    if (factor.degree() == 1) {
        // factor = x - r (monic linear). r = -coeffs[0].
        Rational r = Rational(0) - factor.coeffs[0];
        Rational A = (numerator.coeffs.size() > 0) ? numerator.coeffs[0] : Rational(0);
        if (A == Rational(0)) return SymbolicExpr::number(0);

        if (power == 1) {
            // A * ln|x - r|. We use ln(x - r) here (matching existing strategies
            // which omit the absolute value).
            auto inner = rd_var_minus(var, r);
            return SymbolicExpr::multiply(rd_num_rat(A), SymbolicExpr::ln(inner));
        }
        // power >= 2: -A / ((power-1) * (x - r)^(power-1))
        Rational coeff = Rational(0) - A / Rational(BigInt(power - 1));
        auto base = rd_var_minus(var, r);
        auto exp = rd_num_int(power - 1);
        auto pw = SymbolicExpr::power(base, exp);
        auto inv_pw = SymbolicExpr::power(pw, rd_num_int(-1));
        return SymbolicExpr::multiply(rd_num_rat(coeff), inv_pw);
    }

    if (factor.degree() == 2) {
        // factor = x^2 + p*x + q (monic quadratic). discriminant = p^2 - 4q.
        Rational q = (factor.coeffs.size() > 0) ? factor.coeffs[0] : Rational(0);
        Rational p = (factor.coeffs.size() > 1) ? factor.coeffs[1] : Rational(0);

        // Verify it's irreducible: discriminant must be < 0. (If it's >= 0
        // we should have factored it; but be defensive.)
        Rational disc = p * p - Rational(4) * q;

        // numerator = B*x + C
        Rational B = (numerator.coeffs.size() > 1) ? numerator.coeffs[1] : Rational(0);
        Rational C = (numerator.coeffs.size() > 0) ? numerator.coeffs[0] : Rational(0);

        if (power >= 2) {
            // Higher powers of irreducible quadratics: leave as unevaluated
            // integral; this triggers the safe fallback in try_integrate.
            // Build the symbolic integrand for the unevaluated node.
            auto num_sym = rd_poly_to_sym(numerator, var);
            auto den_sym_base = rd_poly_to_sym(factor, var);
            auto den_pw = SymbolicExpr::power(den_sym_base, rd_num_int(power));
            auto inv_den = SymbolicExpr::power(den_pw, rd_num_int(-1));
            auto integrand = SymbolicExpr::multiply(num_sym, inv_den);
            std::vector<std::shared_ptr<SymbolicNode>> args;
            args.push_back(integrand->root);
            args.push_back(SymbolicExpr::variable(var)->root);
            return std::make_shared<SymbolicExpr>(
                std::make_shared<FunctionNode>(FunctionNode::FuncType::Calculus_Integral, args));
        }

        // power == 1: ∫ (B x + C) / (x^2 + p x + q) dx
        //   Split numerator: (B x + C) = (B/2) * (2x + p) + (C - B p / 2)
        //   ∫ (B/2)(2x+p)/(x^2+px+q) dx = (B/2) * ln(x^2+px+q)
        //   ∫ (C - B p / 2) / (x^2 + p x + q) dx
        //     = (C - B p / 2) * (2 / sqrt(4q - p^2)) * arctan( (2x + p) / sqrt(4q - p^2) )
        std::shared_ptr<SymbolicExpr> result = SymbolicExpr::number(0);

        Rational B_half = B / Rational(2);
        if (!(B_half == Rational(0))) {
            auto den_sym = rd_poly_to_sym(factor, var);
            auto ln_part = SymbolicExpr::ln(den_sym);
            auto term_log = SymbolicExpr::multiply(rd_num_rat(B_half), ln_part);
            result = SymbolicExpr::add(result, term_log);
        }

        Rational atan_coeff = C - (B * p) / Rational(2);
        if (!(atan_coeff == Rational(0))) {
            // s2 = -disc = 4q - p^2 (positive for irreducible factor).
            Rational s2 = Rational(0) - disc;
            // Build sqrt(s2) symbolically (works even when s2 isn't a perfect square).
            auto s2_sym = rd_num_rat(s2);
            auto sqrt_s2 = SymbolicExpr::sqrt(s2_sym);
            // 2x + p
            auto two_x = SymbolicExpr::multiply(rd_num_int(2), SymbolicExpr::variable(var));
            std::shared_ptr<SymbolicExpr> two_x_plus_p;
            if (p == Rational(0)) {
                two_x_plus_p = two_x;
            } else {
                two_x_plus_p = SymbolicExpr::add(two_x, rd_num_rat(p));
            }
            auto inv_sqrt = SymbolicExpr::power(sqrt_s2, rd_num_int(-1));
            auto arctan_arg = SymbolicExpr::multiply(two_x_plus_p, inv_sqrt);
            auto atan_part = make_arctan(arctan_arg);
            auto two = rd_num_int(2);
            auto coeff_part = SymbolicExpr::multiply(
                SymbolicExpr::multiply(rd_num_rat(atan_coeff), two), inv_sqrt);
            auto term_atan = SymbolicExpr::multiply(coeff_part, atan_part);
            result = SymbolicExpr::add(result, term_atan);
        }

        return result;
    }

    // Should not reach here (factor degree > 2 is rejected by factor_denominator).
    return SymbolicExpr::number(0);
}

std::shared_ptr<SymbolicExpr> RationalDecompositionStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth) {
    (void)ctx;
    (void)depth;

    Polynomial<Rational> P, Q;
    try {
        if (!extract_rational(expr, var, P, Q)) {
            return nullptr; // not rational -> let next strategy try
        }
    } catch (...) {
        return nullptr;
    }

    // Defensive: if Q is zero or constant, this is not the right strategy.
    if (rd_is_zero_poly(Q) || Q.degree() < 1) return nullptr;

    // Don't fight with the simpler PartialFractionStrategy on degree <= 2:
    // letting RationalDecomposition handle them too is safe, but yields
    // bulkier outputs for trivial cases. We accept degree >= 3 here; for
    // degrees 1 and 2 we let later strategies (PartialFraction / direct
    // table) take over by returning nullptr.
    if (Q.degree() < 3) return nullptr;
    try {
        // Long division if needed.
        Polynomial<Rational> quot, rem;
        if (P.degree() >= Q.degree()) {
            poly_divide(P, Q, quot, rem);
        } else {
            quot = Polynomial<Rational>(Q.variable_name);
            rem = P;
        }

        // Factor denominator.
        std::vector<std::pair<Polynomial<Rational>, int>> factors;
        if (!factor_denominator(Q, factors)) {
            // Cannot factor over Q -> return unevaluated integral node.
            return Integrator::depends_on(expr, var)
                ? std::make_shared<SymbolicExpr>(
                      std::make_shared<FunctionNode>(
                          FunctionNode::FuncType::Calculus_Integral,
                          std::vector<std::shared_ptr<SymbolicNode>>{
                              expr.root,
                              SymbolicExpr::variable(var)->root}))
                : nullptr;
        }
        if (factors.empty()) {
            // Means Q is constant after factoring, so really there's nothing left;
            // integrate quot only.
            if (rd_is_zero_poly(quot)) return SymbolicExpr::number(0);
            // ∫ quot(x) dx = poly_integral(quot)
            // fall through; handled below
        }

        // Solve coefficients on the proper part rem / Q.
        std::vector<Polynomial<Rational>> numerators;
        if (!rd_is_zero_poly(rem)) {
            if (!solve_coefficients(rem, Q, factors, numerators)) {
                return std::make_shared<SymbolicExpr>(
                    std::make_shared<FunctionNode>(
                        FunctionNode::FuncType::Calculus_Integral,
                        std::vector<std::shared_ptr<SymbolicNode>>{
                            expr.root,
                            SymbolicExpr::variable(var)->root}));
            }
        }

        // Build the result piece by piece.
        std::shared_ptr<SymbolicExpr> result = SymbolicExpr::number(0);

        // Polynomial part from long division.
        if (!rd_is_zero_poly(quot)) {
            // Antiderivative of x^k is x^(k+1) / (k+1).
            for (size_t k = 0; k < quot.coeffs.size(); ++k) {
                if (quot.coeffs[k] == Rational(0)) continue;
                Rational coeff = quot.coeffs[k] / Rational(BigInt(static_cast<long long>(k + 1)));
                std::shared_ptr<SymbolicExpr> term;
                auto v = SymbolicExpr::variable(var);
                auto pw = SymbolicExpr::power(v, rd_num_int(static_cast<long long>(k + 1)));
                if (coeff == Rational(1)) {
                    term = pw;
                } else {
                    term = SymbolicExpr::multiply(rd_num_rat(coeff), pw);
                }
                result = SymbolicExpr::add(result, term);
            }
        }

        // Partial-fraction terms.
        if (!numerators.empty()) {
            size_t idx = 0;
            for (size_t i = 0; i < factors.size(); ++i) {
                const auto& [fpoly, mult] = factors[i];
                for (int l = 1; l <= mult; ++l, ++idx) {
                    if (idx >= numerators.size()) break;
                    auto term = integrate_term(numerators[idx], fpoly, l, var);
                    if (!term) {
                        // Fallback: unevaluated integral over the original.
                        return std::make_shared<SymbolicExpr>(
                            std::make_shared<FunctionNode>(
                                FunctionNode::FuncType::Calculus_Integral,
                                std::vector<std::shared_ptr<SymbolicNode>>{
                                    expr.root,
                                    SymbolicExpr::variable(var)->root}));
                    }
                    result = SymbolicExpr::add(result, term);
                }
            }
        }

        auto simplified = result->simplify();
        if (simplified) {
            // If the simplified form contains an unevaluated integral marker
            // (e.g. high-power irreducible quadratic), prefer keeping the
            // un-simplified result so downstream consumers can still read the
            // partial-fraction form.
            return simplified;
        }
        return result;
    } catch (...) {
        // Any unexpected runtime failure -> return unevaluated integral.
        return std::make_shared<SymbolicExpr>(
            std::make_shared<FunctionNode>(
                FunctionNode::FuncType::Calculus_Integral,
                std::vector<std::shared_ptr<SymbolicNode>>{
                    expr.root,
                    SymbolicExpr::variable(var)->root}));
    }
}

// ---------------------------------------------------------------
// SpecialFunctionStrategy
// ---------------------------------------------------------------
//
// Recognises a small, fixed set of integrand shapes whose antiderivatives are
// not elementary and instead must be expressed via the special functions erf,
// Ei, Si, Ci, Li. The inner argument of each pattern is required to be the
// integration variable itself (or a quadratic c*var^2 with rational c, for the
// erf branch). More general arguments are reduced first by SubstitutionStrategy
// or LinearSubstitutionStrategy upstream and so are not handled here.
//
// Patterns recognised (var = x):
//     exp(-x^2)      ->  (sqrt(pi) / 2) * erf(x)
//     exp(-a*x^2)    ->  (sqrt(pi) / (2*sqrt(a))) * erf(sqrt(a)*x)
//                        where a is a constant w.r.t. var (assumed positive).
//     exp(x) / x     ->  Ei(x)
//     sin(x) / x     ->  Si(x)
//     cos(x) / x     ->  Ci(x)
//     1 / ln(x)      ->  Li(x)
//
// Builders use FunctionNode with the dedicated FuncType enum values added in
// task 1.1, so that print_visitor / differentiation_visitor render and
// differentiate them correctly.

namespace {

// Build a single-argument FunctionNode wrapped in SymbolicExpr.
inline std::shared_ptr<SymbolicExpr> sf_make_fn(
    FunctionNode::FuncType t,
    const std::shared_ptr<SymbolicExpr>& arg) {
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(
            t, std::vector<std::shared_ptr<SymbolicNode>>{arg->root}));
}

// Build sqrt(pi).
inline std::shared_ptr<SymbolicExpr> sf_sqrt_pi() {
    return SymbolicExpr::sqrt(SymbolicExpr::variable("pi"));
}

// Test whether `node` is a single-argument FunctionNode of the given type
// whose argument is exactly the integration variable.
bool sf_is_fn_of_var(const std::shared_ptr<SymbolicNode>& node,
                     FunctionNode::FuncType t,
                     const std::string& var) {
    auto fn = std::dynamic_pointer_cast<FunctionNode>(node);
    if (!fn || fn->type != t) return false;
    if (fn->arguments.size() != 1) return false;
    auto v = std::dynamic_pointer_cast<VariableNode>(fn->arguments[0]);
    return v && v->name == var;
}

// Detect a 1/x factor: a PowerNode whose base is the integration variable and
// exponent is the integer -1.
bool sf_is_inv_var(const std::shared_ptr<SymbolicNode>& node,
                   const std::string& var) {
    auto pw = std::dynamic_pointer_cast<PowerNode>(node);
    if (!pw) return false;
    auto b = std::dynamic_pointer_cast<VariableNode>(pw->base);
    if (!b || b->name != var) return false;
    auto en = std::dynamic_pointer_cast<NumberNode>(pw->exponent);
    if (!en) return false;
    if (std::holds_alternative<BigInt>(en->value)) {
        return std::get<BigInt>(en->value) == BigInt(-1);
    }
    if (std::holds_alternative<Rational>(en->value)) {
        return std::get<Rational>(en->value) == Rational(-1);
    }
    if (std::holds_alternative<lmmc_real_t>(en->value)) {
        lmmc_real_t d = std::get<lmmc_real_t>(en->value);
        int eq = 0;
        lmmc_double_nearly_equal_tol(d, -1.0, 1e-12, 1e-12, &eq);
        return eq != 0;
    }
    return false;
}

// Detect a 1/ln(var) factor, i.e. PowerNode(ln(var), -1).
bool sf_is_inv_ln_var(const std::shared_ptr<SymbolicNode>& node,
                      const std::string& var) {
    auto pw = std::dynamic_pointer_cast<PowerNode>(node);
    if (!pw) return false;
    if (!sf_is_fn_of_var(pw->base, FunctionNode::FuncType::Ln, var)) return false;
    auto en = std::dynamic_pointer_cast<NumberNode>(pw->exponent);
    if (!en) return false;
    if (std::holds_alternative<BigInt>(en->value)) {
        return std::get<BigInt>(en->value) == BigInt(-1);
    }
    if (std::holds_alternative<Rational>(en->value)) {
        return std::get<Rational>(en->value) == Rational(-1);
    }
    if (std::holds_alternative<lmmc_real_t>(en->value)) {
        lmmc_real_t d = std::get<lmmc_real_t>(en->value);
        int eq = 0;
        lmmc_double_nearly_equal_tol(d, -1.0, 1e-12, 1e-12, &eq);
        return eq != 0;
    }
    return false;
}

// Split a node into (factors, has_inv_var) where the inv-var factor is removed
// from `factors` if present. Returns false if there is more than one inv-var
// factor (which would be 1/x^2 and is not the form we handle here).
bool sf_split_inv_var(const std::shared_ptr<SymbolicNode>& node,
                      const std::string& var,
                      std::vector<std::shared_ptr<SymbolicNode>>& other_factors,
                      bool& has_inv_var) {
    other_factors.clear();
    has_inv_var = false;
    std::vector<std::shared_ptr<SymbolicNode>> factors;
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        factors = mul->operands;
    } else {
        factors.push_back(node);
    }
    for (const auto& f : factors) {
        if (sf_is_inv_var(f, var)) {
            if (has_inv_var) return false; // multiple 1/x factors (1/x^2 form)
            has_inv_var = true;
        } else {
            other_factors.push_back(f);
        }
    }
    return true;
}

// Detect exp(-c*var^2) where c is a constant w.r.t. var that is rational and
// strictly positive. On success, sets `c_out` to that rational coefficient.
//
// We accept exp(arg) where the argument is a polynomial in var of degree
// exactly 2, no constant or linear term, and the leading rational coefficient
// is strictly negative. (We then return its absolute value as c.)
bool sf_match_exp_neg_quad(const std::shared_ptr<SymbolicNode>& node,
                           const std::string& var,
                           Rational& c_out) {
    auto fn = std::dynamic_pointer_cast<FunctionNode>(node);
    if (!fn || fn->type != FunctionNode::FuncType::Exp) return false;
    if (fn->arguments.size() != 1) return false;
    SymbolicExpr arg(fn->arguments[0]);

    Polynomial<Rational> poly;
    try {
        poly = symbolic_to_poly<Rational>(std::make_shared<SymbolicExpr>(arg), var);
    } catch (...) {
        return false;
    }
    if (poly.degree() != 2) return false;
    if (poly.coeffs.size() < 3) return false;
    // Constant and linear terms must be zero.
    if (!(poly.coeffs[0] == Rational(0))) return false;
    if (!(poly.coeffs[1] == Rational(0))) return false;
    Rational a2 = poly.coeffs[2];
    // Need a strictly negative coefficient: exp(-c*x^2) with c > 0.
    if (!(a2 < Rational(0))) return false;
    c_out = Rational(0) - a2; // c = -a2 > 0
    return true;
}

} // anonymous namespace

std::shared_ptr<SymbolicExpr> SpecialFunctionStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx, int depth) {
    (void)ctx;
    (void)depth;

    using FT = FunctionNode::FuncType;

    auto v = SymbolicExpr::variable(var);

    // ---- Pattern 1: 1/ln(x) -> Li(x) -------------------------------
    if (sf_is_inv_ln_var(expr.root, var)) {
        auto li = sf_make_fn(FT::Li, v);
        auto simp = li->simplify();
        return simp ? simp : li;
    }

    // ---- Pattern 2: exp(-x^2) or exp(-a*x^2) -> erf -----------------
    {
        Rational c_rat;
        if (sf_match_exp_neg_quad(expr.root, var, c_rat)) {
            // Result: (sqrt(pi) / (2 * sqrt(c))) * erf(sqrt(c) * x)
            auto sqrt_pi = sf_sqrt_pi();
            auto sqrt_c = SymbolicExpr::sqrt(SymbolicExpr::number(c_rat));
            auto scaled_x = SymbolicExpr::multiply(sqrt_c, v);
            auto erf_node = sf_make_fn(FT::Erf, scaled_x);

            auto two = SymbolicExpr::number(2);
            auto two_sqrt_c = SymbolicExpr::multiply(two, sqrt_c);
            auto coeff = SymbolicExpr::multiply(
                sqrt_pi,
                SymbolicExpr::power(two_sqrt_c, SymbolicExpr::number(-1)));
            auto result = SymbolicExpr::multiply(coeff, erf_node);
            auto simp = result->simplify();
            return simp ? simp : result;
        }
    }

    // ---- Patterns 3-5: exp(x)/x, sin(x)/x, cos(x)/x ----------------
    // These all have shape (something)*1/x, where the "something" is a
    // FunctionNode of var. Other 1/x patterns (e.g. 1/x alone, x*1/x) are
    // not our concern.
    {
        std::vector<std::shared_ptr<SymbolicNode>> others;
        bool has_inv = false;
        if (sf_split_inv_var(expr.root, var, others, has_inv) && has_inv && others.size() == 1) {
            const auto& other = others[0];
            if (sf_is_fn_of_var(other, FT::Exp, var)) {
                auto ei = sf_make_fn(FT::Ei, v);
                auto simp = ei->simplify();
                return simp ? simp : ei;
            }
            if (sf_is_fn_of_var(other, FT::Sin, var)) {
                auto si = sf_make_fn(FT::Si, v);
                auto simp = si->simplify();
                return simp ? simp : si;
            }
            if (sf_is_fn_of_var(other, FT::Cos, var)) {
                auto ci = sf_make_fn(FT::Ci, v);
                auto simp = ci->simplify();
                return simp ? simp : ci;
            }
        }
    }

    return nullptr;
}

// ---------------------------------------------------------------
// MultipleIntegralEngine
// ---------------------------------------------------------------
//
// Evaluates iterated (multiple) integrals by sequentially applying
// `Integrator::integrate` (for indefinite steps) or
// `Integrator::integrate_def` (for definite steps), going from the innermost
// step (index 0) to the outermost step (last index).
//
// If the integrand does not depend on a particular variable that has
// definite bounds [a, b], we short-circuit the integration step to a
// straightforward multiplication by (b - a) so we do not turn a constant
// expression into a fresh formal integral.
//
// Whenever a single step yields an expression that still contains an
// unevaluated `Calculus_Integral` node, we stop further iteration and
// return the partially-integrated result.

bool MultipleIntegralEngine::validate(const std::vector<IntegrationStep>& steps) const {
    if (steps.empty() || steps.size() > 3) return false;

    for (size_t i = 0; i < steps.size(); ++i) {
        const auto& s = steps[i];
        if (s.variable.empty()) return false;

        // Bounds must be either both null (indefinite) or both set (definite).
        const bool has_lower = static_cast<bool>(s.lower);
        const bool has_upper = static_cast<bool>(s.upper);
        if (has_lower != has_upper) return false;

        // Reject duplicate variables.
        for (size_t j = 0; j < i; ++j) {
            if (steps[j].variable == s.variable) return false;
        }
    }
    return true;
}

std::shared_ptr<SymbolicExpr> MultipleIntegralEngine::evaluate(
    const SymbolicExpr& integrand,
    const std::vector<IntegrationStep>& steps,
    Integrator& integrator) {

    if (!validate(steps)) return nullptr;

    auto current = std::make_shared<SymbolicExpr>(integrand);

    for (const auto& step : steps) {
        if (!current) return nullptr;

        // If a previous step already produced an unevaluated integral node,
        // stop and propagate that partial result.
        if (has_integral_node_check(current->root)) {
            return current;
        }

        const bool definite = static_cast<bool>(step.lower) && static_cast<bool>(step.upper);

        if (!definite) {
            // Indefinite single-variable integration.
            SymbolicExpr res = integrator.integrate(*current, step.variable);
            current = std::make_shared<SymbolicExpr>(res);
        } else {
            // Definite single-variable integration. If the integrand is
            // independent of the variable, short-circuit to current*(upper-lower).
            if (!Integrator::depends_on(*current, step.variable)) {
                auto diff = sym_sub(*step.upper, *step.lower);
                auto product = SymbolicExpr::multiply(current, diff);
                auto simp = product->simplify();
                current = simp ? simp : product;
            } else {
                SymbolicExpr res = integrator.integrate_def(
                    *current, step.variable, *step.lower, *step.upper);
                auto res_ptr = std::make_shared<SymbolicExpr>(res);
                auto simp = res_ptr->simplify();
                current = simp ? simp : res_ptr;
            }
        }
    }

    return current;
}

}
