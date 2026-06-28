/**
 * @file limit_visitor.hpp
 * @brief 极限访问器，通过代入和 L'Hôpital 法则计算极限。
 */
#pragma once

#include "../symbolic_ast.hpp"
#include "../assumption_context.hpp"
#include "normalization_visitor.hpp"
#include "differentiation_visitor.hpp"
#include <iostream>
#include <cmath>
#include <optional>
#include <limits>
#include <utility>

#ifndef LAMINA_API
#ifdef _WIN32
#ifdef LAMINA_CORE_EXPORTS
#define LAMINA_API __declspec(dllexport)
#else
#define LAMINA_API __declspec(dllimport)
#endif
#else
#define LAMINA_API
#endif
#endif

class SymbolicExpr;

/** @brief 极限访问器，通过代入求值和 L'Hôpital 法则计算符号表达式的极限 */
class LAMINA_API LimitVisitor : public SymbolicVisitor {
    std::string var;
    std::shared_ptr<SymbolicNode> point;
    std::string direction;
    const lamina::AssumptionContext* assumption_ctx_ = nullptr;
    int lhopital_depth_ = 0;
    static constexpr int max_lhopital_depth = 5;

    enum class IndeterminateForm { None, ZeroTimesInf, InfMinusInf, OnePowInf, ZeroPowZero, InfPowZero };

    // --- Low-level helpers ---
    bool is_inf(const std::shared_ptr<SymbolicNode>& node) const;
    bool is_neg_inf(const std::shared_ptr<SymbolicNode>& node) const;
    std::optional<int> get_node_sign(const std::shared_ptr<SymbolicNode>& node) const;
    double get_numeric_value(const std::shared_ptr<NumberNode>& num) const;
    double get_point_value() const;
    bool is_bounded(const std::shared_ptr<SymbolicNode>& node) const;
    bool is_bounded_expression(const std::shared_ptr<SymbolicNode>& node) const;
    bool tends_to_zero(const std::shared_ptr<SymbolicNode>& node) const;
    std::shared_ptr<SymbolicNode> eval_limit(const std::shared_ptr<SymbolicNode>& expr);

    // --- Classification ---
    IndeterminateForm classify_product_form(const std::vector<std::shared_ptr<SymbolicNode>>& vals);
    IndeterminateForm classify_power_form(const std::shared_ptr<SymbolicNode>& bv, const std::shared_ptr<SymbolicNode>& ev);
    IndeterminateForm classify_add_form(const std::vector<std::shared_ptr<SymbolicNode>>& vals);

    // --- Resolution ---
    std::shared_ptr<SymbolicNode> try_squeeze(const std::shared_ptr<SymbolicNode>& expr);
    std::pair<std::shared_ptr<SymbolicNode>, std::shared_ptr<SymbolicNode>> extract_num_den(const std::shared_ptr<SymbolicNode>& expr);
    std::shared_ptr<SymbolicNode> resolve_zero_times_inf(const std::vector<std::shared_ptr<SymbolicNode>>& factors, const std::vector<std::shared_ptr<SymbolicNode>>& factor_vals);
    std::shared_ptr<SymbolicNode> resolve_inf_minus_inf(AddNode& node, const std::vector<std::shared_ptr<SymbolicNode>>& operand_vals);
    std::shared_ptr<SymbolicNode> resolve_exponential_form(const std::shared_ptr<SymbolicNode>& base, const std::shared_ptr<SymbolicNode>& exponent);
    std::shared_ptr<SymbolicNode> apply_lhopital(const std::shared_ptr<SymbolicNode>& num, const std::shared_ptr<SymbolicNode>& den);

    // --- Taylor fallback (defined in limit_visitor.cpp) ---
    std::shared_ptr<SymbolicNode> taylor_fallback(const std::shared_ptr<SymbolicNode>& num, const std::shared_ptr<SymbolicNode>& den, int max_order = 8);
    std::shared_ptr<SymbolicNode> simplify_and_eval_ratio(const std::shared_ptr<SymbolicNode>& ratio_node);
    std::pair<std::shared_ptr<SymbolicNode>, int> find_leading_term(const std::shared_ptr<SymbolicExpr>& series_expr, const std::string& expand_var, const std::shared_ptr<SymbolicExpr>& expand_point, int max_order);
    static int get_sign(const std::shared_ptr<SymbolicNode>& node);

    // --- Limits at infinity helpers ---
    bool is_limit_at_infinity() const;
    bool is_limit_at_neg_infinity() const;
    int get_polynomial_degree(const std::shared_ptr<SymbolicNode>& node) const;
    std::shared_ptr<SymbolicNode> get_leading_coefficient(const std::shared_ptr<SymbolicNode>& node) const;
    std::shared_ptr<SymbolicNode> limit_rational_at_infinity(const std::shared_ptr<SymbolicNode>& num, const std::shared_ptr<SymbolicNode>& den);
    enum class GrowthClass { Constant, Logarithmic, Polynomial, Exponential, Unknown };
    GrowthClass classify_growth(const std::shared_ptr<SymbolicNode>& node) const;
    int get_growth_polynomial_degree(const std::shared_ptr<SymbolicNode>& node) const;
    std::shared_ptr<SymbolicNode> limit_by_growth_comparison(const std::shared_ptr<SymbolicNode>& num, const std::shared_ptr<SymbolicNode>& den);
    std::shared_ptr<SymbolicNode> handle_neg_infinity_limit(const std::shared_ptr<SymbolicNode>& expr);
    std::shared_ptr<SymbolicNode> substitute_neg_t(const std::shared_ptr<SymbolicNode>& node, const std::string& t_var) const;

    // --- Direction helpers ---
    int determine_sign_near_point(const std::shared_ptr<SymbolicNode>& expr, const std::string& dir);
    std::shared_ptr<SymbolicNode> select_branch_by_direction(PiecewiseNode& node, const std::string& dir);
    bool condition_satisfied_by_direction(const std::shared_ptr<SymbolicNode>& condition, const std::string& dir);
    std::optional<int> evaluate_relational_sign(const std::shared_ptr<RelationalNode>& rel, const std::string& dir);
    std::optional<std::shared_ptr<SymbolicNode>> evaluate_sgn_limit(const std::shared_ptr<SymbolicNode>& arg);
    std::optional<std::shared_ptr<SymbolicNode>> evaluate_abs_limit(const std::shared_ptr<SymbolicNode>& arg);

public:
    std::shared_ptr<SymbolicNode> result;

    LimitVisitor(std::string v, std::shared_ptr<SymbolicNode> p, std::string dir = "",
                 const lamina::AssumptionContext* ctx = nullptr)
        : var(std::move(v)), point(std::move(p)), direction(std::move(dir)), assumption_ctx_(ctx) {}

    std::shared_ptr<SymbolicNode> get_result() const { return result; }

    void visit(NumberNode& node) override { result = node.clone(); }
    void visit(VariableNode& node) override { if (node.name == var) result = point->clone(); else result = node.clone(); }
    void visit(AddNode& node) override;
    void visit(MultiplyNode& node) override;
    void visit(PowerNode& node) override;
    void visit(FunctionNode& node) override;
    void visit(MatrixNode& node) override { result = node.clone(); }
    void visit(PiecewiseNode& node) override;
    void visit(ComplexNode& node) override {
        node.real->accept(*this);
        auto lr = result;
        node.imag->accept(*this);
        auto li = result;
        result = SymbolicFactory::create_complex(lr, li);
    }
};

// ============================================================================
// Inline implementations
// ============================================================================

inline bool LimitVisitor::is_inf(const std::shared_ptr<SymbolicNode>& node) const {
    if (!node) return false;
    if (auto f = std::dynamic_pointer_cast<FunctionNode>(node)) return f->type == FunctionNode::FuncType::Infinity;
    if (auto m = std::dynamic_pointer_cast<MultiplyNode>(node)) for (auto& op : m->operands) if (is_inf(op)) return true;
    if (auto p = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (is_inf(p->base)) return true;
        // 0^(negative) is infinity
        if (p->base && p->base->is_zero()) {
            if (auto e_num = std::dynamic_pointer_cast<NumberNode>(p->exponent)) {
                double ev = 0;
                if (std::holds_alternative<double>(e_num->value)) ev = std::get<double>(e_num->value);
                else if (std::holds_alternative<BigInt>(e_num->value)) ev = std::get<BigInt>(e_num->value).to_double();
                else if (std::holds_alternative<Rational>(e_num->value)) ev = std::get<Rational>(e_num->value).to_double();
                if (ev < 0) return true;
            }
        }
    }
    if (auto a = std::dynamic_pointer_cast<AddNode>(node)) for (auto& op : a->operands) if (is_inf(op)) return true;
    return false;
}

inline bool LimitVisitor::is_neg_inf(const std::shared_ptr<SymbolicNode>& node) const {
    if (!node) return false;
    if (auto m = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        bool has_inf = false, has_neg = false;
        for (auto& op : m->operands) {
            if (is_inf(op)) has_inf = true;
            if (auto n = std::dynamic_pointer_cast<NumberNode>(op)) { auto s = get_node_sign(n); if (s && *s < 0) has_neg = true; }
        }
        return has_inf && has_neg;
    }
    return false;
}

inline std::optional<int> LimitVisitor::get_node_sign(const std::shared_ptr<SymbolicNode>& node) const {
    if (!node) return std::nullopt;
    if (node->is_zero()) return 0;
    if (auto num = std::dynamic_pointer_cast<NumberNode>(node)) {
        if (std::holds_alternative<double>(num->value)) { double v = std::get<double>(num->value); return v > 0 ? 1 : (v < 0 ? -1 : 0); }
        if (std::holds_alternative<BigInt>(num->value)) { BigInt v = std::get<BigInt>(num->value); return v > BigInt(0) ? 1 : (v < BigInt(0) ? -1 : 0); }
        if (std::holds_alternative<Rational>(num->value)) { Rational v = std::get<Rational>(num->value); return v > Rational(0) ? 1 : (v < Rational(0) ? -1 : 0); }
    }
    if (is_inf(node)) {
        if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) { for (auto& op : mul->operands) { if (auto n = std::dynamic_pointer_cast<NumberNode>(op)) { auto s = get_node_sign(n); if (s && *s < 0) return -1; } } }
        return 1;
    }
    return std::nullopt;
}

inline double LimitVisitor::get_numeric_value(const std::shared_ptr<NumberNode>& num) const {
    if (std::holds_alternative<double>(num->value)) return std::get<double>(num->value);
    if (std::holds_alternative<BigInt>(num->value)) return std::get<BigInt>(num->value).to_double();
    if (std::holds_alternative<Rational>(num->value)) return std::get<Rational>(num->value).to_double();
    return std::numeric_limits<double>::quiet_NaN();
}

inline double LimitVisitor::get_point_value() const {
    auto num = std::dynamic_pointer_cast<NumberNode>(point);
    if (!num) return std::numeric_limits<double>::quiet_NaN();
    return get_numeric_value(num);
}

inline bool LimitVisitor::is_bounded(const std::shared_ptr<SymbolicNode>& node) const {
    if (!node) return false;
    // sin(expr) and cos(expr) are bounded in [-1, 1]
    // arctan(expr) is bounded in (-π/2, π/2)
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(node)) {
        return func->type == FunctionNode::FuncType::Sin
            || func->type == FunctionNode::FuncType::Cos
            || func->type == FunctionNode::FuncType::ArcTan;
    }
    // Product of bounded functions is bounded
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        bool all_bounded_or_const = true;
        bool has_bounded = false;
        for (const auto& op : mul->operands) {
            if (is_bounded(op)) { has_bounded = true; }
            else if (std::dynamic_pointer_cast<NumberNode>(op)) { /* constants are fine */ }
            else { all_bounded_or_const = false; break; }
        }
        return all_bounded_or_const && has_bounded;
    }
    // bounded^(even positive integer) is bounded (e.g., sin(x)^2)
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (is_bounded(pow->base)) {
            if (auto num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                double e = 0;
                if (std::holds_alternative<double>(num->value)) e = std::get<double>(num->value);
                else if (std::holds_alternative<BigInt>(num->value)) e = std::get<BigInt>(num->value).to_double();
                else if (std::holds_alternative<Rational>(num->value)) e = std::get<Rational>(num->value).to_double();
                if (e > 0) return true;
            }
        }
    }
    return false;
}

/**
 * @brief 判断表达式是否有界（扩展检测）。
 *
 * 检测形如 constant + bounded 或 constant * bounded 的表达式，
 * 这些表达式虽然不一定在 [-1,1] 内，但仍然有界。
 * 用于夹逼定理的一般情形：若 f 有界且 g→0，则 f·g→0。
 *
 * 也检测 bounded × zero_tending 形式的表达式，这类表达式趋向 0，
 * 因此在极限点附近有界。
 */
inline bool LimitVisitor::is_bounded_expression(const std::shared_ptr<SymbolicNode>& node) const {
    if (!node) return false;
    // Already bounded by the basic check
    if (is_bounded(node)) return true;
    // constant + bounded is bounded (e.g., 2 + sin(1/x))
    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        bool has_bounded = false;
        for (const auto& op : add->operands) {
            if (is_bounded(op) || is_bounded_expression(op)) {
                has_bounded = true;
            } else if (std::dynamic_pointer_cast<NumberNode>(op)) {
                // constants are bounded
            } else {
                return false;
            }
        }
        return has_bounded;
    }
    return false;
}

inline bool LimitVisitor::tends_to_zero(const std::shared_ptr<SymbolicNode>& node) const {
    LimitVisitor sub_vis(var, point, direction, assumption_ctx_);
    auto node_copy = node->clone();
    node_copy->accept(sub_vis);
    auto val = sub_vis.get_result();
    if (!val) return false;
    NormalizationVisitor norm; val->accept(norm); val = norm.get_result();
    return val && val->is_zero();
}

inline std::shared_ptr<SymbolicNode> LimitVisitor::eval_limit(const std::shared_ptr<SymbolicNode>& expr) {
    LimitVisitor sub(var, point, direction, assumption_ctx_);
    sub.lhopital_depth_ = this->lhopital_depth_;
    expr->accept(sub);
    auto r = sub.get_result();
    if (r) { NormalizationVisitor norm; r->accept(norm); return norm.get_result(); }
    return r;
}

inline LimitVisitor::IndeterminateForm LimitVisitor::classify_product_form(const std::vector<std::shared_ptr<SymbolicNode>>& vals) {
    bool has_zero = false, has_inf_flag = false;
    for (auto& v : vals) { if (v && v->is_zero()) has_zero = true; if (is_inf(v)) has_inf_flag = true; }
    return (has_zero && has_inf_flag) ? IndeterminateForm::ZeroTimesInf : IndeterminateForm::None;
}

inline LimitVisitor::IndeterminateForm LimitVisitor::classify_power_form(const std::shared_ptr<SymbolicNode>& bv, const std::shared_ptr<SymbolicNode>& ev) {
    if (bv && bv->is_one() && is_inf(ev)) return IndeterminateForm::OnePowInf;
    if (bv && bv->is_zero() && ev && ev->is_zero()) return IndeterminateForm::ZeroPowZero;
    if (is_inf(bv) && ev && ev->is_zero()) return IndeterminateForm::InfPowZero;
    return IndeterminateForm::None;
}

inline LimitVisitor::IndeterminateForm LimitVisitor::classify_add_form(const std::vector<std::shared_ptr<SymbolicNode>>& vals) {
    bool has_pos_inf = false, has_neg_inf = false;
    for (auto& v : vals) { if (is_inf(v) && !is_neg_inf(v)) has_pos_inf = true; if (is_neg_inf(v)) has_neg_inf = true; }
    return (has_pos_inf && has_neg_inf) ? IndeterminateForm::InfMinusInf : IndeterminateForm::None;
}

inline std::shared_ptr<SymbolicNode> LimitVisitor::try_squeeze(const std::shared_ptr<SymbolicNode>& expr) {
    // Case 1: MultiplyNode — product of bounded × zero-tending → 0
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr)) {
        std::vector<std::shared_ptr<SymbolicNode>> bounded_factors, other_factors;
        for (const auto& op : mul->operands) {
            if (is_bounded(op)) {
                bounded_factors.push_back(op);
            } else if (is_bounded_expression(op)) {
                bounded_factors.push_back(op);
            } else {
                other_factors.push_back(op);
            }
        }
        if (bounded_factors.empty() || other_factors.empty()) return nullptr;

        // Check if the non-bounded factors tend to zero
        auto remaining = other_factors.size() == 1
            ? other_factors[0]
            : std::static_pointer_cast<SymbolicNode>(std::make_shared<MultiplyNode>(other_factors));
        if (tends_to_zero(remaining)) {
            return std::make_shared<NumberNode>(BigInt(0));
        }
        return nullptr;
    }

    // Case 2: AddNode — general squeeze theorem.
    // If f(x) = g(x) + h(x) where h(x) is a product of bounded × zero-tending,
    // then lim f = lim g (since h squeezes to 0).
    // This handles the case: f bounded between g and (g + bounded×zero) where
    // lim(lower) = lim(upper) = lim g = L.
    if (auto add = std::dynamic_pointer_cast<AddNode>(expr)) {
        std::vector<std::shared_ptr<SymbolicNode>> squeeze_to_zero_terms;
        std::vector<std::shared_ptr<SymbolicNode>> other_terms;

        for (const auto& op : add->operands) {
            // Check if this term is a product of bounded × zero-tending
            if (auto term_mul = std::dynamic_pointer_cast<MultiplyNode>(op)) {
                auto term_squeeze = try_squeeze(op);
                if (term_squeeze && term_squeeze->is_zero()) {
                    squeeze_to_zero_terms.push_back(op);
                    continue;
                }
            }
            other_terms.push_back(op);
        }

        // If we found squeeze-to-zero terms and have remaining terms with computable limits
        if (!squeeze_to_zero_terms.empty() && !other_terms.empty()) {
            auto remaining_expr = other_terms.size() == 1
                ? other_terms[0]
                : std::static_pointer_cast<SymbolicNode>(std::make_shared<AddNode>(other_terms));
            // Compute the limit of the remaining terms
            auto remaining_limit = eval_limit(remaining_expr);
            if (remaining_limit && !is_inf(remaining_limit)) {
                return remaining_limit;
            }
        }
    }

    return nullptr;
}

inline std::pair<std::shared_ptr<SymbolicNode>, std::shared_ptr<SymbolicNode>> LimitVisitor::extract_num_den(const std::shared_ptr<SymbolicNode>& expr) {
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr)) {
        std::vector<std::shared_ptr<SymbolicNode>> np, dp;
        for (auto& op : mul->operands) {
            if (auto pow = std::dynamic_pointer_cast<PowerNode>(op)) {
                if (auto nn = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                    double e = 0;
                    if (std::holds_alternative<double>(nn->value)) e = std::get<double>(nn->value);
                    else if (std::holds_alternative<BigInt>(nn->value)) e = std::get<BigInt>(nn->value).to_double();
                    else if (std::holds_alternative<Rational>(nn->value)) e = std::get<Rational>(nn->value).to_double();
                    if (e < 0) { dp.push_back(std::make_shared<PowerNode>(pow->base, std::make_shared<NumberNode>(-e))); continue; }
                }
            }
            np.push_back(op);
        }
        if (!dp.empty()) {
            auto n = np.empty() ? std::static_pointer_cast<SymbolicNode>(std::make_shared<NumberNode>(BigInt(1))) : (np.size() == 1 ? np[0] : std::static_pointer_cast<SymbolicNode>(std::make_shared<MultiplyNode>(np)));
            auto d = dp.size() == 1 ? dp[0] : std::static_pointer_cast<SymbolicNode>(std::make_shared<MultiplyNode>(dp));
            return {n, d};
        }
    }
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(expr)) {
        if (auto nn = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
            double e = 0;
            if (std::holds_alternative<double>(nn->value)) e = std::get<double>(nn->value);
            else if (std::holds_alternative<BigInt>(nn->value)) e = std::get<BigInt>(nn->value).to_double();
            else if (std::holds_alternative<Rational>(nn->value)) e = std::get<Rational>(nn->value).to_double();
            if (e < 0) return {std::make_shared<NumberNode>(BigInt(1)), std::make_shared<PowerNode>(pow->base, std::make_shared<NumberNode>(-e))};
        }
    }
    return {expr, std::make_shared<NumberNode>(BigInt(1))};
}

inline std::shared_ptr<SymbolicNode> LimitVisitor::resolve_zero_times_inf(const std::vector<std::shared_ptr<SymbolicNode>>& factors, const std::vector<std::shared_ptr<SymbolicNode>>& factor_vals) {
    std::vector<std::shared_ptr<SymbolicNode>> zero_f, inf_f, other_f;
    for (size_t i = 0; i < factors.size(); ++i) {
        if (factor_vals[i] && factor_vals[i]->is_zero()) zero_f.push_back(factors[i]);
        else if (is_inf(factor_vals[i])) inf_f.push_back(factors[i]);
        else other_f.push_back(factors[i]);
    }
    if (zero_f.empty() || inf_f.empty()) return nullptr;
    auto f = zero_f.size() == 1 ? zero_f[0] : std::static_pointer_cast<SymbolicNode>(std::make_shared<MultiplyNode>(zero_f));
    auto g = inf_f.size() == 1 ? inf_f[0] : std::static_pointer_cast<SymbolicNode>(std::make_shared<MultiplyNode>(inf_f));

    // Choose between 0/0 form (f/(1/g)) and ∞/∞ form (g/(1/f)).
    // Heuristic: if the zero factor has exponential growth class (e.g., e^(-x)),
    // prefer ∞/∞ form because L'Hôpital on 0/0 won't converge (polynomial grows).
    // Otherwise, try 0/0 first as it's simpler for most cases.
    bool prefer_inf_over_inf = false;
    if (is_limit_at_infinity() || is_limit_at_neg_infinity()) {
        GrowthClass f_growth = classify_growth(f);
        GrowthClass g_growth = classify_growth(g);
        // If zero factor is exponential (e.g., e^(-x)→0) and inf factor is polynomial,
        // the ∞/∞ form (polynomial/exponential) converges in one step.
        if (f_growth == GrowthClass::Exponential && g_growth == GrowthClass::Polynomial) {
            prefer_inf_over_inf = true;
        }
        // If zero factor is polynomial (e.g., 1/x→0) and inf factor is logarithmic,
        // the 0/0 form works well.
    }

    std::shared_ptr<SymbolicNode> res = nullptr;

    if (prefer_inf_over_inf) {
        // Try ∞/∞ form first: g / (1/f)
        auto f_inv = std::make_shared<PowerNode>(f, std::make_shared<NumberNode>(BigInt(-1)));
        res = apply_lhopital(g, f_inv);
        if (!res) {
            auto g_inv = std::make_shared<PowerNode>(g, std::make_shared<NumberNode>(BigInt(-1)));
            res = apply_lhopital(f, g_inv);
        }
    } else {
        // Try 0/0 form first: f / (1/g)
        auto g_inv = std::make_shared<PowerNode>(g, std::make_shared<NumberNode>(BigInt(-1)));
        res = apply_lhopital(f, g_inv);
        if (!res) {
            auto f_inv = std::make_shared<PowerNode>(f, std::make_shared<NumberNode>(BigInt(-1)));
            res = apply_lhopital(g, f_inv);
        }
    }

    if (res && !other_f.empty()) {
        std::vector<std::shared_ptr<SymbolicNode>> ff; ff.push_back(res);
        for (auto& of : other_f) { auto ov = eval_limit(of); if (ov) ff.push_back(ov); }
        if (ff.size() == 1) return ff[0];
        auto prod = std::make_shared<MultiplyNode>(ff);
        NormalizationVisitor norm; prod->accept(norm); return norm.get_result();
    }
    return res;
}

inline std::shared_ptr<SymbolicNode> LimitVisitor::resolve_inf_minus_inf(AddNode& node, const std::vector<std::shared_ptr<SymbolicNode>>&) {
    std::vector<std::shared_ptr<SymbolicNode>> nums, dens;
    bool has_nontrivial_den = false;
    for (auto& op : node.operands) { auto nd = extract_num_den(op); nums.push_back(nd.first); dens.push_back(nd.second); if (!nd.second->is_one()) has_nontrivial_den = true; }

    // If no operand has a denominator, combining into a fraction won't help.
    // Fall back to Taylor expansion at the limit point.
    if (!has_nontrivial_den) {
        // Try Taylor fallback: treat the whole expression as numerator / 1
        auto whole_expr = node.operands.size() == 1 ? node.operands[0] : std::static_pointer_cast<SymbolicNode>(std::make_shared<AddNode>(node.operands));
        auto one_node = std::make_shared<NumberNode>(BigInt(1));
        auto tf_result = taylor_fallback(whole_expr, one_node);
        return tf_result;
    }

    auto common_den = dens.size() == 1 ? dens[0] : std::static_pointer_cast<SymbolicNode>(std::make_shared<MultiplyNode>(dens));
    std::vector<std::shared_ptr<SymbolicNode>> new_num_terms;
    for (size_t i = 0; i < nums.size(); ++i) {
        std::vector<std::shared_ptr<SymbolicNode>> tf; tf.push_back(nums[i]);
        for (size_t j = 0; j < dens.size(); ++j) { if (j != i) tf.push_back(dens[j]); }
        new_num_terms.push_back(tf.size() == 1 ? tf[0] : std::static_pointer_cast<SymbolicNode>(std::make_shared<MultiplyNode>(tf)));
    }
    auto combined_num = new_num_terms.size() == 1 ? new_num_terms[0] : std::static_pointer_cast<SymbolicNode>(std::make_shared<AddNode>(new_num_terms));
    NormalizationVisitor norm;
    combined_num->accept(norm); combined_num = norm.get_result();
    common_den->accept(norm); common_den = norm.get_result();
    auto nv = eval_limit(combined_num); auto dv = eval_limit(common_den);
    bool nz = nv && nv->is_zero(), dz = dv && dv->is_zero();
    if ((nz && dz) || (is_inf(nv) && is_inf(dv))) return apply_lhopital(combined_num, common_den);
    if (!dz && !is_inf(dv) && dv) {
        auto ratio = std::make_shared<MultiplyNode>(std::vector<std::shared_ptr<SymbolicNode>>{combined_num, std::make_shared<PowerNode>(common_den, std::make_shared<NumberNode>(BigInt(-1)))});
        return eval_limit(ratio);
    }
    return nullptr;
}

inline std::shared_ptr<SymbolicNode> LimitVisitor::resolve_exponential_form(const std::shared_ptr<SymbolicNode>& base, const std::shared_ptr<SymbolicNode>& exponent) {
    std::vector<std::shared_ptr<SymbolicNode>> ln_args = {base};
    auto ln_base = std::make_shared<FunctionNode>(FunctionNode::FuncType::Ln, ln_args);
    std::vector<std::shared_ptr<SymbolicNode>> prod_ops = {exponent, ln_base};
    auto product = std::make_shared<MultiplyNode>(prod_ops);
    auto inner_limit = eval_limit(product);
    if (!inner_limit) return nullptr;
    if (!is_inf(inner_limit)) {
        std::vector<std::shared_ptr<SymbolicNode>> exp_args = {inner_limit};
        auto exp_result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Exp, exp_args);
        NormalizationVisitor norm; exp_result->accept(norm); return norm.get_result();
    }
    if (is_inf(inner_limit) && !is_neg_inf(inner_limit)) { std::vector<std::shared_ptr<SymbolicNode>> inf_args; return std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args); }
    if (is_neg_inf(inner_limit)) return std::make_shared<NumberNode>(BigInt(0));
    return nullptr;
}

inline std::shared_ptr<SymbolicNode> LimitVisitor::apply_lhopital(const std::shared_ptr<SymbolicNode>& num, const std::shared_ptr<SymbolicNode>& den) {
    if (lhopital_depth_ >= max_lhopital_depth) return taylor_fallback(num, den);
    DifferentiationVisitor diff_vis(var);
    num->accept(diff_vis); auto dN = diff_vis.get_result();
    den->accept(diff_vis); auto dD = diff_vis.get_result();
    if (!dN || !dD) return taylor_fallback(num, den);
    NormalizationVisitor norm;
    dN->accept(norm); dN = norm.get_result();
    dD->accept(norm); dD = norm.get_result();
    if (dD->is_zero()) return taylor_fallback(num, den);

    // Try to simplify the ratio dN/dD algebraically before evaluating limits.
    // This handles cases like x*ln(1+1/x) where L'Hôpital produces derivatives
    // with common factors (e.g., x^-2) that cancel, avoiding infinite 0/0 loops.
    {
        std::vector<std::shared_ptr<SymbolicNode>> ratio_ops = {
            dN, std::make_shared<PowerNode>(dD, std::make_shared<NumberNode>(BigInt(-1)))};
        auto ratio_node = std::make_shared<MultiplyNode>(ratio_ops);
        // Wrap in SymbolicExpr for simplification (defined in limit_visitor.cpp)
        auto simp_result = simplify_and_eval_ratio(ratio_node);
        if (simp_result) return simp_result;
    }

    // Construct the ratio dN/dD. Instead of creating a MultiplyNode that might
    // be misinterpreted as 0×∞, evaluate as a proper fraction: compute limits
    // of dN and dD separately and check for 0/0 or ∞/∞ forms.
    LimitVisitor sub_n(var, point, direction, assumption_ctx_);
    sub_n.lhopital_depth_ = this->lhopital_depth_ + 1;
    dN->accept(sub_n);
    auto val_n = sub_n.get_result();
    if (val_n) { NormalizationVisitor n2; val_n->accept(n2); val_n = n2.get_result(); }

    LimitVisitor sub_d(var, point, direction, assumption_ctx_);
    sub_d.lhopital_depth_ = this->lhopital_depth_ + 1;
    dD->accept(sub_d);
    auto val_d = sub_d.get_result();
    if (val_d) { NormalizationVisitor n2; val_d->accept(n2); val_d = n2.get_result(); }

    if (!val_n || !val_d) return taylor_fallback(num, den);

    bool n_zero = val_n->is_zero();
    bool d_zero = val_d->is_zero();
    bool n_inf = is_inf(val_n);
    bool d_inf = is_inf(val_d);

    // If still indeterminate (0/0 or ∞/∞), recurse with increased depth
    if ((n_zero && d_zero) || (n_inf && d_inf)) {
        LimitVisitor sub(var, point, direction, assumption_ctx_);
        sub.lhopital_depth_ = this->lhopital_depth_ + 1;
        // Build ratio as MultiplyNode for the sub-visitor
        std::vector<std::shared_ptr<SymbolicNode>> ratio_ops = {dN, std::make_shared<PowerNode>(dD, std::make_shared<NumberNode>(BigInt(-1)))};
        auto ratio = std::make_shared<MultiplyNode>(ratio_ops);
        ratio->accept(sub);
        return sub.get_result();
    }

    // Denominator is zero but numerator isn't → ±∞
    if (d_zero && !n_zero) {
        int sign_n = 1;
        if (auto nn = std::dynamic_pointer_cast<NumberNode>(val_n)) {
            auto s = get_node_sign(nn);
            if (s) sign_n = *s;
        }
        std::vector<std::shared_ptr<SymbolicNode>> inf_args;
        auto inf_node = std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
        if (sign_n < 0) {
            std::vector<std::shared_ptr<SymbolicNode>> m = {std::make_shared<NumberNode>(BigInt(-1)), inf_node};
            return std::make_shared<MultiplyNode>(m);
        }
        return inf_node;
    }

    // Normal case: compute the ratio of limits
    if (!d_zero && !d_inf) {
        auto ratio_result = std::make_shared<MultiplyNode>(std::vector<std::shared_ptr<SymbolicNode>>{
            val_n, std::make_shared<PowerNode>(val_d, std::make_shared<NumberNode>(BigInt(-1)))});
        ratio_result->accept(norm);
        return norm.get_result();
    }

    // Numerator is finite, denominator is ∞ → 0
    if (!n_inf && !n_zero && d_inf) {
        return std::make_shared<NumberNode>(BigInt(0));
    }

    // Fallback: evaluate the ratio expression directly
    std::vector<std::shared_ptr<SymbolicNode>> ratio_ops = {dN, std::make_shared<PowerNode>(dD, std::make_shared<NumberNode>(BigInt(-1)))};
    auto ratio = std::make_shared<MultiplyNode>(ratio_ops);
    LimitVisitor sub(var, point, direction, assumption_ctx_);
    sub.lhopital_depth_ = this->lhopital_depth_ + 1;
    ratio->accept(sub);
    return sub.get_result();
}

inline int LimitVisitor::determine_sign_near_point(const std::shared_ptr<SymbolicNode>& expr, const std::string& dir) {
    if (dir.empty()) return 0;
    DifferentiationVisitor diff_vis(var);
    std::shared_ptr<SymbolicNode> curr = expr;
    for (int i = 1; i <= 3; ++i) {
        curr->accept(diff_vis); auto deriv = diff_vis.get_result(); if (!deriv) break;
        LimitVisitor sv2(var, point, dir, assumption_ctx_); deriv->accept(sv2); auto val = sv2.get_result(); if (!val) break;
        NormalizationVisitor norm; val->accept(norm); val = norm.get_result();
        if (!val->is_zero()) {
            int s = 1;
            if (auto num = std::dynamic_pointer_cast<NumberNode>(val)) {
                if (std::holds_alternative<double>(num->value)) s = std::get<double>(num->value) > 0 ? 1 : -1;
                else if (std::holds_alternative<BigInt>(num->value)) s = std::get<BigInt>(num->value) > BigInt(0) ? 1 : -1;
                else if (std::holds_alternative<Rational>(num->value)) s = std::get<Rational>(num->value) > Rational(0) ? 1 : -1;
            }
            if (i % 2 != 0 && dir == "-") s *= -1;
            return s;
        }
        curr = deriv;
    }
    return 0;
}

inline void LimitVisitor::visit(AddNode& node) {
    // Handle negative infinity: substitute x = -t and evaluate lim(t→+∞)
    if (is_limit_at_neg_infinity()) {
        auto neg_inf_result = handle_neg_infinity_limit(std::make_shared<AddNode>(node.operands));
        if (neg_inf_result) { result = neg_inf_result; return; }
    }
    // General squeeze theorem for AddNode: detect terms that are products of
    // bounded × zero-tending (squeeze to 0) and compute limit of remaining terms.
    auto squeeze_result = try_squeeze(std::make_shared<AddNode>(std::vector<std::shared_ptr<SymbolicNode>>(node.operands.begin(), node.operands.end())));
    if (squeeze_result) { result = squeeze_result; return; }
    std::vector<std::shared_ptr<SymbolicNode>> new_ops;
    for (auto& op : node.operands) { op->accept(*this); new_ops.push_back(result); }
    auto form = classify_add_form(new_ops);
    if (form == IndeterminateForm::InfMinusInf) { auto resolved = resolve_inf_minus_inf(node, new_ops); if (resolved) { result = resolved; return; } }
    NormalizationVisitor norm; std::make_shared<AddNode>(new_ops)->accept(norm); result = norm.get_result();
}

inline void LimitVisitor::visit(MultiplyNode& node) {
    // Handle negative infinity: substitute x = -t and evaluate lim(t→+∞)
    if (is_limit_at_neg_infinity()) {
        auto neg_inf_result = handle_neg_infinity_limit(std::make_shared<MultiplyNode>(std::vector<std::shared_ptr<SymbolicNode>>(node.operands.begin(), node.operands.end())));
        if (neg_inf_result) { result = neg_inf_result; return; }
    }
    auto squeeze_result = try_squeeze(std::make_shared<MultiplyNode>(std::vector<std::shared_ptr<SymbolicNode>>(node.operands.begin(), node.operands.end())));
    if (squeeze_result) { result = squeeze_result; return; }
    std::vector<std::shared_ptr<SymbolicNode>> subs_ops;
    for (auto& op : node.operands) { op->accept(*this); subs_ops.push_back(result); }
    NormalizationVisitor norm;
    auto subs_res = std::make_shared<MultiplyNode>(subs_ops); subs_res->accept(norm); auto final_subs = norm.get_result();
    std::vector<std::shared_ptr<SymbolicNode>> num_nodes, den_nodes;
    for (auto& op : node.operands) {
        if (auto pow = std::dynamic_pointer_cast<PowerNode>(op)) {
            if (auto num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                double e = 0;
                if (std::holds_alternative<double>(num->value)) e = std::get<double>(num->value);
                else if (std::holds_alternative<BigInt>(num->value)) e = std::get<BigInt>(num->value).to_double();
                else if (std::holds_alternative<Rational>(num->value)) e = std::get<Rational>(num->value).to_double();
                if (e < 0) { den_nodes.push_back(std::make_shared<PowerNode>(pow->base, std::make_shared<NumberNode>(-e))); continue; }
            }
        }
        num_nodes.push_back(op);
    }
    if (den_nodes.empty()) {
        auto prod_form = classify_product_form(subs_ops);
        if (prod_form == IndeterminateForm::ZeroTimesInf) { auto resolved = resolve_zero_times_inf(node.operands, subs_ops); if (resolved) { result = resolved; return; } }
        result = final_subs; return;
    }
    auto N = num_nodes.size() == 1 ? num_nodes[0] : std::static_pointer_cast<SymbolicNode>(std::make_shared<MultiplyNode>(num_nodes));
    auto D = den_nodes.size() == 1 ? den_nodes[0] : std::static_pointer_cast<SymbolicNode>(std::make_shared<MultiplyNode>(den_nodes));
    auto val_n = eval_limit(N); auto val_d = eval_limit(D);
    bool n_zero = val_n && val_n->is_zero(), d_zero = val_d && val_d->is_zero();
    bool n_inf = is_inf(val_n), d_inf = is_inf(val_d);
    // Limits at infinity: try rational function degree comparison and growth-rate comparison
    if (is_limit_at_infinity() || is_limit_at_neg_infinity()) {
        if (n_inf && d_inf) {
            auto rat_result = limit_rational_at_infinity(N, D);
            if (rat_result) { result = rat_result; return; }
            auto growth_result = limit_by_growth_comparison(N, D);
            if (growth_result) { result = growth_result; return; }
        }
    }
    if ((n_zero && d_zero) || (n_inf && d_inf)) { auto lhop = apply_lhopital(N, D); if (lhop) { result = lhop; return; } }
    if (!n_zero && d_zero) {
        int sign_n = 1;
        if (auto num = std::dynamic_pointer_cast<NumberNode>(val_n)) {
            if (std::holds_alternative<double>(num->value)) sign_n = std::get<double>(num->value) > 0 ? 1 : -1;
            else if (std::holds_alternative<BigInt>(num->value)) sign_n = std::get<BigInt>(num->value) > BigInt(0) ? 1 : -1;
            else if (std::holds_alternative<Rational>(num->value)) sign_n = std::get<Rational>(num->value) > Rational(0) ? 1 : -1;
        }
        int sign_d = (direction == "-" ? -1 : 1);
        int final_sign = sign_n * sign_d;
        std::vector<std::shared_ptr<SymbolicNode>> inf_args;
        auto inf_node = std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
        if (final_sign < 0) { std::vector<std::shared_ptr<SymbolicNode>> m_args = {std::make_shared<NumberNode>(BigInt(-1)), inf_node}; result = std::make_shared<MultiplyNode>(m_args); }
        else { result = inf_node; }
        return;
    }
    result = final_subs;
}

inline void LimitVisitor::visit(PowerNode& node) {
    // Handle negative infinity: substitute x = -t and evaluate lim(t→+∞)
    if (is_limit_at_neg_infinity()) {
        auto neg_inf_result = handle_neg_infinity_limit(std::make_shared<PowerNode>(node.base, node.exponent));
        if (neg_inf_result) { result = neg_inf_result; return; }
    }
    node.base->accept(*this); auto b = result;
    node.exponent->accept(*this); auto e = result;
    auto form = classify_power_form(b, e);
    if (form != IndeterminateForm::None) { auto resolved = resolve_exponential_form(node.base, node.exponent); if (resolved) { result = resolved; return; } }
    // Handle ∞^(negative) → 0 and ∞^(positive) → ∞
    if (is_inf(b)) {
        if (auto e_num = std::dynamic_pointer_cast<NumberNode>(e)) {
            double ev = 0;
            if (std::holds_alternative<double>(e_num->value)) ev = std::get<double>(e_num->value);
            else if (std::holds_alternative<BigInt>(e_num->value)) ev = std::get<BigInt>(e_num->value).to_double();
            else if (std::holds_alternative<Rational>(e_num->value)) ev = std::get<Rational>(e_num->value).to_double();
            if (ev < 0) { result = std::make_shared<NumberNode>(BigInt(0)); return; }
            if (ev > 0) {
                std::vector<std::shared_ptr<SymbolicNode>> inf_args;
                result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
                return;
            }
        }
    }
    // Handle 0^(negative) → ±∞ with direction-aware sign determination
    if (b && b->is_zero()) {
        if (auto e_num = std::dynamic_pointer_cast<NumberNode>(e)) {
            double ev = 0;
            if (std::holds_alternative<double>(e_num->value)) ev = std::get<double>(e_num->value);
            else if (std::holds_alternative<BigInt>(e_num->value)) ev = std::get<BigInt>(e_num->value).to_double();
            else if (std::holds_alternative<Rational>(e_num->value)) ev = std::get<Rational>(e_num->value).to_double();
            if (ev < 0 && !direction.empty()) {
                // base→0, exponent is negative, direction specified → ±∞
                // Determine sign based on approach direction and exponent parity
                std::vector<std::shared_ptr<SymbolicNode>> inf_args;
                auto inf_node = std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
                int exp_int = static_cast<int>(ev);
                bool odd_exponent = (exp_int % 2 != 0);
                if (odd_exponent) {
                    // For odd negative exponents (e.g., x^(-1), x^(-3)):
                    // x→0+ gives +∞, x→0- gives -∞
                    int base_sign = determine_sign_near_point(node.base, direction);
                    if (base_sign < 0) {
                        std::vector<std::shared_ptr<SymbolicNode>> m = {std::make_shared<NumberNode>(BigInt(-1)), inf_node};
                        result = std::make_shared<MultiplyNode>(m);
                        return;
                    }
                }
                // Even negative exponents always give +∞
                result = inf_node;
                return;
            }
        }
    }
    NormalizationVisitor norm; std::make_shared<PowerNode>(b, e)->accept(norm); result = norm.get_result();
}

inline void LimitVisitor::visit(FunctionNode& node) {
    // Handle negative infinity: substitute x = -t and evaluate lim(t→+∞)
    if (is_limit_at_neg_infinity() && node.type != FunctionNode::FuncType::Infinity) {
        auto neg_inf_result = handle_neg_infinity_limit(std::make_shared<FunctionNode>(node.type, node.arguments));
        if (neg_inf_result) { result = neg_inf_result; return; }
    }
    // 方向感知的符号函数处理
    if (node.type == FunctionNode::FuncType::Sgn && node.arguments.size() == 1) {
        auto r = evaluate_sgn_limit(node.arguments[0]);
        if (r) { result = *r; return; }
    }
    // 方向感知的绝对值处理
    if (node.type == FunctionNode::FuncType::Abs && node.arguments.size() == 1) {
        auto r = evaluate_abs_limit(node.arguments[0]);
        if (r) { result = *r; return; }
    }
    // Composed functions at infinity: evaluate inner limit first, then apply outer function
    if (is_limit_at_infinity() && node.arguments.size() == 1 && node.type != FunctionNode::FuncType::Infinity) {
        auto inner_lim = eval_limit(node.arguments[0]);
        if (inner_lim) {
            // If inner limit is a finite number, evaluate the function at that value
            if (auto num = std::dynamic_pointer_cast<NumberNode>(inner_lim)) {
                std::vector<std::shared_ptr<SymbolicNode>> new_args = {inner_lim};
                NormalizationVisitor norm;
                std::make_shared<FunctionNode>(node.type, new_args)->accept(norm);
                result = norm.get_result();
                return;
            }
            // If inner limit is infinity, handle based on function type
            if (is_inf(inner_lim)) {
                bool neg = is_neg_inf(inner_lim);
                switch (node.type) {
                    case FunctionNode::FuncType::Exp:
                        if (neg) { result = std::make_shared<NumberNode>(BigInt(0)); }
                        else { std::vector<std::shared_ptr<SymbolicNode>> inf_args; result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args); }
                        return;
                    case FunctionNode::FuncType::Ln:
                        if (!neg) { std::vector<std::shared_ptr<SymbolicNode>> inf_args; result = std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args); }
                        return;
                    case FunctionNode::FuncType::ArcTan:
                        if (neg) { result = std::make_shared<NumberNode>(Rational(-1, 2)); result = std::make_shared<MultiplyNode>(std::vector<std::shared_ptr<SymbolicNode>>{result, std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, std::vector<std::shared_ptr<SymbolicNode>>{})}); }
                        else { result = std::make_shared<NumberNode>(Rational(1, 2)); result = std::make_shared<MultiplyNode>(std::vector<std::shared_ptr<SymbolicNode>>{result, std::make_shared<FunctionNode>(FunctionNode::FuncType::Infinity, std::vector<std::shared_ptr<SymbolicNode>>{})}); }
                        return;
                    default:
                        break;
                }
            }
        }
    }
    std::vector<std::shared_ptr<SymbolicNode>> new_args;
    for (auto& a : node.arguments) { a->accept(*this); new_args.push_back(result); }
    NormalizationVisitor norm; std::make_shared<FunctionNode>(node.type, new_args)->accept(norm); result = norm.get_result();
}
