#include "visitors/limit_visitor.hpp"
#include "visitors/normalization_visitor.hpp"
#include "visitors/differentiation_visitor.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <utility>

int& LimitVisitor::active_limit_depth() {
    static thread_local int depth = 0;
    return depth;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::make_product_or_one(
    const std::vector<std::shared_ptr<const SymbolicNode>>& factors) {
    if (factors.empty()) return lamina::detail::make_node<NumberNode>(BigInt(1));
    if (factors.size() == 1) return factors[0];
    return lamina::detail::make_node<MultiplyNode>(factors);
}

LimitVisitor::EvaluationScope::EvaluationScope() {
    if (LimitVisitor::active_limit_depth() < LimitVisitor::max_active_limit_depth) {
        ++LimitVisitor::active_limit_depth();
        entered = true;
    }
}

LimitVisitor::EvaluationScope::~EvaluationScope() {
    if (entered) --LimitVisitor::active_limit_depth();
}

std::shared_ptr<const SymbolicNode> LimitVisitor::get_result() const { return result; }

void LimitVisitor::visit(const NumberNode& node) { result = node.clone(); }

void LimitVisitor::visit(const VariableNode& node) {
    result = node.name() == var ? point->clone() : node.clone();
}

void LimitVisitor::visit(const MatrixNode& node) { result = node.clone(); }

void LimitVisitor::visit(const RelationalNode& node) {
    node.left()->accept(*this);
    auto left = result;
    node.right()->accept(*this);
    result = lamina::detail::make_node<RelationalNode>(left, result, node.op());
}

void LimitVisitor::visit(const LogicalNode& node) {
    node.left()->accept(*this);
    auto left = result;
    std::shared_ptr<const SymbolicNode> right;
    if (node.right()) {
        node.right()->accept(*this);
        right = result;
    }
    result = lamina::detail::make_node<LogicalNode>(left, right, node.op());
}

void LimitVisitor::visit(const SummationNode& node) { result = node.clone(); }
void LimitVisitor::visit(const ProductNode& node) { result = node.clone(); }
void LimitVisitor::visit(const TransformNode& node) { result = node.clone(); }
void LimitVisitor::visit(const QuantifierNode& node) { result = node.clone(); }
void LimitVisitor::visit(const SetBuilderNode& node) { result = node.clone(); }
void LimitVisitor::visit(const FiniteSetNode& node) { result = node.clone(); }
void LimitVisitor::visit(const IntervalNode& node) { result = node.clone(); }
void LimitVisitor::visit(const MembershipNode& node) { result = node.clone(); }
void LimitVisitor::visit(const QuantityNode& node) {
    node.value()->accept(*this);
    result = lamina::detail::make_node<QuantityNode>(
        result, node.dimension(), node.scale_to_base(), node.display_unit());
}

void LimitVisitor::visit(const ComplexNode& node) {
    node.real()->accept(*this);
    auto real = result;
    node.imag()->accept(*this);
    result = SymbolicFactory::create_complex(real, result);
}

bool LimitVisitor::is_inf(const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) return false;
    if (auto f = std::dynamic_pointer_cast<const FunctionNode>(node)) return f->type() == FunctionNode::FuncType::Infinity;
    if (auto m = std::dynamic_pointer_cast<const MultiplyNode>(node)) for (auto& op : m->operands()) if (is_inf(op)) return true;
    if (auto p = std::dynamic_pointer_cast<const PowerNode>(node)) {
        if (is_inf(p->base())) return true;
        /// 0^(negative) is infinity
        if (p->base() && p->base()->is_zero()) {
            if (auto e_num = std::dynamic_pointer_cast<const NumberNode>(p->exponent())) {
                double ev = 0;
                if (std::holds_alternative<double>(e_num->value())) ev = std::get<double>(e_num->value());
                else if (std::holds_alternative<BigInt>(e_num->value())) ev = std::get<BigInt>(e_num->value()).to_double();
                else if (std::holds_alternative<Rational>(e_num->value())) ev = std::get<Rational>(e_num->value()).to_double();
                if (ev < 0) return true;
            }
        }
    }
    if (auto a = std::dynamic_pointer_cast<const AddNode>(node)) for (auto& op : a->operands()) if (is_inf(op)) return true;
    return false;
}

bool LimitVisitor::is_neg_inf(const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) return false;
    if (auto m = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        bool has_inf = false, has_neg = false;
        for (auto& op : m->operands()) {
            if (is_inf(op)) has_inf = true;
            if (auto n = std::dynamic_pointer_cast<const NumberNode>(op)) { auto s = get_node_sign(n); if (s && *s < 0) has_neg = true; }
        }
        return has_inf && has_neg;
    }
    return false;
}

std::optional<int> LimitVisitor::get_node_sign(const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) return std::nullopt;
    if (node->is_zero()) return 0;
    if (auto num = std::dynamic_pointer_cast<const NumberNode>(node)) {
        if (std::holds_alternative<double>(num->value())) { double v = std::get<double>(num->value()); return v > 0 ? 1 : (v < 0 ? -1 : 0); }
        if (std::holds_alternative<BigInt>(num->value())) { BigInt v = std::get<BigInt>(num->value()); return v > BigInt(0) ? 1 : (v < BigInt(0) ? -1 : 0); }
        if (std::holds_alternative<Rational>(num->value())) { Rational v = std::get<Rational>(num->value()); return v > Rational(0) ? 1 : (v < Rational(0) ? -1 : 0); }
    }
    if (is_inf(node)) {
        if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) { for (auto& op : mul->operands()) { if (auto n = std::dynamic_pointer_cast<const NumberNode>(op)) { auto s = get_node_sign(n); if (s && *s < 0) return -1; } } }
        return 1;
    }
    return std::nullopt;
}

double LimitVisitor::get_numeric_value(const std::shared_ptr<const NumberNode>& num) const {
    if (std::holds_alternative<double>(num->value())) return std::get<double>(num->value());
    if (std::holds_alternative<BigInt>(num->value())) return std::get<BigInt>(num->value()).to_double();
    if (std::holds_alternative<Rational>(num->value())) return std::get<Rational>(num->value()).to_double();
    return std::numeric_limits<double>::quiet_NaN();
}

double LimitVisitor::get_point_value() const {
    auto num = std::dynamic_pointer_cast<const NumberNode>(point);
    if (!num) return std::numeric_limits<double>::quiet_NaN();
    return get_numeric_value(num);
}

bool LimitVisitor::is_bounded(const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) return false;
    /// sin(expr) and cos(expr) are bounded in [-1, 1]
    /// arctan(expr) is bounded in (-π/2, π/2)
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return func->type() == FunctionNode::FuncType::Sin
            || func->type() == FunctionNode::FuncType::Cos
            || func->type() == FunctionNode::FuncType::ArcTan;
    }
    /// Product of bounded functions is bounded
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        bool all_bounded_or_const = true;
        bool has_bounded = false;
        for (const auto& op : mul->operands()) {
            if (is_bounded(op)) { has_bounded = true; }
            else if (std::dynamic_pointer_cast<const NumberNode>(op)) { /* constants are fine */ }
            else { all_bounded_or_const = false; break; }
        }
        return all_bounded_or_const && has_bounded;
    }
    /// bounded^(even positive integer) is bounded (e.g., sin(x)^2)
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        if (is_bounded(pow->base())) {
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
                double e = 0;
                if (std::holds_alternative<double>(num->value())) e = std::get<double>(num->value());
                else if (std::holds_alternative<BigInt>(num->value())) e = std::get<BigInt>(num->value()).to_double();
                else if (std::holds_alternative<Rational>(num->value())) e = std::get<Rational>(num->value()).to_double();
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
bool LimitVisitor::is_bounded_expression(const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) return false;
    /// Already bounded by the basic check
    if (is_bounded(node)) return true;
    /// constant + bounded is bounded (e.g., 2 + sin(1/x))
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        bool has_bounded = false;
        for (const auto& op : add->operands()) {
            if (is_bounded(op) || is_bounded_expression(op)) {
                has_bounded = true;
            } else if (std::dynamic_pointer_cast<const NumberNode>(op)) {
                /// constants are bounded
            } else {
                return false;
            }
        }
        return has_bounded;
    }
    return false;
}

bool LimitVisitor::tends_to_zero(const std::shared_ptr<const SymbolicNode>& node) const {
    EvaluationScope scope;
    if (!scope) return false;
    LimitVisitor sub_vis(var, point, direction, assumption_ctx_);
    auto node_copy = node->clone();
    node_copy->accept(sub_vis);
    auto val = sub_vis.get_result();
    if (!val) return false;
    NormalizationVisitor norm; val->accept(norm); val = norm.get_result();
    return val && val->is_zero();
}

std::shared_ptr<const SymbolicNode> LimitVisitor::eval_limit(const std::shared_ptr<const SymbolicNode>& expr) {
    EvaluationScope scope;
    if (!scope) return nullptr;
    LimitVisitor sub(var, point, direction, assumption_ctx_);
    sub.lhopital_depth_ = this->lhopital_depth_;
    expr->accept(sub);
    auto r = sub.get_result();
    if (r) { NormalizationVisitor norm; r->accept(norm); return norm.get_result(); }
    return r;
}

LimitVisitor::IndeterminateForm LimitVisitor::classify_product_form(const std::vector<std::shared_ptr<const SymbolicNode>>& vals) {
    bool has_zero = false, has_inf_flag = false;
    for (auto& v : vals) { if (v && v->is_zero()) has_zero = true; if (is_inf(v)) has_inf_flag = true; }
    return (has_zero && has_inf_flag) ? IndeterminateForm::ZeroTimesInf : IndeterminateForm::None;
}

LimitVisitor::IndeterminateForm LimitVisitor::classify_power_form(const std::shared_ptr<const SymbolicNode>& bv, const std::shared_ptr<const SymbolicNode>& ev) {
    if (bv && bv->is_one() && is_inf(ev)) return IndeterminateForm::OnePowInf;
    if (bv && bv->is_zero() && ev && ev->is_zero()) return IndeterminateForm::ZeroPowZero;
    if (is_inf(bv) && ev && ev->is_zero()) return IndeterminateForm::InfPowZero;
    return IndeterminateForm::None;
}

LimitVisitor::IndeterminateForm LimitVisitor::classify_add_form(const std::vector<std::shared_ptr<const SymbolicNode>>& vals) {
    bool has_pos_inf = false, has_neg_inf = false;
    for (auto& v : vals) { if (is_inf(v) && !is_neg_inf(v)) has_pos_inf = true; if (is_neg_inf(v)) has_neg_inf = true; }
    return (has_pos_inf && has_neg_inf) ? IndeterminateForm::InfMinusInf : IndeterminateForm::None;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::try_squeeze(const std::shared_ptr<const SymbolicNode>& expr) {
    /// Case 1: MultiplyNode — product of bounded × zero-tending → 0
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(expr)) {
        std::vector<std::shared_ptr<const SymbolicNode>> bounded_factors, other_factors;
        for (const auto& op : mul->operands()) {
            if (is_bounded(op)) {
                bounded_factors.push_back(op);
            } else if (is_bounded_expression(op)) {
                bounded_factors.push_back(op);
            } else {
                other_factors.push_back(op);
            }
        }
        if (bounded_factors.empty() || other_factors.empty()) return nullptr;

        /// Check if the non-bounded factors tend to zero
        auto remaining = other_factors.size() == 1
            ? other_factors[0]
            : make_product_or_one(other_factors);
        if (tends_to_zero(remaining)) {
            return lamina::detail::make_node<NumberNode>(BigInt(0));
        }
        return nullptr;
    }

    /// Case 2: AddNode — general squeeze theorem.
    /// If f(x) = g(x) + h(x) where h(x) is a product of bounded × zero-tending,
    /// then lim f = lim g (since h squeezes to 0).
    /// This handles the case: f bounded between g and (g + bounded×zero) where
    /// lim(lower) = lim(upper) = lim g = L.
    if (auto add = std::dynamic_pointer_cast<const AddNode>(expr)) {
        std::vector<std::shared_ptr<const SymbolicNode>> squeeze_to_zero_terms;
        std::vector<std::shared_ptr<const SymbolicNode>> other_terms;

        for (const auto& op : add->operands()) {
            /// Check if this term is a product of bounded × zero-tending
            if (auto term_mul = std::dynamic_pointer_cast<const MultiplyNode>(op)) {
                auto term_squeeze = try_squeeze(op);
                if (term_squeeze && term_squeeze->is_zero()) {
                    squeeze_to_zero_terms.push_back(op);
                    continue;
                }
            }
            other_terms.push_back(op);
        }

        /// If we found squeeze-to-zero terms and have remaining terms with computable limits
        if (!squeeze_to_zero_terms.empty() && !other_terms.empty()) {
            auto remaining_expr = other_terms.size() == 1
                ? other_terms[0]
                : std::static_pointer_cast<const SymbolicNode>(lamina::detail::make_node<AddNode>(other_terms));
            /// Compute the limit of the remaining terms
            auto remaining_limit = eval_limit(remaining_expr);
            if (remaining_limit && !is_inf(remaining_limit)) {
                return remaining_limit;
            }
        }
    }

    return nullptr;
}

std::pair<std::shared_ptr<const SymbolicNode>, std::shared_ptr<const SymbolicNode>> LimitVisitor::extract_num_den(const std::shared_ptr<const SymbolicNode>& expr) {
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(expr)) {
        std::vector<std::shared_ptr<const SymbolicNode>> np, dp;
        for (auto& op : mul->operands()) {
            if (auto pow = std::dynamic_pointer_cast<const PowerNode>(op)) {
                if (auto nn = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
                    double e = 0;
                    if (std::holds_alternative<double>(nn->value())) e = std::get<double>(nn->value());
                    else if (std::holds_alternative<BigInt>(nn->value())) e = std::get<BigInt>(nn->value()).to_double();
                    else if (std::holds_alternative<Rational>(nn->value())) e = std::get<Rational>(nn->value()).to_double();
                    if (e < 0) { dp.push_back(lamina::detail::make_node<PowerNode>(pow->base(), lamina::detail::make_node<NumberNode>(-e))); continue; }
                }
            }
            np.push_back(op);
        }
        if (!dp.empty()) {
            auto n = make_product_or_one(np);
            auto d = make_product_or_one(dp);
            return {n, d};
        }
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(expr)) {
        if (auto nn = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
            double e = 0;
            if (std::holds_alternative<double>(nn->value())) e = std::get<double>(nn->value());
            else if (std::holds_alternative<BigInt>(nn->value())) e = std::get<BigInt>(nn->value()).to_double();
            else if (std::holds_alternative<Rational>(nn->value())) e = std::get<Rational>(nn->value()).to_double();
            if (e < 0) return {lamina::detail::make_node<NumberNode>(BigInt(1)), lamina::detail::make_node<PowerNode>(pow->base(), lamina::detail::make_node<NumberNode>(-e))};
        }
    }
    return {expr, lamina::detail::make_node<NumberNode>(BigInt(1))};
}

std::shared_ptr<const SymbolicNode> LimitVisitor::resolve_zero_times_inf(const std::vector<std::shared_ptr<const SymbolicNode>>& factors, const std::vector<std::shared_ptr<const SymbolicNode>>& factor_vals) {
    std::vector<std::shared_ptr<const SymbolicNode>> zero_f, inf_f, other_f;
    for (size_t i = 0; i < factors.size(); ++i) {
        if (factor_vals[i] && factor_vals[i]->is_zero()) zero_f.push_back(factors[i]);
        else if (is_inf(factor_vals[i])) inf_f.push_back(factors[i]);
        else other_f.push_back(factors[i]);
    }
    if (zero_f.empty() || inf_f.empty()) return nullptr;
    auto f = make_product_or_one(zero_f);
    auto g = make_product_or_one(inf_f);

    /// Choose between 0/0 form (f/(1/g)) and ∞/∞ form (g/(1/f)).
    /// Heuristic: if the zero factor has exponential growth class (e.g., e^(-x)),
    /// prefer ∞/∞ form because L'Hôpital on 0/0 won't converge (polynomial grows).
    /// Otherwise, try 0/0 first as it's simpler for most cases.
    bool prefer_inf_over_inf = false;
    if (is_limit_at_infinity() || is_limit_at_neg_infinity()) {
        GrowthClass f_growth = classify_growth(f);
        GrowthClass g_growth = classify_growth(g);
        /// If zero factor is exponential (e.g., e^(-x)→0) and inf factor is polynomial,
        /// the ∞/∞ form (polynomial/exponential) converges in one step.
        if (f_growth == GrowthClass::Exponential && g_growth == GrowthClass::Polynomial) {
            prefer_inf_over_inf = true;
        }
        /// If zero factor is polynomial (e.g., 1/x→0) and inf factor is logarithmic,
        /// the 0/0 form works well.
    }

    std::shared_ptr<const SymbolicNode> res = nullptr;

    if (prefer_inf_over_inf) {
        /// Try ∞/∞ form first: g / (1/f)
        auto f_inv = lamina::detail::make_node<PowerNode>(f, lamina::detail::make_node<NumberNode>(BigInt(-1)));
        res = apply_lhopital(g, f_inv);
        if (!res) {
            auto g_inv = lamina::detail::make_node<PowerNode>(g, lamina::detail::make_node<NumberNode>(BigInt(-1)));
            res = apply_lhopital(f, g_inv);
        }
    } else {
        /// Try 0/0 form first: f / (1/g)
        auto g_inv = lamina::detail::make_node<PowerNode>(g, lamina::detail::make_node<NumberNode>(BigInt(-1)));
        res = apply_lhopital(f, g_inv);
        if (!res) {
            auto f_inv = lamina::detail::make_node<PowerNode>(f, lamina::detail::make_node<NumberNode>(BigInt(-1)));
            res = apply_lhopital(g, f_inv);
        }
    }

    if (res && !other_f.empty()) {
        std::vector<std::shared_ptr<const SymbolicNode>> ff; ff.push_back(res);
        for (auto& of : other_f) { auto ov = eval_limit(of); if (ov) ff.push_back(ov); }
        if (ff.size() == 1) return ff[0];
        auto prod = lamina::detail::make_node<MultiplyNode>(ff);
        NormalizationVisitor norm; prod->accept(norm); return norm.get_result();
    }
    return res;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::resolve_inf_minus_inf(const AddNode& node, const std::vector<std::shared_ptr<const SymbolicNode>>&) {
    std::vector<std::shared_ptr<const SymbolicNode>> nums, dens;
    bool has_nontrivial_den = false;
    for (auto& op : node.operands()) { auto nd = extract_num_den(op); nums.push_back(nd.first); dens.push_back(nd.second); if (!nd.second->is_one()) has_nontrivial_den = true; }

    /// If no operand has a denominator, combining into a fraction won't help.
    /// Fall back to Taylor expansion at the limit point.
    if (!has_nontrivial_den) {
        /// Try Taylor fallback: treat the whole expression as numerator / 1
        auto whole_expr = node.operands().size() == 1 ? node.operands()[0] : std::static_pointer_cast<const SymbolicNode>(lamina::detail::make_node<AddNode>(node.operands()));
        auto one_node = lamina::detail::make_node<NumberNode>(BigInt(1));
        auto tf_result = taylor_fallback(whole_expr, one_node);
        return tf_result;
    }

    auto common_den = make_product_or_one(dens);
    std::vector<std::shared_ptr<const SymbolicNode>> new_num_terms;
    for (size_t i = 0; i < nums.size(); ++i) {
        std::vector<std::shared_ptr<const SymbolicNode>> tf; tf.push_back(nums[i]);
        for (size_t j = 0; j < dens.size(); ++j) { if (j != i) tf.push_back(dens[j]); }
        new_num_terms.push_back(make_product_or_one(tf));
    }
    auto combined_num = new_num_terms.size() == 1 ? new_num_terms[0] : std::static_pointer_cast<const SymbolicNode>(lamina::detail::make_node<AddNode>(new_num_terms));
    NormalizationVisitor norm;
    combined_num->accept(norm); combined_num = norm.get_result();
    common_den->accept(norm); common_den = norm.get_result();
    auto nv = eval_limit(combined_num); auto dv = eval_limit(common_den);
    bool nz = nv && nv->is_zero(), dz = dv && dv->is_zero();
    if ((nz && dz) || (is_inf(nv) && is_inf(dv))) return apply_lhopital(combined_num, common_den);
    if (!dz && !is_inf(dv) && dv) {
        auto ratio = lamina::detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>{combined_num, lamina::detail::make_node<PowerNode>(common_den, lamina::detail::make_node<NumberNode>(BigInt(-1)))});
        return eval_limit(ratio);
    }
    return nullptr;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::resolve_exponential_form(const std::shared_ptr<const SymbolicNode>& base, const std::shared_ptr<const SymbolicNode>& exponent) {
    std::vector<std::shared_ptr<const SymbolicNode>> ln_args = {base};
    auto ln_base = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln, ln_args);
    std::vector<std::shared_ptr<const SymbolicNode>> prod_ops = {exponent, ln_base};
    auto product = lamina::detail::make_node<MultiplyNode>(prod_ops);
    auto inner_limit = eval_limit(product);
    if (!inner_limit) return nullptr;
    if (!is_inf(inner_limit)) {
        std::vector<std::shared_ptr<const SymbolicNode>> exp_args = {inner_limit};
        auto exp_result = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Exp, exp_args);
        NormalizationVisitor norm; exp_result->accept(norm); return norm.get_result();
    }
    if (is_inf(inner_limit) && !is_neg_inf(inner_limit)) { std::vector<std::shared_ptr<const SymbolicNode>> inf_args; return lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args); }
    if (is_neg_inf(inner_limit)) return lamina::detail::make_node<NumberNode>(BigInt(0));
    return nullptr;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::apply_lhopital(const std::shared_ptr<const SymbolicNode>& num, const std::shared_ptr<const SymbolicNode>& den) {
    EvaluationScope scope;
    if (!scope) return nullptr;
    if (lhopital_depth_ >= max_lhopital_depth) return taylor_fallback(num, den);
    DifferentiationVisitor diff_vis(var);
    num->accept(diff_vis); auto dN = diff_vis.get_result();
    den->accept(diff_vis); auto dD = diff_vis.get_result();
    if (!dN || !dD) return taylor_fallback(num, den);
    NormalizationVisitor norm;
    dN->accept(norm); dN = norm.get_result();
    dD->accept(norm); dD = norm.get_result();
    if (dD->is_zero()) return taylor_fallback(num, den);

    /// Try to simplify the ratio dN/dD algebraically before evaluating limits.
    /// This handles cases like x*ln(1+1/x) where L'Hôpital produces derivatives
    /// with common factors (e.g., x^-2) that cancel, avoiding infinite 0/0 loops.
    {
        std::vector<std::shared_ptr<const SymbolicNode>> ratio_ops = {
            dN, lamina::detail::make_node<PowerNode>(dD, lamina::detail::make_node<NumberNode>(BigInt(-1)))};
        auto ratio_node = lamina::detail::make_node<MultiplyNode>(ratio_ops);
        /// Wrap in SymbolicExpr for simplification (defined in limit_visitor.cpp)
        auto simp_result = simplify_and_eval_ratio(ratio_node);
        if (simp_result) return simp_result;
    }

    /// Construct the ratio dN/dD. Instead of creating a MultiplyNode that might
    /// be misinterpreted as 0×∞, evaluate as a proper fraction: compute limits
    /// of dN and dD separately and check for 0/0 or ∞/∞ forms.
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

    /// If still indeterminate (0/0 or ∞/∞), recurse with increased depth
    if ((n_zero && d_zero) || (n_inf && d_inf)) {
        LimitVisitor sub(var, point, direction, assumption_ctx_);
        sub.lhopital_depth_ = this->lhopital_depth_ + 1;
        /// Build ratio as MultiplyNode for the sub-visitor
        std::vector<std::shared_ptr<const SymbolicNode>> ratio_ops = {dN, lamina::detail::make_node<PowerNode>(dD, lamina::detail::make_node<NumberNode>(BigInt(-1)))};
        auto ratio = lamina::detail::make_node<MultiplyNode>(ratio_ops);
        ratio->accept(sub);
        return sub.get_result();
    }

    /// Denominator is zero but numerator isn't → ±∞
    if (d_zero && !n_zero) {
        int sign_n = 1;
        if (auto nn = std::dynamic_pointer_cast<const NumberNode>(val_n)) {
            auto s = get_node_sign(nn);
            if (s) sign_n = *s;
        }
        std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
        auto inf_node = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
        if (sign_n < 0) {
            std::vector<std::shared_ptr<const SymbolicNode>> m = {lamina::detail::make_node<NumberNode>(BigInt(-1)), inf_node};
            return lamina::detail::make_node<MultiplyNode>(m);
        }
        return inf_node;
    }

    /// Normal case: compute the ratio of limits
    if (!d_zero && !d_inf) {
        auto ratio_result = lamina::detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>{
            val_n, lamina::detail::make_node<PowerNode>(val_d, lamina::detail::make_node<NumberNode>(BigInt(-1)))});
        ratio_result->accept(norm);
        return norm.get_result();
    }

    /// Numerator is finite, denominator is ∞ → 0
    if (!n_inf && !n_zero && d_inf) {
        return lamina::detail::make_node<NumberNode>(BigInt(0));
    }

    /// Fallback: evaluate the ratio expression directly
    std::vector<std::shared_ptr<const SymbolicNode>> ratio_ops = {dN, lamina::detail::make_node<PowerNode>(dD, lamina::detail::make_node<NumberNode>(BigInt(-1)))};
    auto ratio = lamina::detail::make_node<MultiplyNode>(ratio_ops);
    LimitVisitor sub(var, point, direction, assumption_ctx_);
    sub.lhopital_depth_ = this->lhopital_depth_ + 1;
    ratio->accept(sub);
    return sub.get_result();
}

int LimitVisitor::determine_sign_near_point(const std::shared_ptr<const SymbolicNode>& expr, const std::string& dir) {
    if (dir.empty()) return 0;
    DifferentiationVisitor diff_vis(var);
    std::shared_ptr<const SymbolicNode> curr = expr;
    for (int i = 1; i <= 3; ++i) {
        curr->accept(diff_vis); auto deriv = diff_vis.get_result(); if (!deriv) break;
        LimitVisitor sv2(var, point, dir, assumption_ctx_); deriv->accept(sv2); auto val = sv2.get_result(); if (!val) break;
        NormalizationVisitor norm; val->accept(norm); val = norm.get_result();
        if (!val->is_zero()) {
            int s = 1;
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(val)) {
                if (std::holds_alternative<double>(num->value())) s = std::get<double>(num->value()) > 0 ? 1 : -1;
                else if (std::holds_alternative<BigInt>(num->value())) s = std::get<BigInt>(num->value()) > BigInt(0) ? 1 : -1;
                else if (std::holds_alternative<Rational>(num->value())) s = std::get<Rational>(num->value()) > Rational(0) ? 1 : -1;
            }
            if (i % 2 != 0 && dir == "-") s *= -1;
            return s;
        }
        curr = deriv;
    }
    return 0;
}

void LimitVisitor::visit(const AddNode& node) {
    /// Handle negative infinity: substitute x = -t and evaluate lim(t→+∞)
    if (is_limit_at_neg_infinity()) {
        auto neg_inf_result = handle_neg_infinity_limit(lamina::detail::make_node<AddNode>(node.operands()));
        if (neg_inf_result) { result = neg_inf_result; return; }
    }
    /// General squeeze theorem for AddNode: detect terms that are products of
    /// bounded × zero-tending (squeeze to 0) and compute limit of remaining terms.
    auto squeeze_result = try_squeeze(lamina::detail::make_node<AddNode>(std::vector<std::shared_ptr<const SymbolicNode>>(node.operands().begin(), node.operands().end())));
    if (squeeze_result) { result = squeeze_result; return; }
    std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
    for (auto& op : node.operands()) { op->accept(*this); new_ops.push_back(result); }
    auto form = classify_add_form(new_ops);
    if (form == IndeterminateForm::InfMinusInf) { auto resolved = resolve_inf_minus_inf(node, new_ops); if (resolved) { result = resolved; return; } }
    NormalizationVisitor norm; lamina::detail::make_node<AddNode>(new_ops)->accept(norm); result = norm.get_result();
}

void LimitVisitor::visit(const MultiplyNode& node) {
    /// Handle negative infinity: substitute x = -t and evaluate lim(t→+∞)
    if (is_limit_at_neg_infinity()) {
        auto neg_inf_result = handle_neg_infinity_limit(lamina::detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>(node.operands().begin(), node.operands().end())));
        if (neg_inf_result) { result = neg_inf_result; return; }
    }
    {
        auto point_num = std::dynamic_pointer_cast<const NumberNode>(point);
        bool at_zero_from_right = point_num && point_num->is_zero() && direction == "+";
        bool has_positive_var_power = false;
        bool has_log_var = false;

        auto is_positive_power_of_limit_var = [&](const std::shared_ptr<const SymbolicNode>& op) {
            if (auto variable = std::dynamic_pointer_cast<const VariableNode>(op)) {
                return variable->name() == var;
            }
            if (auto power = std::dynamic_pointer_cast<const PowerNode>(op)) {
                auto base_var = std::dynamic_pointer_cast<const VariableNode>(power->base());
                auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
                return base_var && base_var->name() == var && exponent && exponent->is_positive();
            }
            return false;
        };

        auto is_log_of_limit_var = [&](const std::shared_ptr<const SymbolicNode>& op) {
            auto function = std::dynamic_pointer_cast<const FunctionNode>(op);
            if (!function || function->arguments().size() != 1) return false;
            if (function->type() != FunctionNode::FuncType::Ln &&
                function->type() != FunctionNode::FuncType::Log) return false;
            auto argument_var = std::dynamic_pointer_cast<const VariableNode>(function->arguments()[0]);
            return argument_var && argument_var->name() == var;
        };

        if (at_zero_from_right) {
            for (const auto& op : node.operands()) {
                has_positive_var_power = has_positive_var_power || is_positive_power_of_limit_var(op);
                has_log_var = has_log_var || is_log_of_limit_var(op);
            }
            if (has_positive_var_power && has_log_var) {
                result = lamina::detail::make_node<NumberNode>(BigInt(0));
                return;
            }
        }
    }
    if (is_limit_at_infinity()) {
        for (const auto& op : node.operands()) {
            auto function = std::dynamic_pointer_cast<const FunctionNode>(op);
            if (!function || function->arguments().size() != 1 ||
                (function->type() != FunctionNode::FuncType::Ln &&
                 function->type() != FunctionNode::FuncType::Log)) {
                continue;
            }
            auto add = std::dynamic_pointer_cast<const AddNode>(function->arguments()[0]);
            if (!add) continue;

            std::shared_ptr<const SymbolicNode> small_num;
            bool has_one = false;
            for (const auto& add_op : add->operands()) {
                if (add_op->is_one()) {
                    has_one = true;
                    continue;
                }
                small_num = add_op;
            }
            if (!has_one || !small_num) continue;

            std::shared_ptr<const SymbolicNode> numerator = lamina::detail::make_node<NumberNode>(BigInt(1));
            bool denominator_is_var = false;
            if (auto power = std::dynamic_pointer_cast<const PowerNode>(small_num)) {
                auto base_var = std::dynamic_pointer_cast<const VariableNode>(power->base());
                auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
                if (base_var && base_var->name() == var && exponent) {
                    double exponent_value = 0.0;
                    if (std::holds_alternative<double>(exponent->value())) exponent_value = std::get<double>(exponent->value());
                    else if (std::holds_alternative<BigInt>(exponent->value())) exponent_value = std::get<BigInt>(exponent->value()).to_double();
                    else if (std::holds_alternative<Rational>(exponent->value())) exponent_value = std::get<Rational>(exponent->value()).to_double();
                    denominator_is_var = exponent_value == -1.0;
                }
            } else if (auto small_mul = std::dynamic_pointer_cast<const MultiplyNode>(small_num)) {
                std::vector<std::shared_ptr<const SymbolicNode>> numerator_parts;
                for (const auto& small_factor : small_mul->operands()) {
                    auto power = std::dynamic_pointer_cast<const PowerNode>(small_factor);
                    auto base_var = power ? std::dynamic_pointer_cast<const VariableNode>(power->base()) : nullptr;
                    auto exponent = power ? std::dynamic_pointer_cast<const NumberNode>(power->exponent()) : nullptr;
                    if (base_var && base_var->name() == var && exponent) {
                        double exponent_value = 0.0;
                        if (std::holds_alternative<double>(exponent->value())) exponent_value = std::get<double>(exponent->value());
                        else if (std::holds_alternative<BigInt>(exponent->value())) exponent_value = std::get<BigInt>(exponent->value()).to_double();
                        else if (std::holds_alternative<Rational>(exponent->value())) exponent_value = std::get<Rational>(exponent->value()).to_double();
                        if (exponent_value == -1.0) {
                            denominator_is_var = true;
                            continue;
                        }
                    }
                    numerator_parts.push_back(small_factor);
                }
                if (!numerator_parts.empty()) {
                    numerator = numerator_parts.size() == 1
                        ? numerator_parts[0]
                : make_product_or_one(numerator_parts);
                }
            }
            if (!denominator_is_var) continue;

            std::vector<std::shared_ptr<const SymbolicNode>> remaining_limits = {numerator};
            bool has_var_factor = false;
            for (const auto& remaining : node.operands()) {
                if (remaining.get() == op.get()) continue;
                if (auto variable = std::dynamic_pointer_cast<const VariableNode>(remaining)) {
                    if (variable->name() == var) {
                        has_var_factor = true;
                        continue;
                    }
                }
                auto remaining_limit = eval_limit(remaining);
                if (!remaining_limit) {
                    remaining_limits.clear();
                    break;
                }
                remaining_limits.push_back(remaining_limit);
            }
            if (!has_var_factor || remaining_limits.empty()) continue;
            auto reduced = remaining_limits.size() == 1
                ? remaining_limits[0]
                : make_product_or_one(remaining_limits);
            NormalizationVisitor norm;
            reduced->accept(norm);
            result = norm.get_result();
            return;
        }
        for (const auto& op : node.operands()) {
            auto function = std::dynamic_pointer_cast<const FunctionNode>(op);
            if (!function || function->arguments().size() != 1 ||
                function->type() != FunctionNode::FuncType::Exp) {
                continue;
            }

            auto exponent_limit = eval_limit(function->arguments()[0]);
            if (!is_neg_inf(exponent_limit)) continue;

            bool other_growth_is_slower = true;
            for (const auto& other : node.operands()) {
                if (other.get() == op.get()) continue;
                GrowthClass growth = classify_growth(other);
                if (growth == GrowthClass::Exponential ||
                    growth == GrowthClass::Unknown) {
                    other_growth_is_slower = false;
                    break;
                }
            }
            if (other_growth_is_slower) {
                result = lamina::detail::make_node<NumberNode>(BigInt(0));
                return;
            }
        }
    }
    {
        for (const auto& op : node.operands()) {
            auto function = std::dynamic_pointer_cast<const FunctionNode>(op);
            if (!function || function->arguments().size() != 1) continue;
            if (function->type() != FunctionNode::FuncType::Sin &&
                function->type() != FunctionNode::FuncType::Tan) {
                continue;
            }

            const auto& argument = function->arguments()[0];
            std::shared_ptr<const SymbolicNode> reciprocal_argument;
            for (const auto& candidate : node.operands()) {
                if (auto power = std::dynamic_pointer_cast<const PowerNode>(candidate)) {
                    auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
                    if (exponent && power->base() && power->base()->compare(*argument) == 0) {
                        double exponent_value = 0.0;
                        if (std::holds_alternative<double>(exponent->value())) {
                            exponent_value = std::get<double>(exponent->value());
                        } else if (std::holds_alternative<BigInt>(exponent->value())) {
                            exponent_value = std::get<BigInt>(exponent->value()).to_double();
                        } else if (std::holds_alternative<Rational>(exponent->value())) {
                            exponent_value = std::get<Rational>(exponent->value()).to_double();
                        }
                        if (exponent_value == -1.0) {
                            reciprocal_argument = candidate;
                            break;
                        }
                    }
                }
            }
            if (!reciprocal_argument) continue;

            auto argument_limit = eval_limit(argument);
            if (argument_limit && argument_limit->is_zero()) {
                std::vector<std::shared_ptr<const SymbolicNode>> remaining_limits;
                for (const auto& remaining : node.operands()) {
                    if (remaining.get() == op.get() ||
                        remaining.get() == reciprocal_argument.get()) {
                        continue;
                    }
                    auto remaining_limit = eval_limit(remaining);
                    if (!remaining_limit) {
                        remaining_limits.clear();
                        break;
                    }
                    remaining_limits.push_back(remaining_limit);
                }
                remaining_limits.push_back(lamina::detail::make_node<NumberNode>(BigInt(1)));
                auto reduced = remaining_limits.size() == 1
                    ? remaining_limits[0]
                    : make_product_or_one(remaining_limits);
                NormalizationVisitor norm;
                reduced->accept(norm);
                result = norm.get_result();
                return;
            }
        }
    }
    auto squeeze_result = try_squeeze(lamina::detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>(node.operands().begin(), node.operands().end())));
    if (squeeze_result) { result = squeeze_result; return; }
    std::vector<std::shared_ptr<const SymbolicNode>> subs_ops;
    for (auto& op : node.operands()) { op->accept(*this); subs_ops.push_back(result); }
    NormalizationVisitor norm;
    auto subs_res = lamina::detail::make_node<MultiplyNode>(subs_ops); subs_res->accept(norm); auto final_subs = norm.get_result();
    std::vector<std::shared_ptr<const SymbolicNode>> num_nodes, den_nodes;
    for (auto& op : node.operands()) {
        if (auto pow = std::dynamic_pointer_cast<const PowerNode>(op)) {
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
                double e = 0;
                if (std::holds_alternative<double>(num->value())) e = std::get<double>(num->value());
                else if (std::holds_alternative<BigInt>(num->value())) e = std::get<BigInt>(num->value()).to_double();
                else if (std::holds_alternative<Rational>(num->value())) e = std::get<Rational>(num->value()).to_double();
                if (e < 0) { den_nodes.push_back(lamina::detail::make_node<PowerNode>(pow->base(), lamina::detail::make_node<NumberNode>(-e))); continue; }
            }
        }
        num_nodes.push_back(op);
    }
    if (den_nodes.empty()) {
        auto prod_form = classify_product_form(subs_ops);
        if (prod_form == IndeterminateForm::ZeroTimesInf) { auto resolved = resolve_zero_times_inf(node.operands(), subs_ops); if (resolved) { result = resolved; return; } }
        result = final_subs; return;
    }
    auto N = make_product_or_one(num_nodes);
    auto D = make_product_or_one(den_nodes);
    auto val_n = eval_limit(N); auto val_d = eval_limit(D);
    bool n_zero = val_n && val_n->is_zero(), d_zero = val_d && val_d->is_zero();
    bool n_inf = is_inf(val_n), d_inf = is_inf(val_d);
    /// Limits at infinity: try rational function degree comparison and growth-rate comparison
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
        if (auto num = std::dynamic_pointer_cast<const NumberNode>(val_n)) {
            if (std::holds_alternative<double>(num->value())) sign_n = std::get<double>(num->value()) > 0 ? 1 : -1;
            else if (std::holds_alternative<BigInt>(num->value())) sign_n = std::get<BigInt>(num->value()) > BigInt(0) ? 1 : -1;
            else if (std::holds_alternative<Rational>(num->value())) sign_n = std::get<Rational>(num->value()) > Rational(0) ? 1 : -1;
        }
        int sign_d = (direction == "-" ? -1 : 1);
        int final_sign = sign_n * sign_d;
        std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
        auto inf_node = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
        if (final_sign < 0) { std::vector<std::shared_ptr<const SymbolicNode>> m_args = {lamina::detail::make_node<NumberNode>(BigInt(-1)), inf_node}; result = lamina::detail::make_node<MultiplyNode>(m_args); }
        else { result = inf_node; }
        return;
    }
    result = final_subs;
}

void LimitVisitor::visit(const PowerNode& node) {
    /// Handle negative infinity: substitute x = -t and evaluate lim(t→+∞)
    if (is_limit_at_neg_infinity()) {
        auto neg_inf_result = handle_neg_infinity_limit(lamina::detail::make_node<PowerNode>(node.base(), node.exponent()));
        if (neg_inf_result) { result = neg_inf_result; return; }
    }
    if (is_limit_at_infinity()) {
        auto base_add = std::dynamic_pointer_cast<const AddNode>(node.base());
        auto exponent_mul = std::dynamic_pointer_cast<const MultiplyNode>(node.exponent());
        auto exponent_var = std::dynamic_pointer_cast<const VariableNode>(node.exponent());
        if (base_add && (exponent_mul || (exponent_var && exponent_var->name() == var))) {
            std::shared_ptr<const SymbolicNode> small_term;
            bool has_one = false;
            for (const auto& op : base_add->operands()) {
                if (op->is_one()) has_one = true;
                else small_term = op;
            }
            bool exponent_has_var = false;
            std::shared_ptr<const SymbolicNode> exponent_coeff = lamina::detail::make_node<NumberNode>(BigInt(1));
            std::vector<std::shared_ptr<const SymbolicNode>> coeff_terms;
            if (exponent_var && exponent_var->name() == var) {
                exponent_has_var = true;
            } else {
                for (const auto& op : exponent_mul->operands()) {
                    auto variable = std::dynamic_pointer_cast<const VariableNode>(op);
                    if (variable && variable->name() == var) {
                        exponent_has_var = true;
                    } else {
                        coeff_terms.push_back(op);
                    }
                }
            }
            if (!coeff_terms.empty()) {
                exponent_coeff = coeff_terms.size() == 1
                    ? coeff_terms[0]
                    : make_product_or_one(coeff_terms);
            }
            if (has_one && small_term && exponent_has_var) {
                std::shared_ptr<const SymbolicNode> small_num = small_term;
                std::shared_ptr<const SymbolicNode> small_den = lamina::detail::make_node<NumberNode>(BigInt(1));
                if (auto small_mul = std::dynamic_pointer_cast<const MultiplyNode>(small_term)) {
                    std::vector<std::shared_ptr<const SymbolicNode>> num_parts;
                    for (const auto& op : small_mul->operands()) {
                        if (auto power = std::dynamic_pointer_cast<const PowerNode>(op)) {
                            auto base_var = std::dynamic_pointer_cast<const VariableNode>(power->base());
                            auto exp_num = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
                            if (base_var && base_var->name() == var && exp_num) {
                                double exp_value = 0.0;
                                if (std::holds_alternative<double>(exp_num->value())) exp_value = std::get<double>(exp_num->value());
                                else if (std::holds_alternative<BigInt>(exp_num->value())) exp_value = std::get<BigInt>(exp_num->value()).to_double();
                                else if (std::holds_alternative<Rational>(exp_num->value())) exp_value = std::get<Rational>(exp_num->value()).to_double();
                                if (exp_value == -1.0) {
                                    small_den = power->base();
                                    continue;
                                }
                            }
                        }
                        num_parts.push_back(op);
                    }
                    if (small_den && small_den->compare(*lamina::detail::make_node<VariableNode>(var)) == 0) {
                        small_num = num_parts.empty()
                            ? std::static_pointer_cast<const SymbolicNode>(lamina::detail::make_node<NumberNode>(BigInt(1)))
                            : make_product_or_one(num_parts);
                    }
                } else if (auto small_power = std::dynamic_pointer_cast<const PowerNode>(small_term)) {
                    auto base_var = std::dynamic_pointer_cast<const VariableNode>(small_power->base());
                    auto exp_num = std::dynamic_pointer_cast<const NumberNode>(small_power->exponent());
                    if (base_var && base_var->name() == var && exp_num) {
                        double exp_value = 0.0;
                        if (std::holds_alternative<double>(exp_num->value())) exp_value = std::get<double>(exp_num->value());
                        else if (std::holds_alternative<BigInt>(exp_num->value())) exp_value = std::get<BigInt>(exp_num->value()).to_double();
                        else if (std::holds_alternative<Rational>(exp_num->value())) exp_value = std::get<Rational>(exp_num->value()).to_double();
                        if (exp_value == -1.0) {
                            small_num = lamina::detail::make_node<NumberNode>(BigInt(1));
                            small_den = small_power->base();
                        }
                    }
                }

                auto den_var = std::dynamic_pointer_cast<const VariableNode>(small_den);
                if (den_var && den_var->name() == var) {
                    auto inner = lamina::detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>{
                        small_num, exponent_coeff});
                    NormalizationVisitor norm;
                    inner->accept(norm);
                    std::vector<std::shared_ptr<const SymbolicNode>> exp_args = {norm.get_result()};
                    auto exp_result = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Exp, exp_args);
                    exp_result->accept(norm);
                    result = norm.get_result();
                    return;
                }
            }
        }
    }
    node.base()->accept(*this); auto b = result;
    node.exponent()->accept(*this); auto e = result;
    auto form = classify_power_form(b, e);
    if (form != IndeterminateForm::None) { auto resolved = resolve_exponential_form(node.base(), node.exponent()); if (resolved) { result = resolved; return; } }
    /// Handle ∞^(negative) → 0 and ∞^(positive) → ∞
    if (is_inf(b)) {
        if (auto e_num = std::dynamic_pointer_cast<const NumberNode>(e)) {
            double ev = 0;
            if (std::holds_alternative<double>(e_num->value())) ev = std::get<double>(e_num->value());
            else if (std::holds_alternative<BigInt>(e_num->value())) ev = std::get<BigInt>(e_num->value()).to_double();
            else if (std::holds_alternative<Rational>(e_num->value())) ev = std::get<Rational>(e_num->value()).to_double();
            if (ev < 0) { result = lamina::detail::make_node<NumberNode>(BigInt(0)); return; }
            if (ev > 0) {
                std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
                result = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
                return;
            }
        }
    }
    /// Handle 0^(negative) → ±∞ with direction-aware sign determination
    if (b && b->is_zero()) {
        if (auto e_num = std::dynamic_pointer_cast<const NumberNode>(e)) {
            double ev = 0;
            if (std::holds_alternative<double>(e_num->value())) ev = std::get<double>(e_num->value());
            else if (std::holds_alternative<BigInt>(e_num->value())) ev = std::get<BigInt>(e_num->value()).to_double();
            else if (std::holds_alternative<Rational>(e_num->value())) ev = std::get<Rational>(e_num->value()).to_double();
            if (ev < 0 && !direction.empty()) {
                /// base→0, exponent is negative, direction specified → ±∞
                /// Determine sign based on approach direction and exponent parity
                std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
                auto inf_node = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
                int exp_int = static_cast<int>(ev);
                bool odd_exponent = (exp_int % 2 != 0);
                if (odd_exponent) {
                    /// For odd negative exponents (e.g., x^(-1), x^(-3)):
                    /// x→0+ gives +∞, x→0- gives -∞
                    int base_sign = determine_sign_near_point(node.base(), direction);
                    if (base_sign < 0) {
                        std::vector<std::shared_ptr<const SymbolicNode>> m = {lamina::detail::make_node<NumberNode>(BigInt(-1)), inf_node};
                        result = lamina::detail::make_node<MultiplyNode>(m);
                        return;
                    }
                }
                /// Even negative exponents always give +∞
                result = inf_node;
                return;
            }
        }
    }
    NormalizationVisitor norm; lamina::detail::make_node<PowerNode>(b, e)->accept(norm); result = norm.get_result();
}

void LimitVisitor::visit(const FunctionNode& node) {
    /// Handle negative infinity: substitute x = -t and evaluate lim(t→+∞)
    if (is_limit_at_neg_infinity() && node.type() != FunctionNode::FuncType::Infinity) {
        auto neg_inf_result = handle_neg_infinity_limit(lamina::detail::make_node<FunctionNode>(node.type(), node.arguments()));
        if (neg_inf_result) { result = neg_inf_result; return; }
    }
    /// 方向感知的符号函数处理
    if (node.type() == FunctionNode::FuncType::Sgn && node.arguments().size() == 1) {
        auto r = evaluate_sgn_limit(node.arguments()[0]);
        if (r) { result = *r; return; }
    }
    /// 方向感知的绝对值处理
    if (node.type() == FunctionNode::FuncType::Abs && node.arguments().size() == 1) {
        auto r = evaluate_abs_limit(node.arguments()[0]);
        if (r) { result = *r; return; }
    }
    /// Composed functions at infinity: evaluate inner limit first, then apply outer function
    if (is_limit_at_infinity() && node.arguments().size() == 1 && node.type() != FunctionNode::FuncType::Infinity) {
        auto inner_lim = eval_limit(node.arguments()[0]);
        if (inner_lim) {
            /// If inner limit is a finite number, evaluate the function at that value
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(inner_lim)) {
                std::vector<std::shared_ptr<const SymbolicNode>> new_args = {inner_lim};
                NormalizationVisitor norm;
                lamina::detail::make_node<FunctionNode>(node.type(), new_args)->accept(norm);
                result = norm.get_result();
                return;
            }
            /// If inner limit is infinity, handle based on function type
            if (is_inf(inner_lim)) {
                bool neg = is_neg_inf(inner_lim);
                switch (node.type()) {
                    case FunctionNode::FuncType::Exp:
                        if (neg) { result = lamina::detail::make_node<NumberNode>(BigInt(0)); }
                        else { std::vector<std::shared_ptr<const SymbolicNode>> inf_args; result = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args); }
                        return;
                    case FunctionNode::FuncType::Ln:
                        if (!neg) { std::vector<std::shared_ptr<const SymbolicNode>> inf_args; result = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args); }
                        return;
                    case FunctionNode::FuncType::ArcTan:
                        if (neg) { result = lamina::detail::make_node<NumberNode>(Rational(-1, 2)); result = lamina::detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>{result, lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, std::vector<std::shared_ptr<const SymbolicNode>>{})}); }
                        else { result = lamina::detail::make_node<NumberNode>(Rational(1, 2)); result = lamina::detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>{result, lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, std::vector<std::shared_ptr<const SymbolicNode>>{})}); }
                        return;
                    default:
                        break;
                }
            }
        }
    }
    std::vector<std::shared_ptr<const SymbolicNode>> new_args;
    for (auto& a : node.arguments()) { a->accept(*this); new_args.push_back(result); }
    NormalizationVisitor norm; lamina::detail::make_node<FunctionNode>(node.type(), new_args)->accept(norm); result = norm.get_result();
}
