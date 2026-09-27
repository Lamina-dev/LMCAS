#include "internal/visitors/limit_visitor.hpp"

namespace LMCAS {

std::shared_ptr<const SymbolicNode> LimitVisitor::singular_ratio_limit(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den) {
    auto ns = determine_sign_near_point(num, direction);
    auto ds = determine_sign_near_point(den, direction);
    if (ns && ds && *ns != 0 && *ds != 0) { return signed_infinity(*ns * *ds); }
    return nullptr;
}



bool LimitVisitor::tends_to_zero(const std::shared_ptr<const SymbolicNode>& node) const {
    EvaluationScope scope;
    if (!scope) {
        return false;
    }
    LimitVisitor sub_vis(var, point, direction, assumption_ctx_, context_);
    auto node_copy = node->clone();
    node_copy->accept(sub_vis);
    auto val = sub_vis.get_result();
    if (!val) {
        return false;
    }
    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget()); val->accept(norm); val = norm.get_result();
    return val && val->is_zero();
}

std::shared_ptr<const SymbolicNode> LimitVisitor::eval_limit(const std::shared_ptr<const SymbolicNode>& expr) {
    EvaluationScope scope;
    if (!scope) {
        return nullptr;
    }
    LimitVisitor sub(var, point, direction, assumption_ctx_, context_);
    sub.lhopital_depth_ = this->lhopital_depth_;
    expr->accept(sub);
    auto r = sub.get_result();
    if (r) { NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget()); r->accept(norm); return norm.get_result(); }
    return r;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::derivative_limit(
    const std::shared_ptr<const SymbolicNode>& node, bool normalize) {
    LimitVisitor visitor(var, point, direction, assumption_ctx_, context_);
    visitor.lhopital_depth_ = lhopital_depth_ + 1;
    node->accept(visitor);
    auto value = visitor.get_result();
    if (!value || !normalize) {
        return value;
    }
    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
    value->accept(norm);
    return norm.get_result();
}

std::shared_ptr<const SymbolicNode> LimitVisitor::derivative_ratio_limit(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den,
    const std::shared_ptr<const SymbolicNode>& num_value,
    const std::shared_ptr<const SymbolicNode>& den_value) {
    bool numerator_zero = num_value->is_zero();
    bool denominator_zero = den_value->is_zero();
    bool numerator_infinity = is_inf(num_value);
    bool denominator_infinity = is_inf(den_value);
    bool indeterminate = numerator_zero && denominator_zero;
    if (!indeterminate) {
        indeterminate = numerator_infinity && denominator_infinity;
    }
    if (!indeterminate) {
        if (denominator_zero && !numerator_zero) {
            return singular_ratio_limit(num, den);
        }
        if (!denominator_zero && !denominator_infinity) {
            auto ratio = LMCAS::detail::make_node<MultiplyNode>(
                std::vector<std::shared_ptr<const SymbolicNode>>{
                    num_value, LMCAS::detail::make_node<PowerNode>(
                        den_value, LMCAS::detail::make_node<NumberNode>(BigInt(-1)))});
            NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
            ratio->accept(norm);
            return norm.get_result();
        }
        if (!numerator_infinity && !numerator_zero && denominator_infinity) {
            return LMCAS::detail::make_node<NumberNode>(BigInt(0));
        }
    }
    auto ratio = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{
            num, LMCAS::detail::make_node<PowerNode>(
                den, LMCAS::detail::make_node<NumberNode>(BigInt(-1)))});
    return derivative_limit(ratio, false);
}

std::shared_ptr<const SymbolicNode> LimitVisitor::apply_lhopital(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den) {
    EvaluationScope scope;
    if (!scope) {
        return nullptr;
    }
    /**
     * @brief 用解析函数芽的首项判断可去相消，减少嵌套商的递归求导。
     * 同样适用于由异号无穷项合成的分式。
     */
    if (!is_inf(point)) {
        auto leading = taylor_fallback(num, den);
        if (leading) return leading;
    }
    if (lhopital_depth_ >= max_lhopital_depth) {
        return taylor_fallback(num, den);
    }
    DifferentiationVisitor differentiation(var);
    num->accept(differentiation);
    auto numerator = differentiation.get_result();
    den->accept(differentiation);
    auto denominator = differentiation.get_result();
    if (!numerator || !denominator) {
        return taylor_fallback(num, den);
    }
    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
    numerator->accept(norm);
    numerator = norm.get_result();
    denominator->accept(norm);
    denominator = norm.get_result();
    if (denominator->is_zero()) {
        return taylor_fallback(num, den);
    }
    auto ratio = LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{
            numerator, LMCAS::detail::make_node<PowerNode>(
                denominator, LMCAS::detail::make_node<NumberNode>(BigInt(-1)))});
    auto simplified = simplify_and_eval_ratio(ratio);
    if (simplified) {
        return simplified;
    }
    auto numerator_value = derivative_limit(numerator, true);
    auto denominator_value = derivative_limit(denominator, true);
    if (!numerator_value || !denominator_value) {
        return taylor_fallback(num, den);
    }
    return derivative_ratio_limit(numerator, denominator, numerator_value, denominator_value);
}

}
