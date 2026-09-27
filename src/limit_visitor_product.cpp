#include "internal/visitors/limit_visitor.hpp"

namespace LMCAS {

std::optional<std::shared_ptr<const SymbolicNode>> LimitVisitor::exact_ratio_at_point(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den) {
    auto nl = exact_leading_at_point(num);
    auto dl = exact_leading_at_point(den);
    if (nl && dl && dl->first != Rational(0)) {
        if (nl->first == Rational(0) || nl->second > dl->second)
            { return detail::make_node<NumberNode>(BigInt(0)); }
        Rational coefficient = nl->first / dl->first;
        if (nl->second == dl->second) { return detail::make_node<NumberNode>(coefficient); }
        if (direction.empty() && (dl->second - nl->second) % 2) { return std::shared_ptr<const SymbolicNode>{}; }
        int sign = coefficient > Rational(0) ? 1 : -1;
        if (direction == "-" && (dl->second - nl->second) % 2) { sign = -sign; }
        return signed_infinity(sign);
    }
    return std::nullopt;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::evaluated_product_ratio(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den, bool at_infinity) {
    auto numerator = eval_limit(num);
    auto denominator = eval_limit(den);
    bool numerator_zero = is_zero_limit(numerator);
    bool denominator_zero = is_zero_limit(denominator);
    const bool both_infinite = is_inf(numerator) && is_inf(denominator);
    if (both_infinite && at_infinity) {
        auto growth = limit_by_growth_comparison(num, den);
        if (growth) {
            return growth;
        }
    }
    if ((numerator_zero && denominator_zero) || both_infinite) {
        auto value = apply_lhopital(num, den);
        if (value) {
            return value;
        }
    }
    if (!numerator_zero && denominator_zero) {
        return singular_ratio_limit(num, den);
    }
    return nullptr;
}

bool LimitVisitor::positive_tail_product_limit(const MultiplyNode& node) {
    std::vector<std::shared_ptr<const SymbolicNode>> rewritten;
    bool changed = false;
    for (const auto& operand : node.operands()) {
        auto function = std::dynamic_pointer_cast<const FunctionNode>(operand);
        auto replacement = function ? positive_tail_log_power(*function) : nullptr;
        changed = changed || replacement != nullptr;
        rewritten.push_back(replacement ? replacement : operand);
    }
    if (changed) {
        result = eval_limit(detail::make_node<MultiplyNode>(std::move(rewritten)));
        return true;
    }
    return false;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::resolve_product_factors(
    const MultiplyNode& node, const std::vector<std::shared_ptr<const SymbolicNode>>& values) {
    std::vector<std::shared_ptr<const SymbolicNode>> numerator, denominator;
    for (const auto& operand : node.operands()) {
        auto factor = denominator_factor(operand);
        if (factor) { denominator.push_back(factor); }
        else { numerator.push_back(operand); }
    }
    std::shared_ptr<const SymbolicNode> resolved;
    if (denominator.empty()) {
        if (classify_product_form(values) == IndeterminateForm::ZeroTimesInf) {
            resolved = resolve_zero_times_inf(node.operands(), values);
        }
    } else {
        auto num = make_product_or_one(numerator);
        auto den = make_product_or_one(denominator);
        resolved = product_ratio_limit(num, den);
    }
    return resolved;
}



std::shared_ptr<const SymbolicNode> LimitVisitor::product_ratio_limit(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den) {
    const bool at_infinity = is_limit_at_infinity() || is_limit_at_neg_infinity();
    if (!at_infinity) {
        auto exact = exact_ratio_at_point(num, den);
        if (exact) { return std::move(*exact); }
    }
    if (at_infinity) {
        auto rational = limit_rational_at_infinity(num, den);
        if (rational) {
            return rational;
        }
    }
    return evaluated_product_ratio(num, den, at_infinity);
}

void LimitVisitor::visit(const MultiplyNode& node) {
    if (is_limit_at_neg_infinity()) {
        auto value = handle_neg_infinity_limit(LMCAS::detail::make_node<MultiplyNode>(node.operands()));
        if (value) {
            result = value;
            return;
        }
    }
    if (is_limit_at_infinity()) {
        if (positive_tail_product_limit(node)) { return; }
    }
    if (logarithmic_zero_product(node)) {
        result = LMCAS::detail::make_node<NumberNode>(BigInt(0));
        return;
    }
    if (is_limit_at_infinity()) {
        if (logarithmic_product_limit(node)) {
            return;
        }
        if (decaying_exponential_product(node)) {
            result = LMCAS::detail::make_node<NumberNode>(BigInt(0));
            return;
        }
    }
    if (trigonometric_product_limit(node)) {
        return;
    }
    auto squeezed = try_squeeze(LMCAS::detail::make_node<MultiplyNode>(node.operands()));
    if (squeezed) {
        result = squeezed;
        return;
    }
    std::vector<std::shared_ptr<const SymbolicNode>> values;
    for (const auto& operand : node.operands()) {
        operand->accept(*this);
        values.push_back(result);
    }
    std::shared_ptr<const SymbolicNode> substituted;
    if (std::all_of(values.begin(), values.end(), [](const auto& value) { return value != nullptr; })) {
        NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
        LMCAS::detail::make_node<MultiplyNode>(values)->accept(norm);
        substituted = norm.get_result();
    }
    auto resolved = resolve_product_factors(node, values);
    result = resolved ? resolved : substituted;
}

}
