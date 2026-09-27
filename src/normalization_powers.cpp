#include "internal/visitors/normalization_visitor.hpp"
#include "internal/normalization_utils.hpp"

namespace LMCAS {

bool NormalizationVisitor::normalize_power_identity(
    const std::shared_ptr<const SymbolicNode>& base,
    const std::shared_ptr<const SymbolicNode>& exponent) {
    if (exponent->is_zero()) {
        if (normalization_can_discard(base, facts_, domain_, context_) &&
            (domain_ == Domain::Complex ||
             normalization_proved(detail::query_nonzero_value(base, facts_, Domain::Complex, context_)))) {
            set_result(detail::make_node<NumberNode>(BigInt(1)));
        } else {
            normalization_check_children(rewrite_budget(), 1, base, exponent);
            set_result(detail::make_node<PowerNode>(base, exponent));
        }
        return true;
    }
    if (exponent->is_one()) {
        set_result(base);
        return true;
    }
    if (base->is_zero()) {
        if (normalization_can_discard(exponent, facts_, domain_, context_) &&
            normalization_proved(detail::query_positive_value(exponent, facts_, context_))) {
            set_result(detail::make_node<NumberNode>(BigInt(0)));
        } else {
            normalization_check_children(rewrite_budget(), 1, base, exponent);
            set_result(detail::make_node<PowerNode>(base, exponent));
        }
        return true;
    }
    if (!base->is_one() || !normalization_can_discard(exponent, facts_, domain_, context_)) {
        return false;
    }
    set_result(detail::make_node<NumberNode>(BigInt(1)));
    return true;
}

bool NormalizationVisitor::normalize_product_power(
    const MultiplyNode& base, const std::shared_ptr<const SymbolicNode>& exponent) {
    for (const auto& operand : base.operands()) {
        const auto* power = dynamic_cast<const PowerNode*>(operand.get());
        if (dynamic_cast<const MatrixNode*>(operand.get()) ||
            (power && dynamic_cast<const MatrixNode*>(power->base().get()))) {
            return false;
        }
    }
    BigInt integer;
    if (!try_get_integer_value(std::dynamic_pointer_cast<const NumberNode>(exponent), integer)) {
        if (!normalization_proved(detail::query_real_value(exponent, facts_, context_))) {
            return false;
        }
        const bool positive = normalization_proved(detail::query_positive_value(exponent, facts_, context_));
        for (const auto& operand : base.operands()) {
            if (!normalization_proved(positive
                    ? detail::query_nonnegative_value(operand, facts_, context_)
                    : detail::query_positive_value(operand, facts_, context_))) {
                return false;
            }
        }
    } else if (integer.is_zero()) {
        return false;
    }
    distribute_power(base, exponent);
    return true;
}

void NormalizationVisitor::distribute_power(
    const MultiplyNode& base, const std::shared_ptr<const SymbolicNode>& exponent) {
    std::vector<std::shared_ptr<const SymbolicNode>> operands;
    std::size_t nodes = 0;
    normalization_check_count(rewrite_budget(), base.operands().size(), 1, false);
    for (const auto& operand : base.operands()) {
        normalization_check_children(rewrite_budget(), 1, operand, exponent);
        detail::make_node<PowerNode>(operand, exponent)->accept(*this);
        if (const auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(result)) {
            for (const auto& child : multiply->operands()) {
                normalization_append(rewrite_budget(), nodes, operands, child);
            }
        } else {
            normalization_append(rewrite_budget(), nodes, operands, result);
        }
    }
    make_normalized_multiply_node(operands, rewrite_budget())->accept(*this);
}

bool NormalizationVisitor::normalize_nested_imaginary_power(
    const PowerNode& base, const std::shared_ptr<const NumberNode>& exponent) {
    BigInt integer;
    if (!try_get_integer_value(exponent, integer)) {
        return false;
    }
    const auto inner_base = std::dynamic_pointer_cast<const NumberNode>(base.base());
    if (!inner_base || !inner_base->is_negative_one()) {
        return false;
    }
    const auto inner_exponent = std::dynamic_pointer_cast<const NumberNode>(base.exponent());
    if (!inner_exponent) {
        return false;
    }
    const auto* rational = std::get_if<Rational>(&inner_exponent->value());
    if (!rational || *rational != Rational(1, 2)) {
        return false;
    }
    set_result(normalization_imaginary_power(integer, rewrite_budget()));
    return true;
}

bool NormalizationVisitor::normalize_nested_power(
    const PowerNode& base, const std::shared_ptr<const SymbolicNode>& exponent) {
    const auto numeric = std::dynamic_pointer_cast<const NumberNode>(exponent);
    if (normalize_nested_imaginary_power(base, numeric)) {
        return true;
    }
    if (!can_combine_nested_power(base, numeric)) {
        return false;
    }
    make_normalized_multiply_node({base.exponent(), exponent}, rewrite_budget())->accept(*this);
    normalization_check_children(rewrite_budget(), 1, base.base(), result);
    detail::make_node<PowerNode>(base.base(), result)->accept(*this);
    return true;
}

void NormalizationVisitor::visit(const PowerNode& node) {
    std::shared_ptr<const SymbolicNode> base;
    if (try_normalize_squared_norm(node, base)) {
        return;
    }
    if (!base) {
        node.base()->accept(*this);
        base = result;
    }
    node.exponent()->accept(*this);
    auto exponent = result;
    if (normalize_power_identity(base, exponent)) {
        return;
    }
    const auto base_number = std::dynamic_pointer_cast<const NumberNode>(base);
    const auto exponent_number = std::dynamic_pointer_cast<const NumberNode>(exponent);
    if (base_number && exponent_number &&
        normalization_can_discard(base, facts_, domain_, context_) && normalization_can_discard(exponent, facts_, domain_, context_)) {
        if (auto folded = normalization_numeric_power(base_number, exponent_number, rewrite_budget())) {
            set_result(std::move(folded));
            return;
        }
    }
    if (const auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(base)) {
        if (normalize_product_power(*multiply, exponent)) {
            return;
        }
    }
    if (const auto power = std::dynamic_pointer_cast<const PowerNode>(base)) {
        if (normalize_nested_power(*power, exponent)) {
            return;
        }
    }
    normalization_check_children(rewrite_budget(), 1, base, exponent);
    set_result(detail::make_node<PowerNode>(base, exponent));
}

bool NormalizationVisitor::is_positive_integer_number(const std::shared_ptr<const NumberNode>& node) {
        BigInt value;
        return try_get_integer_value(node, value) && value > BigInt(0);
    }
bool NormalizationVisitor::can_combine_nested_power(
    const PowerNode& inner,
    const std::shared_ptr<const NumberNode>& outer_exponent) const {
    if (!is_positive_integer_number(outer_exponent)) {
        return false;
    }
    auto exponent = std::dynamic_pointer_cast<const NumberNode>(inner.exponent());
    if (is_positive_integer_number(exponent)) {
        return true;
    }
    if (!exponent) {
        return false;
    }
    const auto* rational = std::get_if<Rational>(&exponent->value());
    if (!rational || rational->get_numerator() <= BigInt(0)) {
        return false;
    }
    return normalization_proved(detail::query_nonnegative_value(inner.base(), facts_, context_));
}

}
