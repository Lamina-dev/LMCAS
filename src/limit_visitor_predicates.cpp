#include "internal/visitors/limit_visitor.hpp"
#include "internal/normalization_utils.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/facts_query.hpp"
#include "symbolic.hpp"

namespace LMCAS {

bool LimitVisitor::is_infinite_product(const MultiplyNode& product) const {
    bool infinite = false;
    for (const auto& operand : product.operands()) {
        if (is_inf(operand)) { infinite = true; }
        else {
            auto sign = get_node_sign(operand);
            if (!sign || !*sign) { return false; }
        }
    }
    return infinite;
}

bool LimitVisitor::is_infinite_sum(const AddNode& sum) const {
    int sign = 0;
    for (const auto& operand : sum.operands()) {
        if (is_inf(operand)) {
            int next = is_neg_inf(operand) ? -1 : 1;
            if (sign && sign != next) { return false; }
            sign = next;
        } else if (!std::dynamic_pointer_cast<const NumberNode>(operand)) { return false; }
    }
    return sign != 0;
}

std::optional<int> LimitVisitor::assumed_node_sign(
    const std::shared_ptr<const SymbolicNode>& node) const {
    std::optional<detail::AssumptionFacts> assumed;
    if (assumption_ctx_) { assumed.emplace(*assumption_ctx_); }
    const auto& facts = assumed ? static_cast<const FactsQuery&>(*assumed) : detail::no_facts();
    auto nonzero = detail::query_nonzero_value(node, facts, Domain::Real, *context_);
    if (!nonzero) { throw nonzero.error(); }
    if (nonzero.value() == Tribool::False) { return 0; }
    auto positive = detail::query_positive_value(node, facts, *context_);
    if (!positive) { throw positive.error(); }
    if (positive.value() == Tribool::True) { return 1; }
    auto nonnegative = detail::query_nonnegative_value(node, facts, *context_);
    if (!nonnegative) { throw nonnegative.error(); }
    if (nonnegative.value() == Tribool::False) {
        auto real = detail::query_real_value(node, facts, *context_);
        if (!real) { throw real.error(); }
        if (real.value() == Tribool::True) { return -1; }
    }
    return std::nullopt;
}

bool LimitVisitor::is_infinite_power(const PowerNode& node) const {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node.exponent());
    BigInt integer;
    return is_inf(node.base()) && number && try_get_integer_value(number, integer) && integer > BigInt(0);
}

bool LimitVisitor::is_inf(const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) { return false; }
    if (auto f = std::dynamic_pointer_cast<const FunctionNode>(node))
        { return f->type() == FunctionNode::FuncType::Infinity; }
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return is_infinite_product(*product);
    }
    if (auto sum = std::dynamic_pointer_cast<const AddNode>(node)) {
        return is_infinite_sum(*sum);
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) { return is_infinite_power(*power); }
    return false;
}

bool LimitVisitor::is_neg_inf(const std::shared_ptr<const SymbolicNode>& node) const {
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        if (!is_inf(node)) { return false; }
        int sign = 1;
        for (const auto& operand : product->operands())
            { sign *= is_inf(operand) ? (is_neg_inf(operand) ? -1 : 1) : *get_node_sign(operand); }
        return sign < 0;
    }
    if (auto sum = std::dynamic_pointer_cast<const AddNode>(node)) {
        if (!is_inf(node)) { return false; }
        for (const auto& operand : sum->operands()) { if (is_inf(operand)) { return is_neg_inf(operand); } }
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        BigInt integer;
        auto number = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
        return number && try_get_integer_value(number, integer) && integer.is_odd() && is_neg_inf(power->base());
    }
    return false;
}

std::optional<int> LimitVisitor::get_node_sign(
    const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) {
        return std::nullopt;
    }
    if (node->is_zero()) {
        return 0;
    }
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        if (std::holds_alternative<double>(number->value())) {
            if (!std::isfinite(std::get<double>(number->value()))) { return std::nullopt; }
            return get_sign(node);
        }
        if (std::holds_alternative<BigInt>(number->value())) {
            return get_sign(node);
        }
        if (std::holds_alternative<Rational>(number->value())) {
            return get_sign(node);
        }
    }
    if (is_inf(node)) {
        return is_neg_inf(node) ? -1 : 1;
    }
    return assumed_node_sign(node);
}

double LimitVisitor::get_numeric_value(const std::shared_ptr<const NumberNode>& num) const {
    if (std::holds_alternative<double>(num->value())) {
        return std::get<double>(num->value());
    }
    if (std::holds_alternative<BigInt>(num->value())) {
        return std::get<BigInt>(num->value()).to_double();
    }
    if (std::holds_alternative<Rational>(num->value())) {
        return std::get<Rational>(num->value()).to_double();
    }
    return std::numeric_limits<double>::quiet_NaN();
}

bool LimitVisitor::is_bounded_product(const MultiplyNode& node) const {
    bool has_bounded = false;
    for (const auto& operand : node.operands()) {
        if (is_bounded(operand)) {
            has_bounded = true;
        } else if (!std::dynamic_pointer_cast<const NumberNode>(operand)) {
            return false;
        }
    }
    return has_bounded;
}

bool LimitVisitor::defined_near_point(const std::shared_ptr<const SymbolicNode>& node) const {
    std::optional<detail::AssumptionFacts> assumed;
    if (assumption_ctx_) { assumed.emplace(*assumption_ctx_); }
    const auto& facts = assumed ? static_cast<const FactsQuery&>(*assumed) : detail::no_facts();
    auto& context = *context_;
    auto obligations = detail::domain_constraints(node, facts, Domain::Real, context);
    if (!obligations) { throw obligations.error(); }
    if (!obligations.value()) { return false; }
    LimitVisitor witness(var, point, direction, assumption_ctx_, &context);
    for (const auto& condition : *obligations.value())
        { if (!witness.condition_holds_near_point(detail::node(condition))) { return false; } }
    return true;
}

bool LimitVisitor::is_bounded(const std::shared_ptr<const SymbolicNode>& node) const {
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return (function->type() == FunctionNode::FuncType::Sin ||
                function->type() == FunctionNode::FuncType::Cos ||
                function->type() == FunctionNode::FuncType::ArcTan) && defined_near_point(node);
    }
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return is_bounded_product(*product);
    }
    auto power = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!power || !is_bounded(power->base())) {
        return false;
    }
    auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
    return exponent && get_numeric_value(exponent) > 0;
}

bool LimitVisitor::is_bounded_expression(const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) {
        return false;
    }
    if (is_bounded(node)) {
        return true;
    }
    auto sum = std::dynamic_pointer_cast<const AddNode>(node);
    if (!sum) {
        return false;
    }
    bool has_bounded = false;
    for (const auto& operand : sum->operands()) {
        if (is_bounded(operand) || is_bounded_expression(operand)) {
            has_bounded = true;
        } else if (!std::dynamic_pointer_cast<const NumberNode>(operand)) {
            return false;
        }
    }
    return has_bounded;
}

bool LimitVisitor::is_zero_limit(const std::shared_ptr<const SymbolicNode>& node) const {
    return node && node->is_zero();
}

LimitVisitor::IndeterminateForm LimitVisitor::classify_product_form(
    const std::vector<std::shared_ptr<const SymbolicNode>>& values) {
    bool zero = false;
    bool infinity = false;
    for (const auto& value : values) {
        if (is_zero_limit(value)) zero = true;
        if (is_inf(value)) infinity = true;
    }
    return zero && infinity ? IndeterminateForm::ZeroTimesInf : IndeterminateForm::None;
}

LimitVisitor::IndeterminateForm LimitVisitor::classify_power_form(
    const std::shared_ptr<const SymbolicNode>& base,
    const std::shared_ptr<const SymbolicNode>& exponent) {
    if (base && base->is_one() && is_inf(exponent)) {
        return IndeterminateForm::OnePowInf;
    }
    if (!is_zero_limit(exponent)) {
        return IndeterminateForm::None;
    }
    if (is_zero_limit(base)) {
        return IndeterminateForm::ZeroPowZero;
    }
    if (is_inf(base)) {
        return IndeterminateForm::InfPowZero;
    }
    return IndeterminateForm::None;
}

LimitVisitor::IndeterminateForm LimitVisitor::classify_add_form(
    const std::vector<std::shared_ptr<const SymbolicNode>>& values) {
    bool positive = false;
    bool negative = false;
    for (const auto& value : values) {
        if (is_inf(value) && !is_neg_inf(value)) positive = true;
        if (is_neg_inf(value)) negative = true;
    }
    return positive && negative ? IndeterminateForm::InfMinusInf : IndeterminateForm::None;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::signed_infinity(int sign) const {
    auto infinity = LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Infinity, std::vector<std::shared_ptr<const SymbolicNode>>{});
    if (sign >= 0) {
        return infinity;
    }
    return LMCAS::detail::make_node<MultiplyNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{
            LMCAS::detail::make_node<NumberNode>(BigInt(-1)), infinity});
}

}
