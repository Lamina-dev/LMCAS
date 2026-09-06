#include "symbolic_set.hpp"

#include <algorithm>
#include <exception>
#include <utility>

#include "symbolic_ast.hpp"

namespace LMCAS {
namespace {

constexpr const char* kOperation = "symbolic_set";

SymbolicSetResult failure(CasErrc code, std::string message) {
    return SymbolicSetResult::failure(code, std::move(message), kOperation);
}

Result<void> consume(ComputationContext& context) {
    return context.consume_steps(1, kOperation);
}

std::shared_ptr<const FiniteSetNode> finite_set_node(const SymbolicExpr& expression) {
    return std::dynamic_pointer_cast<const FiniteSetNode>(detail::node(expression));
}

bool exact_value(const std::shared_ptr<const SymbolicNode>& node, Rational& value) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return false;
    if (std::holds_alternative<BigInt>(number->value())) {
        value = Rational(std::get<BigInt>(number->value()));
        return true;
    }
    if (std::holds_alternative<Rational>(number->value())) {
        value = std::get<Rational>(number->value());
        return true;
    }
    return false;
}

SymbolicSetResult from_nodes(std::vector<std::shared_ptr<const SymbolicNode>> elements) {
    try {
        return SymbolicSetResult::success(detail::make_expression_ptr(
            detail::make_node<FiniteSetNode>(std::move(elements))));
    } catch (const std::bad_alloc&) {
        return failure(CasErrc::ResourceLimit, "finite set allocation failed");
    } catch (const std::exception& error) {
        return failure(CasErrc::InvalidArgument, error.what());
    }
}

} // namespace

SymbolicSetResult make_finite_set(
    std::vector<std::shared_ptr<SymbolicExpr>> elements,
    ComputationContext& context) {
    auto step = consume(context);
    if (!step) return SymbolicSetResult::failure(step.error());
    std::vector<std::shared_ptr<const SymbolicNode>> nodes;
    nodes.reserve(elements.size());
    for (const auto& element : elements) {
        if (!element || !detail::node(element)) {
            return failure(CasErrc::InvalidArgument,
                           "finite set elements cannot be null");
        }
        nodes.push_back(detail::node(element));
    }
    return from_nodes(std::move(nodes));
}

SymbolicSetResult make_interval(
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    bool lower_closed,
    bool upper_closed,
    ComputationContext& context) {
    auto step = consume(context);
    if (!step) return SymbolicSetResult::failure(step.error());
    if (!lower || !upper || !detail::node(lower) || !detail::node(upper)) {
        return failure(CasErrc::InvalidArgument, "interval endpoints cannot be null");
    }
    Rational lower_value;
    Rational upper_value;
    if (exact_value(detail::node(lower), lower_value) &&
        exact_value(detail::node(upper), upper_value) &&
        upper_value < lower_value) {
        return failure(CasErrc::DomainError,
                       "interval lower endpoint exceeds upper endpoint");
    }
    try {
        return SymbolicSetResult::success(detail::make_expression_ptr(
            detail::make_node<IntervalNode>(detail::node(lower), detail::node(upper),
                                            lower_closed, upper_closed)));
    } catch (const std::bad_alloc&) {
        return failure(CasErrc::ResourceLimit, "interval allocation failed");
    }
}

SymbolicSetResult make_membership(
    const std::shared_ptr<SymbolicExpr>& element,
    const std::shared_ptr<SymbolicExpr>& set,
    ComputationContext& context) {
    auto step = consume(context);
    if (!step) return SymbolicSetResult::failure(step.error());
    if (!element || !set || !detail::node(element) || !detail::node(set)) {
        return failure(CasErrc::InvalidArgument, "membership operands cannot be null");
    }
    if (auto finite = finite_set_node(*set)) {
        return SymbolicSetResult::success(SymbolicExpr::number(
            finite->contains(*detail::node(element)) ? 1 : 0));
    }
    if (auto interval = std::dynamic_pointer_cast<const IntervalNode>(detail::node(set))) {
        Rational value;
        Rational lower;
        Rational upper;
        if (exact_value(detail::node(element), value) &&
            exact_value(interval->lower(), lower) && exact_value(interval->upper(), upper)) {
            const bool above = lower < value || (interval->lower_closed() && lower == value);
            const bool below = value < upper || (interval->upper_closed() && value == upper);
            return SymbolicSetResult::success(SymbolicExpr::number(above && below ? 1 : 0));
        }
    }
    return SymbolicSetResult::success(detail::make_expression_ptr(
        detail::make_node<MembershipNode>(detail::node(element), detail::node(set))));
}

SymbolicSetResult finite_set_union(
    const SymbolicExpr& lhs, const SymbolicExpr& rhs,
    ComputationContext& context) {
    auto step = consume(context);
    if (!step) return SymbolicSetResult::failure(step.error());
    auto left = finite_set_node(lhs);
    auto right = finite_set_node(rhs);
    if (!left || !right) return failure(CasErrc::SetOperandTypeMismatch,
                                        "set union requires two finite sets");
    auto elements = left->elements();
    elements.insert(elements.end(), right->elements().begin(), right->elements().end());
    return from_nodes(std::move(elements));
}

SymbolicSetResult finite_set_intersection(
    const SymbolicExpr& lhs, const SymbolicExpr& rhs,
    ComputationContext& context) {
    auto step = consume(context);
    if (!step) return SymbolicSetResult::failure(step.error());
    auto left = finite_set_node(lhs);
    auto right = finite_set_node(rhs);
    if (!left || !right) return failure(CasErrc::SetOperandTypeMismatch,
                                        "set intersection requires two finite sets");
    std::vector<std::shared_ptr<const SymbolicNode>> elements;
    for (const auto& element : left->elements()) {
        if (right->contains(*element)) elements.push_back(element);
    }
    return from_nodes(std::move(elements));
}

SymbolicSetResult finite_set_difference(
    const SymbolicExpr& lhs, const SymbolicExpr& rhs,
    ComputationContext& context) {
    auto step = consume(context);
    if (!step) return SymbolicSetResult::failure(step.error());
    auto left = finite_set_node(lhs);
    auto right = finite_set_node(rhs);
    if (!left || !right) return failure(CasErrc::SetOperandTypeMismatch,
                                        "set difference requires two finite sets");
    std::vector<std::shared_ptr<const SymbolicNode>> elements;
    for (const auto& element : left->elements()) {
        if (!right->contains(*element)) elements.push_back(element);
    }
    return from_nodes(std::move(elements));
}

SymbolicSetResult finite_set_symmetric_difference(
    const SymbolicExpr& lhs, const SymbolicExpr& rhs,
    ComputationContext& context) {
    auto left = finite_set_difference(lhs, rhs, context);
    if (!left) return left;
    auto right = finite_set_difference(rhs, lhs, context);
    if (!right) return right;
    return finite_set_union(*left.value(), *right.value(), context);
}

} // namespace LMCAS
