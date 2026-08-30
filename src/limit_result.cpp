#include "limit_result.hpp"

#include "internal/expression_analysis.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "visitors/limit_visitor.hpp"

#include <exception>

namespace lamina {
namespace {

bool is_infinity(const std::shared_ptr<const SymbolicNode>& node) {
    auto function = std::dynamic_pointer_cast<const FunctionNode>(node);
    return function && function->type() == FunctionNode::FuncType::Infinity;
}

bool is_negative_infinity(const std::shared_ptr<const SymbolicNode>& node) {
    auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!multiply) return false;
    bool negative = false;
    bool infinity = false;
    for (const auto& operand : multiply->operands()) {
        infinity = infinity || is_infinity(operand);
        auto number = std::dynamic_pointer_cast<const NumberNode>(operand);
        if (number && std::holds_alternative<BigInt>(number->value()) &&
            std::get<BigInt>(number->value()) == BigInt(-1)) {
            negative = true;
        }
    }
    return negative && infinity;
}

bool depends_on(const LimitExprPtr& expression, const std::string& variable) {
    return expression_depends_on_variable(
        lamina::detail::node(expression), variable);
}

} // namespace

LimitResult limit_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction,
    ComputationContext& context) {
    constexpr const char* operation = "limit";
    if (!expression || !point || !lamina::detail::node(expression) ||
        !lamina::detail::node(point) || variable.empty()) {
        return LimitResult::failure(
            CasErrc::InvalidArgument,
            "limit requires an expression, named variable, and point",
            operation);
    }
    auto step = context.consume_steps(1, operation);
    if (!step) return LimitResult::failure(step.error());

    if (is_infinity(lamina::detail::node(point)) ||
        is_negative_infinity(lamina::detail::node(point))) {
        auto function = std::dynamic_pointer_cast<const FunctionNode>(
            lamina::detail::node(expression));
        if (function && function->arguments().size() == 1 &&
            (function->type() == FunctionNode::FuncType::Sin ||
             function->type() == FunctionNode::FuncType::Cos) &&
            depends_on(expression, variable)) {
            return LimitResult::success(Verified<LimitOutcome>{
                LimitDoesNotExist{}, ByConstructionProof{}});
        }
    }

    try {
        const char* direction_token = "";
        if (direction == LimitDirection::FromAbove) direction_token = "+";
        if (direction == LimitDirection::FromBelow) direction_token = "-";
        LimitVisitor visitor(
            variable, lamina::detail::node(point), direction_token,
            context.assumptions().get());
        lamina::detail::node(expression)->accept(visitor);
        auto result_node = visitor.get_result();
        if (!result_node ||
            lamina::detail::contains_node_type<LimitNode>(result_node) ||
            expression_depends_on_variable(result_node, variable)) {
            return LimitResult::failure(
                CasErrc::Inconclusive,
                "limit could not be proved in the supported domain",
                operation);
        }
        auto value = lamina::detail::make_expression_ptr(result_node);
        if (is_infinity(result_node)) {
            return LimitResult::success(Verified<LimitOutcome>{
                PositiveInfinityLimit{}, ByConstructionProof{}});
        }
        if (is_negative_infinity(result_node)) {
            return LimitResult::success(Verified<LimitOutcome>{
                NegativeInfinityLimit{}, ByConstructionProof{}});
        }
        return LimitResult::success(Verified<LimitOutcome>{
            FiniteLimit{std::move(value)}, ByConstructionProof{}});
    } catch (const std::bad_alloc&) {
        return LimitResult::failure(
            CasErrc::ResourceLimit, "limit allocation failed", operation);
    } catch (const std::exception& error) {
        return LimitResult::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

LimitResult limit_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction) {
    ComputationContext context;
    return limit_checked(expression, variable, point, direction, context);
}

LimitExpressionResult limit_expression_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction,
    ComputationContext& context) {
    auto result = limit_checked(
        expression, variable, point, direction, context);
    if (!result) return LimitExpressionResult::failure(result.error());
    const auto& outcome = result.value().value;
    if (const auto* finite = std::get_if<FiniteLimit>(&outcome)) {
        return LimitExpressionResult::success(finite->value);
    }
    if (std::holds_alternative<PositiveInfinityLimit>(outcome)) {
        return LimitExpressionResult::success(SymbolicExpr::infinity(1));
    }
    if (std::holds_alternative<NegativeInfinityLimit>(outcome)) {
        return LimitExpressionResult::success(SymbolicExpr::infinity(-1));
    }
    return LimitExpressionResult::failure(
        CasErrc::Inconclusive, "limit does not exist", "limit");
}

LimitExpressionResult limit_expression_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction) {
    ComputationContext context;
    return limit_expression_checked(
        expression, variable, point, direction, context);
}

} // namespace lamina
