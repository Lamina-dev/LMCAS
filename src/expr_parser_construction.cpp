#include "internal/expr_parser.hpp"
#include "internal/expr_common.hpp"

#include <exception>
#include <utility>

namespace LMCAS::expr_detail {
namespace {

ExprResult invalid_construction(std::string message) {
    return expr_common::expression_failure(
        CasErrc::InvalidArgument, std::move(message), kParseOperation);
}

template <typename Factory>
ExprResult construct_expression(std::string allocation_message,
                                Factory&& factory) {
    try {
        auto node = std::forward<Factory>(factory)();
        return ExprResult::success(
            LMCAS::detail::make_expression_ptr(std::move(node)));
    } catch (const std::bad_alloc&) {
        return expr_common::expression_failure(
            CasErrc::ResourceLimit, std::move(allocation_message),
            kParseOperation);
    } catch (const std::exception& error) {
        return expr_common::expression_failure(
            CasErrc::InternalInvariant, error.what(), kParseOperation);
    }
}

}


ExprResult ExprParser::relational(
    const ExprPtr& lhs, const ExprPtr& rhs, RelationOp op) {
    if (!lhs || !rhs) {
        return invalid_construction("relational operands cannot be null");
    }
    return construct_expression(
        "relational expression allocation failed", [&] {
            return LMCAS::detail::make_node<RelationalNode>(
                LMCAS::detail::node(lhs), LMCAS::detail::node(rhs), op);
        });
}

ExprResult ExprParser::logical(
    const ExprPtr& lhs, const ExprPtr& rhs, LogicalNode::Op op) {
    if (!lhs || (op != LogicalNode::Op::Not && !rhs)) {
        return invalid_construction("logical operands cannot be null");
    }
    return construct_expression(
        "logical expression allocation failed", [&] {
            return LMCAS::detail::make_node<LogicalNode>(
                LMCAS::detail::node(lhs),
                rhs ? LMCAS::detail::node(rhs) : nullptr, op);
        });
}

ExprResult ExprParser::finite_set(const std::vector<ExprPtr>& elements) {
    for (const auto& element : elements) {
        if (!element || !LMCAS::detail::node(element)) {
            return invalid_construction("set elements cannot be null");
        }
    }
    return construct_expression(
        "set expression allocation failed", [&] {
            std::vector<std::shared_ptr<const SymbolicNode>> nodes;
            nodes.reserve(elements.size());
            for (const auto& element : elements) {
                nodes.push_back(LMCAS::detail::node(element));
            }
            return LMCAS::detail::make_node<FiniteSetNode>(std::move(nodes));
        });
}

ExprResult ExprParser::interval(
    const ExprPtr& lower, const ExprPtr& upper,
    bool lower_closed, bool upper_closed) {
    if (!lower || !upper) {
        return invalid_construction("interval bounds cannot be null");
    }
    return construct_expression(
        "interval expression allocation failed", [&] {
            return LMCAS::detail::make_node<IntervalNode>(
                LMCAS::detail::node(lower), LMCAS::detail::node(upper),
                lower_closed, upper_closed);
        });
}

ExprResult ExprParser::membership(
    const ExprPtr& element, const ExprPtr& set, bool negated) {
    if (!element || !set) {
        return invalid_construction("membership operands cannot be null");
    }
    return construct_expression(
        "membership expression allocation failed", [&] {
            std::shared_ptr<const SymbolicNode> node =
                LMCAS::detail::make_node<MembershipNode>(
                    LMCAS::detail::node(element), LMCAS::detail::node(set));
            if (negated) {
                node = LMCAS::detail::make_node<LogicalNode>(
                    std::move(node), nullptr, LogicalNode::Op::Not);
            }
            return node;
        });
}

ExprResult ExprParser::function_node(
    const std::string& name, const std::vector<ExprPtr>& arguments,
    FunctionNode::FuncType type) {
    for (const auto& argument : arguments) {
        if (!argument || !LMCAS::detail::node(argument)) {
            return invalid_construction("function arguments cannot be null");
        }
    }
    return construct_expression(
        name + " expression allocation failed", [&] {
            std::vector<std::shared_ptr<const SymbolicNode>> nodes;
            nodes.reserve(arguments.size());
            for (const auto& argument : arguments) {
                nodes.push_back(LMCAS::detail::node(argument));
            }
            return LMCAS::detail::make_node<FunctionNode>(
                type, std::move(nodes));
        });
}

}
