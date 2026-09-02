#include "lsr_expr.hpp"
#include <exception>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "complex_analysis.hpp"
#include "matcher.hpp"
#include "symbolic_ast.hpp"
#include "internal/lsr_expr_common.hpp"
namespace lamina::lsr {
namespace detail::lsr_expr_common {
ExprResult expression_failure(CasErrc code, std::string message,
                              const char* operation) {
    return ExprResult::failure(code, std::move(message), operation);
}

ExprMatchResult expr_match_failure(CasErrc code, std::string message,
                                   const char* operation) {
    return ExprMatchResult::failure(code, std::move(message), operation);
}

ExprResult make_unary_math_expr(const ExprPtr& expression,
                                ComputationContext& context,
                                const char* function_name,
                                ExprPtr (*factory)(ExprPtr)) {
    auto step = context.consume_steps(1, kMathOperation);
    if (!step) return ExprResult::failure(step.error());
    if (!expression || !lamina::detail::node(expression)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  std::string(function_name) +
                                      " argument cannot be null",
                                  kMathOperation);
    }
    try {
        auto result = factory(expression);
        if (!result || !lamina::detail::node(result)) {
            return expression_failure(CasErrc::InternalInvariant,
                                      std::string(function_name) +
                                          " expression construction failed",
                                      kMathOperation);
        }
        return ExprResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  std::string(function_name) +
                                      " expression allocation failed",
                                  kMathOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kMathOperation);
    }
}


ExprResult make_unary_function_expr(const ExprPtr& expression,
                                    ComputationContext& context,
                                    const char* function_name,
                                    FunctionNode::FuncType type) {
    auto step = context.consume_steps(1, kMathOperation);
    if (!step) return ExprResult::failure(step.error());
    if (!expression || !lamina::detail::node(expression)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  std::string(function_name) +
                                      " argument cannot be null",
                                  kMathOperation);
    }
    try {
        auto node = lamina::detail::make_node<FunctionNode>(
            type, std::vector<std::shared_ptr<const SymbolicNode>>{
                      lamina::detail::node(expression)});
        return ExprResult::success(
            lamina::detail::make_expression_ptr(std::move(node)));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  std::string(function_name) +
                                      " expression allocation failed",
                                  kMathOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kMathOperation);
    }
}

ExprResult expr_from_complex_result(const ExpressionResult& result,
                                    const char* operation) {
    if (!result) {
        return ExprResult::failure(result.error());
    }
    if (!result.value() || !lamina::detail::node(result.value())) {
        return expression_failure(CasErrc::InternalInvariant,
                                  "complex expression result is null",
                                  operation);
    }
    return ExprResult::success(result.value());
}


bool is_reserved_symbol_name(const std::string& name) {
    return name == "I" || name == "pi" || name == "π" ||
           name == "e" || name == "phi";
}

bool is_imaginary_unit_name(const std::string& name) {
    return name == "I";
}

ExprResult require_dimensionless(const ExprPtr& expression,
                                 const char* function_name) {
    if (!expression || !lamina::detail::node(expression)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  std::string(function_name) + " argument cannot be null",
                                  kMathOperation);
    }
    auto dimension = dimension_of(*expression);
    if (!dimension) return ExprResult::failure(dimension.error());
    if (!dimension.value().is_dimensionless()) {
        return expression_failure(CasErrc::DimensionMismatch,
                                  std::string(function_name) +
                                      " requires a dimensionless argument",
                                  kMathOperation);
    }
    return ExprResult::success(expression);
}

ExprResult comparison_value(const ExprPtr& expression,
                            ComputationContext& context) {
    if (std::dynamic_pointer_cast<const QuantityNode>(lamina::detail::node(expression))) {
        return strip_unit(expression, UnitStripMode::BaseValue, context);
    }
    return ExprResult::success(expression);
}

} // namespace detail::lsr_expr_common
using namespace detail::lsr_expr_common;


ExprResult function(const std::string& name, std::vector<ExprPtr> arguments) {
    constexpr const char* operation = "lsr.function";
    if (name.empty()) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "function name cannot be empty", operation);
    }
    std::vector<std::shared_ptr<const SymbolicNode>> nodes;
    nodes.reserve(arguments.size());
    for (const auto& argument : arguments) {
        if (!argument || !lamina::detail::node(argument)) {
            return expression_failure(CasErrc::InvalidArgument,
                                      "function arguments cannot be null", operation);
        }
        nodes.push_back(lamina::detail::node(argument));
    }
    try {
        return ExprResult::success(lamina::detail::make_expression_ptr(
            lamina::detail::make_node<UninterpretedFunctionNode>(
                name, std::move(nodes))));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "function allocation failed", operation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument,
                                  error.what(), operation);
    }
}

ExprResult finite_set(std::vector<ExprPtr> elements) {
    ComputationContext context;
    return finite_set(std::move(elements), context);
}

ExprResult interval(ExprPtr lower, ExprPtr upper,
                    bool lower_closed, bool upper_closed) {
    ComputationContext context;
    return interval(lower, upper, lower_closed, upper_closed, context);
}

ExprResult relation(const ExprPtr& lhs, const ExprPtr& rhs, RelationOp op) {
    constexpr const char* operation = "lsr.relation";
    if (!lhs || !rhs) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "relation operands cannot be null", operation);
    }
    try {
        return ExprResult::success(lamina::detail::make_expression_ptr(
            lamina::detail::make_node<RelationalNode>(
                lamina::detail::node(lhs), lamina::detail::node(rhs), op)));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "relation allocation failed", operation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument,
                                  error.what(), operation);
    }
}

ExprResult ne(const ExprPtr& lhs, const ExprPtr& rhs) {
    return relation(lhs, rhs, RelationOp::NEQ);
}

ExprResult lt(const ExprPtr& lhs, const ExprPtr& rhs) {
    return relation(lhs, rhs, RelationOp::LT);
}

ExprResult le(const ExprPtr& lhs, const ExprPtr& rhs) {
    return relation(lhs, rhs, RelationOp::LEQ);
}

ExprResult gt(const ExprPtr& lhs, const ExprPtr& rhs) {
    return relation(lhs, rhs, RelationOp::GT);
}

ExprResult ge(const ExprPtr& lhs, const ExprPtr& rhs) {
    return relation(lhs, rhs, RelationOp::GEQ);
}

namespace {
ExprResult logical_expression(const ExprPtr& lhs, const ExprPtr& rhs,
                              LogicalNode::Op op, const char* operation) {
    if (!lhs || (op != LogicalNode::Op::Not && !rhs)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "logical operands cannot be null", operation);
    }
    try {
        return ExprResult::success(lamina::detail::make_expression_ptr(
            lamina::detail::make_node<LogicalNode>(
                lamina::detail::node(lhs),
                rhs ? lamina::detail::node(rhs) : nullptr, op)));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "logical expression allocation failed",
                                  operation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument,
                                  error.what(), operation);
    }
}
} // namespace

ExprResult logical_and(const ExprPtr& lhs, const ExprPtr& rhs) {
    return logical_expression(lhs, rhs, LogicalNode::Op::And, "lsr.logical_and");
}

ExprResult logical_or(const ExprPtr& lhs, const ExprPtr& rhs) {
    return logical_expression(lhs, rhs, LogicalNode::Op::Or, "lsr.logical_or");
}

ExprResult logical_not(const ExprPtr& expression) {
    return logical_expression(expression, nullptr, LogicalNode::Op::Not,
                              "lsr.logical_not");
}

ExprResult membership(const ExprPtr& element, const ExprPtr& set, bool negated) {
    ComputationContext context;
    auto result = member(element, set, context);
    if (!result || !negated) return result;
    return logical_not(result.value());
}

} // namespace lamina::lsr
