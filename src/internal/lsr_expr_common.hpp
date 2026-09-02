#pragma once

#include "lsr_expr.hpp"
#include "symbolic_ast.hpp"
#include <exception>
#include <memory>
#include <string>
#include <utility>

namespace lamina::lsr::detail::lsr_expr_common {

inline constexpr const char* kSymOperation = "lsr.sym";
inline constexpr const char* kIntegerOperation = "lsr.integer";
inline constexpr const char* kRationalOperation = "lsr.rational";
inline constexpr const char* kApproxOperation = "lsr.approx_real";
inline constexpr const char* kConstantOperation = "lsr.constant";
inline constexpr const char* kComplexOperation = "lsr.complex";
inline constexpr const char* kExprOperation = "lsr.expr_op";
inline constexpr const char* kMathOperation = "lsr.math";
inline constexpr const char* kRealOperation = "lsr.real";
inline constexpr const char* kImagOperation = "lsr.imag";
inline constexpr const char* kConjOperation = "lsr.conj";
inline constexpr const char* kAbsOperation = "lsr.abs";
inline constexpr const char* kSimplifyOperation = "lsr.simplify";
inline constexpr const char* kExpandOperation = "lsr.expand";
inline constexpr const char* kDifferentiateOperation = "lsr.differentiate";
inline constexpr const char* kSubstituteOperation = "lsr.substitute";
inline constexpr const char* kExprMatchOperation = "lsr.expr_match";

ExprResult expression_failure(CasErrc code, std::string message,
                              const char* operation);
ExprMatchResult expr_match_failure(CasErrc code, std::string message,
                                   const char* operation);

ExprResult make_unary_math_expr(const ExprPtr& expression,
                                ComputationContext& context,
                                const char* function_name,
                                ExprPtr (*factory)(ExprPtr));

ExprResult make_unary_function_expr(const ExprPtr& expression,
                                    ComputationContext& context,
                                    const char* function_name,
                                    FunctionNode::FuncType type);

ExprResult expr_from_complex_result(const ExpressionResult& result,
                                    const char* operation);

bool is_reserved_symbol_name(const std::string& name);
bool is_imaginary_unit_name(const std::string& name);
ExprResult require_dimensionless(const ExprPtr& expression,
                                 const char* function_name);
ExprResult comparison_value(const ExprPtr& expression,
                            ComputationContext& context);

template <typename Factory>
ExprResult make_binary_math_expr(const ExprPtr& lhs,
                                 const ExprPtr& rhs,
                                 ComputationContext& context,
                                 const char* function_name,
                                 Factory&& factory) {
    auto step = context.consume_steps(1, kMathOperation);
    if (!step) return ExprResult::failure(step.error());
    if (!lhs || !lamina::detail::node(lhs) ||
        !rhs || !lamina::detail::node(rhs)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  std::string(function_name) +
                                      " arguments cannot be null",
                                  kMathOperation);
    }
    try {
        auto result = std::forward<Factory>(factory)(lhs, rhs);
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

template <typename Factory>
ExprResult make_binary_expr(const ExprPtr& lhs,
                            const ExprPtr& rhs,
                            ComputationContext& context,
                            const char* operation_name,
                            Factory&& factory) {
    auto step = context.consume_steps(1, kExprOperation);
    if (!step) return ExprResult::failure(step.error());
    if (!lhs || !lamina::detail::node(lhs) ||
        !rhs || !lamina::detail::node(rhs)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  std::string(operation_name) +
                                      " operands cannot be null",
                                  kExprOperation);
    }
    try {
        auto result = std::forward<Factory>(factory)(lhs, rhs);
        if (!result || !lamina::detail::node(result)) {
            return expression_failure(CasErrc::InternalInvariant,
                                      std::string(operation_name) +
                                          " expression construction failed",
                                      kExprOperation);
        }
        return ExprResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  std::string(operation_name) +
                                      " expression allocation failed",
                                  kExprOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kExprOperation);
    }
}

template <typename Transform>
ExprResult checked_transform_expr(const ExprPtr& expression,
                                  ComputationContext& context,
                                  const char* operation,
                                  const char* transform_name,
                                  Transform&& transform) {
    auto step = context.consume_steps(1, operation);
    if (!step) return ExprResult::failure(step.error());
    if (!expression || !lamina::detail::node(expression)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  std::string(transform_name) +
                                      " argument cannot be null",
                                  operation);
    }
    try {
        auto result = std::forward<Transform>(transform)(*expression);
        if (!result || !lamina::detail::node(result)) {
            return expression_failure(CasErrc::InternalInvariant,
                                      std::string(transform_name) +
                                          " returned null",
                                      operation);
        }
        return ExprResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  std::string(transform_name) +
                                      " allocation failed",
                                  operation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  operation);
    }
}

} // namespace lamina::lsr::detail::lsr_expr_common
