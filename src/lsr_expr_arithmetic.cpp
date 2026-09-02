#include "lsr_expr.hpp"
#include "complex_analysis.hpp"
#include "symbolic_ast.hpp"
#include "internal/lsr_expr_common.hpp"
#include <cmath>
#include <exception>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace lamina::lsr {

using namespace detail::lsr_expr_common;

ExprResult add(const ExprPtr& lhs,
               const ExprPtr& rhs,
               ComputationContext& context) {
    return quantity_add(lhs, rhs, context);
}

ExprResult add(const ExprPtr& lhs, const ExprPtr& rhs) {
    ComputationContext context;
    return add(lhs, rhs, context);
}

ExprResult mul(const ExprPtr& lhs,
               const ExprPtr& rhs,
               ComputationContext& context) {
    return quantity_multiply(lhs, rhs, context);
}

ExprResult mul(const ExprPtr& lhs, const ExprPtr& rhs) {
    ComputationContext context;
    return mul(lhs, rhs, context);
}

ExprResult div(const ExprPtr& numerator,
               const ExprPtr& denominator,
               ComputationContext& context) {
    return quantity_divide(numerator, denominator, context);
}

ExprResult div(const ExprPtr& numerator, const ExprPtr& denominator) {
    ComputationContext context;
    return div(numerator, denominator, context);
}

ExprResult neg(const ExprPtr& expression, ComputationContext& context) {
    auto minus_one = integer(-1);
    if (!minus_one) return ExprResult::failure(minus_one.error());
    return mul(minus_one.value(), expression, context);
}

ExprResult neg(const ExprPtr& expression) {
    ComputationContext context;
    return neg(expression, context);
}

ExprResult sub(const ExprPtr& lhs,
               const ExprPtr& rhs,
               ComputationContext& context) {
    return quantity_subtract(lhs, rhs, context);
}

ExprResult sub(const ExprPtr& lhs, const ExprPtr& rhs) {
    ComputationContext context;
    return sub(lhs, rhs, context);
}

ExprResult eq(const ExprPtr& lhs,
              const ExprPtr& rhs,
              ComputationContext& context) {
    if (!lhs || !rhs) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "eq operands cannot be null", kExprOperation);
    }
    auto left_dimension = dimension_of(*lhs);
    if (!left_dimension) return ExprResult::failure(left_dimension.error());
    auto right_dimension = dimension_of(*rhs);
    if (!right_dimension) return ExprResult::failure(right_dimension.error());
    if (left_dimension.value() != right_dimension.value()) return integer(0);
    auto left = comparison_value(lhs, context);
    if (!left) return left;
    auto right = comparison_value(rhs, context);
    if (!right) return right;
    return make_binary_expr(left.value(), right.value(), context, "eq", SymbolicExpr::eq);
}

ExprResult eq(const ExprPtr& lhs, const ExprPtr& rhs) {
    ComputationContext context;
    return eq(lhs, rhs, context);
}

ExprResult sqrt(const ExprPtr& expression, ComputationContext& context) {
    auto exponent = rational(Rational(1, 2));
    if (!exponent) return exponent;
    return quantity_power(expression, exponent.value(), context);
}

ExprResult sqrt(const ExprPtr& expression) {
    ComputationContext context;
    return sqrt(expression, context);
}

ExprResult pow(const ExprPtr& base,
               const ExprPtr& exponent,
               ComputationContext& context) {
    return quantity_power(base, exponent, context);
}

ExprResult pow(const ExprPtr& base, const ExprPtr& exponent) {
    ComputationContext context;
    return pow(base, exponent, context);
}

ExprResult sin(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "sin");
    if (!valid) return valid;
    return make_unary_math_expr(expression, context, "sin",
                                SymbolicExpr::sin);
}

ExprResult sin(const ExprPtr& expression) {
    ComputationContext context;
    return sin(expression, context);
}

ExprResult cos(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "cos");
    if (!valid) return valid;
    return make_unary_math_expr(expression, context, "cos",
                                SymbolicExpr::cos);
}

ExprResult cos(const ExprPtr& expression) {
    ComputationContext context;
    return cos(expression, context);
}

ExprResult tan(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "tan");
    if (!valid) return valid;
    return make_unary_math_expr(expression, context, "tan",
                                SymbolicExpr::tan);
}

ExprResult tan(const ExprPtr& expression) {
    ComputationContext context;
    return tan(expression, context);
}

ExprResult asin(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "asin");
    if (!valid) return valid;
    return make_unary_function_expr(expression, context, "asin",
                                    FunctionNode::FuncType::ArcSin);
}

ExprResult asin(const ExprPtr& expression) {
    ComputationContext context;
    return asin(expression, context);
}

ExprResult acos(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "acos");
    if (!valid) return valid;
    return make_unary_function_expr(expression, context, "acos",
                                    FunctionNode::FuncType::ArcCos);
}

ExprResult acos(const ExprPtr& expression) {
    ComputationContext context;
    return acos(expression, context);
}

ExprResult atan(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "atan");
    if (!valid) return valid;
    return make_unary_function_expr(expression, context, "atan",
                                    FunctionNode::FuncType::ArcTan);
}

ExprResult atan(const ExprPtr& expression) {
    ComputationContext context;
    return atan(expression, context);
}

ExprResult exp(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "exp");
    if (!valid) return valid;
    return make_unary_math_expr(expression, context, "exp",
                                SymbolicExpr::exp);
}

ExprResult exp(const ExprPtr& expression) {
    ComputationContext context;
    return exp(expression, context);
}

ExprResult log(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "log");
    if (!valid) return valid;
    return make_unary_math_expr(expression, context, "log",
                                SymbolicExpr::ln);
}

ExprResult log(const ExprPtr& expression) {
    ComputationContext context;
    return log(expression, context);
}

ExprResult log10(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "log10");
    if (!valid) return valid;
    auto step = context.consume_steps(1, kMathOperation);
    if (!step) return ExprResult::failure(step.error());
    if (!expression || !lamina::detail::node(expression)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "log10 argument cannot be null",
                                  kMathOperation);
    }
    try {
        auto numerator = SymbolicExpr::ln(expression);
        auto denominator = SymbolicExpr::ln(SymbolicExpr::number(10));
        auto result = SymbolicExpr::divide(numerator, denominator);
        if (!result || !lamina::detail::node(result)) {
            return expression_failure(CasErrc::InternalInvariant,
                                      "log10 expression construction failed",
                                      kMathOperation);
        }
        return ExprResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "log10 expression allocation failed",
                                  kMathOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kMathOperation);
    }
}

ExprResult log10(const ExprPtr& expression) {
    ComputationContext context;
    return log10(expression, context);
}

ExprResult floor(const ExprPtr& expression, ComputationContext& context) {
    return make_unary_function_expr(expression, context, "floor",
                                    FunctionNode::FuncType::Floor);
}

ExprResult floor(const ExprPtr& expression) {
    ComputationContext context;
    return floor(expression, context);
}

ExprResult ceil(const ExprPtr& expression, ComputationContext& context) {
    return make_unary_function_expr(expression, context, "ceil",
                                    FunctionNode::FuncType::Ceil);
}

ExprResult ceil(const ExprPtr& expression) {
    ComputationContext context;
    return ceil(expression, context);
}

ExprResult round(const ExprPtr& expression, ComputationContext& context) {
    return make_unary_function_expr(expression, context, "round",
                                    FunctionNode::FuncType::Round);
}

ExprResult round(const ExprPtr& expression) {
    ComputationContext context;
    return round(expression, context);
}

ExprResult clamp(const ExprPtr& expression,
                 const ExprPtr& lower,
                 const ExprPtr& upper,
                 ComputationContext& context) {
    auto step = context.consume_steps(1, kMathOperation);
    if (!step) return ExprResult::failure(step.error());
    if (!expression || !lamina::detail::node(expression) ||
        !lower || !lamina::detail::node(lower) ||
        !upper || !lamina::detail::node(upper)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "clamp arguments cannot be null",
                                  kMathOperation);
    }
    try {
        auto max_node = lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Max,
            std::vector<std::shared_ptr<const SymbolicNode>>{
                lamina::detail::node(expression), lamina::detail::node(lower)});
        auto min_node = lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Min,
            std::vector<std::shared_ptr<const SymbolicNode>>{
                std::move(max_node), lamina::detail::node(upper)});
        return ExprResult::success(
            lamina::detail::make_expression_ptr(std::move(min_node)));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "clamp expression allocation failed",
                                  kMathOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kMathOperation);
    }
}

ExprResult clamp(const ExprPtr& expression,
                 const ExprPtr& lower,
                 const ExprPtr& upper) {
    ComputationContext context;
    return clamp(expression, lower, upper, context);
}

ExprResult real(const ExprPtr& expression, ComputationContext& context) {
    return expr_from_complex_result(real_part_checked(expression, context),
                                    kRealOperation);
}

ExprResult real(const ExprPtr& expression) {
    ComputationContext context;
    return real(expression, context);
}

ExprResult imag(const ExprPtr& expression, ComputationContext& context) {
    return expr_from_complex_result(imag_part_checked(expression, context),
                                    kImagOperation);
}

ExprResult imag(const ExprPtr& expression) {
    ComputationContext context;
    return imag(expression, context);
}

ExprResult conj(const ExprPtr& expression, ComputationContext& context) {
    return expr_from_complex_result(conjugate_checked(expression, context),
                                    kConjOperation);
}

ExprResult conj(const ExprPtr& expression) {
    ComputationContext context;
    return conj(expression, context);
}

ExprResult abs(const ExprPtr& expression, ComputationContext& context) {
    auto step = context.consume_steps(1, kAbsOperation);
    if (!step) return ExprResult::failure(step.error());
    auto re = real(expression, context);
    if (!re) return re;
    auto im = imag(expression, context);
    if (!im) return im;
    try {
        auto re_squared = SymbolicExpr::power(re.value(), SymbolicExpr::number(2));
        auto im_squared = SymbolicExpr::power(im.value(), SymbolicExpr::number(2));
        auto sum = SymbolicExpr::add(re_squared, im_squared);
        auto result = SymbolicExpr::sqrt(sum)->simplify();
        if (!result || !lamina::detail::node(result)) {
            return expression_failure(CasErrc::InternalInvariant,
                                      "complex absolute value construction failed",
                                      kAbsOperation);
        }
        return ExprResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "complex absolute value allocation failed",
                                  kAbsOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kAbsOperation);
    }
}

ExprResult abs(const ExprPtr& expression) {
    ComputationContext context;
    return abs(expression, context);
}
} // namespace lamina::lsr
