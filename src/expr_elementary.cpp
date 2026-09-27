#include "expr.hpp"
#include "complex_analysis.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/expr_common.hpp"
#include <cmath>
#include <exception>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace LMCAS {

using namespace expr_detail::expr_common;

ExprResult sqrt(const ExprPtr& expression, ComputationContext& context) {
    auto exponent = rational(Rational(1, 2));
    if (!exponent) { return exponent; }
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
    if (!valid) { return valid; }
    return make_unary_math_expr(expression, context, "sin",
                                SymbolicExpr::sin);
}

ExprResult sin(const ExprPtr& expression) {
    ComputationContext context;
    return sin(expression, context);
}

ExprResult cos(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "cos");
    if (!valid) { return valid; }
    return make_unary_math_expr(expression, context, "cos",
                                SymbolicExpr::cos);
}

ExprResult cos(const ExprPtr& expression) {
    ComputationContext context;
    return cos(expression, context);
}

ExprResult tan(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "tan");
    if (!valid) { return valid; }
    return make_unary_math_expr(expression, context, "tan",
                                SymbolicExpr::tan);
}

ExprResult tan(const ExprPtr& expression) {
    ComputationContext context;
    return tan(expression, context);
}

ExprResult asin(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "asin");
    if (!valid) { return valid; }
    return make_unary_function_expr(expression, context, "asin",
                                    FunctionNode::FuncType::ArcSin);
}

ExprResult asin(const ExprPtr& expression) {
    ComputationContext context;
    return asin(expression, context);
}

ExprResult acos(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "acos");
    if (!valid) { return valid; }
    return make_unary_function_expr(expression, context, "acos",
                                    FunctionNode::FuncType::ArcCos);
}

ExprResult acos(const ExprPtr& expression) {
    ComputationContext context;
    return acos(expression, context);
}

ExprResult atan(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "atan");
    if (!valid) { return valid; }
    return make_unary_function_expr(expression, context, "atan",
                                    FunctionNode::FuncType::ArcTan);
}

ExprResult atan(const ExprPtr& expression) {
    ComputationContext context;
    return atan(expression, context);
}

ExprResult exp(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "exp");
    if (!valid) { return valid; }
    return make_unary_math_expr(expression, context, "exp",
                                SymbolicExpr::exp);
}

ExprResult exp(const ExprPtr& expression) {
    ComputationContext context;
    return exp(expression, context);
}

ExprResult log(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "log");
    if (!valid) { return valid; }
    return make_unary_math_expr(expression, context, "log",
                                SymbolicExpr::ln);
}

ExprResult log(const ExprPtr& expression) {
    ComputationContext context;
    return log(expression, context);
}

ExprResult log10(const ExprPtr& expression, ComputationContext& context) {
    auto valid = require_dimensionless(expression, "log10");
    if (!valid) { return valid; }
    auto step = context.consume_steps(1, kMathOperation);
    if (!step) { return ExprResult::failure(step.error()); }
    if (!expression || !LMCAS::detail::node(expression)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "log10 argument cannot be null",
                                  kMathOperation);
    }
    try {
        auto numerator = SymbolicExpr::ln(expression);
        auto denominator = SymbolicExpr::ln(SymbolicExpr::number(10));
        auto result = SymbolicExpr::divide(numerator, denominator);
        if (!result || !LMCAS::detail::node(result)) {
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
        return expression_failure(CasErrc::InternalInvariant, error.what(),
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
    if (!step) { return ExprResult::failure(step.error()); }
    if (!expression || !LMCAS::detail::node(expression) ||
        !lower || !LMCAS::detail::node(lower) ||
        !upper || !LMCAS::detail::node(upper)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "clamp arguments cannot be null",
                                  kMathOperation);
    }
    try {
        auto max_node = LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Max,
            std::vector<std::shared_ptr<const SymbolicNode>>{
                LMCAS::detail::node(expression), LMCAS::detail::node(lower)});
        auto min_node = LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Min,
            std::vector<std::shared_ptr<const SymbolicNode>>{
                std::move(max_node), LMCAS::detail::node(upper)});
        return ExprResult::success(
            LMCAS::detail::make_expression_ptr(std::move(min_node)));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "clamp expression allocation failed",
                                  kMathOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InternalInvariant, error.what(),
                                  kMathOperation);
    }
}

ExprResult clamp(const ExprPtr& expression,
                 const ExprPtr& lower,
                 const ExprPtr& upper) {
    ComputationContext context;
    return clamp(expression, lower, upper, context);
}

}
