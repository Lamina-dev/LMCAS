#include "lsr_expr.hpp"
#include <exception>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "complex_analysis.hpp"
#include "matcher.hpp"
#include "symbolic_ast.hpp"
namespace lamina::lsr {
namespace {
constexpr const char* kSymOperation = "lsr.sym";
constexpr const char* kIntegerOperation = "lsr.integer";
constexpr const char* kRationalOperation = "lsr.rational";
constexpr const char* kApproxOperation = "lsr.approx_real";
constexpr const char* kConstantOperation = "lsr.constant";
constexpr const char* kComplexOperation = "lsr.complex";
constexpr const char* kExprOperation = "lsr.expr_op";
constexpr const char* kMathOperation = "lsr.math";
constexpr const char* kRealOperation = "lsr.real";
constexpr const char* kImagOperation = "lsr.imag";
constexpr const char* kConjOperation = "lsr.conj";
constexpr const char* kAbsOperation = "lsr.abs";
constexpr const char* kSimplifyOperation = "lsr.simplify";
constexpr const char* kExpandOperation = "lsr.expand";
constexpr const char* kDifferentiateOperation = "lsr.differentiate";
constexpr const char* kSubstituteOperation = "lsr.substitute";
constexpr const char* kExprMatchOperation = "lsr.expr_match";
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

template <typename Factory>
ExprResult make_binary_math_expr(const ExprPtr& lhs,
                                 const ExprPtr& rhs,
                                 ComputationContext& context,
                                 const char* function_name,
                                 Factory factory) {
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
        auto result = factory(lhs, rhs);
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
                            Factory factory) {
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
        auto result = factory(lhs, rhs);
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

template <typename Transform>
ExprResult checked_transform_expr(const ExprPtr& expression,
                                  ComputationContext& context,
                                  const char* operation,
                                  const char* transform_name,
                                  Transform transform) {
    auto step = context.consume_steps(1, operation);
    if (!step) return ExprResult::failure(step.error());
    if (!expression || !lamina::detail::node(expression)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  std::string(transform_name) +
                                      " argument cannot be null",
                                  operation);
    }
    try {
        auto result = transform(*expression);
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

bool is_reserved_symbol_name(const std::string& name) {
    return name == "i" || name == "I" || name == "pi" || name == "π" ||
           name == "e" || name == "phi";
}

bool is_imaginary_unit_name(const std::string& name) {
    return name == "i" || name == "I";
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

} // namespace
ExprResult sym(const std::string& name) {
    if (name.empty()) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "symbol name cannot be empty", kSymOperation);
    }
    if (is_reserved_symbol_name(name)) {
        if (is_imaginary_unit_name(name)) {
            return expression_failure(CasErrc::InvalidArgument,
                                      "imaginary unit symbol is reserved",
                                      kSymOperation);
        }
        return expression_failure(CasErrc::InvalidArgument,
                                  "reserved mathematical constants cannot be shadowed",
                                  kSymOperation);
    }
    try {
        auto expression = SymbolicExpr::variable(name);
        if (!expression) {
            return expression_failure(CasErrc::InternalInvariant,
                                      "symbol factory returned null", kSymOperation);
        }
        return ExprResult::success(std::move(expression));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "symbol allocation failed", kSymOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kSymOperation);
    }
}

ExprResult integer(long long value) {
    try {
        return ExprResult::success(SymbolicExpr::number(value));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "integer allocation failed", kIntegerOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kIntegerOperation);
    }
}

ExprResult integer(const BigInt& value) {
    try {
        return ExprResult::success(SymbolicExpr::number(value));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "integer allocation failed", kIntegerOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kIntegerOperation);
    }
}

ExprResult rational(const Rational& value) {
    try {
        return ExprResult::success(SymbolicExpr::number(value));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "rational allocation failed", kRationalOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kRationalOperation);
    }
}

ExprResult approx_real(double value) {
    if (!std::isfinite(value)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "approximate real must be finite", kApproxOperation);
    }
    try {
        return ExprResult::success(SymbolicExpr::number(value));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "approximate real allocation failed", kApproxOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(), kApproxOperation);
    }
}

ExprResult constant_symbol(const char* name) {
    try {
        auto expression = SymbolicExpr::variable(name);
        if (!expression || !lamina::detail::node(expression)) {
            return expression_failure(CasErrc::InternalInvariant,
                                      "constant factory returned null",
                                      kConstantOperation);
        }
        return ExprResult::success(std::move(expression));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "constant allocation failed",
                                  kConstantOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kConstantOperation);
    }
}

ExprResult pi() {
    return constant_symbol("pi");
}

ExprResult e() {
    return constant_symbol("e");
}

ExprResult phi() {
    return constant_symbol("phi");
}

ExprResult i() {
    return imaginary_unit();
}

ExprResult I() {
    return imaginary_unit();
}

ExprResult imaginary_unit() {
    auto zero = integer(0);
    if (!zero) return ExprResult::failure(zero.error());
    auto one = integer(1);
    if (!one) return ExprResult::failure(one.error());
    return complex(zero.value(), one.value());
}

ExprResult complex(ExprPtr real, ExprPtr imag) {
    if (!real || !imag) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "complex expression parts cannot be null",
                                  kComplexOperation);
    }
    try {
        auto node = SymbolicFactory::create_complex(lamina::detail::node(real),
                                                    lamina::detail::node(imag));
        return ExprResult::success(lamina::detail::make_expression_ptr(std::move(node)));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "complex expression allocation failed",
                                  kComplexOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kComplexOperation);
    }
}

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

ExprResult simplify(const ExprPtr& expression, ComputationContext& context) {
    return checked_transform_expr(
        expression, context, kSimplifyOperation, "simplify",
        [](const SymbolicExpr& value) { return value.simplify(); });
}

ExprResult simplify(const ExprPtr& expression) {
    ComputationContext context;
    return simplify(expression, context);
}

ExprResult expand(const ExprPtr& expression, ComputationContext& context) {
    return checked_transform_expr(
        expression, context, kExpandOperation, "expand",
        [](const SymbolicExpr& value) { return value.expand(); });
}

ExprResult expand(const ExprPtr& expression) {
    ComputationContext context;
    return expand(expression, context);
}

ExprResult differentiate(const ExprPtr& expression,
                         const std::string& variable,
                         ComputationContext& context) {
    if (variable.empty()) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "differentiate variable cannot be empty",
                                  kDifferentiateOperation);
    }
    return checked_transform_expr(
        expression, context, kDifferentiateOperation, "differentiate",
        [&variable](const SymbolicExpr& value) {
            return value.differentiate(variable);
        });
}

ExprResult differentiate(const ExprPtr& expression,
                         const std::string& variable) {
    ComputationContext context;
    return differentiate(expression, variable, context);
}

ExprResult substitute(const ExprPtr& expression,
                      const std::string& variable,
                      const ExprPtr& value,
                      ComputationContext& context) {
    auto step = context.consume_steps(1, kSubstituteOperation);
    if (!step) return ExprResult::failure(step.error());
    if (!expression) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "expression cannot be null",
                                  kSubstituteOperation);
    }
    if (variable.empty()) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "substitution variable cannot be empty",
                                  kSubstituteOperation);
    }
    if (!value) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "substitution value cannot be null",
                                  kSubstituteOperation);
    }
    try {
        auto result = expression->substitute(variable, value);
        if (!result || !lamina::detail::node(result)) {
            return expression_failure(CasErrc::InternalInvariant,
                                      "substitution returned an empty expression",
                                      kSubstituteOperation);
        }
        return ExprResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "substitution allocation failed",
                                  kSubstituteOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kSubstituteOperation);
    }
}

ExprResult substitute(const ExprPtr& expression,
                      const std::string& variable,
                      const ExprPtr& value) {
    ComputationContext context;
    return substitute(expression, variable, value, context);
}

ExprMatchResult expr_match(const ExprPtr& pattern,
                           const ExprPtr& target,
                           const std::vector<std::string>& wildcards,
                           ComputationContext& context) {
    auto step = context.consume_steps(1, kExprMatchOperation);
    if (!step) return ExprMatchResult::failure(step.error());
    if (!pattern || !lamina::detail::node(pattern)) {
        return expr_match_failure(CasErrc::InvalidArgument,
                                  "match pattern cannot be null",
                                  kExprMatchOperation);
    }
    if (!target || !lamina::detail::node(target)) {
        return expr_match_failure(CasErrc::InvalidArgument,
                                  "match target cannot be null",
                                  kExprMatchOperation);
    }

    std::unordered_set<std::string> wildcard_set;
    wildcard_set.reserve(wildcards.size());
    for (const auto& wildcard : wildcards) {
        if (wildcard.empty()) {
            return expr_match_failure(CasErrc::InvalidArgument,
                                      "wildcard names cannot be empty",
                                      kExprMatchOperation);
        }
        if (!wildcard_set.insert(wildcard).second) {
            return expr_match_failure(CasErrc::InvalidArgument,
                                      "wildcard names must be unique",
                                      kExprMatchOperation);
        }
    }

    try {
        MatchMap raw_bindings;
        const bool matched = Matcher::match(
            *pattern, *target, wildcard_set, raw_bindings);
        if (!matched) {
            return ExprMatchResult::success(ExprMatch{false, {}});
        }

        std::vector<std::string> names;
        names.reserve(raw_bindings.size());
        for (const auto& binding : raw_bindings) {
            names.push_back(binding.first);
        }
        std::sort(names.begin(), names.end());

        ExprMatch result;
        result.matched = true;
        result.bindings.reserve(names.size());
        for (const auto& name : names) {
            const auto& value = raw_bindings.at(name);
            auto node = lamina::detail::node(value);
            if (!node) {
                return expr_match_failure(CasErrc::InternalInvariant,
                                          "matcher produced a null binding",
                                          kExprMatchOperation);
            }
            result.bindings.push_back(
                ExprMatchBinding{name,
                                 lamina::detail::make_expression_ptr(node)});
        }
        return ExprMatchResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return expr_match_failure(CasErrc::ResourceLimit,
                                  "expression matching allocation failed",
                                  kExprMatchOperation);
    } catch (const std::exception& error) {
        return expr_match_failure(CasErrc::InternalInvariant, error.what(),
                                  kExprMatchOperation);
    }
}

ExprMatchResult expr_match(const ExprPtr& pattern,
                           const ExprPtr& target,
                           const std::vector<std::string>& wildcards) {
    ComputationContext context;
    return expr_match(pattern, target, wildcards, context);
}

} // namespace lamina::lsr
