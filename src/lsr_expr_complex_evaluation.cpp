#include "lsr_expr.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <utility>

#include "lsr_expr_internal.hpp"
#include "symbolic_ast.hpp"

namespace lamina::lsr {
namespace {

constexpr const char* kEvalfOperation = "lsr.evalf";

Result<ApproxComplex> complex_failure(CasErrc code, std::string message,
                                      const char* operation) {
    return Result<ApproxComplex>::failure(code, std::move(message), operation);
}

Result<ApproxComplex> eval_complex_failure(const CasError& error) {
    return complex_failure(error.code, error.message, kEvalComplexOperation);
}

ApproxReal approx_part(double value) {
    ApproxReal part;
    part.value = value;
    part.absolute_error = std::numeric_limits<double>::epsilon() *
                          std::max(1.0, std::abs(value)) * 4.0;
    part.status = NumericStatus::Finite;
    return part;
}

ApproxComplex approx_complex(double real, double imag) {
    return ApproxComplex{approx_part(real), approx_part(imag)};
}

Result<ApproxComplex> checked_complex(double real, double imag,
                                      const char* operation) {
    if (!std::isfinite(real) || !std::isfinite(imag)) {
        return complex_failure(CasErrc::NumericFailure,
                               "complex evaluation produced a non-finite component",
                               operation);
    }
    return Result<ApproxComplex>::success(approx_complex(real, imag));
}

Result<ApproxComplex> real_to_complex(const Result<ApproxReal>& real) {
    if (!real) {
        return eval_complex_failure(real.error());
    }
    if (!real.value().is_finite()) {
        return complex_failure(CasErrc::NumericFailure,
                               "complex evaluation requires finite real components",
                               kEvalComplexOperation);
    }
    return Result<ApproxComplex>::success(
        ApproxComplex{real.value(), approx_part(0.0)});
}

Result<ApproxComplex> add_complex(const ApproxComplex& lhs,
                                  const ApproxComplex& rhs) {
    auto result = checked_complex(lhs.real.value + rhs.real.value,
                                  lhs.imag.value + rhs.imag.value,
                                  kEvalComplexOperation);
    if (!result) return result;
    result.value().real.absolute_error +=
        lhs.real.absolute_error + rhs.real.absolute_error;
    result.value().imag.absolute_error +=
        lhs.imag.absolute_error + rhs.imag.absolute_error;
    return result;
}

Result<ApproxComplex> multiply_complex(const ApproxComplex& lhs,
                                       const ApproxComplex& rhs) {
    const double a = lhs.real.value;
    const double b = lhs.imag.value;
    const double c = rhs.real.value;
    const double d = rhs.imag.value;
    auto result = checked_complex(a * c - b * d, a * d + b * c,
                                  kEvalComplexOperation);
    if (!result) return result;
    result.value().real.absolute_error +=
        std::abs(c) * lhs.real.absolute_error +
        std::abs(a) * rhs.real.absolute_error +
        std::abs(d) * lhs.imag.absolute_error +
        std::abs(b) * rhs.imag.absolute_error;
    result.value().imag.absolute_error +=
        std::abs(d) * lhs.real.absolute_error +
        std::abs(a) * rhs.imag.absolute_error +
        std::abs(c) * lhs.imag.absolute_error +
        std::abs(b) * rhs.real.absolute_error;
    return result;
}

Result<ApproxComplex> divide_complex(const ApproxComplex& lhs,
                                     const ApproxComplex& rhs) {
    const double c = rhs.real.value;
    const double d = rhs.imag.value;
    const double denom = c * c + d * d;
    if (denom == 0.0 || !std::isfinite(denom)) {
        return complex_failure(CasErrc::DomainError,
                               "complex division by zero or overflow",
                               kEvalComplexOperation);
    }
    return checked_complex((lhs.real.value * c + lhs.imag.value * d) / denom,
                           (lhs.imag.value * c - lhs.real.value * d) / denom,
                           kEvalComplexOperation);
}

bool is_integer_double(double value) {
    return std::isfinite(value) && std::floor(value) == value;
}

Result<ApproxComplex> evaluate_complex_node(
    const std::shared_ptr<const SymbolicNode>& node,
    const NumericBindings& bindings,
    ComputationContext& context) {
    auto step = context.consume_steps(1, kEvalComplexOperation);
    if (!step) return Result<ApproxComplex>::failure(step.error());
    if (!node) {
        return complex_failure(CasErrc::InvalidArgument,
                               "expression contains a null node",
                               kEvalComplexOperation);
    }

    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real = evaluate_numeric(
            *lamina::detail::make_expression_ptr(complex_node->real()),
            bindings, context);
        if (!real) return eval_complex_failure(real.error());
        auto imag = evaluate_numeric(
            *lamina::detail::make_expression_ptr(complex_node->imag()),
            bindings, context);
        if (!imag) return eval_complex_failure(imag.error());
        if (!real.value().is_finite() || !imag.value().is_finite()) {
            return complex_failure(CasErrc::NumericFailure,
                                   "complex components must be finite",
                                   kEvalComplexOperation);
        }
        return Result<ApproxComplex>::success(
            ApproxComplex{real.value(), imag.value()});
    }

    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (is_imaginary_unit_name(variable->name())) {
            return Result<ApproxComplex>::success(approx_complex(0.0, 1.0));
        }
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        auto sum = Result<ApproxComplex>::success(approx_complex(0.0, 0.0));
        for (const auto& operand : add->operands()) {
            auto term = evaluate_complex_node(operand, bindings, context);
            if (!term) return term;
            sum = add_complex(sum.value(), term.value());
            if (!sum) return sum;
        }
        return sum;
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        auto product = Result<ApproxComplex>::success(approx_complex(1.0, 0.0));
        for (const auto& operand : multiply->operands()) {
            auto factor = evaluate_complex_node(operand, bindings, context);
            if (!factor) return factor;
            product = multiply_complex(product.value(), factor.value());
            if (!product) return product;
        }
        return product;
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base = evaluate_complex_node(power->base(), bindings, context);
        if (!base) return base;
        auto exponent = evaluate_numeric(
            *lamina::detail::make_expression_ptr(power->exponent()),
            bindings, context);
        if (!exponent) return eval_complex_failure(exponent.error());
        if (!exponent.value().is_finite() ||
            !std::isfinite(exponent.value().value)) {
            return complex_failure(CasErrc::NumericFailure,
                                   "complex power exponent must be finite",
                                   kEvalComplexOperation);
        }
        const double exponent_value = exponent.value().value;
        if (!is_integer_double(exponent_value) ||
            std::abs(exponent_value) > 64.0) {
            return complex_failure(CasErrc::UnsupportedExpression,
                                   "complex evaluation only supports integer powers with |n| <= 64",
                                   kEvalComplexOperation);
        }
        int exponent_int = static_cast<int>(exponent_value);
        auto result = Result<ApproxComplex>::success(approx_complex(1.0, 0.0));
        for (int i = 0; i < std::abs(exponent_int); ++i) {
            result = multiply_complex(result.value(), base.value());
            if (!result) return result;
        }
        if (exponent_int < 0) {
            result = divide_complex(approx_complex(1.0, 0.0), result.value());
        }
        return result;
    }

    return real_to_complex(
        evaluate_numeric(*lamina::detail::make_expression_ptr(node),
                         bindings, context));
}

} // namespace

Result<ApproxReal> evalf(const SymbolicExpr& expression,
                         const NumericBindings& bindings,
                         ComputationContext& context) {
    auto evaluated = evaluate_numeric(expression, bindings, context);
    if (!evaluated) return evaluated;
    if (!evaluated.value().is_finite() ||
        !std::isfinite(evaluated.value().value)) {
        return Result<ApproxReal>::failure(
            CasErrc::NumericFailure,
            "LSR evalf produced a non-finite result",
            kEvalfOperation);
    }
    return evaluated;
}

Result<ApproxReal> evalf(const SymbolicExpr& expression,
                         const NumericBindings& bindings) {
    ComputationContext context;
    return evalf(expression, bindings, context);
}

Result<ApproxComplex> eval_complex(const SymbolicExpr& expression,
                                   const NumericBindings& bindings,
                                   ComputationContext& context) {
    try {
        if (!lamina::detail::node(expression)) {
            return complex_failure(CasErrc::InvalidArgument,
                                   "cannot evaluate an empty expression as complex",
                                   kEvalComplexOperation);
        }
        return evaluate_complex_node(lamina::detail::node(expression),
                                     bindings, context);
    } catch (const std::bad_alloc&) {
        return complex_failure(CasErrc::ResourceLimit,
                               "complex evaluation allocation failed",
                               kEvalComplexOperation);
    } catch (const std::exception& error) {
        return complex_failure(CasErrc::InternalInvariant, error.what(),
                               kEvalComplexOperation);
    }
}

Result<ApproxComplex> eval_complex(const SymbolicExpr& expression,
                                   const NumericBindings& bindings) {
    ComputationContext context;
    return eval_complex(expression, bindings, context);
}

} // namespace lamina::lsr
