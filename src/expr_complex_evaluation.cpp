#include "expr.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <utility>

#include "expr_internal.hpp"
#include "symbolic_ast.hpp"
#include "lmmc/complex.h"

namespace LMCAS {
namespace {

constexpr const char* kEvalfOperation = "LMCAS.evalf";

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
    const lmmc_complex_t left{a, b};
    const lmmc_complex_t right{c, d};
    lmmc_complex_t product;
    if (lmmc_complex_mul(&left, &right, &product) != LMMC_STATUS_OK) {
        return complex_failure(CasErrc::NumericFailure,
                               "complex multiplication requires finite components and result",
                               kEvalComplexOperation);
    }
    auto result = checked_complex(product.real, product.imag, kEvalComplexOperation);
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
    const double a = lhs.real.value;
    const double b = lhs.imag.value;
    const double c = rhs.real.value;
    const double d = rhs.imag.value;
    if (c == 0.0 && d == 0.0) {
        return complex_failure(CasErrc::DomainError,
                               "complex division by zero",
                               kEvalComplexOperation);
    }
    const lmmc_complex_t left{a, b};
    const lmmc_complex_t right{c, d};
    lmmc_complex_t quotient;
    if (lmmc_complex_div(&left, &right, &quotient) != LMMC_STATUS_OK) {
        return complex_failure(CasErrc::NumericFailure,
                               "complex division requires finite components and result",
                               kEvalComplexOperation);
    }
    return checked_complex(quotient.real, quotient.imag, kEvalComplexOperation);
}

bool is_integer_double(double value) {
    return std::isfinite(value) && std::floor(value) == value;
}

Result<ApproxComplex> evaluate_complex_node(
    const std::shared_ptr<const SymbolicNode>& node,
    const NumericBindings& bindings,
    ComputationContext& context) {
    auto entered = context.enter_recursion(kEvalComplexOperation);
    if (!entered) return Result<ApproxComplex>::failure(entered.error());
    struct RecursionExit {
        ComputationContext& context;
        ~RecursionExit() { context.leave_recursion(); }
    } recursion_exit{context};
    if (!node) {
        return complex_failure(CasErrc::InvalidArgument,
                               "expression contains a null node",
                               kEvalComplexOperation);
    }

    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real = evaluate_numeric(
            *LMCAS::detail::make_expression_ptr(complex_node->real()),
            bindings, context);
        if (!real) return eval_complex_failure(real.error());
        auto imag = evaluate_numeric(
            *LMCAS::detail::make_expression_ptr(complex_node->imag()),
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
        auto exponent_expression = simplify(
            LMCAS::detail::make_expression_ptr(power->exponent()), context);
        if (!exponent_expression) {
            return eval_complex_failure(exponent_expression.error());
        }
        /**
         * Exact rational exponents retain their integer-domain classification
         * through normalization and conversion to the numeric power kernel.
         * @see David Goldberg, "What Every Computer Scientist Should Know About
         * Floating-Point Arithmetic" (1991), Floating-point Formats.
         * https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
         */
        if (auto number = std::dynamic_pointer_cast<const NumberNode>(
                LMCAS::detail::node(exponent_expression.value()))) {
            if (const auto* rational = std::get_if<Rational>(&number->value());
                rational && rational->get_denominator() != BigInt(1)) {
                return complex_failure(
                    CasErrc::UnsupportedExpression,
                    "complex evaluation requires an integer exponent",
                    kEvalComplexOperation);
            }
        }
        auto exponent = evaluate_numeric(
            *exponent_expression.value(), bindings, context);
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
        const int exponent_int = static_cast<int>(exponent_value);
        unsigned exponent_magnitude =
            static_cast<unsigned>(std::abs(exponent_int));
        auto factor = base;
        if (exponent_int < 0) {
            factor = divide_complex(approx_complex(1.0, 0.0), base.value());
            if (!factor) return factor;
        }
        auto result =
            Result<ApproxComplex>::success(approx_complex(1.0, 0.0));
        while (exponent_magnitude != 0U) {
            if ((exponent_magnitude & 1U) != 0U) {
                result = multiply_complex(result.value(), factor.value());
                if (!result) return result;
            }
            exponent_magnitude >>= 1U;
            if (exponent_magnitude != 0U) {
                factor = multiply_complex(factor.value(), factor.value());
                if (!factor) return factor;
            }
        }
        return result;
    }

    return real_to_complex(
        evaluate_numeric(*LMCAS::detail::make_expression_ptr(node),
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
            "LMCAS evalf produced a non-finite result",
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
        if (!LMCAS::detail::node(expression)) {
            return complex_failure(CasErrc::InvalidArgument,
                                   "cannot evaluate an empty expression as complex",
                                   kEvalComplexOperation);
        }
        return evaluate_complex_node(LMCAS::detail::node(expression),
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

} // namespace LMCAS
