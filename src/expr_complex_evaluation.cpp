#include "expr.hpp"

#include <cmath>
#include <exception>
#include <limits>
#include <utility>

#include "internal/expr_internal.hpp"
#include "internal/expr_common.hpp"
#include "internal/symbolic_ast.hpp"
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

ApproxReal approx_part(double value, double error = 0.0) {
    ApproxReal part;
    part.value = value;
    part.absolute_error = error;
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
    const double unknown = std::numeric_limits<double>::infinity();
    return Result<ApproxComplex>::success(
        ApproxComplex{approx_part(real, unknown), approx_part(imag, unknown)});
}

Result<ApproxComplex> real_to_complex(const Result<ApproxReal>& real) {
    if (!real) {
        return eval_complex_failure(real.error());
    }
    if (!real.value().is_finite() || !std::isfinite(real.value().value)) {
        return complex_failure(CasErrc::NumericFailure,
                               "complex evaluation requires finite real components",
                               kEvalComplexOperation);
    }
    return Result<ApproxComplex>::success(
        ApproxComplex{real.value(), approx_part(0.0)});
}

Result<ApproxComplex> add_complex(const ApproxComplex& lhs,
                                  const ApproxComplex& rhs) {
    return checked_complex(lhs.real.value + rhs.real.value,
                           lhs.imag.value + rhs.imag.value,
                           kEvalComplexOperation);
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
    return checked_complex(product.real, product.imag, kEvalComplexOperation);
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

Result<ApproxComplex> evaluate_complex_node(
    const std::shared_ptr<const SymbolicNode>& node,
    const NumericBindings& bindings,
    ComputationContext& context);

Result<ApproxComplex> evaluate_complex_components(
    const ComplexNode& node, const NumericBindings& bindings,
    ComputationContext& context) {
    auto real = evaluate_numeric(
        *detail::make_expression_ptr(node.real()), bindings, context);
    if (!real) {
        return eval_complex_failure(real.error());
    }
    auto imag = evaluate_numeric(
        *detail::make_expression_ptr(node.imag()), bindings, context);
    if (!imag) {
        return eval_complex_failure(imag.error());
    }
    if (!real.value().is_finite() || !imag.value().is_finite() ||
        !std::isfinite(real.value().value) ||
        !std::isfinite(imag.value().value)) {
        return complex_failure(CasErrc::NumericFailure,
                               "complex components must be finite",
                               kEvalComplexOperation);
    }
    return Result<ApproxComplex>::success(
        ApproxComplex{real.value(), imag.value()});
}

template <typename Combine>
Result<ApproxComplex> evaluate_complex_operands(
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands,
    ApproxComplex identity, Combine combine, const NumericBindings& bindings,
    ComputationContext& context) {
    auto result = Result<ApproxComplex>::success(identity);
    for (const auto& operand : operands) {
        auto value = evaluate_complex_node(operand, bindings, context);
        if (!value) {
            return value;
        }
        result = combine(result.value(), value.value());
        if (!result) {
            return result;
        }
    }
    return result;
}

Result<ApproxComplex> evaluate_complex_exponent(
    const std::shared_ptr<const SymbolicNode>& node,
    const NumericBindings& bindings, ComputationContext& context) {
    auto expression = simplify(detail::make_expression_ptr(node), context);
    if (!expression) {
        return complex_failure(expression.error().code,
                               std::move(expression.error().message),
                               kEvalComplexOperation);
    }
    const auto* number = dynamic_cast<const NumberNode*>(
        detail::node(expression.value()).get());
    const auto* rational = number ? std::get_if<Rational>(&number->value()) : nullptr;
    auto exponent = evaluate_complex_node(detail::node(expression.value()),
                                          bindings, context);
    if (!exponent) return exponent;
    if (rational && rational->get_denominator() != BigInt(1) &&
        std::trunc(exponent.value().real.value) == exponent.value().real.value) {
        return complex_failure(CasErrc::UnsupportedExpression,
                               "exact fractional exponent rounds to an integer",
                               kEvalComplexOperation);
    }
    return exponent;
}

Result<ApproxComplex> integer_complex_power(Result<ApproxComplex> base, int integer) {
    unsigned magnitude = static_cast<unsigned>(std::abs(integer));
    auto factor = base;
    if (integer < 0) {
        factor = divide_complex(approx_complex(1.0, 0.0), base.value());
        if (!factor) {
            return factor;
        }
    }
    auto result = Result<ApproxComplex>::success(approx_complex(1.0, 0.0));
    if (magnitude == 0U) {
        return checked_complex(1.0, 0.0, kEvalComplexOperation);
    }
    while (magnitude != 0U) {
        if ((magnitude & 1U) != 0U) {
            result = multiply_complex(result.value(), factor.value());
            if (!result) {
                return result;
            }
        }
        magnitude >>= 1U;
        if (magnitude != 0U) {
            factor = multiply_complex(factor.value(), factor.value());
            if (!factor) {
                return factor;
            }
        }
    }
    return result;
}

Result<ApproxComplex> evaluate_complex_power(
    const PowerNode& power, const NumericBindings& bindings,
    ComputationContext& context) {
    auto base = evaluate_complex_node(power.base(), bindings, context);
    if (!base) {
        return base;
    }
    auto exponent = evaluate_complex_exponent(power.exponent(), bindings, context);
    if (!exponent) {
        return Result<ApproxComplex>::failure(std::move(exponent.error()));
    }
    if (exponent.value().imag.value != 0.0 ||
        std::trunc(exponent.value().real.value) != exponent.value().real.value) {
        const lmmc_complex_t value{base.value().real.value, base.value().imag.value};
        const lmmc_complex_t power{exponent.value().real.value,
                                   exponent.value().imag.value};
        lmmc_complex_t result;
        const auto status = lmmc_complex_pow(&value, &power, &result);
        if (status == LMMC_STATUS_OUT_OF_RANGE) {
            return complex_failure(CasErrc::DomainError, "complex power of zero",
                                   kEvalComplexOperation);
        }
        if (status != LMMC_STATUS_OK) {
            return complex_failure(CasErrc::NumericFailure,
                                   "complex power produced a non-finite result",
                                   kEvalComplexOperation);
        }
        return checked_complex(result.real, result.imag, kEvalComplexOperation);
    }
    if (std::abs(exponent.value().real.value) > 64.0) {
        return complex_failure(CasErrc::UnsupportedExpression,
                               "complex evaluation only supports integer powers with |n| <= 64",
                               kEvalComplexOperation);
    }
    return integer_complex_power(std::move(base), static_cast<int>(exponent.value().real.value));
}

Result<ApproxComplex> evaluate_complex_node(
    const std::shared_ptr<const SymbolicNode>& node,
    const NumericBindings& bindings,
    ComputationContext& context) {
    auto entered = context.enter_recursion(kEvalComplexOperation);
    if (!entered) {
        return Result<ApproxComplex>::failure(entered.error());
    }
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
        return evaluate_complex_components(*complex_node, bindings, context);
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return evaluate_complex_operands(
            add->operands(), approx_complex(0.0, 0.0), add_complex, bindings, context);
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return evaluate_complex_operands(
            multiply->operands(), approx_complex(1.0, 0.0), multiply_complex,
            bindings, context);
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return evaluate_complex_power(*power, bindings, context);
    }

    return real_to_complex(
        evaluate_numeric(*LMCAS::detail::make_expression_ptr(node),
                         bindings, context));
}

}

Result<ApproxReal> evalf(const SymbolicExpr& expression,
                         const NumericBindings& bindings,
                         ComputationContext& context) {
    auto evaluated = evaluate_numeric(expression, bindings, context);
    if (!evaluated) { return evaluated; }
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

}
