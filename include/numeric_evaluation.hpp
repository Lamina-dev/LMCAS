#pragma once

#include <cmath>
#include <string>
#include <unordered_map>

#include "computation_context.hpp"
#include "lmcas_export.hpp"

namespace LMCAS {

class SymbolicExpr;



enum class NumericStatus { Finite, PositiveInfinity, NegativeInfinity };

struct ApproxReal {
    double value = 0.0;
    double absolute_error = 0.0;
    NumericStatus status = NumericStatus::Finite;

    bool is_finite() const noexcept { return status == NumericStatus::Finite; }
};

struct ApproxComplex {
    ApproxReal real;
    ApproxReal imag;

    bool is_finite() const noexcept {
        return real.is_finite() && imag.is_finite();
    }
};

using NumericBindings = std::unordered_map<std::string, double>;

/**
 * @brief Evaluate a real expression using explicit symbol bindings and a computation context.
 *
 * Square roots of one or two real squares are evaluated as absolute values
 * or hypot, avoiding intermediate overflow and underflow in norm expressions.
 *
 * @see ISO C11 committee draft N1570, 7.12.7.3 (hypot).
 * https://www.open-std.org/jtc1/sc22/wg14/www/docs/n1570.pdf
 *
 * Exact numeric exponents retain their integer classification and parity on
 * negative real bases. Non-integral powers of negative real values report
 * DomainError rather than rounding the exponent into the integer domain.
 *
 * @see ISO C11 committee draft N1570, 7.12.7.4 (pow).
 * https://www.open-std.org/jtc1/sc22/wg14/www/docs/n1570.pdf
 *
 * Complex phase on a real argument embeds the value as x + i(+0) and evaluates
 * its principal argument with atan2(+0, x). Negative reals, including negative
 * zero and negative-infinity bindings, therefore have phase +pi; positive reals
 * and positive zero have phase +0. Argument-evaluation failures propagate through Result.
 *
 * @see ISO C11 committee draft N1570, 7.3.9.1 (carg), principal argument.
 * https://www.open-std.org/jtc1/sc22/wg14/www/docs/n1570.pdf
 */
LMCAS_API Result<ApproxReal> evaluate_numeric(const SymbolicExpr& expression,
                                               const NumericBindings& bindings,
                                               ComputationContext& context);

inline Result<ApproxReal> evaluate_numeric(const SymbolicExpr& expression,
                                           const NumericBindings& bindings = {}) {
    ComputationContext context;
    return evaluate_numeric(expression, bindings, context);
}

} // namespace LMCAS
