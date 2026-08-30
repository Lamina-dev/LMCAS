#pragma once

#include <cmath>
#include <string>
#include <unordered_map>

#include "computation_context.hpp"
#include "lamina_export.hpp"

class SymbolicExpr;

namespace lamina {

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

LAMINA_API Result<ApproxReal> evaluate_numeric(const SymbolicExpr& expression,
                                               const NumericBindings& bindings,
                                               ComputationContext& context);

inline Result<ApproxReal> evaluate_numeric(const SymbolicExpr& expression,
                                           const NumericBindings& bindings = {}) {
    ComputationContext context;
    return evaluate_numeric(expression, bindings, context);
}

} // namespace lamina
