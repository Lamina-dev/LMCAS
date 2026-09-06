#pragma once

#include "numeric_evaluation.hpp"

#include <cmath>
#include <memory>
#include <optional>

namespace LMCAS::detail {

inline std::optional<double> try_finite_numeric(
    const SymbolicExpr& expression,
    ComputationContext* context = nullptr)
{
    ComputationContext local_context;
    ComputationContext& evaluation_context = context ? *context : local_context;
    auto evaluated =
        evaluate_numeric(expression, NumericBindings{}, evaluation_context);
    if (!evaluated || !evaluated.value().is_finite() ||
        !std::isfinite(evaluated.value().value)) {
        return std::nullopt;
    }
    return evaluated.value().value;
}

inline std::optional<double> try_finite_numeric(
    const std::shared_ptr<SymbolicExpr>& expression,
    ComputationContext* context = nullptr)
{
    return expression ? try_finite_numeric(*expression, context) : std::nullopt;
}

} // namespace LMCAS::detail
