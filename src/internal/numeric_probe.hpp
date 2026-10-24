#pragma once

#include "numeric_evaluation.hpp"

#include <cmath>
#include <optional>

namespace LMCAS::detail {

inline std::optional<double> try_finite_numeric(
    const SymbolicExpr& expression,
    ComputationContext& context)
{
    auto evaluated = evaluate_numeric(expression, NumericBindings{}, context);
    if (!evaluated) {
        switch (evaluated.error().code) {
        case CasErrc::UnboundSymbol:
        case CasErrc::UnsupportedExpression:
        case CasErrc::DomainError:
            return std::nullopt;
        default:
            throw evaluated.error();
        }
    }
    if (!evaluated.value().is_finite()) {
        return std::nullopt;
    }
    if (!std::isfinite(evaluated.value().value)) {
        throw CasError{CasErrc::NumericFailure,
                       "numeric evaluation produced a nonfinite result",
                       "evaluate_numeric"};
    }
    return evaluated.value().value;
}


} // namespace LMCAS::detail
