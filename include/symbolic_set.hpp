#pragma once

#include <memory>
#include <vector>

#include "computation_context.hpp"
#include "symbolic.hpp"

namespace LMCAS {

using SymbolicSetResult = Result<std::shared_ptr<SymbolicExpr>>;

/** Constructs an immutable finite set with structural deduplication. */
LMCAS_API SymbolicSetResult make_finite_set(
    std::vector<std::shared_ptr<SymbolicExpr>> elements,
    ComputationContext& context);

/** Constructs a real interval. Endpoint ordering is validated when decidable. */
LMCAS_API SymbolicSetResult make_interval(
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    bool lower_closed,
    bool upper_closed,
    ComputationContext& context);

/** Constructs and simplifies a membership proposition when membership is decidable. */
LMCAS_API SymbolicSetResult make_membership(
    const std::shared_ptr<SymbolicExpr>& element,
    const std::shared_ptr<SymbolicExpr>& set,
    ComputationContext& context);

LMCAS_API SymbolicSetResult finite_set_union(
    const SymbolicExpr& lhs, const SymbolicExpr& rhs,
    ComputationContext& context);
LMCAS_API SymbolicSetResult finite_set_intersection(
    const SymbolicExpr& lhs, const SymbolicExpr& rhs,
    ComputationContext& context);
LMCAS_API SymbolicSetResult finite_set_difference(
    const SymbolicExpr& lhs, const SymbolicExpr& rhs,
    ComputationContext& context);
LMCAS_API SymbolicSetResult finite_set_symmetric_difference(
    const SymbolicExpr& lhs, const SymbolicExpr& rhs,
    ComputationContext& context);

} // namespace LMCAS
