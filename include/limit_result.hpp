#pragma once

#include "result.hpp"
#include "computation_context.hpp"
#include "lmcas_export.hpp"
#include "limit_direction.hpp"
#include "proof_outcome.hpp"

#include <string>

#include <memory>
#include <variant>

namespace LMCAS {

class SymbolicExpr;



using LimitExprPtr = std::shared_ptr<SymbolicExpr>;

struct FiniteLimit {
    LimitExprPtr value;
};
struct PositiveInfinityLimit {};
struct NegativeInfinityLimit {};
struct LimitDoesNotExist {};

using LimitOutcome = std::variant<
    FiniteLimit,
    PositiveInfinityLimit,
    NegativeInfinityLimit,
    LimitDoesNotExist>;
using LimitResult = Result<Verified<LimitOutcome>>;
using LimitExpressionResult = Result<LimitExprPtr>;

LMCAS_API LimitResult limit_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction,
    ComputationContext& context);

LMCAS_API LimitResult limit_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction = LimitDirection::Both);

LMCAS_API LimitExpressionResult limit_expression_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction,
    ComputationContext& context);

LMCAS_API LimitExpressionResult limit_expression_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction = LimitDirection::Both);


} // namespace LMCAS
