#pragma once

#include "result.hpp"
#include "computation_context.hpp"
#include "lamina_export.hpp"
#include "limit_direction.hpp"
#include "proof_outcome.hpp"

#include <string>

#include <memory>
#include <variant>

class SymbolicExpr;

namespace lamina {

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

LAMINA_API LimitResult limit_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction,
    ComputationContext& context);

LAMINA_API LimitResult limit_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction = LimitDirection::Both);

LAMINA_API LimitExpressionResult limit_expression_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction,
    ComputationContext& context);

LAMINA_API LimitExpressionResult limit_expression_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction = LimitDirection::Both);


} // namespace lamina
