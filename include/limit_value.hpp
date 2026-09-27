#pragma once

#include "limit_direction.hpp"
#include "proof_outcome.hpp"
#include "result.hpp"
#include <memory>
#include <variant>

namespace LMCAS {

class SymbolicExpr;
using LimitExprPtr = std::shared_ptr<SymbolicExpr>;

/**
 * @brief 已证在请求的实数或复数域内有定义的有限值。
 * 表达式树中所有极限均已求定，扩展实数运算均有确定结果。
 */
struct FiniteLimit {
    LimitExprPtr value;
};
struct PositiveInfinityLimit {};
struct NegativeInfinityLimit {};
struct LimitDoesNotExist {};

using LimitOutcome = std::variant<
    FiniteLimit, PositiveInfinityLimit, NegativeInfinityLimit, LimitDoesNotExist>;
using LimitResult = Result<Verified<LimitOutcome>>;
using LimitExpressionResult = Result<LimitExprPtr>;

}
