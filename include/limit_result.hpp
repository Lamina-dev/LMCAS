#pragma once

#include "limit_value.hpp"
#include "computation_context.hpp"
#include "assumption.hpp"
#include "lmcas_export.hpp"

#include <string>


namespace LMCAS {

/**
 * @brief 证明单侧极限，并在有限点合并两侧结果。
 * 两侧已证极限不相容时返回 LimitDoesNotExist；证明不足时返回 Inconclusive。
 * 资源与取消错误原样保留。无穷点表示尾部趋近，仅有一种趋近方向。
 * Domain::Real 保留实分支与有定义性约束。
 * Domain::Complex 在支持的分支切割之外用连续性证明有限极限；
 * 已证去心邻域有定义时也支持多项式分母的可去奇点，要求方向为 Both、极限点有限，
 * 并采用复域邻域语义。
 * 成功认证的是上下文假设下的带标签结果，可能为无穷或极限不存在。
 * 有限单侧值按指定域逐点比较；残差非恒等于零仍可能在部分允许参数值处为零。
 * 恒等式认证相等前须证明有定义性；参数比较未知时返回 Inconclusive。
 */

LMCAS_API LimitResult limit_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction,
    ComputationContext& context,
    Domain domain = Domain::Real);

LMCAS_API LimitResult limit_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction = LimitDirection::Both,
    Domain domain = Domain::Real);

/**
 * @brief 将有限值与有符号无穷极限投影为表达式。
 * DoesNotExist 转为 Inconclusive，表达式载荷仅承载极限值。
 */
LMCAS_API LimitExpressionResult limit_expression_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction,
    ComputationContext& context,
    Domain domain = Domain::Real);

LMCAS_API LimitExpressionResult limit_expression_checked(
    const LimitExprPtr& expression,
    const std::string& variable,
    const LimitExprPtr& point,
    LimitDirection direction = LimitDirection::Both,
    Domain domain = Domain::Real);


} // namespace LMCAS
