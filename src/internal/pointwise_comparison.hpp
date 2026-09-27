#pragma once

#include "assumption.hpp"
#include "residual_verification.hpp"

namespace LMCAS::detail {

/**
 * @brief 在给定假设的 Real 或 Complex 域内比较逐点值，并保留错误。
 * True/False 分别证明值处处有定义且相等/不等；定义性未证实也归为 Unknown，
 * 仅有非恒等残差不足以判定 False。
 */
Result<Tribool> compare_pointwise_values(
    const ExprPtr& left, const ExprPtr& right,
    ComputationContext& context, Domain domain);

}
