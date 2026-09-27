#pragma once
#include "computation_context.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include <memory>
#include <string>
#include <vector>

namespace LMCAS {

using VectorExprListResult = Result<std::vector<std::shared_ptr<SymbolicExpr>>>;
using VectorAngleResult = Result<double>;
using VectorStringResult = Result<std::string>;

/**
 * @brief 计算两个符号向量的点积。
 * @param a 向量 a 的各分量。
 * @param b 向量 b 的各分量。
 * @return 点积的符号表达式。
 */
LMCAS_API ExpressionResult vector_dot_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b,
    ComputationContext& context
);

LMCAS_API ExpressionResult vector_dot_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

LMCAS_API std::shared_ptr<SymbolicExpr> vector_dot(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

/**
 * @brief 计算两个三维符号向量的叉积。
 * @param a 向量 a 的各分量。
 * @param b 向量 b 的各分量。
 * @return 叉积向量的各分量。
 */
LMCAS_API VectorExprListResult vector_cross_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b,
    ComputationContext& context
);

LMCAS_API VectorExprListResult vector_cross_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

LMCAS_API std::vector<std::shared_ptr<SymbolicExpr>> vector_cross(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

/**
 * @brief 分别按向量尺度归一化，计算两个有限数值向量的夹角（弧度）。
 * 缩放范数与补偿点积控制中间溢出及相消误差。
 * 零向量返回 `CasErrc::DomainError`；非数值或非有限分量返回数值错误。
 * @param a 向量 a 的各分量。
 * @param b 向量 b 的各分量。
 * @return `[0,pi]` 内的有限夹角。
 */
LMCAS_API VectorAngleResult vector_angle_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b,
    ComputationContext& context
);

LMCAS_API VectorAngleResult vector_angle_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
);

}
