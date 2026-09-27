#pragma once

#include "vector_calculus_types.hpp"

namespace LMCAS {

enum class CriticalPointClassification {
    LocalMinimum,
    LocalMaximum,
    Saddle,
    Degenerate,
    Inconclusive
};

/** @brief 临界点的坐标、分类及成立条件。 */
struct CriticalPoint {
    std::map<std::string, std::shared_ptr<SymbolicExpr>> point;
    CriticalPointClassification classification =
        CriticalPointClassification::Inconclusive;
    std::vector<std::shared_ptr<SymbolicExpr>> conditions;
};

using ExtremaResult = Result<std::vector<CriticalPoint>>;
using LagrangeResult =
    Result<std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>>;

/**
 * @brief 求多元函数的极值点并分类。
 *
 * 求解 ∇f = 0，并按海森矩阵的特征值分类：
 * - 所有特征值为正 → 局部极小值 (minimum)
 * - 所有特征值为负 → 局部极大值 (maximum)
 * - 特征值正负混合 → 鞍点 (saddle)
 * - 海森矩阵奇异 → 退化点 (degenerate)
 *
 * @param[in] f    标量函数表达式
 * @param[in] vars 变量名列表
 * @return 临界点列表，每个包含坐标和分类
 */
LMCAS_API ExtremaResult find_extrema_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context);

/** @brief 使用默认计算上下文求极值点；无效输入和未覆盖域以错误返回。 */
LMCAS_API ExtremaResult find_extrema_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars);

LMCAS_API std::vector<CriticalPoint> find_extrema(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars);

/**
 * @brief 使用拉格朗日乘数法求约束极值。
 *
 * 构造系统 ∇f = λ₁∇g₁ + λ₂∇g₂ + ... 并联合约束方程 gᵢ = 0 求解。
 * 支持单个和多个等式约束。
 *
 * @param[in] f           目标函数表达式
 * @param[in] constraints 约束列表（每个约束表达式等于零）
 * @param[in] vars        变量名列表
 * @return 所有临界点的列表，每个为变量名到值的映射（仅包含原始变量，不含乘数）
 */
LMCAS_API LagrangeResult lagrange_multipliers_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars,
    ComputationContext& context);

/** @brief 使用默认计算上下文求拉格朗日候选；无效输入和未覆盖域以错误返回。 */
LMCAS_API LagrangeResult lagrange_multipliers_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars);

LMCAS_API std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> lagrange_multipliers(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars);

}
