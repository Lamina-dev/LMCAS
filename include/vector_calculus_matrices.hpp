#pragma once

#include "vector_calculus_types.hpp"

namespace LMCAS {

/**
 * @brief 计算向量值函数的雅可比矩阵。
 *
 * 返回 m×n 的 MatrixNode，其中 entry(i,j) = ∂fᵢ/∂xⱼ。
 * 支持非方阵（m 个函数，n 个变量）。
 *
 * @param[in] functions 函数列表 (f₁, f₂, ..., fₘ)
 * @param[in] vars      变量名列表 (x₁, x₂, ..., xₙ)
 * @return 包含 MatrixNode 的 SymbolicExpr
 */
LMCAS_API VectorCalculusExprResult jacobian_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars,
    ComputationContext& context);

/** @brief 使用默认计算上下文计算雅可比矩阵，并报告无效输入。 */
LMCAS_API VectorCalculusExprResult jacobian_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars);
LMCAS_API std::shared_ptr<SymbolicExpr> jacobian(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars);


/**
 * @brief 计算标量函数的海森矩阵。
 *
 * 返回对称的 n×n MatrixNode，其中 entry(i,j) = ∂²f/(∂xᵢ∂xⱼ)。
 *
 * @param[in] f    标量函数表达式
 * @param[in] vars 变量名列表 (x₁, x₂, ..., xₙ)
 * @return 包含 MatrixNode 的 SymbolicExpr
 */
LMCAS_API VectorCalculusExprResult hessian_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context);

/** @brief 使用默认计算上下文计算海森矩阵，并报告无效输入。 */
LMCAS_API VectorCalculusExprResult hessian_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars);
LMCAS_API std::shared_ptr<SymbolicExpr> hessian(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars);

}
