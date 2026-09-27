#pragma once

#include "vector_calculus_types.hpp"

namespace LMCAS {

/**
 * @brief 计算标量函数的梯度 ∇f。
 *
 * @param[in] f    标量函数表达式
 * @param[in] vars 变量名列表
 * @return 梯度向量，各分量为偏导数 ∂f/∂xᵢ
 */
LMCAS_API VectorCalculusFieldResult gradient_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context);

/** @brief 使用默认上下文计算梯度；无效输入返回错误。 */
LMCAS_API VectorCalculusFieldResult gradient_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars);
LMCAS_API VectorField gradient(const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars);


/**
 * @brief 计算向量场的散度 ∇·F = ∑∂Fᵢ/∂xᵢ。
 *
 * @param[in] F    向量场（各分量为标量表达式）
 * @param[in] vars 变量名列表（与 F 的分量一一对应）
 * @return 散度标量表达式
 */
LMCAS_API VectorCalculusExprResult divergence_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    ComputationContext& context);

/** @brief 使用默认上下文计算散度；无效输入返回错误。 */
LMCAS_API VectorCalculusExprResult divergence_checked(
    const VectorField& F,
    const std::vector<std::string>& vars);
LMCAS_API std::shared_ptr<SymbolicExpr> divergence(const VectorField& F,
    const std::vector<std::string>& vars);


/**
 * @brief 计算向量场的旋度 ∇×F。
 *
 * 三维情况返回三分量向量；二维情况返回标量旋度（单分量向量）。
 *
 * @param[in] F    向量场
 * @param[in] vars 变量名列表
 * @return 旋度向量场
 */
LMCAS_API VectorCalculusFieldResult curl_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    ComputationContext& context);

/** @brief 使用默认上下文计算旋度；无效输入返回错误。 */
LMCAS_API VectorCalculusFieldResult curl_checked(
    const VectorField& F,
    const std::vector<std::string>& vars);
LMCAS_API VectorField curl(const VectorField& F,
    const std::vector<std::string>& vars);


/**
 * @brief 计算标量函数的拉普拉斯算子 ∇²f = ∑∂²f/∂xᵢ²。
 *
 * @param[in] f    标量函数表达式
 * @param[in] vars 变量名列表
 * @return 拉普拉斯算子结果表达式
 */
LMCAS_API VectorCalculusExprResult laplacian_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context);

/** @brief 使用默认上下文计算拉普拉斯算子；无效输入返回错误。 */
LMCAS_API VectorCalculusExprResult laplacian_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars);
LMCAS_API std::shared_ptr<SymbolicExpr> laplacian(const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars);


/**
 * @brief 计算方向导数 D_u f = ∇f · û。
 *
 * 将方向向量归一化为单位向量后，计算梯度与单位方向的点积。
 * 支持高阶方向导数（重复应用）。
 *
 * @param[in] f         标量函数表达式
 * @param[in] vars      变量名列表
 * @param[in] direction 方向向量
 * @param[in] order     阶数（默认 1）
 * @return 方向导数表达式；方向向量为零时返回 nullptr
 */
LMCAS_API VectorCalculusExprResult directional_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    const VectorField& direction,
    int order,
    ComputationContext& context);

/** @brief 使用默认上下文计算方向导数；零方向及其他无效输入返回错误。 */
LMCAS_API VectorCalculusExprResult directional_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    const VectorField& direction,
    int order = 1);
LMCAS_API std::shared_ptr<SymbolicExpr> directional_derivative(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars,
    const VectorField& direction, int order = 1);

}
