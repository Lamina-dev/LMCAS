#pragma once

#include "vector_calculus_types.hpp"

namespace LMCAS {

/**
 * @brief 计算第一类曲线积分（标量场沿曲线的积分）。
 *
 * 计算 ∫ₐᵇ f(r(t))·|r'(t)| dt，其中 r(t) 为参数化曲线。
 *
 * @param[in] f              标量函数表达式（以曲线参数化中的坐标变量表示）
 * @param[in] parametrization 参数化曲线 [x(t), y(t)] 或 [x(t), y(t), z(t)]
 * @param[in] t              参数变量名
 * @param[in] a              参数下界
 * @param[in] b              参数上界
 * @return 曲线积分结果表达式
 */
LMCAS_API VectorCalculusExprResult curve_integral_scalar_checked(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context);

/** @brief 使用默认上下文计算第一类曲线积分，报告无效输入与未覆盖域。 */
LMCAS_API VectorCalculusExprResult curve_integral_scalar_checked(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b);

LMCAS_API std::shared_ptr<SymbolicExpr> curve_integral_scalar(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b);

/**
 * @brief 计算第二类曲线积分（向量场沿曲线的积分）。
 *
 * 计算 ∫ₐᵇ F(r(t))·r'(t) dt，其中 r(t) 为参数化曲线。
 *
 * @param[in] F              向量场（各分量以坐标变量表示）
 * @param[in] parametrization 参数化曲线 [x(t), y(t)] 或 [x(t), y(t), z(t)]
 * @param[in] t              参数变量名
 * @param[in] a              参数下界
 * @param[in] b              参数上界
 * @return 曲线积分结果表达式
 */
LMCAS_API VectorCalculusExprResult curve_integral_vector_checked(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context);

/** @brief 使用默认上下文计算第二类曲线积分，报告无效输入与未覆盖域。 */
LMCAS_API VectorCalculusExprResult curve_integral_vector_checked(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b);

LMCAS_API std::shared_ptr<SymbolicExpr> curve_integral_vector(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b);

}
