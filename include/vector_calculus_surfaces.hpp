#pragma once

#include "vector_calculus_types.hpp"

namespace LMCAS {

/**
 * @brief 计算第一类曲面积分（标量场在曲面上的积分）。
 *
 * 计算 ∬ f(r(u,v))·|r_u × r_v| du dv，其中 r(u,v) 为参数化曲面。
 *
 * @param[in] f              标量函数表达式（以坐标变量表示）
 * @param[in] parametrization 参数化曲面 [x(u,v), y(u,v), z(u,v)]
 * @param[in] u              第一参数变量名
 * @param[in] v              第二参数变量名
 * @param[in] u_lower        u 参数下界
 * @param[in] u_upper        u 参数上界
 * @param[in] v_lower        v 参数下界
 * @param[in] v_upper        v 参数上界
 * @return 曲面积分结果表达式
 */
LMCAS_API VectorCalculusExprResult surface_integral_scalar_checked(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context);

/** @brief 用默认上下文计算第一类曲面积分，报告无效输入与未覆盖域。 */
LMCAS_API VectorCalculusExprResult surface_integral_scalar_checked(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper);

LMCAS_API std::shared_ptr<SymbolicExpr> surface_integral_scalar(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper);

/**
 * @brief 计算第二类曲面积分（向量场通过曲面的通量）。
 *
 * 计算 ∬ F·(r_u × r_v) du dv，其中 r(u,v) 为参数化曲面。
 *
 * @param[in] F              向量场（三分量，以坐标变量表示）
 * @param[in] parametrization 参数化曲面 [x(u,v), y(u,v), z(u,v)]
 * @param[in] u              第一参数变量名
 * @param[in] v              第二参数变量名
 * @param[in] u_lower        u 参数下界
 * @param[in] u_upper        u 参数上界
 * @param[in] v_lower        v 参数下界
 * @param[in] v_upper        v 参数上界
 * @return 曲面积分结果表达式
 */
LMCAS_API VectorCalculusExprResult surface_integral_vector_checked(
    const VectorField& F, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context);

/** @brief 用默认上下文计算第二类曲面积分，报告无效输入与未覆盖域。 */
LMCAS_API VectorCalculusExprResult surface_integral_vector_checked(
    const VectorField& F, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper);

LMCAS_API std::shared_ptr<SymbolicExpr> surface_integral_vector(
    const VectorField& F, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper);

}
