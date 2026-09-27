#pragma once

#include "vector_calculus_types.hpp"

namespace LMCAS {

/**
 * @brief 格林定理：计算 ∬_R (∂Q/∂x - ∂P/∂y) dA。
 *
 * 等价于闭合曲线积分 ∮_C P dx + Q dy。
 * region 以积分区间表示：x ∈ [x_lo, x_hi], y ∈ [y_lo(x), y_hi(x)]。
 *
 * @param[in] P       向量场 x 分量
 * @param[in] Q       向量场 y 分量
 * @param[in] vars    变量名列表 {x_var, y_var}
 * @param[in] x_bounds x 的积分区间 {下界, 上界}
 * @param[in] y_bounds y 的积分区间 {下界(可含x), 上界(可含x)}
 * @return 格林定理计算的二重积分结果
 */
LMCAS_API VectorCalculusExprResult greens_theorem_checked(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算格林定理面积分，并显式报告无效输入和未覆盖域。
 */
LMCAS_API VectorCalculusExprResult greens_theorem_checked(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds);

LMCAS_API std::shared_ptr<SymbolicExpr> greens_theorem(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds);

/**
 * @brief 利用格林定理计算封闭曲线围成的面积 A = (1/2)∮(x dy - y dx)。
 *
 * @param[in] parametrization 曲线参数化 r(t) = (x(t), y(t))
 * @param[in] t               参数变量名
 * @param[in] a               参数下界
 * @param[in] b               参数上界
 * @return 面积表达式
 */
LMCAS_API VectorCalculusExprResult greens_theorem_area_checked(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算格林面积公式，并显式报告无效输入和未覆盖域。
 */
LMCAS_API VectorCalculusExprResult greens_theorem_area_checked(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b);

LMCAS_API std::shared_ptr<SymbolicExpr> greens_theorem_area(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b);

/**
 * @brief 散度定理（高斯定理）：计算 ∭_V ∇·F dV。
 *
 * 等价于封闭曲面上的通量积分 ∬_S F·dS。
 * volume_bounds 以迭代积分区间表示。
 *
 * @param[in] F             三维向量场
 * @param[in] vars          变量名列表 {x, y, z}
 * @param[in] x_bounds      x 的积分区间
 * @param[in] y_bounds      y 的积分区间（可含 x）
 * @param[in] z_bounds      z 的积分区间（可含 x, y）
 * @return 散度定理计算的三重积分结果
 */
LMCAS_API VectorCalculusExprResult divergence_theorem_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& z_bounds,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算散度定理体积分，并显式报告无效输入和未覆盖域。
 */
LMCAS_API VectorCalculusExprResult divergence_theorem_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& z_bounds);

LMCAS_API std::shared_ptr<SymbolicExpr> divergence_theorem(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& z_bounds);

/**
 * @brief 斯托克斯定理：计算 ∬_S (∇×F)·dS。
 *
 * 等价于边界曲线上的环量积分 ∮_C F·dr。
 * 通过参数化曲面计算 ∬(∇×F)·(r_u × r_v) du dv。
 *
 * @param[in] F               三维向量场
 * @param[in] vars            变量名列表 {x, y, z}
 * @param[in] parametrization 曲面参数化 r(u,v) = (x(u,v), y(u,v), z(u,v))
 * @param[in] u               第一参数变量名
 * @param[in] v               第二参数变量名
 * @param[in] u_bounds        u 的积分区间
 * @param[in] v_bounds        v 的积分区间
 * @return 斯托克斯定理计算的曲面积分结果
 */
LMCAS_API VectorCalculusExprResult stokes_theorem_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& v_bounds,
    ComputationContext& context);

/**
 * @brief 使用默认计算上下文计算斯托克斯定理曲面积分，并显式报告无效输入和未覆盖域。
 */
LMCAS_API VectorCalculusExprResult stokes_theorem_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& v_bounds);

LMCAS_API std::shared_ptr<SymbolicExpr> stokes_theorem(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& v_bounds);

}
