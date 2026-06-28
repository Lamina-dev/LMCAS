/**
 * @file vector_calculus.hpp
 * @brief 向量微积分模块：梯度、散度、旋度、拉普拉斯算子、方向导数、雅可比矩阵、海森矩阵、极值分析。
 */
#pragma once

#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <utility>

class SymbolicExpr;


namespace lamina {

/// 向量场类型：标量表达式的向量
using VectorField = std::vector<std::shared_ptr<SymbolicExpr>>;

// ============================================================
// 微分算子 (Requirements 8, 9, 45, 46, 47, 87, 88)
// ============================================================

/**
 * @brief 计算标量函数的梯度 ∇f。
 *
 * 返回由各偏导数 ∂f/∂xᵢ 组成的向量。
 *
 * @param[in] f    标量函数表达式
 * @param[in] vars 变量名列表
 * @return 梯度向量，各分量为对应偏导数
 */
LAMINA_API VectorField gradient(const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars);

/**
 * @brief 计算向量场的散度 ∇·F = ∑∂Fᵢ/∂xᵢ。
 *
 * @param[in] F    向量场（各分量为标量表达式）
 * @param[in] vars 变量名列表（与 F 的分量一一对应）
 * @return 散度标量表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> divergence(const VectorField& F,
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
LAMINA_API VectorField curl(const VectorField& F,
    const std::vector<std::string>& vars);

/**
 * @brief 计算标量函数的拉普拉斯算子 ∇²f = ∑∂²f/∂xᵢ²。
 *
 * @param[in] f    标量函数表达式
 * @param[in] vars 变量名列表
 * @return 拉普拉斯算子结果表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> laplacian(const std::shared_ptr<SymbolicExpr>& f,
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
LAMINA_API std::shared_ptr<SymbolicExpr> directional_derivative(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars,
    const VectorField& direction, int order = 1);

// ============================================================
// 雅可比矩阵与海森矩阵 (Requirements 10, 11)
// ============================================================

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
LAMINA_API std::shared_ptr<SymbolicExpr> jacobian(
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
LAMINA_API std::shared_ptr<SymbolicExpr> hessian(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars);

// ============================================================
// 曲线积分与曲面积分 (Requirements 48, 49)
// ============================================================

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
LAMINA_API std::shared_ptr<SymbolicExpr> curve_integral_scalar(
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
LAMINA_API std::shared_ptr<SymbolicExpr> curve_integral_vector(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b);

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
LAMINA_API std::shared_ptr<SymbolicExpr> surface_integral_scalar(
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
LAMINA_API std::shared_ptr<SymbolicExpr> surface_integral_vector(
    const VectorField& F, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper);

// ============================================================
// 积分定理 (Requirements 50, 89, 92, 93)
// ============================================================

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
LAMINA_API std::shared_ptr<SymbolicExpr> greens_theorem(
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
LAMINA_API std::shared_ptr<SymbolicExpr> greens_theorem_area(
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
LAMINA_API std::shared_ptr<SymbolicExpr> divergence_theorem(
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
LAMINA_API std::shared_ptr<SymbolicExpr> stokes_theorem(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& v_bounds);

// ============================================================
// 多元极值与拉格朗日乘数法 (Requirements 44, 51, 91)
// ============================================================

/**
 * @brief 临界点信息，包含坐标和分类。
 */
struct CriticalPoint {
    std::map<std::string, std::shared_ptr<SymbolicExpr>> point; ///< 临界点坐标
    std::string classification; ///< 分类: "minimum", "maximum", "saddle", "degenerate"
};

/**
 * @brief 求多元函数的极值点并分类。
 *
 * 通过求解 ∇f = 0 找到所有临界点，然后计算海森矩阵并通过特征值分析进行分类：
 * - 所有特征值为正 → 局部极小值 (minimum)
 * - 所有特征值为负 → 局部极大值 (maximum)
 * - 特征值正负混合 → 鞍点 (saddle)
 * - 海森矩阵奇异 → 退化点 (degenerate)
 *
 * @param[in] f    标量函数表达式
 * @param[in] vars 变量名列表
 * @return 临界点列表，每个包含坐标和分类
 */
LAMINA_API std::vector<CriticalPoint> find_extrema(
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
LAMINA_API std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> lagrange_multipliers(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars);

// ============================================================
// 向量代数运算 (Requirements 40, 41, 42)
// ============================================================

/**
 * @brief 计算两个向量的点积 a·b = ∑aᵢbᵢ。
 * @param[in] a 向量 a
 * @param[in] b 向量 b（维度需与 a 相同）
 * @return 点积标量表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> dot_product(
    const VectorField& a, const VectorField& b);

/**
 * @brief 计算两个三维向量的叉积 a×b。
 * @param[in] a 三维向量 a
 * @param[in] b 三维向量 b
 * @return 叉积向量（三分量）
 */
LAMINA_API VectorField cross_product(
    const VectorField& a, const VectorField& b);

/**
 * @brief 计算向量 a 在向量 b 上的投影向量 proj_b(a) = (a·b / b·b)·b。
 * @param[in] a 被投影向量
 * @param[in] b 投影方向向量
 * @return 投影向量；b 为零向量时返回零向量
 */
LAMINA_API VectorField vector_project(
    const VectorField& a, const VectorField& b);

/**
 * @brief 计算向量 a 在向量 b 上的标量投影 a·b / |b|。
 * @param[in] a 被投影向量
 * @param[in] b 投影方向向量
 * @return 标量投影表达式；b 为零向量时返回 nullptr
 */
LAMINA_API std::shared_ptr<SymbolicExpr> scalar_project(
    const VectorField& a, const VectorField& b);

/**
 * @brief 计算两个向量夹角 arccos(a·b / (|a|·|b|))。
 * @param[in] a 向量 a
 * @param[in] b 向量 b
 * @return 夹角表达式（弧度）；任一向量为零时返回 nullptr
 */
LAMINA_API std::shared_ptr<SymbolicExpr> vector_angle_symbolic(
    const VectorField& a, const VectorField& b);

/**
 * @brief 计算三个三维向量的混合积 a·(b×c)。
 * @param[in] a 向量 a
 * @param[in] b 向量 b
 * @param[in] c 向量 c
 * @return 混合积标量表达式（等于以 a,b,c 为行的行列式）
 */
LAMINA_API std::shared_ptr<SymbolicExpr> mixed_product(
    const VectorField& a, const VectorField& b, const VectorField& c);

} // namespace lamina
