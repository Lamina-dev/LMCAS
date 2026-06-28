/**
 * @file differential_geometry.hpp
 * @brief 微分几何：度量张量、克里斯托弗符号、黎曼曲率张量。
 */
#pragma once
#include "symbolic_ast.hpp"
#include <memory>
#include <vector>
#include <string>

class SymbolicExpr;


namespace lamina {

/**
 * @brief 计算度量张量的逆
 * @param g_ij 协变度量张量（矩阵表示）
 * @return 逆变度量张量 g^ij
 */
LAMINA_API std::shared_ptr<SymbolicExpr> metric_inverse(
    const std::shared_ptr<SymbolicExpr>& g_ij
);

/**
 * @brief 计算第一类克里斯托弗符号 Γ_{kij}
 * @param g_ij 度量张量
 * @param coords 坐标变量列表
 * @param k 索引
 * @param i 索引
 * @param j 索引
 * @return 符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> christoffel_first_kind(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j
);

/**
 * @brief 计算第二类克里斯托弗符号 Γ^k_{ij}
 * @param g_ij 度量张量
 * @param g_up_ij 逆度量张量
 * @param coords 坐标变量列表
 * @param k 上标索引
 * @param i 下标索引
 * @param j 下标索引
 * @return 符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> christoffel_second_kind(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::shared_ptr<SymbolicExpr>& g_up_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j
);

/**
 * @brief 计算黎曼曲率张量 R^rho_{sigma mu nu}
 * @param g_ij 度量张量
 * @param coords 坐标变量
 * @param rho 上标
 * @param sigma 下标1
 * @param mu 下标2
 * @param nu 下标3
 * @return 符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> riemann_curvature_tensor(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::vector<std::string>& coords,
    int rho, int sigma, int mu, int nu
);

/**
 * @brief 计算标量函数沿向量场 X 的李导数 L_X f = ∑ Xⁱ·∂f/∂xⁱ。
 *
 * 支持迭代李导数（order 次重复应用）。
 *
 * @param f      标量函数
 * @param X      向量场分量（与 vars 一一对应）
 * @param vars   坐标变量名列表
 * @param order  迭代阶数（默认 1）
 * @return 李导数表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> lie_derivative(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& X,
    const std::vector<std::string>& vars,
    int order = 1
);

/**
 * @brief 计算微分形式的外微分 d(form)。
 *
 * - 0-形式（标量 f）：返回梯度分量 [∂f/∂x₁, ..., ∂f/∂xₙ]（即 1-形式系数）。
 * - 1-形式（系数 [ω₁,...,ωₙ]）：返回 2-形式分量，按 (i<j) 顺序的
 *   (∂ωⱼ/∂xᵢ - ∂ωᵢ/∂xⱼ)。
 *
 * @param form_coeffs 形式的系数（0-形式传单元素 [f]）
 * @param degree      形式的次数（0 或 1）
 * @param vars        坐标变量名列表
 * @return 外微分后的系数列表
 */
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> exterior_derivative(
    const std::vector<std::shared_ptr<SymbolicExpr>>& form_coeffs,
    int degree,
    const std::vector<std::string>& vars
);

} // namespace lamina
