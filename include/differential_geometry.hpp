/**
 * @file differential_geometry.hpp
 * @brief 微分几何：度量张量、克里斯托弗符号、黎曼曲率张量。
 */
#pragma once
#include "computation_context.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include <memory>
#include <vector>
#include <string>

namespace lamina {

using DifferentialGeometryExprResult = Result<std::shared_ptr<SymbolicExpr>>;
using DifferentialGeometryVectorResult = Result<std::vector<std::shared_ptr<SymbolicExpr>>>;

/**
 * @brief 计算度量张量的逆
 * @param g_ij 协变度量张量（矩阵表示）
 * @return 逆变度量张量 g^ij
 */
LAMINA_API DifferentialGeometryExprResult metric_inverse_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    ComputationContext& context
);

LAMINA_API DifferentialGeometryExprResult metric_inverse_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij
);

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
LAMINA_API DifferentialGeometryExprResult christoffel_first_kind_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j,
    ComputationContext& context
);

LAMINA_API DifferentialGeometryExprResult christoffel_first_kind_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j
);

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
LAMINA_API DifferentialGeometryExprResult christoffel_second_kind_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::shared_ptr<SymbolicExpr>& g_up_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j,
    ComputationContext& context
);

LAMINA_API DifferentialGeometryExprResult christoffel_second_kind_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::shared_ptr<SymbolicExpr>& g_up_ij,
    const std::vector<std::string>& coords,
    int k, int i, int j
);

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
LAMINA_API DifferentialGeometryExprResult riemann_curvature_tensor_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::vector<std::string>& coords,
    int rho, int sigma, int mu, int nu,
    ComputationContext& context
);

LAMINA_API DifferentialGeometryExprResult riemann_curvature_tensor_checked(
    const std::shared_ptr<SymbolicExpr>& g_ij,
    const std::vector<std::string>& coords,
    int rho, int sigma, int mu, int nu
);

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
LAMINA_API DifferentialGeometryExprResult lie_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& X,
    const std::vector<std::string>& vars,
    int order,
    ComputationContext& context
);

LAMINA_API DifferentialGeometryExprResult lie_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& X,
    const std::vector<std::string>& vars,
    int order = 1
);

LAMINA_API std::shared_ptr<SymbolicExpr> lie_derivative(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& X,
    const std::vector<std::string>& vars,
    int order = 1
);

/**
 * @brief 计算微分形式的外微分 d(form)。
 *
 * k-形式的系数按递增坐标指标组合的字典序排列。对于指标组合
 * J=(j0,...,jk)，输出系数为
 * sum_r (-1)^r * partial(omega[J without jr]) / partial(x[jr])。
 * 0-形式传单个系数；n 维空间中的 n-形式返回空的 (n+1)-形式系数表。
 *
 * @param form_coeffs 按上述顺序排列的形式系数（0-形式传单元素 [f]）
 * @param degree      形式次数，范围为 [0, vars.size()]
 * @param vars        坐标变量名列表
 * @return 外微分后的系数列表
 */
LAMINA_API DifferentialGeometryVectorResult exterior_derivative_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& form_coeffs,
    int degree,
    const std::vector<std::string>& vars,
    ComputationContext& context
);

LAMINA_API DifferentialGeometryVectorResult exterior_derivative_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& form_coeffs,
    int degree,
    const std::vector<std::string>& vars
);

LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> exterior_derivative(
    const std::vector<std::shared_ptr<SymbolicExpr>>& form_coeffs,
    int degree,
    const std::vector<std::string>& vars
);

} // namespace lamina
