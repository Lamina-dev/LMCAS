/**
 * @file matrix_decomposition.hpp
 * @brief 矩阵高级分解算法：LU、QR、Cholesky、SVD、特征分解及矩阵指数。
 */
#pragma once
#include "computation_context.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include <memory>
#include <vector>
#include <string>

namespace lamina {

struct LUDecomposition {
    std::shared_ptr<SymbolicExpr> L;
    std::shared_ptr<SymbolicExpr> U;
};

struct QRDecomposition {
    std::shared_ptr<SymbolicExpr> Q;
    std::shared_ptr<SymbolicExpr> R;
};

struct CholeskyDecomposition {
    std::shared_ptr<SymbolicExpr> L;
};

struct SVDDecomposition {
    std::shared_ptr<SymbolicExpr> U;
    std::shared_ptr<SymbolicExpr> S;
    std::shared_ptr<SymbolicExpr> V;
};

struct JordanDecomposition {
    std::shared_ptr<SymbolicExpr> J;
    std::shared_ptr<SymbolicExpr> P;
};

using LUDecompositionResult = Result<LUDecomposition>;
using QRDecompositionResult = Result<QRDecomposition>;
using CholeskyDecompositionResult = Result<CholeskyDecomposition>;
using SVDDecompositionResult = Result<SVDDecomposition>;
using JordanDecompositionResult = Result<JordanDecomposition>;

/**
 * @brief 计算符号矩阵的 LU 分解 (A = L * U)
 * @param A 输入方阵
 * @param L 输出的下三角矩阵
 * @param U 输出的上三角矩阵
 * @return 成功返回 true
 */
LAMINA_API LUDecompositionResult lu_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context
);

LAMINA_API LUDecompositionResult lu_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A
);

LAMINA_API bool lu_decomposition(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& L,
    std::shared_ptr<SymbolicExpr>& U
);

/**
 * @brief 计算符号矩阵的 QR 分解 (A = Q * R)
 * @param A 输入矩阵
 * @param Q 输出的正交矩阵
 * @param R 输出的上三角矩阵
 * @return 成功返回 true
 */
LAMINA_API QRDecompositionResult qr_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context
);

LAMINA_API QRDecompositionResult qr_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A
);

LAMINA_API bool qr_decomposition(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& Q,
    std::shared_ptr<SymbolicExpr>& R
);

/**
 * @brief 计算对称正定矩阵的 Cholesky 分解 (A = L * L^T)
 * @param A 输入矩阵
 * @param L 输出的下三角矩阵
 * @return 成功返回 true
 */
LAMINA_API CholeskyDecompositionResult cholesky_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context
);

LAMINA_API CholeskyDecompositionResult cholesky_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A
);

LAMINA_API bool cholesky_decomposition(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& L
);

/**
 * @brief 计算符号矩阵的奇异值分解 (A = U * S * V^T)
 * @param A 输入矩阵
 * @param U 输出的左奇异矩阵
 * @param S 输出的奇异值对角矩阵
 * @param V 输出的右奇异矩阵
 * @return 成功返回 true
 */
LAMINA_API SVDDecompositionResult svd_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context
);

LAMINA_API SVDDecompositionResult svd_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A
);

LAMINA_API bool svd_decomposition(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& U,
    std::shared_ptr<SymbolicExpr>& S,
    std::shared_ptr<SymbolicExpr>& V
);

/**
 * @brief 计算矩阵指数 e^A
 * @param A 输入矩阵
 * @return 指数矩阵的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> matrix_exp(
    const std::shared_ptr<SymbolicExpr>& A
);

/**
 * @brief 计算方阵的迹（对角元之和）。
 * @param A 输入方阵
 * @return 迹的符号表达式；非方阵返回 nullptr
 */
LAMINA_API std::shared_ptr<SymbolicExpr> matrix_trace(
    const std::shared_ptr<SymbolicExpr>& A
);

/**
 * @brief 计算矩阵的秩（通过 RREF 主元计数）。
 * @param A 输入矩阵
 * @return 秩（整数）；输入非矩阵返回 -1
 */
LAMINA_API int matrix_rank(
    const std::shared_ptr<SymbolicExpr>& A
);

/**
 * @brief 对一组向量执行 Gram-Schmidt 正交化。
 * @param vectors 输入向量列表（每个向量为分量列表）
 * @param normalize 为 true 时返回标准正交基，否则仅正交
 * @return 正交（或标准正交）向量列表；线性相关的向量被跳过
 */
LAMINA_API std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> gram_schmidt(
    const std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& vectors,
    bool normalize = false
);

/**
 * @brief 计算矩阵自然对数 log(A)（通过特征分解）。
 * @param A 输入方阵
 * @return log(A) 的符号表达式；不可对角化或含非正特征值时返回 nullptr
 */
LAMINA_API std::shared_ptr<SymbolicExpr> matrix_log(
    const std::shared_ptr<SymbolicExpr>& A
);

/**
 * @brief 计算两个矩阵的 Kronecker 积 A⊗B。
 * @param A m×n 矩阵
 * @param B p×q 矩阵
 * @return (mp)×(nq) 矩阵
 */
LAMINA_API std::shared_ptr<SymbolicExpr> kronecker(
    const std::shared_ptr<SymbolicExpr>& A,
    const std::shared_ptr<SymbolicExpr>& B
);

/**
 * @brief 计算矩阵范数。
 * @param A 输入矩阵
 * @param type 范数类型："frobenius"、"1"（最大列和）、"inf"（最大行和）
 * @return 范数的符号表达式；类型无效返回 nullptr
 */
LAMINA_API std::shared_ptr<SymbolicExpr> matrix_norm(
    const std::shared_ptr<SymbolicExpr>& A,
    const std::string& type = "frobenius"
);

/**
 * @brief 从二次型表达式提取对称矩阵 A，使得 expr = xᵀAx。
 * @param expr 二次齐次表达式
 * @param vars 变量名列表
 * @return 对称矩阵；提取失败返回 nullptr
 */
LAMINA_API std::shared_ptr<SymbolicExpr> quadratic_form_matrix(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::vector<std::string>& vars
);

/**
 * @brief 对二次型对称矩阵进行定性分类。
 * @param A 对称矩阵
 * @return 分类字符串："positive_definite"、"negative_definite"、
 *         "positive_semidefinite"、"negative_semidefinite"、"indefinite"、"unknown"
 */
LAMINA_API std::string classify_quadratic_form(
    const std::shared_ptr<SymbolicExpr>& A
);

/**
 * @brief 计算方阵的 Jordan 标准型 J 及变换矩阵 P（A = P·J·P⁻¹）。
 * @param A 输入方阵
 * @param J 输出 Jordan 标准型
 * @param P 输出变换矩阵
 * @return 成功返回 true（可对角化或可求 Jordan 块时）
 */
LAMINA_API JordanDecompositionResult jordan_form_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context
);

LAMINA_API JordanDecompositionResult jordan_form_checked(
    const std::shared_ptr<SymbolicExpr>& A
);

LAMINA_API bool jordan_form(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& J,
    std::shared_ptr<SymbolicExpr>& P
);

} // namespace lamina
