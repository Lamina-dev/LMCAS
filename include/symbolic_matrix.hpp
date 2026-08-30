/**
 * @file symbolic_matrix.hpp
 * @brief 符号矩阵运算：乘法、行列式、逆、旋转/反射/缩放矩阵、特征值。
 */
#pragma once
#include "computation_context.hpp"
#include "result.hpp"
#include "symbolic.hpp"
#include <memory>
#include <vector>
#include <string>
#include <variant>

namespace lamina {

using MatrixEigenvalueResult = Result<std::vector<std::shared_ptr<SymbolicExpr>>>;
using MatrixEigenvectorResult =
    Result<std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>>;

/**
 * @brief 计算两个符号矩阵的乘积
 * @param A 左矩阵
 * @param B 右矩阵
 * @return 乘积矩阵的符号表达式
 */
LAMINA_API ExpressionResult matrix_multiply_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    const std::shared_ptr<SymbolicExpr>& B,
    ComputationContext& context
);

LAMINA_API ExpressionResult matrix_multiply_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    const std::shared_ptr<SymbolicExpr>& B
);


/**
 * @brief 计算符号矩阵的行列式
 * @param A 输入矩阵
 * @return 行列式的符号表达式
 */
LAMINA_API ExpressionResult matrix_determinant_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context
);

LAMINA_API ExpressionResult matrix_determinant_checked(
    const std::shared_ptr<SymbolicExpr>& A
);


/**
 * @brief 计算符号矩阵的逆矩阵
 * @param A 输入矩阵
 * @return 逆矩阵的符号表达式
 */
LAMINA_API ExpressionResult matrix_inverse_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context
);

LAMINA_API ExpressionResult matrix_inverse_checked(
    const std::shared_ptr<SymbolicExpr>& A
);


using MatrixRankResult = Result<std::size_t>;
using MatrixNullspaceResult =
    Result<std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>>;

struct MatrixUniqueLinearSolution {
    std::vector<std::shared_ptr<SymbolicExpr>> values;
};

struct MatrixParametricLinearSolution {
    std::vector<std::shared_ptr<SymbolicExpr>> particular;
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> nullspace_basis;
    std::vector<std::size_t> free_columns;
};

struct MatrixInconsistentLinearSolution {};

using MatrixLinearSolution = std::variant<
    MatrixUniqueLinearSolution,
    MatrixParametricLinearSolution,
    MatrixInconsistentLinearSolution>;
using MatrixLinearSolveResult = Result<MatrixLinearSolution>;

LAMINA_API MatrixRankResult matrix_rank_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context);
LAMINA_API MatrixRankResult matrix_rank_checked(
    const std::shared_ptr<SymbolicExpr>& A);

LAMINA_API ExpressionResult matrix_rref_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context);
LAMINA_API ExpressionResult matrix_rref_checked(
    const std::shared_ptr<SymbolicExpr>& A);

LAMINA_API MatrixNullspaceResult matrix_nullspace_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context);
LAMINA_API MatrixNullspaceResult matrix_nullspace_checked(
    const std::shared_ptr<SymbolicExpr>& A);

LAMINA_API MatrixLinearSolveResult matrix_solve_linear_checked(
    const std::shared_ptr<SymbolicExpr>& coefficients,
    const std::shared_ptr<SymbolicExpr>& right_hand_side,
    ComputationContext& context);
LAMINA_API MatrixLinearSolveResult matrix_solve_linear_checked(
    const std::shared_ptr<SymbolicExpr>& coefficients,
    const std::shared_ptr<SymbolicExpr>& right_hand_side);

/**
 * @brief 生成旋转矩阵
 * @param theta 旋转角度（弧度）
 * @param dim 矩阵维度，默认为 2
 * @return 旋转矩阵的符号表达式
 */
LAMINA_API ExpressionResult matrix_rotation_checked(
    double theta,
    int dim,
    ComputationContext& context
);

LAMINA_API ExpressionResult matrix_rotation_checked(
    double theta,
    int dim = 2
);


/**
 * @brief 生成反射矩阵
 * @param angle 反射轴角度（弧度）
 * @param dim 矩阵维度，默认为 2
 * @return 反射矩阵的符号表达式
 */
LAMINA_API ExpressionResult matrix_reflection_checked(
    double angle,
    int dim,
    ComputationContext& context
);

LAMINA_API ExpressionResult matrix_reflection_checked(
    double angle,
    int dim = 2
);


/**
 * @brief 生成缩放矩阵
 * @param sx x 方向缩放因子
 * @param sy y 方向缩放因子
 * @param dim 矩阵维度，默认为 2
 * @return 缩放矩阵的符号表达式
 */
LAMINA_API ExpressionResult matrix_scaling_checked(
    double sx,
    double sy,
    int dim,
    ComputationContext& context
);

LAMINA_API ExpressionResult matrix_scaling_checked(
    double sx,
    double sy,
    int dim = 2
);


/**
 * @brief 计算符号矩阵的特征值
 * @param A 输入矩阵
 * @return 特征值列表
 */
LAMINA_API MatrixEigenvalueResult matrix_eigenvalues_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context
);

LAMINA_API MatrixEigenvalueResult matrix_eigenvalues_checked(
    const std::shared_ptr<SymbolicExpr>& A
);


/**
 * @brief 计算符号矩阵的特征向量
 * @param A 输入矩阵
 * @return 特征向量列表，每个特征向量为一组分量
 */
LAMINA_API MatrixEigenvectorResult matrix_eigenvectors_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context
);

LAMINA_API MatrixEigenvectorResult matrix_eigenvectors_checked(
    const std::shared_ptr<SymbolicExpr>& A
);


}
