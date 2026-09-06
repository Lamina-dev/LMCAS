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

namespace LMCAS {

using MatrixEigenvalueResult = Result<std::vector<std::shared_ptr<SymbolicExpr>>>;
using MatrixEigenvectorResult =
    Result<std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>>;

/**
 * @brief 计算两个符号矩阵的乘积
 *
 * 归一化后的近似元素保留浮点计算值，精确元素保持精确算术。
 * 误差容限用于显式数值比较，不参与矩阵元素的归一化。
 * 矩阵因子按输入顺序相乘，结果形状为左矩阵行数乘以右矩阵列数。
 *
 * @param A 左矩阵
 * @param B 右矩阵
 * @return 乘积矩阵的符号表达式
 */
LMCAS_API ExpressionResult matrix_multiply_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    const std::shared_ptr<SymbolicExpr>& B,
    ComputationContext& context
);

LMCAS_API ExpressionResult matrix_multiply_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    const std::shared_ptr<SymbolicExpr>& B
);


/**
 * @brief 计算符号矩阵的行列式
 * @param A 输入矩阵
 * @return 行列式的符号表达式
 */
LMCAS_API ExpressionResult matrix_determinant_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context
);

LMCAS_API ExpressionResult matrix_determinant_checked(
    const std::shared_ptr<SymbolicExpr>& A
);


/**
 * @brief 计算符号矩阵的逆矩阵
 * @param A 输入矩阵
 * @return 逆矩阵的符号表达式
 */
LMCAS_API ExpressionResult matrix_inverse_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context
);

LMCAS_API ExpressionResult matrix_inverse_checked(
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

LMCAS_API MatrixRankResult matrix_rank_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context);
LMCAS_API MatrixRankResult matrix_rank_checked(
    const std::shared_ptr<SymbolicExpr>& A);

LMCAS_API ExpressionResult matrix_rref_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context);
LMCAS_API ExpressionResult matrix_rref_checked(
    const std::shared_ptr<SymbolicExpr>& A);

LMCAS_API MatrixNullspaceResult matrix_nullspace_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context);
LMCAS_API MatrixNullspaceResult matrix_nullspace_checked(
    const std::shared_ptr<SymbolicExpr>& A);

LMCAS_API MatrixLinearSolveResult matrix_solve_linear_checked(
    const std::shared_ptr<SymbolicExpr>& coefficients,
    const std::shared_ptr<SymbolicExpr>& right_hand_side,
    ComputationContext& context);
LMCAS_API MatrixLinearSolveResult matrix_solve_linear_checked(
    const std::shared_ptr<SymbolicExpr>& coefficients,
    const std::shared_ptr<SymbolicExpr>& right_hand_side);

/**
 * @brief 生成旋转矩阵
 * @param theta 旋转角度（弧度）
 * @param dim 矩阵维度，默认为 2
 * @return 旋转矩阵的符号表达式
 */
LMCAS_API ExpressionResult matrix_rotation_checked(
    double theta,
    int dim,
    ComputationContext& context
);

LMCAS_API ExpressionResult matrix_rotation_checked(
    double theta,
    int dim = 2
);


/**
 * @brief 生成反射矩阵
 * @param angle 反射轴角度（弧度）
 * @param dim 矩阵维度，默认为 2
 * @return 反射矩阵的符号表达式
 */
LMCAS_API ExpressionResult matrix_reflection_checked(
    double angle,
    int dim,
    ComputationContext& context
);

LMCAS_API ExpressionResult matrix_reflection_checked(
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
LMCAS_API ExpressionResult matrix_scaling_checked(
    double sx,
    double sy,
    int dim,
    ComputationContext& context
);

LMCAS_API ExpressionResult matrix_scaling_checked(
    double sx,
    double sy,
    int dim = 2
);


LMCAS_API ExpressionResult matrix_characteristic_polynomial_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    const std::string& variable,
    ComputationContext& context);

LMCAS_API ExpressionResult matrix_characteristic_polynomial_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    const std::string& variable = "lambda");
/**
 * @brief 计算符号矩阵的特征值
 * @param A 输入矩阵
 * @return 特征值列表
 */
LMCAS_API MatrixEigenvalueResult matrix_eigenvalues_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context
);

LMCAS_API MatrixEigenvalueResult matrix_eigenvalues_checked(
    const std::shared_ptr<SymbolicExpr>& A
);


/**
 * @brief 计算符号矩阵的特征向量
 * @param A 输入矩阵
 * @return 特征向量列表，每个特征向量为一组分量
 */
LMCAS_API MatrixEigenvectorResult matrix_eigenvectors_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context
);

LMCAS_API MatrixEigenvectorResult matrix_eigenvectors_checked(
    const std::shared_ptr<SymbolicExpr>& A
);


}
