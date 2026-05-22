/**
 * @file symbolic_matrix.hpp
 * @brief 符号矩阵运算：乘法、行列式、逆、旋转/反射/缩放矩阵、特征值。
 */
#pragma once
#include "symbolic.hpp"
#include <memory>
#include <vector>
#include <string>

namespace lamina {

/**
 * @brief 计算两个符号矩阵的乘积
 * @param A 左矩阵
 * @param B 右矩阵
 * @return 乘积矩阵的符号表达式
 */
std::shared_ptr<SymbolicExpr> matrix_multiply(
    const std::shared_ptr<SymbolicExpr>& A,
    const std::shared_ptr<SymbolicExpr>& B
);

/**
 * @brief 计算符号矩阵的行列式
 * @param A 输入矩阵
 * @return 行列式的符号表达式
 */
std::shared_ptr<SymbolicExpr> matrix_determinant(
    const std::shared_ptr<SymbolicExpr>& A
);

/**
 * @brief 计算符号矩阵的逆矩阵
 * @param A 输入矩阵
 * @return 逆矩阵的符号表达式
 */
std::shared_ptr<SymbolicExpr> matrix_inverse(
    const std::shared_ptr<SymbolicExpr>& A
);

/**
 * @brief 生成旋转矩阵
 * @param theta 旋转角度（弧度）
 * @param dim 矩阵维度，默认为 2
 * @return 旋转矩阵的符号表达式
 */
std::shared_ptr<SymbolicExpr> matrix_rotation(
    double theta,
    int dim = 2
);

/**
 * @brief 生成反射矩阵
 * @param angle 反射轴角度（弧度）
 * @param dim 矩阵维度，默认为 2
 * @return 反射矩阵的符号表达式
 */
std::shared_ptr<SymbolicExpr> matrix_reflection(
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
std::shared_ptr<SymbolicExpr> matrix_scaling(
    double sx,
    double sy,
    int dim = 2
);

/**
 * @brief 计算符号矩阵的特征值
 * @param A 输入矩阵
 * @return 特征值列表
 */
std::vector<std::shared_ptr<SymbolicExpr>> matrix_eigenvalues(
    const std::shared_ptr<SymbolicExpr>& A
);

/**
 * @brief 计算符号矩阵的特征向量
 * @param A 输入矩阵
 * @return 特征向量列表，每个特征向量为一组分量
 */
std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> matrix_eigenvectors(
    const std::shared_ptr<SymbolicExpr>& A
);

}
