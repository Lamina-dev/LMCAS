#pragma once
#include "symbolic.hpp"
#include <memory>
#include <vector>
#include <string>

namespace lamina {
// 矩阵乘法
std::shared_ptr<SymbolicExpr> matrix_multiply(
    const std::shared_ptr<SymbolicExpr>& A,
    const std::shared_ptr<SymbolicExpr>& B
);

// 行列式
std::shared_ptr<SymbolicExpr> matrix_determinant(
    const std::shared_ptr<SymbolicExpr>& A
);

// 逆矩阵
std::shared_ptr<SymbolicExpr> matrix_inverse(
    const std::shared_ptr<SymbolicExpr>& A
);

// 线性变换：旋转
std::shared_ptr<SymbolicExpr> matrix_rotation(
    double theta,
    int dim = 2
);

// 线性变换：反射
std::shared_ptr<SymbolicExpr> matrix_reflection(
    double angle,
    int dim = 2
);

// 线性变换：伸缩
std::shared_ptr<SymbolicExpr> matrix_scaling(
    double sx,
    double sy,
    int dim = 2
);

// 特征值/特征向量
std::vector<std::shared_ptr<SymbolicExpr>> matrix_eigenvalues(
    const std::shared_ptr<SymbolicExpr>& A
);
std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> matrix_eigenvectors(
    const std::shared_ptr<SymbolicExpr>& A
);

}
