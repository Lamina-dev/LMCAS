#include "../include/symbolic_matrix.hpp"
#include "../include/symbolic.hpp"
#include <vector>
#include <cmath>

namespace lamina {
// 矩阵乘法
std::shared_ptr<SymbolicExpr> matrix_multiply(const std::shared_ptr<SymbolicExpr>& A, const std::shared_ptr<SymbolicExpr>& B) {
    return SymbolicExpr::multiply(A, B);
}

// 行列式
std::shared_ptr<SymbolicExpr> matrix_determinant(const std::shared_ptr<SymbolicExpr>& A) {
    return SymbolicExpr::determinant(A);
}

// 逆矩阵
std::shared_ptr<SymbolicExpr> matrix_inverse(const std::shared_ptr<SymbolicExpr>& A) {
    return SymbolicExpr::inverse(A);
}

// 旋转矩阵（二维）
std::shared_ptr<SymbolicExpr> matrix_rotation(double theta, int dim) {
    if (dim == 2) {
        auto c = SymbolicExpr::number(std::cos(theta));
        auto s = SymbolicExpr::number(std::sin(theta));
        return SymbolicExpr::matrix({{c, SymbolicExpr::multiply(SymbolicExpr::number(-1), s)}, {s, c}});
    }
    // 三维旋转需指定轴，略
    return nullptr;
}

// 反射矩阵（二维，关于 y=tan(angle)x）
std::shared_ptr<SymbolicExpr> matrix_reflection(double angle, int dim) {
    if (dim == 2) {
        auto c = SymbolicExpr::number(std::cos(angle));
        auto s = SymbolicExpr::number(std::sin(angle));
        return SymbolicExpr::matrix({{SymbolicExpr::add(c, SymbolicExpr::multiply(s, s)), SymbolicExpr::multiply(SymbolicExpr::number(-1), s)}, {s, SymbolicExpr::add(c, SymbolicExpr::multiply(s, s))}});
    }
    return nullptr;
}

// 伸缩矩阵
std::shared_ptr<SymbolicExpr> matrix_scaling(double sx, double sy, int dim) {
    if (dim == 2) {
        return SymbolicExpr::matrix({{SymbolicExpr::number(sx), SymbolicExpr::number(0)}, {SymbolicExpr::number(0), SymbolicExpr::number(sy)}});
    }
    return nullptr;
}

// 特征值（仅支持2x2数值矩阵）
std::vector<std::shared_ptr<SymbolicExpr>> matrix_eigenvalues(const std::shared_ptr<SymbolicExpr>& A) {
    // 仅支持2x2数值矩阵
    // 形式：A = [[a, b], [c, d]]
    // λ^2 - (a+d)λ + (ad-bc) = 0
    // 返回 λ1, λ2
    // 需用户传入具体数值
    return {};
}

// 特征向量（仅支持2x2数值矩阵）
std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> matrix_eigenvectors(const std::shared_ptr<SymbolicExpr>& A) {
    // 仅支持2x2数值矩阵
    return {};
}

} // namespace lamina
