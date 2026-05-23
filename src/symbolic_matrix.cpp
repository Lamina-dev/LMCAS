#include "../include/symbolic_matrix.hpp"
#include "../include/symbolic.hpp"
#include <vector>
#include <cmath>

namespace lamina {

std::shared_ptr<SymbolicExpr> matrix_multiply(const std::shared_ptr<SymbolicExpr>& A, const std::shared_ptr<SymbolicExpr>& B) {
    return SymbolicExpr::multiply(A, B);
}

std::shared_ptr<SymbolicExpr> matrix_determinant(const std::shared_ptr<SymbolicExpr>& A) {
    return SymbolicExpr::determinant(A);
}

std::shared_ptr<SymbolicExpr> matrix_inverse(const std::shared_ptr<SymbolicExpr>& A) {
    return SymbolicExpr::inverse(A);
}

std::shared_ptr<SymbolicExpr> matrix_rotation(double theta, int dim) {
    if (dim == 2) {
        auto c = SymbolicExpr::number(std::cos(theta));
        auto s = SymbolicExpr::number(std::sin(theta));
        return SymbolicExpr::matrix({{c, SymbolicExpr::multiply(SymbolicExpr::number(-1), s)}, {s, c}});
    }

    return nullptr;
}

std::shared_ptr<SymbolicExpr> matrix_reflection(double angle, int dim) {
    if (dim == 2) {
        auto c = SymbolicExpr::number(std::cos(angle));
        auto s = SymbolicExpr::number(std::sin(angle));
        return SymbolicExpr::matrix({{SymbolicExpr::add(c, SymbolicExpr::multiply(s, s)), SymbolicExpr::multiply(SymbolicExpr::number(-1), s)}, {s, SymbolicExpr::add(c, SymbolicExpr::multiply(s, s))}});
    }
    return nullptr;
}

std::shared_ptr<SymbolicExpr> matrix_scaling(double sx, double sy, int dim) {
    if (dim == 2) {
        return SymbolicExpr::matrix({{SymbolicExpr::number(sx), SymbolicExpr::number(0)}, {SymbolicExpr::number(0), SymbolicExpr::number(sy)}});
    }
    return nullptr;
}

std::vector<std::shared_ptr<SymbolicExpr>> matrix_eigenvalues(const std::shared_ptr<SymbolicExpr>& A) {

    return {};
}

std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> matrix_eigenvectors(const std::shared_ptr<SymbolicExpr>& A) {

    return {};
}

}
