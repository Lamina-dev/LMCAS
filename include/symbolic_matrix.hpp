#pragma once
#include "symbolic.hpp"
#include <memory>
#include <vector>
#include <string>

namespace lamina {

std::shared_ptr<SymbolicExpr> matrix_multiply(
    const std::shared_ptr<SymbolicExpr>& A,
    const std::shared_ptr<SymbolicExpr>& B
);

std::shared_ptr<SymbolicExpr> matrix_determinant(
    const std::shared_ptr<SymbolicExpr>& A
);

std::shared_ptr<SymbolicExpr> matrix_inverse(
    const std::shared_ptr<SymbolicExpr>& A
);

std::shared_ptr<SymbolicExpr> matrix_rotation(
    double theta,
    int dim = 2
);

std::shared_ptr<SymbolicExpr> matrix_reflection(
    double angle,
    int dim = 2
);

std::shared_ptr<SymbolicExpr> matrix_scaling(
    double sx,
    double sy,
    int dim = 2
);

std::vector<std::shared_ptr<SymbolicExpr>> matrix_eigenvalues(
    const std::shared_ptr<SymbolicExpr>& A
);
std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> matrix_eigenvectors(
    const std::shared_ptr<SymbolicExpr>& A
);

}
