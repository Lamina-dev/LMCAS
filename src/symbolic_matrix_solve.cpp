#include "internal/symbolic_matrix_support.hpp"

namespace LMCAS {
using namespace detail::matrix_api;
MatrixNullspaceResult matrix_nullspace_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context) {
    const std::string operation = "matrix_nullspace";
    auto matrix = require_matrix(A, context, operation);
    if (!matrix) return MatrixNullspaceResult::failure(matrix.error());
    auto nullspace = detail::nullspace_exact(
        exact_matrix_data(*matrix.value()), context, operation);
    if (!nullspace) return MatrixNullspaceResult::failure(nullspace.error());
    return MatrixNullspaceResult::success(std::move(nullspace.value()));
}

MatrixNullspaceResult matrix_nullspace_checked(
    const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_nullspace_checked(A, context);
}

MatrixLinearSolveResult matrix_solve_linear_checked(
    const std::shared_ptr<SymbolicExpr>& coefficients,
    const std::shared_ptr<SymbolicExpr>& right_hand_side,
    ComputationContext& context) {
    const std::string operation = "matrix_solve_linear";
    auto matrix = require_matrix(coefficients, context, operation);
    if (!matrix) return MatrixLinearSolveResult::failure(matrix.error());
    auto rhs = require_matrix(right_hand_side, context, operation);
    if (!rhs) return MatrixLinearSolveResult::failure(rhs.error());
    if (rhs.value()->rows() != matrix.value()->rows() ||
        rhs.value()->cols() != 1) {
        return MatrixLinearSolveResult::failure(
            CasErrc::InvalidArgument,
            "right-hand side must be a matching column matrix", operation);
    }

    detail::ExactMatrixData augmented{
        matrix.value()->rows(), matrix.value()->cols() + 1, {}};
    augmented.entries.reserve(augmented.rows * augmented.cols);
    for (std::size_t row = 0; row < augmented.rows; ++row) {
        for (std::size_t column = 0;
             column < matrix.value()->cols(); ++column) {
            augmented.entries.push_back(
                LMCAS::detail::make_expression_ptr(
                    matrix.value()->get(row, column)));
        }
        augmented.entries.push_back(
            LMCAS::detail::make_expression_ptr(rhs.value()->get(row, 0)));
    }
    auto solved = detail::solve_linear_exact(
        std::move(augmented), matrix.value()->cols(), context, operation);
    if (!solved) return MatrixLinearSolveResult::failure(solved.error());
    if (auto* unique =
            std::get_if<detail::UniqueLinearSolution>(&solved.value())) {
        return MatrixLinearSolveResult::success(
            MatrixUniqueLinearSolution{std::move(unique->values)});
    }
    if (auto* parametric =
            std::get_if<detail::ParametricLinearSolution>(&solved.value())) {
        return MatrixLinearSolveResult::success(
            MatrixParametricLinearSolution{
                std::move(parametric->particular),
                std::move(parametric->nullspace_basis),
                std::move(parametric->free_columns)});
    }
    return MatrixLinearSolveResult::success(
        MatrixInconsistentLinearSolution{});
}

MatrixLinearSolveResult matrix_solve_linear_checked(
    const std::shared_ptr<SymbolicExpr>& coefficients,
    const std::shared_ptr<SymbolicExpr>& right_hand_side) {
    ComputationContext context;
    return matrix_solve_linear_checked(
        coefficients, right_hand_side, context);
}

}
