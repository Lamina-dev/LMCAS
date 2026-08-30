#include "../include/symbolic_matrix.hpp"
#include "symbolic_ast.hpp"
#include "../include/symbolic.hpp"
#include "internal/exact_matrix.hpp"
#include <vector>
#include <cmath>

namespace lamina {

namespace {

Result<void> matrix_consume_step(ComputationContext& context,
                                 const std::string& operation)
{
    return context.consume_steps(1, operation);
}

Result<std::shared_ptr<const MatrixNode>> require_matrix(
    const std::shared_ptr<SymbolicExpr>& expr,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = matrix_consume_step(context, operation);
    if (!step) return Result<std::shared_ptr<const MatrixNode>>::failure(step.error());
    if (!expr || !lamina::detail::node(expr)) {
        return Result<std::shared_ptr<const MatrixNode>>::failure(
            CasErrc::InvalidArgument, "matrix expression cannot be null", operation);
    }
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(expr));
    if (!mat) {
        return Result<std::shared_ptr<const MatrixNode>>::failure(
            CasErrc::InvalidArgument, "expression is not a matrix", operation);
    }
    return Result<std::shared_ptr<const MatrixNode>>::success(mat);
}

Result<std::shared_ptr<const MatrixNode>> require_square_matrix(
    const std::shared_ptr<SymbolicExpr>& expr,
    ComputationContext& context,
    const std::string& operation)
{
    auto mat = require_matrix(expr, context, operation);
    if (!mat) return mat;
    if (mat.value()->rows() != mat.value()->cols()) {
        return Result<std::shared_ptr<const MatrixNode>>::failure(
            CasErrc::InvalidArgument, "matrix must be square", operation);
    }
    return mat;
}

ExpressionResult require_matrix_result(std::shared_ptr<SymbolicExpr> result,
                                       const std::string& operation)
{
    if (!result || !lamina::detail::node(result)) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                         "matrix operation returned an empty expression",
                                         operation);
    }
    return ExpressionResult::success(std::move(result));
}

detail::ExactMatrixData exact_matrix_data(const MatrixNode& matrix) {
    detail::ExactMatrixData data{matrix.rows(), matrix.cols(), {}};
    data.entries.reserve(matrix.rows() * matrix.cols());
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t column = 0; column < matrix.cols(); ++column) {
            data.entries.push_back(
                lamina::detail::make_expression_ptr(matrix.get(row, column)));
        }
    }
    return data;
}

std::shared_ptr<SymbolicExpr> exact_matrix_expression(
    const detail::ExactMatrixData& matrix) {
    MatrixNode::DenseStorage entries;
    entries.reserve(matrix.entries.size());
    for (const auto& entry : matrix.entries) {
        entries.push_back(lamina::detail::node(entry));
    }
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<MatrixNode>(
            matrix.rows, matrix.cols, std::move(entries)));
}

ExpressionResult unsupported_matrix_dimension(const std::string& operation)
{
    return ExpressionResult::failure(
        CasErrc::UnsupportedExpression,
        "only two-dimensional transformation matrices are supported", operation);
}


} // namespace

ExpressionResult matrix_multiply_checked(const std::shared_ptr<SymbolicExpr>& A,
                                         const std::shared_ptr<SymbolicExpr>& B,
                                         ComputationContext& context) {
    const std::string operation = "matrix_multiply";
    auto left = require_matrix(A, context, operation);
    if (!left) return ExpressionResult::failure(left.error());
    auto right = require_matrix(B, context, operation);
    if (!right) return ExpressionResult::failure(right.error());
    if (left.value()->cols() != right.value()->rows()) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument, "matrix dimensions are not compatible", operation);
    }
    return require_matrix_result(SymbolicExpr::multiply(A, B), operation);
}

ExpressionResult matrix_multiply_checked(const std::shared_ptr<SymbolicExpr>& A,
                                         const std::shared_ptr<SymbolicExpr>& B) {
    ComputationContext context;
    return matrix_multiply_checked(A, B, context);
}


ExpressionResult matrix_determinant_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context) {
    const std::string operation = "matrix_determinant";
    auto matrix = require_square_matrix(A, context, operation);
    if (!matrix) return ExpressionResult::failure(matrix.error());
    auto determinant = detail::determinant_exact(
        exact_matrix_data(*matrix.value()), context, operation);
    if (!determinant) return ExpressionResult::failure(determinant.error());
    return ExpressionResult::success(std::move(determinant.value()));
}

ExpressionResult matrix_determinant_checked(const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_determinant_checked(A, context);
}


ExpressionResult matrix_inverse_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context) {
    const std::string operation = "matrix_inverse";
    auto matrix = require_square_matrix(A, context, operation);
    if (!matrix) return ExpressionResult::failure(matrix.error());
    auto inverse = detail::inverse_exact(
        exact_matrix_data(*matrix.value()), context, operation);
    if (!inverse) return ExpressionResult::failure(inverse.error());
    return ExpressionResult::success(
        exact_matrix_expression(inverse.value()));
}

ExpressionResult matrix_inverse_checked(const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_inverse_checked(A, context);
}


MatrixRankResult matrix_rank_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context) {
    const std::string operation = "matrix_rank";
    auto matrix = require_matrix(A, context, operation);
    if (!matrix) return MatrixRankResult::failure(matrix.error());
    auto rank = detail::rank_exact(
        exact_matrix_data(*matrix.value()), matrix.value()->cols(),
        context, operation);
    if (!rank) return MatrixRankResult::failure(rank.error());
    return MatrixRankResult::success(rank.value());
}

MatrixRankResult matrix_rank_checked(
    const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_rank_checked(A, context);
}

ExpressionResult matrix_rref_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context) {
    const std::string operation = "matrix_rref";
    auto matrix = require_matrix(A, context, operation);
    if (!matrix) return ExpressionResult::failure(matrix.error());
    auto rref = detail::rref_exact(
        exact_matrix_data(*matrix.value()), matrix.value()->cols(),
        context, operation);
    if (!rref) return ExpressionResult::failure(rref.error());
    return ExpressionResult::success(exact_matrix_expression(rref.value()));
}

ExpressionResult matrix_rref_checked(
    const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_rref_checked(A, context);
}

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
                lamina::detail::make_expression_ptr(
                    matrix.value()->get(row, column)));
        }
        augmented.entries.push_back(
            lamina::detail::make_expression_ptr(rhs.value()->get(row, 0)));
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

ExpressionResult matrix_rotation_checked(double theta, int dim,
                                         ComputationContext& context) {
    const std::string operation = "matrix_rotation";
    auto step = matrix_consume_step(context, operation);
    if (!step) return ExpressionResult::failure(step.error());
    if (dim == 2) {
        auto c = SymbolicExpr::number(std::cos(theta));
        auto s = SymbolicExpr::number(std::sin(theta));
        return require_matrix_result(
            SymbolicExpr::matrix({{c, SymbolicExpr::multiply(SymbolicExpr::number(-1), s)}, {s, c}}),
            operation);
    }

    return unsupported_matrix_dimension(operation);
}

ExpressionResult matrix_rotation_checked(double theta, int dim) {
    ComputationContext context;
    return matrix_rotation_checked(theta, dim, context);
}


ExpressionResult matrix_reflection_checked(double angle, int dim,
                                           ComputationContext& context) {
    const std::string operation = "matrix_reflection";
    auto step = matrix_consume_step(context, operation);
    if (!step) return ExpressionResult::failure(step.error());
    if (dim == 2) {
        auto c = SymbolicExpr::number(std::cos(angle));
        auto s = SymbolicExpr::number(std::sin(angle));
        return require_matrix_result(
            SymbolicExpr::matrix({{SymbolicExpr::add(c, SymbolicExpr::multiply(s, s)), SymbolicExpr::multiply(SymbolicExpr::number(-1), s)}, {s, SymbolicExpr::add(c, SymbolicExpr::multiply(s, s))}}),
            operation);
    }
    return unsupported_matrix_dimension(operation);
}

ExpressionResult matrix_reflection_checked(double angle, int dim) {
    ComputationContext context;
    return matrix_reflection_checked(angle, dim, context);
}


ExpressionResult matrix_scaling_checked(double sx, double sy, int dim,
                                        ComputationContext& context) {
    const std::string operation = "matrix_scaling";
    auto step = matrix_consume_step(context, operation);
    if (!step) return ExpressionResult::failure(step.error());
    if (dim == 2) {
        return require_matrix_result(
            SymbolicExpr::matrix({{SymbolicExpr::number(sx), SymbolicExpr::number(0)}, {SymbolicExpr::number(0), SymbolicExpr::number(sy)}}),
            operation);
    }
    return unsupported_matrix_dimension(operation);
}

ExpressionResult matrix_scaling_checked(double sx, double sy, int dim) {
    ComputationContext context;
    return matrix_scaling_checked(sx, sy, dim, context);
}


MatrixEigenvalueResult matrix_eigenvalues_checked(const std::shared_ptr<SymbolicExpr>& A,
                                                  ComputationContext& context) {
    const std::string operation = "matrix_eigenvalues";
    auto mat_result = require_square_matrix(A, context, operation);
    if (!mat_result) return MatrixEigenvalueResult::failure(mat_result.error());
    auto mat = mat_result.value();

    bool triangular = true;
    for (size_t r = 0; r < mat->rows() && triangular; ++r) {
        for (size_t c = 0; c < r; ++c) {
            auto below = lamina::detail::make_expression_ptr(mat->get(r, c))->simplify();
            if (!lamina::detail::node(below)->is_zero()) {
                triangular = false;
                break;
            }
        }
    }
    if (triangular) {
        std::vector<std::shared_ptr<SymbolicExpr>> values;
        values.reserve(mat->rows());
        for (size_t i = 0; i < mat->rows(); ++i) {
            values.push_back(lamina::detail::make_expression_ptr(mat->get(i, i))->simplify());
        }
        return MatrixEigenvalueResult::success(std::move(values));
    }

    auto values_matrix = SymbolicExpr::eigenvalues(A);
    auto values_node = values_matrix ? std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(values_matrix)) : nullptr;
    if (!values_node || values_node->rows() == 0) {
        return MatrixEigenvalueResult::failure(
            CasErrc::UnsupportedExpression,
            "eigenvalue computation did not produce a matrix result",
            operation);
    }

    std::vector<std::shared_ptr<SymbolicExpr>> values;
    for (size_t c = 0; c < values_node->cols(); ++c) {
        values.push_back(lamina::detail::make_expression_ptr(values_node->get(0, c))->simplify());
    }
    return MatrixEigenvalueResult::success(std::move(values));
}

MatrixEigenvalueResult matrix_eigenvalues_checked(const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_eigenvalues_checked(A, context);
}


MatrixEigenvectorResult matrix_eigenvectors_checked(const std::shared_ptr<SymbolicExpr>& A,
                                                    ComputationContext& context) {
    const std::string operation = "matrix_eigenvectors";
    auto mat_result = require_square_matrix(A, context, operation);
    if (!mat_result) return MatrixEigenvectorResult::failure(mat_result.error());
    auto mat = mat_result.value();

    if (mat->rows() == 2) {
        auto a = lamina::detail::make_expression_ptr(mat->get(0, 0))->simplify();
        auto b = lamina::detail::make_expression_ptr(mat->get(0, 1))->simplify();
        auto lower = lamina::detail::make_expression_ptr(mat->get(1, 0))->simplify();
        auto d = lamina::detail::make_expression_ptr(mat->get(1, 1))->simplify();
        auto delta = SymbolicExpr::add(d, SymbolicExpr::multiply(SymbolicExpr::number(-1), a))->simplify();
        if (lamina::detail::node(lower)->is_zero() && !lamina::detail::node(delta)->is_zero()) {
            return MatrixEigenvectorResult::success({
                {SymbolicExpr::number(1), SymbolicExpr::number(0)},
                {b, delta}
            });
        }
    }

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> vectors;
    auto pairs = SymbolicExpr::eigenvectors(A);
    for (const auto& [lambda, vector_matrices] : pairs) {
        (void)lambda;
        for (const auto& vector_matrix : vector_matrices) {
            auto vec_node = vector_matrix ? std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(vector_matrix)) : nullptr;
            if (!vec_node || vec_node->cols() != 1) continue;
            std::vector<std::shared_ptr<SymbolicExpr>> vector;
            for (size_t r = 0; r < vec_node->rows(); ++r) {
                vector.push_back(lamina::detail::make_expression_ptr(vec_node->get(r, 0))->simplify());
            }
            vectors.push_back(std::move(vector));
        }
    }
    if (vectors.empty()) {
        return MatrixEigenvectorResult::failure(
            CasErrc::UnsupportedExpression,
            "eigenvector computation did not produce supported vectors",
            operation);
    }
    return MatrixEigenvectorResult::success(std::move(vectors));
}

MatrixEigenvectorResult matrix_eigenvectors_checked(const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_eigenvectors_checked(A, context);
}


}
