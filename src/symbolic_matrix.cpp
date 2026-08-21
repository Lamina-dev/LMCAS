#include "../include/symbolic_matrix.hpp"
#include "symbolic_ast.hpp"
#include "../include/symbolic.hpp"
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

ExpressionResult unsupported_matrix_dimension(const std::string& operation)
{
    return ExpressionResult::failure(
        CasErrc::UnsupportedExpression,
        "only two-dimensional transformation matrices are supported", operation);
}

bool matrix_is_provably_nonzero_number(const std::shared_ptr<SymbolicExpr>& expr)
{
    if (!expr || !lamina::detail::node(expr)) return false;
    auto simplified = expr->simplify();
    return simplified && lamina::detail::node(simplified) &&
           lamina::detail::node(simplified)->is_number() &&
           !lamina::detail::node(simplified)->is_zero();
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

std::shared_ptr<SymbolicExpr> matrix_multiply(const std::shared_ptr<SymbolicExpr>& A,
                                              const std::shared_ptr<SymbolicExpr>& B) {
    auto checked = matrix_multiply_checked(A, B);
    if (!checked) return nullptr;
    return checked.value();
}

ExpressionResult matrix_determinant_checked(const std::shared_ptr<SymbolicExpr>& A,
                                            ComputationContext& context) {
    const std::string operation = "matrix_determinant";
    auto mat = require_square_matrix(A, context, operation);
    if (!mat) return ExpressionResult::failure(mat.error());
    return require_matrix_result(SymbolicExpr::determinant(A), operation);
}

ExpressionResult matrix_determinant_checked(const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_determinant_checked(A, context);
}

std::shared_ptr<SymbolicExpr> matrix_determinant(const std::shared_ptr<SymbolicExpr>& A) {
    auto checked = matrix_determinant_checked(A);
    if (!checked) return nullptr;
    return checked.value();
}

ExpressionResult matrix_inverse_checked(const std::shared_ptr<SymbolicExpr>& A,
                                        ComputationContext& context) {
    const std::string operation = "matrix_inverse";
    auto mat = require_square_matrix(A, context, operation);
    if (!mat) return ExpressionResult::failure(mat.error());
    auto determinant_budget =
        context.consume_steps(mat.value()->rows() * mat.value()->cols() + 1, operation);
    if (!determinant_budget) return ExpressionResult::failure(determinant_budget.error());

    auto det = SymbolicExpr::determinant(A);
    if (!det || !lamina::detail::node(det)) {
        return ExpressionResult::failure(
            CasErrc::InternalInvariant,
            "matrix determinant could not be constructed",
            operation);
    }
    auto det_simplified = det->simplify();
    if (!det_simplified || !lamina::detail::node(det_simplified)) {
        return ExpressionResult::failure(
            CasErrc::InternalInvariant,
            "matrix determinant simplification produced an empty expression",
            operation);
    }
    if (lamina::detail::node(det_simplified)->is_zero()) {
        return ExpressionResult::failure(
            CasErrc::DomainError,
            "matrix is singular and has no inverse",
            operation);
    }
    if (!matrix_is_provably_nonzero_number(det_simplified)) {
        return ExpressionResult::failure(
            CasErrc::Inconclusive,
            "matrix determinant is not provably non-zero in the current support domain",
            operation);
    }

    auto inverse_budget =
        context.consume_steps(mat.value()->rows() * mat.value()->cols() * 4 + 1, operation);
    if (!inverse_budget) return ExpressionResult::failure(inverse_budget.error());
    return require_matrix_result(SymbolicExpr::inverse(A), operation);
}

ExpressionResult matrix_inverse_checked(const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_inverse_checked(A, context);
}

std::shared_ptr<SymbolicExpr> matrix_inverse(const std::shared_ptr<SymbolicExpr>& A) {
    auto checked = matrix_inverse_checked(A);
    if (!checked) return nullptr;
    return checked.value();
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

std::shared_ptr<SymbolicExpr> matrix_rotation(double theta, int dim) {
    auto checked = matrix_rotation_checked(theta, dim);
    if (!checked) return nullptr;
    return checked.value();
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

std::shared_ptr<SymbolicExpr> matrix_reflection(double angle, int dim) {
    auto checked = matrix_reflection_checked(angle, dim);
    if (!checked) return nullptr;
    return checked.value();
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

std::shared_ptr<SymbolicExpr> matrix_scaling(double sx, double sy, int dim) {
    auto checked = matrix_scaling_checked(sx, sy, dim);
    if (!checked) return nullptr;
    return checked.value();
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

std::vector<std::shared_ptr<SymbolicExpr>> matrix_eigenvalues(const std::shared_ptr<SymbolicExpr>& A) {
    auto checked = matrix_eigenvalues_checked(A);
    if (!checked) return {};
    return std::move(checked.value());
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

std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> matrix_eigenvectors(const std::shared_ptr<SymbolicExpr>& A) {
    auto checked = matrix_eigenvectors_checked(A);
    if (!checked) return {};
    return std::move(checked.value());
}

}
