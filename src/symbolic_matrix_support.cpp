#include "internal/symbolic_matrix_support.hpp"

namespace LMCAS::detail::matrix_api {

Result<std::shared_ptr<const MatrixNode>> require_matrix(
    const std::shared_ptr<SymbolicExpr>& expr,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return Result<std::shared_ptr<const MatrixNode>>::failure(step.error());
    if (!expr || !LMCAS::detail::node(expr)) {
        return Result<std::shared_ptr<const MatrixNode>>::failure(
            CasErrc::InvalidArgument, "matrix expression cannot be null", operation);
    }
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(expr));
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
    if (!result || !LMCAS::detail::node(result)) {
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
                LMCAS::detail::make_expression_ptr(matrix.get(row, column)));
        }
    }
    return data;
}

std::shared_ptr<SymbolicExpr> exact_matrix_expression(
    const detail::ExactMatrixData& matrix) {
    MatrixNode::DenseStorage entries;
    entries.reserve(matrix.entries.size());
    for (const auto& entry : matrix.entries) {
        entries.push_back(LMCAS::detail::node(entry));
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<MatrixNode>(
            matrix.rows, matrix.cols, std::move(entries)));
}

}
