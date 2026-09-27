#pragma once

#include "symbolic_matrix.hpp"
#include "internal/exact_matrix.hpp"
#include "internal/symbolic_ast.hpp"
#include <utility>

namespace LMCAS::detail::matrix_api {

Result<std::shared_ptr<const MatrixNode>> require_matrix(
    const std::shared_ptr<SymbolicExpr>& expr, ComputationContext& context,
    const std::string& operation);
Result<std::shared_ptr<const MatrixNode>> require_square_matrix(
    const std::shared_ptr<SymbolicExpr>& expr, ComputationContext& context,
    const std::string& operation);
ExpressionResult require_matrix_result(std::shared_ptr<SymbolicExpr> result,
                                       const std::string& operation);
ExactMatrixData exact_matrix_data(const MatrixNode& matrix);
std::shared_ptr<SymbolicExpr> exact_matrix_expression(const ExactMatrixData& matrix);

}

namespace LMCAS::detail {
Result<std::vector<std::pair<ExprPtr, std::vector<std::vector<ExprPtr>>>>>
matrix_eigenspaces_checked(const ExprPtr& matrix, ComputationContext& context);
}
