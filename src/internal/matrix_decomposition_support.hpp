#pragma once

#include "matrix_decomposition.hpp"
#include "internal/symbolic_ast.hpp"

#include <optional>

namespace LMCAS::detail::matrix_decomposition {

Result<std::shared_ptr<const MatrixNode>> validate_decomposition_matrix(
    const std::shared_ptr<SymbolicExpr>& A, bool require_square,
    ComputationContext& context, const std::string& operation);
Result<void> validate_decomposition_output(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& name, const std::string& operation);
std::optional<Rational> exact_rational_expr(
    const std::shared_ptr<const SymbolicNode>& node);
std::optional<Rational> exact_rational_expr(
    const std::shared_ptr<SymbolicExpr>& expr);
std::shared_ptr<SymbolicExpr> exact_number_expr(const Rational& value);
std::optional<std::vector<std::vector<Rational>>> exact_rational_matrix(
    const std::shared_ptr<const MatrixNode>& matrix);
std::shared_ptr<SymbolicExpr> identity_matrix_expr(std::size_t n);
std::optional<std::vector<Rational>> exact_rectangular_diagonal_entries(
    const std::shared_ptr<const MatrixNode>& matrix);
std::optional<std::vector<Rational>> exact_diagonal_entries(
    const std::shared_ptr<const MatrixNode>& matrix);
std::shared_ptr<SymbolicExpr> rectangular_diagonal_expr(
    std::size_t rows, std::size_t cols, const std::vector<Rational>& diagonal);

}
