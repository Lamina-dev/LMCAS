#pragma once

#include "computation_context.hpp"
#include "result.hpp"
#include "symbolic.hpp"

#include <cstddef>
#include <string>
#include <variant>
#include <vector>

namespace lamina::detail {

enum class ZeroProof {
    Zero,
    NonZero,
    Unknown
};

Result<ZeroProof> classify_exact_zero(
    const ExprPtr& expression,
    ComputationContext& context,
    const std::string& operation);

struct ExactMatrixData {
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::vector<ExprPtr> entries;

    const ExprPtr& at(std::size_t row, std::size_t col) const {
        return entries.at(row * cols + col);
    }
    ExprPtr& at(std::size_t row, std::size_t col) {
        return entries.at(row * cols + col);
    }
};

enum class EliminationForm {
    Echelon,
    ReducedRowEchelon
};

struct EliminationResult {
    ExactMatrixData matrix;
    std::vector<std::size_t> pivot_columns;
    std::size_t rank = 0;
    std::size_t row_swaps = 0;
    std::vector<ExprPtr> pivot_values;
};

Result<EliminationResult> eliminate_exact(
    ExactMatrixData input,
    std::size_t coefficient_columns,
    EliminationForm form,
    ComputationContext& context,
    const std::string& operation);

Result<ExprPtr> determinant_exact(
    const ExactMatrixData& input,
    ComputationContext& context,
    const std::string& operation);

Result<ExactMatrixData> rref_exact(
    ExactMatrixData input,
    std::size_t coefficient_columns,
    ComputationContext& context,
    const std::string& operation);
Result<std::vector<Rational>> solve_rational_unique(
    std::size_t rows,
    std::size_t coefficient_columns,
    std::vector<Rational> augmented_entries,
    ComputationContext& context,
    const std::string& operation);


Result<std::size_t> rank_exact(
    ExactMatrixData input,
    std::size_t coefficient_columns,
    ComputationContext& context,
    const std::string& operation);

Result<ExactMatrixData> inverse_exact(
    const ExactMatrixData& input,
    ComputationContext& context,
    const std::string& operation);

struct UniqueLinearSolution {
    std::vector<ExprPtr> values;
};

struct ParametricLinearSolution {
    std::vector<ExprPtr> particular;
    std::vector<std::vector<ExprPtr>> nullspace_basis;
    std::vector<std::size_t> free_columns;
};

struct InconsistentLinearSolution {};

using ExactLinearSolution = std::variant<
    UniqueLinearSolution,
    ParametricLinearSolution,
    InconsistentLinearSolution>;

Result<ExactLinearSolution> solve_linear_exact(
    ExactMatrixData augmented,
    std::size_t coefficient_columns,
    ComputationContext& context,
    const std::string& operation);

Result<std::vector<std::vector<ExprPtr>>> nullspace_exact(
    ExactMatrixData coefficients,
    ComputationContext& context,
    const std::string& operation);

} // namespace lamina::detail
