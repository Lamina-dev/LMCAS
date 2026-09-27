#include "internal/exact_matrix_support.hpp"
#include "internal/symbolic_ast.hpp"

#include <algorithm>
#include <limits>
#include <new>
#include <set>

namespace LMCAS::detail {
using namespace matrix_kernel;
namespace {
void subtract_pivot_row(ExactMatrixData& matrix, std::size_t row,
                        std::size_t pivot_row, const ExprPtr& factor) {
    for (std::size_t target = 0; target < matrix.cols; ++target) {
        matrix.at(row, target) = exact_subtract(
            matrix.at(row, target),
            exact_multiply(factor, matrix.at(pivot_row, target)));
    }
}

Result<void> eliminate_below(ExactMatrixData& matrix, std::size_t pivot_row,
                            std::size_t column, const ExprPtr& pivot_value,
                            ComputationContext& context,
                            const std::string& operation) {
    for (std::size_t row = pivot_row + 1; row < matrix.rows; ++row) {
        auto proof = classify_exact_zero(matrix.at(row, column), context, operation);
        if (!proof) return Result<void>::failure(proof.error());
        if (proof.value() == ZeroProof::Zero) continue;
        auto factor = exact_divide(matrix.at(row, column), pivot_value);
        subtract_pivot_row(matrix, row, pivot_row, factor);
    }
    return Result<void>::success();
}

Result<void> reduce_above(EliminationResult& result, ComputationContext& context,
                         const std::string& operation) {
    for (std::size_t pivot_index = result.rank; pivot_index-- > 0;) {
        const std::size_t row = pivot_index;
        const std::size_t column = result.pivot_columns[pivot_index];
        auto pivot_value = simplify_expr(result.matrix.at(row, column));
        for (std::size_t target = 0; target < result.matrix.cols; ++target) {
            result.matrix.at(row, target) = exact_divide(
                result.matrix.at(row, target), pivot_value);
        }
        for (std::size_t upper = 0; upper < row; ++upper) {
            auto proof = classify_exact_zero(
                result.matrix.at(upper, column), context, operation);
            if (!proof) return Result<void>::failure(proof.error());
            if (proof.value() == ZeroProof::Zero) continue;
            auto factor = result.matrix.at(upper, column);
            subtract_pivot_row(result.matrix, upper, row, factor);
        }
    }
    return Result<void>::success();
}
}

namespace {
Result<void> eliminate_forward(
    EliminationResult& result, std::size_t coefficient_columns,
    ComputationContext& context, const std::string& operation) {
    std::size_t pivot_row = 0;
    for (std::size_t column = 0;
         column < coefficient_columns && pivot_row < result.matrix.rows;
         ++column) {
        auto pivot = choose_pivot(
            result.matrix, pivot_row, column, context, operation);
        if (!pivot) return Result<void>::failure(pivot.error());
        if (!pivot.value()) continue;
        if (*pivot.value() != pivot_row) {
            swap_rows(result.matrix, *pivot.value(), pivot_row);
            ++result.row_swaps;
        }
        auto pivot_value = simplify_expr(result.matrix.at(pivot_row, column));
        result.pivot_columns.push_back(column);
        result.pivot_values.push_back(pivot_value);
        auto reduced = eliminate_below(
            result.matrix, pivot_row, column, pivot_value, context, operation);
        if (!reduced) return Result<void>::failure(reduced.error());
        ++pivot_row;
    }
    result.rank = result.pivot_columns.size();
    return Result<void>::success();
}
}

Result<EliminationResult> eliminate_exact(
    ExactMatrixData input,
    std::size_t coefficient_columns,
    EliminationForm form,
    ComputationContext& context,
    const std::string& operation) {
    auto valid = validate_matrix(input, operation);
    if (!valid) return Result<EliminationResult>::failure(valid.error());
    if (coefficient_columns == 0 || coefficient_columns > input.cols) {
        return Result<EliminationResult>::failure(
            CasErrc::InvalidArgument,
            "coefficient column count is outside matrix shape", operation);
    }
    try {
        EliminationResult result;
        result.matrix = std::move(input);
        auto forward = eliminate_forward(result, coefficient_columns, context, operation);
        if (!forward) return Result<EliminationResult>::failure(forward.error());

        if (form == EliminationForm::ReducedRowEchelon) {
            auto reduced = reduce_above(result, context, operation);
            if (!reduced) return Result<EliminationResult>::failure(reduced.error());
        }
        return Result<EliminationResult>::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return Result<EliminationResult>::failure(
            CasErrc::ResourceLimit,
            "exact elimination allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<EliminationResult>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

Result<ExactMatrixData> rref_exact(
    ExactMatrixData input,
    std::size_t coefficient_columns,
    ComputationContext& context,
    const std::string& operation) {
    auto eliminated = eliminate_exact(
        std::move(input), coefficient_columns,
        EliminationForm::ReducedRowEchelon, context, operation);
    if (!eliminated) return Result<ExactMatrixData>::failure(eliminated.error());
    return Result<ExactMatrixData>::success(
        std::move(eliminated.value().matrix));
}

Result<std::size_t> rank_exact(
    ExactMatrixData input,
    std::size_t coefficient_columns,
    ComputationContext& context,
    const std::string& operation) {
    auto eliminated = eliminate_exact(
        std::move(input), coefficient_columns,
        EliminationForm::Echelon, context, operation);
    if (!eliminated) return Result<std::size_t>::failure(eliminated.error());
    return Result<std::size_t>::success(eliminated.value().rank);
}

Result<ExactMatrixData> inverse_exact(
    const ExactMatrixData& input,
    ComputationContext& context,
    const std::string& operation) {
    auto valid = validate_matrix(input, operation);
    if (!valid) return Result<ExactMatrixData>::failure(valid.error());
    if (input.rows != input.cols) {
        return Result<ExactMatrixData>::failure(
            CasErrc::InvalidArgument,
            "inverse requires a square matrix", operation);
    }
    const std::size_t n = input.rows;
    ExactMatrixData augmented{n, 2 * n, {}};
    augmented.entries.reserve(n * 2 * n);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t column = 0; column < n; ++column) {
            augmented.entries.push_back(input.at(row, column));
        }
        for (std::size_t column = 0; column < n; ++column) {
            augmented.entries.push_back(
                SymbolicExpr::number(row == column ? 1 : 0));
        }
    }
    auto eliminated = eliminate_exact(
        std::move(augmented), n,
        EliminationForm::ReducedRowEchelon, context, operation);
    if (!eliminated) return Result<ExactMatrixData>::failure(eliminated.error());
    if (eliminated.value().rank != n) {
        return Result<ExactMatrixData>::failure(
            CasErrc::DomainError,
            "matrix is singular and has no inverse", operation);
    }
    ExactMatrixData inverse{n, n, {}};
    inverse.entries.reserve(n * n);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t column = 0; column < n; ++column) {
            inverse.entries.push_back(
                eliminated.value().matrix.at(row, n + column));
        }
    }
    return Result<ExactMatrixData>::success(std::move(inverse));
}

}
