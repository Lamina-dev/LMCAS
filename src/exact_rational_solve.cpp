#include "internal/exact_rational_matrix.hpp"

#include <algorithm>

namespace LMCAS::detail {
namespace {

void reduce_rational_pivot(
    std::vector<Rational>& matrix, std::size_t rows, std::size_t columns,
    std::size_t pivot_row, std::size_t column, std::size_t first_target) {
    const Rational pivot = matrix[pivot_row * columns + column];
    for (std::size_t target = first_target; target < columns; ++target) {
        matrix[pivot_row * columns + target] =
            matrix[pivot_row * columns + target] / pivot;
    }
    for (std::size_t row = 0; row < rows; ++row) {
        if (row == pivot_row) continue;
        const Rational factor = matrix[row * columns + column];
        if (factor == Rational(0)) continue;
        for (std::size_t target = first_target; target < columns; ++target) {
            matrix[row * columns + target] =
                matrix[row * columns + target] -
                factor * matrix[pivot_row * columns + target];
        }
    }
}

}

std::vector<std::optional<std::size_t>> reduce_rational_rows(
    std::size_t rows, std::size_t coefficient_columns,
    std::vector<Rational>& matrix, bool start_at_pivot) {
    const std::size_t columns = coefficient_columns + 1;
    std::vector<std::optional<std::size_t>> pivot_rows(coefficient_columns);
    std::size_t pivot_row = 0;
    for (std::size_t column = 0;
         column < coefficient_columns && pivot_row < rows; ++column) {
        std::size_t selected = pivot_row;
        while (selected < rows &&
               matrix[selected * columns + column] == Rational(0)) {
            ++selected;
        }
        if (selected == rows) continue;
        if (selected != pivot_row) {
            for (std::size_t target = 0; target < columns; ++target) {
                std::swap(matrix[selected * columns + target],
                          matrix[pivot_row * columns + target]);
            }
        }
        reduce_rational_pivot(matrix, rows, columns, pivot_row, column,
                              start_at_pivot ? column : 0);
        pivot_rows[column] = pivot_row++;
    }
    return pivot_rows;
}

bool rational_system_inconsistent(
    std::size_t rows, std::size_t coefficient_columns,
    const std::vector<Rational>& matrix) {
    const std::size_t columns = coefficient_columns + 1;
    for (std::size_t row = 0; row < rows; ++row) {
        bool zero = true;
        for (std::size_t column = 0; column < coefficient_columns; ++column) {
            zero = zero && matrix[row * columns + column] == Rational(0);
        }
        if (zero && matrix[row * columns + coefficient_columns] != Rational(0)) {
            return true;
        }
    }
    return false;
}

Result<std::vector<Rational>> solve_rational_unique(
    std::size_t rows,
    std::size_t coefficient_columns,
    std::vector<Rational> matrix,
    ComputationContext& context,
    const std::string& operation) {
    const std::size_t columns = coefficient_columns + 1;
    if (rows == 0 || coefficient_columns == 0 ||
        matrix.size() != rows * columns) {
        return Result<std::vector<Rational>>::failure(
            CasErrc::InvalidArgument,
            "rational linear-system storage is invalid", operation);
    }
    auto budget = context.consume_steps(rows * columns, operation);
    if (!budget) return Result<std::vector<Rational>>::failure(budget.error());
    auto pivot_rows = reduce_rational_rows(rows, coefficient_columns, matrix, true);
    if (rational_system_inconsistent(rows, coefficient_columns, matrix)) {
        return Result<std::vector<Rational>>::failure(
            CasErrc::DomainError,
            "rational linear system is inconsistent", operation);
    }
    std::vector<Rational> solution(coefficient_columns, Rational(0));
    for (std::size_t column = 0; column < coefficient_columns; ++column) {
        if (!pivot_rows[column]) {
            return Result<std::vector<Rational>>::failure(
                CasErrc::Inconclusive,
                "rational linear system does not have a unique solution",
                operation);
        }
        solution[column] =
            matrix[*pivot_rows[column] * columns + coefficient_columns];
    }
    return Result<std::vector<Rational>>::success(std::move(solution));
}

}
