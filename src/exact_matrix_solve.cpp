#include "internal/exact_matrix_support.hpp"
#include "internal/symbolic_ast.hpp"

#include <algorithm>
#include <limits>
#include <new>
#include <set>

namespace LMCAS::detail {
using namespace matrix_kernel;
namespace {
template <typename Entry>
Result<ExactLinearSolution> assemble_linear_solution(
    std::size_t coefficient_columns,
    const std::vector<std::optional<std::size_t>>& pivot_row, Entry entry) {
    std::vector<std::size_t> free_columns;
    for (std::size_t column = 0; column < coefficient_columns; ++column) {
        if (!pivot_row[column]) free_columns.push_back(column);
    }
    std::vector<ExprPtr> particular(
        coefficient_columns, SymbolicExpr::number(0));
    for (std::size_t column = 0; column < coefficient_columns; ++column) {
        if (pivot_row[column]) {
            particular[column] = entry(*pivot_row[column], coefficient_columns, false);
        }
    }
    if (free_columns.empty()) {
        return Result<ExactLinearSolution>::success(
            UniqueLinearSolution{std::move(particular)});
    }

    std::vector<std::vector<ExprPtr>> basis;
    basis.reserve(free_columns.size());
    for (const auto free_column : free_columns) {
        std::vector<ExprPtr> vector(
            coefficient_columns, SymbolicExpr::number(0));
        vector[free_column] = SymbolicExpr::number(1);
        for (std::size_t column = 0; column < coefficient_columns; ++column) {
            if (pivot_row[column]) {
                vector[column] = entry(*pivot_row[column], free_column, true);
            }
        }
        basis.push_back(std::move(vector));
    }
    return Result<ExactLinearSolution>::success(ParametricLinearSolution{
        std::move(particular), std::move(basis), std::move(free_columns)});
}

Result<bool> inconsistent_symbolic_system(
    const ExactMatrixData& rref, std::size_t coefficient_columns,
    ComputationContext& context, const std::string& operation) {
    for (std::size_t row = 0; row < rref.rows; ++row) {
        bool coefficients_zero = true;
        for (std::size_t column = 0; column < coefficient_columns; ++column) {
            auto proof = classify_exact_zero(rref.at(row, column), context, operation);
            if (!proof) return Result<bool>::failure(proof.error());
            if (proof.value() == ZeroProof::Unknown) {
                return Result<bool>::failure(
                    CasErrc::Inconclusive,
                    "linear-system consistency is not provable", operation);
            }
            coefficients_zero = coefficients_zero && proof.value() == ZeroProof::Zero;
        }
        if (!coefficients_zero) continue;
        auto rhs = classify_exact_zero(
            rref.at(row, coefficient_columns), context, operation);
        if (!rhs) return Result<bool>::failure(rhs.error());
        if (rhs.value() == ZeroProof::Unknown) {
            return Result<bool>::failure(
                CasErrc::Inconclusive,
                "linear-system right-hand side is not provably zero", operation);
        }
        if (rhs.value() == ZeroProof::NonZero) {
            return Result<bool>::success(true);
        }
    }
    return Result<bool>::success(false);
}

Result<ExactLinearSolution> rational_linear_solve(
    const ExactMatrixData& augmented,
    std::size_t coefficient_columns,
    ComputationContext& context,
    const std::string& operation) {
    auto budget = context.consume_steps(
        augmented.rows * augmented.cols, operation);
    if (!budget) {
        return Result<ExactLinearSolution>::failure(budget.error());
    }
    std::vector<Rational> matrix;
    matrix.reserve(augmented.entries.size());
    for (const auto& entry : augmented.entries) {
        auto value = exact_rational(entry);
        if (!value) {
            return Result<ExactLinearSolution>::failure(
                CasErrc::Inconclusive,
                "rational specialization received a symbolic entry",
                operation);
        }
        matrix.push_back(*value);
    }
    auto pivot_rows = reduce_rational_rows(
        augmented.rows, coefficient_columns, matrix, false);
    if (rational_system_inconsistent(augmented.rows, coefficient_columns, matrix)) {
        return Result<ExactLinearSolution>::success(InconsistentLinearSolution{});
    }
    return assemble_linear_solution(coefficient_columns, pivot_rows,
        [&](std::size_t row, std::size_t column, bool negate) {
            const auto& value = matrix[row * augmented.cols + column];
            return SymbolicExpr::number(negate ? -value : value);
        });
}

}

Result<ExactLinearSolution> solve_linear_exact(
    ExactMatrixData augmented,
    std::size_t coefficient_columns,
    ComputationContext& context,
    const std::string& operation) {
    auto valid = validate_matrix(augmented, operation);

    if (!valid) return Result<ExactLinearSolution>::failure(valid.error());
    if (augmented.cols != coefficient_columns + 1) {
        return Result<ExactLinearSolution>::failure(
            CasErrc::InvalidArgument,
            "linear solve requires one augmented right-hand side", operation);
    }
    bool all_rational = true;
    for (const auto& entry : augmented.entries) {
        all_rational = all_rational && exact_rational(entry).has_value();
    }
    if (all_rational) {
        return rational_linear_solve(
            augmented, coefficient_columns, context, operation);
    }
    auto eliminated = eliminate_exact(
        std::move(augmented), coefficient_columns,
        EliminationForm::ReducedRowEchelon, context, operation);
    if (!eliminated) return Result<ExactLinearSolution>::failure(eliminated.error());
    auto& rref = eliminated.value().matrix;

    auto inconsistent = inconsistent_symbolic_system(
        rref, coefficient_columns, context, operation);
    if (!inconsistent) return Result<ExactLinearSolution>::failure(inconsistent.error());
    if (inconsistent.value()) {
        return Result<ExactLinearSolution>::success(InconsistentLinearSolution{});
    }
    std::vector<std::optional<std::size_t>> pivot_row(coefficient_columns);
    for (std::size_t index = 0;
         index < eliminated.value().pivot_columns.size(); ++index) {
        const auto column = eliminated.value().pivot_columns[index];
        if (column < coefficient_columns) pivot_row[column] = index;
    }
    return assemble_linear_solution(coefficient_columns, pivot_row,
        [&](std::size_t row, std::size_t column, bool negate) {
            const auto& value = rref.at(row, column);
            return negate ? exact_negate(value) : value;
        });
}

Result<std::vector<std::vector<ExprPtr>>> nullspace_exact(
    ExactMatrixData coefficients,
    ComputationContext& context,
    const std::string& operation) {
    auto valid = validate_matrix(coefficients, operation);
    if (!valid) {
        return Result<std::vector<std::vector<ExprPtr>>>::failure(valid.error());
    }
    const std::size_t columns = coefficients.cols;
    auto eliminated = eliminate_exact(
        std::move(coefficients), columns,
        EliminationForm::ReducedRowEchelon, context, operation);
    if (!eliminated) {
        return Result<std::vector<std::vector<ExprPtr>>>::failure(
            eliminated.error());
    }
    std::set<std::size_t> pivots(
        eliminated.value().pivot_columns.begin(),
        eliminated.value().pivot_columns.end());
    std::vector<std::vector<ExprPtr>> basis;
    for (std::size_t free_column = 0; free_column < columns; ++free_column) {
        if (pivots.count(free_column) != 0) continue;
        std::vector<ExprPtr> vector(columns, SymbolicExpr::number(0));
        vector[free_column] = SymbolicExpr::number(1);
        for (std::size_t row = 0;
             row < eliminated.value().pivot_columns.size(); ++row) {
            vector[eliminated.value().pivot_columns[row]] = exact_negate(
                eliminated.value().matrix.at(row, free_column));
        }
        basis.push_back(std::move(vector));
    }
    return Result<std::vector<std::vector<ExprPtr>>>::success(
        std::move(basis));
}

}
