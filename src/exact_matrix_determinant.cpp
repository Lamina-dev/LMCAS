#include "internal/exact_matrix_support.hpp"
#include "internal/symbolic_ast.hpp"

#include <algorithm>
#include <limits>
#include <new>
#include <set>

namespace LMCAS::detail {
using namespace matrix_kernel;
namespace {
Result<ExprPtr> rational_determinant(
    const ExactMatrixData& input,
    ComputationContext& context,
    const std::string& operation) {
    std::vector<Rational> values;
    values.reserve(input.entries.size());
    for (const auto& entry : input.entries) {
        auto value = exact_rational(entry);
        if (!value) {
            return Result<ExprPtr>::failure(
                CasErrc::Inconclusive,
                "matrix is not entirely rational",
                operation);
        }
        values.push_back(*value);
    }
    auto determinant = rational_determinant_value_impl(
        input.rows, std::move(values), context, operation);
    if (!determinant) return Result<ExprPtr>::failure(determinant.error());
    return Result<ExprPtr>::success(
        SymbolicExpr::number(determinant.value()));
}

Result<ExprPtr> determinant_cofactor(
    const ExactMatrixData& matrix,
    ComputationContext& context,
    const std::string& operation) {
    auto budget = context.consume_steps(1, operation);
    if (!budget) return Result<ExprPtr>::failure(budget.error());
    if (matrix.rows == 1) {
        return Result<ExprPtr>::success(matrix.entries[0]);
    }
    auto sum = SymbolicExpr::number(0);
    for (std::size_t column = 0; column < matrix.cols; ++column) {
        if (auto coefficient = exact_rational(matrix.at(0, column));
            coefficient && *coefficient == Rational(0)) {
            continue;
        }
        ExactMatrixData minor{
            matrix.rows - 1, matrix.cols - 1, {}};
        minor.entries.reserve(minor.rows * minor.cols);
        for (std::size_t row = 1; row < matrix.rows; ++row) {
            for (std::size_t source = 0; source < matrix.cols; ++source) {
                if (source != column) {
                    minor.entries.push_back(matrix.at(row, source));
                }
            }
        }
        auto minor_determinant = determinant_cofactor(
            minor, context, operation);
        if (!minor_determinant) {
            return Result<ExprPtr>::failure(minor_determinant.error());
        }
        auto term = exact_multiply(
            matrix.at(0, column), minor_determinant.value());
        if (column % 2 != 0) term = exact_negate(term);
        sum = exact_add(sum, term);
    }
    return Result<ExprPtr>::success(simplify_expr(sum));
}

Result<ExprPtr> symbolic_determinant_value(
    const ExactMatrixData& input, ComputationContext& context,
    const std::string& operation) {
    ExactMatrixData matrix = input;
    const std::size_t n = matrix.rows;
    if (n == 1) return Result<ExprPtr>::success(matrix.entries[0]);
    auto previous = SymbolicExpr::number(1);
    std::size_t swaps = 0;
    for (std::size_t column = 0; column + 1 < n; ++column) {
        auto pivot = choose_pivot(matrix, column, column, context, operation);
        if (!pivot) {
            if (pivot.error().code == CasErrc::Inconclusive) {
                return determinant_cofactor(input, context, operation);
            }
            return Result<ExprPtr>::failure(pivot.error());
        }
        if (!pivot.value()) {
            return Result<ExprPtr>::success(SymbolicExpr::number(0));
        }
        if (*pivot.value() != column) {
            swap_rows(matrix, *pivot.value(), column);
            ++swaps;
        }
        auto pivot_value = matrix.at(column, column);
        for (std::size_t row = column + 1; row < n; ++row) {
            for (std::size_t target = column + 1; target < n; ++target) {
                auto numerator = exact_subtract(
                    exact_multiply(matrix.at(row, target), pivot_value),
                    exact_multiply(matrix.at(row, column), matrix.at(column, target)));
                matrix.at(row, target) = exact_divide(numerator, previous);
            }
            matrix.at(row, column) = SymbolicExpr::number(0);
        }
        previous = pivot_value;
    }
    auto determinant = matrix.at(n - 1, n - 1);
    if (swaps % 2 != 0) determinant = exact_negate(determinant);
    return Result<ExprPtr>::success(simplify_expr(determinant));
}

}

Result<ExprPtr> determinant_exact(
    const ExactMatrixData& input,
    ComputationContext& context,
    const std::string& operation) {
    auto valid = validate_matrix(input, operation);
    if (!valid) return Result<ExprPtr>::failure(valid.error());
    if (input.rows != input.cols) {
        return Result<ExprPtr>::failure(
            CasErrc::InvalidArgument,
            "determinant requires a square matrix", operation);
    }
    for (const auto& entry : input.entries) {
        if (contains_approximate_number(entry)) {
            return Result<ExprPtr>::failure(
                CasErrc::Inconclusive,
                "exact determinant rejects approximate entries", operation);
        }
    }
    bool all_rational = true;
    for (const auto& entry : input.entries) {
        all_rational = all_rational && exact_rational(entry).has_value();
    }
    if (all_rational) return rational_determinant(input, context, operation);

    try {
        return symbolic_determinant_value(input, context, operation);
    } catch (const std::bad_alloc&) {
        return Result<ExprPtr>::failure(
            CasErrc::ResourceLimit,
            "exact determinant allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<ExprPtr>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

}
