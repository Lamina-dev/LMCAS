#include "internal/matrix_decomposition_support.hpp"

#include <algorithm>
#include <new>

namespace LMCAS {
using namespace detail::matrix_decomposition;
namespace {
using RationalRows = std::vector<std::vector<Rational>>;

void eliminate_lu_column(RationalRows& lower, RationalRows& upper,
                         std::size_t column) {
    const std::size_t n = upper.size();
            for (std::size_t row = column + 1; row < n; ++row) {
                const Rational multiplier =
                    upper[row][column] / upper[column][column];
                lower[row][column] = multiplier;
                for (std::size_t col = column; col < n; ++col) {
                    upper[row][col] =
                        upper[row][col] -
                        multiplier * upper[column][col];
                }
            }
}

Result<void> factor_lu(RationalRows& lower, RationalRows& upper,
                       RationalRows& permutation, const std::string& operation) {
    const std::size_t n = upper.size();
        for (std::size_t i = 0; i < n; ++i) {
            lower[i][i] = Rational(1);
            permutation[i][i] = Rational(1);
        }
        for (std::size_t column = 0; column < n; ++column) {
            std::size_t pivot = column;
            while (pivot < n && upper[pivot][column] == Rational(0)) ++pivot;
            if (pivot == n) {
                return Result<void>::failure(
                    CasErrc::DomainError,
                    "PLU decomposition requires a nonsingular matrix",
                    operation);
            }
            if (pivot != column) {
                std::swap(upper[pivot], upper[column]);
                std::swap(permutation[pivot], permutation[column]);
                for (std::size_t prior = 0; prior < column; ++prior) {
                    std::swap(lower[pivot][prior], lower[column][prior]);
                }
            }
            eliminate_lu_column(lower, upper, column);
        }
    return Result<void>::success();
}

std::shared_ptr<SymbolicExpr> rational_rows_expression(const RationalRows& input) {
    const std::size_t n = input.size();
            std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(
                n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
            for (std::size_t row = 0; row < n; ++row) {
                for (std::size_t col = 0; col < n; ++col) {
                    grid[row][col] = exact_number_expr(input[row][col]);
                }
            }
            return SymbolicExpr::matrix(grid);
}
}

LUDecompositionResult lu_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context)
{
    const std::string operation = "lu_decomposition";
    auto matrix = validate_decomposition_matrix(A, true, context, operation);
    if (!matrix) return LUDecompositionResult::failure(matrix.error());
    auto budget = context.consume_steps(matrix.value()->rows() * matrix.value()->cols() * 8 + 8,
                                        operation);
    if (!budget) return LUDecompositionResult::failure(budget.error());
    auto values = exact_rational_matrix(matrix.value());
    if (!values) {
        return LUDecompositionResult::failure(
            CasErrc::Inconclusive,
            "exact PLU requires a rational matrix or proved symbolic pivots",
            operation);
    }
    try {
        const std::size_t n = values->size();
        auto U_values = *values;
        std::vector<std::vector<Rational>> L_values(
            n, std::vector<Rational>(n, Rational(0)));
        std::vector<std::vector<Rational>> P_values(
            n, std::vector<Rational>(n, Rational(0)));
        auto factored = factor_lu(L_values, U_values, P_values, operation);
        if (!factored) return LUDecompositionResult::failure(factored.error());
        return LUDecompositionResult::success(LUDecomposition{
            rational_rows_expression(P_values),
            rational_rows_expression(L_values),
            rational_rows_expression(U_values)});
    } catch (const std::bad_alloc&) {
        return LUDecompositionResult::failure(CasErrc::ResourceLimit,
                                             "allocation failed while calculating LU decomposition",
                                             operation);
    } catch (const std::exception& ex) {
        return LUDecompositionResult::failure(CasErrc::InternalInvariant,
                                             ex.what(),
                                             operation);
    }
}
LUDecompositionResult lu_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A)
{
    ComputationContext context;
    return lu_decomposition_checked(A, context);
}

}
