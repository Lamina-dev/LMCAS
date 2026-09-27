#include "internal/exact_rational_matrix.hpp"

#include <algorithm>
#include <limits>
#include <new>

namespace LMCAS::detail {
Result<Rational> rational_determinant_value_impl(
    std::size_t dimension,
    std::vector<Rational> values,
    ComputationContext& context,
    const std::string& operation) {
    if (dimension == 0 ||
        dimension > std::numeric_limits<std::size_t>::max() / dimension ||
        values.size() != dimension * dimension) {
        return Result<Rational>::failure(
            CasErrc::InvalidArgument,
            "rational determinant storage does not match its dimension",
            operation);
    }
    Rational sign(1);
    Rational previous(1);
    if (dimension == 1) return Result<Rational>::success(values[0]);
    for (std::size_t pivot_column = 0;
         pivot_column + 1 < dimension; ++pivot_column) {
        auto budget = context.consume_steps(1, operation);
        if (!budget) return Result<Rational>::failure(budget.error());
        std::size_t pivot_row = pivot_column;
        while (pivot_row < dimension &&
               values[pivot_row * dimension + pivot_column] == Rational(0)) {
            ++pivot_row;
        }
        if (pivot_row == dimension) {
            return Result<Rational>::success(Rational(0));
        }
        if (pivot_row != pivot_column) {
            for (std::size_t column = 0; column < dimension; ++column) {
                std::swap(values[pivot_row * dimension + column],
                          values[pivot_column * dimension + column]);
            }
            sign = -sign;
        }
        const Rational pivot =
            values[pivot_column * dimension + pivot_column];
        for (std::size_t row = pivot_column + 1; row < dimension; ++row) {
            for (std::size_t column = pivot_column + 1;
                 column < dimension; ++column) {
                values[row * dimension + column] =
                    (values[row * dimension + column] * pivot -
                     values[row * dimension + pivot_column] *
                         values[pivot_column * dimension + column]) /
                    previous;
            }
            values[row * dimension + pivot_column] = Rational(0);
        }
        previous = pivot;
    }
    return Result<Rational>::success(
        sign * values[(dimension - 1) * dimension + (dimension - 1)]);
}

Result<Rational> rational_determinant_exact(
    std::size_t dimension,
    std::vector<Rational> values,
    ComputationContext& context,
    const std::string& operation) {
    try {
        return rational_determinant_value_impl(
            dimension, std::move(values), context, operation);
    } catch (const std::bad_alloc&) {
        return Result<Rational>::failure(
            CasErrc::ResourceLimit,
            "rational determinant allocation failed",
            operation);
    } catch (const std::exception& error) {
        return Result<Rational>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

}
