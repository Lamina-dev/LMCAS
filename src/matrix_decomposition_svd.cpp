#include "internal/matrix_decomposition_support.hpp"

#include <algorithm>
#include <new>

namespace LMCAS {
using namespace detail::matrix_decomposition;
SVDDecompositionResult svd_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context)
{
    const std::string operation = "svd_decomposition";
    auto matrix = validate_decomposition_matrix(A, false, context, operation);
    if (!matrix) return SVDDecompositionResult::failure(matrix.error());
    auto budget = context.consume_steps(matrix.value()->rows() * matrix.value()->cols() * 32 + 32,
                                        operation);
    if (!budget) return SVDDecompositionResult::failure(budget.error());
    auto diagonal = exact_rectangular_diagonal_entries(matrix.value());
    if (diagonal) {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> signs(
            matrix.value()->rows(),
            std::vector<std::shared_ptr<SymbolicExpr>>(
                matrix.value()->rows(), SymbolicExpr::number(0)));
        std::vector<Rational> magnitudes = *diagonal;
        for (std::size_t row = 0; row < matrix.value()->rows(); ++row) {
            signs[row][row] = SymbolicExpr::number(1);
        }
        for (std::size_t index = 0; index < magnitudes.size(); ++index) {
            if (magnitudes[index] < Rational(0)) {
                signs[index][index] = SymbolicExpr::number(-1);
                magnitudes[index] = Rational(0) - magnitudes[index];
            }
        }
        auto U = SymbolicExpr::matrix(signs);
        auto S = rectangular_diagonal_expr(
            matrix.value()->rows(), matrix.value()->cols(), magnitudes);
        auto V = identity_matrix_expr(matrix.value()->cols());
        return SVDDecompositionResult::success(SVDDecomposition{U, S, V});
    }
    return SVDDecompositionResult::failure(
        CasErrc::Inconclusive,
        "exact non-diagonal SVD requires a complete proved singular basis",
        operation);
}

SVDDecompositionResult svd_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A)
{
    ComputationContext context;
    return svd_decomposition_checked(A, context);
}

}
