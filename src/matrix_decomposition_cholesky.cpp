#include "internal/matrix_decomposition_support.hpp"
#include "internal/exact_rational_matrix.hpp"

#include <algorithm>
#include <new>

namespace LMCAS {
using namespace detail::matrix_decomposition;
namespace {
Result<void> prove_exact_spd(
    const std::shared_ptr<const MatrixNode>& matrix,
    ComputationContext& context,
    const std::string& operation) {
    auto values = exact_rational_matrix(matrix);
    if (!values) {
        return Result<void>::failure(
            CasErrc::Inconclusive,
            "Cholesky input is outside the exact rational SPD support domain",
            operation);
    }
    const size_t n = values->size();
    for (size_t row = 0; row < n; ++row) {
        for (size_t col = row + 1; col < n; ++col) {
            if ((*values)[row][col] != (*values)[col][row]) {
                return Result<void>::failure(
                    CasErrc::DomainError,
                    "Cholesky input must be symmetric",
                    operation);
            }
        }
    }
    for (size_t order = 1; order <= n; ++order) {
        std::vector<Rational> leading;
        leading.reserve(order * order);
        for (size_t row = 0; row < order; ++row) {
            for (size_t col = 0; col < order; ++col) {
                leading.push_back((*values)[row][col]);
            }
        }
        auto determinant = detail::rational_determinant_exact(
            order, std::move(leading), context, operation);
        if (!determinant) return Result<void>::failure(determinant.error());
        if (determinant.value() <= Rational(0)) {
            return Result<void>::failure(
                CasErrc::DomainError,
                "Cholesky input is not positive definite",
                operation);
        }
    }
    return Result<void>::success();
}
static bool cholesky_decomposition_impl(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& L) {
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(A));
    if (!mat || mat->rows() != mat->cols()) return false;
    size_t n = mat->rows();
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> L_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            auto sum = SymbolicExpr::number(0);
            for (size_t k = 0; k < j; k++) sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(L_grid[i][k], L_grid[j][k]));
            auto a_ij = LMCAS::detail::make_expression_ptr(mat->get(i, j));
            auto diff = SymbolicExpr::add(a_ij, SymbolicExpr::multiply(SymbolicExpr::number(-1), sum));
            if (i == j) L_grid[i][j] = SymbolicExpr::sqrt(diff);
            else L_grid[i][j] = SymbolicExpr::divide(diff, L_grid[j][j]);
        }
    }
    L = SymbolicExpr::matrix(L_grid);
    return true;
}
}

CholeskyDecompositionResult cholesky_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context)
{
    const std::string operation = "cholesky_decomposition";
    auto matrix = validate_decomposition_matrix(A, true, context, operation);
    if (!matrix) return CholeskyDecompositionResult::failure(matrix.error());
    auto budget = context.consume_steps(matrix.value()->rows() * matrix.value()->cols() * 10 + 8,
                                        operation);
    if (!budget) return CholeskyDecompositionResult::failure(budget.error());
    auto spd = prove_exact_spd(matrix.value(), context, operation);
    if (!spd) return CholeskyDecompositionResult::failure(spd.error());
    try {
        std::shared_ptr<SymbolicExpr> L;
        if (!cholesky_decomposition_impl(A, L)) {
            return CholeskyDecompositionResult::failure(
                CasErrc::Inconclusive,
                "Cholesky decomposition could not be constructed in the supported symbolic domain",
                operation);
        }
        auto l_check = validate_decomposition_output(L, "L", operation);
        if (!l_check) return CholeskyDecompositionResult::failure(l_check.error());
        return CholeskyDecompositionResult::success(CholeskyDecomposition{L});
    } catch (const std::bad_alloc&) {
        return CholeskyDecompositionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while calculating Cholesky decomposition",
            operation);
    } catch (const std::exception& ex) {
        return CholeskyDecompositionResult::failure(CasErrc::InternalInvariant,
                                                   ex.what(),
                                                   operation);
    }
}

CholeskyDecompositionResult cholesky_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A)
{
    ComputationContext context;
    return cholesky_decomposition_checked(A, context);
}

}
