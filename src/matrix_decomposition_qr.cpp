#include "internal/matrix_decomposition_support.hpp"
#include "internal/exact_rational_matrix.hpp"

#include <algorithm>
#include <new>

namespace LMCAS {
using namespace detail::matrix_decomposition;
namespace {
Result<void> prove_exact_full_column_rank_qr_support(
    const std::shared_ptr<const MatrixNode>& matrix,
    ComputationContext& context,
    const std::string& operation) {
    auto values = exact_rational_matrix(matrix);
    if (!values) {
        return Result<void>::failure(
            CasErrc::Inconclusive,
            "QR input is outside the exact rational full-column-rank support domain",
            operation);
    }
    const size_t rows = values->size();
    const size_t cols = rows == 0 ? 0 : (*values)[0].size();
    if (rows < cols) {
        return Result<void>::failure(
            CasErrc::Inconclusive,
            "QR decomposition requires at least as many rows as columns in the supported domain",
            operation);
    }
    std::vector<std::vector<Rational>> gram(
        cols, std::vector<Rational>(cols, Rational(0)));
    for (size_t i = 0; i < cols; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            Rational sum(0);
            for (size_t row = 0; row < rows; ++row) {
                sum = sum + (*values)[row][i] * (*values)[row][j];
            }
            gram[i][j] = sum;
        }
    }
    for (size_t order = 1; order <= cols; ++order) {
        std::vector<Rational> leading;
        leading.reserve(order * order);
        for (size_t row = 0; row < order; ++row) {
            for (size_t col = 0; col < order; ++col) {
                leading.push_back(gram[row][col]);
            }
        }
        auto determinant = detail::rational_determinant_exact(
            order, std::move(leading), context, operation);
        if (!determinant) return Result<void>::failure(determinant.error());
        if (determinant.value() <= Rational(0)) {
            return Result<void>::failure(
                CasErrc::Inconclusive,
                "QR decomposition requires proven full column rank",
                operation);
        }
    }
    return Result<void>::success();
}
static bool qr_decomposition_impl(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& Q,
    std::shared_ptr<SymbolicExpr>& R) {
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(A));
    if (!mat) return false;
    size_t m = mat->rows();
    size_t n = mat->cols();

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> Q_grid(m, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> R_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A_cols(n, std::vector<std::shared_ptr<SymbolicExpr>>(m));
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < m; i++) A_cols[j][i] = LMCAS::detail::make_expression_ptr(mat->get(i, j));
    }

    for (size_t j = 0; j < n; j++) {
        std::vector<std::shared_ptr<SymbolicExpr>> u_j = A_cols[j];
        for (size_t i = 0; i < j; i++) {
            auto dot = SymbolicExpr::number(0);
            for (size_t k = 0; k < m; k++) dot = SymbolicExpr::add(dot, SymbolicExpr::multiply(Q_grid[k][i], A_cols[j][k]));
            R_grid[i][j] = dot;

            for (size_t k = 0; k < m; k++) {
                u_j[k] = SymbolicExpr::add(u_j[k], SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(R_grid[i][j], Q_grid[k][i])));
            }
        }

        auto norm_sq = SymbolicExpr::number(0);
        for (size_t k = 0; k < m; k++) norm_sq = SymbolicExpr::add(norm_sq, SymbolicExpr::multiply(u_j[k], u_j[k]));
        R_grid[j][j] = SymbolicExpr::sqrt(norm_sq);

        for (size_t k = 0; k < m; k++) Q_grid[k][j] = SymbolicExpr::divide(u_j[k], R_grid[j][j]);
    }

    Q = SymbolicExpr::matrix(Q_grid);
    R = SymbolicExpr::matrix(R_grid);
    return true;
}
}

QRDecompositionResult qr_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context)
{
    const std::string operation = "qr_decomposition";
    auto matrix = validate_decomposition_matrix(A, false, context, operation);
    if (!matrix) return QRDecompositionResult::failure(matrix.error());
    auto budget = context.consume_steps(matrix.value()->rows() * matrix.value()->cols() * 12 + 8,
                                        operation);
    if (!budget) return QRDecompositionResult::failure(budget.error());
    auto full_rank = prove_exact_full_column_rank_qr_support(
        matrix.value(), context, operation);
    if (!full_rank) return QRDecompositionResult::failure(full_rank.error());
    try {
        std::shared_ptr<SymbolicExpr> Q;
        std::shared_ptr<SymbolicExpr> R;
        if (!qr_decomposition_impl(A, Q, R)) {
            return QRDecompositionResult::failure(
                CasErrc::Inconclusive,
                "QR decomposition could not be constructed in the supported symbolic domain",
                operation);
        }
        auto q_check = validate_decomposition_output(Q, "Q", operation);
        if (!q_check) return QRDecompositionResult::failure(q_check.error());
        auto r_check = validate_decomposition_output(R, "R", operation);
        if (!r_check) return QRDecompositionResult::failure(r_check.error());
        return QRDecompositionResult::success(QRDecomposition{Q, R});
    } catch (const std::bad_alloc&) {
        return QRDecompositionResult::failure(CasErrc::ResourceLimit,
                                             "allocation failed while calculating QR decomposition",
                                             operation);
    } catch (const std::exception& ex) {
        return QRDecompositionResult::failure(CasErrc::InternalInvariant,
                                             ex.what(),
                                             operation);
    }
}

QRDecompositionResult qr_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A)
{
    ComputationContext context;
    return qr_decomposition_checked(A, context);
}

std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> gram_schmidt(
    const std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& vectors,
    bool normalize) {
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> basis;
    for (const auto& v : vectors) {
        std::vector<std::shared_ptr<SymbolicExpr>> u = v;
        for (const auto& b : basis) {
            auto vb = SymbolicExpr::number(0);
            auto bb = SymbolicExpr::number(0);
            for (size_t k = 0; k < u.size(); ++k) {
                vb = SymbolicExpr::add(vb, SymbolicExpr::multiply(v[k], b[k]));
                bb = SymbolicExpr::add(bb, SymbolicExpr::multiply(b[k], b[k]));
            }
            if (LMCAS::detail::node(bb) && LMCAS::detail::node(bb)->is_zero()) continue;
            auto coeff = SymbolicExpr::divide(vb, bb);
            for (size_t k = 0; k < u.size(); ++k) {
                u[k] = SymbolicExpr::add(u[k],
                    SymbolicExpr::multiply(SymbolicExpr::number(-1),
                        SymbolicExpr::multiply(coeff, b[k])))->simplify();
            }
        }
        auto norm_sq = SymbolicExpr::number(0);
        for (auto& x : u) norm_sq = SymbolicExpr::add(norm_sq, SymbolicExpr::multiply(x, x));
        norm_sq = norm_sq->simplify();
        if (LMCAS::detail::node(norm_sq) && LMCAS::detail::node(norm_sq)->is_zero()) continue;
        basis.push_back(u);
    }
    if (normalize) {
        for (auto& u : basis) {
            auto norm_sq = SymbolicExpr::number(0);
            for (auto& x : u) norm_sq = SymbolicExpr::add(norm_sq, SymbolicExpr::multiply(x, x));
            auto norm = SymbolicExpr::sqrt(norm_sq);
            for (auto& x : u) x = SymbolicExpr::divide(x, norm)->simplify();
        }
    }
    return basis;
}

}
