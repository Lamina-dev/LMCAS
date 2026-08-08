/**
 * @file matrix_decomposition.cpp
 * @brief 矩阵高级分解算法实现。
 */
#include "matrix_decomposition.hpp"
#include "numeric_evaluation.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <string>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <optional>

namespace lamina {

namespace {

Result<std::shared_ptr<const MatrixNode>> validate_decomposition_matrix(
    const std::shared_ptr<SymbolicExpr>& A,
    bool require_square,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return Result<std::shared_ptr<const MatrixNode>>::failure(step.error());
    if (!A || !lamina::detail::node(A)) {
        return Result<std::shared_ptr<const MatrixNode>>::failure(
            CasErrc::InvalidArgument,
            "matrix input cannot be null",
            operation);
    }
    auto matrix = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(A));
    if (!matrix) {
        return Result<std::shared_ptr<const MatrixNode>>::failure(
            CasErrc::InvalidArgument,
            "input must be a matrix",
            operation);
    }
    if (require_square && matrix->rows() != matrix->cols()) {
        return Result<std::shared_ptr<const MatrixNode>>::failure(
            CasErrc::InvalidArgument,
            "matrix must be square",
            operation);
    }
    return Result<std::shared_ptr<const MatrixNode>>::success(matrix);
}

Result<void> validate_decomposition_output(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& name,
    const std::string& operation)
{
    if (!expr || !lamina::detail::node(expr) ||
        !std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(expr))) {
        return Result<void>::failure(
            CasErrc::InternalInvariant,
            name + " output is not a matrix",
            operation);
    }
    return Result<void>::success();
}

std::optional<double> checked_finite_numeric(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !lamina::detail::node(expr)) return std::nullopt;
    ComputationContext context;
    auto evaluated = evaluate_numeric(*expr, NumericBindings{}, context);
    if (!evaluated || !evaluated.value().is_finite() ||
        !std::isfinite(evaluated.value().value)) {
        return std::nullopt;
    }
    return evaluated.value().value;
}

std::optional<Rational> exact_rational_expr(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return std::nullopt;
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        if (std::holds_alternative<BigInt>(number->value())) {
            return Rational(std::get<BigInt>(number->value()));
        }
        if (std::holds_alternative<Rational>(number->value())) {
            return std::get<Rational>(number->value());
        }
        return std::nullopt;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        Rational sum(0);
        for (const auto& operand : add->operands()) {
            auto value = exact_rational_expr(operand);
            if (!value) return std::nullopt;
            sum = sum + *value;
        }
        return sum;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        Rational product(1);
        for (const auto& operand : mul->operands()) {
            auto value = exact_rational_expr(operand);
            if (!value) return std::nullopt;
            product = product * *value;
        }
        return product;
    }
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (fn->type() == FunctionNode::FuncType::Abs && fn->arguments().size() == 1) {
            auto value = exact_rational_expr(fn->arguments()[0]);
            if (!value) return std::nullopt;
            return (*value < Rational(0)) ? (Rational(0) - *value) : *value;
        }
    }
    return std::nullopt;
}

std::optional<Rational> exact_rational_expr(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !lamina::detail::node(expr)) return std::nullopt;
    auto direct = exact_rational_expr(lamina::detail::node(expr));
    if (direct) return direct;
    auto simplified = expr->simplify();
    return simplified ? exact_rational_expr(lamina::detail::node(simplified)) : std::nullopt;
}

std::shared_ptr<SymbolicExpr> exact_number_expr(const Rational& value) {
    if (value.get_denominator() == BigInt(1)) {
        return lamina::detail::make_expression_ptr(
            lamina::detail::make_node<NumberNode>(
                std::variant<BigInt, Rational, lmmc_real_t>{
                    std::in_place_type<BigInt>, value.get_numerator()}));
    }
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(
        std::variant<BigInt, Rational, lmmc_real_t>{
            std::in_place_type<Rational>, value}));
}

std::optional<std::vector<std::vector<Rational>>> exact_rational_matrix(
    const std::shared_ptr<const MatrixNode>& matrix) {
    if (!matrix) return std::nullopt;
    std::vector<std::vector<Rational>> values(
        matrix->rows(), std::vector<Rational>(matrix->cols(), Rational(0)));
    for (size_t row = 0; row < matrix->rows(); ++row) {
        for (size_t col = 0; col < matrix->cols(); ++col) {
            auto value = exact_rational_expr(matrix->get(row, col));
            if (!value) return std::nullopt;
            values[row][col] = *value;
        }
    }
    return values;
}

Rational exact_det(std::vector<std::vector<Rational>> matrix) {
    const size_t n = matrix.size();
    Rational det(1);
    int sign = 1;
    for (size_t col = 0; col < n; ++col) {
        size_t pivot = col;
        while (pivot < n && matrix[pivot][col] == Rational(0)) {
            ++pivot;
        }
        if (pivot == n) return Rational(0);
        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
            sign = -sign;
        }
        const Rational pivot_value = matrix[col][col];
        det = det * pivot_value;
        for (size_t row = col + 1; row < n; ++row) {
            const Rational factor = matrix[row][col] / pivot_value;
            for (size_t k = col; k < n; ++k) {
                matrix[row][k] = matrix[row][k] - factor * matrix[col][k];
            }
        }
    }
    return sign < 0 ? (Rational(0) - det) : det;
}

Result<void> prove_exact_spd(
    const std::shared_ptr<const MatrixNode>& matrix,
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
        std::vector<std::vector<Rational>> leading(
            order, std::vector<Rational>(order, Rational(0)));
        for (size_t row = 0; row < order; ++row) {
            for (size_t col = 0; col < order; ++col) {
                leading[row][col] = (*values)[row][col];
            }
        }
        if (exact_det(std::move(leading)) <= Rational(0)) {
            return Result<void>::failure(
                CasErrc::DomainError,
                "Cholesky input is not positive definite",
                operation);
        }
    }
    return Result<void>::success();
}

Result<void> prove_exact_no_pivot_lu_support(
    const std::shared_ptr<const MatrixNode>& matrix,
    const std::string& operation) {
    auto values = exact_rational_matrix(matrix);
    if (!values) {
        return Result<void>::failure(
            CasErrc::Inconclusive,
            "LU input is outside the exact rational no-pivot support domain",
            operation);
    }
    const size_t n = values->size();
    for (size_t order = 1; order <= n; ++order) {
        std::vector<std::vector<Rational>> leading(
            order, std::vector<Rational>(order, Rational(0)));
        for (size_t row = 0; row < order; ++row) {
            for (size_t col = 0; col < order; ++col) {
                leading[row][col] = (*values)[row][col];
            }
        }
        if (exact_det(std::move(leading)) == Rational(0)) {
            return Result<void>::failure(
                CasErrc::Inconclusive,
                "LU decomposition requires nonzero no-pivot leading minors",
                operation);
        }
    }
    return Result<void>::success();
}

Result<void> prove_exact_full_column_rank_qr_support(
    const std::shared_ptr<const MatrixNode>& matrix,
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
        std::vector<std::vector<Rational>> leading(
            order, std::vector<Rational>(order, Rational(0)));
        for (size_t row = 0; row < order; ++row) {
            for (size_t col = 0; col < order; ++col) {
                leading[row][col] = gram[row][col];
            }
        }
        if (exact_det(std::move(leading)) <= Rational(0)) {
            return Result<void>::failure(
                CasErrc::Inconclusive,
                "QR decomposition requires proven full column rank",
                operation);
        }
    }
    return Result<void>::success();
}

std::shared_ptr<SymbolicExpr> identity_matrix_expr(size_t n) {
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(
        n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < n; ++i) {
        grid[i][i] = SymbolicExpr::number(1);
    }
    return SymbolicExpr::matrix(grid);
}

std::optional<std::vector<Rational>> exact_nonnegative_diagonal_entries(
    const std::shared_ptr<const MatrixNode>& matrix) {
    auto values = exact_rational_matrix(matrix);
    if (!values) return std::nullopt;
    const size_t rows = values->size();
    const size_t cols = rows == 0 ? 0 : (*values)[0].size();
    const size_t diag_count = std::min(rows, cols);
    std::vector<Rational> diagonal(diag_count, Rational(0));
    for (size_t row = 0; row < rows; ++row) {
        for (size_t col = 0; col < cols; ++col) {
            const bool is_diag = row == col && row < diag_count;
            if (!is_diag && (*values)[row][col] != Rational(0)) {
                return std::nullopt;
            }
            if (is_diag) {
                if ((*values)[row][col] < Rational(0)) {
                    return std::nullopt;
                }
                diagonal[row] = (*values)[row][col];
            }
        }
    }
    return diagonal;
}

std::optional<std::vector<Rational>> exact_diagonal_entries(
    const std::shared_ptr<const MatrixNode>& matrix) {
    auto values = exact_rational_matrix(matrix);
    if (!values) return std::nullopt;
    const size_t rows = values->size();
    const size_t cols = rows == 0 ? 0 : (*values)[0].size();
    if (rows != cols) return std::nullopt;
    std::vector<Rational> diagonal(rows, Rational(0));
    for (size_t row = 0; row < rows; ++row) {
        for (size_t col = 0; col < cols; ++col) {
            if (row == col) {
                diagonal[row] = (*values)[row][col];
            } else if ((*values)[row][col] != Rational(0)) {
                return std::nullopt;
            }
        }
    }
    return diagonal;
}

std::shared_ptr<SymbolicExpr> rectangular_diagonal_expr(
    size_t rows,
    size_t cols,
    const std::vector<Rational>& diagonal) {
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(
        rows, std::vector<std::shared_ptr<SymbolicExpr>>(cols, SymbolicExpr::number(0)));
    for (size_t i = 0; i < diagonal.size(); ++i) {
        grid[i][i] = exact_number_expr(diagonal[i]);
    }
    return SymbolicExpr::matrix(grid);
}

} // namespace

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
    auto no_pivot = prove_exact_no_pivot_lu_support(matrix.value(), operation);
    if (!no_pivot) return LUDecompositionResult::failure(no_pivot.error());
    try {
        std::shared_ptr<SymbolicExpr> L;
        std::shared_ptr<SymbolicExpr> U;
        if (!lu_decomposition(A, L, U)) {
            return LUDecompositionResult::failure(
                CasErrc::Inconclusive,
                "LU decomposition is outside the current no-pivot symbolic support domain",
                operation);
        }
        auto l_check = validate_decomposition_output(L, "L", operation);
        if (!l_check) return LUDecompositionResult::failure(l_check.error());
        auto u_check = validate_decomposition_output(U, "U", operation);
        if (!u_check) return LUDecompositionResult::failure(u_check.error());
        return LUDecompositionResult::success(LUDecomposition{L, U});
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

bool lu_decomposition(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& L,
    std::shared_ptr<SymbolicExpr>& U) {
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(A));
    if (!mat || mat->rows() != mat->cols()) return false;
    size_t n = mat->rows();
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> L_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> U_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    
    for (size_t i = 0; i < n; i++) {
        for (size_t k = i; k < n; k++) {
            auto sum = SymbolicExpr::number(0);
            for (size_t j = 0; j < i; j++) sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(L_grid[i][j], U_grid[j][k]));
            auto a_ik = lamina::detail::make_expression_ptr(mat->get(i, k));
            U_grid[i][k] = SymbolicExpr::add(a_ik, SymbolicExpr::multiply(SymbolicExpr::number(-1), sum));
        }

        for (size_t k = i; k < n; k++) {
            if (i == k) {
                L_grid[i][i] = SymbolicExpr::number(1);
            } else {
                auto sum = SymbolicExpr::number(0);
                for (size_t j = 0; j < i; j++) sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(L_grid[k][j], U_grid[j][i]));
                auto a_ki = lamina::detail::make_expression_ptr(mat->get(k, i));
                L_grid[k][i] = SymbolicExpr::divide(SymbolicExpr::add(a_ki, SymbolicExpr::multiply(SymbolicExpr::number(-1), sum)), U_grid[i][i]);
            }
        }
    }
    
    L = SymbolicExpr::matrix(L_grid);
    U = SymbolicExpr::matrix(U_grid);
    return true;
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
    auto full_rank = prove_exact_full_column_rank_qr_support(matrix.value(), operation);
    if (!full_rank) return QRDecompositionResult::failure(full_rank.error());
    try {
        std::shared_ptr<SymbolicExpr> Q;
        std::shared_ptr<SymbolicExpr> R;
        if (!qr_decomposition(A, Q, R)) {
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

bool qr_decomposition(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& Q,
    std::shared_ptr<SymbolicExpr>& R) {
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(A));
    if (!mat) return false;
    size_t m = mat->rows();
    size_t n = mat->cols();
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> Q_grid(m, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> R_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A_cols(n, std::vector<std::shared_ptr<SymbolicExpr>>(m));
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < m; i++) A_cols[j][i] = lamina::detail::make_expression_ptr(mat->get(i, j));
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
    auto spd = prove_exact_spd(matrix.value(), operation);
    if (!spd) return CholeskyDecompositionResult::failure(spd.error());
    try {
        std::shared_ptr<SymbolicExpr> L;
        if (!cholesky_decomposition(A, L)) {
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

bool cholesky_decomposition(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& L) {
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(A));
    if (!mat || mat->rows() != mat->cols()) return false;
    size_t n = mat->rows();
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> L_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            auto sum = SymbolicExpr::number(0);
            for (size_t k = 0; k < j; k++) sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(L_grid[i][k], L_grid[j][k]));
            auto a_ij = lamina::detail::make_expression_ptr(mat->get(i, j));
            auto diff = SymbolicExpr::add(a_ij, SymbolicExpr::multiply(SymbolicExpr::number(-1), sum));
            if (i == j) L_grid[i][j] = SymbolicExpr::sqrt(diff);
            else L_grid[i][j] = SymbolicExpr::divide(diff, L_grid[j][j]);
        }
    }
    L = SymbolicExpr::matrix(L_grid);
    return true;
}

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
    auto diagonal = exact_nonnegative_diagonal_entries(matrix.value());
    if (!diagonal) {
        return SVDDecompositionResult::failure(
            CasErrc::Inconclusive,
            "SVD checked support is currently limited to exact rational nonnegative diagonal matrices",
            operation);
    }
    auto U = identity_matrix_expr(matrix.value()->rows());
    auto S = rectangular_diagonal_expr(matrix.value()->rows(), matrix.value()->cols(), *diagonal);
    auto V = identity_matrix_expr(matrix.value()->cols());
    return SVDDecompositionResult::success(SVDDecomposition{U, S, V});
}

SVDDecompositionResult svd_decomposition_checked(
    const std::shared_ptr<SymbolicExpr>& A)
{
    ComputationContext context;
    return svd_decomposition_checked(A, context);
}

static std::vector<std::shared_ptr<SymbolicExpr>> normalize_vec(std::vector<std::shared_ptr<SymbolicExpr>> v) {
    auto norm_sq = SymbolicExpr::number(0);
    for (auto& x : v) norm_sq = SymbolicExpr::add(norm_sq, SymbolicExpr::multiply(x, x));
    auto norm = SymbolicExpr::sqrt(norm_sq);
    std::vector<std::shared_ptr<SymbolicExpr>> res;
    for (auto& x : v) res.push_back(SymbolicExpr::divide(x, norm));
    return res;
}

bool svd_decomposition(
    const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& U,
    std::shared_ptr<SymbolicExpr>& S,
    std::shared_ptr<SymbolicExpr>& V) {
    
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(A));
    if (!mat) return false;
    
    auto AT = SymbolicExpr::transpose(A);
    auto ATA = SymbolicExpr::multiply(AT, A);
    
    auto eigen_V = SymbolicExpr::eigenvectors(ATA);
    if (eigen_V.empty()) return false; // Fallback fails if matrix > 4x4 or no symbolic roots
    
    size_t n = mat->cols();
    size_t m = mat->rows();
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> V_cols;
    std::vector<std::shared_ptr<SymbolicExpr>> sigmas;
    
    for (auto& pair : eigen_V) {
        auto sigma = SymbolicExpr::sqrt(pair.first);
        for (auto& vec_mat : pair.second) {
            auto v_mat_node = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(vec_mat));
            std::vector<std::shared_ptr<SymbolicExpr>> v_col;
            for (size_t i = 0; i < n; i++) v_col.push_back(lamina::detail::make_expression_ptr(v_mat_node->get(i, 0)));
            V_cols.push_back(normalize_vec(v_col));
            sigmas.push_back(sigma);
        }
    }
    
    if (V_cols.size() != n) return false; // Incomplete basis
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> U_cols;
    for (size_t i = 0; i < n; i++) {
        if (!lamina::detail::node(sigmas[i])->is_zero()) {
            /// u_i = A v_i / sigma_i
            std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> v_col_mat(n, std::vector<std::shared_ptr<SymbolicExpr>>(1));
            for (size_t j = 0; j < n; j++) v_col_mat[j][0] = V_cols[i][j];
            auto Avi = SymbolicExpr::multiply(A, SymbolicExpr::matrix(v_col_mat));
            auto Avi_node = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(Avi));
            std::vector<std::shared_ptr<SymbolicExpr>> u_col;
            for (size_t j = 0; j < m; j++) {
                u_col.push_back(SymbolicExpr::divide(lamina::detail::make_expression_ptr(Avi_node->get(j, 0)), sigmas[i]));
            }
            U_cols.push_back(u_col);
        }
    }
    
    /// Complete U basis with standard basis via Gram-Schmidt if needed (simplified for stub replacement)
    for (size_t i = 0; i < m && U_cols.size() < m; i++) {
        std::vector<std::shared_ptr<SymbolicExpr>> ei(m, SymbolicExpr::number(0));
        ei[i] = SymbolicExpr::number(1);
        auto ui = ei;
        for (auto& uj : U_cols) {
            auto dot = SymbolicExpr::number(0);
            for (size_t k = 0; k < m; k++) dot = SymbolicExpr::add(dot, SymbolicExpr::multiply(uj[k], ei[k]));
            for (size_t k = 0; k < m; k++) ui[k] = SymbolicExpr::add(ui[k], SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(dot, uj[k])));
        }
        auto norm_sq = SymbolicExpr::number(0);
        for (auto& x : ui) norm_sq = SymbolicExpr::add(norm_sq, SymbolicExpr::multiply(x, x));
        if (!lamina::detail::node(norm_sq)->is_zero()) {
            U_cols.push_back(normalize_vec(ui));
        }
    }
    
    if (U_cols.size() != m) return false;
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> U_grid(m, std::vector<std::shared_ptr<SymbolicExpr>>(m));
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) U_grid[i][j] = U_cols[j][i];
    }
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> V_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) V_grid[i][j] = V_cols[j][i];
    }
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> S_grid(m, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < std::min(m, n); i++) S_grid[i][i] = sigmas[i];
    
    U = SymbolicExpr::matrix(U_grid);
    S = SymbolicExpr::matrix(S_grid);
    V = SymbolicExpr::matrix(V_grid);
    return true;
}

std::shared_ptr<SymbolicExpr> matrix_exp(
    const std::shared_ptr<SymbolicExpr>& A) {
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(A));
    if (!mat || mat->rows() != mat->cols()) return SymbolicExpr::exp(A);
    size_t n = mat->rows();

    bool is_zero_matrix = true;
    for (size_t r = 0; r < n && is_zero_matrix; ++r) {
        for (size_t c = 0; c < n; ++c) {
            auto entry = lamina::detail::make_expression_ptr(mat->get(r, c))->simplify();
            if (!lamina::detail::node(entry)->is_zero()) {
                is_zero_matrix = false;
                break;
            }
        }
    }
    if (is_zero_matrix) {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> identity(
            n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
        for (size_t i = 0; i < n; ++i) {
            identity[i][i] = SymbolicExpr::number(1);
        }
        return SymbolicExpr::matrix(identity);
    }
    
    auto eigen_V = SymbolicExpr::eigenvectors(A);
    if (eigen_V.empty()) return SymbolicExpr::exp(A); // Fallback to AST node
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> P_cols;
    std::vector<std::shared_ptr<SymbolicExpr>> evals;
    for (auto& pair : eigen_V) {
        for (auto& vec : pair.second) {
            auto v_mat_node = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(vec));
            std::vector<std::shared_ptr<SymbolicExpr>> v_col;
            for (size_t i = 0; i < n; i++) v_col.push_back(lamina::detail::make_expression_ptr(v_mat_node->get(i, 0)));
            P_cols.push_back(v_col);
            evals.push_back(pair.first);
        }
    }
    
    if (P_cols.size() != n) return SymbolicExpr::exp(A); // Not diagonalizable
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> P_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) P_grid[i][j] = P_cols[j][i];
    }
    auto P = SymbolicExpr::matrix(P_grid);
    auto P_inv = SymbolicExpr::inverse(P);
    if (!P_inv) return SymbolicExpr::exp(A);
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> expD_grid(n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < n; i++) expD_grid[i][i] = SymbolicExpr::exp(evals[i]);
    auto expD = SymbolicExpr::matrix(expD_grid);
    
    return SymbolicExpr::multiply(P, SymbolicExpr::multiply(expD, P_inv));
}


std::shared_ptr<SymbolicExpr> matrix_trace(const std::shared_ptr<SymbolicExpr>& A) {
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(A));
    if (!mat || mat->rows() != mat->cols()) return nullptr;
    auto sum = SymbolicExpr::number(0);
    for (size_t i = 0; i < mat->rows(); ++i) {
        sum = SymbolicExpr::add(sum, lamina::detail::make_expression_ptr(mat->get(i, i)));
    }
    return sum->simplify();
}

int matrix_rank(const std::shared_ptr<SymbolicExpr>& A) {
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(A));
    if (!mat) return -1;
    size_t m = mat->rows(), n = mat->cols();
    if (m == 0 || n == 0) return 0;

    /// 构造可变副本（数值化为 double，符号项无法判零时视为非零主元）
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> g(m,
        std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            g[i][j] = lamina::detail::make_expression_ptr(mat->get(i, j))->simplify();

    auto is_zero_expr = [](const std::shared_ptr<SymbolicExpr>& e) {
        return e && lamina::detail::node(e) && lamina::detail::node(e)->is_zero();
    };

    int rank = 0;
    size_t row = 0;
    for (size_t col = 0; col < n && row < m; ++col) {
        /// 寻找主元
        size_t pivot = row;
        while (pivot < m && is_zero_expr(g[pivot][col])) ++pivot;
        if (pivot == m) continue;
        std::swap(g[row], g[pivot]);

        auto pval = g[row][col];
        /// 消元
        for (size_t i = 0; i < m; ++i) {
            if (i == row || is_zero_expr(g[i][col])) continue;
            auto factor = SymbolicExpr::divide(g[i][col], pval);
            for (size_t j = col; j < n; ++j) {
                auto sub = SymbolicExpr::multiply(factor, g[row][j]);
                g[i][j] = SymbolicExpr::add(g[i][j],
                    SymbolicExpr::multiply(SymbolicExpr::number(-1), sub))->simplify();
            }
        }
        ++row;
        ++rank;
    }
    return rank;
}


std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> gram_schmidt(
    const std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& vectors,
    bool normalize) {
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> basis;
    for (const auto& v : vectors) {
        std::vector<std::shared_ptr<SymbolicExpr>> u = v;
        for (const auto& b : basis) {
            /// 投影系数 <v,b>/<b,b>
            auto vb = SymbolicExpr::number(0);
            auto bb = SymbolicExpr::number(0);
            for (size_t k = 0; k < u.size(); ++k) {
                vb = SymbolicExpr::add(vb, SymbolicExpr::multiply(v[k], b[k]));
                bb = SymbolicExpr::add(bb, SymbolicExpr::multiply(b[k], b[k]));
            }
            if (lamina::detail::node(bb) && lamina::detail::node(bb)->is_zero()) continue;
            auto coeff = SymbolicExpr::divide(vb, bb);
            for (size_t k = 0; k < u.size(); ++k) {
                u[k] = SymbolicExpr::add(u[k],
                    SymbolicExpr::multiply(SymbolicExpr::number(-1),
                        SymbolicExpr::multiply(coeff, b[k])))->simplify();
            }
        }
        /// 检查 u 是否为零向量（线性相关）
        auto norm_sq = SymbolicExpr::number(0);
        for (auto& x : u) norm_sq = SymbolicExpr::add(norm_sq, SymbolicExpr::multiply(x, x));
        norm_sq = norm_sq->simplify();
        if (lamina::detail::node(norm_sq) && lamina::detail::node(norm_sq)->is_zero()) continue;
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


std::shared_ptr<SymbolicExpr> matrix_log(const std::shared_ptr<SymbolicExpr>& A) {
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(A));
    if (!mat || mat->rows() != mat->cols()) return nullptr;
    size_t n = mat->rows();

    auto eigen_V = SymbolicExpr::eigenvectors(A);
    if (eigen_V.empty()) return nullptr;

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> P_cols;
    std::vector<std::shared_ptr<SymbolicExpr>> evals;
    for (auto& pr : eigen_V) {
        for (auto& vec : pr.second) {
            auto vnode = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(vec));
            if (!vnode) return nullptr;
            std::vector<std::shared_ptr<SymbolicExpr>> col;
            for (size_t i = 0; i < n; ++i) col.push_back(lamina::detail::make_expression_ptr(vnode->get(i, 0)));
            P_cols.push_back(col);
            evals.push_back(pr.first);
        }
    }
    if (P_cols.size() != n) return nullptr;

    /// 检查特征值为正（实数对数存在性）
    for (auto& ev : evals) {
        auto se = ev->simplify();
        if (se->is_number()) {
            auto d = checked_finite_numeric(se);
            if (!d || *d <= 0) return nullptr;
        }
    }

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> P_grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) P_grid[i][j] = P_cols[j][i];
    auto P = SymbolicExpr::matrix(P_grid);
    auto P_inv = SymbolicExpr::inverse(P);
    if (!P_inv) return nullptr;

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> logD(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < n; ++i) logD[i][i] = SymbolicExpr::ln(evals[i]);
    auto logDm = SymbolicExpr::matrix(logD);

    return SymbolicExpr::multiply(P, SymbolicExpr::multiply(logDm, P_inv));
}


std::shared_ptr<SymbolicExpr> kronecker(const std::shared_ptr<SymbolicExpr>& A,
    const std::shared_ptr<SymbolicExpr>& B) {
    auto a = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(A));
    auto b = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(B));
    if (!a || !b) return nullptr;
    size_t m = a->rows(), n = a->cols(), p = b->rows(), q = b->cols();
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(m * p,
        std::vector<std::shared_ptr<SymbolicExpr>>(n * q));
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j) {
            auto aij = lamina::detail::make_expression_ptr(a->get(i, j));
            for (size_t k = 0; k < p; ++k)
                for (size_t l = 0; l < q; ++l) {
                    auto bkl = lamina::detail::make_expression_ptr(b->get(k, l));
                    grid[i * p + k][j * q + l] = SymbolicExpr::multiply(aij, bkl)->simplify();
                }
        }
    return SymbolicExpr::matrix(grid);
}


std::shared_ptr<SymbolicExpr> matrix_norm(const std::shared_ptr<SymbolicExpr>& A,
    const std::string& type) {
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(A));
    if (!mat) return nullptr;
    size_t m = mat->rows(), n = mat->cols();

    auto abs_of = [](const std::shared_ptr<SymbolicExpr>& e) {
        auto node = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(e)});
        return lamina::detail::make_expression_ptr(node);
    };

    if (type == "frobenius") {
        auto sum = SymbolicExpr::number(0);
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j) {
                auto e = lamina::detail::make_expression_ptr(mat->get(i, j));
                sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(e, e));
            }
        return SymbolicExpr::sqrt(sum)->simplify();
    }
    if (type == "1") {
        /// 最大列和
        std::shared_ptr<SymbolicExpr> best = nullptr;
        std::optional<Rational> best_exact;
        double best_val = -1;
        for (size_t j = 0; j < n; ++j) {
            auto colsum = SymbolicExpr::number(0);
            std::optional<Rational> exact_sum = Rational(0);
            for (size_t i = 0; i < m; ++i) {
                auto entry = lamina::detail::make_expression_ptr(mat->get(i, j));
                colsum = SymbolicExpr::add(colsum, abs_of(entry));
                if (exact_sum) {
                    auto exact_entry = exact_rational_expr(entry);
                    if (exact_entry) {
                        *exact_sum = *exact_sum +
                            ((*exact_entry < Rational(0)) ? (Rational(0) - *exact_entry) : *exact_entry);
                    } else {
                        exact_sum.reset();
                    }
                }
            }
            colsum = colsum->simplify();
            auto exact = exact_sum ? exact_sum : exact_rational_expr(colsum);
            if (exact && (!best || (best_exact && *exact > *best_exact))) {
                best = exact_number_expr(*exact);
                best_exact = *exact;
                best_val = -1;
                continue;
            }
            if (exact && !best_exact && best) {
                auto maybe_v = checked_finite_numeric(colsum);
                if (maybe_v && *maybe_v > best_val) {
                    best = exact_number_expr(*exact);
                    best_val = *maybe_v;
                }
                continue;
            }
            if (!exact && !best_exact) {
                auto maybe_v = checked_finite_numeric(colsum);
                double v = maybe_v ? *maybe_v : -std::numeric_limits<double>::infinity();
                if (!best || v > best_val) { best = colsum; best_val = v; }
            }
        }
        return best;
    }
    if (type == "inf") {
        /// 最大行和
        std::shared_ptr<SymbolicExpr> best = nullptr;
        std::optional<Rational> best_exact;
        double best_val = -1;
        for (size_t i = 0; i < m; ++i) {
            auto rowsum = SymbolicExpr::number(0);
            std::optional<Rational> exact_sum = Rational(0);
            for (size_t j = 0; j < n; ++j) {
                auto entry = lamina::detail::make_expression_ptr(mat->get(i, j));
                rowsum = SymbolicExpr::add(rowsum, abs_of(entry));
                if (exact_sum) {
                    auto exact_entry = exact_rational_expr(entry);
                    if (exact_entry) {
                        *exact_sum = *exact_sum +
                            ((*exact_entry < Rational(0)) ? (Rational(0) - *exact_entry) : *exact_entry);
                    } else {
                        exact_sum.reset();
                    }
                }
            }
            rowsum = rowsum->simplify();
            auto exact = exact_sum ? exact_sum : exact_rational_expr(rowsum);
            if (exact && (!best || (best_exact && *exact > *best_exact))) {
                best = exact_number_expr(*exact);
                best_exact = *exact;
                best_val = -1;
                continue;
            }
            if (exact && !best_exact && best) {
                auto maybe_v = checked_finite_numeric(rowsum);
                if (maybe_v && *maybe_v > best_val) {
                    best = exact_number_expr(*exact);
                    best_val = *maybe_v;
                }
                continue;
            }
            if (!exact && !best_exact) {
                auto maybe_v = checked_finite_numeric(rowsum);
                double v = maybe_v ? *maybe_v : -std::numeric_limits<double>::infinity();
                if (!best || v > best_val) { best = rowsum; best_val = v; }
            }
        }
        return best;
    }
    return nullptr;
}


std::shared_ptr<SymbolicExpr> quadratic_form_matrix(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::vector<std::string>& vars) {
    if (!expr) return nullptr;
    size_t n = vars.size();
    auto e = expr->expand();
    if (!e) e = expr;

    /// A[i][j] = (1/2) ∂²(expr)/∂xᵢ∂xⱼ （对称矩阵，对角为 ∂²/2 即系数本身的一半*2）
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < n; ++i) {
        auto di = e->differentiate(vars[i]);
        for (size_t j = 0; j < n; ++j) {
            auto dij = di->differentiate(vars[j]);
            /// 对二次型，二阶偏导为常数；A_ij = (1/2)·∂²/∂xi∂xj
            grid[i][j] = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)), dij)->simplify();
        }
    }
    return SymbolicExpr::matrix(grid);
}

std::string classify_quadratic_form(const std::shared_ptr<SymbolicExpr>& A) {
    /// 使用特征值列表（来自特征多项式求根），避免脆弱的特征向量求解。
    auto evals_expr = SymbolicExpr::eigenvalues(A);
    if (!evals_expr) return "unknown";
    auto mat_node = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(evals_expr));
    if (!mat_node || mat_node->cols() == 0) return "unknown";

    std::vector<std::shared_ptr<SymbolicExpr>> evals;
    for (size_t i = 0; i < mat_node->cols(); ++i) {
        evals.push_back(lamina::detail::make_expression_ptr(mat_node->get(0, i)));
    }
    if (evals.empty()) return "unknown";

    bool all_pos = true, all_neg = true, any_zero = false;
    bool has_pos = false, has_neg = false;
    const double eps = 1e-9;
    for (auto& ev : evals) {
        auto se = ev->simplify();
        auto maybe_d = checked_finite_numeric(se);
        if (!maybe_d) return "unknown";
        double d = *maybe_d;
        if (d > eps) { all_neg = false; has_pos = true; }
        else if (d < -eps) { all_pos = false; has_neg = true; }
        else { any_zero = true; all_pos = false; all_neg = false; }
    }
    if (all_pos) return "positive_definite";
    if (all_neg) return "negative_definite";
    if (has_pos && has_neg) return "indefinite";
    if (has_pos && any_zero) return "positive_semidefinite";
    if (has_neg && any_zero) return "negative_semidefinite";
    return "indefinite";
}


JordanDecompositionResult jordan_form_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context)
{
    const std::string operation = "jordan_form";
    auto matrix = validate_decomposition_matrix(A, true, context, operation);
    if (!matrix) return JordanDecompositionResult::failure(matrix.error());
    auto budget = context.consume_steps(matrix.value()->rows() * matrix.value()->cols() * 48 + 32,
                                        operation);
    if (!budget) return JordanDecompositionResult::failure(budget.error());
    try {
        auto diagonal = exact_diagonal_entries(matrix.value());
        if (!diagonal) {
            return JordanDecompositionResult::failure(
                CasErrc::Inconclusive,
                "Jordan checked support is currently limited to exact rational diagonal matrices",
                operation);
        }
        auto J = rectangular_diagonal_expr(matrix.value()->rows(), matrix.value()->cols(), *diagonal);
        auto P = identity_matrix_expr(matrix.value()->rows());
        auto j_check = validate_decomposition_output(J, "J", operation);
        if (!j_check) return JordanDecompositionResult::failure(j_check.error());
        auto p_check = validate_decomposition_output(P, "P", operation);
        if (!p_check) return JordanDecompositionResult::failure(p_check.error());
        return JordanDecompositionResult::success(JordanDecomposition{J, P});
    } catch (const std::bad_alloc&) {
        return JordanDecompositionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while calculating Jordan form",
            operation);
    } catch (const std::exception& ex) {
        return JordanDecompositionResult::failure(CasErrc::InternalInvariant,
                                                 ex.what(),
                                                 operation);
    }
}

JordanDecompositionResult jordan_form_checked(
    const std::shared_ptr<SymbolicExpr>& A)
{
    ComputationContext context;
    return jordan_form_checked(A, context);
}

bool jordan_form(const std::shared_ptr<SymbolicExpr>& A,
    std::shared_ptr<SymbolicExpr>& J, std::shared_ptr<SymbolicExpr>& P) {
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(A));
    if (!mat || mat->rows() != mat->cols()) return false;
    size_t n = mat->rows();

    bool is_diagonal = true;
    for (size_t r = 0; r < n && is_diagonal; ++r) {
        for (size_t c = 0; c < n; ++c) {
            if (r == c) continue;
            auto entry = lamina::detail::make_expression_ptr(mat->get(r, c))->simplify();
            if (!lamina::detail::node(entry)->is_zero()) {
                is_diagonal = false;
                break;
            }
        }
    }
    if (is_diagonal) {
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> J_grid(
            n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> P_grid(
            n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
        for (size_t i = 0; i < n; ++i) {
            J_grid[i][i] = lamina::detail::make_expression_ptr(mat->get(i, i))->simplify();
            P_grid[i][i] = SymbolicExpr::number(1);
        }
        J = SymbolicExpr::matrix(J_grid);
        P = SymbolicExpr::matrix(P_grid);
        return true;
    }

    auto eigen_V = SymbolicExpr::eigenvectors(A);
    if (eigen_V.empty()) return false;

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> P_cols;
    std::vector<std::shared_ptr<SymbolicExpr>> diag;
    for (auto& pr : eigen_V) {
        for (auto& vec : pr.second) {
            auto vnode = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(vec));
            if (!vnode) return false;
            std::vector<std::shared_ptr<SymbolicExpr>> col;
            for (size_t i = 0; i < n; ++i) col.push_back(lamina::detail::make_expression_ptr(vnode->get(i, 0)));
            P_cols.push_back(col);
            diag.push_back(pr.first);
        }
    }
    /// 当特征向量数等于 n 时，矩阵可对角化，Jordan 型即对角阵。
    /// 缺陷情形（重根但特征向量不足）当前不支持广义特征向量链，返回 false。
    if (P_cols.size() != n) return false;

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> P_grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) P_grid[i][j] = P_cols[j][i];
    P = SymbolicExpr::matrix(P_grid);

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> J_grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < n; ++i) J_grid[i][i] = diag[i];
    J = SymbolicExpr::matrix(J_grid);
    return true;
}

} // namespace lamina
