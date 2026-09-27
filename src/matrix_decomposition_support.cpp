#include "internal/matrix_decomposition_support.hpp"
#include "internal/expression_analysis.hpp"

namespace LMCAS::detail::matrix_decomposition {
Result<std::shared_ptr<const MatrixNode>> validate_decomposition_matrix(
    const std::shared_ptr<SymbolicExpr>& A,
    bool require_square,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return Result<std::shared_ptr<const MatrixNode>>::failure(step.error());
    if (!A || !LMCAS::detail::node(A)) {
        return Result<std::shared_ptr<const MatrixNode>>::failure(
            CasErrc::InvalidArgument,
            "matrix input cannot be null",
            operation);
    }
    auto matrix = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(A));
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
    if (!expr || !LMCAS::detail::node(expr) ||
        !std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(expr))) {
        return Result<void>::failure(
            CasErrc::InternalInvariant,
            name + " output is not a matrix",
            operation);
    }
    return Result<void>::success();
}
namespace {
std::optional<Rational> exact_rational_sum(const AddNode& add) {
    Rational sum(0);
    for (const auto& operand : add.operands()) {
        auto value = exact_rational_expr(operand);
        if (!value) return std::nullopt;
        sum = sum + *value;
    }
    return sum;
}

std::optional<Rational> exact_rational_product(const MultiplyNode& mul) {
    Rational product(1);
    for (const auto& operand : mul.operands()) {
        auto value = exact_rational_expr(operand);
        if (!value) return std::nullopt;
        product = product * *value;
    }
    return product;
}
}

std::optional<Rational> exact_rational_expr(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return std::nullopt;
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return LMCAS::detail::exact_rational_value(*number);
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return exact_rational_sum(*add);
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return exact_rational_product(*mul);
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

namespace {
std::optional<Rational> simplify_exact_rational_expr(
    const std::shared_ptr<const SymbolicNode>& node) {
    auto expression = LMCAS::detail::make_expression_ptr(node);
    auto simplified = expression->simplify();
    return simplified
        ? exact_rational_expr(LMCAS::detail::node(simplified))
        : std::nullopt;
}
}

std::optional<Rational> exact_rational_expr(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) return std::nullopt;
    const auto& node = LMCAS::detail::node(expr);
    if (auto direct = exact_rational_expr(node)) return direct;
    return simplify_exact_rational_expr(node);
}

std::shared_ptr<SymbolicExpr> exact_number_expr(const Rational& value) {
    if (value.get_denominator() == BigInt(1)) {
        return LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<NumberNode>(
                std::variant<BigInt, Rational, lmmc_real_t>{
                    std::in_place_type<BigInt>, value.get_numerator()}));
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(
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
            const auto& node = matrix->get(row, col);
            auto value = exact_rational_expr(node);
            if (!value) {
                value = simplify_exact_rational_expr(node);
            }
            if (!value) return std::nullopt;
            values[row][col] = *value;
        }
    }
    return values;
}
std::shared_ptr<SymbolicExpr> identity_matrix_expr(size_t n) {
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(
        n, std::vector<std::shared_ptr<SymbolicExpr>>(n, SymbolicExpr::number(0)));
    for (size_t i = 0; i < n; ++i) {
        grid[i][i] = SymbolicExpr::number(1);
    }
    return SymbolicExpr::matrix(grid);
}

std::optional<std::vector<Rational>> exact_rectangular_diagonal_entries(
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

}
