#include "internal/matrix_decomposition_support.hpp"

#include <algorithm>
#include <new>

namespace LMCAS {
using namespace detail::matrix_decomposition;
namespace {
std::shared_ptr<SymbolicExpr> abs_expression(const std::shared_ptr<SymbolicExpr>& e) {
    auto node = detail::make_node<FunctionNode>(FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<const SymbolicNode>>{detail::node(e)});
    return detail::make_expression_ptr(node);
}

struct NormSum {
    std::shared_ptr<SymbolicExpr> expression;
    std::optional<Rational> exact;
};

NormSum matrix_axis_sum(const MatrixNode& matrix, std::size_t axis, bool columns) {
    auto sum = SymbolicExpr::number(0);
    std::optional<Rational> exact_sum = Rational(0);
    const auto count = columns ? matrix.rows() : matrix.cols();
    for (std::size_t index = 0; index < count; ++index) {
        auto entry = detail::make_expression_ptr(columns
            ? matrix.get(index, axis) : matrix.get(axis, index));
        sum = SymbolicExpr::add(sum, abs_expression(entry));
        if (exact_sum) {
            auto exact_entry = exact_rational_expr(entry);
            if (exact_entry) {
                *exact_sum = *exact_sum + ((*exact_entry < Rational(0))
                    ? (Rational(0) - *exact_entry) : *exact_entry);
            } else {
                exact_sum.reset();
            }
        }
    }
    sum = sum->simplify();
    return {sum, exact_sum ? exact_sum : exact_rational_expr(sum)};
}

std::shared_ptr<SymbolicExpr> maximum_axis_sum(const MatrixNode& matrix, bool columns) {
    std::optional<Rational> maximum_exact;
    std::vector<std::shared_ptr<const SymbolicNode>> candidates;
    const auto count = columns ? matrix.cols() : matrix.rows();
    for (std::size_t axis = 0; axis < count; ++axis) {
        auto candidate = matrix_axis_sum(matrix, axis, columns);
        if (candidate.exact) {
            if (!maximum_exact || *candidate.exact > *maximum_exact)
                maximum_exact = *candidate.exact;
        } else {
            candidates.push_back(detail::node(candidate.expression));
        }
    }
    if (maximum_exact)
        candidates.push_back(detail::node(exact_number_expr(*maximum_exact)));
    if (candidates.empty()) return nullptr;
    if (candidates.size() == 1) return detail::make_expression_ptr(candidates.front());
    return detail::make_expression_ptr(detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Max, std::move(candidates)));
}
}

std::shared_ptr<SymbolicExpr> matrix_trace(const std::shared_ptr<SymbolicExpr>& A) {
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(A));
    if (!mat || mat->rows() != mat->cols()) {
        return nullptr;
    }
    auto sum = SymbolicExpr::number(0);
    for (size_t i = 0; i < mat->rows(); ++i) {
        sum = SymbolicExpr::add(sum, LMCAS::detail::make_expression_ptr(mat->get(i, i)));
    }
    return sum->simplify();
}
std::shared_ptr<SymbolicExpr> kronecker(const std::shared_ptr<SymbolicExpr>& A,
    const std::shared_ptr<SymbolicExpr>& B) {
    auto a = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(A));
    auto b = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(B));
    if (!a || !b) {
        return nullptr;
    }
    size_t m = a->rows(), n = a->cols(), p = b->rows(), q = b->cols();
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(m * p,
        std::vector<std::shared_ptr<SymbolicExpr>>(n * q));
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j) {
            auto aij = LMCAS::detail::make_expression_ptr(a->get(i, j));
            for (size_t k = 0; k < p; ++k)
                for (size_t l = 0; l < q; ++l) {
                    auto bkl = LMCAS::detail::make_expression_ptr(b->get(k, l));
                    grid[i * p + k][j * q + l] = SymbolicExpr::multiply(aij, bkl)->simplify();
                }
        }
    return SymbolicExpr::matrix(grid);
}
std::shared_ptr<SymbolicExpr> matrix_norm(const std::shared_ptr<SymbolicExpr>& A,
    const std::string& type) {
    auto mat = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(A));
    if (!mat) {
        return nullptr;
    }
    size_t m = mat->rows(), n = mat->cols();

    if (type == "frobenius") {
        auto sum = SymbolicExpr::number(0);
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j) {
                auto e = LMCAS::detail::make_expression_ptr(mat->get(i, j));
                sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(e, e));
            }
        return SymbolicExpr::sqrt(sum)->simplify();
    }
    if (type == "1") {
        return maximum_axis_sum(*mat, true);
    }
    if (type == "inf") {
        return maximum_axis_sum(*mat, false);
    }
    return nullptr;
}

}
