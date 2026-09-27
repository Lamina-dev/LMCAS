#include "internal/matrix_decomposition_support.hpp"

#include "internal/symbolic_matrix_support.hpp"
#include "internal/assumption_facts.hpp"
#include "assumption_context.hpp"
#include <new>
#include <optional>

namespace LMCAS {
using namespace detail::matrix_decomposition;
namespace {
Result<void> verify_jordan_identity(const ExprPtr& A,
    const JordanDecomposition& value, const MatrixNode& matrix,
    ComputationContext& context, const std::string& operation) {
    auto ap = matrix_multiply_checked(A, value.P, context);
    if (!ap) { return Result<void>::failure(ap.error()); }
    auto pj = matrix_multiply_checked(value.P, value.J, context);
    if (!pj) { return Result<void>::failure(pj.error()); }
    auto left = std::dynamic_pointer_cast<const MatrixNode>(detail::node(ap.value()));
    auto right = std::dynamic_pointer_cast<const MatrixNode>(detail::node(pj.value()));
    if (!left || !right) {
        return Result<void>::failure(CasErrc::InternalInvariant,
            "Jordan matrix product did not produce a matrix", operation);
    }
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t col = 0; col < matrix.cols(); ++col) {
            auto difference = SymbolicExpr::add(detail::make_expression_ptr(left->get(row, col)),
                SymbolicExpr::multiply(SymbolicExpr::number(-1),
                    detail::make_expression_ptr(right->get(row, col))));
            auto proof = detail::classify_exact_zero(difference, context, operation);
            if (!proof) { return Result<void>::failure(proof.error()); }
            if (proof.value() != detail::ZeroProof::Zero) {
                return Result<void>::failure(proof.value() == detail::ZeroProof::Unknown
                    ? CasErrc::Inconclusive : CasErrc::InternalInvariant,
                    "Jordan intertwining identity is not proved at (" +
                        std::to_string(row) + "," + std::to_string(col) +
                        "): " + difference->to_string(), operation);
            }
        }
    }
    return Result<void>::success();
}

bool jordan_dimensions_match(const MatrixNode& a, const MatrixNode& p,
    const MatrixNode& j) {
    return a.rows() == a.cols() && p.rows() == a.rows() && p.cols() == a.cols() &&
        j.rows() == a.rows() && j.cols() == a.cols();
}

Result<void> verify_jordan_decomposition(const ExprPtr& A,
    const JordanDecomposition& value, ComputationContext& context) {
    const std::string operation = "jordan_form.verify";
    auto a = std::dynamic_pointer_cast<const MatrixNode>(detail::node(A));
    auto p = std::dynamic_pointer_cast<const MatrixNode>(detail::node(value.P));
    auto j = std::dynamic_pointer_cast<const MatrixNode>(detail::node(value.J));
    if (!a || !p || !j || !jordan_dimensions_match(*a, *p, *j)) {
        return Result<void>::failure(CasErrc::InternalInvariant,
            "Jordan certificate has invalid dimensions", operation);
    }
    auto rank = matrix_rank_checked(value.P, context);
    if (!rank) { return Result<void>::failure(rank.error()); }
    if (rank.value() != a->rows()) {
        return Result<void>::failure(CasErrc::InternalInvariant,
            "Jordan change of basis is singular", operation);
    }
    return verify_jordan_identity(A, value, *a, context, operation);
}

JordanDecompositionResult certified_jordan(const ExprPtr& A,
    JordanDecomposition value, ComputationContext& context) {
    auto proof = verify_jordan_decomposition(A, value, context);
    if (!proof) return JordanDecompositionResult::failure(proof.error());
    return JordanDecompositionResult::success(std::move(value));
}

JordanDecompositionResult diagonalization_checked(const ExprPtr& A,
    ComputationContext& context) {
    auto matrix = validate_decomposition_matrix(A, true, context, "matrix.diagonalize");
    if (!matrix) { return JordanDecompositionResult::failure(matrix.error()); }
    if (auto diagonal = exact_diagonal_entries(matrix.value())) {
        return certified_jordan(A, {rectangular_diagonal_expr(matrix.value()->rows(),
            matrix.value()->cols(), *diagonal), identity_matrix_expr(matrix.value()->rows())}, context);
    }
    auto spaces = detail::matrix_eigenspaces_checked(A, context);
    if (!spaces) { return JordanDecompositionResult::failure(spaces.error()); }
    const auto n = matrix.value()->rows();
    std::size_t count = 0;
    for (const auto& space : spaces.value()) { count += space.second.size(); }
    if (count != n) {
        return JordanDecompositionResult::failure(CasErrc::Inconclusive,
            "ordinary eigenspaces do not span the matrix", "matrix.diagonalize");
    }
    std::vector<std::vector<ExprPtr>> p(n, std::vector<ExprPtr>(n));
    std::vector<std::vector<ExprPtr>> j(n,
        std::vector<ExprPtr>(n, SymbolicExpr::number(0)));
    std::size_t column = 0;
    for (const auto& space : spaces.value()) {
        for (const auto& vector : space.second) {
            if (vector.size() != n) {
                return JordanDecompositionResult::failure(CasErrc::InternalInvariant,
                    "eigenspace vector has invalid dimensions", "matrix.diagonalize");
            }
            for (std::size_t row = 0; row < n; ++row) { p[row][column] = vector[row]; }
            j[column][column] = space.first;
            ++column;
        }
    }
    return certified_jordan(A,
        {SymbolicExpr::matrix(j), SymbolicExpr::matrix(p)}, context);
}

struct QuadraticSigns {
    bool has_pos = false;
    bool has_neg = false;
    bool any_zero = false;
    void include(int proved_sign) {
        if (proved_sign > 0) has_pos = true;
        else if (proved_sign < 0) has_neg = true;
        else any_zero = true;
    }
    std::string classification() const {
        if (has_pos && has_neg) { return "indefinite"; }
        if (has_pos) { return any_zero ? "positive_semidefinite" : "positive_definite"; }
        if (has_neg) { return any_zero ? "negative_semidefinite" : "negative_definite"; }
        return any_zero ? "positive_semidefinite" : "unknown";
    }
};

std::optional<int> symbolic_quadratic_sign(const ExprPtr& value,
    const detail::AssumptionFacts& facts, ComputationContext& context) {
    auto real = detail::query_real_value(detail::node(value), facts, context);
    if (!real || real.value() != Tribool::True) { return std::nullopt; }
    auto zero = detail::classify_exact_zero(value, context, "matrix.quadratic_sign");
    if (!zero) { return std::nullopt; }
    if (zero.value() == detail::ZeroProof::Zero) { return 0; }
    auto positive = detail::query_positive_value(detail::node(value), facts, context);
    if (!positive) { return std::nullopt; }
    if (positive.value() == Tribool::True) { return 1; }
    auto negative_value = SymbolicExpr::multiply(SymbolicExpr::number(-1), value);
    auto negative = detail::query_positive_value(detail::node(negative_value), facts, context);
    if (!negative || negative.value() != Tribool::True) { return std::nullopt; }
    return -1;
}

std::optional<int> quadratic_eigenvalue_sign(const ExprPtr& eigenvalue,
    const detail::AssumptionFacts& facts, ComputationContext& context) {
    auto value = eigenvalue->simplify();
    if (auto rational = exact_rational_expr(value)) {
        return *rational > Rational(0) ? 1 : (*rational < Rational(0) ? -1 : 0);
    }
    if (const auto* number = dynamic_cast<const NumberNode*>(detail::node(value).get())) {
        if (const auto* approximate = std::get_if<lmmc_real_t>(&number->value())) {
            if (!std::isfinite(*approximate)) { return std::nullopt; }
            return *approximate > 0 ? 1 : (*approximate < 0 ? -1 : 0);
        }
    }
    return symbolic_quadratic_sign(value, facts, context);
}

bool is_defective_upper_pair(const std::vector<std::vector<Rational>>& matrix) {
    return matrix.size() == 2 && matrix[0].size() == 2 &&
        matrix[0][0] == matrix[1][1] && matrix[1][0] == Rational(0) &&
        matrix[0][1] != Rational(0);
}

}

JordanDecompositionResult jordan_form_checked(const ExprPtr& A,
    ComputationContext& context) {
    const std::string operation = "jordan_form";
    auto matrix = validate_decomposition_matrix(A, true, context, operation);
    if (!matrix) { return JordanDecompositionResult::failure(matrix.error()); }
    auto budget = context.consume_steps(matrix.value()->rows() * matrix.value()->cols() * 48 + 32,
        operation);
    if (!budget) { return JordanDecompositionResult::failure(budget.error()); }
    try {
        if (auto diagonal = exact_diagonal_entries(matrix.value())) {
            return certified_jordan(A, {rectangular_diagonal_expr(
                matrix.value()->rows(), matrix.value()->cols(), *diagonal),
                identity_matrix_expr(matrix.value()->rows())}, context);
        }
        auto exact = exact_rational_matrix(matrix.value());
        if (exact && is_defective_upper_pair(*exact)) {
            auto zero = SymbolicExpr::number(0);
            auto one = SymbolicExpr::number(1);
            auto lambda = exact_number_expr((*exact)[0][0]);
            auto j = SymbolicExpr::matrix({{lambda, one}, {zero, lambda}});
            auto p = SymbolicExpr::matrix({{exact_number_expr((*exact)[0][1]), zero}, {zero, one}});
            return certified_jordan(A, {std::move(j), std::move(p)}, context);
        }
        return diagonalization_checked(A, context);
    } catch (const std::bad_alloc&) {
        return JordanDecompositionResult::failure(CasErrc::ResourceLimit,
            "allocation failed while calculating Jordan form", operation);
    } catch (const std::exception& ex) {
        return JordanDecompositionResult::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
}

JordanDecompositionResult jordan_form_checked(const ExprPtr& A) {
    ComputationContext context;
    return jordan_form_checked(A, context);
}

ExprPtr matrix_exp(const ExprPtr& A) {
    ComputationContext context;
    auto decomposition = diagonalization_checked(A, context);
    if (!decomposition) return SymbolicExpr::exp(A);
    auto inverse = matrix_inverse_checked(decomposition.value().P, context);
    if (!inverse) return SymbolicExpr::exp(A);
    auto diagonal = std::dynamic_pointer_cast<const MatrixNode>(detail::node(decomposition.value().J));
    const auto n = diagonal->rows();
    std::vector<std::vector<ExprPtr>> values(n, std::vector<ExprPtr>(n, SymbolicExpr::number(0)));
    for (std::size_t i = 0; i < n; ++i)
        values[i][i] = SymbolicExpr::exp(detail::make_expression_ptr(diagonal->get(i, i)));
    return SymbolicExpr::multiply(decomposition.value().P,
        SymbolicExpr::multiply(SymbolicExpr::matrix(values), inverse.value()))->simplify();
}

ExprPtr matrix_log(const ExprPtr& A) {
    ComputationContext context;
    auto decomposition = diagonalization_checked(A, context);
    if (!decomposition) return nullptr;
    auto inverse = matrix_inverse_checked(decomposition.value().P, context);
    if (!inverse) return nullptr;
    auto diagonal = std::dynamic_pointer_cast<const MatrixNode>(detail::node(decomposition.value().J));
    const auto n = diagonal->rows();
    AssumptionContext assumptions;
    detail::AssumptionFacts facts(assumptions);
    std::vector<std::vector<ExprPtr>> values(n, std::vector<ExprPtr>(n, SymbolicExpr::number(0)));
    for (std::size_t i = 0; i < n; ++i) {
        auto positive = detail::query_positive_value(diagonal->get(i, i), facts, context);
        if (!positive || positive.value() != Tribool::True) return nullptr;
        values[i][i] = SymbolicExpr::ln(detail::make_expression_ptr(diagonal->get(i, i)));
    }
    return SymbolicExpr::multiply(decomposition.value().P,
        SymbolicExpr::multiply(SymbolicExpr::matrix(values), inverse.value()))->simplify();
}

std::string classify_quadratic_form(const ExprPtr& A) {
    ComputationContext context;
    auto eigenvalues = matrix_eigenvalues_checked(A, context);
    if (!eigenvalues || eigenvalues.value().empty()) { return "unknown"; }
    AssumptionContext assumptions;
    detail::AssumptionFacts facts(assumptions);
    QuadraticSigns signs;
    for (const auto& eigenvalue : eigenvalues.value()) {
        auto sign = quadratic_eigenvalue_sign(eigenvalue, facts, context);
        if (!sign) { return "unknown"; }
        signs.include(*sign);
    }
    return signs.classification();
}
}
