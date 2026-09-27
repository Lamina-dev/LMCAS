#include "internal/symbolic_matrix_support.hpp"
#include "internal/normalization_utils.hpp"
#include "internal/visitors/normalization_visitor.hpp"

#include <new>

namespace LMCAS {
using namespace detail::matrix_api;
ExpressionResult matrix_multiply_checked(const std::shared_ptr<SymbolicExpr>& A,
                                         const std::shared_ptr<SymbolicExpr>& B,
                                         ComputationContext& context) {
    const std::string operation = "matrix_multiply";
    auto left = require_matrix(A, context, operation);
    if (!left) return ExpressionResult::failure(left.error());
    auto right = require_matrix(B, context, operation);
    if (!right) return ExpressionResult::failure(right.error());
    if (left.value()->cols() != right.value()->rows()) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument, "matrix dimensions are not compatible", operation);
    }
    try {
        detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
                                     context.limits().max_ast_nodes, "matrix_multiply");
        const auto count = normalization_check_count(
            &budget, left.value()->rows(), right.value()->cols());
        auto output_terms = context.require_expansion_terms(count, operation);
        if (!output_terms) { return ExpressionResult::failure(output_terms.error()); }
        auto inner_terms = context.require_expansion_terms(left.value()->cols(), operation);
        if (!inner_terms) { return ExpressionResult::failure(inner_terms.error()); }
        normalization_check_arithmetic<MultiplyNode>(&budget, 0, left.value(), right.value());
        auto product = detail::make_node<MultiplyNode>(
            std::vector<std::shared_ptr<const SymbolicNode>>{left.value(), right.value()});
        NormalizationVisitor visitor(context, detail::no_facts(), Domain::Real, &budget);
        product->accept(visitor);
        auto result = visitor.get_result();
        budget.measure(result);
        if (!std::dynamic_pointer_cast<const MatrixNode>(result)) {
            return ExpressionResult::failure(CasErrc::InternalInvariant,
                "matrix multiplication did not produce an explicit matrix", operation);
        }
        return detail::make_expression_ptr(std::move(result));
    } catch (const CasError& error) {
        return ExpressionResult::failure(error);
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
            "allocation failed while multiplying matrices", operation);
    } catch (const std::exception& ex) {
        return ExpressionResult::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
}

ExpressionResult matrix_multiply_checked(const std::shared_ptr<SymbolicExpr>& A,
                                         const std::shared_ptr<SymbolicExpr>& B) {
    ComputationContext context;
    return matrix_multiply_checked(A, B, context);
}


ExpressionResult matrix_determinant_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context) {
    const std::string operation = "matrix_determinant";
    auto matrix = require_square_matrix(A, context, operation);
    if (!matrix) return ExpressionResult::failure(matrix.error());
    auto determinant = detail::determinant_exact(
        exact_matrix_data(*matrix.value()), context, operation);
    if (!determinant) return ExpressionResult::failure(determinant.error());
    return ExpressionResult::success(std::move(determinant.value()));
}

ExpressionResult matrix_determinant_checked(const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_determinant_checked(A, context);
}


ExpressionResult matrix_inverse_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context) {
    const std::string operation = "matrix_inverse";
    auto matrix = require_square_matrix(A, context, operation);
    if (!matrix) return ExpressionResult::failure(matrix.error());
    auto inverse = detail::inverse_exact(
        exact_matrix_data(*matrix.value()), context, operation);
    if (!inverse) return ExpressionResult::failure(inverse.error());
    return ExpressionResult::success(
        exact_matrix_expression(inverse.value()));
}

ExpressionResult matrix_inverse_checked(const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_inverse_checked(A, context);
}


MatrixRankResult matrix_rank_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context) {
    const std::string operation = "matrix_rank";
    auto matrix = require_matrix(A, context, operation);
    if (!matrix) return MatrixRankResult::failure(matrix.error());
    auto rank = detail::rank_exact(
        exact_matrix_data(*matrix.value()), matrix.value()->cols(),
        context, operation);
    if (!rank) return MatrixRankResult::failure(rank.error());
    return MatrixRankResult::success(rank.value());
}

MatrixRankResult matrix_rank_checked(
    const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_rank_checked(A, context);
}

ExpressionResult matrix_rref_checked(
    const std::shared_ptr<SymbolicExpr>& A,
    ComputationContext& context) {
    const std::string operation = "matrix_rref";
    auto matrix = require_matrix(A, context, operation);
    if (!matrix) return ExpressionResult::failure(matrix.error());
    auto rref = detail::rref_exact(
        exact_matrix_data(*matrix.value()), matrix.value()->cols(),
        context, operation);
    if (!rref) return ExpressionResult::failure(rref.error());
    return ExpressionResult::success(exact_matrix_expression(rref.value()));
}

ExpressionResult matrix_rref_checked(
    const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_rref_checked(A, context);
}

}
