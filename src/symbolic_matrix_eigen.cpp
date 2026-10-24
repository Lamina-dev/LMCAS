#include <stdexcept>
#include "internal/expression_transform.hpp"
#include "internal/symbolic_matrix_support.hpp"

#include "solve_strategies.hpp"
#include "internal/assumption_facts.hpp"
#include "assumption_context.hpp"
#include <new>
#include <set>

namespace LMCAS {
using namespace detail::matrix_api;
ExpressionResult matrix_characteristic_polynomial_checked(
    const std::shared_ptr<SymbolicExpr>& A, const std::string& variable,
    ComputationContext& context) {
    constexpr const char* operation = "matrix_characteristic_polynomial";
    auto matrix_result = require_square_matrix(A, context, operation);
    if (!matrix_result) return ExpressionResult::failure(matrix_result.error());
    if (variable.empty()) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument,
            "characteristic polynomial variable cannot be empty", operation);
    }
    try {
        const auto& matrix = *matrix_result.value();
        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> entries(
            matrix.rows(),
            std::vector<std::shared_ptr<SymbolicExpr>>(matrix.cols()));
        auto lambda = SymbolicExpr::variable(variable);
        for (std::size_t row = 0; row < matrix.rows(); ++row) {
            for (std::size_t column = 0; column < matrix.cols(); ++column) {
                auto value = detail::make_expression_ptr(
                    matrix.get(row, column));
                entries[row][column] = row == column
                    ? SymbolicExpr::add(
                          value,
                          SymbolicExpr::multiply(
                              lambda, SymbolicExpr::number(-1)))
                    : std::move(value);
            }
        }
        auto determinant = matrix_determinant_checked(
            SymbolicExpr::matrix(entries), context);
        if (!determinant) return ExpressionResult::failure(determinant.error());
        return ExpressionResult::success(
            std::move(determinant.value()));
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while constructing characteristic polynomial",
            operation);
    } catch (const std::exception& ex) {
        return ExpressionResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

ExpressionResult matrix_characteristic_polynomial_checked(
    const std::shared_ptr<SymbolicExpr>& A, const std::string& variable) {
    ComputationContext context;
    return matrix_characteristic_polynomial_checked(A, variable, context);
}

MatrixEigenvalueResult matrix_eigenvalues_checked(
    const std::shared_ptr<SymbolicExpr>& A, ComputationContext& context) try {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    constexpr const char* operation = "matrix_eigenvalues";
    auto matrix_result = require_square_matrix(A, context, operation);
    if (!matrix_result) {
        return MatrixEigenvalueResult::failure(matrix_result.error());
    }
    const auto& matrix = *matrix_result.value();
    bool triangular = true;
    for (std::size_t row = 0; row < matrix.rows() && triangular; ++row) {
        for (std::size_t column = 0; column < row; ++column) {
            auto value =
                detail::make_expression_ptr(detail::normalize_node(matrix.get(row, column), transform_budget));
            if (!detail::node(value)->is_zero()) {
                triangular = false;
                break;
            }
        }
    }
    if (triangular) {
        std::vector<std::shared_ptr<SymbolicExpr>> values;
        values.reserve(matrix.rows());
        for (std::size_t index = 0; index < matrix.rows(); ++index) {
            values.push_back(
                detail::make_expression_ptr(detail::normalize_node(matrix.get(index, index), transform_budget)));
        }
        return MatrixEigenvalueResult::success(std::move(values));
    }

    auto polynomial = matrix_characteristic_polynomial_checked(
        A, "lambda", context);
    if (!polynomial) {
        return MatrixEigenvalueResult::failure(polynomial.error());
    }
    auto roots = solve_finite_checked(
        polynomial.value(), "lambda", context, SolveOptions{});
    if (!roots) return MatrixEigenvalueResult::failure(roots.error());
    std::vector<std::shared_ptr<SymbolicExpr>> values;
    std::set<std::string> seen;
    for (auto& root : roots.value()) {
        const auto key = root->to_string();
        if (seen.insert(key).second) values.push_back(std::move(root));
    }
    if (values.empty()) {
        return MatrixEigenvalueResult::failure(
            CasErrc::Inconclusive,
            "eigenvalue solver produced no finite roots", operation);
    }
    return MatrixEigenvalueResult::success(std::move(values));
} catch (const CasError& error) {
    return MatrixEigenvalueResult::failure(error);
} catch (const std::bad_alloc&) {
    return MatrixEigenvalueResult::failure(CasErrc::ResourceLimit, "allocation failed", "matrix_eigenvalues");
} catch (const std::length_error& error) {
    return MatrixEigenvalueResult::failure(CasErrc::ResourceLimit, error.what(), "matrix_eigenvalues");
} catch (const std::exception& error) {
    return MatrixEigenvalueResult::failure(CasErrc::InternalInvariant, error.what(), "matrix_eigenvalues");
}

MatrixEigenvalueResult matrix_eigenvalues_checked(
    const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_eigenvalues_checked(A, context);
}

namespace {
using Eigenspaces = std::vector<std::pair<ExprPtr, std::vector<std::vector<ExprPtr>>>>;

Result<bool> eigenvalue_already_collected(const ExprPtr& lambda,
    const Eigenspaces& spaces, ComputationContext& context,
    const detail::AssumptionFacts& facts) try {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    constexpr const char* operation = "matrix_eigenspaces";
    for (const auto& space : spaces) {
        if (detail::node(lambda)->equals(*detail::node(space.first))) { return true; }
        auto difference = SymbolicExpr::add(lambda,
            SymbolicExpr::multiply(SymbolicExpr::number(-1), space.first));
        auto proof = detail::classify_exact_zero(difference, context, operation);
        if (!proof) { return Result<bool>::failure(proof.error()); }
        if (proof.value() == detail::ZeroProof::Zero) { return true; }
        auto nonzero = detail::query_nonzero_value(detail::normalize_node(detail::node(difference), transform_budget), facts, Domain::Complex, context);
        if (!nonzero) { return Result<bool>::failure(nonzero.error()); }
        if (nonzero.value() != Tribool::True) {
            return Result<bool>::failure(CasErrc::Inconclusive,
                "eigenvalues cannot be proved pointwise distinct", operation);
        }
    }
    return false;
} catch (const CasError& error) {
    return Result<bool>::failure(error);
} catch (const std::bad_alloc&) {
    return Result<bool>::failure(CasErrc::ResourceLimit, "allocation failed", "matrix_eigenspaces");
} catch (const std::length_error& error) {
    return Result<bool>::failure(CasErrc::ResourceLimit, error.what(), "matrix_eigenspaces");
} catch (const std::exception& error) {
    return Result<bool>::failure(CasErrc::InternalInvariant, error.what(), "matrix_eigenspaces");
}

MatrixEigenvectorResult eigenvalue_nullspace(const MatrixNode& matrix,
    const ExprPtr& lambda, ComputationContext& context) try {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> shifted(
        matrix.rows(),
        std::vector<std::shared_ptr<SymbolicExpr>>(matrix.cols()));
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t column = 0; column < matrix.cols(); ++column) {
            auto value = detail::make_expression_ptr(matrix.get(row, column));
            shifted[row][column] = row == column
                ? detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::add(value, SymbolicExpr::multiply(
                      SymbolicExpr::number(-1), lambda))), transform_budget))
                : std::move(value);
        }
    }
    auto nullspace = matrix_nullspace_checked(SymbolicExpr::matrix(shifted), context);
    if (!nullspace) { return MatrixEigenvectorResult::failure(nullspace.error()); }
    if (nullspace.value().empty()) {
        return MatrixEigenvectorResult::failure(CasErrc::Inconclusive,
            "eigenvalue has no proved nonzero nullspace basis", "matrix_eigenspaces");
    }
    return nullspace;
} catch (const CasError& error) {
    return MatrixEigenvectorResult::failure(error);
} catch (const std::bad_alloc&) {
    return MatrixEigenvectorResult::failure(CasErrc::ResourceLimit, "allocation failed", "eigenvalue_nullspace");
} catch (const std::length_error& error) {
    return MatrixEigenvectorResult::failure(CasErrc::ResourceLimit, error.what(), "eigenvalue_nullspace");
} catch (const std::exception& error) {
    return MatrixEigenvectorResult::failure(CasErrc::InternalInvariant, error.what(), "eigenvalue_nullspace");
}
}

Result<std::vector<std::pair<ExprPtr, std::vector<std::vector<ExprPtr>>>>>
detail::matrix_eigenspaces_checked(const ExprPtr& A, ComputationContext& context) try {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    using SpacesResult = Result<Eigenspaces>;
    constexpr const char* operation = "matrix_eigenspaces";
    auto matrix_result = require_square_matrix(A, context, operation);
    if (!matrix_result) { return SpacesResult::failure(matrix_result.error()); }
    auto eigenvalues = matrix_eigenvalues_checked(A, context);
    if (!eigenvalues) { return SpacesResult::failure(eigenvalues.error()); }
    const auto& matrix = *matrix_result.value();
    Eigenspaces spaces;
    AssumptionContext assumptions;
    detail::AssumptionFacts facts(context.assumptions() ? *context.assumptions() : assumptions);
    try {
        for (const auto& eigenvalue : eigenvalues.value()) {
            auto lambda = detail::make_expression_ptr(detail::normalize_node(detail::node(eigenvalue), transform_budget));
            auto duplicate = eigenvalue_already_collected(lambda, spaces, context, facts);
            if (!duplicate) { return SpacesResult::failure(duplicate.error()); }
            if (duplicate.value()) { continue; }
            auto nullspace = eigenvalue_nullspace(matrix, lambda, context);
            if (!nullspace) { return SpacesResult::failure(nullspace.error()); }
            spaces.emplace_back(std::move(lambda), std::move(nullspace.value()));
        }
        return SpacesResult::success(std::move(spaces));
    } catch (const std::bad_alloc&) {
        return SpacesResult::failure(CasErrc::ResourceLimit,
            "allocation failed while computing eigenspaces", operation);
    } catch (const std::length_error& error) { return Result<std::vector<std::pair<ExprPtr, std::vector<std::vector<ExprPtr>>>>>::failure(CasErrc::ResourceLimit, error.what(), "matrix_eigenspaces"); } catch (const std::exception& ex) {
        return SpacesResult::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
} catch (const CasError& error) {
    return Result<std::vector<std::pair<ExprPtr, std::vector<std::vector<ExprPtr>>>>>::failure(error);
} catch (const std::bad_alloc&) {
    return Result<std::vector<std::pair<ExprPtr, std::vector<std::vector<ExprPtr>>>>>::failure(CasErrc::ResourceLimit, "allocation failed", "matrix_eigenspaces");
} catch (const std::length_error& error) {
    return Result<std::vector<std::pair<ExprPtr, std::vector<std::vector<ExprPtr>>>>>::failure(CasErrc::ResourceLimit, error.what(), "matrix_eigenspaces");
} catch (const std::exception& error) {
    return Result<std::vector<std::pair<ExprPtr, std::vector<std::vector<ExprPtr>>>>>::failure(CasErrc::InternalInvariant, error.what(), "matrix_eigenspaces");
}

MatrixEigenvectorResult matrix_eigenvectors_checked(
    const std::shared_ptr<SymbolicExpr>& A, ComputationContext& context) {
    auto spaces = detail::matrix_eigenspaces_checked(A, context);
    if (!spaces) return MatrixEigenvectorResult::failure(spaces.error());
    try {
        std::vector<std::vector<ExprPtr>> vectors;
        for (auto& space : spaces.value())
            for (auto& vector : space.second) vectors.push_back(std::move(vector));
        return MatrixEigenvectorResult::success(std::move(vectors));
    } catch (const std::bad_alloc&) {
        return MatrixEigenvectorResult::failure(CasErrc::ResourceLimit,
            "allocation failed while collecting eigenvectors", "matrix_eigenvectors");
    }
}

MatrixEigenvectorResult matrix_eigenvectors_checked(
    const std::shared_ptr<SymbolicExpr>& A) {
    ComputationContext context;
    return matrix_eigenvectors_checked(A, context);
}

}
