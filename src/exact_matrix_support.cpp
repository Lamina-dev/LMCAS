#include "internal/expression_transform.hpp"
#include "internal/exact_matrix_support.hpp"
#include "internal/normalization_utils.hpp"
#include "assumption_context.hpp"
#include "internal/expression_analysis.hpp"

#include <limits>
#include <new>

namespace LMCAS::detail {
namespace matrix_kernel {
Result<void> validate_matrix(
    const ExactMatrixData& matrix,
    const std::string& operation) {
    if (matrix.rows == 0 || matrix.cols == 0) {
        return Result<void>::failure(
            CasErrc::InvalidArgument,
            "exact matrix dimensions must be non-zero", operation);
    }
    if (matrix.rows > std::numeric_limits<std::size_t>::max() / matrix.cols ||
        matrix.entries.size() != matrix.rows * matrix.cols) {
        return Result<void>::failure(
            CasErrc::InvalidArgument,
            "exact matrix storage does not match its dimensions", operation);
    }
    for (const auto& entry : matrix.entries) {
        if (!entry || !node(entry)) {
            return Result<void>::failure(
                CasErrc::InvalidArgument,
                "exact matrix entries cannot be null", operation);
        }
    }
    return Result<void>::success();
}


bool contains_approximate_number(const ExprPtr& expression) {
    return expression &&
        contains_inexact_number(node(expression));
}

std::optional<Rational> exact_rational(const ExprPtr& expression) {
    if (!expression || !node(expression)) {
        return std::nullopt;
    }
    return exact_rational_value(node(expression));
}

ExprPtr simplify_expr(ExprPtr expression, ComputationContext& context) {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    if (!expression) {
        return expression;
    }
    auto simplified = detail::make_expression_ptr(detail::normalize_node(detail::node(expression), transform_budget));
    return simplified ? simplified : expression;
}

ExprPtr exact_add(const ExprPtr& left, const ExprPtr& right, ComputationContext& context) {
    return simplify_expr(SymbolicExpr::add(left, right), context);
}

ExprPtr exact_negate(const ExprPtr& expression, ComputationContext& context) {
    return simplify_expr(SymbolicExpr::multiply(
        SymbolicExpr::number(-1), expression), context);
}

ExprPtr exact_subtract(const ExprPtr& left, const ExprPtr& right, ComputationContext& context) {
    return exact_add(left, exact_negate(right, context), context);
}

ExprPtr exact_multiply(const ExprPtr& left, const ExprPtr& right, ComputationContext& context) {
    return simplify_expr(SymbolicExpr::multiply(left, right), context);
}

ExprPtr exact_divide(const ExprPtr& numerator, const ExprPtr& denominator, ComputationContext& context) {
    return simplify_expr(SymbolicExpr::divide(numerator, denominator), context);
}

void swap_rows(ExactMatrixData& matrix, std::size_t first, std::size_t second) {
    if (first == second) {
        return;
    }
    for (std::size_t column = 0; column < matrix.cols; ++column) {
        std::swap(matrix.at(first, column), matrix.at(second, column));
    }
}

Result<std::optional<std::size_t>> choose_pivot(
    const ExactMatrixData& matrix,
    std::size_t first_row,
    std::size_t column,
    ComputationContext& context,
    const std::string& operation) {
    bool saw_unknown = false;
    for (std::size_t row = first_row; row < matrix.rows; ++row) {
        auto proof = classify_exact_zero(matrix.at(row, column), context, operation);
        if (!proof) {
            return Result<std::optional<std::size_t>>::failure(proof.error());
        }
        if (proof.value() == ZeroProof::NonZero) {
            return Result<std::optional<std::size_t>>::success(row);
        }
        saw_unknown = saw_unknown || proof.value() == ZeroProof::Unknown;
    }
    if (saw_unknown) {
        return Result<std::optional<std::size_t>>::failure(
            CasErrc::Inconclusive,
            "exact pivot is not provably zero or non-zero", operation);
    }
    return Result<std::optional<std::size_t>>::success(std::nullopt);
}

}
using namespace matrix_kernel;
namespace {
std::optional<ZeroProof> classify_numeric_power(const ExprPtr& simplified) {
        if (auto power =
                std::dynamic_pointer_cast<const PowerNode>(node(simplified))) {
            auto base = exact_rational(
                make_expression_ptr(power->base()));
            auto exponent = exact_rational(
                make_expression_ptr(power->exponent()));
            if (base && exponent) {
                if (*base > Rational(0)) {
                    return ZeroProof::NonZero;
                }
                if (*base == Rational(0) && *exponent > Rational(0)) {
                    return ZeroProof::Zero;
                }
            }
        }
    return std::nullopt;
}

std::optional<ZeroProof> classify_numeric_sqrt(const ExprPtr& simplified) {
    auto function =
        std::dynamic_pointer_cast<const FunctionNode>(node(simplified));
    if (!function) {
        return std::nullopt;
    }
    if (function->type() != FunctionNode::FuncType::Sqrt ||
        function->arguments().size() != 1) {
        return std::nullopt;
    }
    auto argument = exact_rational(
        make_expression_ptr(function->arguments()[0]));
    if (!argument || *argument < Rational(0)) {
        return std::nullopt;
    }
    return *argument == Rational(0) ? ZeroProof::Zero : ZeroProof::NonZero;
}

std::optional<ZeroProof> classify_numeric_value(const ExprPtr& expression) {
    if (auto value = exact_rational(expression)) {
        return *value == Rational(0) ? ZeroProof::Zero : ZeroProof::NonZero;
    }
    if (auto proof = classify_numeric_power(expression)) {
        return proof;
    }
    if (auto proof = classify_numeric_sqrt(expression)) {
        return proof;
    }
    if (auto root = std::dynamic_pointer_cast<const RootOfNode>(node(expression));
        root && !root->exact_id().polynomial.coeffs.empty() &&
        root->exact_id().polynomial.coeffs[0] != Rational(0)) {
        return ZeroProof::NonZero;
    }
    return std::nullopt;
}
}

Result<ZeroProof> classify_exact_zero(
    const ExprPtr& expression,
    ComputationContext& context,
    const std::string& operation) {
    auto access = context.consume_steps(1, operation);
    if (!access) {
        return Result<ZeroProof>::failure(access.error());
    }
    if (!expression || !node(expression)) {
        return Result<ZeroProof>::failure(
            CasErrc::InvalidArgument,
            "zero classification requires an expression", operation);
    }
    try {
        if (contains_approximate_number(expression)) {
            return Result<ZeroProof>::success(ZeroProof::Unknown);
        }
        auto simplified = simplify_expr(expression, context);
        if (contains_approximate_number(simplified)) {
            return Result<ZeroProof>::success(ZeroProof::Unknown);
        }
        if (auto proof = classify_numeric_value(simplified)) {
            return Result<ZeroProof>::success(*proof);
        }
        if (!context.assumptions()) {
            return Result<ZeroProof>::success(ZeroProof::Unknown);
        }
        auto assumption = context.assumptions()->is_nonzero_checked(*simplified, context);
        if (!assumption) {
            return Result<ZeroProof>::failure(assumption.error());
        }
        if (assumption.value() == Tribool::True) {
            return Result<ZeroProof>::success(ZeroProof::NonZero);
        }
        if (assumption.value() == Tribool::False) {
            return Result<ZeroProof>::success(ZeroProof::Zero);
        }
        return Result<ZeroProof>::success(ZeroProof::Unknown);
    } catch (const CasError& error) {
        return Result<ZeroProof>::failure(error);
    } catch (const std::length_error& error) {
        return Result<ZeroProof>::failure(CasErrc::ResourceLimit, error.what(), operation);
    } catch (const std::bad_alloc&) {
        return Result<ZeroProof>::failure(
            CasErrc::ResourceLimit,
            "zero classification allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<ZeroProof>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

}
