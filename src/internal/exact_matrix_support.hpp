#include "computation_context.hpp"
#pragma once

#include "internal/exact_matrix.hpp"

#include <optional>

namespace LMCAS::detail::matrix_kernel {

Result<void> validate_matrix(const ExactMatrixData& matrix,
                             const std::string& operation);
std::optional<Rational> exact_rational(const ExprPtr& expression);
ExprPtr simplify_expr(ExprPtr expression, ComputationContext& context);
ExprPtr exact_add(const ExprPtr& left, const ExprPtr& right, ComputationContext& context);
ExprPtr exact_negate(const ExprPtr& expression, ComputationContext& context);
ExprPtr exact_subtract(const ExprPtr& left, const ExprPtr& right, ComputationContext& context);
ExprPtr exact_multiply(const ExprPtr& left, const ExprPtr& right, ComputationContext& context);
ExprPtr exact_divide(const ExprPtr& numerator, const ExprPtr& denominator, ComputationContext& context);
void swap_rows(ExactMatrixData& matrix, std::size_t first, std::size_t second);
Result<std::optional<std::size_t>> choose_pivot(
    const ExactMatrixData& matrix, std::size_t first_row, std::size_t column,
    ComputationContext& context, const std::string& operation);

}
