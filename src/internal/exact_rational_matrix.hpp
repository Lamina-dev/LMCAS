#pragma once

#include "computation_context.hpp"
#include "rational.hpp"
#include "result.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace LMCAS::detail {
Result<Rational> rational_determinant_value_impl(
    std::size_t dimension,
    std::vector<Rational> values,
    ComputationContext& context,
    const std::string& operation);


Result<Rational> rational_determinant_exact(
    std::size_t dimension,
    std::vector<Rational> values,
    ComputationContext& context,
    const std::string& operation);

Result<std::vector<Rational>> solve_rational_unique(
    std::size_t rows,
    std::size_t coefficient_columns,
    std::vector<Rational> augmented_entries,
    ComputationContext& context,
    const std::string& operation);
std::vector<std::optional<std::size_t>> reduce_rational_rows(
    std::size_t rows, std::size_t coefficient_columns,
    std::vector<Rational>& matrix, bool start_at_pivot);

bool rational_system_inconsistent(
    std::size_t rows, std::size_t coefficient_columns,
    const std::vector<Rational>& matrix);

}
