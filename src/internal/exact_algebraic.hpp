#pragma once

#include "computation_context.hpp"
#include "polynomial.hpp"
#include "result.hpp"

#include <cstddef>

namespace LMCAS::detail {

struct ExactRealAlgebraic {
    Polynomial<Rational> polynomial;
    Rational lower;
    Rational upper;
    std::size_t root_index = 0;
    std::size_t multiplicity = 1;

    bool is_rational() const noexcept { return lower == upper; }
};

using ExactRealAlgebraicResult = Result<ExactRealAlgebraic>;

LMCAS_API ExactRealAlgebraicResult make_exact_real_algebraic(
    Polynomial<Rational> polynomial,
    std::size_t root_index,
    std::size_t multiplicity,
    ComputationContext& context);

LMCAS_API Result<void> refine_exact_real_algebraic(
    ExactRealAlgebraic& value,
    ComputationContext& context,
    const std::string& operation = "exact_algebraic.refine");

LMCAS_API Result<void> refine_exact_real_algebraic_to_tolerance(
    ExactRealAlgebraic& value,
    double absolute_tolerance,
    double relative_tolerance,
    ComputationContext& context,
    const std::string& operation = "exact_algebraic.refine");

LMCAS_API Result<bool> equal_exact_real_algebraic(
    ExactRealAlgebraic lhs,
    ExactRealAlgebraic rhs,
    ComputationContext& context);

LMCAS_API Result<int> compare_exact_real_algebraic(
    ExactRealAlgebraic lhs,
    ExactRealAlgebraic rhs,
    ComputationContext& context);

} // namespace LMCAS::detail
