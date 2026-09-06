#pragma once

#include "computation_context.hpp"
#include "polynomial.hpp"
#include "result.hpp"

#include <string>
#include <utility>
#include <vector>

namespace LMCAS::detail {

using RationalInterval = std::pair<Rational, Rational>;

Result<std::vector<RationalInterval>> isolate_real_roots_exact(
    const Polynomial<Rational>& polynomial,
    ComputationContext& context,
    const std::string& operation);

Result<std::size_t> count_real_roots_exact(
    const Polynomial<Rational>& polynomial,
    const Rational& lower,
    const Rational& upper,
    ComputationContext& context,
    const std::string& operation);

} // namespace LMCAS::detail
