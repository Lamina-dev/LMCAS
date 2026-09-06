#pragma once

#include "computation_context.hpp"
#include "polynomial.hpp"
#include "result.hpp"

namespace LMCAS {

using UnivariateFactorResult =
    Result<MathResult<std::vector<Polynomial<Rational>>>>;

UnivariateFactorResult factor_univariate_bridge_checked(
    const Polynomial<Rational>& polynomial,
    ComputationContext& context);

} // namespace LMCAS
