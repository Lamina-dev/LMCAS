#pragma once
#include "multivariate_factor.hpp"
#include <optional>

namespace LMCAS::multivariate_checked_detail {
MultiFactorCheckedResult assemble_checked_factorization(
    const MultiPoly& original, const std::vector<MultiPoly>& factors,
    const std::vector<int>& multiplicities, Completeness completeness,
    std::string reason);
std::optional<MultiFactorCheckedResult> factor_evaluated_linear(
    const MultiPoly& original, const MultiPoly& primitive,
    const std::string& main_variable, ComputationContext& context);
}
