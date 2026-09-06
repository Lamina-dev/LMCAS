#pragma once

#include "computation_context.hpp"
#include "expr.hpp"
#include "proof_outcome.hpp"
#include "result.hpp"
#include "symbolic.hpp"

#include <variant>

namespace LMCAS {

struct ProvedZeroResidual {
    ProofCertificate certificate;
};

struct ProvedNonzeroResidual {
    ExprPtr normalized_residual;
};

struct UnprovedResidual {
    ExprPtr normalized_residual;
};

using ResidualCheck = std::variant<
    ProvedZeroResidual,
    ProvedNonzeroResidual,
    UnprovedResidual>;

using ResidualCheckResult = Result<ResidualCheck>;

LMCAS_API ResidualCheckResult check_zero_residual(
    const ExprPtr& residual,
    ComputationContext& context,
    const LMCAS::EqvOptions& options = {});

LMCAS_API ResidualCheckResult check_equivalent(
    const ExprPtr& left,
    const ExprPtr& right,
    ComputationContext& context,
    const LMCAS::EqvOptions& options = {});

} // namespace LMCAS
