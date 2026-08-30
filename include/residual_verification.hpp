#pragma once

#include "computation_context.hpp"
#include "lsr_expr.hpp"
#include "proof_outcome.hpp"
#include "result.hpp"
#include "symbolic.hpp"

#include <variant>

namespace lamina {

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

LAMINA_API ResidualCheckResult check_zero_residual(
    const ExprPtr& residual,
    ComputationContext& context,
    const lsr::EqvOptions& options = {});

LAMINA_API ResidualCheckResult check_equivalent(
    const ExprPtr& left,
    const ExprPtr& right,
    ComputationContext& context,
    const lsr::EqvOptions& options = {});

} // namespace lamina
