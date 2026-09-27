#pragma once

#include "computation_context.hpp"
#include "equivalence_options.hpp"
#include "proof_outcome.hpp"
#include "result.hpp"
#include "symbolic.hpp"

#include <variant>

namespace LMCAS {


/** @brief 在表达式的原公共定义域上认证值相等；点值定义域与分母条件需另行证明。 */
struct ProvedZeroResidual {
    ProofCertificate certificate;
};

/**
 * @brief 证明残差在原公共定义域上非恒零，例如非零有理多项式。
 * 残差仍可能在特定参数值处为零；点值不等性、定义域及上下文假设下的非零性均需另行证明。
 */
struct ProvedNonIdentityResidual {
    ExprPtr normalized_residual;
};

struct UnprovedResidual {
    ExprPtr normalized_residual;
};

using ResidualCheck = std::variant<
    ProvedZeroResidual,
    ProvedNonIdentityResidual,
    UnprovedResidual>;

using ResidualCheckResult = Result<ResidualCheck>;

/**
 * @brief 成功时返回恒等、非恒等或未知分类，计算错误由外层 Result 返回。
 * 假设可辅助恒等规则；形式非恒等与点值事实分别判定，未知分类仍属成功返回。
 */

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
