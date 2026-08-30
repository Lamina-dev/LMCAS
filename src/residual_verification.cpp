#include "residual_verification.hpp"

#include "internal/exact_matrix.hpp"
#include "internal/expression_analysis.hpp"
#include "poly_utils.hpp"

#include <array>
#include <new>

namespace lamina {
namespace {

ExprPtr normalized(ExprPtr expression) {
    if (!expression) return expression;
    auto simplified = expression->simplify();
    return simplified ? simplified : expression;
}

ResidualCheckResult classify_residual(
    const ExprPtr& expression,
    ComputationContext& context,
    const std::string& operation,
    bool final_stage) {
    auto proof = detail::classify_exact_zero(
        expression, context, operation);
    if (!proof) return ResidualCheckResult::failure(proof.error());
    if (proof.value() == detail::ZeroProof::Zero) {
        return ResidualCheckResult::success(ProvedZeroResidual{
            ExactNormalizationProof{expression}});
    }
    const auto variables = free_variables(detail::node(expression));
    if (variables.size() == 1) {
        auto polynomial = recognize_rational_polynomial(
            *expression, *variables.begin(), context);
        if (!polynomial) {
            return ResidualCheckResult::failure(polynomial.error());
        }
        if (polynomial.value() && !polynomial.value()->is_zero()) {
            return ResidualCheckResult::success(
                ProvedNonzeroResidual{expression});
        }
    }
    if (proof.value() == detail::ZeroProof::NonZero) {
        return ResidualCheckResult::success(
            ProvedNonzeroResidual{expression});
    }
    if (final_stage) {
        return ResidualCheckResult::success(UnprovedResidual{expression});
    }
    return ResidualCheckResult::success(UnprovedResidual{expression});
}

bool is_decided(const ResidualCheck& result) {
    return !std::holds_alternative<UnprovedResidual>(result);
}

} // namespace

ResidualCheckResult check_zero_residual(
    const ExprPtr& residual,
    ComputationContext& context,
    const lsr::EqvOptions& options) {
    constexpr const char* operation = "residual.check_zero";
    if (!residual || !detail::node(residual)) {
        return ResidualCheckResult::failure(
            CasErrc::InvalidArgument,
            "residual expression cannot be null", operation);
    }
    auto access = context.consume_steps(1, operation);
    if (!access) return ResidualCheckResult::failure(access.error());

    try {
        auto current = normalized(residual);
        auto first = classify_residual(current, context, operation, false);
        if (!first || is_decided(first.value())) return first;

        current = normalized(current->expand());
        auto expanded = classify_residual(current, context, operation, false);
        if (!expanded || is_decided(expanded.value())) return expanded;

        current = normalized(current->cancel());
        auto cancelled = classify_residual(current, context, operation, false);
        if (!cancelled || is_decided(cancelled.value())) return cancelled;

        current = normalized(current->expand());
        current = normalized(current->cancel());
        auto combined = classify_residual(current, context, operation, false);
        if (!combined || is_decided(combined.value())) return combined;

        const auto zero = SymbolicExpr::number(0);
        for (const auto profile : {
                 lsr::EqvProfile::Core,
                 lsr::EqvProfile::TrigBasic,
                 lsr::EqvProfile::ExpLogBasic}) {
            auto profile_options = options;
            profile_options.profile = profile;
            auto equivalent = lsr::equivalent_core(
                *current, *zero, context, profile_options);
            if (!equivalent) {
                if (equivalent.error().code == CasErrc::UnsupportedExpression ||
                    equivalent.error().code == CasErrc::Inconclusive) {
                    continue;
                }
                return ResidualCheckResult::failure(equivalent.error());
            }
            if (equivalent.value()) {
                return ResidualCheckResult::success(ProvedZeroResidual{
                    ExactResidualProof{current}});
            }
        }
        return ResidualCheckResult::success(UnprovedResidual{current});
    } catch (const std::bad_alloc&) {
        return ResidualCheckResult::failure(
            CasErrc::ResourceLimit,
            "residual verification allocation failed", operation);
    } catch (const std::exception& error) {
        return ResidualCheckResult::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

ResidualCheckResult check_equivalent(
    const ExprPtr& left,
    const ExprPtr& right,
    ComputationContext& context,
    const lsr::EqvOptions& options) {
    if (!left || !right || !detail::node(left) || !detail::node(right)) {
        return ResidualCheckResult::failure(
            CasErrc::InvalidArgument,
            "equivalence check requires two expressions",
            "residual.check_equivalent");
    }
    try {
        auto residual = SymbolicExpr::add(
            left,
            SymbolicExpr::multiply(SymbolicExpr::number(-1), right));
        return check_zero_residual(residual, context, options);
    } catch (const std::bad_alloc&) {
        return ResidualCheckResult::failure(
            CasErrc::ResourceLimit,
            "equivalence residual allocation failed",
            "residual.check_equivalent");
    } catch (const std::exception& error) {
        return ResidualCheckResult::failure(
            CasErrc::InternalInvariant, error.what(),
            "residual.check_equivalent");
    }
}

} // namespace lamina
