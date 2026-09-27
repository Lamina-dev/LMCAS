#include "residual_verification.hpp"

#include "internal/assumption_facts.hpp"
#include "internal/facts_query.hpp"
#include "internal/pointwise_comparison.hpp"
#include "internal/symbolic_ast.hpp"

#include <exception>
#include <new>
#include <optional>

namespace LMCAS {

ResidualCheckResult check_equivalent(
    const ExprPtr& left,
    const ExprPtr& right,
    ComputationContext& context,
    const LMCAS::EqvOptions& options) {
    if (!left || !right || !detail::node(left) || !detail::node(right)) {
        return ResidualCheckResult::failure(
            CasErrc::InvalidArgument,
            "equivalence check requires two expressions",
            "residual.check_equivalent");
    }
    auto access = context.consume_steps(1, "residual.check_equivalent");
    if (!access) return ResidualCheckResult::failure(access.error());
    if (detail::node(left)->equals(*detail::node(right)))
        return ResidualCheckResult::success(ProvedZeroResidual{
            ExactNormalizationProof{left}});
    try {
        auto residual = SymbolicExpr::add(
            left,
            SymbolicExpr::multiply(SymbolicExpr::number(-1), right));
        return check_zero_residual(residual, context, options);
    } catch (const CasError& error) {
        return ResidualCheckResult::failure(error);
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

namespace {

Result<bool> pointwise_values_defined(
    const ExprPtr& left, const ExprPtr& right, const FactsQuery& facts,
    Domain domain, ComputationContext& context) {
    auto left_defined = detail::query_definedness(
        detail::node(left), facts, domain, context);
    if (!left_defined) {
        return Result<bool>::failure(left_defined.error());
    }
    auto right_defined = detail::query_definedness(
        detail::node(right), facts, domain, context);
    if (!right_defined) {
        return Result<bool>::failure(right_defined.error());
    }
    return left_defined.value() == Tribool::True &&
        right_defined.value() == Tribool::True;
}

Result<Tribool> compare_defined_values(
    const ExprPtr& left, const ExprPtr& right, const FactsQuery& facts,
    ComputationContext& context, Domain domain) {
    if (detail::node(left)->equals(*detail::node(right))) {
        return Tribool::True;
    }
    auto residual = SymbolicExpr::add(
        left,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), right))->simplify();
    if (!residual) {
        return Result<Tribool>::failure(
            CasErrc::InternalInvariant,
            "value comparison produced no residual",
            "residual.compare_values");
    }
    auto nonzero = detail::query_nonzero_value(
        detail::node(residual), facts, domain, context);
    if (!nonzero) {
        return Result<Tribool>::failure(nonzero.error());
    }
    if (nonzero.value() == Tribool::True) {
        return Tribool::False;
    }
    if (nonzero.value() == Tribool::False) {
        return Tribool::True;
    }
    auto identity = check_equivalent(left, right, context);
    if (!identity) {
        return Result<Tribool>::failure(identity.error());
    }
    return std::holds_alternative<ProvedZeroResidual>(identity.value())
        ? Tribool::True : Tribool::Unknown;
}

}

Result<Tribool> detail::compare_pointwise_values(
    const ExprPtr& left, const ExprPtr& right,
    ComputationContext& context, Domain domain) {
    constexpr const char* operation = "residual.compare_values";
    auto access = context.consume_steps(1, operation);
    if (!access) {
        return Result<Tribool>::failure(access.error());
    }
    if (!left || !right || !detail::node(left) || !detail::node(right)) {
        return Result<Tribool>::failure(
            CasErrc::InvalidArgument,
            "value comparison requires two expressions", operation);
    }
    try {
        std::optional<detail::AssumptionFacts> assumed;
        if (context.assumptions()) {
            assumed.emplace(*context.assumptions());
        }
        const auto& facts = assumed
            ? static_cast<const FactsQuery&>(*assumed)
            : detail::no_facts();
        auto defined = pointwise_values_defined(
            left, right, facts, domain, context);
        if (!defined) {
            return Result<Tribool>::failure(defined.error());
        }
        if (!defined.value()) {
            return Tribool::Unknown;
        }
        return compare_defined_values(
            left, right, facts, context, domain);
    } catch (const CasError& error) {
        return Result<Tribool>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<Tribool>::failure(
            CasErrc::ResourceLimit,
            "value comparison allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<Tribool>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

}
