#include "limit_result.hpp"
#include "internal/calculus_utils_support.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/facts_query.hpp"
#include "internal/pointwise_comparison.hpp"
#include "assumption_context.hpp"

namespace LMCAS {
using namespace calculus_utils_detail;

namespace {
Result<detail::SymbolicNodePtr> select_continuity_branch(
    const PiecewiseNode& piecewise, const FactsQuery& facts,
    const AssumptionContext& assumptions, ComputationContext& context) {
    for (const auto& branch : piecewise.branches()) {
        auto access = context.consume_steps(1, "continuity_at");
        if (!access) { return Result<detail::SymbolicNodePtr>::failure(access.error()); }
        auto defined = detail::query_definedness(branch.condition, facts, Domain::Real, context);
        if (!defined) { return Result<detail::SymbolicNodePtr>::failure(defined.error()); }
        const auto truth = defined.value() == Tribool::True
            ? assumptions.evaluate_condition(detail::expression_from_node(branch.condition))
            : Tribool::Unknown;
        if (truth == Tribool::Unknown) {
            return Result<detail::SymbolicNodePtr>::failure(CasErrc::Inconclusive,
                "point-value branch is undecided", "continuity_at");
        }
        if (truth == Tribool::True) { return branch.expression; }
    }
    return piecewise.default_expr();
}

Result<ExprPtr> continuity_point_value(
    ExprPtr value, const FactsQuery& facts, ComputationContext& context) {
    const AssumptionContext empty_assumptions;
    const auto& assumptions = context.assumptions() ? *context.assumptions() : empty_assumptions;
    while (const auto piecewise = std::dynamic_pointer_cast<const PiecewiseNode>(detail::node(value))) {
        auto selected = select_continuity_branch(*piecewise, facts, assumptions, context);
        if (!selected) { return Result<ExprPtr>::failure(selected.error()); }
        if (!selected.value()) { return ExprPtr{}; }
        value = detail::make_expression_ptr(selected.value());
    }
    return value;
}

ContinuityResult classify_continuity_point_value(
    const ExprPtr& limit, const ExprPtr& value,
    const FactsQuery& facts, ComputationContext& context) {
    auto defined = detail::query_definedness(detail::node(value), facts, Domain::Real, context);
    if (!defined) { return ContinuityResult::failure(defined.error()); }
    if (defined.value() == Tribool::False) {
        return ContinuityResult::success(ContinuityType::Removable);
    }
    if (defined.value() != Tribool::True) {
        return ContinuityResult::failure(CasErrc::Inconclusive,
            "point-value definedness is unknown", "continuity_at");
    }
    auto same = detail::compare_pointwise_values(limit, value, context, Domain::Real);
    if (!same) { return ContinuityResult::failure(same.error()); }
    if (same.value() == Tribool::True) {
        return ContinuityResult::success(ContinuityType::Continuous);
    }
    if (same.value() == Tribool::False) {
        return ContinuityResult::success(ContinuityType::Removable);
    }
    return ContinuityResult::failure(CasErrc::Inconclusive,
        "point value could not be compared with the limit", "continuity_at");
}

ContinuityResult continuity_with_limit(
    const ExprPtr& f, const std::string& var, const ExprPtr& point,
    const ExprPtr& limit, ComputationContext& context) {
    auto value = f->substitute(var, point);
    if (!value) {
        return ContinuityResult::failure(CasErrc::InternalInvariant,
            "substitution produced no value", "continuity_at");
    }
    std::optional<detail::AssumptionFacts> assumed;
    if (context.assumptions()) { assumed.emplace(*context.assumptions()); }
    const auto& facts = assumed ? static_cast<const FactsQuery&>(*assumed) : detail::no_facts();
    auto selected = continuity_point_value(value, facts, context);
    if (!selected) { return ContinuityResult::failure(selected.error()); }
    if (!selected.value()) { return ContinuityResult::success(ContinuityType::Removable); }
    return classify_continuity_point_value(limit, selected.value(), facts, context);
}

ContinuityResult continuity_from_limits(
    const ExprPtr& f, const std::string& var, const ExprPtr& point,
    const LimitOutcome& left, const LimitOutcome& right, ComputationContext& context) {
    const auto* left_finite = std::get_if<FiniteLimit>(&left);
    const auto* right_finite = std::get_if<FiniteLimit>(&right);
    if (!left_finite || !right_finite) {
        return ContinuityResult::success(ContinuityType::Essential);
    }
    auto sides = detail::compare_pointwise_values(
        left_finite->value, right_finite->value, context, Domain::Real);
    if (!sides) { return ContinuityResult::failure(sides.error()); }
    if (sides.value() == Tribool::False) {
        return ContinuityResult::success(ContinuityType::Jump);
    }
    if (sides.value() != Tribool::True) {
        return ContinuityResult::failure(CasErrc::Inconclusive,
            "one-sided values could not be compared", "continuity_at");
    }
    return continuity_with_limit(f, var, point, left_finite->value, context);
}
}

ContinuityResult continuity_at_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point, ComputationContext& context) {
    const std::string operation = "continuity_at";
    try {
        if (!f || !point || var.empty()) {
            return ContinuityResult::failure(CasErrc::InvalidArgument,
                "continuity analysis requires an expression, point, and variable", operation);
        }
        auto input = calculus_utils_validate_expr(f, var, context, operation);
        if (!input) { return ContinuityResult::failure(input.error()); }
        auto left = limit_checked(f, var, point, LimitDirection::FromBelow, context);
        if (!left) { return ContinuityResult::failure(left.error()); }
        auto right = limit_checked(f, var, point, LimitDirection::FromAbove, context);
        if (!right) { return ContinuityResult::failure(right.error()); }
        return continuity_from_limits(f, var, point, left.value().value, right.value().value, context);
    } catch (const std::bad_alloc&) {
        return ContinuityResult::failure(CasErrc::ResourceLimit, "continuity analysis allocation failed", operation);
    } catch (const std::exception& ex) {
        return ContinuityResult::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
}

ContinuityResult continuity_at_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point) {
    ComputationContext context;
    return continuity_at_checked(f, var, point, context);
}
}
