#include "solve_transcendental.hpp"
#include "internal/assumption_simplification.hpp"
#include "internal/transcendental_solver_support.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/assumption_facts.hpp"
#include "poly_utils.hpp"
#include "solve_polynomial.hpp"
#include "solve_strategies.hpp"
#include "residual_verification.hpp"
#include <algorithm>

namespace LMCAS {
namespace {

using namespace detail::transcendental;
constexpr int maximum_transcendental_depth = 5;

SolveResult inconclusive(const char* message) {
    return SolveResult::failure(CasErrc::Inconclusive, message, "solve.transcendental");
}

SolveResult solve_transcendental_impl(
    const ExprPtr& original, const ExprPtr& expr, const std::string& var,
    ComputationContext& context, const SolveOptions& options, int depth);

SolveResult solve_inner_equation(
    const ExprPtr& original, const ExprPtr& inner, const ExprPtr& value,
    const std::string& var, ComputationContext& context,
    const SolveOptions& options, int depth) {
    auto variable = std::dynamic_pointer_cast<const VariableNode>(detail::node(inner));
    if (variable && !variable->is_constant() && variable->name() == var) {
        return SolutionSet{FiniteSolutions{{FiniteSolution{value, 1, {}}}}};
    }
    const auto product =
        std::dynamic_pointer_cast<const MultiplyNode>(detail::node(inner));
    if (product && product->operands().size() == 2) {
        const auto& first = product->operands()[0];
        const auto& second = product->operands()[1];
        const auto coefficient =
            std::dynamic_pointer_cast<const NumberNode>(first);
        const auto product_variable =
            std::dynamic_pointer_cast<const VariableNode>(second);
        if (coefficient && coefficient->is_negative_one() &&
            product_variable && !product_variable->is_constant() &&
            product_variable->name() == var) {
            auto inverted = SymbolicExpr::multiply(
                SymbolicExpr::number(-1), value)->simplify();
            return SolutionSet{FiniteSolutions{{
                FiniteSolution{std::move(inverted), 1, {}}}}};
        }
    }
    auto equation = SymbolicExpr::add(inner,
        SymbolicExpr::multiply(value, SymbolicExpr::number(-1)));
    auto polynomial = symbolic_to_poly<SymbolicPolyCoeff>(equation, var);
    if (polynomial) return solve_equation(equation, var, context, options);
    if (polynomial.error().code != CasErrc::UnsupportedExpression) {
        return SolveResult::failure(polynomial.error());
    }
    return solve_transcendental_impl(original, equation, var, context, options, depth + 1);
}

void append_finite_union(FiniteSolutions& finite, FiniteSolutions& prior) {
    for (auto& candidate : finite.values) {
        auto duplicate = std::find_if(prior.values.begin(), prior.values.end(),
            [&](const FiniteSolution& previous) {
                return candidate.conditions.empty() && previous.conditions.empty() &&
                    detail::node(candidate.value)->equals(*detail::node(previous.value));
            });
        if (duplicate == prior.values.end()) { prior.values.push_back(std::move(candidate)); }
    }
}

void append_parametric_union(FiniteSolutions* finite, ParametricSolutions* families,
    FiniteSolutions* prior, SolutionSet& output) {
    if (prior) {
        ParametricSolutions promoted;
        for (auto& point : prior->values) {
            promoted.values.push_back({std::move(point.value), {}, std::move(point.conditions)});
        }
        output = std::move(promoted);
    }
    auto& result = std::get<ParametricSolutions>(output).values;
    if (finite) {
        for (auto& point : finite->values) {
            result.push_back({std::move(point.value), {}, std::move(point.conditions)});
        }
    } else {
        for (auto& family : families->values) { result.push_back(std::move(family)); }
    }
}
void append_solutions(SolutionSet solved, SolutionSet& output,
                      const ExprPtr& original, const std::string& var) {
    if (std::holds_alternative<EmptySolutions>(solved) ||
        std::holds_alternative<UniversalSolutions>(output)) { return; }
    if (std::holds_alternative<EmptySolutions>(output) ||
        std::holds_alternative<UniversalSolutions>(solved)) {
        output = std::move(solved);
        return;
    }
    auto* finite = std::get_if<FiniteSolutions>(&solved);
    auto* prior = std::get_if<FiniteSolutions>(&output);
    if (finite && prior) {
        append_finite_union(*finite, *prior);
        return;
    }
    auto* families = std::get_if<ParametricSolutions>(&solved);
    if ((finite || families) && (prior || std::holds_alternative<ParametricSolutions>(output))) {
        append_parametric_union(finite, families, prior, output);
        return;
    }
    output = conditional_preimage(original, var);
}

Result<bool> zero_difference(const ExprPtr& left, const ExprPtr& right, ComputationContext& context) {
    auto difference = detail::simplify_expression(SymbolicExpr::add(left,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), right)), context);
    if (!difference) { return Result<bool>::failure(difference.error()); }
    auto proof = check_zero_residual(difference.value(), context);
    if (!proof) { return Result<bool>::failure(proof.error()); }
    return std::holds_alternative<ProvedZeroResidual>(proof.value());
}

bool compatible_periodic_conditions(const ParametricSolution& left, const ParametricSolution& right) {
    if (right.integer_parameters.size() != 1 || left.conditions.size() != right.conditions.size()) {
        return false;
    }
    for (std::size_t c = 0; c < left.conditions.size(); ++c) {
        if (expression_depends_on_variable(detail::node(left.conditions[c]),
                                          left.integer_parameters.front()) ||
            expression_depends_on_variable(detail::node(right.conditions[c]),
                                          right.integer_parameters.front()) ||
            !detail::node(left.conditions[c])->equals(*detail::node(right.conditions[c]))) {
            return false;
        }
    }
    return true;
}

Result<std::optional<ExprPtr>> merged_half_period_value(const detail::AffineForm& first,
    const detail::AffineForm& second, const std::string& parameter, ComputationContext& context) {
    using MergeResult = Result<std::optional<ExprPtr>>;
    auto half = detail::simplify_expression(SymbolicExpr::multiply(
        SymbolicExpr::number(Rational(1, 2)), first.slope), context);
    if (!half) { return MergeResult::failure(half.error()); }
    const auto& half_period = half.value();
    auto next_offset = SymbolicExpr::add(first.offset, half_period);
    auto prior_offset = SymbolicExpr::add(second.offset, half_period);
    ExprPtr offset;
    auto forward = zero_difference(next_offset, second.offset, context);
    if (!forward) { return MergeResult::failure(forward.error()); }
    auto backward = zero_difference(prior_offset, first.offset, context);
    if (!backward) { return MergeResult::failure(backward.error()); }
    if (forward.value()) { offset = second.offset; }
    else if (backward.value()) { offset = first.offset; }
    else { return std::optional<ExprPtr>{}; }
    auto merged = detail::simplify_expression(SymbolicExpr::add(offset,
        SymbolicExpr::multiply(half_period, SymbolicExpr::variable(parameter))), context);
    if (!merged) { return MergeResult::failure(merged.error()); }
    return std::optional<ExprPtr>{std::move(merged.value())};
}

Result<bool> coalesce_periodic_pair(ParametricSolution& left, const ParametricSolution& right,
    ComputationContext& context) {
    if (!compatible_periodic_conditions(left, right)) { return false; }
    auto first = detail::recognize_affine(*left.value, left.integer_parameters.front(), context);
    if (!first) { return Result<bool>::failure(first.error()); }
    auto second = detail::recognize_affine(*right.value, right.integer_parameters.front(), context);
    if (!second) { return Result<bool>::failure(second.error()); }
    if (!first.value() || !second.value()) { return false; }
    for (auto* coefficient : {&first.value()->slope, &first.value()->offset,
                              &second.value()->slope, &second.value()->offset}) {
        auto normalized = detail::simplify_expression(*coefficient, context);
        if (!normalized) { return Result<bool>::failure(normalized.error()); }
        *coefficient = std::move(normalized.value());
    }
    auto equal_slopes = zero_difference(first.value()->slope, second.value()->slope, context);
    if (!equal_slopes) { return Result<bool>::failure(equal_slopes.error()); }
    if (!equal_slopes.value()) { return false; }
    auto merged = merged_half_period_value(
        *first.value(), *second.value(), left.integer_parameters.front(), context);
    if (!merged) { return Result<bool>::failure(merged.error()); }
    if (!merged.value()) { return false; }
    left.value = std::move(*merged.value());
    return true;
}

Result<void> coalesce_periodic_families(SolutionSet& solutions, ComputationContext& context) {
    auto* parametric = std::get_if<ParametricSolutions>(&solutions);
    if (!parametric) { return Result<void>::success(); }
    auto& families = parametric->values;
    for (std::size_t i = 0; i < families.size(); ++i) {
        if (families[i].integer_parameters.size() != 1) { continue; }
        for (std::size_t j = i + 1; j < families.size();) {
            auto merged = coalesce_periodic_pair(families[i], families[j], context);
            if (!merged) { return Result<void>::failure(merged.error()); }
            if (!merged.value()) { ++j; continue; }
            families.erase(families.begin() + j);
            j = i + 1;
        }
    }
    return Result<void>::success();
}

SolveResult solve_inverted_values(
    const ExprPtr& original, const ExprPtr& inner,
    const std::vector<ExprPtr>& values, const std::string& var,
    ComputationContext& context, const SolveOptions& options, int depth) {
    if (values.empty()) return inconclusive("Inverse did not establish a complete value set");
    SolutionSet results = EmptySolutions{};
    for (const auto& value : values) {
        auto solved = solve_inner_equation(original, inner, value, var, context, options, depth);
        if (!solved) return solved;
        append_solutions(std::move(solved).value(), results, original, var);
    }
    return results;
}

SolveResult solve_substitution(
    const ExprPtr& original, const SubstitutionResult& substitution, const std::string& var,
    ComputationContext& context, const SolveOptions& options, int depth) {
    auto polynomial = symbolic_to_poly<SymbolicPolyCoeff>(
        substitution.poly_in_u, substitution.u_var);
    if (!polynomial) {
        if (polynomial.error().code == CasErrc::UnsupportedExpression) {
            return inconclusive("Substitution is not a supported polynomial");
        }
        return SolveResult::failure(polynomial.error());
    }
    auto solved = solve_equation(substitution.poly_in_u, substitution.u_var, context, options);
    if (!solved) return solved;
    if (std::holds_alternative<EmptySolutions>(solved.value())) return solved;
    auto finite = std::get_if<FiniteSolutions>(&solved.value());
    if (!finite) return conditional_preimage(original, var);
    SolutionSet results = EmptySolutions{};
    for (const auto& value : finite->values) {
        if (!value.conditions.empty()) return conditional_preimage(original, var);
        auto inverted = solve_inner_equation(
            original, substitution.u_expr, value.value, var, context, options, depth);
        if (!inverted) return inverted;
        append_solutions(std::move(inverted).value(), results, original, var);
    }
    return results;
}

SolveResult solve_transcendental_constant(const ExprPtr& original, const ExprPtr& expr,
    const std::string& var, ComputationContext& context) {
    std::optional<detail::AssumptionFacts> assumptions;
    if (context.assumptions()) { assumptions.emplace(*context.assumptions()); }
    const LMCAS::FactsQuery& facts = assumptions ?
        static_cast<const LMCAS::FactsQuery&>(*assumptions) : detail::no_facts();
    auto nonzero = detail::query_nonzero_value(detail::node(expr->simplify()), facts, Domain::Real, context);
    if (!nonzero) { return SolveResult::failure(nonzero.error()); }
    if (nonzero.value() == Tribool::True) { return SolutionSet{EmptySolutions{}}; }
    if (nonzero.value() == Tribool::False) { return SolutionSet{UniversalSolutions{}}; }
    return conditional_preimage(original, var);
}

SolveResult solve_function_pattern(const InversePattern& pattern, const ExprPtr& original,
    const std::string& var, ComputationContext& context, const SolveOptions& options, int depth) {
    auto inverted = invert_function(pattern, original, var, context, options);
    if (!inverted) { return inverted; }
    auto finite = std::get_if<FiniteSolutions>(&inverted.value());
    if (!finite) { return inverted; }
    SolutionSet results = EmptySolutions{};
    for (const auto& value : finite->values) {
        if (!value.conditions.empty()) { return conditional_preimage(original, var); }
        auto solved = solve_inner_equation(original, pattern.inner, value.value,
            var, context, options, depth);
        if (!solved) { return solved; }
        append_solutions(std::move(solved).value(), results, original, var);
    }
    return results;
}

SolveResult solve_transcendental_impl(
    const ExprPtr& original, const ExprPtr& expr, const std::string& var,
    ComputationContext& context, const SolveOptions& options, int depth) {
    auto step = context.consume_steps(1, "solve.transcendental");
    if (!step) { return SolveResult::failure(step.error()); }
    if (depth > maximum_transcendental_depth) { return inconclusive("Transcendental inversion depth exceeded"); }
    if (!expression_depends_on_variable(detail::node(expr), var)) {
        return solve_transcendental_constant(original, expr, var, context);
    }
    if (auto pattern = decompose_lambert_w_pattern(expr, var)) {
        return solve_inverted_values(original, pattern->inner, invert_lambert_w(pattern->rhs),
            var, context, options, depth);
    }
    if (auto pattern = decompose_exp_base_pattern(expr, var)) {
        return solve_inverted_values(original, pattern->inner,
            invert_exp_base(pattern->base, pattern->rhs), var, context, options, depth);
    }
    if (auto pattern = decompose_trig_exp_pattern(expr, var)) {
        return solve_function_pattern(*pattern, original, var, context, options, depth);
    }
    auto substitution = detect_substitution(expr, var);
    if (!substitution) {
        return SolveResult::failure(substitution.error());
    }
    if (substitution.value()) {
        return solve_substitution(original, *substitution.value(), var, context, options, depth);
    }
    return inconclusive("No complete transcendental inversion matched");
}

}

SolveResult solve_transcendental(
    const ExprPtr& expr, const std::string& var,
    ComputationContext& context, const SolveOptions& options) {
    if (!expr || !detail::node(expr) || var.empty()) {
        return SolveResult::failure(CasErrc::InvalidArgument,
            "Expression and variable must not be empty", "solve.transcendental");
    }
    auto equation = expr;
    if (auto relation = std::dynamic_pointer_cast<const RelationalNode>(detail::node(expr))) {
        if (relation->op() != RelationalNode::Op::EQ) {
            return SolveResult::failure(CasErrc::InvalidArgument,
                "Expected an equality", "solve.transcendental");
        }
        equation = SymbolicExpr::add(detail::make_expression_ptr(relation->left()),
            SymbolicExpr::multiply(SymbolicExpr::number(-1),
                detail::make_expression_ptr(relation->right())));
    }
    equation = equation->simplify();
    auto solved = std::dynamic_pointer_cast<const MultiplyNode>(detail::node(equation)) ?
        solve_equation(equation, var, context, options) :
        solve_transcendental_impl(expr, equation, var, context, options, 0);
    if (!solved) return solved;
    auto coalesced = coalesce_periodic_families(solved.value(), context);
    if (!coalesced) return SolveResult::failure(coalesced.error());
    return detail::finalize_solution_set(expr, var, std::move(solved).value(), context, Domain::Real, options);
}

SolveResult solve_transcendental(
    const ExprPtr& expr, const std::string& var, const SolveOptions& options) {
    ComputationContext context;
    return solve_transcendental(expr, var, context, options);
}

}
