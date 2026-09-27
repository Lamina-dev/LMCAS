#include "symbolic.hpp"
#include "parametric_solver.hpp"
#include "assumption_context.hpp"
#include "solver.hpp"
#include <iostream>
#include <map>
#include <set>
#include <cmath>
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "internal/visitors/normalization_visitor.hpp"
#include "poly_utils.hpp"
#include "solve_strategies.hpp"
#include "solve_polynomial.hpp"
#include "solve_transcendental.hpp"
#include "solve_mixed_transcendental.hpp"
#include "newton_raphson.hpp"
#include "root_of_utils.hpp"
#include "rational_polynomial.hpp"
#include "residual_verification.hpp"
#include "internal/pointwise_comparison.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/complex_quadratic.hpp"
#include "internal/transcendental_solver_support.hpp"
#include "internal/normalization_utils.hpp"
#include <limits>
#include "internal/polynomial_solver_support.hpp"

namespace LMCAS {

using SolveVectorResult =
    LMCAS::Result<std::vector<std::shared_ptr<SymbolicExpr>>>;
using detail::equation_predicate;
using detail::substitute_raw;

std::shared_ptr<SymbolicExpr> get_coeff(const Polynomial<SymbolicPolyCoeff>& p, int deg) {
    if (deg < 0 || deg > p.degree()) {
        return SymbolicExpr::number(0);
    }
    return p.coeffs[deg].val;
}

static SolveResult solve_set_core(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts);


static SolveResult conditional_equation(const ExprPtr& original, const std::string& variable) {
    return SolutionSet{ConditionalSolutions{ConditionSet{variable, equation_predicate(original), {}}}};
}

class ApproximateSolverInput final : public detail::RecursiveSymbolicVisitor {
public:
    bool found = false;
    void visit(const NumberNode& number) override {
        found = found || std::holds_alternative<lmmc_real_t>(number.value());
    }
};

static Result<void> validate_symbolic_polynomial_input(
    const ExprPtr& original, const SolveOptions& options) {
    if (options.allow_numeric) {
        return Result<void>::success();
    }
    ApproximateSolverInput approximate;
    detail::node(original)->accept(approximate);
    if (!approximate.found) {
        return Result<void>::success();
    }
    return Result<void>::failure(
        CasErrc::Inconclusive,
        "approximate polynomial input requires allow_numeric",
        "solve.polynomial");
}

static ExprPtr polynomial_relation(
    const ExprPtr& left, const ExprPtr& right, RelationOp op) {
    return detail::make_expression_ptr(
        detail::make_node<RelationalNode>(
            detail::node(left), detail::node(right), op));
}

static ExprPtr polynomial_logic(
    const ExprPtr& left, const ExprPtr& right, LogicalNode::Op op) {
    return detail::make_expression_ptr(
        detail::make_node<LogicalNode>(
            detail::node(left), detail::node(right), op));
}

static SolveResult conditional_linear_polynomial(
    const Polynomial<SymbolicPolyCoeff>& polynomial,
    const ExprPtr& leading, const std::string& variable) {
    auto constant = get_coeff(polynomial, 0);
    auto zero = SymbolicExpr::number(0);
    auto generic = polynomial_logic(
        polynomial_relation(leading, zero, RelationOp::NEQ),
        polynomial_relation(
            SymbolicExpr::multiply(
                leading, SymbolicExpr::variable(variable)),
            SymbolicExpr::multiply(
                SymbolicExpr::number(-1), constant),
            RelationOp::EQ),
        LogicalNode::Op::And);
    auto degenerate = polynomial_logic(
        polynomial_relation(leading, zero, RelationOp::EQ),
        polynomial_relation(constant, zero, RelationOp::EQ),
        LogicalNode::Op::And);
    return SolutionSet{ConditionalSolutions{ConditionSet{
        variable,
        polynomial_logic(generic, degenerate, LogicalNode::Op::Or),
        {}}}};
}

static SolveResult solve_symbolic_quadratic(
    const Polynomial<SymbolicPolyCoeff>& polynomial,
    const ExprPtr& leading, const FactsQuery& facts,
    ComputationContext& context) {
    auto linear = get_coeff(polynomial, 1);
    auto constant = get_coeff(polynomial, 0);
    auto discriminant = SymbolicExpr::add(
        SymbolicExpr::power(linear, SymbolicExpr::number(2)),
        SymbolicExpr::multiply(
            SymbolicExpr::number(-4),
            SymbolicExpr::multiply(leading, constant)))->simplify();
    auto nonzero = detail::query_nonzero_value(
        detail::node(discriminant), facts, Domain::Complex, context);
    if (!nonzero) {
        return SolveResult::failure(nonzero.error());
    }

    auto radical = detail::make_expression_ptr(
        detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sqrt,
            std::vector<std::shared_ptr<const SymbolicNode>>{
                detail::node(discriminant)}));
    auto negative_linear =
        SymbolicExpr::multiply(SymbolicExpr::number(-1), linear);
    auto denominator =
        SymbolicExpr::multiply(SymbolicExpr::number(2), leading);
    auto positive = SymbolicExpr::divide(
        SymbolicExpr::add(negative_linear, radical),
        denominator)->simplify();
    auto negative = SymbolicExpr::divide(
        SymbolicExpr::add(
            negative_linear,
            SymbolicExpr::multiply(
                SymbolicExpr::number(-1), radical)),
        denominator)->simplify();
    auto repeated =
        SymbolicExpr::divide(negative_linear, denominator)->simplify();

    if (nonzero.value() == Tribool::False) {
        return SolutionSet{FiniteSolutions{{
            FiniteSolution{repeated, 2, {}}}}};
    }
    if (nonzero.value() == Tribool::True) {
        return SolutionSet{FiniteSolutions{{
            FiniteSolution{positive, 1, {}},
            FiniteSolution{negative, 1, {}}}}};
    }

    auto zero = SymbolicExpr::number(0);
    auto distinct =
        polynomial_relation(discriminant, zero, RelationOp::NEQ);
    auto coincident =
        polynomial_relation(discriminant, zero, RelationOp::EQ);
    return SolutionSet{FiniteSolutions{{
        FiniteSolution{positive, 1, {distinct}},
        FiniteSolution{negative, 1, {distinct}},
        FiniteSolution{positive, 2, {coincident}}}}};
}

static SolveResult solve_symbolic_polynomial_set(
    const Polynomial<SymbolicPolyCoeff>& polynomial,
    const ExprPtr& original, const std::string& variable,
    ComputationContext& context, const SolveOptions& options) {
    auto valid = validate_symbolic_polynomial_input(original, options);
    if (!valid) {
        return SolveResult::failure(valid.error());
    }

    std::optional<detail::AssumptionFacts> adapter;
    if (context.assumptions()) {
        adapter.emplace(*context.assumptions());
    }
    const FactsQuery& facts = adapter
        ? static_cast<const FactsQuery&>(*adapter)
        : detail::no_facts();
    int degree = polynomial.degree();
    while (degree >= 0) {
        auto leading = get_coeff(polynomial, degree)->simplify();
        auto nonzero = detail::query_nonzero_value(
            detail::node(leading), facts, Domain::Complex, context);
        if (!nonzero) {
            return SolveResult::failure(nonzero.error());
        }
        if (nonzero.value() == Tribool::False) {
            --degree;
            continue;
        }
        if (nonzero.value() == Tribool::Unknown || degree >= 3) {
            return degree == 1
                ? conditional_linear_polynomial(
                    polynomial, leading, variable)
                : conditional_equation(original, variable);
        }
        if (degree == 0) {
            return SolutionSet{EmptySolutions{}};
        }
        if (degree == 2) {
            return solve_symbolic_quadratic(
                polynomial, leading, facts, context);
        }
        auto value = SymbolicExpr::divide(
            SymbolicExpr::multiply(
                SymbolicExpr::number(-1),
                get_coeff(polynomial, 0)),
            leading)->simplify();
        return SolutionSet{FiniteSolutions{{
            FiniteSolution{value, 1, {}}}}};
    }
    return SolutionSet{UniversalSolutions{}};
}

static SolveResult solve_product_factor(const detail::SymbolicNodePtr& node,
    const std::string& variable, ComputationContext& context, const SolveOptions& options,
    std::size_t& multiplicity) {
    auto factor = detail::make_expression_ptr(node);
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        BigInt exponent;
        if (try_get_integer_value(
                std::dynamic_pointer_cast<const NumberNode>(power->exponent()), exponent)) {
            if (exponent <= BigInt(0)) { return SolutionSet{EmptySolutions{}}; }
            if (exponent > BigInt(static_cast<std::int64_t>(std::numeric_limits<int>::max()))) {
                return SolveResult::failure(CasErrc::ResourceLimit, "root multiplicity overflow", "solve.product");
            }
            multiplicity = static_cast<std::size_t>(*exponent.try_to_int64());
            factor = detail::make_expression_ptr(power->base());
        }
    }
    return solve_equation(factor, variable, context, options);
}

static Result<void> append_product_roots(FiniteSolutions& values, std::size_t multiplicity,
    std::vector<FiniteSolution>& finite) {
    for (auto& value : values.values) {
        if (value.multiplicity > std::numeric_limits<std::size_t>::max() / multiplicity) {
            return Result<void>::failure(CasErrc::ResourceLimit, "root multiplicity overflow", "solve.product");
        }
        value.multiplicity *= multiplicity;
        auto existing = std::find_if(finite.begin(), finite.end(), [&](const FiniteSolution& candidate) {
            return detail::node(candidate.value)->equals(*detail::node(value.value)) &&
                candidate.conditions.empty() && value.conditions.empty();
        });
        if (existing == finite.end()) { finite.push_back(std::move(value)); }
        else {
            if (existing->multiplicity > std::numeric_limits<std::size_t>::max() - value.multiplicity) {
                return Result<void>::failure(CasErrc::ResourceLimit, "root multiplicity overflow", "solve.product");
            }
            existing->multiplicity += value.multiplicity;
        }
    }
    return Result<void>::success();
}

static void append_product_families(ParametricSolutions& values,
    const std::set<std::string>& original_names, std::vector<ParametricSolution>& families) {
    for (auto& family : values.values) {
        auto occupied = original_names;
        const auto value_names = detail::all_variable_names(detail::node(family.value));
        occupied.insert(value_names.begin(), value_names.end());
        for (auto& parameter : family.integer_parameters) {
            if (!original_names.count(parameter)) { continue; }
            std::string fresh = "_k";
            for (std::size_t suffix = 1; occupied.count(fresh); ++suffix) {
                fresh = "_k" + std::to_string(suffix);
            }
            auto symbol = SymbolicExpr::variable(fresh);
            family.value = substitute_raw(family.value, parameter, symbol);
            for (auto& condition : family.conditions) {
                condition = substitute_raw(condition, parameter, symbol);
            }
            parameter = fresh;
            occupied.insert(fresh);
        }
        families.push_back(std::move(family));
    }
}

static std::optional<SolveResult> solve_product_set(const ExprPtr& expression,
    const std::string& variable, ComputationContext& context, const SolveOptions& options) {
    auto product = std::dynamic_pointer_cast<const MultiplyNode>(detail::node(expression));
    if (!product) { return std::nullopt; }
    std::vector<FiniteSolution> finite;
    std::vector<ParametricSolution> families;
    const auto original_names = detail::all_variable_names(detail::node(expression));
    for (const auto& node : product->operands()) {
        std::size_t multiplicity = 1;
        auto solved = solve_product_factor(node, variable, context, options, multiplicity);
        if (!solved) {
            if (solved.error().code == CasErrc::Inconclusive) { return conditional_equation(expression, variable); }
            return SolveResult::failure(solved.error());
        }
        if (std::holds_alternative<EmptySolutions>(solved.value())) { continue; }
        if (std::holds_alternative<UniversalSolutions>(solved.value())) { return conditional_equation(expression, variable); }
        if (auto values = std::get_if<FiniteSolutions>(&solved.value())) {
            auto appended = append_product_roots(*values, multiplicity, finite);
            if (!appended) { return SolveResult::failure(appended.error()); }
        } else if (auto values = std::get_if<ParametricSolutions>(&solved.value())) {
            append_product_families(*values, original_names, families);
        } else { return conditional_equation(expression, variable); }
    }
    if (!finite.empty() && !families.empty()) { return conditional_equation(expression, variable); }
    if (!families.empty()) { return SolveResult::success(SolutionSet{ParametricSolutions{std::move(families)}}); }
    if (!finite.empty()) { return SolveResult::success(SolutionSet{FiniteSolutions{std::move(finite)}}); }
    return SolveResult::success(SolutionSet{EmptySolutions{}});
}

static SolveResult collect_exact_solutions(
    std::vector<std::shared_ptr<SymbolicExpr>> roots, ComputationContext& context) {
    constexpr const char* operation = "solve_equation";
    std::vector<FiniteSolution> solutions;
    for (auto& root : roots) {
        auto root_step = context.consume_steps(1, operation);
        if (!root_step) {
            return SolveResult::failure(root_step.error());
        }
        if (!root || !LMCAS::detail::node(root)) {
            return SolveResult::failure(
                CasErrc::InternalInvariant,
                "polynomial solver produced a null root", operation);
        }
        auto existing = std::find_if(
            solutions.begin(), solutions.end(),
            [&](const FiniteSolution& solution) {
                return LMCAS::detail::node(solution.value)->equals(
                    *LMCAS::detail::node(root));
            });
        if (existing != solutions.end()) {
            ++existing->multiplicity;
        } else {
            solutions.push_back(
                FiniteSolution{std::move(root), 1, {}});
        }
    }
    return SolveResult::success(
        SolutionSet{FiniteSolutions{std::move(solutions)}});
}

static SolveVectorResult certified_quartic_roots(const Polynomial<Rational>& remainder,
    const ExprPtr& expression, const std::string& var, ComputationContext& context) {
    auto candidates = solve_quartic(
        SymbolicExpr::number(remainder.coeffs[4]),
        SymbolicExpr::number(remainder.coeffs[3]),
        SymbolicExpr::number(remainder.coeffs[2]),
        SymbolicExpr::number(remainder.coeffs[1]),
        SymbolicExpr::number(remainder.coeffs[0]), var);
    bool certified = candidates.size() == 4;
    for (std::size_t i = 0; certified && i < candidates.size(); ++i) {
        auto proof = check_zero_residual(substitute_raw(expression, var, candidates[i]), context);
        if (!proof) { return SolveVectorResult::failure(proof.error()); }
        certified = std::holds_alternative<ProvedZeroResidual>(proof.value());
        for (std::size_t j = 0; certified && j < i; ++j) {
            auto equal = detail::compare_pointwise_values(
                candidates[i], candidates[j], context, Domain::Complex);
            if (!equal) { return SolveVectorResult::failure(equal.error()); }
            certified = equal.value() == Tribool::False;
        }
    }
    if (certified) { return candidates; }
    return std::vector<ExprPtr>{};
}

static SolveVectorResult solve_higher_exact_roots(const Polynomial<Rational>& remainder,
    bool bounded_candidates, const std::string& var, ComputationContext& context) {
    auto expression = poly_to_symbolic(remainder);
    std::vector<ExprPtr> represented_roots;
    if (remainder.degree() == 4 && bounded_candidates) {
        auto candidates = certified_quartic_roots(remainder, expression, var, context);
        if (!candidates) { return SolveVectorResult::failure(candidates.error()); }
        represented_roots = std::move(candidates.value());
    }
    if (represented_roots.empty()) {
        for (int index = 0; index < remainder.degree(); ++index) {
            auto root = make_rootof_checked(expression, var, static_cast<std::size_t>(index), context);
            if (!root) { return SolveVectorResult::failure(root.error()); }
            represented_roots.push_back(std::move(root.value()));
        }
    }
    return represented_roots;
}

static SolveVectorResult solve_exact_factor(Polynomial<Rational> remainder,
    const std::string& var, ComputationContext& context) {
    std::vector<ExprPtr> factor_roots;
    const bool bounded_candidates =
        polynomial_solver_detail::bounded_rational_root_search(remainder);
    if (bounded_candidates) {
        for (const auto& root : find_rational_roots(remainder)) {
            Polynomial<Rational> divisor({-root, Rational(1)}, var);
            auto division = remainder.div_mod(divisor);
            if (!division.second.is_zero()) { continue; }
            factor_roots.push_back(SymbolicExpr::number(root));
            remainder = std::move(division.first);
        }
    }
    if (remainder.degree() == 1) {
        factor_roots.push_back(SymbolicExpr::number(-remainder.coeffs[0] / remainder.coeffs[1]));
    } else if (remainder.degree() == 2) {
        auto quadratic = detail::real_quadratic_roots_checked(
            SymbolicExpr::number(remainder.coeffs[2]),
            SymbolicExpr::number(remainder.coeffs[1]),
            SymbolicExpr::number(remainder.coeffs[0]), context);
        if (!quadratic) { return SolveVectorResult::failure(quadratic.error()); }
        for (const auto& root : quadratic.value()) {
            factor_roots.push_back(detail::make_expression_ptr(
                SymbolicFactory::create_complex(detail::node(root.real), detail::node(root.imag))));
        }
    } else if (remainder.degree() > 2) {
        auto represented_roots = solve_higher_exact_roots(remainder, bounded_candidates, var, context);
        if (!represented_roots) { return SolveVectorResult::failure(represented_roots.error()); }
        factor_roots.insert(factor_roots.end(), represented_roots.value().begin(), represented_roots.value().end());
    }
    return factor_roots;
}

static SolveResult solve_exact_polynomial_set(const Polynomial<Rational>& exact,
                                            const std::string& var,
                                            ComputationContext& context) {
    if (exact.is_zero()) { return SolutionSet{UniversalSolutions{}}; }
    if (exact.degree() == 0) { return SolutionSet{EmptySolutions{}}; }
    std::vector<ExprPtr> roots;
    for (auto& factor : square_free_factorization(exact)) {
        auto step = context.consume_steps(1, "solve.exact");
        if (!step) { return SolveResult::failure(step.error()); }
        auto solved = solve_exact_factor(factor.first, var, context);
        if (!solved) { return SolveResult::failure(solved.error()); }
        const auto& factor_roots = solved.value();
        if (factor_roots.size() != static_cast<std::size_t>(factor.first.degree()))
            { return SolveResult::failure(CasErrc::InternalInvariant,
                "square-free factor root count mismatch", "solve.exact"); }
        for (int copy = 0; copy < factor.second; ++copy)
            { roots.insert(roots.end(), factor_roots.begin(), factor_roots.end()); }
    }
    if (roots.size() != static_cast<std::size_t>(exact.degree()))
        { return SolveResult::failure(CasErrc::InternalInvariant,
            "polynomial root multiplicities do not cover its degree", "solve.exact"); }
    return collect_exact_solutions(std::move(roots), context);
}

static SolveResult dispatch_equation(const ExprPtr& expr, const ExprPtr& f_expr,
    const std::string& var, ComputationContext& context, const SolveOptions& opts) {
    auto polynomial = recognize_rational_polynomial(*f_expr, var, context);
    if (!polynomial) { return SolveResult::failure(polynomial.error()); }
    if (polynomial.value()) {
        auto solved = solve_exact_polynomial_set(*polynomial.value(), var, context);
        if (!solved) { return solved; }
        return detail::finalize_solution_set(expr, var, std::move(solved.value()), context, Domain::Complex, opts);
    }
    auto symbolic = symbolic_to_poly<SymbolicPolyCoeff>(f_expr, var);
    if (symbolic) {
        auto solved = solve_symbolic_polynomial_set(symbolic.value(), f_expr, var, context, opts);
        if (!solved) { return solved; }
        return detail::finalize_solution_set(expr, var, std::move(solved.value()), context, Domain::Complex, opts);
    }
    if (symbolic.error().code != CasErrc::UnsupportedExpression) {
        return SolveResult::failure(symbolic.error());
    }
    if (auto product = solve_product_set(f_expr, var, context, opts)) {
        if (!*product) { return std::move(*product); }
        return detail::finalize_solution_set(expr, var, std::move(product->value()), context, Domain::Real, opts);
    }
    auto solved = solve_set_core(f_expr, var, context, opts);
    if (!solved) { return solved; }
    return detail::finalize_solution_set(expr, var, std::move(solved.value()), context, Domain::Real, opts);
}

LMCAS::SolveResult solve_equation(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts) {
    constexpr const char* operation = "solve_equation";
    if (!expr || !LMCAS::detail::node(expr)) {
        return SolveResult::failure(
            CasErrc::InvalidArgument, "equation expression cannot be null", operation);
    }
    if (var.empty()) {
        return SolveResult::failure(
            CasErrc::InvalidArgument, "solve variable cannot be empty", operation);
    }

    auto step = context.consume_steps(1, operation);
    if (!step) {
        return SolveResult::failure(step.error());
    }

    std::shared_ptr<SymbolicExpr> f_expr = expr;
    if (auto relation = std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(expr))) {
        if (relation->op() != RelationalNode::Op::EQ) {
            return SolveResult::failure(
                CasErrc::InvalidArgument,
                "solve_equation accepts equations, not inequalities",
                operation);
        }
        auto left = LMCAS::detail::make_expression_ptr(relation->left());
        auto right = LMCAS::detail::make_expression_ptr(relation->right());
        f_expr = SymbolicExpr::add(
            left, SymbolicExpr::multiply(SymbolicExpr::number(-1), right));
    }

    return dispatch_equation(expr, f_expr, var, context, opts);
}

LMCAS::SolveResult solve_equation(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts) {
    ComputationContext context;
    return solve_equation(expr, var, context, opts);
}


static std::optional<SolveVectorResult> solve_finite_mixed(
    const std::shared_ptr<SymbolicExpr>& f_expr, const std::string& var,
    ComputationContext& context, const SolveOptions& opts) {
    constexpr const char* operation = "solve_dispatch_vector";
    if (!contains_transcendental_of_var(f_expr, var)) {
        return std::nullopt;
    }
    auto substitution = detect_trans_substitutions(f_expr, var);
    if (!substitution.mappings.empty() && is_polynomial_after_substitution(substitution)) {
        return std::nullopt;
    }
    if (!opts.allow_numeric) {
        return SolveVectorResult::failure(
            CasErrc::Inconclusive,
            "mixed transcendental equation is outside the finite-vector support domain",
            operation);
    }
    auto mixed = solve_mixed_transcendental_checked(f_expr, var, opts, context);
    if (!mixed) {
        return SolveVectorResult::failure(mixed.error());
    }
    if (mixed.value().completeness == Completeness::Complete) {
        return SolveVectorResult::success(std::move(mixed.value().value));
    }
    return SolveVectorResult::failure(
        CasErrc::Inconclusive,
        mixed.value().reason.empty() ? "mixed transcendental search is incomplete" : mixed.value().reason,
        operation);
}

static SolveResult finite_values_to_set(SolveVectorResult result) {
    if (!result) return SolveResult::failure(result.error());
    if (result.value().empty()) {
        return SolveResult::success(SolutionSet{EmptySolutions{}});
    }
    std::vector<FiniteSolution> solutions;
    solutions.reserve(result.value().size());
    for (auto& value : result.value()) {
        solutions.push_back(FiniteSolution{std::move(value), 1, {}});
    }
    return SolveResult::success(SolutionSet{FiniteSolutions{std::move(solutions)}});
}

static SolveResult solve_set_core(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts) {
    constexpr const char* operation = "solve_dispatch_set";
    if (!expr || !LMCAS::detail::node(expr)) {
        return SolveResult::failure(CasErrc::InvalidArgument,
                                          "expression cannot be null",
                                          operation);
    }
    if (var.empty()) {
        return SolveResult::failure(CasErrc::InvalidArgument,
                                          "solve variable cannot be empty",
                                          operation);
    }

    auto simplified = expr->simplify();
    if (!simplified || !LMCAS::detail::node(simplified)) {
        return SolveResult::failure(CasErrc::InvalidArgument,
                                          "expression simplification failed",
                                          operation);
    }

    std::shared_ptr<SymbolicExpr> f_expr = simplified;
    if (auto rel = std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(simplified))) {
        if (rel->op() == RelationalNode::Op::EQ) {
            auto left = LMCAS::detail::make_expression_ptr(rel->left());
            auto right = LMCAS::detail::make_expression_ptr(rel->right());
            f_expr = SymbolicExpr::add(
                left, SymbolicExpr::multiply(right, SymbolicExpr::number(-1)));
        }
    }

    auto polynomial = symbolic_to_poly<SymbolicPolyCoeff>(f_expr, var);
    if (polynomial) {
        return solve_symbolic_polynomial_set(polynomial.value(), expr, var, context, opts);
    }
    if (polynomial.error().code != CasErrc::UnsupportedExpression) {
        return SolveResult::failure(polynomial.error());
    }

    auto transcendental = solve_transcendental(f_expr, var, context, opts);
    if (transcendental) return transcendental;
    if (transcendental.error().code != CasErrc::Inconclusive) {
        return SolveResult::failure(transcendental.error());
    }

    if (auto mixed = solve_finite_mixed(f_expr, var, context, opts)) {
        return finite_values_to_set(std::move(*mixed));
    }


    return SolveResult::failure(CasErrc::Inconclusive,
        "expression is outside the supported finite exact solve domain", operation);
}

LMCAS::FiniteSolveResult solve_finite_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts) {
    constexpr const char* operation = "solve_finite_projection";
    auto solved = solve_equation(expr, var, context, opts);
    if (!solved) {
        return SolveVectorResult::failure(solved.error());
    }
    auto solution_set = std::move(solved.value());
    if (std::holds_alternative<EmptySolutions>(solution_set)) {
        return SolveVectorResult::success({});
    }
    auto* finite = std::get_if<FiniteSolutions>(&solution_set);
    if (!finite) {
        return SolveVectorResult::failure(
            CasErrc::Inconclusive,
            "solution set is not finitely enumerable", operation);
    }
    std::size_t value_count = 0;
    for (const auto& solution : finite->values) {
        if (!solution.conditions.empty()) {
            return SolveVectorResult::failure(CasErrc::Inconclusive,
                "Conditional roots cannot be projected to an unconditional vector", operation);
        }
        value_count += solution.multiplicity;
    }
    std::vector<std::shared_ptr<SymbolicExpr>> values;
    values.reserve(value_count);
    for (const auto& solution : finite->values) {
        for (std::size_t copy = 0; copy < solution.multiplicity; ++copy) {
            values.push_back(solution.value);
        }
    }
    return SolveVectorResult::success(std::move(values));
}

LMCAS::FiniteSolveResult solve_finite_checked(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable,
    const SolveOptions& options) {
    ComputationContext context;
    return solve_finite_checked(
        expression, variable, context, options);
}



std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
SymbolicExpr::solve_system(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& vars) {
    std::vector<SymbolicExpr> values;
    values.reserve(equations.size());
    for (const auto& equation : equations) {
        if (!equation) {
            return {};
        }
        if (auto relation =
                std::dynamic_pointer_cast<const RelationalNode>(
                    LMCAS::detail::node(equation))) {
            if (relation->op() != RelationalNode::Op::EQ) {
                return {};
            }
            auto normalized = SymbolicExpr::add(
                LMCAS::detail::make_expression_ptr(relation->left()),
                SymbolicExpr::multiply(
                    SymbolicExpr::number(-1),
                    LMCAS::detail::make_expression_ptr(relation->right())));
            values.push_back(*normalized->simplify());
        } else {
            values.push_back(*equation);
        }
    }
    auto solved = LMCAS::Solver::solve_linear_system(values, vars);
    if (solved.empty()) {
        return {};
    }
    std::map<std::string, std::shared_ptr<SymbolicExpr>> projected;
    for (auto& [variable, value] : solved) {
        projected.emplace(
            variable, std::make_shared<SymbolicExpr>(std::move(value)));
    }
    return {std::move(projected)};
}

std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> SymbolicExpr::solve_system(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>& parameters) {
    return ParametricSolver::solve_system(equations, unknowns, parameters);
}


std::shared_ptr<SymbolicExpr> SymbolicExpr::transpose(const std::shared_ptr<SymbolicExpr>& mat) {
    if (!mat) {
        return mat;
    }
    auto m_node = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(mat));
    if (!m_node) {
        return mat;
    }
    size_t r = m_node->rows();
    size_t c = m_node->cols();

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> new_data(c, std::vector<std::shared_ptr<SymbolicExpr>>(r));

    for(size_t i=0; i<r; ++i) {
        for(size_t j=0; j<c; ++j) {
            new_data[j][i] = LMCAS::detail::make_expression_ptr(m_node->get(i,j));
        }
    }

    return SymbolicExpr::matrix(new_data);
}

}
