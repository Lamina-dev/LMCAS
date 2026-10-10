#include "solver.hpp"
#include "solve_strategies.hpp"
#include "fglm.hpp"
#include "internal/solver_groebner_builder.hpp"
#include "poly_utils.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/solver_support.hpp"
#include "internal/assumption_facts.hpp"
#include "residual_verification.hpp"
#include <iterator>
#include <limits>
#include <optional>

namespace LMCAS {
namespace {
using namespace solver_detail;
using PolynomialSolution = std::map<std::string, SymbolicExpr>;
using PolynomialSolutions = std::vector<PolynomialSolution>;
using Expressions = std::vector<std::shared_ptr<SymbolicExpr>>;

std::shared_ptr<SymbolicExpr> clear_denominators(
    const SymbolicExpr& equation, Expressions& constraints) {
    if (!detail::node(equation)) return nullptr;
    std::vector<std::shared_ptr<const SymbolicNode>> factors;
    Expressions local_constraints;
    if (!collect_denominator_factors(detail::node(equation), factors, local_constraints)) {
        return nullptr;
    }
    auto cleared = to_ptr(equation);
    if (factors.empty()) {
        cleared = cleared->simplify();
    } else if (auto sum = std::dynamic_pointer_cast<const AddNode>(detail::node(cleared))) {
        std::vector<std::shared_ptr<const SymbolicNode>> terms;
        terms.reserve(sum->operands().size());
        for (const auto& term : sum->operands()) {
            auto product = multiply_no_expand(term, factors);
            terms.push_back(detail::node(product));
        }
        cleared = detail::make_expression_ptr(detail::make_node<AddNode>(terms));
    } else {
        cleared = multiply_no_expand(detail::node(cleared), factors);
    }
    if (!cleared || !detail::node(cleared) || !is_polynomial_node(detail::node(cleared))) {
        return nullptr;
    }
    for (const auto& constraint : local_constraints) constraints.push_back(constraint);
    return cleared;
}

std::shared_ptr<SymbolicExpr> substitute_solution(
    const std::shared_ptr<SymbolicExpr>& expression, const PolynomialSolution& solution) {
    auto result = expression;
    for (const auto& [name, value] : solution) {
        result = result->substitute(name, detail::make_expression_ptr(value));
        if (!result) return nullptr;
    }
    return result->simplify();
}

PolynomialSystemResult filter_denominators(PolynomialSolutions candidates,
    const Expressions& constraints, ComputationContext& context) {
    std::optional<detail::AssumptionFacts> adapter;
    if (context.assumptions()) adapter.emplace(*context.assumptions());
    const FactsQuery& facts = adapter ? static_cast<const FactsQuery&>(*adapter) : detail::no_facts();
    PolynomialSolutions filtered;
    filtered.reserve(candidates.size());
    for (auto& solution : candidates) {
        bool valid = true;
        for (const auto& denominator : constraints) {
            auto bound = substitute_solution(denominator, solution);
            if (!bound) return PolynomialSystemResult::failure(CasErrc::InternalInvariant,
                "denominator substitution failed", "solve_polynomial_system");
            auto nonzero = detail::query_nonzero_value(detail::node(bound), facts, Domain::Complex, context);
            if (!nonzero) return PolynomialSystemResult::failure(nonzero.error());
            if (nonzero.value() == Tribool::False) { valid = false; break; }
            if (nonzero.value() != Tribool::True)
                return PolynomialSystemResult::failure(CasErrc::Inconclusive,
                    "candidate denominator is not proved nonzero", "solve_polynomial_system");
        }
        if (valid) filtered.push_back(std::move(solution));
    }
    return filtered;
}

PolynomialSystemResult solve_single_equation(const SymbolicExpr& equation,
                                            const std::string& variable,
                                            ComputationContext& context) {
    auto solved = solve_finite_checked(detail::make_expression_ptr(equation),
                                       variable, context, SolveOptions{});
    if (!solved) return PolynomialSystemResult::failure(solved.error());
    PolynomialSolutions solutions;
    for (const auto& root : solved.value()) solutions.push_back({{variable, *root}});
    return solutions;
}

Result<void> validate_basis_input(const std::vector<SymbolicExpr>& equations,
                                  const std::vector<std::string>& variables,
                                  groebner_detail::PolyContext& conversion,
                                  ComputationContext& context) {
    constexpr const char* operation = "solve_polynomial_system";
    for (const auto& equation : equations) {
        auto budget = context.consume_steps(1, operation);
        if (!budget) return Result<void>::failure(budget.error());
        groebner_detail::PolyBuilder builder(variables, conversion, true);
        detail::node(equation)->accept(builder);
        if (builder.failed) return Result<void>::failure(
            CasErrc::Inconclusive, "equation is not an exact polynomial in the system variables", operation);
    }
    return Result<void>::success();
}

Result<groebner_detail::Poly> convert_basis_polynomial(
    const SymbolicExpr& polynomial, const std::vector<std::string>& variables,
    groebner_detail::PolyContext& conversion, ComputationContext& context,
    std::shared_ptr<SymbolicExpr>& simplified) {
    constexpr const char* operation = "solve_polynomial_system";
    auto budget = context.consume_steps(1, operation);
    if (!budget) return Result<groebner_detail::Poly>::failure(budget.error());
    simplified = detail::make_expression_ptr(polynomial)->simplify();
    if (!simplified) return Result<groebner_detail::Poly>::failure(
        CasErrc::Inconclusive, "basis polynomial cannot be converted", operation);
    groebner_detail::PolyBuilder builder(variables, conversion, true);
    detail::node(simplified)->accept(builder);
    if (builder.failed) return Result<groebner_detail::Poly>::failure(
        CasErrc::Inconclusive, "basis polynomial cannot be converted", operation);
    return builder.get_result();
}

Result<Expressions> simplified_basis(const std::vector<SymbolicExpr>& equations,
                                     const std::vector<std::string>& variables,
                                     ComputationContext& context) {
    constexpr const char* operation = "solve_polynomial_system";
    groebner_detail::PolyContext conversion(variables);
    auto validated = validate_basis_input(equations, variables, conversion, context);
    if (!validated) return Result<Expressions>::failure(validated.error());
    auto groebner = Solver::reduced_groebner_basis(equations, variables);
    Expressions basis;
    std::vector<FGLMPoly> leading;
    basis.reserve(groebner.size());
    leading.reserve(groebner.size());
    for (const auto& polynomial : groebner) {
        std::shared_ptr<SymbolicExpr> simplified;
        auto converted = convert_basis_polynomial(polynomial, variables, conversion, context, simplified);
        if (!converted) return Result<Expressions>::failure(converted.error());
        if (simplified->is_number() && !simplified->is_zero()) return Expressions{simplified};
        if (converted.value().is_zero()) continue;
        FGLMPoly term(variables.size());
        term.add_term(converted.value().lead_monomial(), converted.value().lead_coeff());
        leading.push_back(std::move(term));
        basis.push_back(std::move(simplified));
    }
    if (!is_zero_dimensional(leading, variables.size())) {
        return Result<Expressions>::failure(CasErrc::Inconclusive,
            "polynomial system has no certified finite solution set", operation);
    }
    return basis;
}
class PolynomialSystemSearch {
public:
    PolynomialSystemSearch(const Expressions& basis,
                           const std::vector<std::string>& variables,
                           ComputationContext& context)
        : basis_(basis), variables_(variables), context_(context) {}

    PolynomialSystemResult solve(int variable, const PolynomialSolution& partial) const;

private:
    const Expressions& basis_;
    const std::vector<std::string>& variables_;
    ComputationContext& context_;

    bool depends_on_remaining(const SymbolicExpr& expression, int last) const;
    Result<std::optional<Expressions>> substitute_basis(int variable,
                                                const PolynomialSolution& partial) const;
    Result<std::shared_ptr<SymbolicExpr>> select_equation(const Expressions& reduced,
                                                         int variable, bool& appears) const;
};

bool PolynomialSystemSearch::depends_on_remaining(const SymbolicExpr& expression,
                                                   int last) const {
    for (int index = 0; index <= last && index < static_cast<int>(variables_.size()); ++index) {
        if (contains(expression, variables_[index])) return true;
    }
    return false;
}

Result<std::optional<Expressions>> PolynomialSystemSearch::substitute_basis(
    int variable, const PolynomialSolution& partial) const {
    Expressions reduced;
    reduced.reserve(basis_.size());
    for (const auto& polynomial : basis_) {
        auto substituted = substitute_solution(polynomial, partial);
        if (!substituted) { return Result<std::optional<Expressions>>::failure(
            CasErrc::InternalInvariant, "basis substitution failed", "solve_polynomial_system"); }
        if (substituted->is_zero()) { continue; }
        if (!depends_on_remaining(*substituted, variable)) {
            auto proof = check_zero_residual(substituted, context_);
            if (!proof) { return Result<std::optional<Expressions>>::failure(proof.error()); }
            if (std::holds_alternative<ProvedZeroResidual>(proof.value())) { continue; }
            std::optional<detail::AssumptionFacts> adapter;
            if (context_.assumptions()) { adapter.emplace(*context_.assumptions()); }
            const FactsQuery& facts = adapter ? static_cast<const FactsQuery&>(*adapter) : detail::no_facts();
            auto nonzero = detail::query_nonzero_value(detail::node(substituted), facts, Domain::Complex, context_);
            if (!nonzero) { return Result<std::optional<Expressions>>::failure(nonzero.error()); }
            if (nonzero.value() == Tribool::True) { return std::optional<Expressions>{}; }
            return Result<std::optional<Expressions>>::failure(CasErrc::Inconclusive,
                "remaining system constraint is not proved", "solve_polynomial_system");
        }
        reduced.push_back(substituted);
    }
    return std::optional<Expressions>{std::move(reduced)};
}

Result<std::shared_ptr<SymbolicExpr>> PolynomialSystemSearch::select_equation(
    const Expressions& reduced, int variable, bool& appears) const {
    std::shared_ptr<SymbolicExpr> target;
    int best_degree = std::numeric_limits<int>::max();
    const auto& name = variables_[variable];
    for (const auto& polynomial : reduced) {
        if (!contains(*polynomial, name)) continue;
        appears = true;
        if (depends_on_remaining(*polynomial, variable - 1)) continue;
        auto univariate = symbolic_to_poly<SymbolicPolyCoeff>(polynomial, name);
        if (!univariate) {
            if (univariate.error().code != CasErrc::UnsupportedExpression) {
                return Result<std::shared_ptr<SymbolicExpr>>::failure(univariate.error());
            }
            continue;
        }
        int degree = univariate.value().degree();
        if (degree >= 1 && degree < best_degree) {
            best_degree = degree;
            target = polynomial;
        }
    }
    return target;
}

PolynomialSystemResult PolynomialSystemSearch::solve(
    int variable, const PolynomialSolution& partial) const {
    auto reduced = substitute_basis(variable, partial);
    if (!reduced) return PolynomialSystemResult::failure(reduced.error());
    if (!reduced.value()) return PolynomialSolutions{};
    if (variable < 0) return PolynomialSolutions{partial};
    const auto& name = variables_[variable];
    bool appears = false;
    auto target = select_equation(*reduced.value(), variable, appears);
    if (!target) return PolynomialSystemResult::failure(target.error());
    if (!appears) {
        return PolynomialSystemResult::failure(CasErrc::Inconclusive,
            "remaining variable is not constrained", "solve_polynomial_system");
    }
    if (!target.value()) {
        return PolynomialSystemResult::failure(CasErrc::Inconclusive,
            "No convertible polynomial target for remaining variable", "solve_polynomial_system");
    }
    auto solved = solve_finite_checked(target.value(), name, context_, SolveOptions{});
    if (!solved) return PolynomialSystemResult::failure(solved.error());
    PolynomialSolutions results;
    for (const auto& root : solved.value()) {
        auto next = partial;
        next.insert_or_assign(name, *root);
        auto solved_branch = solve(variable - 1, next);
        if (!solved_branch) return solved_branch;
        auto& values = solved_branch.value();
        results.insert(results.end(), std::make_move_iterator(values.begin()),
                       std::make_move_iterator(values.end()));
    }
    return results;
}

PolynomialSystemResult solve_polynomial_system_impl(
    const std::vector<SymbolicExpr>& equations,
    const std::vector<std::string>& variables, ComputationContext& context) {
    std::vector<SymbolicExpr> cleared_equations;
    Expressions constraints;
    cleared_equations.reserve(equations.size());
    for (const auto& equation : equations) {
        auto cleared = clear_denominators(equation, constraints);
        if (!cleared) return PolynomialSystemResult::failure(CasErrc::Inconclusive,
            "equation cannot be converted to a polynomial without losing its domain", "solve_polynomial_system");
        cleared_equations.push_back(*cleared);
    }
    if (cleared_equations.size() == 1 && variables.size() == 1) {
        auto solved = solve_single_equation(cleared_equations[0], variables[0], context);
        if (!solved) return solved;
        return filter_denominators(std::move(solved.value()), constraints, context);
    }
    auto basis = simplified_basis(cleared_equations, variables, context);
    if (!basis) return PolynomialSystemResult::failure(basis.error());
    if (basis.value().size() == 1 && basis.value().front()->is_number() &&
        !basis.value().front()->is_zero()) return PolynomialSolutions{};
    PolynomialSystemSearch search(basis.value(), variables, context);
    auto solved = search.solve(static_cast<int>(variables.size()) - 1, PolynomialSolution{});
    if (!solved) return solved;
    return filter_denominators(std::move(solved.value()), constraints, context);
}
}

PolynomialSystemResult Solver::solve_polynomial_system_checked(
    const std::vector<SymbolicExpr>& equations,
    const std::vector<std::string>& variables, ComputationContext& context) {
    constexpr const char* operation = "solve_polynomial_system";
    if (equations.empty() || variables.empty()) {
        return PolynomialSystemResult::failure(
            CasErrc::InvalidArgument,
            "polynomial system requires equations and variables", operation);
    }
    auto budget = context.consume_steps(equations.size() * variables.size() + 1, operation);
    if (!budget) return PolynomialSystemResult::failure(budget.error());
    try {
        return solve_polynomial_system_impl(equations, variables, context);
    } catch (const std::bad_alloc&) {
        return PolynomialSystemResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while solving polynomial system", operation);
    } catch (const std::exception& error) {
        return PolynomialSystemResult::failure(CasErrc::InternalInvariant, error.what(), operation);
    }
}

PolynomialSystemResult Solver::solve_polynomial_system_checked(
    const std::vector<SymbolicExpr>& equations,
    const std::vector<std::string>& variables) {
    ComputationContext context;
    return solve_polynomial_system_checked(equations, variables, context);
}
}
