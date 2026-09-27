#include "parametric_solver.hpp"
#include "solver.hpp"
#include "solve_strategies.hpp"
#include "poly_utils.hpp"
#include "internal/symbolic_ast.hpp"
#include "residual_verification.hpp"
#include "internal/assumption_facts.hpp"
#include <algorithm>
#include <limits>

namespace LMCAS {
namespace {

using ParametricAssignment = std::map<std::string, std::shared_ptr<SymbolicExpr>>;
using ParametricBasis = std::vector<std::shared_ptr<SymbolicExpr>>;

Result<Tribool> parameter_constraint(const ExprPtr& expression, ComputationContext& context) {
    std::optional<detail::AssumptionFacts> assumptions;
    if (context.assumptions()) { assumptions.emplace(*context.assumptions()); }
    const FactsQuery& facts = assumptions ? static_cast<const FactsQuery&>(*assumptions) : detail::no_facts();
    auto nonzero = detail::query_nonzero_value(detail::node(expression), facts, Domain::Complex, context);
    if (!nonzero) { return Result<Tribool>::failure(nonzero.error()); }
    if (nonzero.value() == Tribool::True) { return Tribool::False; }
    if (nonzero.value() == Tribool::False) { return Tribool::True; }
    auto residual = check_zero_residual(expression, context);
    if (!residual) { return Result<Tribool>::failure(residual.error()); }
    return std::holds_alternative<ProvedZeroResidual>(residual.value()) ?
        Tribool::True : Tribool::Unknown;
}

ParametricSolutionList free_unknowns(const std::vector<std::string>& unknowns) {
    ParametricAssignment solution;
    for (const auto& variable : unknowns) solution[variable] = SymbolicExpr::variable(variable);
    return {std::move(solution)};
}

std::vector<SymbolicExpr> parametric_input_polynomials(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations) {
    std::vector<SymbolicExpr> input_polys;
    input_polys.reserve(equations.size());
    for (const auto& eq : equations) {
        if (eq && LMCAS::detail::node(eq)) {
            auto simplified = eq->simplify();
            if (simplified && !simplified->is_zero()) {
                input_polys.push_back(*simplified);
            }
        }
    }
    return input_polys;
}

ParametricBasis parametric_basis(const std::vector<SymbolicExpr>& input_polys,
                                 const std::vector<std::string>& unknowns) {
    auto G_basis = Solver::groebner_basis(input_polys, unknowns);

    std::vector<std::shared_ptr<SymbolicExpr>> basis;
    basis.reserve(G_basis.size());
    for (const auto& g : G_basis) {
        auto g_ptr = LMCAS::detail::make_expression_ptr(g);
        auto simp = g_ptr->simplify();
        if (simp && !simp->is_zero()) {
            basis.push_back(simp);
        }
    }
    return basis;
}

bool depends_on_remaining_unknown(const SymbolicExpr& expression,
                                  const std::vector<std::string>& unknowns,
                                  int last) {
    for (int index = 0; index <= last && index < static_cast<int>(unknowns.size()); ++index) {
        if (contains(expression, unknowns[index])) {
            return true;
        }
    }
    return false;
}

std::shared_ptr<SymbolicExpr> substitute_assignment(
    const std::shared_ptr<SymbolicExpr>& expression, const ParametricAssignment& assignment) {
    auto result = expression;
    for (const auto& [name, value] : assignment) {
        result = result->substitute(name, value);
        if (!result) {
            return nullptr;
        }
    }
    return result->simplify();
}

class ParametricBranchSolver {
    const ParametricBasis& basis_;
    const std::vector<std::string>& unknowns_;
    ComputationContext& context_;

    Result<bool> reduce_basis(int position, const ParametricAssignment& partial,
                      ParametricBasis& reduced) const {
        reduced.reserve(basis_.size());
        for (const auto& polynomial : basis_) {
            auto expression = substitute_assignment(polynomial, partial);
            if (!expression) {
                return Result<bool>::failure(CasErrc::InternalInvariant,
                    "assignment substitution failed", "solve.parametric");
            }
            if (expression->is_zero()) {
                continue;
            }
            if (!depends_on_remaining_unknown(*expression, unknowns_, position)) {
                auto residual = parameter_constraint(expression, context_);
                if (!residual) return Result<bool>::failure(residual.error());
                if (residual.value() == Tribool::False) return false;
                if (residual.value() != Tribool::True)
                    return Result<bool>::failure(CasErrc::Inconclusive,
                        "remaining parameter constraint is not proved", "solve.parametric");
                continue;
            }
            reduced.push_back(expression);
        }
        return true;
    }

    Result<std::shared_ptr<SymbolicExpr>> select_target(const ParametricBasis& reduced,
                                                       int position) const {
        const auto& variable = unknowns_[position];
        std::shared_ptr<SymbolicExpr> target;
        int best_degree = std::numeric_limits<int>::max();
        bool appears = false;
        for (const auto& expression : reduced) {
            if (!contains(*expression, variable)) {
                continue;
            }
            appears = true;
            if (depends_on_remaining_unknown(*expression, unknowns_, position - 1)) {
                continue;
            }
            auto polynomial = symbolic_to_poly<SymbolicPolyCoeff>(expression, variable);
            if (!polynomial) {
                if (polynomial.error().code != CasErrc::UnsupportedExpression) {
                    return Result<std::shared_ptr<SymbolicExpr>>::failure(polynomial.error());
                }
                continue;
            }
            int degree = polynomial.value().degree();
            if (degree >= 1 && degree < best_degree) {
                best_degree = degree;
                target = expression;
            }
        }
        if (!target && appears) {
            return Result<std::shared_ptr<SymbolicExpr>>::failure(CasErrc::Inconclusive,
                "No convertible polynomial target for remaining variable", "solve.parametric");
        }
        return target;
    }

    ParametricSolutionsResult solve_free_unknown(int position,
                                                 const ParametricAssignment& partial) const {
        auto next = partial;
        const auto& variable = unknowns_[position];
        next[variable] = SymbolicExpr::variable(variable);
        return solve(position - 1, next);
    }

public:
    ParametricBranchSolver(const ParametricBasis& basis,
                           const std::vector<std::string>& unknowns,
                           ComputationContext& context)
        : basis_(basis), unknowns_(unknowns), context_(context) {}

    ParametricSolutionsResult solve(int position, const ParametricAssignment& partial) const {
        ParametricBasis reduced;
        auto consistent = reduce_basis(position, partial, reduced);
        if (!consistent) return ParametricSolutionsResult::failure(consistent.error());
        if (!consistent.value()) {
            return ParametricSolutionList{};
        }
        if (position < 0) {
            return ParametricSolutionList{partial};
        }
        auto selected = select_target(reduced, position);
        if (!selected) return ParametricSolutionsResult::failure(selected.error());
        if (!selected.value()) {
            return solve_free_unknown(position, partial);
        }
        const auto& variable = unknowns_[position];
        auto solved = solve_finite_checked(selected.value(), variable, context_);
        if (!solved) {
            return ParametricSolutionsResult::failure(solved.error());
        }
        auto roots = std::move(solved.value());
        if (roots.empty()) {
            return ParametricSolutionList{};
        }
        ParametricSolutionList results;
        for (const auto& root : roots) {
            auto next = partial;
            next[variable] = root;
            auto branch = solve(position - 1, next);
            if (!branch) {
                return branch;
            }
            auto& values = branch.value();
            results.insert(results.end(), std::make_move_iterator(values.begin()),
                           std::make_move_iterator(values.end()));
        }
        return results;
    }
};

}

ParametricSolutionsResult ParametricSolver::solve_polynomial_parametric_impl(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    ComputationContext& context)
{
    if (equations.empty() || unknowns.empty()) {
        return free_unknowns(unknowns);
    }
    auto input_polys = parametric_input_polynomials(equations);
    if (input_polys.empty()) {
        return free_unknowns(unknowns);
    }
    ParametricBasis basis;
    if (unknowns.size() == 1) {
        basis.reserve(input_polys.size());
        for (const auto& polynomial : input_polys)
            basis.push_back(detail::make_expression_ptr(polynomial));
    } else {
        basis = parametric_basis(input_polys, unknowns);
    }
    for (const auto& polynomial : basis) {
        if (!depends_on_remaining_unknown(*polynomial, unknowns, static_cast<int>(unknowns.size()) - 1)) {
            auto residual = parameter_constraint(polynomial, context);
            if (!residual) return ParametricSolutionsResult::failure(residual.error());
            if (residual.value() == Tribool::False)
                return ParametricSolutionList{};
            if (residual.value() != Tribool::True)
                return ParametricSolutionsResult::failure(CasErrc::Inconclusive,
                    "parameter-only constraint is not proved", "solve.parametric");
        }
    }
    ParametricBranchSolver branches(basis, unknowns, context);
    ParametricAssignment empty;
    return branches.solve(static_cast<int>(unknowns.size()) - 1, empty);
}

}
