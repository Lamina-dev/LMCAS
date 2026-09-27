#include "expr.hpp"

#include <utility>
#include <vector>

#include "internal/expr_internal.hpp"

namespace LMCAS {

SolveResult solve_set(const ExprPtr& equation,
                      const std::string& variable,
                      ComputationContext& context,
                      const SolveOptions& options) {
    if (!equation) {
        return SolveResult::failure(CasErrc::InvalidArgument,
                                    "equation cannot be null", "LMCAS.solve_set");
    }
    if (variable.empty()) {
        return SolveResult::failure(CasErrc::InvalidArgument,
                                    "solve variable cannot be empty",
                                    "LMCAS.solve_set");
    }
    return solve_equation(equation, variable, context, options);
}

SolveResult solve_set(const ExprPtr& equation,
                      const std::string& variable,
                      const SolveOptions& options) {
    ComputationContext context;
    return solve_set(equation, variable, context, options);
}

static ExprSetResult project_finite_expr_set(const ExprPtr& original,
                                        const std::string& variable,
                                        ComputationContext& context,
                                        const SolveOptions& options) {
    auto solved = solve_set(original, variable, context, options);
    if (!solved) {
        return ExprSetResult::failure(solved.error());
    }

    const auto& solution_set = solved.value();
    if (std::holds_alternative<EmptySolutions>(solution_set)) {
        return expr_set({});
    }
    const auto* finite = std::get_if<FiniteSolutions>(&solution_set);
    if (!finite) {
        return expr_set_failure(
            CasErrc::Inconclusive,
            "solution set is not a finite enumerable set<Expr>",
            kSolveExprSetOperation);
    }

    std::vector<ExprPtr> elements;
    elements.reserve(finite->values.size());
    for (const auto& solution : finite->values) {
        if (!solution.conditions.empty()) {
            return expr_set_failure(CasErrc::Inconclusive,
                                    "conditional finite solutions cannot be lowered to set<Expr>",
                                    kSolveExprSetOperation);
        }
        elements.push_back(solution.value);
    }
    return expr_set(std::move(elements));
}

ExprSetResult solve_expr_set(const ExprPtr& equation,
                             const std::string& variable,
                             ComputationContext& context,
                             const SolveOptions& options) {
    return project_finite_expr_set(equation, variable, context, options);
}
ExprSetResult solve_expr_set(const ExprPtr& equation,
                             const std::string& variable,
                             const SolveOptions& options) {
    ComputationContext context;
    return solve_expr_set(equation, variable, context, options);
}

ExprSetResult roots(const ExprPtr& expression,
                    const std::string& variable,
                    ComputationContext& context,
                    const SolveOptions& options) {
    return project_finite_expr_set(expression, variable, context, options);
}

ExprSetResult roots(const ExprPtr& expression,
                    const std::string& variable,
                    const SolveOptions& options) {
    ComputationContext context;
    return roots(expression, variable, context, options);
}

ExprSetResult solve(const ExprPtr& equation,
                    const std::string& variable,
                    ComputationContext& context,
                    const SolveOptions& options) {
    return solve_expr_set(equation, variable, context, options);
}

ExprSetResult solve(const ExprPtr& equation,
                    const std::string& variable,
                    const SolveOptions& options) {
    ComputationContext context;
    return solve(equation, variable, context, options);
}

}
