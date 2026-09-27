#include "internal/vector_calculus_support.hpp"
#include "vector_calculus_extrema.hpp"
#include "integration.hpp"
#include "numeric_evaluation.hpp"
#include "solver.hpp"
#include "solve_strategies.hpp"
#include "symbolic_matrix.hpp"
#include "internal/symbolic_ast.hpp"
#include "residual_verification.hpp"

#include <cmath>
#include <exception>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace LMCAS {

using namespace vector_calculus_detail;

static VectorCalculusFieldResult lagrange_gradient_strict(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars, const std::string& operation)
{
    VectorField gradient;
    gradient.reserve(vars.size());
    for (const auto& var : vars) {
        auto partial = vector_calculus_differentiate_strict(f, var, operation);
        if (!partial) {
            return VectorCalculusFieldResult::failure(partial.error());
        }
        gradient.push_back(std::move(partial.value()));
    }
    return VectorCalculusFieldResult::success(std::move(gradient));
}


static std::vector<std::string> lagrange_multiplier_names(size_t count)
{
    std::vector<std::string> names;
    names.reserve(count);
    for (size_t k = 0; k < count; ++k) {
        names.push_back("lambda_" + std::to_string(k + 1));
    }
    return names;
}

static std::shared_ptr<SymbolicExpr> lagrange_stationarity_equation(
    const VectorField& gradient, const std::vector<VectorField>& constraints,
    const std::vector<std::string>& multipliers, size_t index)
{
    auto equation = gradient[index];
    for (size_t k = 0; k < multipliers.size(); ++k) {
        auto multiplier = SymbolicExpr::variable(multipliers[k]);
        auto term = SymbolicExpr::multiply(multiplier, constraints[k][index]);
        equation = SymbolicExpr::add(equation,
            SymbolicExpr::multiply(SymbolicExpr::number(-1), term));
    }
    return equation;
}

static LagrangeResult solve_lagrange_strict(
    const VectorField& equations, const std::vector<std::string>& all_vars,
    ComputationContext& context, const std::string& operation)
{
    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> full_solutions;
    std::vector<SymbolicExpr> poly_eqs;
    poly_eqs.reserve(equations.size());
    for (const auto& eq : equations) {
        if (!eq || !LMCAS::detail::node(eq)) {
            return LagrangeResult::failure(
                CasErrc::Inconclusive,
                "Lagrange equation construction failed in the supported domain",
                operation);
        }
        poly_eqs.push_back(*eq);
    }

    auto checked_solutions =
        Solver::solve_polynomial_system_checked(poly_eqs, all_vars, context);
    if (!checked_solutions) return LagrangeResult::failure(checked_solutions.error());
    const auto& poly_solutions = checked_solutions.value();
    for (const auto& sol : poly_solutions) {
        std::map<std::string, std::shared_ptr<SymbolicExpr>> pt;
        for (const auto& [name, val] : sol) {
            pt[name] = LMCAS::detail::make_expression_ptr(val);
        }
        full_solutions.push_back(std::move(pt));
    }

    return LagrangeResult::success(std::move(full_solutions));
}

static LagrangeResult verify_lagrange_solutions(
    const std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>& full_solutions,
    const VectorField& equations, const std::vector<std::string>& all_vars,
    const std::vector<std::string>& vars, ComputationContext& context, const std::string& operation)
{
    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> verified;
    std::set<std::string> var_set(vars.begin(), vars.end());
    for (const auto& full_solution : full_solutions) {
        if (!vector_calculus_point_has_vars(full_solution, all_vars)) {
            return LagrangeResult::failure(
                CasErrc::Inconclusive,
                "Lagrange candidate omits one or more variables or multipliers",
                operation);
        }
        for (const auto& eq : equations) {
            auto residual = eq;
            for (const auto& [variable, value] : full_solution)
                residual = residual->substitute(variable, value);
            auto proof = check_zero_residual(residual, context);
            if (!proof) return LagrangeResult::failure(proof.error());
            if (!std::holds_alternative<ProvedZeroResidual>(proof.value())) {
                return LagrangeResult::failure(
                    CasErrc::Inconclusive,
                    "Lagrange candidate does not verify against the full stationarity system",
                    operation);
            }
        }

        std::map<std::string, std::shared_ptr<SymbolicExpr>> filtered;
        for (const auto& [name, val] : full_solution) {
            if (var_set.count(name)) {
                filtered[name] = val;
            }
        }
        if (!vector_calculus_point_has_vars(filtered, vars)) {
            return LagrangeResult::failure(
                CasErrc::Inconclusive,
                "Lagrange candidate omits one or more variables",
                operation);
        }
        verified.push_back(std::move(filtered));
    }

    return LagrangeResult::success(std::move(verified));
}


static LagrangeResult lagrange_multipliers_strict(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars, ComputationContext& context, const std::string& operation)
{
    const size_t n = vars.size();
    const size_t m = constraints.size();
    auto gradient = lagrange_gradient_strict(f, vars, operation);
    if (!gradient) {
        return LagrangeResult::failure(gradient.error());
    }
    std::vector<VectorField> grad_constraints;
    grad_constraints.reserve(m);
    for (const auto& constraint : constraints) {
        auto partials = lagrange_gradient_strict(constraint, vars, operation);
        if (!partials) {
            return LagrangeResult::failure(partials.error());
        }
        grad_constraints.push_back(std::move(partials.value()));
    }
    auto lambda_names = lagrange_multiplier_names(m);
    VectorField equations;
    equations.reserve(n + m);
    for (size_t i = 0; i < n; ++i) {
        auto equation = lagrange_stationarity_equation(
            gradient.value(), grad_constraints, lambda_names, i);
        auto checked = vector_calculus_simplify_strict(equation, operation,
            "Lagrange stationarity equation is outside the supported domain");
        if (!checked) {
            return LagrangeResult::failure(checked.error());
        }
        equations.push_back(std::move(checked.value()));
    }
    for (const auto& constraint : constraints) equations.push_back(constraint);
    std::vector<std::string> all_vars;
    all_vars.reserve(n + m);
    all_vars.insert(all_vars.end(), vars.begin(), vars.end());
    all_vars.insert(all_vars.end(), lambda_names.begin(), lambda_names.end());
    auto solutions = solve_lagrange_strict(equations, all_vars, context, operation);
    if (!solutions) {
        return solutions;
    }
    return verify_lagrange_solutions(solutions.value(), equations, all_vars, vars, context, operation);
}

LagrangeResult lagrange_multipliers_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    const std::string operation = "lagrange_multipliers";
    auto valid = vector_calculus_validate_expr_vars(f, vars, context, operation);
    if (!valid) {
        return LagrangeResult::failure(valid.error());
    }
    auto distinct = vector_calculus_validate_distinct_vars(vars, vars.size(), context, operation);
    if (!distinct) {
        return LagrangeResult::failure(distinct.error());
    }
    if (constraints.empty()) {
        return LagrangeResult::failure(CasErrc::InvalidArgument,
                                       "constraint list cannot be empty",
                                       operation);
    }
    for (const auto& constraint : constraints) {
        if (!constraint || !LMCAS::detail::node(constraint)) {
            return LagrangeResult::failure(CasErrc::InvalidArgument,
                                           "constraint expressions cannot be null",
                                           operation);
        }
    }
    auto step = context.consume_steps(vars.size() * constraints.size() * 6 +
                                      vars.size() * 6 + constraints.size() * 4 + 8,
                                      operation);
    if (!step) {
        return LagrangeResult::failure(step.error());
    }

    try {
        return lagrange_multipliers_strict(f, constraints, vars, context, operation);
    } catch (const std::bad_alloc&) {
        return LagrangeResult::failure(CasErrc::ResourceLimit,
                                       "Lagrange multiplier allocation failed",
                                       operation);
    } catch (const std::exception& e) {
        return LagrangeResult::failure(CasErrc::InternalInvariant,
                                       e.what(), operation);
    }
}

LagrangeResult lagrange_multipliers_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return lagrange_multipliers_checked(f, constraints, vars, context);
}
std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> lagrange_multipliers(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars)
{
    auto result = lagrange_multipliers_checked(f, constraints, vars);
    return result ? std::move(result.value()) :
        std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>{};
}

}
