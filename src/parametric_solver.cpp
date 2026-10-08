#include "parametric_solver.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/assumption_facts.hpp"
#include "residual_verification.hpp"
#include <algorithm>
#include <set>
#include <optional>

namespace LMCAS {

static std::shared_ptr<SymbolicExpr> extract_coefficient(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var)
{

    return expr->differentiate(var)->simplify();
}

static std::shared_ptr<SymbolicExpr> extract_constant_term(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::vector<std::string>& unknowns)
{


    auto constant = expr;
    for (const auto& var : unknowns) {
        constant = constant->substitute(var, SymbolicExpr::number(0));
    }
    return constant->simplify();
}

bool ParametricSolver::is_linear_in_unknowns(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns)
{
    for (const auto& eq : equations) {
        for (const auto& var : unknowns) {

            auto first_deriv = eq->differentiate(var)->simplify();
            auto second_deriv = first_deriv->differentiate(var)->simplify();
            if (!second_deriv->is_zero()) {
                return false;
            }

            for (const auto& other_var : unknowns) {
                if (other_var == var) {
                    continue;
                }
                auto mixed = first_deriv->differentiate(other_var)->simplify();
                if (!mixed->is_zero()) {
                    return false;
                }
            }
        }
    }
    return true;
}

using ParametricMatrix = std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>;
using ParametricVector = std::vector<std::shared_ptr<SymbolicExpr>>;

static bool parametric_residual_zero(const std::shared_ptr<SymbolicExpr>& expression,
                                     ComputationContext& context) {
    auto residual = check_zero_residual(expression, context);
    if (!residual) throw residual.error();
    if (std::holds_alternative<ProvedZeroResidual>(residual.value())) return true;
    throw CasError{CasErrc::Inconclusive,
        "parameter coefficient or constraint is undecided", "solve.parametric"};
}

static Tribool parametric_nonzero_value(const std::shared_ptr<SymbolicExpr>& expression,
                                        ComputationContext& context) {
    std::optional<detail::AssumptionFacts> assumptions;
    if (context.assumptions()) assumptions.emplace(*context.assumptions());
    const FactsQuery& facts = assumptions ? static_cast<const FactsQuery&>(*assumptions) : detail::no_facts();
    auto nonzero = detail::query_nonzero_value(detail::node(expression), facts, Domain::Complex, context);
    if (!nonzero) throw nonzero.error();
    return nonzero.value();
}

static bool parametric_zero(const std::shared_ptr<SymbolicExpr>& expression,
                            ComputationContext& context) {
    if (expression->is_zero()) return true;
    const auto nonzero = parametric_nonzero_value(expression, context);
    if (nonzero == Tribool::True) return false;
    if (nonzero == Tribool::False) return true;
    return parametric_residual_zero(expression, context);
}

static size_t find_parametric_pivot(ParametricMatrix& matrix, size_t row, size_t column,
                                    ComputationContext& context) {
    while (row < matrix.size()) {
        auto simplified = matrix[row][column]->simplify();
        matrix[row][column] = simplified;
        if (!parametric_zero(simplified, context)) {
            break;
        }
        ++row;
    }
    return row;
}

static void eliminate_parametric_column(ParametricMatrix& A, ParametricVector& b,
                                        size_t current_row, size_t col, size_t n) {
    size_t m = A.size();
    auto pivot = A[current_row][col];

    for (size_t r = 0; r < m; ++r) {
        if (r == current_row) {
            continue;
        }
        auto entry = A[r][col]->simplify();
        if (entry->is_zero()) {
            continue;
        }

        auto factor = SymbolicExpr::divide(entry, pivot)->simplify();

        for (size_t k = col; k < n; ++k) {
            auto term = SymbolicExpr::multiply(factor, A[current_row][k]);
            A[r][k] = SymbolicExpr::add(A[r][k], SymbolicExpr::multiply(term, SymbolicExpr::number(-1)))->simplify();
        }

        auto b_term = SymbolicExpr::multiply(factor, b[current_row]);
        b[r] = SymbolicExpr::add(b[r], SymbolicExpr::multiply(b_term, SymbolicExpr::number(-1)))->simplify();
    }
}

static std::map<std::string, std::shared_ptr<SymbolicExpr>> parametric_linear_solution(
    const ParametricMatrix& A, const ParametricVector& b,
    const std::vector<size_t>& pivot_cols, const std::vector<size_t>& pivot_rows,
    const std::vector<std::string>& unknowns, ComputationContext& context) {
    size_t n = unknowns.size();
    std::map<std::string, std::shared_ptr<SymbolicExpr>> solution;

    std::set<size_t> pivot_col_set(pivot_cols.begin(), pivot_cols.end());
    for (size_t col = 0; col < n; ++col) {
        if (pivot_col_set.find(col) == pivot_col_set.end()) {

            solution[unknowns[col]] = SymbolicExpr::variable(unknowns[col]);
        }
    }

    for (size_t k = 0; k < pivot_cols.size(); ++k) {
        size_t col = pivot_cols[k];
        size_t row = pivot_rows[k];

        auto value = b[row];

        for (size_t j = col + 1; j < n; ++j) {
            if (pivot_col_set.find(j) != pivot_col_set.end()) {
                continue;
            }
            auto coeff = A[row][j]->simplify();
            if (!parametric_zero(coeff, context)) {

                auto term = SymbolicExpr::multiply(coeff, SymbolicExpr::variable(unknowns[j]));
                value = SymbolicExpr::add(value, SymbolicExpr::multiply(term, SymbolicExpr::number(-1)));
            }
        }

        auto pivot_val = A[row][col];
        value = SymbolicExpr::divide(value, pivot_val)->simplify();

        solution[unknowns[col]] = value;
    }
    return solution;
}

std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
ParametricSolver::solve_linear_parametric(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    ComputationContext& context)
{
    size_t m = equations.size();
    size_t n = unknowns.size();

    if (m == 0 || n == 0) {

        std::map<std::string, std::shared_ptr<SymbolicExpr>> solution;
        for (const auto& var : unknowns) {
            solution[var] = SymbolicExpr::variable(var);
        }
        return { solution };
    }

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A(m, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    std::vector<std::shared_ptr<SymbolicExpr>> b(m);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            A[i][j] = extract_coefficient(equations[i], unknowns[j]);
        }

        auto constant = extract_constant_term(equations[i], unknowns);
        b[i] = SymbolicExpr::multiply(constant, SymbolicExpr::number(-1))->simplify();
    }

    std::vector<size_t> pivot_cols;
    std::vector<size_t> pivot_rows;
    size_t current_row = 0;

    for (size_t col = 0; col < n && current_row < m; ++col) {

        size_t pivot_row = find_parametric_pivot(A, current_row, col, context);
        if (pivot_row == m) {
            continue;
        }

        if (pivot_row != current_row) {
            std::swap(A[pivot_row], A[current_row]);
            std::swap(b[pivot_row], b[current_row]);
        }

        eliminate_parametric_column(A, b, current_row, col, n);

        pivot_cols.push_back(col);
        pivot_rows.push_back(current_row);
        current_row++;
    }

    for (size_t r = current_row; r < m; ++r) {
        auto b_simplified = b[r]->simplify();
        if (!parametric_zero(b_simplified, context)) {

            return {};
        }
    }

    return {parametric_linear_solution(A, b, pivot_cols, pivot_rows, unknowns, context)};
}


ParametricSolutionsResult
ParametricSolver::solve_polynomial_parametric_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>& parameters,
    ComputationContext& context)
{
    constexpr const char* operation = "solve_polynomial_parametric";
    if (equations.empty() || unknowns.empty()) {
        return ParametricSolutionsResult::failure(
            CasErrc::InvalidArgument,
            "parametric polynomial solve requires equations and unknowns",
            operation);
    }
    auto budget = context.consume_steps(
        equations.size() * unknowns.size() + 1, operation);
    if (!budget) {
        return ParametricSolutionsResult::failure(budget.error());
    }
    try {
        if (is_linear_in_unknowns(equations, unknowns))
            return solve_linear_parametric(equations, unknowns, context);
        return solve_polynomial_parametric_impl(
            equations, unknowns, parameters, context);
    } catch (const CasError& error) {
        return ParametricSolutionsResult::failure(error);
    } catch (const std::bad_alloc&) {
        return ParametricSolutionsResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while solving parametric polynomial system",
            operation);
    } catch (const std::exception& ex) {
        return ParametricSolutionsResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

ParametricSolutionsResult
ParametricSolver::solve_polynomial_parametric_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>& parameters)
{
    ComputationContext context;
    return solve_polynomial_parametric_checked(
        equations, unknowns, parameters, context);
}

std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
ParametricSolver::solve_system(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>& parameters)
{
    if (equations.empty() || unknowns.empty()) {
        return {};
    }


    ComputationContext context;
    if (is_linear_in_unknowns(equations, unknowns)) {
        return solve_linear_parametric(equations, unknowns, context);
    }
    auto solved = solve_polynomial_parametric_checked(
        equations, unknowns, parameters, context);
    if (!solved) throw solved.error();
    return std::move(solved.value());
}


PiecewiseSolution ParametricSolver::solve_system_piecewise(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>& parameters)
{
    PiecewiseSolution result;
    if (equations.empty() || unknowns.empty()) return result;
    auto solutions = solve_system(equations, unknowns, parameters);
    if (!solutions.empty()) {
        PiecewiseSolution::Case only_case;
        only_case.condition = SymbolicExpr::number(1);
        only_case.solutions = std::move(solutions);
        result.cases.push_back(std::move(only_case));
    }
    return result;
}

}
