#include "parametric_solver.hpp"
#include "internal/exact_matrix.hpp"
#include "solver.hpp"
#include "solve_strategies.hpp"
#include "poly_utils.hpp"
#include "internal/symbolic_ast.hpp"
#include <algorithm>
#include <set>

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

static size_t find_parametric_pivot(ParametricMatrix& matrix, size_t row, size_t column) {
    while (row < matrix.size()) {
        auto simplified = matrix[row][column]->simplify();
        matrix[row][column] = simplified;
        if (!simplified->is_zero()) {
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
    const std::vector<std::string>& unknowns) {
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
            if (!coeff->is_zero()) {

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
    const std::vector<std::string>& unknowns)
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

        size_t pivot_row = find_parametric_pivot(A, current_row, col);
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
        if (!b_simplified->is_zero()) {

            return {};
        }
    }

    return {parametric_linear_solution(A, b, pivot_cols, pivot_rows, unknowns)};
}


ParametricSolutionsResult
ParametricSolver::solve_polynomial_parametric_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>&,
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
        return solve_polynomial_parametric_impl(
            equations, unknowns, context);
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
    const std::vector<std::string>&)
{
    if (equations.empty() || unknowns.empty()) {
        return {};
    }


    if (is_linear_in_unknowns(equations, unknowns)) {
        return solve_linear_parametric(equations, unknowns);
    }

    ComputationContext context;
    auto solved = solve_polynomial_parametric_impl(
        equations, unknowns, context);
    return solved ? std::move(solved.value()) : ParametricSolutionList{};
}

static std::shared_ptr<SymbolicExpr> compute_determinant(
    const std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& matrix,
    size_t dimension)
{
    if (dimension == 0 || matrix.size() != dimension) {
        return nullptr;
    }
    detail::ExactMatrixData exact{dimension, dimension, {}};
    exact.entries.reserve(dimension * dimension);
    for (const auto& row : matrix) {
        if (row.size() != dimension) {
            return nullptr;
        }
        for (const auto& entry : row) {
            if (!entry) {
                return nullptr;
            }
            exact.entries.push_back(entry);
        }
    }
    ComputationContext context;
    auto determinant = detail::determinant_exact(
        exact, context, "parametric_determinant");
    return determinant ? std::move(determinant.value()) : nullptr;
}

static bool depends_on_parameters(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::vector<std::string>& parameters)
{
    if (!expr) {
        return false;
    }
    for (const auto& p : parameters) {
        if (contains(*expr, p)) {
            return true;
        }
    }
    return false;
}

PiecewiseSolution ParametricSolver::solve_system_piecewise(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>& parameters)
{
    PiecewiseSolution result;
    size_t m = equations.size();
    size_t n = unknowns.size();

    if (m == 0 || n == 0) {
        return result;
    }

    std::vector<std::string> effective_params;
    std::set<std::string> unknown_set(unknowns.begin(), unknowns.end());
    for (const auto& p : parameters) {
        if (unknown_set.find(p) == unknown_set.end()) {
            effective_params.push_back(p);
        }
    }

    if (m != n) {
        auto solutions = solve_system(equations, unknowns, parameters);
        if (!solutions.empty()) {
            PiecewiseSolution::Case generic_case;
            generic_case.condition = SymbolicExpr::number(1);
            generic_case.solutions = solutions;
            result.cases.push_back(generic_case);
        }
        return result;
    }

    if (!is_linear_in_unknowns(equations, unknowns)) {
        auto solutions = solve_system(equations, unknowns, parameters);
        if (!solutions.empty()) {
            PiecewiseSolution::Case generic_case;
            generic_case.condition = SymbolicExpr::number(1);
            generic_case.solutions = solutions;
            result.cases.push_back(generic_case);
        }
        return result;
    }

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            A[i][j] = extract_coefficient(equations[i], unknowns[j]);
        }
    }

    auto det = compute_determinant(A, n);
    if (!det) {
        auto solutions = solve_system(equations, unknowns, parameters);
        if (!solutions.empty()) {
            PiecewiseSolution::Case generic_case;
            generic_case.condition = SymbolicExpr::number(1);
            generic_case.solutions = std::move(solutions);
            result.cases.push_back(std::move(generic_case));
        }
        return result;
    }

    if (det->is_zero()) {

        auto solutions = solve_system(equations, unknowns, parameters);
        PiecewiseSolution::Case degenerate_case;

        degenerate_case.condition = SymbolicExpr::number(0);
        degenerate_case.solutions = solutions;
        result.cases.push_back(degenerate_case);
    } else if (!depends_on_parameters(det, effective_params)) {

        auto solutions = solve_system(equations, unknowns, parameters);
        PiecewiseSolution::Case generic_case;

        generic_case.condition = det;
        generic_case.solutions = solutions;
        result.cases.push_back(generic_case);
    } else {

        auto generic_solutions = solve_system(equations, unknowns, parameters);
        PiecewiseSolution::Case generic_case;
        generic_case.condition = det;
        generic_case.solutions = generic_solutions;
        result.cases.push_back(generic_case);

        PiecewiseSolution::Case degenerate_case;
        degenerate_case.condition = SymbolicExpr::number(0);

        degenerate_case.solutions = generic_solutions;
        result.cases.push_back(degenerate_case);
    }

    return result;
}

}
