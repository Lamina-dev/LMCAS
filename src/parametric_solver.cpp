#include "parametric_solver.hpp"
#include "solver.hpp"
#include "poly_utils.hpp"
#include "symbolic_ast.hpp"
#include <algorithm>
#include <set>

namespace lamina {

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

    auto result = expr;
    for (const auto& var : unknowns) {
        auto coeff = extract_coefficient(expr, var);
        if (!coeff->is_zero()) {

        }
    }

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
                if (other_var == var) continue;
                auto mixed = first_deriv->differentiate(other_var)->simplify();
                if (!mixed->is_zero()) {
                    return false;
                }
            }
        }
    }
    return true;
}

std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
ParametricSolver::solve_linear_parametric(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>&)
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

        size_t pivot_row = current_row;
        bool found_pivot = false;
        while (pivot_row < m) {
            auto simplified = A[pivot_row][col]->simplify();
            A[pivot_row][col] = simplified;
            if (!simplified->is_zero()) {
                found_pivot = true;
                break;
            }
            pivot_row++;
        }

        if (!found_pivot) {

            continue;
        }

        if (pivot_row != current_row) {
            std::swap(A[pivot_row], A[current_row]);
            std::swap(b[pivot_row], b[current_row]);
        }

        auto pivot = A[current_row][col];

        for (size_t r = 0; r < m; ++r) {
            if (r == current_row) continue;
            auto entry = A[r][col]->simplify();
            if (entry->is_zero()) continue;

            auto factor = SymbolicExpr::divide(entry, pivot)->simplify();

            for (size_t k = col; k < n; ++k) {
                auto term = SymbolicExpr::multiply(factor, A[current_row][k]);
                A[r][k] = SymbolicExpr::add(A[r][k], SymbolicExpr::multiply(term, SymbolicExpr::number(-1)))->simplify();
            }

            auto b_term = SymbolicExpr::multiply(factor, b[current_row]);
            b[r] = SymbolicExpr::add(b[r], SymbolicExpr::multiply(b_term, SymbolicExpr::number(-1)))->simplify();
        }

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
            if (pivot_col_set.find(j) != pivot_col_set.end()) continue;
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

    return { solution };
}

std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
ParametricSolver::solve_polynomial_parametric(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>&)
{
    if (equations.empty() || unknowns.empty()) {

        std::map<std::string, std::shared_ptr<SymbolicExpr>> solution;
        for (const auto& var : unknowns) {
            solution[var] = SymbolicExpr::variable(var);
        }
        return { solution };
    }

    std::vector<SymbolicExpr> input_polys;
    input_polys.reserve(equations.size());
    for (const auto& eq : equations) {
        if (eq && lamina::detail::node(eq)) {
            auto simplified = eq->simplify();
            if (simplified && !simplified->is_zero()) {
                input_polys.push_back(*simplified);
            }
        }
    }

    if (input_polys.empty()) {

        std::map<std::string, std::shared_ptr<SymbolicExpr>> solution;
        for (const auto& var : unknowns) {
            solution[var] = SymbolicExpr::variable(var);
        }
        return { solution };
    }

    auto G_basis = Solver::groebner_basis(input_polys, unknowns);

    std::vector<std::shared_ptr<SymbolicExpr>> basis;
    basis.reserve(G_basis.size());
    for (const auto& g : G_basis) {
        auto g_ptr = lamina::detail::make_expression_ptr(g);
        auto simp = g_ptr->simplify();
        if (simp && !simp->is_zero()) {
            basis.push_back(simp);
        }
    }

    if (unknowns.empty()) {
        for (const auto& p : basis) {
            if (p->is_number() && !p->is_zero()) return {};
        }
        return { {} };
    }

    for (const auto& p : basis) {
        bool depends_on_unknown = false;
        for (const auto& var : unknowns) {
            if (contains(*p, var)) {
                depends_on_unknown = true;
                break;
            }
        }
        if (!depends_on_unknown) {

            if (!p->is_zero()) {
                return {};
            }
        }
    }

    auto substitute_all = [&](const std::shared_ptr<SymbolicExpr>& expr,
                              const std::map<std::string, std::shared_ptr<SymbolicExpr>>& subs)
        -> std::shared_ptr<SymbolicExpr> {
        auto res = expr;
        for (const auto& [name, val] : subs) {
            res = res->substitute(name, val);
            if (!res) return nullptr;
        }
        return res->simplify();
    };

    auto solve_rec = [&](auto&& self, int var_pos,
                         const std::map<std::string, std::shared_ptr<SymbolicExpr>>& partial)
        -> std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> {

        std::vector<std::shared_ptr<SymbolicExpr>> reduced;
        reduced.reserve(basis.size());

        for (const auto& p : basis) {
            auto r = substitute_all(p, partial);
            if (!r) continue;
            if (r->is_zero()) continue;

            bool depends = false;
            for (int i = 0; i <= var_pos && i < (int)unknowns.size(); ++i) {
                if (contains(*r, unknowns[i])) {
                    depends = true;
                    break;
                }
            }

            if (!depends) {

                if (r->is_number() && !r->is_zero()) {
                    return {};
                }

                continue;
            }

            reduced.push_back(r);
        }

        if (var_pos < 0) {
            return { partial };
        }

        const auto& curr_var = unknowns[var_pos];
        bool curr_var_appears = false;
        std::shared_ptr<SymbolicExpr> target = nullptr;
        int best_deg = std::numeric_limits<int>::max();

        for (const auto& r : reduced) {
            if (!contains(*r, curr_var)) continue;
            curr_var_appears = true;

            bool has_other = false;
            for (int i = 0; i < var_pos; ++i) {
                if (contains(*r, unknowns[i])) {
                    has_other = true;
                    break;
                }
            }
            if (has_other) continue;

            auto poly = symbolic_to_poly<SymbolicPolyCoeff>(r, curr_var);
            int deg = poly.degree();
            if (deg >= 1 && deg < best_deg) {
                best_deg = deg;
                target = r;
            }
        }

        if (!curr_var_appears) {
            auto next_partial = partial;
            next_partial[curr_var] = SymbolicExpr::variable(curr_var);
            return self(self, var_pos - 1, next_partial);
        }

        if (!target) {
            auto next_partial = partial;
            next_partial[curr_var] = SymbolicExpr::variable(curr_var);
            return self(self, var_pos - 1, next_partial);
        }

        auto roots = SymbolicExpr::solve(target, curr_var);
        if (roots.empty()) {

            auto next_partial = partial;
            next_partial[curr_var] = SymbolicExpr::variable(curr_var);
            return self(self, var_pos - 1, next_partial);
        }

        std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> results;
        for (const auto& r : roots) {
            auto next_partial = partial;
            next_partial[curr_var] = r;
            auto sub_res = self(self, var_pos - 1, next_partial);
            results.insert(results.end(), sub_res.begin(), sub_res.end());
        }
        return results;
    };

    std::map<std::string, std::shared_ptr<SymbolicExpr>> empty;
    return solve_rec(solve_rec, static_cast<int>(unknowns.size()) - 1, empty);
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

    std::vector<std::string> effective_params;
    std::set<std::string> unknown_set(unknowns.begin(), unknowns.end());
    for (const auto& p : parameters) {
        if (unknown_set.find(p) == unknown_set.end()) {
            effective_params.push_back(p);
        }
    }

    if (is_linear_in_unknowns(equations, unknowns)) {
        return solve_linear_parametric(equations, unknowns, effective_params);
    }

    return solve_polynomial_parametric(equations, unknowns, effective_params);
}

static std::shared_ptr<SymbolicExpr> compute_determinant(
    const std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& matrix,
    size_t n)
{
    if (n == 0) return SymbolicExpr::number(0);
    if (n == 1) return matrix[0][0];

    if (n == 2) {
        auto ad = SymbolicExpr::multiply(matrix[0][0], matrix[1][1]);
        auto bc = SymbolicExpr::multiply(matrix[0][1], matrix[1][0]);

        auto neg_bc = SymbolicExpr::multiply(bc, SymbolicExpr::number(-1));
        return SymbolicExpr::add(ad, neg_bc)->simplify();
    }

    if (n == 3) {

        auto a = matrix[0][0], b = matrix[0][1], c = matrix[0][2];
        auto d = matrix[1][0], e = matrix[1][1], f = matrix[1][2];
        auto g = matrix[2][0], h = matrix[2][1], i = matrix[2][2];

        auto ei = SymbolicExpr::multiply(e, i);
        auto fh = SymbolicExpr::multiply(f, h);
        auto ei_fh = SymbolicExpr::add(ei, SymbolicExpr::multiply(fh, SymbolicExpr::number(-1)));

        auto di = SymbolicExpr::multiply(d, i);
        auto fg = SymbolicExpr::multiply(f, g);
        auto di_fg = SymbolicExpr::add(di, SymbolicExpr::multiply(fg, SymbolicExpr::number(-1)));

        auto dh = SymbolicExpr::multiply(d, h);
        auto eg = SymbolicExpr::multiply(e, g);
        auto dh_eg = SymbolicExpr::add(dh, SymbolicExpr::multiply(eg, SymbolicExpr::number(-1)));

        auto term1 = SymbolicExpr::multiply(a, ei_fh);
        auto term2 = SymbolicExpr::multiply(b, di_fg);
        auto term3 = SymbolicExpr::multiply(c, dh_eg);

        auto neg_term2 = SymbolicExpr::multiply(term2, SymbolicExpr::number(-1));
        auto sum12 = SymbolicExpr::add(term1, neg_term2);
        return SymbolicExpr::add(sum12, term3)->simplify();
    }

    auto det = SymbolicExpr::number(0);
    for (size_t j = 0; j < n; ++j) {

        std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> minor(n - 1, std::vector<std::shared_ptr<SymbolicExpr>>(n - 1));
        for (size_t row = 1; row < n; ++row) {
            size_t minor_col = 0;
            for (size_t col = 0; col < n; ++col) {
                if (col == j) continue;
                minor[row - 1][minor_col] = matrix[row][col];
                minor_col++;
            }
        }

        auto cofactor = compute_determinant(minor, n - 1);
        auto term = SymbolicExpr::multiply(matrix[0][j], cofactor);

        if (j % 2 == 0) {
            det = SymbolicExpr::add(det, term);
        } else {
            det = SymbolicExpr::add(det, SymbolicExpr::multiply(term, SymbolicExpr::number(-1)));
        }
    }
    return det->simplify();
}

static bool depends_on_parameters(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::vector<std::string>& parameters)
{
    if (!expr) return false;
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
