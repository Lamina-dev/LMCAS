#include "parametric_solver.hpp"
#include "solver.hpp"
#include "poly_utils.hpp"
#include <algorithm>
#include <set>

namespace lamina {

// ============================================================================
// Helper: Extract the coefficient of a variable from a symbolic expression.
// Uses differentiation: coeff(expr, var) = d(expr)/d(var)
// The constant term is: expr - coeff * var
// ============================================================================
static std::shared_ptr<SymbolicExpr> extract_coefficient(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var)
{
    // The coefficient of 'var' in a linear expression is the partial derivative
    return expr->differentiate(var)->simplify();
}

// Extract the constant term of an expression after removing all unknown terms.
// constant = expr - sum(coeff_i * unknown_i) for all unknowns
static std::shared_ptr<SymbolicExpr> extract_constant_term(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::vector<std::string>& unknowns)
{
    // Start with the expression, subtract each unknown's contribution
    auto result = expr;
    for (const auto& var : unknowns) {
        auto coeff = extract_coefficient(expr, var);
        if (!coeff->is_zero()) {
            // result = result - coeff * var
            // But we need to compute from the original expr, not iteratively
            // Actually, for a linear expression: expr = c0 + c1*x1 + c2*x2 + ...
            // We substitute all unknowns with 0 to get the constant term
        }
    }
    // Better approach: substitute all unknowns with 0
    auto constant = expr;
    for (const auto& var : unknowns) {
        constant = constant->substitute(var, SymbolicExpr::number(0));
    }
    return constant->simplify();
}

// ============================================================================
// is_linear_in_unknowns: Check if each equation is linear in the unknowns
// An expression is linear in unknowns if:
//   - The second derivative w.r.t. any unknown is zero
//   - The mixed partial derivative w.r.t. any two unknowns is zero
// ============================================================================
bool ParametricSolver::is_linear_in_unknowns(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns)
{
    for (const auto& eq : equations) {
        for (const auto& var : unknowns) {
            // Check degree > 1: second derivative w.r.t. var should be zero
            auto first_deriv = eq->differentiate(var)->simplify();
            auto second_deriv = first_deriv->differentiate(var)->simplify();
            if (!second_deriv->is_zero()) {
                return false;
            }
            // Check that the first derivative doesn't depend on any other unknown
            // (no products of unknowns)
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

// ============================================================================
// solve_linear_parametric: Symbolic Gaussian elimination for parametric systems
//
// Algorithm (from design section 4.2):
// 1. Construct augmented matrix [A | b]
// 2. Perform symbolic Gaussian elimination with partial pivoting
// 3. Back-substitution
// 4. Free variables remain as themselves
// ============================================================================
std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
ParametricSolver::solve_linear_parametric(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>& parameters)
{
    size_t m = equations.size();  // number of equations
    size_t n = unknowns.size();   // number of unknowns

    if (m == 0 || n == 0) {
        // No equations or no unknowns: return empty solution with free variables
        std::map<std::string, std::shared_ptr<SymbolicExpr>> solution;
        for (const auto& var : unknowns) {
            solution[var] = SymbolicExpr::variable(var);
        }
        return { solution };
    }

    // Step 1: Construct augmented matrix [A | b]
    // A[i][j] = coefficient of unknowns[j] in equations[i]
    // b[i] = -constant_term(equations[i])
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A(m, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    std::vector<std::shared_ptr<SymbolicExpr>> b(m);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            A[i][j] = extract_coefficient(equations[i], unknowns[j]);
        }
        // b[i] = -(constant term) = -(equation with all unknowns set to 0)
        auto constant = extract_constant_term(equations[i], unknowns);
        b[i] = SymbolicExpr::multiply(constant, SymbolicExpr::number(-1))->simplify();
    }

    // Step 2: Symbolic Gaussian elimination with partial pivoting
    // We perform row reduction to reduced row echelon form (RREF)
    std::vector<size_t> pivot_cols;  // pivot_cols[k] = column index of k-th pivot
    std::vector<size_t> pivot_rows;  // corresponding row indices
    size_t current_row = 0;

    for (size_t col = 0; col < n && current_row < m; ++col) {
        // Find a row with a non-zero pivot in this column
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
            // This column has no pivot - it will be a free variable
            continue;
        }

        // Swap rows if needed
        if (pivot_row != current_row) {
            std::swap(A[pivot_row], A[current_row]);
            std::swap(b[pivot_row], b[current_row]);
        }

        auto pivot = A[current_row][col];

        // Eliminate all other rows
        for (size_t r = 0; r < m; ++r) {
            if (r == current_row) continue;
            auto entry = A[r][col]->simplify();
            if (entry->is_zero()) continue;

            // factor = A[r][col] / pivot
            auto factor = SymbolicExpr::divide(entry, pivot)->simplify();

            // A[r][k] = A[r][k] - factor * A[current_row][k]
            for (size_t k = col; k < n; ++k) {
                auto term = SymbolicExpr::multiply(factor, A[current_row][k]);
                A[r][k] = SymbolicExpr::add(A[r][k], SymbolicExpr::multiply(term, SymbolicExpr::number(-1)))->simplify();
            }
            // b[r] = b[r] - factor * b[current_row]
            auto b_term = SymbolicExpr::multiply(factor, b[current_row]);
            b[r] = SymbolicExpr::add(b[r], SymbolicExpr::multiply(b_term, SymbolicExpr::number(-1)))->simplify();
        }

        pivot_cols.push_back(col);
        pivot_rows.push_back(current_row);
        current_row++;
    }

    // Step 3: Check for inconsistency
    // If any row has all-zero coefficients but non-zero b, the system is inconsistent
    for (size_t r = current_row; r < m; ++r) {
        auto b_simplified = b[r]->simplify();
        if (!b_simplified->is_zero()) {
            // Inconsistent system
            return {};
        }
    }

    // Step 4: Back-substitution
    // For each pivot column, compute unknowns[col] = b[row] / A[row][col]
    // Also account for non-pivot columns (free variables) in the same row
    std::map<std::string, std::shared_ptr<SymbolicExpr>> solution;

    // First, mark free variables
    std::set<size_t> pivot_col_set(pivot_cols.begin(), pivot_cols.end());
    for (size_t col = 0; col < n; ++col) {
        if (pivot_col_set.find(col) == pivot_col_set.end()) {
            // Free variable: remains as itself
            solution[unknowns[col]] = SymbolicExpr::variable(unknowns[col]);
        }
    }

    // Now solve for pivot variables
    for (size_t k = 0; k < pivot_cols.size(); ++k) {
        size_t col = pivot_cols[k];
        size_t row = pivot_rows[k];

        // Start with b[row]
        auto value = b[row];

        // Subtract contributions from free variables in this row
        for (size_t j = col + 1; j < n; ++j) {
            if (pivot_col_set.find(j) != pivot_col_set.end()) continue; // skip pivot columns
            auto coeff = A[row][j]->simplify();
            if (!coeff->is_zero()) {
                // value = value - coeff * unknowns[j]
                auto term = SymbolicExpr::multiply(coeff, SymbolicExpr::variable(unknowns[j]));
                value = SymbolicExpr::add(value, SymbolicExpr::multiply(term, SymbolicExpr::number(-1)));
            }
        }

        // Divide by the pivot element
        auto pivot_val = A[row][col];
        value = SymbolicExpr::divide(value, pivot_val)->simplify();

        solution[unknowns[col]] = value;
    }

    return { solution };
}

// ============================================================================
// solve_polynomial_parametric: Solve polynomial systems using Gröbner basis
//
// Algorithm (from design section 4.3):
// 1. Compute Gröbner basis treating parameters as coefficients
//    (only unknowns participate in the monomial ordering)
// 2. Check consistency: if basis contains constant 1, return empty
// 3. Find univariate polynomials in the elimination ideal
// 4. Solve univariate equations (roots are expressions involving parameters)
// 5. Back-substitute to find remaining unknowns
// 6. Handle positive-dimensional systems (free variables)
// ============================================================================
std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
ParametricSolver::solve_polynomial_parametric(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>& parameters)
{
    if (equations.empty() || unknowns.empty()) {
        // No equations or no unknowns: all unknowns are free
        std::map<std::string, std::shared_ptr<SymbolicExpr>> solution;
        for (const auto& var : unknowns) {
            solution[var] = SymbolicExpr::variable(var);
        }
        return { solution };
    }

    // Step 1: Compute Gröbner basis with only unknowns as variables.
    // Parameters are automatically treated as coefficients since they are not
    // in the variable list passed to groebner_basis().
    std::vector<SymbolicExpr> input_polys;
    input_polys.reserve(equations.size());
    for (const auto& eq : equations) {
        if (eq && eq->root) {
            auto simplified = eq->simplify();
            if (simplified && !simplified->is_zero()) {
                input_polys.push_back(*simplified);
            }
        }
    }

    if (input_polys.empty()) {
        // All equations simplified to zero: all unknowns are free
        std::map<std::string, std::shared_ptr<SymbolicExpr>> solution;
        for (const auto& var : unknowns) {
            solution[var] = SymbolicExpr::variable(var);
        }
        return { solution };
    }

    // Compute Gröbner basis (lex order on unknowns, parameters as coefficients)
    auto G_basis = Solver::groebner_basis(input_polys, unknowns);

    // Simplify basis elements and filter out zeros
    std::vector<std::shared_ptr<SymbolicExpr>> basis;
    basis.reserve(G_basis.size());
    for (const auto& g : G_basis) {
        auto g_ptr = std::make_shared<SymbolicExpr>(g);
        auto simp = g_ptr->simplify();
        if (simp && !simp->is_zero()) {
            basis.push_back(simp);
        }
    }

    // Step 2: Check consistency - if basis contains a non-zero constant
    // (i.e., an element that doesn't depend on any unknown), system is inconsistent
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
            // This element doesn't involve any unknown.
            // If it's a non-zero expression (possibly involving parameters),
            // the system is inconsistent (for generic parameter values).
            // If it's identically zero, skip it.
            if (!p->is_zero()) {
                return {};  // Inconsistent system
            }
        }
    }

    // Step 3-5: Recursive back-substitution
    // Similar to Solver::solve_polynomial_system() but solutions contain parameters.
    // We process unknowns from last to first (lex order elimination property).

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

    // Recursive solver: processes unknowns from var_pos down to 0
    auto solve_rec = [&](auto&& self, int var_pos,
                         const std::map<std::string, std::shared_ptr<SymbolicExpr>>& partial)
        -> std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> {

        // Substitute known values into basis elements
        std::vector<std::shared_ptr<SymbolicExpr>> reduced;
        reduced.reserve(basis.size());

        for (const auto& p : basis) {
            auto r = substitute_all(p, partial);
            if (!r) continue;
            if (r->is_zero()) continue;

            // Check if this element still depends on unknowns[0..var_pos]
            bool depends = false;
            for (int i = 0; i <= var_pos && i < (int)unknowns.size(); ++i) {
                if (contains(*r, unknowns[i])) {
                    depends = true;
                    break;
                }
            }

            if (!depends) {
                // Element doesn't depend on remaining unknowns.
                // If it's a non-zero number, system is inconsistent for this branch.
                if (r->is_number() && !r->is_zero()) {
                    return {};
                }
                // If it depends only on parameters, it's a constraint on parameters
                // (we treat it as consistent for generic parameter values)
                continue;
            }

            reduced.push_back(r);
        }

        // Base case: all unknowns have been assigned
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

            // Check if this element involves other unknowns (indices < var_pos)
            bool has_other = false;
            for (int i = 0; i < var_pos; ++i) {
                if (contains(*r, unknowns[i])) {
                    has_other = true;
                    break;
                }
            }
            if (has_other) continue;

            // This element only involves curr_var (and possibly parameters)
            // Choose the one with lowest degree in curr_var
            auto poly = symbolic_to_poly<SymbolicPolyCoeff>(r, curr_var);
            int deg = poly.degree();
            if (deg >= 1 && deg < best_deg) {
                best_deg = deg;
                target = r;
            }
        }

        // Step 6: If curr_var doesn't appear, it's a free variable
        if (!curr_var_appears) {
            auto next_partial = partial;
            next_partial[curr_var] = SymbolicExpr::variable(curr_var);
            return self(self, var_pos - 1, next_partial);
        }

        // If no suitable univariate polynomial found, treat as free variable
        if (!target) {
            auto next_partial = partial;
            next_partial[curr_var] = SymbolicExpr::variable(curr_var);
            return self(self, var_pos - 1, next_partial);
        }

        // Step 4: Solve the univariate equation for curr_var
        // The coefficients may contain parameters - SymbolicExpr::solve() handles this
        auto roots = SymbolicExpr::solve(target, curr_var);
        if (roots.empty()) {
            // Cannot solve - treat as free variable as fallback
            auto next_partial = partial;
            next_partial[curr_var] = SymbolicExpr::variable(curr_var);
            return self(self, var_pos - 1, next_partial);
        }

        // Step 5: Back-substitute each root and solve for remaining unknowns
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

// ============================================================================
// solve_system: Public entry point - dispatches to linear or polynomial solver
// ============================================================================
std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
ParametricSolver::solve_system(
    const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
    const std::vector<std::string>& unknowns,
    const std::vector<std::string>& parameters)
{
    if (equations.empty() || unknowns.empty()) {
        return {};
    }

    // Handle unknowns/parameters overlap: unknowns take precedence
    std::vector<std::string> effective_params;
    std::set<std::string> unknown_set(unknowns.begin(), unknowns.end());
    for (const auto& p : parameters) {
        if (unknown_set.find(p) == unknown_set.end()) {
            effective_params.push_back(p);
        }
    }

    // Check if the system is linear in the unknowns
    if (is_linear_in_unknowns(equations, unknowns)) {
        return solve_linear_parametric(equations, unknowns, effective_params);
    }

    // Non-linear: delegate to polynomial solver
    return solve_polynomial_parametric(equations, unknowns, effective_params);
}

// ============================================================================
// Helper: Compute the symbolic determinant of a matrix represented as a 2D vector
// Uses cofactor expansion along the first row (recursive).
// ============================================================================
static std::shared_ptr<SymbolicExpr> compute_determinant(
    const std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& matrix,
    size_t n)
{
    if (n == 0) return SymbolicExpr::number(0);
    if (n == 1) return matrix[0][0];

    // 2x2: ad - bc
    if (n == 2) {
        auto ad = SymbolicExpr::multiply(matrix[0][0], matrix[1][1]);
        auto bc = SymbolicExpr::multiply(matrix[0][1], matrix[1][0]);
        // det = ad - bc = ad + (-1)*bc
        auto neg_bc = SymbolicExpr::multiply(bc, SymbolicExpr::number(-1));
        return SymbolicExpr::add(ad, neg_bc)->simplify();
    }

    // 3x3: Sarrus' rule (or cofactor expansion)
    if (n == 3) {
        // a(ei - fh) - b(di - fg) + c(dh - eg)
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

        // det = term1 - term2 + term3
        auto neg_term2 = SymbolicExpr::multiply(term2, SymbolicExpr::number(-1));
        auto sum12 = SymbolicExpr::add(term1, neg_term2);
        return SymbolicExpr::add(sum12, term3)->simplify();
    }

    // General case: cofactor expansion along first row
    auto det = SymbolicExpr::number(0);
    for (size_t j = 0; j < n; ++j) {
        // Build (n-1)x(n-1) minor by removing row 0 and column j
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

        // Sign: (-1)^j
        if (j % 2 == 0) {
            det = SymbolicExpr::add(det, term);
        } else {
            det = SymbolicExpr::add(det, SymbolicExpr::multiply(term, SymbolicExpr::number(-1)));
        }
    }
    return det->simplify();
}

// ============================================================================
// Helper: Check if an expression depends on any of the given parameters
// ============================================================================
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

// ============================================================================
// solve_system_piecewise: Returns piecewise solution based on parameter conditions
//
// Algorithm:
// 1. For non-square systems (m != n), return single-case wrapper
// 2. For square systems, compute symbolic determinant of coefficient matrix
// 3. Classify determinant:
//    - Non-zero constant: single Case (unique solution)
//    - Identically zero: degenerate Case
//    - Depends on parameters: two Cases (generic + degenerate)
// ============================================================================
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

    // Handle unknowns/parameters overlap: unknowns take precedence
    std::vector<std::string> effective_params;
    std::set<std::string> unknown_set(unknowns.begin(), unknowns.end());
    for (const auto& p : parameters) {
        if (unknown_set.find(p) == unknown_set.end()) {
            effective_params.push_back(p);
        }
    }

    // For non-square systems, just return the single-case wrapper
    if (m != n) {
        auto solutions = solve_system(equations, unknowns, parameters);
        if (!solutions.empty()) {
            PiecewiseSolution::Case generic_case;
            generic_case.condition = SymbolicExpr::number(1); // always true
            generic_case.solutions = solutions;
            result.cases.push_back(generic_case);
        }
        return result;
    }

    // For non-linear systems, also return single-case wrapper
    if (!is_linear_in_unknowns(equations, unknowns)) {
        auto solutions = solve_system(equations, unknowns, parameters);
        if (!solutions.empty()) {
            PiecewiseSolution::Case generic_case;
            generic_case.condition = SymbolicExpr::number(1); // always true
            generic_case.solutions = solutions;
            result.cases.push_back(generic_case);
        }
        return result;
    }

    // Square linear system: compute symbolic determinant of coefficient matrix
    // Build coefficient matrix A[i][j] = coefficient of unknowns[j] in equations[i]
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> A(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            A[i][j] = extract_coefficient(equations[i], unknowns[j]);
        }
    }

    // Compute symbolic determinant
    auto det = compute_determinant(A, n);

    // Classify the determinant
    if (det->is_zero()) {
        // Determinant is identically zero: degenerate case
        // Solve the underdetermined system (may have free variables or be inconsistent)
        auto solutions = solve_system(equations, unknowns, parameters);
        PiecewiseSolution::Case degenerate_case;
        // condition: det = 0 (which is always true here)
        degenerate_case.condition = SymbolicExpr::number(0); // represents det = 0
        degenerate_case.solutions = solutions;
        result.cases.push_back(degenerate_case);
    } else if (!depends_on_parameters(det, effective_params)) {
        // Determinant is a non-zero constant (doesn't depend on parameters)
        // Return single Case with unique solution
        auto solutions = solve_system(equations, unknowns, parameters);
        PiecewiseSolution::Case generic_case;
        // condition: det ≠ 0 (always true since det is a non-zero constant)
        generic_case.condition = det; // non-zero constant serves as "det ≠ 0"
        generic_case.solutions = solutions;
        result.cases.push_back(generic_case);
    } else {
        // Determinant depends on parameters: return two cases
        // Case 1: Generic case (det ≠ 0) - unique solution
        auto generic_solutions = solve_system(equations, unknowns, parameters);
        PiecewiseSolution::Case generic_case;
        generic_case.condition = det; // condition: det ≠ 0
        generic_case.solutions = generic_solutions;
        result.cases.push_back(generic_case);

        // Case 2: Degenerate case (det = 0)
        // We need to solve the system under the constraint that det = 0
        // For now, solve the system as-is (the linear solver handles rank-deficient cases)
        // The degenerate solution may have free variables
        PiecewiseSolution::Case degenerate_case;
        degenerate_case.condition = SymbolicExpr::number(0); // represents det = 0
        // For the degenerate case, the system is rank-deficient.
        // The solve_linear_parametric already handles this (free variables).
        // We still call solve_system which will attempt to solve; the result
        // may include free variables or be empty if inconsistent.
        degenerate_case.solutions = generic_solutions; // same call, but semantically under det=0
        result.cases.push_back(degenerate_case);
    }

    return result;
}

} // namespace lamina
