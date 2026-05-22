// root_of_utils.cpp - RootOf symbolic representation utilities

#include "root_of_utils.hpp"
#include "newton_raphson.hpp"
#include "solve_polynomial.hpp"
#include "symbolic_ast.hpp"
#include <algorithm>
#include <cmath>
#include <complex>

namespace lamina {

// Helper: convert Polynomial<SymbolicPolyCoeff> to a SymbolicExpr
static std::shared_ptr<SymbolicExpr> symbolic_poly_to_expr(
    const Polynomial<SymbolicPolyCoeff>& poly) {
    if (poly.is_zero()) return SymbolicExpr::number(0);

    std::vector<std::shared_ptr<SymbolicExpr>> terms;
    for (int i = poly.degree(); i >= 0; --i) {
        const auto& coeff = poly.coeffs[i];
        if (coeff == SymbolicPolyCoeff(0)) continue;

        auto coeff_expr = coeff.val;
        if (i == 0) {
            terms.push_back(coeff_expr);
        } else {
            auto var_expr = SymbolicExpr::variable(poly.variable_name);
            std::shared_ptr<SymbolicExpr> var_part;
            if (i == 1) {
                var_part = var_expr;
            } else {
                var_part = SymbolicExpr::power(var_expr, SymbolicExpr::number(i));
            }

            // Check if coefficient is 1 or -1 to simplify output
            if (coeff_expr->is_one()) {
                terms.push_back(var_part);
            } else {
                terms.push_back(SymbolicExpr::multiply(coeff_expr, var_part));
            }
        }
    }

    if (terms.empty()) return SymbolicExpr::number(0);
    if (terms.size() == 1) return terms[0];

    auto result = terms[0];
    for (size_t i = 1; i < terms.size(); ++i) {
        result = SymbolicExpr::add(result, terms[i]);
    }
    return result;
}

// Generate RootOf expression list for an irreducible polynomial
std::vector<std::shared_ptr<SymbolicExpr>> make_rootof_solutions(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var) {

    int n = poly.degree();
    if (n <= 0) return {};

    // Convert the polynomial to a symbolic expression for embedding in RootOf
    auto poly_expr = symbolic_poly_to_expr(poly);

    std::vector<std::shared_ptr<SymbolicExpr>> solutions;
    solutions.reserve(n);

    for (int k = 0; k < n; ++k) {
        solutions.push_back(SymbolicExpr::root_of(poly_expr, var, k));
    }

    return solutions;
}

// ============================================================================
// rootof_evaluate: Numerically evaluate a RootOf(poly, var, k) expression
// Returns the k-th root's numeric value for real roots, nullopt for complex
// roots or out-of-range indices.
//
// Index ordering convention:
//   - Real roots sorted ascending: index 0 = smallest real root
//   - Complex roots follow real roots in lexicographic (Re, Im) order
//   - Conjugate pairs at adjacent indices (positive Im first)
//
// For polynomials with parametric coefficients that block numerical ordering,
// indices 0..n-1 are assigned without ordering guarantees.
// ============================================================================

// Helper: check if a SymbolicExpr contains any free variables
static bool expr_has_free_variables(const std::shared_ptr<SymbolicNode>& node) {
    if (!node) return false;

    struct VarCheckVisitor : public SymbolicVisitor {
        bool found = false;
        void visit(NumberNode&) override {}
        void visit(VariableNode&) override { found = true; }
        void visit(AddNode& n) override { for (auto& op : n.operands) { if (found) return; op->accept(*this); } }
        void visit(MultiplyNode& n) override { for (auto& op : n.operands) { if (found) return; op->accept(*this); } }
        void visit(PowerNode& n) override { n.base->accept(*this); if (!found) n.exponent->accept(*this); }
        void visit(FunctionNode& n) override { for (auto& arg : n.arguments) { if (found) return; arg->accept(*this); } }
        void visit(MatrixNode& n) override {
            if (std::holds_alternative<MatrixNode::DenseStorage>(n.storage)) {
                for (auto& item : std::get<MatrixNode::DenseStorage>(n.storage)) {
                    if (item) { item->accept(*this); if (found) return; }
                }
            } else {
                for (auto& [idx, item] : std::get<MatrixNode::SparseStorage>(n.storage)) {
                    item->accept(*this); if (found) return;
                }
            }
        }
    } visitor;
    node->accept(visitor);
    return visitor.found;
}

// Helper: try to convert a SymbolicPolyCoeff polynomial to Rational polynomial
// Returns true if all coefficients are numeric (no free variables besides the poly variable)
static bool try_convert_to_rational_poly(
    const Polynomial<SymbolicPolyCoeff>& sym_poly,
    Polynomial<Rational>& out_poly)
{
    out_poly = Polynomial<Rational>(sym_poly.variable_name);
    out_poly.coeffs.resize(sym_poly.coeffs.size(), Rational(0));

    for (size_t i = 0; i < sym_poly.coeffs.size(); ++i) {
        auto coeff_expr = sym_poly.coeffs[i].val;
        if (!coeff_expr) {
            out_poly.coeffs[i] = Rational(0);
            continue;
        }

        auto simplified = coeff_expr->simplify();

        // Check if the coefficient has free variables (parametric)
        if (expr_has_free_variables(simplified->root)) {
            return false;
        }

        // Try to extract a Rational value
        if (auto num_node = std::dynamic_pointer_cast<NumberNode>(simplified->root)) {
            if (std::holds_alternative<Rational>(num_node->value)) {
                out_poly.coeffs[i] = std::get<Rational>(num_node->value);
            } else if (std::holds_alternative<BigInt>(num_node->value)) {
                out_poly.coeffs[i] = Rational(std::get<BigInt>(num_node->value));
            } else {
                // lmmc_real_t (double) - convert to Rational
                out_poly.coeffs[i] = Rational::from_double(std::get<lmmc_real_t>(num_node->value));
            }
        } else {
            // Try numeric evaluation
            lmmc_real_t val = simplified->to_numeric();
            if (std::isnan(val) || std::isinf(val)) {
                return false;
            }
            out_poly.coeffs[i] = Rational::from_double(val);
        }
    }

    out_poly.trim();
    return true;
}

std::optional<lmmc_real_t> rootof_evaluate(
    const std::shared_ptr<SymbolicExpr>& rootof_expr)
{
    if (!rootof_expr || !rootof_expr->root) return std::nullopt;

    // Extract the FunctionNode with type RootOf
    auto func_node = std::dynamic_pointer_cast<FunctionNode>(rootof_expr->root);
    if (!func_node || func_node->type != FunctionNode::FuncType::RootOf) {
        return std::nullopt;
    }

    // RootOf has 3 arguments: [poly_expr, var_node, index_node]
    if (func_node->arguments.size() != 3) {
        return std::nullopt;
    }

    auto poly_node = func_node->arguments[0];
    auto var_node = func_node->arguments[1];
    auto index_node = func_node->arguments[2];

    // Extract variable name
    auto var_sym = std::dynamic_pointer_cast<VariableNode>(var_node);
    if (!var_sym) return std::nullopt;
    std::string var = var_sym->name;

    // Extract index k
    auto idx_num = std::dynamic_pointer_cast<NumberNode>(index_node);
    if (!idx_num) return std::nullopt;

    int k = 0;
    if (std::holds_alternative<BigInt>(idx_num->value)) {
        k = std::get<BigInt>(idx_num->value).to_int();
    } else if (std::holds_alternative<Rational>(idx_num->value)) {
        k = static_cast<int>(std::get<Rational>(idx_num->value).to_double());
    } else {
        k = static_cast<int>(std::get<lmmc_real_t>(idx_num->value));
    }

    // Convert the polynomial expression to a Polynomial<SymbolicPolyCoeff>
    auto poly_expr = std::make_shared<SymbolicExpr>(poly_node);
    auto sym_poly = symbolic_to_poly<SymbolicPolyCoeff>(poly_expr, var);

    int degree = sym_poly.degree();
    if (degree <= 0) return std::nullopt;

    // Check index range: must be in [0, degree-1]
    if (k < 0 || k >= degree) {
        return std::nullopt;
    }

    // Try to convert to a Rational polynomial for numerical evaluation
    Polynomial<Rational> rat_poly;
    if (!try_convert_to_rational_poly(sym_poly, rat_poly)) {
        // Parametric coefficients: cannot numerically order roots.
        // We cannot evaluate numerically, return nullopt.
        return std::nullopt;
    }

    // Use Sturm sequence to isolate real roots
    auto intervals = isolate_real_roots(rat_poly);

    // The number of real roots found
    int num_real_roots = static_cast<int>(intervals.size());

    // Index convention:
    //   indices [0, num_real_roots-1] → real roots sorted ascending
    //   indices [num_real_roots, degree-1] → complex roots in lexicographic (Re, Im) order
    //
    // Since this function returns lmmc_real_t, we can only evaluate real roots.
    // For complex root indices, return nullopt.
    if (k >= num_real_roots) {
        return std::nullopt;
    }

    // Refine the k-th real root interval using Newton-Raphson
    const auto& interval = intervals[k]; // Already sorted ascending by isolate_real_roots

    // Compute midpoint as starting point
    lmmc_real_t x0 = (interval.first.to_double() + interval.second.to_double()) / 2.0;

    // Build the symbolic polynomial expression and its derivative for Newton-Raphson
    auto f_expr = poly_expr->simplify();
    auto df_expr = f_expr->differentiate(var);
    if (!df_expr) {
        // Fallback: just return midpoint if we can't differentiate
        return x0;
    }
    df_expr = df_expr->simplify();

    // Use Newton-Raphson to refine
    SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = newton_raphson(f_expr, df_expr, var, x0, opts);
    if (result.has_value()) {
        return result->value;
    }

    // If Newton-Raphson didn't converge, try bisection as fallback
    lmmc_real_t lo = interval.first.to_double();
    lmmc_real_t hi = interval.second.to_double();

    // Simple bisection fallback
    for (int iter = 0; iter < 200; ++iter) {
        lmmc_real_t mid = (lo + hi) / 2.0;
        auto fmid_expr = f_expr->substitute(var, SymbolicExpr::number(mid));
        lmmc_real_t fmid = fmid_expr->to_numeric();

        if (std::abs(fmid) < 1e-12) {
            return mid;
        }

        auto flo_expr = f_expr->substitute(var, SymbolicExpr::number(lo));
        lmmc_real_t flo = flo_expr->to_numeric();

        if ((flo > 0) != (fmid > 0)) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    // Return best approximation (midpoint of final interval)
    return (lo + hi) / 2.0;
}

// ============================================================================
// rootof_simplify: Simplify a RootOf(poly, var, k) expression
//
// 1. If poly.degree() <= 4, solve using closed-form formulas, sort roots by
//    the index ordering convention (real ascending, complex in lexicographic
//    (Re, Im) order), and return the k-th root.
// 2. If poly factors (stretch goal), recurse into each factor with
//    redistributed indices.
// 3. Otherwise, return the original RootOf node unchanged.
// ============================================================================

// Helper: try to evaluate a symbolic expression to a complex number
// Returns (value, success)
static std::pair<std::complex<double>, bool> try_eval_complex(
    const std::shared_ptr<SymbolicExpr>& expr)
{
    if (!expr || !expr->root) return {{0.0, 0.0}, false};

    // Try direct numeric evaluation
    lmmc_real_t val = expr->to_numeric();
    if (!std::isnan(val) && !std::isinf(val)) {
        return {{val, 0.0}, true};
    }

    // Check if expression contains 'i' (imaginary unit) - try to decompose
    // For now, if to_numeric fails, we check if it's a purely imaginary or complex expression
    // by looking for patterns like a + b*i
    auto simplified = expr->simplify();
    val = simplified->to_numeric();
    if (!std::isnan(val) && !std::isinf(val)) {
        return {{val, 0.0}, true};
    }

    // Cannot evaluate numerically (may contain symbolic parameters or complex parts)
    return {{0.0, 0.0}, false};
}

// Helper struct for sorting roots by the ordering convention
struct RootWithIndex {
    std::shared_ptr<SymbolicExpr> expr;
    double real_part;
    double imag_part;
    bool is_real;
    bool eval_success;
};

// Helper: solve a polynomial of degree 1-4 using closed-form formulas
static std::vector<std::shared_ptr<SymbolicExpr>> solve_closed_form_from_poly(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var)
{
    int deg = poly.degree();
    if (deg <= 0) return {};

    // Extract coefficients as SymbolicExpr
    auto get_coeff = [&](int d) -> std::shared_ptr<SymbolicExpr> {
        if (d < 0 || d > deg) return SymbolicExpr::number(0);
        return poly.coeffs[d].val ? poly.coeffs[d].val : SymbolicExpr::number(0);
    };

    if (deg == 1) {
        // ax + b = 0 => x = -b/a
        auto a = get_coeff(1);
        auto b = get_coeff(0);
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        return { SymbolicExpr::divide(neg_b, a)->simplify() };
    } else if (deg == 2) {
        auto a = get_coeff(2);
        auto b = get_coeff(1);
        auto c = get_coeff(0);

        // Use quadratic formula: (-b ± sqrt(b²-4ac)) / (2a)
        auto b2 = SymbolicExpr::power(b, SymbolicExpr::number(2));
        auto four_ac = SymbolicExpr::multiply(SymbolicExpr::number(4), SymbolicExpr::multiply(a, c));
        auto delta = SymbolicExpr::add(b2, SymbolicExpr::multiply(four_ac, SymbolicExpr::number(-1)));
        auto sqrt_delta = SymbolicExpr::sqrt(delta);
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        auto two_a = SymbolicExpr::multiply(SymbolicExpr::number(2), a);

        auto x1 = SymbolicExpr::divide(SymbolicExpr::add(neg_b, sqrt_delta), two_a)->simplify();
        auto x2 = SymbolicExpr::divide(SymbolicExpr::add(neg_b, SymbolicExpr::multiply(sqrt_delta, SymbolicExpr::number(-1))), two_a)->simplify();
        return { x1, x2 };
    } else if (deg == 3) {
        auto a = get_coeff(3);
        auto b = get_coeff(2);
        auto c = get_coeff(1);
        auto d = get_coeff(0);
        return solve_cubic(a, b, c, d, var);
    } else if (deg == 4) {
        auto a = get_coeff(4);
        auto b = get_coeff(3);
        auto c = get_coeff(2);
        auto d = get_coeff(1);
        auto e = get_coeff(0);
        return solve_quartic(a, b, c, d, e, var);
    }

    return {};
}

// Sort roots by the ordering convention:
//   - Real roots ascending (index 0 = smallest real root)
//   - Complex roots follow in lexicographic (Re, Im) order
//   - Conjugate pairs at adjacent indices (positive Im first)
static std::vector<std::shared_ptr<SymbolicExpr>> sort_roots_by_convention(
    const std::vector<std::shared_ptr<SymbolicExpr>>& roots)
{
    if (roots.empty()) return roots;

    // Evaluate each root numerically
    std::vector<RootWithIndex> evaluated;
    evaluated.reserve(roots.size());

    bool all_evaluated = true;
    for (const auto& root : roots) {
        RootWithIndex ri;
        ri.expr = root;
        auto [val, success] = try_eval_complex(root);
        ri.eval_success = success;
        if (success) {
            ri.real_part = val.real();
            ri.imag_part = val.imag();
            // Consider a root "real" if imaginary part is negligible
            ri.is_real = (std::abs(ri.imag_part) < 1e-10);
        } else {
            all_evaluated = false;
            ri.real_part = 0.0;
            ri.imag_part = 0.0;
            ri.is_real = true; // Assume real if we can't evaluate
        }
        evaluated.push_back(ri);
    }

    // If we can't evaluate all roots numerically, return them in original order
    if (!all_evaluated) {
        return roots;
    }

    // Separate real and complex roots
    std::vector<RootWithIndex> real_roots;
    std::vector<RootWithIndex> complex_roots;

    for (auto& ri : evaluated) {
        if (ri.is_real) {
            real_roots.push_back(ri);
        } else {
            complex_roots.push_back(ri);
        }
    }

    // Sort real roots ascending
    std::sort(real_roots.begin(), real_roots.end(),
        [](const RootWithIndex& a, const RootWithIndex& b) {
            return a.real_part < b.real_part;
        });

    // Sort complex roots by lexicographic (Re, Im) order
    // Conjugate pairs should be adjacent with positive Im first
    std::sort(complex_roots.begin(), complex_roots.end(),
        [](const RootWithIndex& a, const RootWithIndex& b) {
            if (std::abs(a.real_part - b.real_part) > 1e-10) {
                return a.real_part < b.real_part;
            }
            // Same real part: positive imaginary first (descending Im)
            return a.imag_part > b.imag_part;
        });

    // Combine: real roots first, then complex roots
    std::vector<std::shared_ptr<SymbolicExpr>> sorted;
    sorted.reserve(roots.size());
    for (const auto& ri : real_roots) {
        sorted.push_back(ri.expr);
    }
    for (const auto& ri : complex_roots) {
        sorted.push_back(ri.expr);
    }

    return sorted;
}

std::shared_ptr<SymbolicExpr> rootof_simplify(
    const std::shared_ptr<SymbolicExpr>& rootof_expr)
{
    if (!rootof_expr || !rootof_expr->root) return rootof_expr;

    // Extract the FunctionNode with type RootOf
    auto func_node = std::dynamic_pointer_cast<FunctionNode>(rootof_expr->root);
    if (!func_node || func_node->type != FunctionNode::FuncType::RootOf) {
        return rootof_expr;
    }

    // RootOf has 3 arguments: [poly_expr, var_node, index_node]
    if (func_node->arguments.size() != 3) {
        return rootof_expr;
    }

    auto poly_node = func_node->arguments[0];
    auto var_node = func_node->arguments[1];
    auto index_node = func_node->arguments[2];

    // Extract variable name
    auto var_sym = std::dynamic_pointer_cast<VariableNode>(var_node);
    if (!var_sym) return rootof_expr;
    std::string var = var_sym->name;

    // Extract index k
    auto idx_num = std::dynamic_pointer_cast<NumberNode>(index_node);
    if (!idx_num) return rootof_expr;

    int k = 0;
    if (std::holds_alternative<BigInt>(idx_num->value)) {
        k = std::get<BigInt>(idx_num->value).to_int();
    } else if (std::holds_alternative<Rational>(idx_num->value)) {
        k = static_cast<int>(std::get<Rational>(idx_num->value).to_double());
    } else {
        k = static_cast<int>(std::get<lmmc_real_t>(idx_num->value));
    }

    // Convert the polynomial expression to a Polynomial<SymbolicPolyCoeff>
    auto poly_expr = std::make_shared<SymbolicExpr>(poly_node);
    auto sym_poly = symbolic_to_poly<SymbolicPolyCoeff>(poly_expr, var);

    int degree = sym_poly.degree();
    if (degree <= 0) return rootof_expr;

    // Check index range
    if (k < 0 || k >= degree) {
        return rootof_expr;
    }

    // Case 1: degree <= 4 - use closed-form formulas
    if (degree <= 4) {
        auto roots = solve_closed_form_from_poly(sym_poly, var);
        if (roots.empty()) {
            return rootof_expr;
        }

        // Sort roots by the ordering convention
        auto sorted_roots = sort_roots_by_convention(roots);

        // Return the k-th root (bounds check)
        if (k >= 0 && k < static_cast<int>(sorted_roots.size())) {
            return sorted_roots[k];
        }
        // If we got fewer roots than expected, return original
        return rootof_expr;
    }

    // Case 2: degree > 4 - factoring (stretch goal, not implemented)
    // For now, return the original RootOf node unchanged
    return rootof_expr;
}

} // namespace lamina
