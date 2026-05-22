#include "../include/inequality_solver.hpp"
#include "../include/poly_utils.hpp"
#include "../include/solve_polynomial.hpp"
#include "../include/solve_strategies.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <limits>
#include <functional>

namespace lamina {

// ============================================================================
// Internal helpers
// ============================================================================

// Determine the sign of the leading coefficient of a polynomial expression.
// Returns +1 or -1. If the sign cannot be determined, defaults to +1.
static int determine_leading_sign(const Polynomial<SymbolicPolyCoeff>& poly) {
    if (poly.is_zero()) return 0;
    auto lc = poly.lead_coeff().val;
    if (!lc) return 1;
    auto simplified = lc->simplify();
    if (!simplified) return 1;

    try {
        double val = simplified->to_numeric();
        if (val > 0) return 1;
        if (val < 0) return -1;
    } catch (...) {}

    // Try checking if it's a NumberNode
    if (auto num = std::dynamic_pointer_cast<NumberNode>(simplified->root)) {
        if (std::holds_alternative<BigInt>(num->value)) {
            return std::get<BigInt>(num->value).IsNegative() ? -1 : 1;
        }
        if (std::holds_alternative<Rational>(num->value)) {
            return std::get<Rational>(num->value).get_numerator().IsNegative() ? -1 : 1;
        }
        if (std::holds_alternative<lmmc_real_t>(num->value)) {
            return std::get<lmmc_real_t>(num->value) < 0 ? -1 : 1;
        }
    }
    return 1;
}

// Find all real roots of a polynomial with their multiplicities.
// Uses square-free factorization to determine multiplicities, then solves each factor.
// Returns pairs of (root_expr, multiplicity).
static std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>> find_roots_with_multiplicity(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& variable) {

    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>> result;

    // Convert to Rational polynomial for square-free factorization
    auto poly_rat = symbolic_to_poly<Rational>(expr, variable);
    if (poly_rat.is_zero() || poly_rat.degree() <= 0) {
        return result;
    }

    // Use square-free factorization to get factors with multiplicities
    auto factors = square_free_factorization(poly_rat);
    if (factors.empty()) {
        // Fallback: try solving the whole polynomial directly
        auto poly_spc = symbolic_to_poly<SymbolicPolyCoeff>(expr, variable);
        if (!poly_spc.is_zero() && poly_spc.degree() >= 1) {
            auto roots = solve_by_factoring(poly_spc, variable);
            for (const auto& root : roots) {
                if (!root) continue;
                try {
                    double val = root->to_numeric();
                    if (std::isfinite(val)) {
                        result.push_back({root, 1});
                    }
                } catch (...) {}
            }
        }
        return result;
    }

    // For each factor, solve it using solve_by_factoring on the SymbolicPolyCoeff version
    // This avoids the double-simplification issue in SymbolicExpr::solve()
    for (const auto& [factor, mult] : factors) {
        if (factor.degree() <= 0) continue;
        
        // For degree 1: solve directly (linear: ax + b = 0 → x = -b/a)
        if (factor.degree() == 1) {
            // coeffs[0] = constant, coeffs[1] = leading
            Rational a = factor.coeffs[1];
            Rational b = factor.coeffs[0];
            if (a != Rational(0)) {
                Rational root_val = Rational(0) - b / a;
                auto root_expr = SymbolicExpr::number(root_val);
                try {
                    double val = root_expr->to_numeric();
                    if (std::isfinite(val)) {
                        result.push_back({root_expr, mult});
                    }
                } catch (...) {}
            }
            continue;
        }
        
        // For degree 2: use quadratic formula directly
        if (factor.degree() == 2) {
            Rational a = factor.coeffs[2];
            Rational b = factor.coeffs[1];
            Rational c = factor.coeffs[0];
            // discriminant = b² - 4ac
            Rational disc = b * b - Rational(4) * a * c;
            double disc_val = disc.to_double();
            if (disc_val < -1e-10) continue;  // No real roots
            if (disc_val < 0) disc_val = 0;  // Snap near-zero discriminant
            
            double a_val = a.to_double();
            double b_val = b.to_double();
            double sqrt_disc = std::sqrt(disc_val);
            
            double r1 = (-b_val + sqrt_disc) / (2.0 * a_val);
            double r2 = (-b_val - sqrt_disc) / (2.0 * a_val);
            
            // Use numeric values directly for the roots
            if (std::isfinite(r1)) {
                result.push_back({SymbolicExpr::number(r1), mult});
            }
            if (std::isfinite(r2) && std::abs(r1 - r2) > 1e-10) {
                result.push_back({SymbolicExpr::number(r2), mult});
            }
            continue;
        }
        
        // For degree >= 3: convert to SymbolicPolyCoeff and use solve_by_factoring
        std::vector<SymbolicPolyCoeff> spc_coeffs;
        for (int i = 0; i <= factor.degree(); ++i) {
            spc_coeffs.push_back(SymbolicPolyCoeff(SymbolicExpr::number(factor.coeffs[i])));
        }
        Polynomial<SymbolicPolyCoeff> factor_spc(spc_coeffs, variable);
        
        auto factor_roots = solve_by_factoring(factor_spc, variable);
        for (const auto& root : factor_roots) {
            if (!root) continue;
            try {
                double val = root->to_numeric();
                if (std::isfinite(val)) {
                    result.push_back({root, mult});
                }
            } catch (...) {}
        }
    }

    return result;
}

// Compare two root expressions numerically for sorting
static bool root_less_than(const std::shared_ptr<SymbolicExpr>& a,
                           const std::shared_ptr<SymbolicExpr>& b) {
    try {
        double va = a->to_numeric();
        double vb = b->to_numeric();
        return va < vb;
    } catch (...) {
        return false;
    }
}

// Check if two root expressions are numerically equal
static bool roots_equal(const std::shared_ptr<SymbolicExpr>& a,
                        const std::shared_ptr<SymbolicExpr>& b) {
    try {
        double va = a->to_numeric();
        double vb = b->to_numeric();
        return std::abs(va - vb) < 1e-10;
    } catch (...) {
        return false;
    }
}

// ============================================================================
// InequalitySolver::build_sign_chart
// ============================================================================

std::vector<SignChartEntry> InequalitySolver::build_sign_chart(
    const std::shared_ptr<SymbolicExpr>& poly,
    const std::string& variable,
    const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
    const std::vector<int>& multiplicities) {

    std::vector<SignChartEntry> chart;

    if (roots.empty()) {
        // No roots: the polynomial has constant sign everywhere
        auto p = symbolic_to_poly<SymbolicPolyCoeff>(poly, variable);
        int sign = determine_leading_sign(p);
        chart.push_back({Interval::entire_line(), sign});
        return chart;
    }

    // Determine leading sign (sign at +∞)
    auto p = symbolic_to_poly<SymbolicPolyCoeff>(poly, variable);
    int leading_sign = determine_leading_sign(p);

    // Build chart from right to left
    // Start with the rightmost interval (last_root, +∞) which has sign = leading_sign
    int current_sign = leading_sign;

    // We'll build intervals from left to right but track sign from right to left
    // First, compute signs for each interval
    size_t n = roots.size();
    std::vector<int> interval_signs(n + 1);

    // Sign in (roots[n-1], +∞) = leading_sign
    interval_signs[n] = leading_sign;

    // Going from right to left, flip sign at odd-multiplicity roots
    for (int i = (int)n - 1; i >= 0; --i) {
        interval_signs[i] = interval_signs[i + 1];
        if (multiplicities[i] % 2 != 0) {
            interval_signs[i] = -interval_signs[i];
        }
    }

    // Now build the chart entries from left to right
    // Interval 0: (-∞, roots[0])
    {
        Interval iv;
        iv.lower = Endpoint::neg_inf();
        iv.upper = Endpoint::open(roots[0]);
        chart.push_back({iv, interval_signs[0]});
    }

    // Intervals between consecutive roots
    for (size_t i = 0; i + 1 < n; ++i) {
        Interval iv;
        iv.lower = Endpoint::open(roots[i]);
        iv.upper = Endpoint::open(roots[i + 1]);
        chart.push_back({iv, interval_signs[i + 1]});
    }

    // Last interval: (roots[n-1], +∞)
    {
        Interval iv;
        iv.lower = Endpoint::open(roots[n - 1]);
        iv.upper = Endpoint::pos_inf();
        chart.push_back({iv, interval_signs[n]});
    }

    return chart;
}

// ============================================================================
// InequalitySolver::select_intervals
// ============================================================================

IntervalUnion InequalitySolver::select_intervals(
    const std::vector<SignChartEntry>& chart,
    InequalityType type,
    const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
    const std::vector<int>& multiplicities) {

    std::vector<Interval> result_intervals;

    // Determine which sign we're looking for
    bool want_positive = (type == InequalityType::GreaterThan || type == InequalityType::GreaterEqual);
    bool is_strict = (type == InequalityType::GreaterThan || type == InequalityType::LessThan);
    int target_sign = want_positive ? 1 : -1;

    for (const auto& entry : chart) {
        if (entry.sign == target_sign) {
            result_intervals.push_back(entry.interval);
        }
    }

    // For non-strict inequalities, include roots where the expression equals zero
    if (!is_strict) {
        // Merge root points into adjacent intervals by closing endpoints
        for (size_t i = 0; i < roots.size(); ++i) {
            // For non-strict, include the root as a closed point
            // We need to close the endpoints of adjacent intervals that touch this root
            bool merged = false;
            for (auto& iv : result_intervals) {
                // Check if this root is the upper bound of this interval
                if (!iv.upper.is_pos_infinity && iv.upper.value) {
                    if (roots_equal(iv.upper.value, roots[i])) {
                        iv.upper.is_open = false;
                        merged = true;
                    }
                }
                // Check if this root is the lower bound of this interval
                if (!iv.lower.is_neg_infinity && iv.lower.value) {
                    if (roots_equal(iv.lower.value, roots[i])) {
                        iv.lower.is_open = false;
                        merged = true;
                    }
                }
            }
            // If the root wasn't adjacent to any selected interval, add it as a point
            if (!merged) {
                result_intervals.push_back(Interval::point(roots[i]));
            }
        }
    }

    return IntervalUnion(result_intervals);
}

// ============================================================================
// InequalitySolver::solve_inequality
// ============================================================================

IntervalUnion InequalitySolver::solve_inequality(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable) {

    if (!expr) return IntervalUnion::empty();

    // Check if expression is polynomial in the variable
    auto poly = symbolic_to_poly<SymbolicPolyCoeff>(expr, variable);
    if (poly.is_zero()) {
        // Check if the expression actually depends on the variable
        // If it does, it's a non-polynomial expression (like sin(x)) → return empty
        if (depends_on_var(expr->root, variable)) {
            return IntervalUnion::empty();
        }
        // Truly zero polynomial: 0 > 0 → empty, 0 >= 0 → entire line
        if (type == InequalityType::GreaterEqual || type == InequalityType::LessEqual) {
            return IntervalUnion::entire_line();
        }
        return IntervalUnion::empty();
    }

    if (poly.degree() <= 0) {
        // Constant polynomial
        auto lc = poly.lead_coeff().val;
        if (!lc) return IntervalUnion::empty();
        try {
            double val = lc->simplify()->to_numeric();
            bool satisfies = false;
            switch (type) {
                case InequalityType::GreaterThan: satisfies = (val > 0); break;
                case InequalityType::GreaterEqual: satisfies = (val >= 0); break;
                case InequalityType::LessThan: satisfies = (val < 0); break;
                case InequalityType::LessEqual: satisfies = (val <= 0); break;
            }
            return satisfies ? IntervalUnion::entire_line() : IntervalUnion::empty();
        } catch (...) {
            return IntervalUnion::empty();
        }
    }

    // Verify the expression is truly polynomial by checking with Rational coefficients
    // If the expression contains non-polynomial functions of the variable, return empty
    {
        auto poly_rat = symbolic_to_poly<Rational>(expr, variable);
        if (poly_rat.is_zero() && depends_on_var(expr->root, variable)) {
            // The expression depends on the variable but couldn't be converted to polynomial
            return IntervalUnion::empty();
        }
        // If poly_rat conversion succeeded with degree >= 1, the expression is polynomial.
        // Skip the expensive reconstruction + simplify check which can hang on certain expressions.
    }

    // Find roots with multiplicities
    auto roots_with_mult = find_roots_with_multiplicity(expr, variable);

    // Sort roots numerically
    std::sort(roots_with_mult.begin(), roots_with_mult.end(),
        [](const auto& a, const auto& b) {
            return root_less_than(a.first, b.first);
        });

    // Remove duplicate roots (keep highest multiplicity)
    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>> unique_roots;
    for (const auto& [root, mult] : roots_with_mult) {
        if (!unique_roots.empty() && roots_equal(unique_roots.back().first, root)) {
            // Keep the higher multiplicity
            unique_roots.back().second = std::max(unique_roots.back().second, mult);
        } else {
            unique_roots.push_back({root, mult});
        }
    }

    // Separate into roots and multiplicities vectors
    std::vector<std::shared_ptr<SymbolicExpr>> roots;
    std::vector<int> multiplicities;
    for (const auto& [root, mult] : unique_roots) {
        roots.push_back(root);
        multiplicities.push_back(mult);
    }

    // Build sign chart
    auto chart = build_sign_chart(expr, variable, roots, multiplicities);

    // Select intervals based on inequality type
    return select_intervals(chart, type, roots, multiplicities);
}

// ============================================================================
// InequalitySolver::solve_rational_inequality
// ============================================================================

IntervalUnion InequalitySolver::solve_rational_inequality(
    const std::shared_ptr<SymbolicExpr>& numerator,
    const std::shared_ptr<SymbolicExpr>& denominator,
    InequalityType type,
    const std::string& variable) {

    if (!numerator || !denominator) return IntervalUnion::empty();

    // Check that both numerator and denominator are polynomial in the variable
    auto num_poly = symbolic_to_poly<SymbolicPolyCoeff>(numerator, variable);
    auto den_poly = symbolic_to_poly<SymbolicPolyCoeff>(denominator, variable);

    // If denominator is zero polynomial, undefined - return empty
    if (den_poly.is_zero()) return IntervalUnion::empty();

    // If numerator is zero polynomial: 0/q(x) = 0
    // 0 > 0 → empty, 0 >= 0 → entire line minus denominator roots, 0 < 0 → empty, 0 <= 0 → entire line minus den roots
    if (num_poly.is_zero()) {
        if (type == InequalityType::GreaterThan || type == InequalityType::LessThan) {
            return IntervalUnion::empty();
        }
        // Non-strict: 0 >= 0 or 0 <= 0 is true everywhere except where denominator is zero
        auto den_roots_with_mult = find_roots_with_multiplicity(denominator, variable);
        if (den_roots_with_mult.empty()) {
            return IntervalUnion::entire_line();
        }
        // Start with entire line and exclude denominator roots
        // Build intervals excluding denominator roots
        std::vector<std::shared_ptr<SymbolicExpr>> den_roots;
        for (const auto& [root, mult] : den_roots_with_mult) {
            den_roots.push_back(root);
        }
        std::sort(den_roots.begin(), den_roots.end(), root_less_than);

        std::vector<Interval> intervals;
        // (-∞, den_roots[0])
        {
            Interval iv;
            iv.lower = Endpoint::neg_inf();
            iv.upper = Endpoint::open(den_roots[0]);
            intervals.push_back(iv);
        }
        for (size_t i = 0; i + 1 < den_roots.size(); ++i) {
            Interval iv;
            iv.lower = Endpoint::open(den_roots[i]);
            iv.upper = Endpoint::open(den_roots[i + 1]);
            intervals.push_back(iv);
        }
        // (den_roots.back(), +∞)
        {
            Interval iv;
            iv.lower = Endpoint::open(den_roots.back());
            iv.upper = Endpoint::pos_inf();
            intervals.push_back(iv);
        }
        return IntervalUnion(intervals);
    }

    // Verify both are actually polynomial (not transcendental)
    // Check numerator using Rational polynomial conversion
    {
        auto num_poly_rat = symbolic_to_poly<Rational>(numerator, variable);
        if (num_poly_rat.is_zero() && depends_on_var(numerator->root, variable)) {
            return IntervalUnion::empty();  // Numerator not polynomial
        }
        // If conversion succeeded with degree >= 1, it's polynomial. Skip expensive simplify check.
    }
    // Check denominator using Rational polynomial conversion
    {
        auto den_poly_rat = symbolic_to_poly<Rational>(denominator, variable);
        if (den_poly_rat.is_zero() && depends_on_var(denominator->root, variable)) {
            return IntervalUnion::empty();  // Denominator not polynomial
        }
        // If conversion succeeded with degree >= 1, it's polynomial. Skip expensive simplify check.
    }

    // Step 1: Find all real roots of numerator with multiplicities
    auto num_roots_with_mult = find_roots_with_multiplicity(numerator, variable);

    // Step 2: Find all real roots of denominator with multiplicities
    auto den_roots_with_mult = find_roots_with_multiplicity(denominator, variable);

    // Step 3: Merge all critical points and sort
    // We need to track which roots come from numerator vs denominator
    struct CriticalPoint {
        std::shared_ptr<SymbolicExpr> value;
        int num_multiplicity;  // multiplicity in numerator (0 if not a num root)
        int den_multiplicity;  // multiplicity in denominator (0 if not a den root)
    };

    std::vector<CriticalPoint> critical_points;

    for (const auto& [root, mult] : num_roots_with_mult) {
        critical_points.push_back({root, mult, 0});
    }
    for (const auto& [root, mult] : den_roots_with_mult) {
        // Check if this root is already in the list (shared root)
        bool found = false;
        for (auto& cp : critical_points) {
            if (roots_equal(cp.value, root)) {
                cp.den_multiplicity = mult;
                found = true;
                break;
            }
        }
        if (!found) {
            critical_points.push_back({root, 0, mult});
        }
    }

    // Sort critical points numerically
    std::sort(critical_points.begin(), critical_points.end(),
        [](const CriticalPoint& a, const CriticalPoint& b) {
            return root_less_than(a.value, b.value);
        });

    // Remove duplicates (shouldn't happen after merge, but be safe)
    // Already handled above

    // Step 4: Build combined sign chart
    // The sign of p(x)/q(x) = sign(p(x)) * sign(q(x))
    // We determine the sign at +∞ from leading coefficients, then flip at odd-multiplicity roots

    int num_leading_sign = determine_leading_sign(num_poly);
    int den_leading_sign = determine_leading_sign(den_poly);
    int combined_leading_sign = num_leading_sign * den_leading_sign;

    // Compute combined multiplicity for sign-change purposes:
    // At each critical point, the combined expression changes sign if the total
    // odd-multiplicity contribution causes a flip.
    // For the rational function p/q, the sign flips at a critical point if
    // the sum of multiplicities (from both num and den) at that point is odd.
    size_t n_cp = critical_points.size();
    std::vector<int> interval_signs(n_cp + 1);

    // Sign in the rightmost interval (last_cp, +∞) = combined_leading_sign
    interval_signs[n_cp] = combined_leading_sign;

    // Going from right to left, flip sign at odd total multiplicity
    for (int i = (int)n_cp - 1; i >= 0; --i) {
        int total_mult = critical_points[i].num_multiplicity + critical_points[i].den_multiplicity;
        interval_signs[i] = interval_signs[i + 1];
        if (total_mult % 2 != 0) {
            interval_signs[i] = -interval_signs[i];
        }
    }

    // Step 5: Select intervals where the sign satisfies the inequality type
    bool want_positive = (type == InequalityType::GreaterThan || type == InequalityType::GreaterEqual);
    bool is_strict = (type == InequalityType::GreaterThan || type == InequalityType::LessThan);
    int target_sign = want_positive ? 1 : -1;

    std::vector<Interval> result_intervals;

    if (n_cp == 0) {
        // No critical points: constant sign everywhere
        if (interval_signs[0] == target_sign) {
            result_intervals.push_back(Interval::entire_line());
        }
    } else {
        // First interval: (-∞, cp[0])
        if (interval_signs[0] == target_sign) {
            Interval iv;
            iv.lower = Endpoint::neg_inf();
            iv.upper = Endpoint::open(critical_points[0].value);
            result_intervals.push_back(iv);
        }

        // Middle intervals: (cp[i], cp[i+1])
        for (size_t i = 0; i + 1 < n_cp; ++i) {
            if (interval_signs[i + 1] == target_sign) {
                Interval iv;
                iv.lower = Endpoint::open(critical_points[i].value);
                iv.upper = Endpoint::open(critical_points[i + 1].value);
                result_intervals.push_back(iv);
            }
        }

        // Last interval: (cp[n-1], +∞)
        if (interval_signs[n_cp] == target_sign) {
            Interval iv;
            iv.lower = Endpoint::open(critical_points[n_cp - 1].value);
            iv.upper = Endpoint::pos_inf();
            result_intervals.push_back(iv);
        }
    }

    // Step 6: Endpoint handling
    // - Denominator roots: ALWAYS excluded (open endpoints) regardless of strict/non-strict
    // - Numerator-only roots (where expression = 0):
    //   - Strict (> or <): excluded (open endpoints)
    //   - Non-strict (≥ or ≤): included (closed endpoints)
    if (!is_strict) {
        for (const auto& cp : critical_points) {
            // Only include numerator-only roots (where den_multiplicity == 0)
            // At these points, the expression equals zero, which satisfies >= or <=
            if (cp.num_multiplicity > 0 && cp.den_multiplicity == 0) {
                // Try to close the endpoint in adjacent intervals
                bool merged = false;
                for (auto& iv : result_intervals) {
                    // Check if this root is the upper bound of this interval
                    if (!iv.upper.is_pos_infinity && iv.upper.value) {
                        if (roots_equal(iv.upper.value, cp.value)) {
                            iv.upper.is_open = false;
                            merged = true;
                        }
                    }
                    // Check if this root is the lower bound of this interval
                    if (!iv.lower.is_neg_infinity && iv.lower.value) {
                        if (roots_equal(iv.lower.value, cp.value)) {
                            iv.lower.is_open = false;
                            merged = true;
                        }
                    }
                }
                // If the root wasn't adjacent to any selected interval, add it as a point
                if (!merged) {
                    result_intervals.push_back(Interval::point(cp.value));
                }
            }
        }
    }

    return IntervalUnion(result_intervals);
}

// ============================================================================
// InequalitySolver::solve_inequalities
// ============================================================================

IntervalUnion InequalitySolver::solve_inequalities(
    const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                 InequalityType>>& inequalities,
    const std::string& variable) {

    if (inequalities.empty()) return IntervalUnion::entire_line();

    // Solve each inequality individually
    IntervalUnion result = solve_inequality(inequalities[0].first, inequalities[0].second, variable);

    // Intersect all results
    for (size_t i = 1; i < inequalities.size(); ++i) {
        auto solution = solve_inequality(inequalities[i].first, inequalities[i].second, variable);
        result = result.intersect(solution);
        if (result.is_empty()) break;  // Early exit
    }

    return result;
}

// ============================================================================
// InequalitySolver::build_parametric_solution
// ============================================================================

IntervalUnion InequalitySolver::build_parametric_solution(
    const std::vector<std::shared_ptr<SymbolicExpr>>& symbolic_roots,
    const std::vector<int>& multiplicities,
    int leading_sign,
    InequalityType type) {

    if (symbolic_roots.empty()) {
        // No roots: constant sign everywhere = leading_sign
        bool want_positive = (type == InequalityType::GreaterThan || type == InequalityType::GreaterEqual);
        int target_sign = want_positive ? 1 : -1;
        if (leading_sign == target_sign) {
            return IntervalUnion::entire_line();
        }
        return IntervalUnion::empty();
    }

    // Build sign chart from right to left using leading_sign and multiplicities
    size_t n = symbolic_roots.size();
    std::vector<int> interval_signs(n + 1);

    // Sign in (roots[n-1], +∞) = leading_sign
    interval_signs[n] = leading_sign;

    // Going from right to left, flip sign at odd-multiplicity roots
    for (int i = (int)n - 1; i >= 0; --i) {
        interval_signs[i] = interval_signs[i + 1];
        if (multiplicities[i] % 2 != 0) {
            interval_signs[i] = -interval_signs[i];
        }
    }

    // Select intervals based on inequality type
    bool want_positive = (type == InequalityType::GreaterThan || type == InequalityType::GreaterEqual);
    bool is_strict = (type == InequalityType::GreaterThan || type == InequalityType::LessThan);
    int target_sign = want_positive ? 1 : -1;

    std::vector<Interval> result_intervals;

    // First interval: (-∞, roots[0])
    if (interval_signs[0] == target_sign) {
        Interval iv;
        iv.lower = Endpoint::neg_inf();
        iv.upper = Endpoint::open(symbolic_roots[0]);
        result_intervals.push_back(iv);
    }

    // Middle intervals: (roots[i], roots[i+1])
    for (size_t i = 0; i + 1 < n; ++i) {
        if (interval_signs[i + 1] == target_sign) {
            Interval iv;
            iv.lower = Endpoint::open(symbolic_roots[i]);
            iv.upper = Endpoint::open(symbolic_roots[i + 1]);
            result_intervals.push_back(iv);
        }
    }

    // Last interval: (roots[n-1], +∞)
    if (interval_signs[n] == target_sign) {
        Interval iv;
        iv.lower = Endpoint::open(symbolic_roots[n - 1]);
        iv.upper = Endpoint::pos_inf();
        result_intervals.push_back(iv);
    }

    // For non-strict inequalities, include roots (close endpoints)
    if (!is_strict) {
        for (size_t i = 0; i < symbolic_roots.size(); ++i) {
            bool merged = false;
            for (auto& iv : result_intervals) {
                // Check if this root is the upper bound of this interval
                if (!iv.upper.is_pos_infinity && iv.upper.value) {
                    // Compare symbolically
                    auto diff = SymbolicExpr::add(iv.upper.value,
                        SymbolicExpr::multiply(symbolic_roots[i], SymbolicExpr::number(-1)));
                    if (diff->simplify()->is_zero()) {
                        iv.upper.is_open = false;
                        merged = true;
                    }
                }
                // Check if this root is the lower bound of this interval
                if (!iv.lower.is_neg_infinity && iv.lower.value) {
                    auto diff = SymbolicExpr::add(iv.lower.value,
                        SymbolicExpr::multiply(symbolic_roots[i], SymbolicExpr::number(-1)));
                    if (diff->simplify()->is_zero()) {
                        iv.lower.is_open = false;
                        merged = true;
                    }
                }
            }
            // If the root wasn't adjacent to any selected interval, add it as a point
            if (!merged) {
                result_intervals.push_back(Interval::point(symbolic_roots[i]));
            }
        }
    }

    return IntervalUnion(result_intervals);
}

// ============================================================================
// InequalitySolver::solve_parametric_inequality
// ============================================================================

// Helper: check if a symbolic expression depends on any of the given parameters
static bool depends_on_any_param(const std::shared_ptr<SymbolicExpr>& expr,
                                  const std::vector<std::string>& parameters) {
    if (!expr || !expr->root) return false;
    for (const auto& param : parameters) {
        if (depends_on_var(expr->root, param)) return true;
    }
    return false;
}

// Helper: solve a polynomial with symbolic coefficients, returning symbolic roots
// Uses the existing solve_by_factoring / closed-form solvers which already handle
// symbolic coefficients and return symbolic expressions.
static std::vector<std::shared_ptr<SymbolicExpr>> solve_symbolic_poly(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& variable) {

    if (poly.is_zero() || poly.degree() < 1) return {};

    // For degree <= 4, use closed-form solvers which handle symbolic coefficients
    int deg = poly.degree();
    auto get_coeff = [&](int d) -> std::shared_ptr<SymbolicExpr> {
        if (d < 0 || d > deg) return SymbolicExpr::number(0);
        return poly.coeffs[d].val ? poly.coeffs[d].val : SymbolicExpr::number(0);
    };

    if (deg == 1) {
        // Linear: ax + b = 0 → x = -b/a
        auto a = get_coeff(1);
        auto b = get_coeff(0);
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        auto root = SymbolicExpr::divide(neg_b, a)->simplify();
        return { root };
    }

    // For degree >= 2, use solve_by_factoring which handles symbolic coefficients
    auto results = solve_by_factoring(poly, variable);
    return results;
}

PiecewiseIntervalResult InequalitySolver::solve_parametric_inequality(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable,
    const std::vector<std::string>& parameters) {

    PiecewiseIntervalResult result;

    if (!expr) return result;

    // If no parameters specified, fall back to the non-parametric solver
    if (parameters.empty()) {
        auto solution = solve_inequality(expr, type, variable);
        PiecewiseIntervalResult::Case single_case;
        single_case.condition = nullptr;  // unconditional
        single_case.solution = solution;
        result.cases.push_back(single_case);
        return result;
    }

    // Convert expression to polynomial in the variable (coefficients may be symbolic)
    auto poly = symbolic_to_poly<SymbolicPolyCoeff>(expr, variable);

    if (poly.is_zero()) {
        // Check if the expression depends on the variable but isn't polynomial
        if (depends_on_var(expr->root, variable)) {
            // Non-polynomial expression → return empty
            return result;
        }
        // Truly zero: 0 > 0 → empty, 0 >= 0 → entire line
        PiecewiseIntervalResult::Case zero_case;
        zero_case.condition = nullptr;
        if (type == InequalityType::GreaterEqual || type == InequalityType::LessEqual) {
            zero_case.solution = IntervalUnion::entire_line();
        } else {
            zero_case.solution = IntervalUnion::empty();
        }
        result.cases.push_back(zero_case);
        return result;
    }

    // Get the leading coefficient
    int deg = poly.degree();
    auto leading_coeff = poly.coeffs[deg].val;
    if (!leading_coeff) leading_coeff = SymbolicExpr::number(0);
    leading_coeff = leading_coeff->simplify();

    // Check if the leading coefficient depends on parameters
    bool lc_depends_on_params = depends_on_any_param(leading_coeff, parameters);

    if (!lc_depends_on_params) {
        // Leading coefficient is a constant (doesn't depend on parameters)
        // Determine its sign
        int leading_sign = 1;
        try {
            double val = leading_coeff->to_numeric();
            leading_sign = (val > 0) ? 1 : -1;
        } catch (...) {
            // Try NumberNode check
            if (auto num = std::dynamic_pointer_cast<NumberNode>(leading_coeff->root)) {
                if (std::holds_alternative<BigInt>(num->value)) {
                    leading_sign = std::get<BigInt>(num->value).IsNegative() ? -1 : 1;
                } else if (std::holds_alternative<Rational>(num->value)) {
                    leading_sign = std::get<Rational>(num->value).get_numerator().IsNegative() ? -1 : 1;
                } else if (std::holds_alternative<lmmc_real_t>(num->value)) {
                    leading_sign = std::get<lmmc_real_t>(num->value) < 0 ? -1 : 1;
                }
            }
        }

        // Solve the polynomial symbolically
        auto symbolic_roots = solve_symbolic_poly(poly, variable);

        // All roots have multiplicity 1 from solve_by_factoring (it returns repeated roots)
        // For simplicity, assign multiplicity 1 to each root
        std::vector<int> multiplicities(symbolic_roots.size(), 1);

        // Build the parametric solution
        auto solution = build_parametric_solution(symbolic_roots, multiplicities, leading_sign, type);

        PiecewiseIntervalResult::Case single_case;
        single_case.condition = nullptr;  // unconditional
        single_case.solution = solution;
        result.cases.push_back(single_case);
    } else {
        // Leading coefficient depends on parameters → piecewise solution
        // Case 1: leading_coeff > 0
        {
            PiecewiseIntervalResult::Case pos_case;
            pos_case.condition = std::make_shared<SymbolicExpr>(
                std::make_shared<RelationalNode>(
                    leading_coeff->root,
                    SymbolicExpr::number(0)->root,
                    RelationalNode::Op::GT));

            // Solve with leading_sign = +1
            auto symbolic_roots = solve_symbolic_poly(poly, variable);
            std::vector<int> multiplicities(symbolic_roots.size(), 1);
            pos_case.solution = build_parametric_solution(symbolic_roots, multiplicities, 1, type);
            result.cases.push_back(pos_case);
        }

        // Case 2: leading_coeff < 0
        {
            PiecewiseIntervalResult::Case neg_case;
            neg_case.condition = std::make_shared<SymbolicExpr>(
                std::make_shared<RelationalNode>(
                    leading_coeff->root,
                    SymbolicExpr::number(0)->root,
                    RelationalNode::Op::LT));

            // Solve with leading_sign = -1
            auto symbolic_roots = solve_symbolic_poly(poly, variable);
            std::vector<int> multiplicities(symbolic_roots.size(), 1);
            neg_case.solution = build_parametric_solution(symbolic_roots, multiplicities, -1, type);
            result.cases.push_back(neg_case);
        }

        // Case 3: Degenerate case - leading_coeff == 0 (degree drops)
        {
            PiecewiseIntervalResult::Case degen_case;
            degen_case.condition = std::make_shared<SymbolicExpr>(
                std::make_shared<RelationalNode>(
                    leading_coeff->root,
                    SymbolicExpr::number(0)->root,
                    RelationalNode::Op::EQ));

            // When leading coefficient is zero, the polynomial has lower degree
            // Build the reduced polynomial (remove leading term)
            if (deg >= 1) {
                std::vector<SymbolicPolyCoeff> reduced_coeffs;
                for (int i = 0; i < deg; ++i) {
                    reduced_coeffs.push_back(poly.coeffs[i]);
                }
                Polynomial<SymbolicPolyCoeff> reduced_poly(reduced_coeffs, variable);

                if (reduced_poly.is_zero()) {
                    // All coefficients are zero → handle as zero polynomial
                    if (type == InequalityType::GreaterEqual || type == InequalityType::LessEqual) {
                        degen_case.solution = IntervalUnion::entire_line();
                    } else {
                        degen_case.solution = IntervalUnion::empty();
                    }
                } else {
                    // Recursively solve the reduced polynomial inequality
                    // Check if the new leading coefficient also depends on parameters
                    auto new_lc = reduced_poly.lead_coeff().val;
                    if (new_lc) new_lc = new_lc->simplify();

                    if (new_lc && depends_on_any_param(new_lc, parameters)) {
                        // Still parametric - recursively call parametric solver
                        // Convert reduced poly back to expression
                        auto reduced_expr = SymbolicExpr::number(0);
                        auto var_expr = SymbolicExpr::variable(variable);
                        for (int i = reduced_poly.degree(); i >= 0; --i) {
                            auto coeff_val = reduced_poly.coeffs[i].val;
                            if (!coeff_val) continue;
                            if (i == 0) {
                                reduced_expr = SymbolicExpr::add(reduced_expr, coeff_val);
                            } else if (i == 1) {
                                reduced_expr = SymbolicExpr::add(reduced_expr,
                                    SymbolicExpr::multiply(coeff_val, var_expr));
                            } else {
                                reduced_expr = SymbolicExpr::add(reduced_expr,
                                    SymbolicExpr::multiply(coeff_val,
                                        SymbolicExpr::power(var_expr, SymbolicExpr::number(i))));
                            }
                        }
                        reduced_expr = reduced_expr->simplify();
                        auto sub_result = solve_parametric_inequality(reduced_expr, type, variable, parameters);
                        // Take the first case from the sub-result as the degenerate solution
                        if (!sub_result.cases.empty()) {
                            degen_case.solution = sub_result.cases[0].solution;
                        } else {
                            degen_case.solution = IntervalUnion::empty();
                        }
                    } else {
                        // Reduced polynomial has constant leading coefficient
                        int reduced_leading_sign = 1;
                        if (new_lc) {
                            try {
                                double val = new_lc->to_numeric();
                                reduced_leading_sign = (val > 0) ? 1 : -1;
                            } catch (...) {
                                if (auto num = std::dynamic_pointer_cast<NumberNode>(new_lc->root)) {
                                    if (std::holds_alternative<BigInt>(num->value)) {
                                        reduced_leading_sign = std::get<BigInt>(num->value).IsNegative() ? -1 : 1;
                                    } else if (std::holds_alternative<Rational>(num->value)) {
                                        reduced_leading_sign = std::get<Rational>(num->value).get_numerator().IsNegative() ? -1 : 1;
                                    } else if (std::holds_alternative<lmmc_real_t>(num->value)) {
                                        reduced_leading_sign = std::get<lmmc_real_t>(num->value) < 0 ? -1 : 1;
                                    }
                                }
                            }
                        }

                        auto symbolic_roots = solve_symbolic_poly(reduced_poly, variable);
                        std::vector<int> multiplicities(symbolic_roots.size(), 1);
                        degen_case.solution = build_parametric_solution(
                            symbolic_roots, multiplicities, reduced_leading_sign, type);
                    }
                }
            } else {
                degen_case.solution = IntervalUnion::empty();
            }

            result.cases.push_back(degen_case);
        }
    }

    return result;
}

} // namespace lamina
