// newton_raphson.cpp - Newton-Raphson numerical solver with Sturm isolation
// Implements Sturm sequence real root isolation for polynomials with Rational coefficients

#include "newton_raphson.hpp"
#include "poly_utils.hpp"
#include <algorithm>
#include <cmath>

namespace lamina {

// ============================================================================
// Helper: Count sign changes in a sequence of Rational values (skipping zeros)
// ============================================================================
static int count_sign_changes(const std::vector<Rational>& values) {
    int changes = 0;
    int last_sign = 0; // 0 = no sign yet, 1 = positive, -1 = negative

    for (const auto& v : values) {
        if (v.is_zero()) continue;

        int current_sign = (v > Rational(0)) ? 1 : -1;
        if (last_sign != 0 && current_sign != last_sign) {
            ++changes;
        }
        last_sign = current_sign;
    }
    return changes;
}

// ============================================================================
// Helper: Evaluate Sturm sequence at a point and count sign changes
// ============================================================================
static int sturm_sign_changes_at(
    const std::vector<Polynomial<Rational>>& sturm,
    const Rational& x)
{
    std::vector<Rational> values;
    values.reserve(sturm.size());
    for (const auto& p : sturm) {
        values.push_back(p.eval(x));
    }
    return count_sign_changes(values);
}

// ============================================================================
// Helper: Count sign changes at +infinity using leading coefficients' signs
// For Sturm[i] of degree d_i, sign at +inf is sign of leading coefficient
// ============================================================================
static int sturm_sign_changes_at_pos_inf(
    const std::vector<Polynomial<Rational>>& sturm)
{
    std::vector<Rational> signs;
    signs.reserve(sturm.size());
    for (const auto& p : sturm) {
        if (p.is_zero()) continue;
        signs.push_back(p.lead_coeff());
    }
    return count_sign_changes(signs);
}

// ============================================================================
// Helper: Count sign changes at -infinity
// For Sturm[i] of degree d_i, sign at -inf is sign(lead_coeff) * (-1)^d_i
// ============================================================================
static int sturm_sign_changes_at_neg_inf(
    const std::vector<Polynomial<Rational>>& sturm)
{
    std::vector<Rational> signs;
    signs.reserve(sturm.size());
    for (const auto& p : sturm) {
        if (p.is_zero()) continue;
        int deg = p.degree();
        Rational lc = p.lead_coeff();
        // At -inf, sign is lc * (-1)^deg
        if (deg % 2 == 1) {
            lc = -lc;
        }
        signs.push_back(lc);
    }
    return count_sign_changes(signs);
}

// ============================================================================
// Helper: Compute Cauchy bound for real roots of a polynomial
// For polynomial a_n*x^n + ... + a_0, bound = 1 + max(|a_i/a_n|) for i < n
// ============================================================================
static Rational cauchy_bound(const Polynomial<Rational>& poly) {
    if (poly.degree() <= 0) return Rational(1);

    Rational lc = poly.lead_coeff();
    Rational max_ratio(0);

    for (int i = 0; i < poly.degree(); ++i) {
        Rational coeff = poly.coeffs[i];
        if (coeff.is_zero()) continue;
        Rational ratio = coeff.abs() / lc.abs();
        if (ratio > max_ratio) {
            max_ratio = ratio;
        }
    }

    return Rational(1) + max_ratio;
}

// ============================================================================
// Helper: Number of distinct real roots in interval (a, b] using Sturm theorem
// V(a) - V(b) gives the number of distinct real roots in (a, b]
// ============================================================================
static int roots_in_interval(
    const std::vector<Polynomial<Rational>>& sturm,
    const Rational& a,
    const Rational& b)
{
    int va = sturm_sign_changes_at(sturm, a);
    int vb = sturm_sign_changes_at(sturm, b);
    return va - vb;
}

// ============================================================================
// isolate_real_roots: Sturm sequence real root isolation
// Returns intervals [lo, hi] each containing exactly one real root
// ============================================================================
std::vector<std::pair<Rational, Rational>> isolate_real_roots(
    const Polynomial<Rational>& poly)
{
    std::vector<std::pair<Rational, Rational>> result;

    if (poly.is_zero() || poly.degree() <= 0) {
        return result;
    }

    // Work with the square-free part to handle only distinct roots
    Polynomial<Rational> sqfree = poly.square_free_part();
    if (sqfree.is_zero() || sqfree.degree() <= 0) {
        return result;
    }

    // Step 1: Build the Sturm sequence
    // sturm[0] = sqfree
    // sturm[1] = sqfree'
    // sturm[i] = -remainder(sturm[i-2], sturm[i-1])
    std::vector<Polynomial<Rational>> sturm;
    sturm.push_back(sqfree);
    sturm.push_back(sqfree.differentiate());

    while (true) {
        size_t n = sturm.size();
        const auto& prev2 = sturm[n - 2];
        const auto& prev1 = sturm[n - 1];

        if (prev1.is_zero()) break;

        // Compute remainder of prev2 / prev1
        auto [quotient, remainder] = prev2.div_mod(prev1);

        if (remainder.is_zero()) break;

        // Negate the remainder: sturm[i] = -remainder(sturm[i-2], sturm[i-1])
        Polynomial<Rational> neg_rem(remainder.variable_name);
        neg_rem.coeffs.resize(remainder.coeffs.size());
        for (size_t i = 0; i < remainder.coeffs.size(); ++i) {
            neg_rem.coeffs[i] = -remainder.coeffs[i];
        }
        neg_rem.trim();

        sturm.push_back(neg_rem);
    }

    // Step 2: Compute Cauchy bound to get initial interval
    Rational bound = cauchy_bound(sqfree);
    Rational lo = -bound;
    Rational hi = bound;

    // Check total number of real roots
    int total_roots = roots_in_interval(sturm, lo, hi);
    if (total_roots <= 0) {
        return result;
    }

    // Step 3: Bisect intervals until each contains exactly one root
    // Use a work queue of intervals
    struct Interval {
        Rational lo, hi;
        int root_count;
    };

    std::vector<Interval> work_queue;
    work_queue.push_back({lo, hi, total_roots});

    while (!work_queue.empty()) {
        Interval current = work_queue.back();
        work_queue.pop_back();

        if (current.root_count == 0) {
            continue;
        }

        if (current.root_count == 1) {
            result.push_back({current.lo, current.hi});
            continue;
        }

        // Split at midpoint
        Rational mid = (current.lo + current.hi) / Rational(2);

        int left_count = roots_in_interval(sturm, current.lo, mid);
        int right_count = roots_in_interval(sturm, mid, current.hi);

        // Check if mid itself is a root (left_count + right_count < current.root_count)
        // This can happen if mid is exactly a root
        if (left_count + right_count < current.root_count) {
            // mid is a root - add it as a degenerate interval [mid, mid]
            result.push_back({mid, mid});
        }

        if (right_count > 0) {
            work_queue.push_back({mid, current.hi, right_count});
        }
        if (left_count > 0) {
            work_queue.push_back({current.lo, mid, left_count});
        }
    }

    // Sort intervals by their lower bound
    std::sort(result.begin(), result.end(),
        [](const std::pair<Rational, Rational>& a,
           const std::pair<Rational, Rational>& b) {
            return a.first < b.first;
        });

    return result;
}

// ============================================================================
// bisection: Bisection method on a bracket [lo, hi]
// Used as fallback when Newton-Raphson encounters near-zero derivative
// ============================================================================
std::optional<NumericRoot> bisection(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    lmmc_real_t lo,
    lmmc_real_t hi,
    const SolveOptions& opts)
{
    // Evaluate f at endpoints
    lmmc_real_t f_lo = f->substitute(var, SymbolicExpr::number(lo))->to_numeric();
    lmmc_real_t f_hi = f->substitute(var, SymbolicExpr::number(hi))->to_numeric();

    // Check if either endpoint is already a root
    if (std::abs(f_lo) < opts.tolerance) {
        return NumericRoot{lo, std::abs(f_lo), 0};
    }
    if (std::abs(f_hi) < opts.tolerance) {
        return NumericRoot{hi, std::abs(f_hi), 0};
    }

    // Ensure sign change exists; if not, cannot bisect
    if (f_lo * f_hi > 0) {
        return std::nullopt;
    }

    // Bisection iterations (use max_newton_iterations as the iteration cap)
    int max_iter = opts.max_newton_iterations * 3; // bisection converges slower
    for (int i = 1; i <= max_iter; ++i) {
        lmmc_real_t mid = (lo + hi) * 0.5;
        lmmc_real_t f_mid = f->substitute(var, SymbolicExpr::number(mid))->to_numeric();

        if (std::abs(f_mid) < opts.tolerance) {
            return NumericRoot{mid, std::abs(f_mid), i};
        }

        // Also check if interval is tiny enough
        if (std::abs(hi - lo) < opts.tolerance) {
            return NumericRoot{mid, std::abs(f_mid), i};
        }

        if (f_lo * f_mid < 0) {
            hi = mid;
            f_hi = f_mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }

    // Return best estimate even if not fully converged
    lmmc_real_t mid = (lo + hi) * 0.5;
    lmmc_real_t f_mid = f->substitute(var, SymbolicExpr::number(mid))->to_numeric();
    if (std::abs(f_mid) < opts.tolerance * 1000) {
        return NumericRoot{mid, std::abs(f_mid), max_iter};
    }
    return std::nullopt;
}

// ============================================================================
// newton_raphson: Newton-Raphson with bracket for bisection fallback
// When |f'(x)| < 1e-15, falls back to bisection on [bracket_lo, bracket_hi]
// ============================================================================
std::optional<NumericRoot> newton_raphson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    lmmc_real_t bracket_lo,
    lmmc_real_t bracket_hi,
    const SolveOptions& opts)
{
    lmmc_real_t x = x0;
    for (int i = 1; i <= opts.max_newton_iterations; ++i) {
        // Evaluate f(x) and f'(x)
        auto fx_expr = f->substitute(var, SymbolicExpr::number(x));
        auto dfx_expr = df->substitute(var, SymbolicExpr::number(x));

        lmmc_real_t fx = fx_expr->to_numeric();
        lmmc_real_t dfx = dfx_expr->to_numeric();

        if (std::abs(fx) < opts.tolerance) {
            return NumericRoot{x, std::abs(fx), i};
        }

        if (std::abs(dfx) < 1e-15) {
            // Derivative near zero — fall back to bisection on the original bracket
            return bisection(f, var, bracket_lo, bracket_hi, opts);
        }

        lmmc_real_t x_new = x - fx / dfx;

        // Damping: if step is too large after first iteration, halve the step
        if (i > 1 && std::abs(x_new - x) > 2.0 * std::abs(x - x0)) {
            x_new = x - 0.5 * fx / dfx;
        }

        x = x_new;
    }

    return std::nullopt; // Did not converge
}

// ============================================================================
// newton_raphson: Single Newton-Raphson root refinement (no bracket version)
// Without a bracket, derivative-near-zero returns nullopt
// ============================================================================
std::optional<NumericRoot> newton_raphson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    const SolveOptions& opts)
{
    lmmc_real_t x = x0;
    for (int i = 1; i <= opts.max_newton_iterations; ++i) {
        // Evaluate f(x) and f'(x)
        auto fx_expr = f->substitute(var, SymbolicExpr::number(x));
        auto dfx_expr = df->substitute(var, SymbolicExpr::number(x));

        lmmc_real_t fx = fx_expr->to_numeric();
        lmmc_real_t dfx = dfx_expr->to_numeric();

        if (std::abs(fx) < opts.tolerance) {
            return NumericRoot{x, std::abs(fx), i};
        }

        if (std::abs(dfx) < 1e-15) {
            // Derivative near zero, no bracket available — cannot continue
            return std::nullopt;
        }

        lmmc_real_t x_new = x - fx / dfx;

        // Damping: if step is too large after first iteration, halve the step
        if (i > 1 && std::abs(x_new - x) > 2.0 * std::abs(x - x0)) {
            x_new = x - 0.5 * fx / dfx;
        }

        x = x_new;
    }

    return std::nullopt; // Did not converge
}

// ============================================================================
// solve_numeric: Numerical root finding
// Polynomial inputs: isolate real roots via Sturm, run Newton on each interval,
// deflate (x - r) after each success, continue until all intervals are processed
// Non-polynomial inputs: require a user-supplied initial guess x0 (no Sturm),
// return at most one numeric root per call
// ============================================================================
std::vector<NumericRoot> solve_numeric(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts)
{
    std::vector<NumericRoot> results;

    // Try to convert the expression to a polynomial with Rational coefficients
    auto poly = symbolic_to_poly<Rational>(expr, var);

    if (!poly.is_zero() && poly.degree() >= 1) {
        // ====================================================================
        // Polynomial path: Sturm isolation + Newton on each interval + deflation
        // ====================================================================
        Polynomial<Rational> current_poly = poly;

        // Compute the symbolic derivative of the expression for Newton-Raphson
        auto df_expr = expr->differentiate(var);

        while (current_poly.degree() >= 1) {
            // Isolate real roots of the current polynomial
            auto intervals = isolate_real_roots(current_poly);

            if (intervals.empty()) {
                break; // No more real roots to find
            }

            bool found_any = false;
            for (const auto& [lo_rat, hi_rat] : intervals) {
                // Check max_roots limit
                if (opts.max_roots > 0 && (int)results.size() >= opts.max_roots) {
                    return results;
                }

                lmmc_real_t lo = lo_rat.to_double();
                lmmc_real_t hi = hi_rat.to_double();
                lmmc_real_t x0 = (lo + hi) * 0.5;

                // Convert current polynomial to symbolic for Newton evaluation
                auto current_expr = poly_to_symbolic(current_poly);
                auto current_df = current_expr->differentiate(var);

                // Use Newton-Raphson with bracket fallback
                auto root_opt = newton_raphson(current_expr, current_df, var,
                                              x0, lo, hi, opts);

                if (root_opt.has_value()) {
                    NumericRoot root = root_opt.value();

                    // Verify residual against the original expression
                    auto check = expr->substitute(var, SymbolicExpr::number(root.value));
                    lmmc_real_t residual = std::abs(check->to_numeric());
                    root.residual = residual;

                    if (residual < opts.tolerance * 1000) {
                        results.push_back(root);
                        found_any = true;

                        // Deflate: divide current_poly by (x - r)
                        // Construct (x - r) as a polynomial
                        Rational r_rat = Rational::from_double(root.value);
                        Polynomial<Rational> linear_factor({-r_rat, Rational(1)},
                                                          current_poly.variable_name);
                        auto [quotient, remainder] = current_poly.div_mod(linear_factor);

                        // If division isn't exact (due to floating-point root),
                        // use synthetic division with the approximate root
                        if (!remainder.is_zero()) {
                            // Perform synthetic division manually with the numeric root
                            // This is acceptable since we're in numeric mode
                            int deg = current_poly.degree();
                            std::vector<Rational> new_coeffs(deg);
                            // Synthetic division: coefficients from high to low
                            new_coeffs[deg - 1] = current_poly.coeffs[deg];
                            for (int i = deg - 2; i >= 0; --i) {
                                new_coeffs[i] = current_poly.coeffs[i + 1] + new_coeffs[i + 1] * r_rat;
                            }
                            quotient = Polynomial<Rational>(new_coeffs, current_poly.variable_name);
                        }

                        current_poly = quotient;
                        break; // Restart isolation with the deflated polynomial
                    }
                }
            }

            if (!found_any) {
                break; // Could not find any more roots in remaining intervals
            }
        }
    } else {
        // ====================================================================
        // Non-polynomial path: use initial guess, return at most one root
        // ====================================================================
        lmmc_real_t x0 = opts.has_initial_guess ? opts.initial_guess : 0.0;

        // Compute derivative symbolically
        auto df_expr = expr->differentiate(var);

        // Run Newton-Raphson without bracket
        auto root_opt = newton_raphson(expr, df_expr, var, x0, opts);

        if (root_opt.has_value()) {
            NumericRoot root = root_opt.value();
            // Verify residual
            if (root.residual < opts.tolerance) {
                results.push_back(root);
            }
        }
    }

    return results;
}

} // namespace lamina
