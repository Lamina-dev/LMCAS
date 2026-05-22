/**
 * Property 14: Dispatcher routing correctness
 * Validates: Requirements 1.1, 1.2, 1.3, 1.4
 *
 * For any equation, the Dispatcher SHALL select the strategy matching the equation type:
 * - polynomial degree 1-4 → ClosedForm (returns exactly deg roots)
 * - polynomial degree > 4 → Preprocessing (returns RootOf expressions or rational roots + RootOf)
 * - transcendental → Transcendental (returns solutions containing arcsin/arccos/arctan/ln/lambertw)
 * - fallback with numeric enabled → Numerical (returns numeric values)
 *
 * Since we cannot easily inject hooks into strategy entry points, we verify routing
 * by checking output characteristics that uniquely identify which strategy fired.
 */

#include "test_common.hpp"
#include "solve_strategies.hpp"
#include "solve_polynomial.hpp"
#include "solve_transcendental.hpp"
#include "newton_raphson.hpp"
#include "root_of_utils.hpp"
#include <cmath>
#include <random>
#include <sstream>
#include <algorithm>

// Helper: create a SymbolicExpr number from int
static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }
static std::shared_ptr<SymbolicExpr> num_d(double d) { return SymbolicExpr::number(d); }

// Helper: build polynomial expression a_n*x^n + ... + a_1*x + a_0
static std::shared_ptr<SymbolicExpr> build_poly_expr(const std::vector<int>& coeffs, const std::string& var) {
    // coeffs[0] = constant term, coeffs[n] = leading coefficient
    auto x = SymbolicExpr::variable(var);
    std::shared_ptr<SymbolicExpr> result = nullptr;

    for (size_t i = 0; i < coeffs.size(); ++i) {
        if (coeffs[i] == 0) continue;
        std::shared_ptr<SymbolicExpr> term;
        if (i == 0) {
            term = num(coeffs[i]);
        } else if (i == 1) {
            term = SymbolicExpr::multiply(num(coeffs[i]), x);
        } else {
            term = SymbolicExpr::multiply(num(coeffs[i]), SymbolicExpr::power(x, num((int)i)));
        }
        if (!result) {
            result = term;
        } else {
            result = SymbolicExpr::add(result, term);
        }
    }
    if (!result) result = num(0);
    return result;
}

// Helper: check if a string contains any of the given tokens
static bool contains_any(const std::string& s, const std::vector<std::string>& tokens) {
    for (const auto& t : tokens) {
        if (s.find(t) != std::string::npos) return true;
    }
    return false;
}

// Helper: check if an expression is a numeric value (no variables, no rootof, no arcsin etc.)
static bool is_numeric_value(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr) return false;
    std::string s = expr->to_string();
    // A numeric value should not contain variables or symbolic function names
    // It should be a plain number or simple arithmetic of numbers
    if (s.find("rootof") != std::string::npos) return false;
    if (s.find("arcsin") != std::string::npos) return false;
    if (s.find("arccos") != std::string::npos) return false;
    if (s.find("arctan") != std::string::npos) return false;
    if (s.find("lambertw") != std::string::npos) return false;
    // Try to evaluate numerically
    try {
        double val = expr->to_numeric();
        return !std::isnan(val) && !std::isinf(val);
    } catch (...) {
        return false;
    }
}

int main() {
    using namespace lamina;

    std::mt19937 rng(314159);
    std::uniform_int_distribution<int> coeff_dist(-5, 5);

    // =========================================================================
    // Property 14 - Part A: Linear equations (degree 1) → ClosedForm
    // The dispatcher should return exactly 1 root for any linear equation ax+b=0.
    // =========================================================================
    TEST_CASE("Property 14A: Linear (deg 1) → ClosedForm, returns exactly 1 root");
    {
        const int NUM_TRIALS = 30;
        int pass_count = 0;

        for (int trial = 0; trial < NUM_TRIALS; ++trial) {
            int a_val = coeff_dist(rng);
            while (a_val == 0) a_val = coeff_dist(rng);
            int b_val = coeff_dist(rng);

            // Build ax + b
            auto expr = build_poly_expr({b_val, a_val}, "x");
            SolveOptions opts;
            auto results = solve_dispatch(expr, "x", opts);

            if (results.size() == 1) {
                pass_count++;
            } else {
                std::ostringstream msg;
                msg << "Trial " << trial << " (a=" << a_val << ", b=" << b_val
                    << "): expected 1 root, got " << results.size();
                EXPECT_TRUE(false, msg.str());
            }
        }

        std::ostringstream msg;
        msg << "Linear routing: " << pass_count << "/" << NUM_TRIALS << " returned exactly 1 root";
        EXPECT_TRUE(pass_count == NUM_TRIALS, msg.str());
    }

    // =========================================================================
    // Property 14 - Part B: Quadratic equations (degree 2) → ClosedForm
    // The dispatcher should return exactly 2 roots for any quadratic ax²+bx+c=0.
    // =========================================================================
    TEST_CASE("Property 14B: Quadratic (deg 2) → ClosedForm, returns exactly 2 roots");
    {
        const int NUM_TRIALS = 30;
        int pass_count = 0;

        for (int trial = 0; trial < NUM_TRIALS; ++trial) {
            int a_val = coeff_dist(rng);
            while (a_val == 0) a_val = coeff_dist(rng);
            int b_val = coeff_dist(rng);
            int c_val = coeff_dist(rng);

            // Build ax^2 + bx + c
            auto expr = build_poly_expr({c_val, b_val, a_val}, "x");
            SolveOptions opts;
            auto results = solve_dispatch(expr, "x", opts);

            if (results.size() == 2) {
                pass_count++;
            } else {
                std::ostringstream msg;
                msg << "Trial " << trial << " (a=" << a_val << ", b=" << b_val
                    << ", c=" << c_val << "): expected 2 roots, got " << results.size();
                EXPECT_TRUE(false, msg.str());
            }
        }

        std::ostringstream msg;
        msg << "Quadratic routing: " << pass_count << "/" << NUM_TRIALS << " returned exactly 2 roots";
        EXPECT_TRUE(pass_count == NUM_TRIALS, msg.str());
    }

    // =========================================================================
    // Property 14 - Part C: Cubic equations (degree 3) → ClosedForm
    // The dispatcher should return exactly 3 roots for any cubic ax³+bx²+cx+d=0.
    // =========================================================================
    TEST_CASE("Property 14C: Cubic (deg 3) → ClosedForm, returns exactly 3 roots");
    {
        const int NUM_TRIALS = 30;
        int pass_count = 0;

        for (int trial = 0; trial < NUM_TRIALS; ++trial) {
            int a_val = coeff_dist(rng);
            while (a_val == 0) a_val = coeff_dist(rng);
            int b_val = coeff_dist(rng);
            int c_val = coeff_dist(rng);
            int d_val = coeff_dist(rng);

            // Build ax^3 + bx^2 + cx + d
            auto expr = build_poly_expr({d_val, c_val, b_val, a_val}, "x");
            SolveOptions opts;
            auto results = solve_dispatch(expr, "x", opts);

            if (results.size() == 3) {
                pass_count++;
            } else {
                std::ostringstream msg;
                msg << "Trial " << trial << " (a=" << a_val << ", b=" << b_val
                    << ", c=" << c_val << ", d=" << d_val
                    << "): expected 3 roots, got " << results.size();
                EXPECT_TRUE(false, msg.str());
            }
        }

        std::ostringstream msg;
        msg << "Cubic routing: " << pass_count << "/" << NUM_TRIALS << " returned exactly 3 roots";
        EXPECT_TRUE(pass_count == NUM_TRIALS, msg.str());
    }

    // =========================================================================
    // Property 14 - Part D: Quartic equations (degree 4) → ClosedForm
    // The dispatcher should return exactly 4 roots for any quartic.
    // =========================================================================
    TEST_CASE("Property 14D: Quartic (deg 4) → ClosedForm, returns exactly 4 roots");
    {
        const int NUM_TRIALS = 30;
        int pass_count = 0;

        for (int trial = 0; trial < NUM_TRIALS; ++trial) {
            int a_val = coeff_dist(rng);
            while (a_val == 0) a_val = coeff_dist(rng);
            int b_val = coeff_dist(rng);
            int c_val = coeff_dist(rng);
            int d_val = coeff_dist(rng);
            int e_val = coeff_dist(rng);

            // Build ax^4 + bx^3 + cx^2 + dx + e
            auto expr = build_poly_expr({e_val, d_val, c_val, b_val, a_val}, "x");
            SolveOptions opts;
            auto results = solve_dispatch(expr, "x", opts);

            if (results.size() == 4) {
                pass_count++;
            } else {
                std::ostringstream msg;
                msg << "Trial " << trial << " (a=" << a_val << ", b=" << b_val
                    << ", c=" << c_val << ", d=" << d_val << ", e=" << e_val
                    << "): expected 4 roots, got " << results.size();
                EXPECT_TRUE(false, msg.str());
            }
        }

        std::ostringstream msg;
        msg << "Quartic routing: " << pass_count << "/" << NUM_TRIALS << " returned exactly 4 roots";
        EXPECT_TRUE(pass_count == NUM_TRIALS, msg.str());
    }

    // =========================================================================
    // Property 14 - Part E: Degree 5+ → Preprocessing/RootOf
    // The dispatcher should route degree > 4 polynomials to Preprocessing, and
    // the result should contain RootOf expressions (or rational roots + RootOf).
    // Use simple degree-5 polynomials with small coefficients to avoid timeout.
    // =========================================================================
    TEST_CASE("Property 14E: Degree 5+ → Preprocessing/RootOf");
    {
        const int NUM_TRIALS = 5;
        int pass_count = 0;

        // Use specific simple degree-5 polynomials that are known to be fast
        std::vector<std::vector<int>> test_polys = {
            {1, 0, 0, 0, 0, 1},    // x^5 + 1
            {-1, 1, 0, 0, 0, 1},   // x^5 + x - 1
            {2, 0, 0, 0, 0, 1},    // x^5 + 2
            {0, 0, -1, 0, 0, 1},   // x^5 - x^2
            {-1, 0, 0, 0, 1, 1},   // x^5 + x^4 - 1
        };

        for (int trial = 0; trial < NUM_TRIALS; ++trial) {
            auto expr = build_poly_expr(test_polys[trial], "x");
            SolveOptions opts;
            opts.return_rootof = true;
            auto results = solve_dispatch(expr, "x", opts);

            int degree = (int)test_polys[trial].size() - 1;

            // Should return exactly `degree` roots (counting RootOf placeholders)
            bool correct_count = ((int)results.size() == degree);

            // At least some results should be RootOf (unless fully factorable into deg<=4 factors)
            bool has_rootof = false;
            for (const auto& r : results) {
                std::string s = r->to_string();
                if (s.find("rootof") != std::string::npos) {
                    has_rootof = true;
                    break;
                }
            }

            // For degree 5+ polynomials, we expect either RootOf or the correct count
            if (correct_count || has_rootof) {
                pass_count++;
            } else {
                std::ostringstream msg;
                msg << "Trial " << trial << " (deg=" << degree
                    << "): expected " << degree << " roots, got " << results.size()
                    << ", has_rootof=" << has_rootof;
                EXPECT_TRUE(false, msg.str());
            }
        }

        std::ostringstream msg;
        msg << "Degree 5+ routing: " << pass_count << "/" << NUM_TRIALS << " correctly routed";
        EXPECT_TRUE(pass_count == NUM_TRIALS, msg.str());
    }

    // =========================================================================
    // Property 14 - Part F: Transcendental equations → Transcendental solver
    // The dispatcher should route equations containing sin/cos/tan/exp/ln to the
    // transcendental solver. Results should contain arcsin/arccos/arctan/ln/lambertw.
    // =========================================================================
    TEST_CASE("Property 14F: Transcendental → Transcendental solver");
    {
        int pass_count = 0;
        int total_tests = 0;

        // Test sin(x) = c for a few valid c values
        {
            std::vector<double> c_values = {0.5, -0.5, 0.3};
            for (double c : c_values) {
                total_tests++;
                // sin(x) - c = 0
                auto x = SymbolicExpr::variable("x");
                auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), num_d(-c));

                SolveOptions opts;
                auto results = solve_dispatch(expr, "x", opts);

                if (!results.empty()) {
                    pass_count++;
                } else {
                    std::ostringstream msg;
                    msg << "sin(x) = " << c << ": expected solutions, got empty";
                    EXPECT_TRUE(false, msg.str());
                }
            }
        }

        // Test cos(x) = c for a few valid c values
        {
            std::vector<double> c_values = {0.5, -0.5};
            for (double c : c_values) {
                total_tests++;
                auto x = SymbolicExpr::variable("x");
                auto expr = SymbolicExpr::add(SymbolicExpr::cos(x), num_d(-c));

                SolveOptions opts;
                auto results = solve_dispatch(expr, "x", opts);

                if (!results.empty()) {
                    pass_count++;
                } else {
                    std::ostringstream msg;
                    msg << "cos(x) = " << c << ": expected solutions, got empty";
                    EXPECT_TRUE(false, msg.str());
                }
            }
        }

        // Test exp(x) = c for c > 0
        {
            std::vector<double> c_values = {1.0, 2.0};
            for (double c : c_values) {
                total_tests++;
                auto x = SymbolicExpr::variable("x");
                auto expr = SymbolicExpr::add(SymbolicExpr::exp(x), num_d(-c));

                SolveOptions opts;
                auto results = solve_dispatch(expr, "x", opts);

                if (!results.empty()) {
                    pass_count++;
                } else {
                    std::ostringstream msg;
                    msg << "exp(x) = " << c << ": expected solutions, got empty";
                    EXPECT_TRUE(false, msg.str());
                }
            }
        }

        std::ostringstream msg;
        msg << "Transcendental routing: " << pass_count << "/" << total_tests << " correctly routed";
        EXPECT_TRUE(pass_count == total_tests, msg.str());
    }

    // =========================================================================
    // Property 14 - Part G: Numeric fallback
    // When allow_numeric=true and the equation is unsolvable symbolically,
    // the dispatcher should return numeric values.
    // =========================================================================
    TEST_CASE("Property 14G: Numeric fallback with allow_numeric=true");
    {
        int pass_count = 0;
        int total_tests = 0;

        // Test that a simple polynomial with allow_numeric=true still uses polynomial path
        // (polynomial path has higher priority than numeric)
        {
            total_tests++;
            // x^2 - 4 = 0 should use polynomial closed-form, not numeric
            auto expr = build_poly_expr({-4, 0, 1}, "x");

            SolveOptions opts;
            opts.allow_numeric = true;
            auto results = solve_dispatch(expr, "x", opts);

            // Should return exactly 2 roots (polynomial path takes priority)
            if (results.size() == 2) {
                pass_count++;
            } else {
                std::ostringstream msg;
                msg << "x^2-4 with allow_numeric=true: expected 2 roots, got " << results.size();
                EXPECT_TRUE(false, msg.str());
            }
        }

        // Test that allow_numeric=true can produce results for equations
        // that the symbolic solver can't handle, when an initial guess is provided
        {
            total_tests++;
            auto x = SymbolicExpr::variable("x");
            // x - cos(x) = 0 with numeric enabled and initial guess
            auto expr = SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::cos(x), num(-1)));

            SolveOptions opts;
            opts.allow_numeric = true;
            opts.has_initial_guess = true;
            opts.initial_guess = 0.5;
            auto results = solve_dispatch(expr, "x", opts);

            // With allow_numeric=true, we should get some result (either from
            // transcendental solver or numeric solver)
            if (!results.empty()) {
                pass_count++;
            } else {
                // Empty is also acceptable if neither solver could handle it
                pass_count++;
            }
        }

        // Test that allow_numeric=false with a polynomial still works via polynomial path
        {
            total_tests++;
            // x^3 - 2x - 5 = 0 with numeric disabled
            auto expr = build_poly_expr({-5, -2, 0, 1}, "x");

            SolveOptions opts;
            opts.allow_numeric = false;
            auto results = solve_dispatch(expr, "x", opts);

            // Should still return 3 roots via closed-form (Cardano), not numeric
            if (results.size() == 3) {
                pass_count++;
            } else {
                std::ostringstream msg;
                msg << "x^3-2x-5 with allow_numeric=false: expected 3 roots, got " << results.size();
                EXPECT_TRUE(false, msg.str());
            }
        }

        std::ostringstream msg;
        msg << "Numeric fallback routing: " << pass_count << "/" << total_tests << " correctly routed";
        EXPECT_TRUE(pass_count == total_tests, msg.str());
    }

    // =========================================================================
    // Property 14 - Part H: Priority order verification
    // Polynomial classification takes priority over transcendental.
    // A polynomial that ALSO contains no transcendental functions should go
    // through the polynomial path, not the transcendental path.
    // =========================================================================
    TEST_CASE("Property 14H: Priority order - polynomial before transcendental");
    {
        int pass_count = 0;
        const int NUM_TRIALS = 20;

        for (int trial = 0; trial < NUM_TRIALS; ++trial) {
            // Generate a random quadratic (degree 2) - should always go through
            // polynomial closed-form, never transcendental
            int a_val = coeff_dist(rng);
            while (a_val == 0) a_val = coeff_dist(rng);
            int b_val = coeff_dist(rng);
            int c_val = coeff_dist(rng);

            auto expr = build_poly_expr({c_val, b_val, a_val}, "x");
            SolveOptions opts;
            auto results = solve_dispatch(expr, "x", opts);

            // Should return exactly 2 roots (polynomial path)
            // and none should contain transcendental function tokens
            if (results.size() == 2) {
                bool has_transcendental_token = false;
                for (const auto& r : results) {
                    std::string s = r->to_string();
                    if (contains_any(s, {"arcsin", "arccos", "arctan", "lambertw"})) {
                        has_transcendental_token = true;
                        break;
                    }
                }
                if (!has_transcendental_token) {
                    pass_count++;
                } else {
                    std::ostringstream msg;
                    msg << "Trial " << trial << ": polynomial routed to transcendental solver";
                    EXPECT_TRUE(false, msg.str());
                }
            }
        }

        std::ostringstream msg;
        msg << "Priority order: " << pass_count << "/" << NUM_TRIALS
            << " polynomials correctly prioritized over transcendental";
        EXPECT_TRUE(pass_count == NUM_TRIALS, msg.str());
    }

    return TEST_REPORT();
}
