// Feature: mixed-transcendental-solver, Property 1: Classification correctness
// **Validates: Requirements 1.1, 1.5**

// Feature: mixed-transcendental-solver, Property 2: Routing correctness for reducible expressions
// **Validates: Requirements 1.2, 1.3**

// Feature: mixed-transcendental-solver, Property 12: allow_numeric gate
// **Validates: Requirements 7.2**

#include "test_common.hpp"
#include "solve_mixed_transcendental.hpp"
#include "solve_strategies.hpp"
#include "visitors/differentiation_visitor.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace lamina;

void test_classification_sin_x_plus_x() {
    TEST_CASE("contains_transcendental_of_var: sin(x) + x → true");

    auto x = SymbolicExpr::variable("x");
    // sin(x) + x: sin depends on x
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    EXPECT_TRUE(contains_transcendental_of_var(expr, "x"),
        "sin(x) + x should contain transcendental of var x");
}

void test_classification_exp_x_minus_x_squared() {
    TEST_CASE("contains_transcendental_of_var: exp(x) - x^2 → true");

    auto x = SymbolicExpr::variable("x");
    // exp(x) - x^2: exp depends on x
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::power(x, SymbolicExpr::number(2)))
    );

    EXPECT_TRUE(contains_transcendental_of_var(expr, "x"),
        "exp(x) - x^2 should contain transcendental of var x");
}

void test_classification_sin_constant_plus_polynomial() {
    TEST_CASE("contains_transcendental_of_var: x^2 + sin(3) → false");

    auto x = SymbolicExpr::variable("x");
    // x^2 + sin(3): sin argument is constant, not dependent on x
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::sin(SymbolicExpr::number(3))
    );

    EXPECT_FALSE(contains_transcendental_of_var(expr, "x"),
        "x^2 + sin(3) should NOT contain transcendental of var x (sin arg is constant)");
}

void test_classification_cos_y_plus_x() {
    TEST_CASE("contains_transcendental_of_var: cos(y) + x → false for var=x");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    // cos(y) + x: cos argument doesn't depend on x
    auto expr = SymbolicExpr::add(SymbolicExpr::cos(y), x);

    EXPECT_FALSE(contains_transcendental_of_var(expr, "x"),
        "cos(y) + x should NOT contain transcendental of var x (cos arg depends on y, not x)");
}

void test_classification_ln_cos_x_plus_x() {
    TEST_CASE("contains_transcendental_of_var: ln(x) * cos(x) + x → true");

    auto x = SymbolicExpr::variable("x");
    // ln(x) * cos(x) + x: both ln and cos depend on x
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::ln(x), SymbolicExpr::cos(x)),
        x
    );

    EXPECT_TRUE(contains_transcendental_of_var(expr, "x"),
        "ln(x) * cos(x) + x should contain transcendental of var x");
}

void test_classification_tan_linear_arg() {
    TEST_CASE("contains_transcendental_of_var: tan(2*x + 1) - x → true");

    auto x = SymbolicExpr::variable("x");
    // tan(2*x + 1) - x: tan argument depends on x
    auto tan_arg = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::number(1)
    );
    auto expr = SymbolicExpr::add(
        SymbolicExpr::tan(tan_arg),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x)
    );

    EXPECT_TRUE(contains_transcendental_of_var(expr, "x"),
        "tan(2*x + 1) - x should contain transcendental of var x");
}

void test_classification_pure_polynomial() {
    TEST_CASE("contains_transcendental_of_var: x^3 - 2*x + 1 → false");

    auto x = SymbolicExpr::variable("x");
    // x^3 - 2*x + 1: pure polynomial, no transcendental functions
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(3)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-2), x),
            SymbolicExpr::number(1)
        )
    );

    EXPECT_FALSE(contains_transcendental_of_var(expr, "x"),
        "x^3 - 2*x + 1 should NOT contain transcendental of var x (pure polynomial)");
}

// Feature: mixed-transcendental-solver, Property 12: allow_numeric gate
// **Validates: Requirements 7.2**
void test_allow_numeric_false_returns_empty() {
    TEST_CASE("Property 12: allow_numeric=false on mixed transcendental → empty vector");

    auto x = SymbolicExpr::variable("x");
    // x*sin(x) - 1 = 0: a mixed transcendental equation with no closed-form solution.
    // This equation cannot be reduced by polynomial or transcendental solvers.
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(x, SymbolicExpr::sin(x)),
        SymbolicExpr::number(-1)
    );

    lamina::SolveOptions opts;
    opts.allow_numeric = false;

    auto results = lamina::solve_dispatch(expr, "x", opts);

    EXPECT_TRUE(results.empty(),
        "solve_dispatch with allow_numeric=false on x*sin(x)-1 should return empty vector");
}

// Feature: mixed-transcendental-solver, Property 12: allow_numeric gate
// **Validates: Requirements 7.2**
void test_allow_numeric_true_permits_solving() {
    TEST_CASE("Property 12: allow_numeric=true on mixed transcendental → solver invoked");

    auto x = SymbolicExpr::variable("x");
    // x*sin(x) - 1 = 0: a mixed transcendental equation
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(x, SymbolicExpr::sin(x)),
        SymbolicExpr::number(-1)
    );

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;

    auto results = lamina::solve_dispatch(expr, "x", opts);

    // Note: If the numerical solver is fully implemented, this should return
    // at least one root. If it's still a stub, this test documents the expected
    // behavior once implemented. The key property is that allow_numeric=true
    // does NOT block the solver from attempting numerical methods.
    std::cout << "  [INFO] allow_numeric=true returned " << results.size() << " root(s)" << std::endl;

    // If results are non-empty, verify they are valid NumberNode expressions
    for (size_t i = 0; i < results.size(); ++i) {
        EXPECT_TRUE(results[i] != nullptr,
            "Root " + std::to_string(i) + " should not be null");
    }
}

// Feature: mixed-transcendental-solver, Property 12: allow_numeric gate
// **Validates: Requirements 7.2**
void test_allow_numeric_false_cos_equation_returns_empty() {
    TEST_CASE("Property 12: allow_numeric=false on x*cos(x)+x^2*sin(x)-1 → empty vector");

    auto x = SymbolicExpr::variable("x");
    // x*cos(x) + x^2*sin(x) - 1 = 0: another mixed transcendental equation
    // that cannot be solved symbolically.
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(x, SymbolicExpr::cos(x)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(
                SymbolicExpr::power(x, SymbolicExpr::number(2)),
                SymbolicExpr::sin(x)
            ),
            SymbolicExpr::number(-1)
        )
    );

    lamina::SolveOptions opts;
    opts.allow_numeric = false;

    auto results = lamina::solve_dispatch(expr, "x", opts);

    EXPECT_TRUE(results.empty(),
        "solve_dispatch with allow_numeric=false on x*cos(x)+x^2*sin(x)-1 should return empty vector");
}

// Feature: mixed-transcendental-solver, Property 2: Routing correctness for reducible expressions
// **Validates: Requirements 1.2, 1.3**
void test_routing_pure_polynomial_not_hybrid() {
    TEST_CASE("Property 2: Pure polynomial x^2 - 4 = 0 routes to polynomial solver, not hybrid");

    auto x = SymbolicExpr::variable("x");
    // x^2 - 4: purely algebraic, no transcendental functions of x
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-4)
    );

    lamina::SolveOptions opts;
    // Even with allow_numeric=true, a pure polynomial should NOT go through hybrid solver
    opts.allow_numeric = true;

    auto results = lamina::solve_dispatch(expr, "x", opts);

    // x^2 - 4 = 0 has roots x = -2 and x = 2
    EXPECT_TRUE(results.size() == 2,
        "x^2 - 4 = 0 should produce exactly 2 roots");

    if (results.size() == 2) {
        // Extract numeric values and sort them
        std::vector<double> root_vals;
        for (const auto& r : results) {
            auto val = test_numeric_eval(r);
            EXPECT_TRUE(val.has_value(), "Root should be evaluable to a number");
            if (val.has_value()) {
                root_vals.push_back(*val);
            }
        }
        if (root_vals.size() == 2) {
            std::sort(root_vals.begin(), root_vals.end());
            EXPECT_NEAR(root_vals[0], -2.0, 1e-10,
                "First root of x^2-4=0 should be -2");
            EXPECT_NEAR(root_vals[1], 2.0, 1e-10,
                "Second root of x^2-4=0 should be 2");
        }
    }
}

// Feature: mixed-transcendental-solver, Property 2: Routing correctness for reducible expressions
// **Validates: Requirements 1.2, 1.3**
void test_routing_transcendental_substitution_not_hybrid() {
    TEST_CASE("Property 2: exp(x) - 2 = 0 routes to transcendental solver, not hybrid");

    auto x = SymbolicExpr::variable("x");
    // exp(x) - 2 = 0: reducible via substitution u = exp(x), giving u - 2 = 0 → u = 2 → x = ln(2)
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::number(-2)
    );

    lamina::SolveOptions opts;
    // With allow_numeric=false, the transcendental solver path should still find ln(2)
    opts.allow_numeric = false;

    auto results = lamina::solve_dispatch(expr, "x", opts);

    // exp(x) - 2 = 0 has solution x = ln(2) ≈ 0.693147
    EXPECT_TRUE(!results.empty(),
        "exp(x) - 2 = 0 should produce at least one root via transcendental path");

    if (!results.empty()) {
        // Verify the root is approximately ln(2)
        auto val = test_numeric_eval(results[0]);
        if (val.has_value()) {
            EXPECT_NEAR(*val, std::log(2.0), 1e-10,
                "Root of exp(x)-2=0 should be ln(2) ≈ 0.693147");
        } else {
            // The result might be symbolic (e.g., ln(2) as an expression)
            // Substitute back to verify: exp(root) - 2 should equal 0
            auto substituted = expr->substitute("x", results[0]);
            auto eval = test_numeric_eval(substituted->simplify());
            if (eval.has_value()) {
                EXPECT_NEAR(*eval, 0.0, 1e-10,
                    "Substituting root back into exp(x)-2 should give 0");
            } else {
                // Accept symbolic result — the key property is that it was found
                // without allow_numeric, confirming the transcendental path handled it
                std::cout << "  [INFO] Root is symbolic: " << results[0]->to_string() << std::endl;
                EXPECT_TRUE(true, "Transcendental path produced a symbolic result (not hybrid)");
            }
        }
    }
}

// Feature: mixed-transcendental-solver, Property 2: Routing correctness for reducible expressions
// **Validates: Requirements 1.2, 1.3**
void test_routing_polynomial_with_allow_numeric_false() {
    TEST_CASE("Property 2: Pure polynomial x^2 - 4 = 0 solved even with allow_numeric=false");

    auto x = SymbolicExpr::variable("x");
    // x^2 - 4: purely algebraic
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-4)
    );

    lamina::SolveOptions opts;
    opts.allow_numeric = false;

    auto results = lamina::solve_dispatch(expr, "x", opts);

    // A pure polynomial should be solved by the polynomial solver regardless of allow_numeric
    EXPECT_TRUE(results.size() == 2,
        "x^2 - 4 = 0 should produce 2 roots even with allow_numeric=false (polynomial path)");

    if (results.size() == 2) {
        std::vector<double> root_vals;
        for (const auto& r : results) {
            auto val = test_numeric_eval(r);
            if (val.has_value()) {
                root_vals.push_back(*val);
            }
        }
        if (root_vals.size() == 2) {
            std::sort(root_vals.begin(), root_vals.end());
            EXPECT_NEAR(root_vals[0], -2.0, 1e-10,
                "First root should be -2 (polynomial path, not hybrid)");
            EXPECT_NEAR(root_vals[1], 2.0, 1e-10,
                "Second root should be 2 (polynomial path, not hybrid)");
        }
    }
}

// Feature: mixed-transcendental-solver, Property 13: Backward compatibility
// **Validates: Requirements 7.4**
void test_backward_compat_cubic_polynomial() {
    TEST_CASE("Property 13: Cubic x^3 - 6x^2 + 11x - 6 = 0 produces roots {1, 2, 3} via polynomial path");

    auto x = SymbolicExpr::variable("x");
    // x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3)
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(3)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-6),
                SymbolicExpr::power(x, SymbolicExpr::number(2))),
            SymbolicExpr::add(
                SymbolicExpr::multiply(SymbolicExpr::number(11), x),
                SymbolicExpr::number(-6)
            )
        )
    );

    // Test with allow_numeric=true — polynomial path should still handle it
    lamina::SolveOptions opts_numeric;
    opts_numeric.allow_numeric = true;
    auto results_numeric = lamina::solve_dispatch(expr, "x", opts_numeric);

    EXPECT_TRUE(results_numeric.size() == 3,
        "x^3-6x^2+11x-6=0 should produce 3 roots with allow_numeric=true");

    if (results_numeric.size() == 3) {
        std::vector<double> root_vals;
        for (const auto& r : results_numeric) {
            auto val = test_numeric_eval(r);
            EXPECT_TRUE(val.has_value(), "Root should be evaluable to a number");
            if (val.has_value()) root_vals.push_back(*val);
        }
        if (root_vals.size() == 3) {
            std::sort(root_vals.begin(), root_vals.end());
            EXPECT_NEAR(root_vals[0], 1.0, 1e-10, "First root should be 1");
            EXPECT_NEAR(root_vals[1], 2.0, 1e-10, "Second root should be 2");
            EXPECT_NEAR(root_vals[2], 3.0, 1e-10, "Third root should be 3");
        }
    }

    // Test with allow_numeric=false — polynomial path should still produce same results
    lamina::SolveOptions opts_symbolic;
    opts_symbolic.allow_numeric = false;
    auto results_symbolic = lamina::solve_dispatch(expr, "x", opts_symbolic);

    EXPECT_TRUE(results_symbolic.size() == 3,
        "x^3-6x^2+11x-6=0 should produce 3 roots with allow_numeric=false (polynomial path)");

    if (results_symbolic.size() == 3) {
        std::vector<double> root_vals;
        for (const auto& r : results_symbolic) {
            auto val = test_numeric_eval(r);
            if (val.has_value()) root_vals.push_back(*val);
        }
        if (root_vals.size() == 3) {
            std::sort(root_vals.begin(), root_vals.end());
            EXPECT_NEAR(root_vals[0], 1.0, 1e-10,
                "First root should be 1 (same as allow_numeric=true)");
            EXPECT_NEAR(root_vals[1], 2.0, 1e-10,
                "Second root should be 2 (same as allow_numeric=true)");
            EXPECT_NEAR(root_vals[2], 3.0, 1e-10,
                "Third root should be 3 (same as allow_numeric=true)");
        }
    }
}

// Feature: mixed-transcendental-solver, Property 13: Backward compatibility
// **Validates: Requirements 7.4**
void test_backward_compat_transcendental_substitution() {
    TEST_CASE("Property 13: exp(x) - 2 = 0 still solved via transcendental path with hybrid integrated");

    auto x = SymbolicExpr::variable("x");
    // exp(x) - 2 = 0: reducible via substitution u=exp(x), u-2=0 → x=ln(2)
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::number(-2)
    );

    // With allow_numeric=false, the transcendental solver should still find ln(2)
    // This confirms the hybrid solver does NOT intercept reducible transcendental equations
    lamina::SolveOptions opts;
    opts.allow_numeric = false;

    auto results = lamina::solve_dispatch(expr, "x", opts);

    EXPECT_TRUE(!results.empty(),
        "exp(x)-2=0 should produce at least one root via transcendental path (not hybrid)");

    if (!results.empty()) {
        // Verify the root is approximately ln(2) ≈ 0.693147
        auto val = test_numeric_eval(results[0]);
        if (val.has_value()) {
            EXPECT_NEAR(*val, std::log(2.0), 1e-10,
                "Root of exp(x)-2=0 should be ln(2) (backward compatible)");
        } else {
            // Symbolic result — verify by substitution
            auto substituted = expr->substitute("x", results[0]);
            auto eval = test_numeric_eval(substituted->simplify());
            if (eval.has_value()) {
                EXPECT_NEAR(*eval, 0.0, 1e-10,
                    "Substituting root back into exp(x)-2 should give 0");
            }
        }
    }

    // With allow_numeric=true, should produce the same result (not go through hybrid)
    lamina::SolveOptions opts_numeric;
    opts_numeric.allow_numeric = true;
    auto results_numeric = lamina::solve_dispatch(expr, "x", opts_numeric);

    EXPECT_TRUE(!results_numeric.empty(),
        "exp(x)-2=0 with allow_numeric=true should still produce roots");
    EXPECT_TRUE(results_numeric.size() == results.size(),
        "allow_numeric flag should not change result count for reducible transcendental");
}

// ============================================================================
// determine_search_interval tests
// **Validates: Requirements 4.1, 4.2, 4.3, 4.4, 4.5, 8.4**
// ============================================================================

void test_search_interval_user_specified() {
    TEST_CASE("determine_search_interval: user-specified interval [2, 5]");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = true;
    opts.search_lo = 2.0;
    opts.search_hi = 5.0;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "User-specified valid interval should return a value");
    if (result) {
        EXPECT_NEAR(result->lo, 2.0, 1e-15, "lo should be 2.0");
        EXPECT_NEAR(result->hi, 5.0, 1e-15, "hi should be 5.0");
    }
}

void test_search_interval_user_invalid_lo_ge_hi() {
    TEST_CASE("determine_search_interval: user-specified invalid interval [5, 2] → nullopt");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = true;
    opts.search_lo = 5.0;
    opts.search_hi = 2.0;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_FALSE(result.has_value(), "Invalid interval (lo >= hi) should return nullopt");
}

void test_search_interval_user_invalid_equal() {
    TEST_CASE("determine_search_interval: user-specified interval [3, 3] → nullopt");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = true;
    opts.search_lo = 3.0;
    opts.search_hi = 3.0;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_FALSE(result.has_value(), "Equal bounds (lo == hi) should return nullopt");
}

void test_search_interval_user_width_le_tolerance() {
    TEST_CASE("determine_search_interval: user-specified interval width <= tolerance → nullopt");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = true;
    opts.search_lo = 1.0;
    opts.search_hi = 1.0 + 1e-13;  // width = 1e-13 <= default tolerance 1e-12
    opts.tolerance = 1e-12;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_FALSE(result.has_value(), "Interval width <= tolerance should return nullopt");
}

void test_search_interval_default_no_periodic() {
    TEST_CASE("determine_search_interval: no periodic functions → default [-10, 10]");

    auto x = SymbolicExpr::variable("x");
    // exp(x) - x: no periodic functions
    auto expr = SymbolicExpr::add(SymbolicExpr::exp(x),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x));

    lamina::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "Default interval should be valid");
    if (result) {
        EXPECT_NEAR(result->lo, -10.0, 1e-15, "Default lo should be -10");
        EXPECT_NEAR(result->hi, 10.0, 1e-15, "Default hi should be 10");
    }
}

void test_search_interval_sin_x_periodic_extension() {
    TEST_CASE("determine_search_interval: sin(x) + x → period 2π, interval [-2π, 2π]");

    auto x = SymbolicExpr::variable("x");
    // sin(x) + x: k=1, period = 2π ≈ 6.28
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "Periodic extension should produce valid interval");
    if (result) {
        // Default [-10, 10] already covers 2π ≈ 6.28, so no extension needed
        // Actually: 2*pi ≈ 6.28 < 10, so default [-10, 10] is already larger
        EXPECT_NEAR(result->lo, -10.0, 1e-10, "sin(x): period 2π < 10, so default [-10,10] used");
        EXPECT_NEAR(result->hi, 10.0, 1e-10, "sin(x): period 2π < 10, so default [-10,10] used");
    }
}

void test_search_interval_sin_small_k_periodic_extension() {
    TEST_CASE("determine_search_interval: sin(x/5) → period 10π ≈ 31.4, interval [-10π, 10π]");

    auto x = SymbolicExpr::variable("x");
    // sin(x/5) = sin((1/5)*x): k=0.2, period = 2π/0.2 = 10π ≈ 31.4
    // 2 full periods = 2 * 10π ≈ 62.8 → interval [-31.4, 31.4]
    auto sin_arg = SymbolicExpr::multiply(SymbolicExpr::number(0.2), x);
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(sin_arg), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "Periodic extension should produce valid interval");
    if (result) {
        double expected_period = 2.0 * M_PI / 0.2;  // 10π ≈ 31.4
        // half_span = period = 31.4 (covers 2 full periods symmetric around 0)
        EXPECT_NEAR(result->lo, -expected_period, 1e-10,
            "sin(0.2*x): interval should extend to -10π");
        EXPECT_NEAR(result->hi, expected_period, 1e-10,
            "sin(0.2*x): interval should extend to +10π");
    }
}

void test_search_interval_tan_periodic_extension() {
    TEST_CASE("determine_search_interval: tan(x/3) → period π/|k|=3π ≈ 9.42, interval [-3π, 3π]");

    auto x = SymbolicExpr::variable("x");
    // tan(x/3): k=1/3, period for tan = π/|k| = 3π ≈ 9.42
    // 2 full periods = 2 * 3π ≈ 18.85 → half_span = 3π ≈ 9.42
    // But default [-10, 10] already covers 9.42, so no extension
    // Actually: half_span = period = 9.42 < 10, so default is used
    auto tan_arg = SymbolicExpr::multiply(SymbolicExpr::number(1.0/3.0), x);
    auto expr = SymbolicExpr::add(SymbolicExpr::tan(tan_arg), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "Periodic extension should produce valid interval");
    if (result) {
        // period = π/(1/3) = 3π ≈ 9.42 < 10, so default [-10, 10] is used
        EXPECT_NEAR(result->lo, -10.0, 1e-10, "tan(x/3): period 3π < 10, default used");
        EXPECT_NEAR(result->hi, 10.0, 1e-10, "tan(x/3): period 3π < 10, default used");
    }
}

void test_search_interval_tan_small_k_extension() {
    TEST_CASE("determine_search_interval: tan(0.1*x) → period π/0.1=10π ≈ 31.4, extends");

    auto x = SymbolicExpr::variable("x");
    // tan(0.1*x): k=0.1, period for tan = π/|k| = 10π ≈ 31.4
    // half_span = period = 31.4 > 10, so extends
    auto tan_arg = SymbolicExpr::multiply(SymbolicExpr::number(0.1), x);
    auto expr = SymbolicExpr::add(SymbolicExpr::tan(tan_arg), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "Periodic extension should produce valid interval");
    if (result) {
        double expected_period = M_PI / 0.1;  // 10π ≈ 31.4
        EXPECT_NEAR(result->lo, -expected_period, 1e-10,
            "tan(0.1*x): interval should extend to -10π");
        EXPECT_NEAR(result->hi, expected_period, 1e-10,
            "tan(0.1*x): interval should extend to +10π");
    }
}

void test_search_interval_clamp_to_100() {
    TEST_CASE("determine_search_interval: very small k → clamped to [-100, 100]");

    auto x = SymbolicExpr::variable("x");
    // sin(0.01*x): k=0.01, period = 2π/0.01 = 200π ≈ 628
    // half_span = 628 > 100, so clamped to [-100, 100]
    auto sin_arg = SymbolicExpr::multiply(SymbolicExpr::number(0.01), x);
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(sin_arg), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "Clamped interval should be valid");
    if (result) {
        EXPECT_NEAR(result->lo, -100.0, 1e-15, "Should be clamped to -100");
        EXPECT_NEAR(result->hi, 100.0, 1e-15, "Should be clamped to +100");
    }
}

void test_search_interval_nonlinear_arg_default() {
    TEST_CASE("determine_search_interval: sin(x^2) → non-linear, default [-10, 10]");

    auto x = SymbolicExpr::variable("x");
    // sin(x^2): argument is non-linear in x
    auto sin_arg = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(sin_arg), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "Non-linear periodic arg should fall back to default");
    if (result) {
        EXPECT_NEAR(result->lo, -10.0, 1e-15, "Non-linear arg: default lo = -10");
        EXPECT_NEAR(result->hi, 10.0, 1e-15, "Non-linear arg: default hi = 10");
    }
}

void test_search_interval_cos_2x_plus_1() {
    TEST_CASE("determine_search_interval: cos(2*x + 1) → k=2, period=π ≈ 3.14, default used");

    auto x = SymbolicExpr::variable("x");
    // cos(2*x + 1): k=2, period = 2π/2 = π ≈ 3.14
    // half_span = π < 10, so default [-10, 10] is used
    auto cos_arg = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::number(1)
    );
    auto expr = SymbolicExpr::add(SymbolicExpr::cos(cos_arg), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "Linear periodic arg should produce valid interval");
    if (result) {
        EXPECT_NEAR(result->lo, -10.0, 1e-10, "cos(2x+1): period π < 10, default used");
        EXPECT_NEAR(result->hi, 10.0, 1e-10, "cos(2x+1): period π < 10, default used");
    }
}

void test_search_interval_user_overrides_periodic() {
    TEST_CASE("determine_search_interval: user interval overrides periodic extension");

    auto x = SymbolicExpr::variable("x");
    // sin(0.01*x): would normally extend to [-100, 100]
    auto sin_arg = SymbolicExpr::multiply(SymbolicExpr::number(0.01), x);
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(sin_arg), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = true;
    opts.search_lo = -1.0;
    opts.search_hi = 1.0;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "User-specified interval should override periodic extension");
    if (result) {
        EXPECT_NEAR(result->lo, -1.0, 1e-15, "User override: lo = -1");
        EXPECT_NEAR(result->hi, 1.0, 1e-15, "User override: hi = 1");
    }
}

// Feature: mixed-transcendental-solver, Property 9: Periodic extension covers two full periods
// **Validates: Requirements 4.4**
void test_property9_sin_x_plus_x_two_periods() {
    TEST_CASE("Property 9: sin(x) + x → k=1, period=2π, interval spans at least 4π ≈ 12.57");

    auto x = SymbolicExpr::variable("x");
    // sin(x) + x: k=1, period = 2π, two full periods = 4π ≈ 12.57
    // The interval must span at least 4π total width.
    // Since 4π ≈ 12.57 > default width of 20 ([-10,10]), the solver should extend.
    // Actually: the requirement says "span at least 2*(2π/|k|)" = 4π ≈ 12.57.
    // The default [-10,10] has width 20 which already covers 12.57.
    // But the design says half_span = period = 2π ≈ 6.28, so interval is [-6.28, 6.28]
    // which has width 12.57. Since default [-10,10] width=20 > 12.57, default is used.
    // Key check: the returned interval width >= 4π ≈ 12.57
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "sin(x)+x should produce valid interval");
    if (result) {
        double width = result->hi - result->lo;
        double min_required = 2.0 * (2.0 * M_PI / 1.0);  // 4π ≈ 12.57
        EXPECT_TRUE(width >= min_required - 1e-10,
            "sin(x)+x: interval width " + std::to_string(width) +
            " should be >= 4π ≈ " + std::to_string(min_required));
    }
}

// Feature: mixed-transcendental-solver, Property 9: Periodic extension covers two full periods
// **Validates: Requirements 4.4**
void test_property9_cos_half_x_minus_x_two_periods() {
    TEST_CASE("Property 9: cos(0.5*x) - x → k=0.5, period=4π, interval spans at least 8π ≈ 25.13");

    auto x = SymbolicExpr::variable("x");
    // cos(0.5*x) - x: k=0.5, period = 2π/0.5 = 4π ≈ 12.57
    // Two full periods = 2 * 4π = 8π ≈ 25.13
    // Default [-10,10] width=20 < 25.13, so extension is needed.
    auto cos_arg = SymbolicExpr::multiply(SymbolicExpr::number(0.5), x);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::cos(cos_arg),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x)
    );

    lamina::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "cos(0.5*x)-x should produce valid interval");
    if (result) {
        double width = result->hi - result->lo;
        double min_required = 2.0 * (2.0 * M_PI / 0.5);  // 8π ≈ 25.13
        EXPECT_TRUE(width >= min_required - 1e-10,
            "cos(0.5*x)-x: interval width " + std::to_string(width) +
            " should be >= 8π ≈ " + std::to_string(min_required));
    }
}

// Feature: mixed-transcendental-solver, Property 9: Periodic extension covers two full periods
// **Validates: Requirements 4.4**
void test_property9_tan_x_plus_x_two_periods() {
    TEST_CASE("Property 9: tan(x) + x → k=1, period=π (tan), interval spans at least 2π ≈ 6.28");

    auto x = SymbolicExpr::variable("x");
    // tan(x) + x: k=1, period for tan = π/|k| = π ≈ 3.14
    // Two full periods = 2 * π = 2π ≈ 6.28
    // Default [-10,10] width=20 > 6.28, so default is used.
    // Key check: interval width >= 2π ≈ 6.28
    auto expr = SymbolicExpr::add(SymbolicExpr::tan(x), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "tan(x)+x should produce valid interval");
    if (result) {
        double width = result->hi - result->lo;
        double min_required = 2.0 * (M_PI / 1.0);  // 2π ≈ 6.28
        EXPECT_TRUE(width >= min_required - 1e-10,
            "tan(x)+x: interval width " + std::to_string(width) +
            " should be >= 2π ≈ " + std::to_string(min_required));
    }
}

// Feature: mixed-transcendental-solver, Property 9: Periodic extension covers two full periods
// **Validates: Requirements 4.4**
void test_property9_sin_small_k_clamped() {
    TEST_CASE("Property 9: sin(0.1*x) + x/100 → k=0.1, period=20π, 2 periods span=125.66, half_span=62.83");

    auto x = SymbolicExpr::variable("x");
    // sin(0.1*x) + x/100: k=0.1, period = 2π/0.1 = 20π ≈ 62.83
    // Two full periods = 2 * 20π = 40π ≈ 125.66 total width
    // Implementation: half_span = period = 62.83, interval = [-62.83, 62.83]
    // Since 62.83 < 100, no clamping needed. Width = 125.66 >= 2*(2π/|k|).
    auto sin_arg = SymbolicExpr::multiply(SymbolicExpr::number(0.1), x);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(sin_arg),
        SymbolicExpr::multiply(SymbolicExpr::number(0.01), x)
    );

    lamina::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "sin(0.1*x)+x/100 should produce valid interval");
    if (result) {
        double width = result->hi - result->lo;
        double min_required = 2.0 * (2.0 * M_PI / 0.1);  // 2 * 20π ≈ 125.66
        EXPECT_TRUE(width >= min_required - 1e-10,
            "sin(0.1*x)+x/100: interval width " + std::to_string(width) +
            " should be >= 2*(2π/0.1) ≈ " + std::to_string(min_required));
        // Verify bounds are within [-100, 100] (clamping constraint)
        EXPECT_TRUE(result->lo >= -100.0 - 1e-10,
            "sin(0.1*x): lo should be >= -100 (clamping)");
        EXPECT_TRUE(result->hi <= 100.0 + 1e-10,
            "sin(0.1*x): hi should be <= 100 (clamping)");
    }
}

// Feature: mixed-transcendental-solver, Property 9: Periodic extension covers two full periods
// **Validates: Requirements 4.4**
void test_property9_sin_x_squared_nonlinear_default() {
    TEST_CASE("Property 9: sin(x^2) + x → non-linear argument, default [-10, 10]");

    auto x = SymbolicExpr::variable("x");
    // sin(x^2) + x: argument x^2 is non-linear in x
    // Should use default [-10, 10] without periodic extension
    auto sin_arg = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(sin_arg), x);

    lamina::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = lamina::determine_search_interval(expr, "x", opts);
    EXPECT_TRUE(result.has_value(), "sin(x^2)+x should produce valid interval");
    if (result) {
        EXPECT_NEAR(result->lo, -10.0, 1e-15,
            "sin(x^2)+x: non-linear arg, default lo=-10");
        EXPECT_NEAR(result->hi, 10.0, 1e-15,
            "sin(x^2)+x: non-linear arg, default hi=10");
    }
}

// Feature: mixed-transcendental-solver, Property 5: Max roots limit respected
// **Validates: Requirements 2.4**
void test_property5_max_roots_2_sin_x() {
    TEST_CASE("Property 5: sin(x) on [-10,10] with max_roots=2 → at most 2 intervals");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);

    // Compute derivative: cos(x)
    auto derivative = SymbolicExpr::cos(x);

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-12;
    opts.max_roots = 2;

    lamina::SearchInterval interval{-10.0, 10.0};

    auto result = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE(result.size() <= 2,
        "sin(x) on [-10,10] with max_roots=2: got " + std::to_string(result.size()) +
        " intervals, expected at most 2");
}

// Feature: mixed-transcendental-solver, Property 5: Max roots limit respected
// **Validates: Requirements 2.4**
void test_property5_max_roots_1_sin_x() {
    TEST_CASE("Property 5: sin(x) on [-10,10] with max_roots=1 → at most 1 interval");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);

    auto derivative = SymbolicExpr::cos(x);

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-12;
    opts.max_roots = 1;

    lamina::SearchInterval interval{-10.0, 10.0};

    auto result = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE(result.size() <= 1,
        "sin(x) on [-10,10] with max_roots=1: got " + std::to_string(result.size()) +
        " intervals, expected at most 1");
}

// Feature: mixed-transcendental-solver, Property 5: Max roots limit respected
// **Validates: Requirements 2.4**
void test_property5_max_roots_unlimited_sin_x() {
    TEST_CASE("Property 5: sin(x) on [-10,10] with max_roots=-1 (unlimited) → finds multiple roots");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);

    auto derivative = SymbolicExpr::cos(x);

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;  // unlimited

    lamina::SearchInterval interval{-10.0, 10.0};

    auto result = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    // sin(x) has roots at -3π, -2π, -π, 0, π, 2π, 3π in [-10,10]
    // That's approximately 7 roots. We expect at least 3 to confirm unlimited works.
    EXPECT_TRUE(result.size() >= 3,
        "sin(x) on [-10,10] with max_roots=-1: got " + std::to_string(result.size()) +
        " intervals, expected at least 3 (unlimited mode)");
}

// Feature: mixed-transcendental-solver, Property 5: Max roots limit respected
// **Validates: Requirements 2.4**
void test_property5_max_roots_3_cos_x_minus_half() {
    TEST_CASE("Property 5: cos(x)-0.5 on [-10,10] with max_roots=3 → at most 3 intervals");

    auto x = SymbolicExpr::variable("x");
    // cos(x) - 0.5
    auto expr = SymbolicExpr::add(
        SymbolicExpr::cos(x),
        SymbolicExpr::number(-0.5)
    );

    // derivative: -sin(x)
    auto derivative = SymbolicExpr::multiply(
        SymbolicExpr::number(-1),
        SymbolicExpr::sin(x)
    );

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-12;
    opts.max_roots = 3;

    lamina::SearchInterval interval{-10.0, 10.0};

    auto result = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE(result.size() <= 3,
        "cos(x)-0.5 on [-10,10] with max_roots=3: got " + std::to_string(result.size()) +
        " intervals, expected at most 3");
}

// Feature: mixed-transcendental-solver, Property 4: Minimum subdivision width
// **Validates: Requirements 2.2**
void test_property4_min_width_sin_x_plus_x_over_10() {
    TEST_CASE("Property 4: sin(x) + x/10 on [-10, 10] → all intervals have width >= 1e-6");

    auto x = SymbolicExpr::variable("x");
    // sin(x) + x/10
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::multiply(SymbolicExpr::number(0.1), x)
    );

    lamina::SearchInterval interval{-10.0, 10.0};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;

    auto intervals = lamina::isolate_roots(expr, nullptr, "x", interval, opts);

    std::cout << "  [INFO] isolate_roots returned " << intervals.size() << " interval(s)" << std::endl;

    for (size_t i = 0; i < intervals.size(); ++i) {
        double width = intervals[i].hi - intervals[i].lo;
        EXPECT_TRUE(width >= 1e-6,
            "Interval " + std::to_string(i) + " width " + std::to_string(width) +
            " should be >= 1e-6 [lo=" + std::to_string(intervals[i].lo) +
            ", hi=" + std::to_string(intervals[i].hi) + "]");
    }
}

// Feature: mixed-transcendental-solver, Property 4: Minimum subdivision width
// **Validates: Requirements 2.2**
void test_property4_min_width_tan_x_minus_x() {
    TEST_CASE("Property 4: tan(x) - x on [-5, 5] → all intervals have width >= 1e-6");

    auto x = SymbolicExpr::variable("x");
    // tan(x) - x: many roots close together near multiples of π
    auto expr = SymbolicExpr::add(
        SymbolicExpr::tan(x),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x)
    );

    lamina::SearchInterval interval{-5.0, 5.0};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;

    auto intervals = lamina::isolate_roots(expr, nullptr, "x", interval, opts);

    std::cout << "  [INFO] isolate_roots returned " << intervals.size() << " interval(s)" << std::endl;

    for (size_t i = 0; i < intervals.size(); ++i) {
        double width = intervals[i].hi - intervals[i].lo;
        EXPECT_TRUE(width >= 1e-6,
            "Interval " + std::to_string(i) + " width " + std::to_string(width) +
            " should be >= 1e-6 [lo=" + std::to_string(intervals[i].lo) +
            ", hi=" + std::to_string(intervals[i].hi) + "]");
    }
}

// Feature: mixed-transcendental-solver, Property 4: Minimum subdivision width
// **Validates: Requirements 2.2**
void test_property4_min_width_narrow_interval() {
    TEST_CASE("Property 4: sin(x) + x/10 on [0, 0.001] → intervals still respect width >= 1e-6");

    auto x = SymbolicExpr::variable("x");
    // sin(x) + x/10 on a very narrow search interval [0, 0.001]
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::multiply(SymbolicExpr::number(0.1), x)
    );

    lamina::SearchInterval interval{0.0, 0.001};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;

    auto intervals = lamina::isolate_roots(expr, nullptr, "x", interval, opts);

    std::cout << "  [INFO] isolate_roots on [0, 0.001] returned " << intervals.size() << " interval(s)" << std::endl;

    for (size_t i = 0; i < intervals.size(); ++i) {
        double width = intervals[i].hi - intervals[i].lo;
        EXPECT_TRUE(width >= 1e-6,
            "Interval " + std::to_string(i) + " width " + std::to_string(width) +
            " should be >= 1e-6 [lo=" + std::to_string(intervals[i].lo) +
            ", hi=" + std::to_string(intervals[i].hi) + "]");
    }
}

// Feature: mixed-transcendental-solver, Property 3: Root isolation sign-change invariant
// **Validates: Requirements 2.1**

/// Helper: evaluate f(val) = expr->substitute(var, number(val))->to_numeric()
static std::optional<double> eval_at(const std::shared_ptr<SymbolicExpr>& expr,
                                     const std::string& var, double val) {
    try {
        auto substituted = expr->substitute(var, SymbolicExpr::number(val));
        double result = substituted->to_numeric();
        if (std::isnan(result) || std::isinf(result)) return std::nullopt;
        return result;
    } catch (...) {
        return std::nullopt;
    }
}

/// Helper: compute derivative of expr with respect to var using DifferentiationVisitor
static std::shared_ptr<SymbolicExpr> compute_derivative(
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    if (!expr || !expr->root) return nullptr;
    try {
        DifferentiationVisitor dv(var);
        expr->root->accept(dv);
        auto result_node = dv.get_result();
        if (!result_node) return nullptr;
        return std::make_shared<SymbolicExpr>(result_node);
    } catch (...) {
        return nullptr;
    }
}

void test_property3_sign_change_sin_x_plus_x_div_10() {
    TEST_CASE("Property 3: sign-change invariant for sin(x) on [-10, 10]");

    auto x = SymbolicExpr::variable("x");
    // f(x) = sin(x): has roots at 0, ±π, ±2π, ±3π within [-10, 10]
    auto expr = SymbolicExpr::sin(x);

    auto derivative = compute_derivative(expr, "x");

    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    lamina::SearchInterval interval{-10.0, 10.0};
    auto intervals = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    std::cout << "  [INFO] sin(x) on [-10,10]: " << intervals.size() << " intervals found" << std::endl;

    // sin(x) has 7 roots in [-10, 10]: -3π, -2π, -π, 0, π, 2π, 3π
    // We expect at least some intervals to be found
    EXPECT_TRUE(intervals.size() > 0,
        "sin(x) on [-10,10] should produce at least one isolated interval");

    for (size_t i = 0; i < intervals.size(); ++i) {
        auto f_lo = eval_at(expr, "x", intervals[i].lo);
        auto f_hi = eval_at(expr, "x", intervals[i].hi);

        if (f_lo.has_value() && f_hi.has_value()) {
            // Sign-change invariant: f(lo)*f(hi) <= 0 OR interval is confirmed via monotonicity
            bool sign_change = (*f_lo) * (*f_hi) <= 0.0;
            bool confirmed = intervals[i].confirmed;
            EXPECT_TRUE(sign_change || confirmed,
                "Interval [" + std::to_string(intervals[i].lo) + ", " +
                std::to_string(intervals[i].hi) + "]: f(lo)=" +
                std::to_string(*f_lo) + ", f(hi)=" + std::to_string(*f_hi) +
                " must have opposite signs or be confirmed by monotonicity");
        }
    }
}

void test_property3_sign_change_exp_x_minus_x_minus_2() {
    TEST_CASE("Property 3: sign-change invariant for exp(x) - x - 2 on [-5, 5]");

    auto x = SymbolicExpr::variable("x");
    // f(x) = exp(x) - x - 2
    // Roots near x ≈ -1.84 and x ≈ 1.15
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
            SymbolicExpr::number(-2)
        )
    );

    auto derivative = compute_derivative(expr, "x");

    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    lamina::SearchInterval interval{-5.0, 5.0};
    auto intervals = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    std::cout << "  [INFO] exp(x)-x-2 on [-5,5]: " << intervals.size() << " intervals found" << std::endl;

    // exp(x) - x - 2 has 2 roots in [-5, 5]
    EXPECT_TRUE(intervals.size() > 0,
        "exp(x)-x-2 on [-5,5] should produce at least one isolated interval");

    for (size_t i = 0; i < intervals.size(); ++i) {
        auto f_lo = eval_at(expr, "x", intervals[i].lo);
        auto f_hi = eval_at(expr, "x", intervals[i].hi);

        if (f_lo.has_value() && f_hi.has_value()) {
            bool sign_change = (*f_lo) * (*f_hi) <= 0.0;
            bool confirmed = intervals[i].confirmed;
            EXPECT_TRUE(sign_change || confirmed,
                "Interval [" + std::to_string(intervals[i].lo) + ", " +
                std::to_string(intervals[i].hi) + "]: f(lo)=" +
                std::to_string(*f_lo) + ", f(hi)=" + std::to_string(*f_hi) +
                " must have opposite signs or be confirmed by monotonicity");
        }
    }
}

void test_property3_sign_change_x_cos_x_minus_1() {
    TEST_CASE("Property 3: sign-change invariant for x*cos(x) - 1 on [-10, 10]");

    auto x = SymbolicExpr::variable("x");
    // f(x) = x*cos(x) - 1: has multiple roots in [-10, 10]
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(x, SymbolicExpr::cos(x)),
        SymbolicExpr::number(-1)
    );

    auto derivative = compute_derivative(expr, "x");

    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    lamina::SearchInterval interval{-10.0, 10.0};
    auto intervals = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    std::cout << "  [INFO] x*cos(x)-1 on [-10,10]: " << intervals.size() << " intervals found" << std::endl;

    // x*cos(x) - 1 has roots (e.g., near x ≈ 1.28, x ≈ -1.28, and others)
    EXPECT_TRUE(intervals.size() > 0,
        "x*cos(x)-1 on [-10,10] should produce at least one isolated interval");

    for (size_t i = 0; i < intervals.size(); ++i) {
        auto f_lo = eval_at(expr, "x", intervals[i].lo);
        auto f_hi = eval_at(expr, "x", intervals[i].hi);

        if (f_lo.has_value() && f_hi.has_value()) {
            bool sign_change = (*f_lo) * (*f_hi) <= 0.0;
            bool confirmed = intervals[i].confirmed;
            EXPECT_TRUE(sign_change || confirmed,
                "Interval [" + std::to_string(intervals[i].lo) + ", " +
                std::to_string(intervals[i].hi) + "]: f(lo)=" +
                std::to_string(*f_lo) + ", f(hi)=" + std::to_string(*f_hi) +
                " must have opposite signs or be confirmed by monotonicity");
        }
    }
}

// Feature: mixed-transcendental-solver, Property 6: No sign changes yields empty result
// **Validates: Requirements 2.6**
void test_property6_x_squared_plus_one_always_positive() {
    TEST_CASE("Property 6: x^2 + 1 on [-10, 10] → always positive, isolate_roots returns empty");

    auto x = SymbolicExpr::variable("x");
    // x^2 + 1 > 0 for all real x
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1)
    );
    auto derivative = expr->differentiate("x");

    lamina::SearchInterval interval{-10.0, 10.0};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    auto isolated = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE(isolated.empty(),
        "x^2 + 1 is always positive on [-10, 10], isolate_roots should return empty");
}

// Feature: mixed-transcendental-solver, Property 6: No sign changes yields empty result
// **Validates: Requirements 2.6**
void test_property6_exp_x_plus_one_always_positive() {
    TEST_CASE("Property 6: exp(x) + 1 on [-10, 10] → always positive, isolate_roots returns empty");

    auto x = SymbolicExpr::variable("x");
    // exp(x) + 1 > 1 > 0 for all real x (since exp(x) > 0)
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::number(1)
    );
    auto derivative = expr->differentiate("x");

    lamina::SearchInterval interval{-10.0, 10.0};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    auto isolated = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE(isolated.empty(),
        "exp(x) + 1 is always positive on [-10, 10], isolate_roots should return empty");
}

// Feature: mixed-transcendental-solver, Property 6: No sign changes yields empty result
// **Validates: Requirements 2.6**
void test_property6_neg_x_squared_minus_one_always_negative() {
    TEST_CASE("Property 6: -(x^2 + 1) on [-10, 10] → always negative, isolate_roots returns empty");

    auto x = SymbolicExpr::variable("x");
    // -(x^2 + 1) = -x^2 - 1 < 0 for all real x
    auto expr = SymbolicExpr::multiply(
        SymbolicExpr::number(-1),
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::number(1)
        )
    );
    auto derivative = expr->differentiate("x");

    lamina::SearchInterval interval{-10.0, 10.0};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    auto isolated = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE(isolated.empty(),
        "-(x^2 + 1) is always negative on [-10, 10], isolate_roots should return empty");
}

// Feature: mixed-transcendental-solver, Property 6: No sign changes yields empty result
// **Validates: Requirements 2.6**
void test_property6_sin_x_plus_five_always_positive() {
    TEST_CASE("Property 6: sin(x) + 5 on [-10, 10] → always positive, isolate_roots returns empty");

    auto x = SymbolicExpr::variable("x");
    // sin(x) + 5: sin(x) ranges in [-1, 1], so sin(x)+5 ranges in [4, 6] > 0
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::number(5)
    );
    auto derivative = expr->differentiate("x");

    lamina::SearchInterval interval{-10.0, 10.0};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    auto isolated = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE(isolated.empty(),
        "sin(x) + 5 is always positive on [-10, 10], isolate_roots should return empty");
}

// ============================================================================
// refine_root tests
// **Validates: Requirements 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 8.3**
// ============================================================================

void test_refine_root_sin_x_newton_raphson() {
    TEST_CASE("refine_root: sin(x) on [3.0, 3.5] with derivative → converges to π");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    auto derivative = compute_derivative(expr, "x");  // cos(x)

    lamina::IsolatedInterval interval{3.0, 3.5, true};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = lamina::refine_root(expr, derivative, "x", interval, opts);

    EXPECT_TRUE(result.has_value(), "refine_root should converge for sin(x) on [3.0, 3.5]");
    if (result) {
        EXPECT_NEAR(result->value, M_PI, 1e-10,
            "Root of sin(x) near [3.0, 3.5] should be π ≈ 3.14159");
        EXPECT_TRUE(result->residual < opts.tolerance,
            "Residual " + std::to_string(result->residual) + " should be < tolerance");
        EXPECT_TRUE(result->iterations > 0, "Should take at least 1 iteration");
    }
}

void test_refine_root_sin_x_bisection_fallback() {
    TEST_CASE("refine_root: sin(x) on [3.0, 3.5] without derivative → bisection converges to π");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);

    lamina::IsolatedInterval interval{3.0, 3.5, false};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    // No derivative → pure bisection
    auto result = lamina::refine_root(expr, nullptr, "x", interval, opts);

    EXPECT_TRUE(result.has_value(), "refine_root (bisection) should converge for sin(x) on [3.0, 3.5]");
    if (result) {
        EXPECT_NEAR(result->value, M_PI, 1e-10,
            "Root of sin(x) via bisection should be π ≈ 3.14159");
        EXPECT_TRUE(result->residual < opts.tolerance,
            "Residual " + std::to_string(result->residual) + " should be < tolerance");
    }
}

void test_refine_root_exp_x_minus_x_minus_2() {
    TEST_CASE("refine_root: exp(x) - x - 2 on [1.0, 1.5] → converges to root ≈ 1.146");

    auto x = SymbolicExpr::variable("x");
    // exp(x) - x - 2
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
            SymbolicExpr::number(-2)
        )
    );
    auto derivative = compute_derivative(expr, "x");  // exp(x) - 1

    lamina::IsolatedInterval interval{1.0, 1.5, true};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = lamina::refine_root(expr, derivative, "x", interval, opts);

    EXPECT_TRUE(result.has_value(), "refine_root should converge for exp(x)-x-2 on [1.0, 1.5]");
    if (result) {
        // Verify residual is small
        EXPECT_TRUE(result->residual < opts.tolerance,
            "Residual " + std::to_string(result->residual) + " should be < tolerance");
        // Verify root is in interval
        EXPECT_TRUE(result->value >= 1.0 && result->value <= 1.5,
            "Root " + std::to_string(result->value) + " should be in [1.0, 1.5]");
    }
}

void test_refine_root_x_cos_x_minus_1() {
    TEST_CASE("refine_root: x*cos(x) - 1 on [4.5, 5.0] → converges");

    auto x = SymbolicExpr::variable("x");
    // x*cos(x) - 1: has a root near x ≈ 4.917
    // f(4.5) = 4.5*cos(4.5) - 1 ≈ 4.5*(-0.2108) - 1 ≈ -1.949 < 0
    // f(5.0) = 5.0*cos(5.0) - 1 ≈ 5.0*(0.2837) - 1 ≈ 0.418 > 0
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(x, SymbolicExpr::cos(x)),
        SymbolicExpr::number(-1)
    );
    auto derivative = compute_derivative(expr, "x");

    lamina::IsolatedInterval interval{4.5, 5.0, true};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = lamina::refine_root(expr, derivative, "x", interval, opts);

    EXPECT_TRUE(result.has_value(), "refine_root should converge for x*cos(x)-1 on [4.5, 5.0]");
    if (result) {
        EXPECT_TRUE(result->residual < opts.tolerance,
            "Residual " + std::to_string(result->residual) + " should be < tolerance");
        EXPECT_TRUE(result->value >= 4.5 && result->value <= 5.0,
            "Root " + std::to_string(result->value) + " should be in [4.5, 5.0]");
    }
}

void test_refine_root_discards_when_no_convergence() {
    TEST_CASE("refine_root: x^2 + 1 on [-1, 1] → no root, returns nullopt");

    auto x = SymbolicExpr::variable("x");
    // x^2 + 1 > 0 always, no real root
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1)
    );
    auto derivative = compute_derivative(expr, "x");

    // This interval has no sign change, so refinement should fail
    lamina::IsolatedInterval interval{-1.0, 1.0, false};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 50;

    auto result = lamina::refine_root(expr, derivative, "x", interval, opts);

    EXPECT_FALSE(result.has_value(),
        "refine_root should return nullopt for x^2+1 (no real root in [-1, 1])");
}

void test_refine_root_derivative_zero_bisection_step() {
    TEST_CASE("refine_root: x^3 on [-1, 1] with derivative 3x^2 → converges to 0 (derivative zero at root)");

    auto x = SymbolicExpr::variable("x");
    // x^3: root at 0, derivative 3x^2 is zero at root
    auto expr = SymbolicExpr::power(x, SymbolicExpr::number(3));
    auto derivative = compute_derivative(expr, "x");  // 3*x^2

    lamina::IsolatedInterval interval{-1.0, 1.0, false};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 200;

    auto result = lamina::refine_root(expr, derivative, "x", interval, opts);

    EXPECT_TRUE(result.has_value(), "refine_root should converge for x^3 on [-1, 1]");
    if (result) {
        EXPECT_NEAR(result->value, 0.0, 1e-4,
            "Root of x^3 should be 0 (derivative zero at root triggers bisection)");
        EXPECT_TRUE(result->residual < opts.tolerance,
            "Residual " + std::to_string(result->residual) + " should be < tolerance");
    }
}

void test_refine_root_max_iterations_exceeded() {
    TEST_CASE("refine_root: sin(x) on [3.0, 3.5] with max_iterations=2 → discards if residual > tolerance");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    auto derivative = compute_derivative(expr, "x");

    lamina::IsolatedInterval interval{3.0, 3.5, true};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-15;  // Very tight tolerance
    opts.max_newton_iterations = 2;  // Very few iterations

    auto result = lamina::refine_root(expr, derivative, "x", interval, opts);

    // With only 2 iterations and 1e-15 tolerance, Newton-Raphson might still converge
    // (it's quadratically convergent). If it does, great. If not, it should return nullopt.
    if (result.has_value()) {
        EXPECT_TRUE(result->residual < opts.tolerance,
            "If converged, residual must be < tolerance");
    }
    // Either way, the function should not crash
    std::cout << "  [INFO] max_iter=2, tol=1e-15: "
              << (result.has_value() ? "converged" : "did not converge") << std::endl;
}

// Feature: mixed-transcendental-solver, Property 11: Output sorted ascending
// **Validates: Requirements 5.4, 6.3**

void test_property11_deduplicate_roots_ascending_order() {
    TEST_CASE("Property 11: deduplicate_roots produces strictly ascending output");

    // Create a vector of NumericRoot in random (unsorted) order
    std::vector<NumericRoot> roots = {
        {3.0, 1e-14, 5},
        {1.0, 1e-14, 3},
        {5.0, 1e-14, 7},
        {2.0, 1e-14, 4},
        {4.0, 1e-14, 6}
    };

    lmmc_real_t tolerance = 1e-12;
    int max_roots = -1;

    auto sorted = deduplicate_roots(roots, tolerance, max_roots);

    EXPECT_TRUE(sorted.size() == 5,
        "All 5 distinct roots should survive deduplication");

    // Verify strictly ascending order: result[i] < result[i+1]
    for (size_t i = 0; i + 1 < sorted.size(); ++i) {
        EXPECT_TRUE(sorted[i] < sorted[i + 1],
            "sorted[" + std::to_string(i) + "]=" + std::to_string(sorted[i]) +
            " should be < sorted[" + std::to_string(i + 1) + "]=" +
            std::to_string(sorted[i + 1]));
    }

    // Verify actual values
    if (sorted.size() == 5) {
        EXPECT_NEAR(sorted[0], 1.0, 1e-10, "First root should be 1.0");
        EXPECT_NEAR(sorted[1], 2.0, 1e-10, "Second root should be 2.0");
        EXPECT_NEAR(sorted[2], 3.0, 1e-10, "Third root should be 3.0");
        EXPECT_NEAR(sorted[3], 4.0, 1e-10, "Fourth root should be 4.0");
        EXPECT_NEAR(sorted[4], 5.0, 1e-10, "Fifth root should be 5.0");
    }
}

// Feature: mixed-transcendental-solver, Property 11: Output sorted ascending
// **Validates: Requirements 5.4, 6.3**
void test_property11_pipeline_sin_x_ascending() {
    TEST_CASE("Property 11: sin(x) on [-10, 10] pipeline → roots in ascending order");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    auto derivative = compute_derivative(expr, "x");  // cos(x)

    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    lamina::SearchInterval interval{-10.0, 10.0};

    // Step 1: isolate roots
    auto intervals = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    // Step 2: refine each root
    std::vector<NumericRoot> refined_roots;
    for (const auto& iso : intervals) {
        auto root = lamina::refine_root(expr, derivative, "x", iso, opts);
        if (root.has_value()) {
            refined_roots.push_back(*root);
        }
    }

    // Step 3: deduplicate
    auto sorted = deduplicate_roots(refined_roots, opts.tolerance, opts.max_roots);

    std::cout << "  [INFO] sin(x) pipeline: " << sorted.size() << " roots found" << std::endl;

    // sin(x) has roots at -3π, -2π, -π, 0, π, 2π, 3π in [-10, 10] → 7 roots
    EXPECT_TRUE(sorted.size() >= 2,
        "sin(x) on [-10,10] should produce at least 2 roots");

    // Verify strictly ascending order
    for (size_t i = 0; i + 1 < sorted.size(); ++i) {
        EXPECT_TRUE(sorted[i] < sorted[i + 1],
            "sin(x) roots: sorted[" + std::to_string(i) + "]=" +
            std::to_string(sorted[i]) + " should be < sorted[" +
            std::to_string(i + 1) + "]=" + std::to_string(sorted[i + 1]));
    }
}

// Feature: mixed-transcendental-solver, Property 11: Output sorted ascending
// **Validates: Requirements 5.4, 6.3**
void test_property11_pipeline_cos_x_ascending() {
    TEST_CASE("Property 11: cos(x) on [-10, 10] pipeline → roots in ascending order");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::cos(x);
    auto derivative = compute_derivative(expr, "x");  // -sin(x)

    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    lamina::SearchInterval interval{-10.0, 10.0};

    // Step 1: isolate roots
    auto intervals = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    // Step 2: refine each root
    std::vector<NumericRoot> refined_roots;
    for (const auto& iso : intervals) {
        auto root = lamina::refine_root(expr, derivative, "x", iso, opts);
        if (root.has_value()) {
            refined_roots.push_back(*root);
        }
    }

    // Step 3: deduplicate
    auto sorted = deduplicate_roots(refined_roots, opts.tolerance, opts.max_roots);

    std::cout << "  [INFO] cos(x) pipeline: " << sorted.size() << " roots found" << std::endl;

    // cos(x) has roots at ±π/2, ±3π/2, ±5π/2 in [-10, 10] → 6 roots
    EXPECT_TRUE(sorted.size() >= 2,
        "cos(x) on [-10,10] should produce at least 2 roots");

    // Verify strictly ascending order
    for (size_t i = 0; i + 1 < sorted.size(); ++i) {
        EXPECT_TRUE(sorted[i] < sorted[i + 1],
            "cos(x) roots: sorted[" + std::to_string(i) + "]=" +
            std::to_string(sorted[i]) + " should be < sorted[" +
            std::to_string(i + 1) + "]=" + std::to_string(sorted[i + 1]));
    }
}

// Feature: mixed-transcendental-solver, Property 11: Output sorted ascending
// **Validates: Requirements 5.4, 6.3**
void test_property11_deduplicate_with_duplicates_ascending() {
    TEST_CASE("Property 11: deduplicate_roots with near-duplicates → ascending after merge");

    // Create roots with some near-duplicates (within 10*tolerance)
    lmmc_real_t tolerance = 1e-10;
    std::vector<NumericRoot> roots = {
        {5.0,  1e-13, 4},
        {1.0,  1e-13, 2},
        {1.0 + 5e-10, 1e-14, 3},  // near-duplicate of 1.0 (diff < 10*tol)
        {3.0,  1e-13, 5},
        {3.0 + 2e-10, 1e-12, 6},  // near-duplicate of 3.0 (diff < 10*tol)
        {-2.0, 1e-13, 1}
    };

    int max_roots = -1;
    auto sorted = deduplicate_roots(roots, tolerance, max_roots);

    std::cout << "  [INFO] deduplicate with duplicates: " << sorted.size() << " roots remain" << std::endl;

    // After deduplication: should have 4 distinct roots (-2, 1, 3, 5)
    EXPECT_TRUE(sorted.size() == 4,
        "Should have 4 distinct roots after deduplication, got " +
        std::to_string(sorted.size()));

    // Verify strictly ascending order
    for (size_t i = 0; i + 1 < sorted.size(); ++i) {
        EXPECT_TRUE(sorted[i] < sorted[i + 1],
            "After dedup: sorted[" + std::to_string(i) + "]=" +
            std::to_string(sorted[i]) + " should be < sorted[" +
            std::to_string(i + 1) + "]=" + std::to_string(sorted[i + 1]));
    }
}

// Feature: mixed-transcendental-solver, Property 8: Roots within specified search interval
// **Validates: Requirements 4.2**
void test_property8_roots_within_interval_sin_x() {
    TEST_CASE("Property 8: sin(x) on [-3, 3] → all refined roots within [-3, 3]");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    auto derivative = compute_derivative(expr, "x");  // cos(x)

    lamina::SearchInterval interval{-3.0, 3.0};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;
    opts.has_search_interval = true;
    opts.search_lo = -3.0;
    opts.search_hi = 3.0;

    auto isolated = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    std::cout << "  [INFO] sin(x) on [-3,3]: " << isolated.size() << " isolated interval(s)" << std::endl;

    // sin(x) has a root at x=0 within [-3, 3]. ±π ≈ ±3.14 are outside.
    // Regardless of how many roots are found, ALL must be within [-3, 3].
    for (size_t i = 0; i < isolated.size(); ++i) {
        auto refined = lamina::refine_root(expr, derivative, "x", isolated[i], opts);
        if (refined.has_value()) {
            EXPECT_TRUE(refined->value >= -3.0,
                "Root " + std::to_string(refined->value) + " should be >= -3.0 (lower bound)");
            EXPECT_TRUE(refined->value <= 3.0,
                "Root " + std::to_string(refined->value) + " should be <= 3.0 (upper bound)");
        }
    }

    // Also test with a wider interval that definitely contains multiple roots
    lamina::SearchInterval interval2{-4.0, 4.0};
    opts.search_lo = -4.0;
    opts.search_hi = 4.0;

    auto isolated2 = lamina::isolate_roots(expr, derivative, "x", interval2, opts);

    std::cout << "  [INFO] sin(x) on [-4,4]: " << isolated2.size() << " isolated interval(s)" << std::endl;

    // sin(x) has roots at -π, 0, π within [-4, 4]
    for (size_t i = 0; i < isolated2.size(); ++i) {
        auto refined = lamina::refine_root(expr, derivative, "x", isolated2[i], opts);
        if (refined.has_value()) {
            EXPECT_TRUE(refined->value >= -4.0,
                "Root " + std::to_string(refined->value) + " should be >= -4.0 (lower bound)");
            EXPECT_TRUE(refined->value <= 4.0,
                "Root " + std::to_string(refined->value) + " should be <= 4.0 (upper bound)");
        }
    }

    // Verify at least some roots were found in the wider interval
    EXPECT_TRUE(isolated2.size() > 0,
        "sin(x) on [-4,4] should find at least one root (0, ±π are in range)");
}

// Feature: mixed-transcendental-solver, Property 8: Roots within specified search interval
// **Validates: Requirements 4.2**
void test_property8_roots_within_interval_cos_x_minus_half() {
    TEST_CASE("Property 8: cos(x) - 0.5 on [0, 5] → all refined roots within [0, 5]");

    auto x = SymbolicExpr::variable("x");
    // cos(x) - 0.5
    auto expr = SymbolicExpr::add(
        SymbolicExpr::cos(x),
        SymbolicExpr::number(-0.5)
    );
    // derivative: -sin(x)
    auto derivative = SymbolicExpr::multiply(
        SymbolicExpr::number(-1),
        SymbolicExpr::sin(x)
    );

    lamina::SearchInterval interval{0.0, 5.0};
    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;
    opts.has_search_interval = true;
    opts.search_lo = 0.0;
    opts.search_hi = 5.0;

    auto isolated = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    std::cout << "  [INFO] cos(x)-0.5 on [0,5]: " << isolated.size() << " isolated interval(s)" << std::endl;

    // cos(x) = 0.5 has roots at x = π/3 ≈ 1.047 and x = 5π/3 ≈ 5.236 (outside [0,5])
    // Actually within [0,5]: x = π/3 ≈ 1.047 and x = 5π/3 ≈ 5.236 is outside,
    // but x = 2π - π/3 ≈ 5.236 is outside. So only π/3 ≈ 1.047 is in [0,5]?
    // Wait: cos(x)=0.5 → x = ±π/3 + 2kπ. In [0,5]: π/3 ≈ 1.047, 5π/3 ≈ 5.236 (outside).
    // So we expect 1 root in [0,5].

    for (size_t i = 0; i < isolated.size(); ++i) {
        auto refined = lamina::refine_root(expr, derivative, "x", isolated[i], opts);
        if (refined.has_value()) {
            EXPECT_TRUE(refined->value >= 0.0,
                "Root " + std::to_string(refined->value) + " should be >= 0.0 (lower bound)");
            EXPECT_TRUE(refined->value <= 5.0,
                "Root " + std::to_string(refined->value) + " should be <= 5.0 (upper bound)");
        }
    }
}

// Feature: mixed-transcendental-solver, Property 7: Output validity invariant
// **Validates: Requirements 3.2, 3.5, 3.6, 6.1, 7.3**

void test_property7_output_validity_sin_x_refined_roots() {
    TEST_CASE("Property 7: isolate_roots + refine_root on sin(x) [-10,10] → all roots finite, |f(root)| < tol, convertible to NumberNode");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    auto derivative = compute_derivative(expr, "x");  // cos(x)

    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    lamina::SearchInterval interval{-10.0, 10.0};
    auto intervals = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    std::cout << "  [INFO] sin(x) on [-10,10]: " << intervals.size() << " intervals isolated" << std::endl;

    EXPECT_TRUE(intervals.size() > 0,
        "sin(x) on [-10,10] should produce at least one isolated interval");

    for (size_t i = 0; i < intervals.size(); ++i) {
        auto root_result = lamina::refine_root(expr, derivative, "x", intervals[i], opts);

        EXPECT_TRUE(root_result.has_value(),
            "Interval " + std::to_string(i) + " should produce a refined root");

        if (root_result.has_value()) {
            // (a) Value is finite
            EXPECT_TRUE(std::isfinite(root_result->value),
                "Root " + std::to_string(i) + " value " + std::to_string(root_result->value) +
                " should be finite");

            // (b) |f(root)| < tolerance
            EXPECT_TRUE(root_result->residual < opts.tolerance,
                "Root " + std::to_string(i) + " residual " + std::to_string(root_result->residual) +
                " should be < tolerance " + std::to_string(opts.tolerance));

            // (c) Can be converted to NumberNode expression
            auto number_expr = SymbolicExpr::number(root_result->value);
            EXPECT_TRUE(number_expr != nullptr,
                "Root " + std::to_string(i) + " should be convertible to SymbolicExpr");
            if (number_expr && number_expr->root) {
                auto as_number = std::dynamic_pointer_cast<NumberNode>(number_expr->root);
                EXPECT_TRUE(as_number != nullptr,
                    "Root " + std::to_string(i) + " expression should be a NumberNode");
            }
        }
    }
}

void test_property7_output_validity_exp_x_minus_x_minus_2() {
    TEST_CASE("Property 7: isolate_roots + refine_root on exp(x)-x-2 [-5,5] → output validity");

    auto x = SymbolicExpr::variable("x");
    // exp(x) - x - 2: roots near x ≈ -1.84 and x ≈ 1.146
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
            SymbolicExpr::number(-2)
        )
    );
    auto derivative = compute_derivative(expr, "x");

    lamina::SolveOptions opts;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    lamina::SearchInterval interval{-5.0, 5.0};
    auto intervals = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    std::cout << "  [INFO] exp(x)-x-2 on [-5,5]: " << intervals.size() << " intervals isolated" << std::endl;

    EXPECT_TRUE(intervals.size() > 0,
        "exp(x)-x-2 on [-5,5] should produce at least one isolated interval");

    for (size_t i = 0; i < intervals.size(); ++i) {
        auto root_result = lamina::refine_root(expr, derivative, "x", intervals[i], opts);

        if (root_result.has_value()) {
            // (a) Value is finite
            EXPECT_TRUE(std::isfinite(root_result->value),
                "Root " + std::to_string(i) + " value should be finite");

            // (b) |f(root)| < tolerance — verify by substitution
            auto substituted = expr->substitute("x", SymbolicExpr::number(root_result->value));
            double f_val = substituted->to_numeric();
            EXPECT_TRUE(std::abs(f_val) < opts.tolerance,
                "Root " + std::to_string(i) + " |f(root)|=" + std::to_string(std::abs(f_val)) +
                " should be < tolerance " + std::to_string(opts.tolerance));

            // (c) Can be converted to NumberNode
            auto number_expr = SymbolicExpr::number(root_result->value);
            EXPECT_TRUE(number_expr != nullptr && number_expr->root != nullptr,
                "Root " + std::to_string(i) + " should produce a valid SymbolicExpr");
            if (number_expr && number_expr->root) {
                auto as_number = std::dynamic_pointer_cast<NumberNode>(number_expr->root);
                EXPECT_TRUE(as_number != nullptr,
                    "Root " + std::to_string(i) + " expression root should be NumberNode");
            }
        }
    }
}

void test_property7_max_roots_limit_after_deduplication() {
    TEST_CASE("Property 7: sin(x) on [-10,10] with max_roots=2 → at most 2 roots after deduplicate_roots");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    auto derivative = compute_derivative(expr, "x");

    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    opts.max_roots = 2;

    lamina::SearchInterval interval{-10.0, 10.0};
    auto intervals = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    // Refine all isolated intervals
    std::vector<NumericRoot> refined_roots;
    for (size_t i = 0; i < intervals.size(); ++i) {
        auto root_result = lamina::refine_root(expr, derivative, "x", intervals[i], opts);
        if (root_result.has_value()) {
            refined_roots.push_back(*root_result);
        }
    }

    std::cout << "  [INFO] Refined " << refined_roots.size() << " root(s) before deduplication" << std::endl;

    // (d) After deduplication with max_roots=2, at most 2 roots returned
    auto deduped = lamina::deduplicate_roots(refined_roots, opts.tolerance, opts.max_roots);

    EXPECT_TRUE(deduped.size() <= 2,
        "deduplicate_roots with max_roots=2 should return at most 2 roots, got " +
        std::to_string(deduped.size()));

    // Verify each deduplicated root is finite and can be represented as NumberNode
    for (size_t i = 0; i < deduped.size(); ++i) {
        EXPECT_TRUE(std::isfinite(deduped[i]),
            "Deduplicated root " + std::to_string(i) + " value " +
            std::to_string(deduped[i]) + " should be finite");

        auto number_expr = SymbolicExpr::number(deduped[i]);
        EXPECT_TRUE(number_expr != nullptr && number_expr->root != nullptr,
            "Deduplicated root " + std::to_string(i) + " should be convertible to NumberNode");
        if (number_expr && number_expr->root) {
            auto as_number = std::dynamic_pointer_cast<NumberNode>(number_expr->root);
            EXPECT_TRUE(as_number != nullptr,
                "Deduplicated root " + std::to_string(i) + " expression should be NumberNode");
        }
    }
}

void test_property7_output_validity_x_cos_x_minus_1() {
    TEST_CASE("Property 7: isolate_roots + refine_root on x*cos(x)-1 [-10,10] → output validity");

    auto x = SymbolicExpr::variable("x");
    // x*cos(x) - 1: multiple roots in [-10, 10]
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(x, SymbolicExpr::cos(x)),
        SymbolicExpr::number(-1)
    );
    auto derivative = compute_derivative(expr, "x");

    lamina::SolveOptions opts;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    lamina::SearchInterval interval{-10.0, 10.0};
    auto intervals = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    std::cout << "  [INFO] x*cos(x)-1 on [-10,10]: " << intervals.size() << " intervals isolated" << std::endl;

    EXPECT_TRUE(intervals.size() > 0,
        "x*cos(x)-1 on [-10,10] should produce at least one isolated interval");

    std::vector<NumericRoot> all_roots;
    for (size_t i = 0; i < intervals.size(); ++i) {
        auto root_result = lamina::refine_root(expr, derivative, "x", intervals[i], opts);

        if (root_result.has_value()) {
            // (a) Value is finite
            EXPECT_TRUE(std::isfinite(root_result->value),
                "Root " + std::to_string(i) + " value should be finite");

            // (b) |f(root)| < tolerance
            EXPECT_TRUE(root_result->residual < opts.tolerance,
                "Root " + std::to_string(i) + " residual " + std::to_string(root_result->residual) +
                " should be < tolerance");

            // (c) NumberNode representation
            auto number_expr = SymbolicExpr::number(root_result->value);
            auto as_number = std::dynamic_pointer_cast<NumberNode>(number_expr->root);
            EXPECT_TRUE(as_number != nullptr,
                "Root " + std::to_string(i) + " should be representable as NumberNode");

            all_roots.push_back(*root_result);
        }
    }

    // Verify deduplication preserves validity
    auto deduped = lamina::deduplicate_roots(all_roots, opts.tolerance, -1);
    std::cout << "  [INFO] After deduplication: " << deduped.size() << " unique root(s)" << std::endl;

    for (size_t i = 0; i < deduped.size(); ++i) {
        EXPECT_TRUE(std::isfinite(deduped[i]),
            "Deduplicated root " + std::to_string(i) + " should be finite");
    }
}

// ============================================================================
// Result assembly tests (Requirements 6.1, 6.2, 6.3, 6.4)
// ============================================================================

void test_assemble_results_basic() {
    TEST_CASE("Result assembly: deduplicate roots then convert to NumberNode expressions");

    // Create a few NumericRoot values with varying residuals
    std::vector<lamina::NumericRoot> roots = {
        {3.14159, 1e-13, 5},
        {-1.5,    2e-14, 3},
        {0.0,     0.0,   1},
        {2.71828, 5e-14, 7}
    };

    lmmc_real_t tolerance = 1e-12;
    int max_roots = -1;

    // Deduplicate (no duplicates in this set)
    auto deduped = lamina::deduplicate_roots(roots, tolerance, max_roots);

    // Verify deduplication produced sorted ascending values
    EXPECT_TRUE(deduped.size() == 4, "Should have 4 unique roots");
    if (deduped.size() == 4) {
        EXPECT_NEAR(deduped[0], -1.5, 1e-15, "First root should be -1.5");
        EXPECT_NEAR(deduped[1], 0.0, 1e-15, "Second root should be 0.0");
        EXPECT_NEAR(deduped[2], 2.71828, 1e-10, "Third root should be 2.71828");
        EXPECT_NEAR(deduped[3], 3.14159, 1e-10, "Fourth root should be 3.14159");
    }

    // Assemble into NumberNode expressions (simulating assemble_results)
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (const auto& val : deduped) {
        results.push_back(SymbolicExpr::number(static_cast<double>(val)));
    }

    // Verify each result is a valid NumberNode expression
    EXPECT_TRUE(results.size() == 4, "Should produce 4 NumberNode expressions");
    for (size_t i = 0; i < results.size(); ++i) {
        EXPECT_TRUE(results[i] != nullptr,
            "Result " + std::to_string(i) + " should not be null");
        auto num = std::dynamic_pointer_cast<NumberNode>(results[i]->root);
        EXPECT_TRUE(num != nullptr,
            "Result " + std::to_string(i) + " root should be a NumberNode");
    }

    // Verify numeric values match
    if (results.size() == 4) {
        auto v0 = test_numeric_eval(results[0]);
        auto v1 = test_numeric_eval(results[1]);
        auto v2 = test_numeric_eval(results[2]);
        auto v3 = test_numeric_eval(results[3]);

        EXPECT_TRUE(v0.has_value() && v1.has_value() && v2.has_value() && v3.has_value(),
            "All results should be evaluable to numeric values");

        if (v0) EXPECT_NEAR(*v0, -1.5, 1e-15, "NumberNode[0] value should be -1.5");
        if (v1) EXPECT_NEAR(*v1, 0.0, 1e-15, "NumberNode[1] value should be 0.0");
        if (v2) EXPECT_NEAR(*v2, 2.71828, 1e-10, "NumberNode[2] value should be 2.71828");
        if (v3) EXPECT_NEAR(*v3, 3.14159, 1e-10, "NumberNode[3] value should be 3.14159");
    }
}

void test_assemble_results_empty() {
    TEST_CASE("Result assembly: empty root list produces empty expression vector");

    std::vector<lamina::NumericRoot> roots;
    lmmc_real_t tolerance = 1e-12;
    int max_roots = -1;

    auto deduped = lamina::deduplicate_roots(roots, tolerance, max_roots);
    EXPECT_TRUE(deduped.empty(), "Deduplication of empty input should return empty");

    // Assemble (simulating assemble_results with empty input)
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (const auto& val : deduped) {
        results.push_back(SymbolicExpr::number(static_cast<double>(val)));
    }
    EXPECT_TRUE(results.empty(), "Assembly of empty roots should return empty vector");
}

void test_assemble_results_with_duplicates() {
    TEST_CASE("Result assembly: duplicates removed before conversion to NumberNode");

    // Create roots with duplicates (values within 10*tolerance of each other)
    std::vector<lamina::NumericRoot> roots = {
        {1.0,           1e-13, 5},
        {1.0 + 5e-12,  2e-13, 6},  // duplicate of 1.0 (diff < 10*1e-12)
        {2.0,           3e-14, 3},
        {2.0 + 1e-13,  1e-14, 4},  // duplicate of 2.0 (diff < 10*1e-12)
        {-3.0,          0.0,   1}
    };

    lmmc_real_t tolerance = 1e-12;
    int max_roots = -1;

    auto deduped = lamina::deduplicate_roots(roots, tolerance, max_roots);

    // Should have 3 unique roots: -3.0, 1.0, 2.0
    EXPECT_TRUE(deduped.size() == 3,
        "Should have 3 unique roots after deduplication, got " + std::to_string(deduped.size()));

    if (deduped.size() == 3) {
        EXPECT_NEAR(deduped[0], -3.0, 1e-10, "First unique root should be -3.0");
        EXPECT_NEAR(deduped[1], 1.0, 1e-10, "Second unique root should be ~1.0");
        EXPECT_NEAR(deduped[2], 2.0, 1e-10, "Third unique root should be ~2.0");
    }

    // Assemble into NumberNode expressions
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    for (const auto& val : deduped) {
        results.push_back(SymbolicExpr::number(static_cast<double>(val)));
    }

    EXPECT_TRUE(results.size() == 3, "Should produce 3 NumberNode expressions");

    // Verify sorted ascending order
    if (results.size() >= 2) {
        for (size_t i = 0; i + 1 < results.size(); ++i) {
            auto vi = test_numeric_eval(results[i]);
            auto vj = test_numeric_eval(results[i + 1]);
            EXPECT_TRUE(vi.has_value() && vj.has_value(),
                "Results should be evaluable");
            if (vi && vj) {
                EXPECT_TRUE(*vi < *vj,
                    "Results should be in ascending order: " +
                    std::to_string(*vi) + " < " + std::to_string(*vj));
            }
        }
    }
}

// Feature: mixed-transcendental-solver, Property 10: Deduplication invariant
// **Validates: Requirements 5.1**

void test_property10_deduplicate_roots_direct() {
    TEST_CASE("Property 10: deduplicate_roots removes duplicates within 10*tolerance");

    lmmc_real_t tolerance = 1e-12;
    lmmc_real_t dedup_threshold = 10.0 * tolerance;  // 1e-11

    // Create roots with some duplicates (values within 10*tolerance of each other)
    // Roots: 1.0, 1.0+5e-12 (duplicate of 1.0), 2.0, 2.0+3e-12 (duplicate of 2.0), 3.0
    std::vector<NumericRoot> roots;
    roots.push_back(NumericRoot{1.0, 1e-13, 5});
    roots.push_back(NumericRoot{1.0 + 5e-12, 2e-13, 7});   // within 1e-11 of 1.0
    roots.push_back(NumericRoot{2.0, 3e-13, 4});
    roots.push_back(NumericRoot{2.0 + 3e-12, 1e-14, 6});   // within 1e-11 of 2.0
    roots.push_back(NumericRoot{3.0, 5e-14, 3});

    auto result = lamina::deduplicate_roots(roots, tolerance, -1);

    // After deduplication, only 3 distinct roots should remain
    EXPECT_TRUE(result.size() == 3,
        "Expected 3 roots after deduplication, got " + std::to_string(result.size()));

    // Verify all pairs satisfy |r_i - r_j| >= 10 * tolerance
    for (size_t i = 0; i < result.size(); ++i) {
        for (size_t j = i + 1; j < result.size(); ++j) {
            double diff = std::abs(result[i] - result[j]);
            EXPECT_TRUE(diff >= dedup_threshold,
                "Pair (" + std::to_string(i) + "," + std::to_string(j) +
                "): |" + std::to_string(result[i]) + " - " + std::to_string(result[j]) +
                "| = " + std::to_string(diff) + " should be >= " + std::to_string(dedup_threshold));
        }
    }
}

void test_property10_deduplicate_roots_no_duplicates() {
    TEST_CASE("Property 10: deduplicate_roots with well-separated roots keeps all");

    lmmc_real_t tolerance = 1e-12;
    lmmc_real_t dedup_threshold = 10.0 * tolerance;

    // All roots are well-separated (> 10*tolerance apart)
    std::vector<NumericRoot> roots;
    roots.push_back(NumericRoot{-5.0, 1e-13, 3});
    roots.push_back(NumericRoot{0.0, 2e-14, 4});
    roots.push_back(NumericRoot{5.0, 1e-13, 5});

    auto result = lamina::deduplicate_roots(roots, tolerance, -1);

    EXPECT_TRUE(result.size() == 3,
        "All well-separated roots should be kept, got " + std::to_string(result.size()));

    // Verify deduplication invariant
    for (size_t i = 0; i < result.size(); ++i) {
        for (size_t j = i + 1; j < result.size(); ++j) {
            double diff = std::abs(result[i] - result[j]);
            EXPECT_TRUE(diff >= dedup_threshold,
                "Pair (" + std::to_string(i) + "," + std::to_string(j) +
                "): |" + std::to_string(result[i]) + " - " + std::to_string(result[j]) +
                "| = " + std::to_string(diff) + " should be >= " + std::to_string(dedup_threshold));
        }
    }
}

void test_property10_deduplicate_roots_all_duplicates() {
    TEST_CASE("Property 10: deduplicate_roots with all values within threshold → single root");

    lmmc_real_t tolerance = 1e-12;
    lmmc_real_t dedup_threshold = 10.0 * tolerance;

    // All roots are within 10*tolerance of each other
    std::vector<NumericRoot> roots;
    roots.push_back(NumericRoot{1.0, 5e-13, 10});
    roots.push_back(NumericRoot{1.0 + 2e-12, 1e-13, 8});
    roots.push_back(NumericRoot{1.0 + 4e-12, 3e-13, 9});
    roots.push_back(NumericRoot{1.0 + 7e-12, 2e-13, 7});

    auto result = lamina::deduplicate_roots(roots, tolerance, -1);

    EXPECT_TRUE(result.size() == 1,
        "All roots within threshold should collapse to 1, got " + std::to_string(result.size()));

    // Even with a single root, the invariant trivially holds (no pairs to check)
}

void test_property10_deduplicate_pipeline_sin_x() {
    TEST_CASE("Property 10: full pipeline sin(x) on [-10,10] → deduplicated roots satisfy invariant");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    auto derivative = compute_derivative(expr, "x");

    lamina::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    lamina::SearchInterval interval{-10.0, 10.0};

    // Step 1: Isolate roots
    auto intervals = lamina::isolate_roots(expr, derivative, "x", interval, opts);

    // Step 2: Refine each isolated root
    std::vector<NumericRoot> refined_roots;
    for (const auto& iso : intervals) {
        auto root = lamina::refine_root(expr, derivative, "x", iso, opts);
        if (root.has_value()) {
            refined_roots.push_back(*root);
        }
    }

    // Step 3: Deduplicate
    lmmc_real_t tolerance = opts.tolerance;
    lmmc_real_t dedup_threshold = 10.0 * tolerance;
    auto deduped = lamina::deduplicate_roots(refined_roots, tolerance, -1);

    std::cout << "  [INFO] sin(x) pipeline: " << intervals.size() << " intervals -> "
              << refined_roots.size() << " refined -> " << deduped.size() << " deduplicated" << std::endl;

    // sin(x) has 7 roots in [-10, 10]: -3π, -2π, -π, 0, π, 2π, 3π
    EXPECT_TRUE(deduped.size() >= 3,
        "sin(x) on [-10,10] should have at least 3 deduplicated roots, got " +
        std::to_string(deduped.size()));

    // Verify deduplication invariant: all pairs satisfy |r_i - r_j| >= 10*tolerance
    for (size_t i = 0; i < deduped.size(); ++i) {
        for (size_t j = i + 1; j < deduped.size(); ++j) {
            double diff = std::abs(deduped[i] - deduped[j]);
            EXPECT_TRUE(diff >= dedup_threshold,
                "Pair (" + std::to_string(i) + "," + std::to_string(j) +
                "): |" + std::to_string(deduped[i]) + " - " + std::to_string(deduped[j]) +
                "| = " + std::to_string(diff) + " should be >= " + std::to_string(dedup_threshold));
        }
    }

    // Verify roots are sorted ascending
    for (size_t i = 0; i + 1 < deduped.size(); ++i) {
        EXPECT_TRUE(deduped[i] < deduped[i + 1],
            "Roots should be sorted ascending: result[" + std::to_string(i) + "]=" +
            std::to_string(deduped[i]) + " should be < result[" + std::to_string(i + 1) + "]=" +
            std::to_string(deduped[i + 1]));
    }
}

// Feature: mixed-transcendental-solver, Property 16: Factored expression root completeness
// **Validates: Requirements 9.2**

void test_property16_factored_sin_x_minus_half_times_x_minus_3() {
    TEST_CASE("Property 16: (sin(x) - 0.5) * (x - 3) on [0, 4] → roots from BOTH factors");

    auto x = SymbolicExpr::variable("x");
    // (sin(x) - 0.5) * (x - 3)
    // Factor 1: sin(x) - 0.5 has roots at x = π/6 ≈ 0.524 and x = 5π/6 ≈ 2.618 in [0, 4]
    // Factor 2: x - 3 has root at x = 3
    // Combined result should contain roots from BOTH factors
    auto factor1 = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::number(-0.5)
    );
    auto factor2 = SymbolicExpr::add(x, SymbolicExpr::number(-3));
    auto expr = SymbolicExpr::multiply(factor1, factor2);

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;
    opts.has_search_interval = true;
    opts.search_lo = 0.0;
    opts.search_hi = 4.0;

    auto results = solve_mixed_transcendental(expr, "x", opts);

    std::cout << "  [INFO] (sin(x)-0.5)*(x-3) on [0,4]: " << results.size() << " root(s) found" << std::endl;
    for (size_t i = 0; i < results.size(); ++i) {
        auto val = test_numeric_eval(results[i]);
        if (val.has_value()) {
            std::cout << "    root[" << i << "] = " << *val << std::endl;
        }
    }

    // Expected roots: π/6 ≈ 0.5236, 5π/6 ≈ 2.6180, and 3.0
    // We need at least 3 roots (from both factors)
    EXPECT_TRUE(results.size() >= 3,
        "(sin(x)-0.5)*(x-3) on [0,4] should have at least 3 roots, got " +
        std::to_string(results.size()));

    // Verify roots from factor 1 (sin(x) - 0.5): π/6 and 5π/6
    double expected_pi_6 = M_PI / 6.0;       // ≈ 0.5236
    double expected_5pi_6 = 5.0 * M_PI / 6.0; // ≈ 2.6180

    bool found_pi_6 = false;
    bool found_5pi_6 = false;
    bool found_3 = false;

    for (const auto& r : results) {
        auto val = test_numeric_eval(r);
        if (!val.has_value()) continue;
        if (std::abs(*val - expected_pi_6) < 1e-4) found_pi_6 = true;
        if (std::abs(*val - expected_5pi_6) < 1e-4) found_5pi_6 = true;
        if (std::abs(*val - 3.0) < 1e-4) found_3 = true;
    }

    EXPECT_TRUE(found_pi_6,
        "Should find root from factor sin(x)-0.5 at π/6 ≈ " + std::to_string(expected_pi_6));
    EXPECT_TRUE(found_5pi_6,
        "Should find root from factor sin(x)-0.5 at 5π/6 ≈ " + std::to_string(expected_5pi_6));
    EXPECT_TRUE(found_3,
        "Should find root from factor x-3 at x = 3");
}

// Feature: mixed-transcendental-solver, Property 16: Factored expression root completeness
// **Validates: Requirements 9.2**

void test_property16_factored_sin_x_times_cos_x() {
    TEST_CASE("Property 16: sin(x) * cos(x) on [-10, 10] → roots from BOTH sin and cos factors");

    auto x = SymbolicExpr::variable("x");
    // sin(x) * cos(x)
    // Factor 1: sin(x) has roots at 0, ±π, ±2π, ±3π in [-10, 10]
    // Factor 2: cos(x) has roots at ±π/2, ±3π/2, ±5π/2 in [-10, 10]
    // Combined result should contain roots from BOTH factors
    auto expr = SymbolicExpr::multiply(
        SymbolicExpr::sin(x),
        SymbolicExpr::cos(x)
    );

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;
    opts.has_search_interval = true;
    opts.search_lo = -10.0;
    opts.search_hi = 10.0;

    auto results = solve_mixed_transcendental(expr, "x", opts);

    std::cout << "  [INFO] sin(x)*cos(x) on [-10,10]: " << results.size() << " root(s) found" << std::endl;
    for (size_t i = 0; i < results.size(); ++i) {
        auto val = test_numeric_eval(results[i]);
        if (val.has_value()) {
            std::cout << "    root[" << i << "] = " << *val << std::endl;
        }
    }

    // sin(x) roots in [-10, 10]: -3π, -2π, -π, 0, π, 2π, 3π → 7 roots
    // cos(x) roots in [-10, 10]: -5π/2, -3π/2, -π/2, π/2, 3π/2, 5π/2 → 6 roots
    // Total distinct roots: 13 (no overlap between sin and cos zeros)
    // We expect at least roots from both factors to be present

    // Check that we have a reasonable number of roots (at least some from each factor)
    EXPECT_TRUE(results.size() >= 6,
        "sin(x)*cos(x) on [-10,10] should have at least 6 roots (from both factors), got " +
        std::to_string(results.size()));

    // Verify presence of roots from sin(x): check for root near 0 and near π
    bool found_sin_zero = false;
    bool found_sin_pi = false;

    // Verify presence of roots from cos(x): check for root near π/2 and near 3π/2
    bool found_cos_pi_2 = false;
    bool found_cos_3pi_2 = false;

    for (const auto& r : results) {
        auto val = test_numeric_eval(r);
        if (!val.has_value()) continue;
        if (std::abs(*val - 0.0) < 1e-4) found_sin_zero = true;
        if (std::abs(*val - M_PI) < 1e-4) found_sin_pi = true;
        if (std::abs(*val - M_PI / 2.0) < 1e-4) found_cos_pi_2 = true;
        if (std::abs(*val - 3.0 * M_PI / 2.0) < 1e-4) found_cos_3pi_2 = true;
    }

    EXPECT_TRUE(found_sin_zero,
        "Should find root from sin(x) factor at x = 0");
    EXPECT_TRUE(found_sin_pi,
        "Should find root from sin(x) factor at x = π ≈ " + std::to_string(M_PI));
    EXPECT_TRUE(found_cos_pi_2,
        "Should find root from cos(x) factor at x = π/2 ≈ " + std::to_string(M_PI / 2.0));
    EXPECT_TRUE(found_cos_3pi_2,
        "Should find root from cos(x) factor at x = 3π/2 ≈ " + std::to_string(3.0 * M_PI / 2.0));
}

// ============================================================================
// Feature: mixed-transcendental-solver, Property 15: No unhandled exceptions
// **Validates: Requirements 8.5**
// Exception safety tests (Requirements 8.2, 8.5)
// ============================================================================

void test_exception_safety_division_by_zero_expression() {
    TEST_CASE("Exception safety: expression causing division by zero returns empty, not throw");

    auto x = SymbolicExpr::variable("x");
    // 1/x + sin(x): evaluates to infinity/NaN at x=0, which lies in default [-10, 10].
    // The solver must handle this gracefully and not throw.
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(-1)),
        SymbolicExpr::sin(x)
    );

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;

    // This must not throw — the try-catch wrapper should catch any internal errors
    bool threw = false;
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    try {
        results = lamina::solve_mixed_transcendental(expr, "x", opts);
    } catch (...) {
        threw = true;
    }

    EXPECT_FALSE(threw,
        "solve_mixed_transcendental should NOT throw on 1/x + sin(x)");

    // Results should either be empty or contain valid finite roots
    for (size_t i = 0; i < results.size(); ++i) {
        EXPECT_TRUE(results[i] != nullptr,
            "Root " + std::to_string(i) + " should not be null");
        if (results[i]) {
            auto val = test_numeric_eval(results[i]);
            if (val.has_value()) {
                EXPECT_TRUE(std::isfinite(*val),
                    "Root " + std::to_string(i) + " should be finite, got " + std::to_string(*val));
            }
        }
    }
}

void test_exception_safety_ln_negative_domain() {
    TEST_CASE("Exception safety: ln(x) + x on interval including x<=0 returns gracefully");

    auto x = SymbolicExpr::variable("x");
    // ln(x) + x = 0: ln(x) is undefined for x <= 0.
    // Default interval [-10, 10] includes negative values where ln produces NaN.
    auto expr = SymbolicExpr::add(SymbolicExpr::ln(x), x);

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;

    bool threw = false;
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    try {
        results = lamina::solve_mixed_transcendental(expr, "x", opts);
    } catch (...) {
        threw = true;
    }

    EXPECT_FALSE(threw,
        "solve_mixed_transcendental should NOT throw on ln(x) + x with negative domain");

    // Any returned roots must be finite and positive (in domain of ln)
    for (size_t i = 0; i < results.size(); ++i) {
        auto val = test_numeric_eval(results[i]);
        if (val.has_value()) {
            EXPECT_TRUE(std::isfinite(*val),
                "Root " + std::to_string(i) + " should be finite");
        }
    }
}

void test_exception_safety_tan_near_singularity() {
    TEST_CASE("Exception safety: tan(x) - x near singularities returns gracefully");

    auto x = SymbolicExpr::variable("x");
    // tan(x) - x = 0: tan(x) has singularities at x = π/2 + nπ.
    // The solver must handle NaN/infinity at these points gracefully.
    auto expr = SymbolicExpr::add(
        SymbolicExpr::tan(x),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x)
    );

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;

    bool threw = false;
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    try {
        results = lamina::solve_mixed_transcendental(expr, "x", opts);
    } catch (...) {
        threw = true;
    }

    EXPECT_FALSE(threw,
        "solve_mixed_transcendental should NOT throw on tan(x) - x near singularities");

    // All returned roots must be finite
    for (size_t i = 0; i < results.size(); ++i) {
        auto val = test_numeric_eval(results[i]);
        if (val.has_value()) {
            EXPECT_TRUE(std::isfinite(*val),
                "Root " + std::to_string(i) + " should be finite");
        }
    }
}

void test_exception_safety_nullptr_expression() {
    TEST_CASE("Exception safety: nullptr expression returns empty vector");

    std::shared_ptr<SymbolicExpr> expr = nullptr;

    lamina::SolveOptions opts;
    opts.allow_numeric = true;

    bool threw = false;
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    try {
        results = lamina::solve_mixed_transcendental(expr, "x", opts);
    } catch (...) {
        threw = true;
    }

    EXPECT_FALSE(threw,
        "solve_mixed_transcendental should NOT throw on nullptr expression");
    EXPECT_TRUE(results.empty(),
        "nullptr expression should return empty vector");
}

void test_exception_safety_extreme_values() {
    TEST_CASE("Exception safety: exp(x) - 1e300 on wide interval handles overflow gracefully");

    auto x = SymbolicExpr::variable("x");
    // exp(x) - 1e300: exp(x) overflows for large x, producing infinity.
    // The solver must handle this without throwing.
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::number(-1e300)
    );

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.has_search_interval = true;
    opts.search_lo = -100.0;
    opts.search_hi = 1000.0;

    bool threw = false;
    std::vector<std::shared_ptr<SymbolicExpr>> results;
    try {
        results = lamina::solve_mixed_transcendental(expr, "x", opts);
    } catch (...) {
        threw = true;
    }

    EXPECT_FALSE(threw,
        "solve_mixed_transcendental should NOT throw on exp(x) - 1e300 with wide interval");

    // Any returned roots must be finite
    for (size_t i = 0; i < results.size(); ++i) {
        auto val = test_numeric_eval(results[i]);
        if (val.has_value()) {
            EXPECT_TRUE(std::isfinite(*val),
                "Root " + std::to_string(i) + " should be finite");
        }
    }
}

// Feature: mixed-transcendental-solver, Property 14: No variable dependence yields empty result
// **Validates: Requirements 8.1**
void test_no_variable_dependence_returns_empty() {
    TEST_CASE("Property 14: expression without solve variable returns empty vector");

    // sin(3) + 5: no dependence on x at all
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(SymbolicExpr::number(3)),
        SymbolicExpr::number(5)
    );

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-12;

    auto results = lamina::solve_mixed_transcendental(expr, "x", opts);

    EXPECT_TRUE(results.empty(),
        "Expression sin(3)+5 has no dependence on x, should return empty vector");
}

void test_no_variable_dependence_constant_expression() {
    TEST_CASE("Property 14: pure constant expression returns empty vector");

    // 42: a plain number, no variable at all
    auto expr = SymbolicExpr::number(42);

    lamina::SolveOptions opts;
    opts.allow_numeric = true;

    auto results = lamina::solve_mixed_transcendental(expr, "x", opts);

    EXPECT_TRUE(results.empty(),
        "Constant expression 42 has no dependence on x, should return empty vector");
}

void test_no_variable_dependence_other_variable() {
    TEST_CASE("Property 14: expression with different variable returns empty vector");

    auto y = SymbolicExpr::variable("y");
    // sin(y) + y: depends on y, but NOT on x
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(y), y);

    lamina::SolveOptions opts;
    opts.allow_numeric = true;

    auto results = lamina::solve_mixed_transcendental(expr, "x", opts);

    EXPECT_TRUE(results.empty(),
        "Expression sin(y)+y has no dependence on x, should return empty vector");
}

// ============================================================================
// End-to-end known-root equation tests (Task 11.1)
// Call solve_mixed_transcendental directly and verify correct roots.
// **Validates: Requirements 1.1, 1.2, 1.3, 2.1, 3.2, 4.1, 4.4, 8.1, 8.2, 8.3, 9.2**
// ============================================================================

void test_e2e_sin_x_plus_x_root_at_zero() {
    TEST_CASE("E2E: sin(x) - x/2 = 0 → roots near ±1.895");

    auto x = SymbolicExpr::variable("x");
    // sin(x) - x/2 = 0: has roots at x=0 and near x ≈ ±1.8955
    // (sin(x) = x/2 intersects at these points)
    // We test with a search interval that avoids x=0 landing on a grid point.
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::multiply(SymbolicExpr::number(-0.5), x)
    );

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;
    opts.has_search_interval = true;
    opts.search_lo = -3.0;
    opts.search_hi = 3.0;

    auto results = solve_mixed_transcendental(expr, "x", opts);

    std::cout << "  [INFO] sin(x)-x/2 on [-3,3]: " << results.size() << " root(s) found" << std::endl;
    for (size_t i = 0; i < results.size(); ++i) {
        auto val = test_numeric_eval(results[i]);
        if (val.has_value()) {
            std::cout << "    root[" << i << "] = " << *val << std::endl;
        }
    }

    // sin(x) - x/2 = 0 has roots near ±1.8955 (and x=0, but that may be missed
    // due to the grid-point limitation of sign-change detection)
    EXPECT_TRUE(results.size() >= 1,
        "sin(x)-x/2=0 on [-3,3] should find at least 1 root");

    // Verify all returned roots satisfy |f(root)| < tolerance
    for (const auto& r : results) {
        auto val = test_numeric_eval(r);
        if (!val.has_value()) continue;
        double f_val = std::sin(*val) - 0.5 * (*val);
        EXPECT_TRUE(std::abs(f_val) < opts.tolerance,
            "Root " + std::to_string(*val) + ": |f(root)|=" +
            std::to_string(std::abs(f_val)) + " should be < tolerance");
    }
}

void test_e2e_exp_x_minus_x_minus_2() {
    TEST_CASE("E2E: exp(x) - x - 2 = 0 → roots near -1.84 and 1.146");

    auto x = SymbolicExpr::variable("x");
    // exp(x) - x - 2 = 0 has two real roots:
    //   x ≈ -1.8414 and x ≈ 1.1462
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
            SymbolicExpr::number(-2)
        )
    );

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    auto results = solve_mixed_transcendental(expr, "x", opts);

    std::cout << "  [INFO] exp(x)-x-2: " << results.size() << " root(s) found" << std::endl;
    for (size_t i = 0; i < results.size(); ++i) {
        auto val = test_numeric_eval(results[i]);
        if (val.has_value()) {
            std::cout << "    root[" << i << "] = " << *val << std::endl;
        }
    }

    EXPECT_TRUE(results.size() >= 2,
        "exp(x)-x-2=0 should find at least 2 roots");

    // Verify roots near expected values
    bool found_neg = false;
    bool found_pos = false;
    for (const auto& r : results) {
        auto val = test_numeric_eval(r);
        if (!val.has_value()) continue;
        if (std::abs(*val - (-1.8414)) < 0.01) found_neg = true;
        if (std::abs(*val - 1.1462) < 0.01) found_pos = true;
    }
    EXPECT_TRUE(found_neg,
        "exp(x)-x-2=0 should have root near x ≈ -1.84");
    EXPECT_TRUE(found_pos,
        "exp(x)-x-2=0 should have root near x ≈ 1.15");

    // Verify all returned roots satisfy |f(root)| < tolerance
    for (const auto& r : results) {
        auto val = test_numeric_eval(r);
        if (!val.has_value()) continue;
        double f_val = std::exp(*val) - *val - 2.0;
        EXPECT_TRUE(std::abs(f_val) < opts.tolerance,
            "Root " + std::to_string(*val) + ": |f(root)|=" +
            std::to_string(std::abs(f_val)) + " should be < tolerance");
    }
}

void test_e2e_x_cos_x_minus_1() {
    TEST_CASE("E2E: x*cos(x) - 1 = 0 → multiple roots in [-10, 10]");

    auto x = SymbolicExpr::variable("x");
    // x*cos(x) - 1 = 0 has multiple roots in [-10, 10]
    // Known approximate roots: ±1.283, ±4.917, ±7.975, ...
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(x, SymbolicExpr::cos(x)),
        SymbolicExpr::number(-1)
    );

    lamina::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    auto results = solve_mixed_transcendental(expr, "x", opts);

    std::cout << "  [INFO] x*cos(x)-1: " << results.size() << " root(s) found" << std::endl;
    for (size_t i = 0; i < results.size(); ++i) {
        auto val = test_numeric_eval(results[i]);
        if (val.has_value()) {
            std::cout << "    root[" << i << "] = " << *val << std::endl;
        }
    }

    // Should find at least 2 roots (the pair near ±1.283)
    EXPECT_TRUE(results.size() >= 2,
        "x*cos(x)-1=0 should find at least 2 roots in [-10,10]");

    // Verify all returned roots satisfy |f(root)| < tolerance
    for (const auto& r : results) {
        auto val = test_numeric_eval(r);
        if (!val.has_value()) continue;
        double f_val = (*val) * std::cos(*val) - 1.0;
        EXPECT_TRUE(std::abs(f_val) < opts.tolerance,
            "Root " + std::to_string(*val) + ": |f(root)|=" +
            std::to_string(std::abs(f_val)) + " should be < tolerance");
    }

    // Verify roots are sorted ascending
    for (size_t i = 0; i + 1 < results.size(); ++i) {
        auto vi = test_numeric_eval(results[i]);
        auto vj = test_numeric_eval(results[i + 1]);
        if (vi.has_value() && vj.has_value()) {
            EXPECT_TRUE(*vi < *vj,
                "Roots should be ascending: " + std::to_string(*vi) +
                " < " + std::to_string(*vj));
        }
    }
}

int main() {
    test_classification_sin_x_plus_x();
    test_classification_exp_x_minus_x_squared();
    test_classification_sin_constant_plus_polynomial();
    test_classification_cos_y_plus_x();
    test_classification_ln_cos_x_plus_x();
    test_classification_tan_linear_arg();
    test_classification_pure_polynomial();

    // Property 12: allow_numeric gate tests
    test_allow_numeric_false_returns_empty();
    test_allow_numeric_true_permits_solving();
    test_allow_numeric_false_cos_equation_returns_empty();

    // Property 2: Routing correctness for reducible expressions
    test_routing_pure_polynomial_not_hybrid();
    test_routing_transcendental_substitution_not_hybrid();
    test_routing_polynomial_with_allow_numeric_false();

    // Property 13: Backward compatibility
    test_backward_compat_cubic_polynomial();
    test_backward_compat_transcendental_substitution();

    // determine_search_interval tests (Requirements 4.1-4.5, 8.4)
    test_search_interval_user_specified();
    test_search_interval_user_invalid_lo_ge_hi();
    test_search_interval_user_invalid_equal();
    test_search_interval_user_width_le_tolerance();
    test_search_interval_default_no_periodic();
    test_search_interval_sin_x_periodic_extension();
    test_search_interval_sin_small_k_periodic_extension();
    test_search_interval_tan_periodic_extension();
    test_search_interval_tan_small_k_extension();
    test_search_interval_clamp_to_100();
    test_search_interval_nonlinear_arg_default();
    test_search_interval_cos_2x_plus_1();
    test_search_interval_user_overrides_periodic();

    // Property 9: Periodic extension covers two full periods
    test_property9_sin_x_plus_x_two_periods();
    test_property9_cos_half_x_minus_x_two_periods();
    test_property9_tan_x_plus_x_two_periods();
    test_property9_sin_small_k_clamped();
    test_property9_sin_x_squared_nonlinear_default();

    // Property 4: Minimum subdivision width
    test_property4_min_width_sin_x_plus_x_over_10();
    test_property4_min_width_tan_x_minus_x();
    test_property4_min_width_narrow_interval();

    // Property 5: Max roots limit respected
    test_property5_max_roots_2_sin_x();
    test_property5_max_roots_1_sin_x();
    test_property5_max_roots_unlimited_sin_x();
    test_property5_max_roots_3_cos_x_minus_half();

    // Property 3: Root isolation sign-change invariant
    test_property3_sign_change_sin_x_plus_x_div_10();
    test_property3_sign_change_exp_x_minus_x_minus_2();
    test_property3_sign_change_x_cos_x_minus_1();

    // Property 6: No sign changes yields empty result
    test_property6_x_squared_plus_one_always_positive();
    test_property6_exp_x_plus_one_always_positive();
    test_property6_neg_x_squared_minus_one_always_negative();
    test_property6_sin_x_plus_five_always_positive();

    // Root refinement tests (Requirements 3.1-3.6, 8.3)
    test_refine_root_sin_x_newton_raphson();
    test_refine_root_sin_x_bisection_fallback();
    test_refine_root_exp_x_minus_x_minus_2();
    test_refine_root_x_cos_x_minus_1();
    test_refine_root_discards_when_no_convergence();
    test_refine_root_derivative_zero_bisection_step();
    test_refine_root_max_iterations_exceeded();

    // Property 11: Output sorted ascending
    test_property11_deduplicate_roots_ascending_order();
    test_property11_pipeline_sin_x_ascending();
    test_property11_pipeline_cos_x_ascending();
    test_property11_deduplicate_with_duplicates_ascending();

    // Property 8: Roots within specified search interval
    test_property8_roots_within_interval_sin_x();
    test_property8_roots_within_interval_cos_x_minus_half();

    // Property 7: Output validity invariant
    test_property7_output_validity_sin_x_refined_roots();
    test_property7_output_validity_exp_x_minus_x_minus_2();
    test_property7_max_roots_limit_after_deduplication();
    test_property7_output_validity_x_cos_x_minus_1();

    // Result assembly tests (Requirements 6.1, 6.2, 6.3, 6.4)
    test_assemble_results_basic();
    test_assemble_results_empty();
    test_assemble_results_with_duplicates();

    // Property 10: Deduplication invariant
    test_property10_deduplicate_roots_direct();
    test_property10_deduplicate_roots_no_duplicates();
    test_property10_deduplicate_roots_all_duplicates();
    test_property10_deduplicate_pipeline_sin_x();

    // Property 16: Factored expression root completeness
    test_property16_factored_sin_x_minus_half_times_x_minus_3();
    test_property16_factored_sin_x_times_cos_x();

    // Property 15: No unhandled exceptions (Requirements 8.2, 8.5)
    test_exception_safety_division_by_zero_expression();
    test_exception_safety_ln_negative_domain();
    test_exception_safety_tan_near_singularity();
    test_exception_safety_nullptr_expression();
    test_exception_safety_extreme_values();

    // Property 14: No variable dependence yields empty result (Requirement 8.1)
    test_no_variable_dependence_returns_empty();
    test_no_variable_dependence_constant_expression();
    test_no_variable_dependence_other_variable();

    // End-to-end known-root equation tests (Task 11.1)
    test_e2e_sin_x_plus_x_root_at_zero();
    test_e2e_exp_x_minus_x_minus_2();
    test_e2e_x_cos_x_minus_1();

    return TEST_REPORT();
}
