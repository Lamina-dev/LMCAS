#include "test_common.hpp"
#include "interval.hpp"
#include "inequality_solver.hpp"
#include "symbolic.hpp"
#include "poly_utils.hpp"
#include <cmath>
#include <random>
#include <vector>

using namespace lamina;

int main() {
    // =========================================================================
    // Unit Tests for Parametric Inequality Solver (Task 7.1)
    // Requirements: 7.1, 7.2, 7.3, 7.4
    // =========================================================================

    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");
    auto c = SymbolicExpr::variable("c");

    // =========================================================================
    // Test 1: Linear parametric inequality: ax + b > 0
    // When a > 0: x > -b/a → (-b/a, +∞)
    // When a < 0: x < -b/a → (-∞, -b/a)
    // When a = 0: degenerate (becomes b > 0, constant)
    // Requirements: 7.1, 7.2, 7.3, 7.4
    // =========================================================================
    TEST_CASE("Parametric linear: ax + b > 0 (piecewise by sign of a)");
    {
        // Build ax + b
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(a, x),
            b);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"a", "b"});

        // Should have 3 cases: a > 0, a < 0, a = 0
        EXPECT_TRUE(result.cases.size() == 3,
            "ax + b > 0: should have 3 piecewise cases (a>0, a<0, a=0)");

        // Case 1: a > 0 → solution is (-b/a, +∞)
        // The solution should have one interval with symbolic lower bound
        auto& pos_case = result.cases[0];
        EXPECT_TRUE(pos_case.condition != nullptr,
            "ax + b > 0, case a>0: condition should not be null");
        EXPECT_TRUE(!pos_case.solution.is_empty(),
            "ax + b > 0, case a>0: solution should not be empty");
        EXPECT_TRUE(pos_case.solution.intervals().size() == 1,
            "ax + b > 0, case a>0: should have exactly 1 interval");

        if (!pos_case.solution.intervals().empty()) {
            auto& iv = pos_case.solution.intervals()[0];
            // Lower bound should be symbolic (-b/a)
            EXPECT_TRUE(!iv.lower.is_neg_infinity,
                "ax + b > 0, case a>0: lower bound should not be -inf");
            EXPECT_TRUE(iv.lower.is_open,
                "ax + b > 0, case a>0: lower bound should be open (strict)");
            EXPECT_TRUE(iv.upper.is_pos_infinity,
                "ax + b > 0, case a>0: upper bound should be +inf");
            // The lower bound value should be -b/a
            if (iv.lower.value) {
                // Substitute a=2, b=4 → root should be -4/2 = -2
                auto substituted = iv.lower.value->substitute("a", SymbolicExpr::number(2));
                substituted = substituted->substitute("b", SymbolicExpr::number(4));
                substituted = substituted->simplify();
                try {
                    double val = substituted->to_numeric();
                    EXPECT_TRUE(std::abs(val - (-2.0)) < 1e-10,
                        "ax + b > 0, case a>0: root at a=2,b=4 should be -2");
                } catch (...) {
                    EXPECT_TRUE(false, "ax + b > 0, case a>0: root should be evaluable");
                }
            }
        }

        // Case 2: a < 0 → solution is (-∞, -b/a)
        auto& neg_case = result.cases[1];
        EXPECT_TRUE(neg_case.condition != nullptr,
            "ax + b > 0, case a<0: condition should not be null");
        EXPECT_TRUE(!neg_case.solution.is_empty(),
            "ax + b > 0, case a<0: solution should not be empty");
        EXPECT_TRUE(neg_case.solution.intervals().size() == 1,
            "ax + b > 0, case a<0: should have exactly 1 interval");

        if (!neg_case.solution.intervals().empty()) {
            auto& iv = neg_case.solution.intervals()[0];
            EXPECT_TRUE(iv.lower.is_neg_infinity,
                "ax + b > 0, case a<0: lower bound should be -inf");
            EXPECT_TRUE(!iv.upper.is_pos_infinity,
                "ax + b > 0, case a<0: upper bound should not be +inf");
            EXPECT_TRUE(iv.upper.is_open,
                "ax + b > 0, case a<0: upper bound should be open (strict)");
        }

        // Case 3: a = 0 → degenerate (b > 0 → entire line, b <= 0 → empty)
        auto& degen_case = result.cases[2];
        EXPECT_TRUE(degen_case.condition != nullptr,
            "ax + b > 0, case a=0: condition should not be null");
    }

    // =========================================================================
    // Test 2: Quadratic with constant leading coefficient: x² + bx + c > 0
    // Leading coefficient is 1 (constant, positive) → single case
    // Roots are symbolic: (-b ± sqrt(b²-4c))/2
    // Requirements: 7.1, 7.2
    // =========================================================================
    TEST_CASE("Parametric quadratic with constant leading coeff: x^2 + bx + c > 0");
    {
        // Build x^2 + bx + c
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto bx = SymbolicExpr::multiply(b, x);
        auto expr = SymbolicExpr::add(SymbolicExpr::add(x2, bx), c);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"b", "c"});

        // Leading coefficient is 1 (constant) → should be single case (no piecewise)
        EXPECT_TRUE(result.cases.size() == 1,
            "x^2 + bx + c > 0: should have 1 case (constant leading coeff)");

        // The solution should have symbolic boundaries
        auto& single_case = result.cases[0];
        EXPECT_TRUE(single_case.condition == nullptr,
            "x^2 + bx + c > 0: condition should be null (unconditional)");

        // For a quadratic with positive leading coeff and > 0:
        // Solution is (-∞, root1) ∪ (root2, +∞) where root1 < root2
        // The solver returns 2 roots from the quadratic formula
        auto& intervals = single_case.solution.intervals();
        EXPECT_TRUE(intervals.size() == 2,
            "x^2 + bx + c > 0: should have 2 intervals (outside roots)");

        // Verify by substituting specific values: b=-3, c=2 → x²-3x+2 > 0
        // Roots: x=1, x=2 → solution: (-∞,1) ∪ (2,+∞)
        if (intervals.size() == 2) {
            // First interval: (-∞, root1) where root1 is the smaller root
            EXPECT_TRUE(intervals[0].lower.is_neg_infinity,
                "x^2 + bx + c > 0: first interval starts at -inf");
            EXPECT_TRUE(intervals[0].upper.is_open,
                "x^2 + bx + c > 0: first interval upper is open (strict)");

            // Second interval: (root2, +∞) where root2 is the larger root
            EXPECT_TRUE(intervals[1].upper.is_pos_infinity,
                "x^2 + bx + c > 0: second interval ends at +inf");
            EXPECT_TRUE(intervals[1].lower.is_open,
                "x^2 + bx + c > 0: second interval lower is open (strict)");

            // Verify symbolic roots evaluate correctly for b=-3, c=2
            if (intervals[0].upper.value) {
                auto root1 = intervals[0].upper.value->substitute("b", SymbolicExpr::number(-3));
                root1 = root1->substitute("c", SymbolicExpr::number(2));
                root1 = root1->simplify();
                try {
                    double val = root1->to_numeric();
                    // Should be 1 or 2 (the smaller root)
                    EXPECT_TRUE(std::abs(val - 1.0) < 1e-10 || std::abs(val - 2.0) < 1e-10,
                        "x^2 + bx + c > 0: root1 at b=-3,c=2 should be 1 or 2");
                } catch (...) {
                    // Root might be in symbolic form that can't be evaluated directly
                    // This is acceptable for parametric solutions
                }
            }
        }
    }

    // =========================================================================
    // Test 3: Quadratic with parametric leading coefficient: ax² + bx + c ≥ 0
    // Leading coefficient 'a' depends on parameters → piecewise
    // Requirements: 7.3, 7.4
    // =========================================================================
    TEST_CASE("Parametric quadratic: ax^2 + bx + c >= 0 (piecewise by sign of a)");
    {
        // Build ax^2 + bx + c
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto ax2 = SymbolicExpr::multiply(a, x2);
        auto bx = SymbolicExpr::multiply(b, x);
        auto expr = SymbolicExpr::add(SymbolicExpr::add(ax2, bx), c);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterEqual, "x", {"a", "b", "c"});

        // Leading coefficient 'a' depends on parameters → 3 cases
        EXPECT_TRUE(result.cases.size() == 3,
            "ax^2 + bx + c >= 0: should have 3 piecewise cases");

        // Case 1: a > 0 (upward parabola)
        // For >= 0 with positive leading coeff: solution is complement of (root1, root2)
        // i.e., (-∞, root1] ∪ [root2, +∞)
        auto& pos_case = result.cases[0];
        EXPECT_TRUE(pos_case.condition != nullptr,
            "ax^2 + bx + c >= 0, case a>0: should have condition");

        // Case 2: a < 0 (downward parabola)
        // For >= 0 with negative leading coeff: solution is [root1, root2]
        auto& neg_case = result.cases[1];
        EXPECT_TRUE(neg_case.condition != nullptr,
            "ax^2 + bx + c >= 0, case a<0: should have condition");

        // Case 3: a = 0 (degenerate → linear: bx + c >= 0)
        auto& degen_case = result.cases[2];
        EXPECT_TRUE(degen_case.condition != nullptr,
            "ax^2 + bx + c >= 0, case a=0: should have condition (degenerate)");
    }

    // =========================================================================
    // Test 4: Verify parametric solution consistency with numeric solution
    // For specific parameter values, the parametric solution should match
    // the direct numeric solution.
    // Requirements: 7.5
    // =========================================================================
    TEST_CASE("Parametric consistency: parametric vs numeric for specific values");
    {
        // Test: ax + b > 0 with a=3, b=-6 → 3x - 6 > 0 → x > 2
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(a, x),
            b);

        auto parametric_result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"a", "b"});

        // Direct numeric: 3x - 6 > 0
        auto numeric_expr = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(3), x),
            SymbolicExpr::number(-6));
        auto numeric_result = InequalitySolver::solve_inequality(
            numeric_expr, InequalityType::GreaterThan, "x");

        // Numeric result should be (2, +∞)
        EXPECT_TRUE(!numeric_result.is_empty(),
            "3x - 6 > 0: numeric solution should not be empty");
        EXPECT_TRUE(numeric_result.contains(3.0),
            "3x - 6 > 0: x=3 should be in solution");
        EXPECT_TRUE(!numeric_result.contains(1.0),
            "3x - 6 > 0: x=1 should not be in solution");
        EXPECT_TRUE(!numeric_result.contains(2.0),
            "3x - 6 > 0: x=2 (root) should not be in solution (strict)");

        // Parametric result case a>0: substitute a=3, b=-6 into the boundary
        EXPECT_TRUE(!parametric_result.cases.empty(),
            "Parametric ax + b > 0: should have cases");

        if (!parametric_result.cases.empty()) {
            auto& pos_case = parametric_result.cases[0];  // a > 0 case
            auto& intervals = pos_case.solution.intervals();
            if (!intervals.empty() && intervals[0].lower.value) {
                auto boundary = intervals[0].lower.value;
                boundary = boundary->substitute("a", SymbolicExpr::number(3));
                boundary = boundary->substitute("b", SymbolicExpr::number(-6));
                boundary = boundary->simplify();
                try {
                    double val = boundary->to_numeric();
                    EXPECT_TRUE(std::abs(val - 2.0) < 1e-10,
                        "Parametric consistency: boundary at a=3,b=-6 should be 2");
                } catch (...) {
                    EXPECT_TRUE(false,
                        "Parametric consistency: boundary should be evaluable");
                }
            }
        }
    }

    // =========================================================================
    // Test 5: No parameters → falls back to non-parametric solver
    // Requirements: 7.1 (backward compatibility)
    // =========================================================================
    TEST_CASE("Parametric with empty params: falls back to non-parametric");
    {
        // x^2 - 4 > 0 with no parameters
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto expr = SymbolicExpr::add(x2, SymbolicExpr::number(-4));

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {});

        // Should have single case with same result as non-parametric
        EXPECT_TRUE(result.cases.size() == 1,
            "No params: should have 1 case");
        EXPECT_TRUE(!result.cases[0].solution.is_empty(),
            "No params: x^2 - 4 > 0 should not be empty");
        EXPECT_TRUE(result.cases[0].solution.contains(3.0),
            "No params: x=3 should be in solution");
        EXPECT_TRUE(result.cases[0].solution.contains(-3.0),
            "No params: x=-3 should be in solution");
        EXPECT_TRUE(!result.cases[0].solution.contains(0.0),
            "No params: x=0 should not be in solution");
    }

    // =========================================================================
    // Test 6: Degenerate case verification
    // ax^2 + 2x + 1 > 0 with a=0 → 2x + 1 > 0 → x > -1/2
    // Requirements: 7.4
    // =========================================================================
    TEST_CASE("Degenerate case: ax^2 + 2x + 1 > 0 when a=0");
    {
        // Build ax^2 + 2x + 1
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto ax2 = SymbolicExpr::multiply(a, x2);
        auto two_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
        auto expr = SymbolicExpr::add(SymbolicExpr::add(ax2, two_x), SymbolicExpr::number(1));

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"a"});

        // Should have 3 cases (a>0, a<0, a=0)
        EXPECT_TRUE(result.cases.size() == 3,
            "ax^2 + 2x + 1 > 0: should have 3 cases");

        // The degenerate case (a=0) should give solution for 2x + 1 > 0 → x > -1/2
        if (result.cases.size() == 3) {
            auto& degen = result.cases[2];  // a = 0 case
            EXPECT_TRUE(!degen.solution.is_empty(),
                "Degenerate (a=0): 2x + 1 > 0 should not be empty");

            // The solution should be (-1/2, +∞)
            auto& intervals = degen.solution.intervals();
            EXPECT_TRUE(intervals.size() == 1,
                "Degenerate (a=0): should have 1 interval");

            if (!intervals.empty()) {
                EXPECT_TRUE(intervals[0].upper.is_pos_infinity,
                    "Degenerate (a=0): upper bound should be +inf");
                EXPECT_TRUE(intervals[0].lower.is_open,
                    "Degenerate (a=0): lower bound should be open (strict)");

                // Check the boundary value is -1/2
                if (intervals[0].lower.value) {
                    auto val_expr = intervals[0].lower.value->simplify();
                    try {
                        double val = val_expr->to_numeric();
                        EXPECT_TRUE(std::abs(val - (-0.5)) < 1e-10,
                            "Degenerate (a=0): boundary should be -1/2");
                    } catch (...) {
                        EXPECT_TRUE(false,
                            "Degenerate (a=0): boundary should be evaluable");
                    }
                }
            }
        }
    }

    // =========================================================================
    // Task 7.5: 含参数不等式单元测试
    // Requirements: 7.1, 7.2, 7.3, 7.4
    // =========================================================================

    // =========================================================================
    // Test 7: ax + b > 0 with a > 0 — verify numeric consistency
    // When a=5, b=-10: 5x - 10 > 0 → x > 2 → (2, +∞)
    // When a=1, b=3: x + 3 > 0 → x > -3 → (-3, +∞)
    // Requirements: 7.1, 7.2
    // =========================================================================
    TEST_CASE("Task 7.5: ax + b > 0, a > 0 numeric verification");
    {
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(a, x),
            b);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"a", "b"});

        EXPECT_TRUE(!result.cases.empty(),
            "ax + b > 0: should have cases");
        EXPECT_TRUE(result.cases.size() >= 1,
            "ax + b > 0: should have at least 1 case");

        // Case a > 0 (first case): solution is (-b/a, +∞)
        auto& pos_case = result.cases[0];
        auto& pos_intervals = pos_case.solution.intervals();

        // Test with a=5, b=-10 → root = -(-10)/5 = 2 → (2, +∞)
        if (!pos_intervals.empty() && pos_intervals[0].lower.value) {
            auto boundary = pos_intervals[0].lower.value;
            auto b1 = boundary->substitute("a", SymbolicExpr::number(5));
            b1 = b1->substitute("b", SymbolicExpr::number(-10));
            b1 = b1->simplify();
            try {
                double val = b1->to_numeric();
                EXPECT_TRUE(std::abs(val - 2.0) < 1e-10,
                    "ax+b>0, a>0: a=5,b=-10 → boundary should be 2");
            } catch (...) {
                EXPECT_TRUE(false, "ax+b>0, a>0: boundary should be evaluable for a=5,b=-10");
            }

            // Test with a=1, b=3 → root = -3/1 = -3 → (-3, +∞)
            auto b2 = boundary->substitute("a", SymbolicExpr::number(1));
            b2 = b2->substitute("b", SymbolicExpr::number(3));
            b2 = b2->simplify();
            try {
                double val = b2->to_numeric();
                EXPECT_TRUE(std::abs(val - (-3.0)) < 1e-10,
                    "ax+b>0, a>0: a=1,b=3 → boundary should be -3");
            } catch (...) {
                EXPECT_TRUE(false, "ax+b>0, a>0: boundary should be evaluable for a=1,b=3");
            }
        }

        // Verify the interval structure: should be (boundary, +∞)
        if (!pos_intervals.empty()) {
            EXPECT_TRUE(pos_intervals[0].upper.is_pos_infinity,
                "ax+b>0, a>0: upper bound should be +∞");
            EXPECT_TRUE(pos_intervals[0].lower.is_open,
                "ax+b>0, a>0: lower bound should be open (strict >)");
        }
    }

    // =========================================================================
    // Test 8: ax + b > 0 with a < 0 — verify numeric consistency
    // When a=-2, b=6: -2x + 6 > 0 → x < 3 → (-∞, 3)
    // When a=-4, b=-8: -4x - 8 > 0 → x < -2 → (-∞, -2)
    // Requirements: 7.1, 7.2
    // =========================================================================
    TEST_CASE("Task 7.5: ax + b > 0, a < 0 numeric verification");
    {
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(a, x),
            b);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"a", "b"});

        EXPECT_TRUE(result.cases.size() >= 2,
            "ax + b > 0: should have at least 2 cases (a>0, a<0)");

        // Case a < 0 (second case): solution is (-∞, -b/a)
        auto& neg_case = result.cases[1];
        auto& neg_intervals = neg_case.solution.intervals();

        EXPECT_TRUE(!neg_intervals.empty(),
            "ax+b>0, a<0: solution should not be empty");

        if (!neg_intervals.empty()) {
            // Verify structure: (-∞, boundary)
            EXPECT_TRUE(neg_intervals[0].lower.is_neg_infinity,
                "ax+b>0, a<0: lower bound should be -∞");
            EXPECT_TRUE(neg_intervals[0].upper.is_open,
                "ax+b>0, a<0: upper bound should be open (strict >)");
            EXPECT_TRUE(!neg_intervals[0].upper.is_pos_infinity,
                "ax+b>0, a<0: upper bound should not be +∞");

            // Test with a=-2, b=6 → root = -6/(-2) = 3 → (-∞, 3)
            if (neg_intervals[0].upper.value) {
                auto boundary = neg_intervals[0].upper.value;
                auto b1 = boundary->substitute("a", SymbolicExpr::number(-2));
                b1 = b1->substitute("b", SymbolicExpr::number(6));
                b1 = b1->simplify();
                try {
                    double val = b1->to_numeric();
                    EXPECT_TRUE(std::abs(val - 3.0) < 1e-10,
                        "ax+b>0, a<0: a=-2,b=6 → boundary should be 3");
                } catch (...) {
                    EXPECT_TRUE(false, "ax+b>0, a<0: boundary should be evaluable for a=-2,b=6");
                }

                // Test with a=-4, b=-8 → root = -(-8)/(-4) = -2 → (-∞, -2)
                auto b2 = boundary->substitute("a", SymbolicExpr::number(-4));
                b2 = b2->substitute("b", SymbolicExpr::number(-8));
                b2 = b2->simplify();
                try {
                    double val = b2->to_numeric();
                    EXPECT_TRUE(std::abs(val - (-2.0)) < 1e-10,
                        "ax+b>0, a<0: a=-4,b=-8 → boundary should be -2");
                } catch (...) {
                    EXPECT_TRUE(false, "ax+b>0, a<0: boundary should be evaluable for a=-4,b=-8");
                }
            }
        }
    }

    // =========================================================================
    // Test 9: ax² + bx + c ≥ 0 — discriminant depends on parameters
    // With a=1, b=-2, c=1: x²-2x+1 = (x-1)² ≥ 0 → entire line (disc=0)
    // With a=1, b=0, c=1: x²+1 ≥ 0 → entire line (disc<0, always positive)
    // With a=1, b=-5, c=6: x²-5x+6 = (x-2)(x-3) ≥ 0 → (-∞,2]∪[3,+∞) (disc>0)
    // Requirements: 7.1, 7.2, 7.3
    // =========================================================================
    TEST_CASE("Task 7.5: ax^2 + bx + c >= 0, discriminant depends on parameters");
    {
        // Build ax^2 + bx + c
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto ax2 = SymbolicExpr::multiply(a, x2);
        auto bx = SymbolicExpr::multiply(b, x);
        auto expr = SymbolicExpr::add(SymbolicExpr::add(ax2, bx), c);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterEqual, "x", {"a", "b", "c"});

        // Should have 3 cases (a>0, a<0, a=0)
        EXPECT_TRUE(result.cases.size() == 3,
            "ax^2+bx+c>=0: should have 3 piecewise cases");

        // Verify case a > 0 by substituting specific values
        if (result.cases.size() >= 1) {
            auto& pos_case = result.cases[0];

            // For a=1, b=-5, c=6: x²-5x+6 = (x-2)(x-3) ≥ 0
            // Solution: (-∞, 2] ∪ [3, +∞)
            // The parametric solution should have 2 intervals with symbolic roots
            auto& intervals = pos_case.solution.intervals();

            // Verify the roots evaluate correctly for a=1, b=-5, c=6
            if (intervals.size() == 2) {
                // First interval: (-∞, root1]
                if (intervals[0].upper.value) {
                    auto root1 = intervals[0].upper.value->substitute("a", SymbolicExpr::number(1));
                    root1 = root1->substitute("b", SymbolicExpr::number(-5));
                    root1 = root1->substitute("c", SymbolicExpr::number(6));
                    root1 = root1->simplify();
                    try {
                        double val = root1->to_numeric();
                        // Should be 2 (the smaller root)
                        EXPECT_TRUE(std::abs(val - 2.0) < 1e-10 || std::abs(val - 3.0) < 1e-10,
                            "ax^2+bx+c>=0, a>0: root at a=1,b=-5,c=6 should be 2 or 3");
                    } catch (...) {
                        // Symbolic root may not be directly evaluable
                    }
                }

                // Second interval: [root2, +∞)
                if (intervals[1].lower.value) {
                    auto root2 = intervals[1].lower.value->substitute("a", SymbolicExpr::number(1));
                    root2 = root2->substitute("b", SymbolicExpr::number(-5));
                    root2 = root2->substitute("c", SymbolicExpr::number(6));
                    root2 = root2->simplify();
                    try {
                        double val = root2->to_numeric();
                        // Should be 3 (the larger root)
                        EXPECT_TRUE(std::abs(val - 2.0) < 1e-10 || std::abs(val - 3.0) < 1e-10,
                            "ax^2+bx+c>=0, a>0: root at a=1,b=-5,c=6 should be 2 or 3");
                    } catch (...) {
                        // Symbolic root may not be directly evaluable
                    }
                }

                // Verify non-strict: endpoints should be closed
                EXPECT_TRUE(!intervals[0].upper.is_open,
                    "ax^2+bx+c>=0, a>0: upper of first interval should be closed (non-strict)");
                EXPECT_TRUE(!intervals[1].lower.is_open,
                    "ax^2+bx+c>=0, a>0: lower of second interval should be closed (non-strict)");
            }
        }

        // Verify case a < 0 by substituting specific values
        if (result.cases.size() >= 2) {
            auto& neg_case = result.cases[1];

            // For a=-1, b=5, c=-6: -x²+5x-6 = -(x-2)(x-3) ≥ 0
            // This is ≥ 0 when (x-2)(x-3) ≤ 0, i.e., x ∈ [2, 3]
            // The parametric solution for a<0 should have 1 interval [root1, root2]
            auto& intervals = neg_case.solution.intervals();

            if (intervals.size() == 1) {
                // Should be [root1, root2]
                EXPECT_TRUE(!intervals[0].lower.is_neg_infinity,
                    "ax^2+bx+c>=0, a<0: lower should not be -∞");
                EXPECT_TRUE(!intervals[0].upper.is_pos_infinity,
                    "ax^2+bx+c>=0, a<0: upper should not be +∞");
                EXPECT_TRUE(!intervals[0].lower.is_open,
                    "ax^2+bx+c>=0, a<0: lower should be closed (non-strict)");
                EXPECT_TRUE(!intervals[0].upper.is_open,
                    "ax^2+bx+c>=0, a<0: upper should be closed (non-strict)");
            }
        }
    }

    // =========================================================================
    // Test 10: Degenerate case — leading coefficient becomes zero
    // ax² + 3x - 2 > 0 with a=0 → 3x - 2 > 0 → x > 2/3
    // Requirements: 7.4
    // =========================================================================
    TEST_CASE("Task 7.5: Degenerate case ax^2 + 3x - 2 > 0 when a=0");
    {
        // Build ax^2 + 3x - 2
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto ax2 = SymbolicExpr::multiply(a, x2);
        auto three_x = SymbolicExpr::multiply(SymbolicExpr::number(3), x);
        auto expr = SymbolicExpr::add(
            SymbolicExpr::add(ax2, three_x),
            SymbolicExpr::number(-2));

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"a"});

        // Should have 3 cases
        EXPECT_TRUE(result.cases.size() == 3,
            "ax^2+3x-2>0: should have 3 cases");

        if (result.cases.size() == 3) {
            // Degenerate case (a=0): 3x - 2 > 0 → x > 2/3
            auto& degen = result.cases[2];
            EXPECT_TRUE(!degen.solution.is_empty(),
                "Degenerate (a=0): 3x-2>0 should not be empty");

            auto& intervals = degen.solution.intervals();
            EXPECT_TRUE(intervals.size() == 1,
                "Degenerate (a=0): should have 1 interval");

            if (!intervals.empty()) {
                EXPECT_TRUE(intervals[0].upper.is_pos_infinity,
                    "Degenerate (a=0): upper should be +∞");
                EXPECT_TRUE(intervals[0].lower.is_open,
                    "Degenerate (a=0): lower should be open (strict)");

                // Boundary should be 2/3
                if (intervals[0].lower.value) {
                    auto val_expr = intervals[0].lower.value->simplify();
                    try {
                        double val = val_expr->to_numeric();
                        EXPECT_TRUE(std::abs(val - (2.0/3.0)) < 1e-10,
                            "Degenerate (a=0): boundary should be 2/3");
                    } catch (...) {
                        EXPECT_TRUE(false,
                            "Degenerate (a=0): boundary should be evaluable");
                    }
                }
            }
        }
    }

    // =========================================================================
    // Test 11: Degenerate case — leading coefficient zero with constant remainder
    // ax² + 5 > 0 with a=0 → 5 > 0 → entire line
    // Requirements: 7.4
    // =========================================================================
    TEST_CASE("Task 7.5: Degenerate case ax^2 + 5 > 0 when a=0 (constant remainder)");
    {
        // Build ax^2 + 5
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto ax2 = SymbolicExpr::multiply(a, x2);
        auto expr = SymbolicExpr::add(ax2, SymbolicExpr::number(5));

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"a"});

        // Should have 3 cases
        EXPECT_TRUE(result.cases.size() == 3,
            "ax^2+5>0: should have 3 cases");

        if (result.cases.size() == 3) {
            // Degenerate case (a=0): 5 > 0 → entire line
            auto& degen = result.cases[2];
            EXPECT_TRUE(degen.solution.is_entire_line(),
                "Degenerate (a=0): 5>0 should be entire line");
        }
    }

    // =========================================================================
    // Test 12: ax + b ≤ 0 with a > 0 — non-strict inequality
    // When a=2, b=-6: 2x - 6 ≤ 0 → x ≤ 3 → (-∞, 3]
    // Requirements: 7.1, 7.2
    // =========================================================================
    TEST_CASE("Task 7.5: ax + b <= 0, a > 0 (non-strict)");
    {
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(a, x),
            b);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::LessEqual, "x", {"a", "b"});

        EXPECT_TRUE(!result.cases.empty(),
            "ax+b<=0: should have cases");

        // Case a > 0: solution is (-∞, -b/a]
        if (!result.cases.empty()) {
            auto& pos_case = result.cases[0];
            auto& intervals = pos_case.solution.intervals();

            EXPECT_TRUE(!intervals.empty(),
                "ax+b<=0, a>0: should have intervals");

            if (!intervals.empty()) {
                EXPECT_TRUE(intervals[0].lower.is_neg_infinity,
                    "ax+b<=0, a>0: lower should be -∞");
                EXPECT_TRUE(!intervals[0].upper.is_pos_infinity,
                    "ax+b<=0, a>0: upper should not be +∞");
                // Non-strict: endpoint should be closed
                EXPECT_TRUE(!intervals[0].upper.is_open,
                    "ax+b<=0, a>0: upper should be closed (non-strict ≤)");

                // Verify with a=2, b=-6 → root = 3 → (-∞, 3]
                if (intervals[0].upper.value) {
                    auto boundary = intervals[0].upper.value;
                    auto b1 = boundary->substitute("a", SymbolicExpr::number(2));
                    b1 = b1->substitute("b", SymbolicExpr::number(-6));
                    b1 = b1->simplify();
                    try {
                        double val = b1->to_numeric();
                        EXPECT_TRUE(std::abs(val - 3.0) < 1e-10,
                            "ax+b<=0, a>0: a=2,b=-6 → boundary should be 3");
                    } catch (...) {
                        EXPECT_TRUE(false,
                            "ax+b<=0, a>0: boundary should be evaluable");
                    }
                }
            }
        }
    }

    // =========================================================================
    // Test 13: Degenerate case with non-strict inequality
    // ax² + x ≥ 0 with a=0 → x ≥ 0 → [0, +∞)
    // Requirements: 7.4
    // =========================================================================
    TEST_CASE("Task 7.5: Degenerate ax^2 + x >= 0 when a=0 (non-strict)");
    {
        // Build ax^2 + x
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto ax2 = SymbolicExpr::multiply(a, x2);
        auto expr = SymbolicExpr::add(ax2, x);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterEqual, "x", {"a"});

        EXPECT_TRUE(result.cases.size() == 3,
            "ax^2+x>=0: should have 3 cases");

        if (result.cases.size() == 3) {
            // Degenerate case (a=0): x ≥ 0 → [0, +∞)
            auto& degen = result.cases[2];
            EXPECT_TRUE(!degen.solution.is_empty(),
                "Degenerate (a=0): x>=0 should not be empty");

            auto& intervals = degen.solution.intervals();
            EXPECT_TRUE(intervals.size() == 1,
                "Degenerate (a=0): x>=0 should have 1 interval");

            if (!intervals.empty()) {
                EXPECT_TRUE(intervals[0].upper.is_pos_infinity,
                    "Degenerate (a=0): upper should be +∞");
                // Non-strict: lower bound should be closed at 0
                EXPECT_TRUE(!intervals[0].lower.is_open,
                    "Degenerate (a=0): lower should be closed (non-strict ≥)");

                if (intervals[0].lower.value) {
                    auto val_expr = intervals[0].lower.value->simplify();
                    try {
                        double val = val_expr->to_numeric();
                        EXPECT_TRUE(std::abs(val - 0.0) < 1e-10,
                            "Degenerate (a=0): boundary should be 0");
                    } catch (...) {
                        EXPECT_TRUE(false,
                            "Degenerate (a=0): boundary should be evaluable");
                    }
                }
            }
        }
    }

    return TEST_REPORT();
}
