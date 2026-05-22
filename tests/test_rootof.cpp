#include "test_common.hpp"
#include "root_of_utils.hpp"
#include "solve_polynomial.hpp"
#include "poly_utils.hpp"
#include <cmath>
#include <algorithm>
#include <vector>

// Helper: create a SymbolicExpr number from int
static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

// Helper: evaluate a symbolic expression numerically (handles common node types)
static double eval_numeric(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !expr->root) return std::nan("");
    return expr->to_numeric();
}

int main() {
    // =========================================================================
    // Test 1: Out-of-range index returns nullopt (k >= degree)
    // Validates: Requirements 7.7
    // =========================================================================
    TEST_CASE("RootOf - Out-of-range index k >= degree returns nullopt");
    {
        // Build a degree-3 polynomial: x^3 + x + 1
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(SymbolicExpr::variable("x"), num(3)),
                SymbolicExpr::variable("x")),
            num(1));

        // Create RootOf with k = 3 (out of range for degree 3, valid range is 0..2)
        auto rootof_k3 = SymbolicExpr::root_of(poly_expr, "x", 3);
        auto result = lamina::rootof_evaluate(rootof_k3);
        EXPECT_TRUE(!result.has_value(),
            "rootof_evaluate with k=3 on degree-3 poly returns nullopt");

        // Create RootOf with k = 5 (well beyond range)
        auto rootof_k5 = SymbolicExpr::root_of(poly_expr, "x", 5);
        result = lamina::rootof_evaluate(rootof_k5);
        EXPECT_TRUE(!result.has_value(),
            "rootof_evaluate with k=5 on degree-3 poly returns nullopt");

        // Create RootOf with k = 100 (extreme out of range)
        auto rootof_k100 = SymbolicExpr::root_of(poly_expr, "x", 100);
        result = lamina::rootof_evaluate(rootof_k100);
        EXPECT_TRUE(!result.has_value(),
            "rootof_evaluate with k=100 on degree-3 poly returns nullopt");
    }

    // =========================================================================
    // Test 2: Out-of-range index returns nullopt (k < 0)
    // Validates: Requirements 7.7
    // =========================================================================
    TEST_CASE("RootOf - Out-of-range index k < 0 returns nullopt");
    {
        // Build a degree-3 polynomial: x^3 - 2x + 1
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(SymbolicExpr::variable("x"), num(3)),
                SymbolicExpr::multiply(num(-2), SymbolicExpr::variable("x"))),
            num(1));

        // Create RootOf with k = -1 (negative index)
        auto rootof_neg1 = SymbolicExpr::root_of(poly_expr, "x", -1);
        auto result = lamina::rootof_evaluate(rootof_neg1);
        EXPECT_TRUE(!result.has_value(),
            "rootof_evaluate with k=-1 returns nullopt");

        // Create RootOf with k = -10 (very negative)
        auto rootof_neg10 = SymbolicExpr::root_of(poly_expr, "x", -10);
        result = lamina::rootof_evaluate(rootof_neg10);
        EXPECT_TRUE(!result.has_value(),
            "rootof_evaluate with k=-10 returns nullopt");
    }

    // =========================================================================
    // Test 3: Parametric coefficients - rootof_evaluate returns nullopt
    // Validates: Requirements 7.3
    // =========================================================================
    TEST_CASE("RootOf - Parametric coefficients returns nullopt");
    {
        // Build a polynomial with symbolic parameter: x^3 + a*x + 1
        // The coefficient 'a' is a free variable, so numerical ordering is impossible
        auto x = SymbolicExpr::variable("x");
        auto a = SymbolicExpr::variable("a");
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(x, num(3)),
                SymbolicExpr::multiply(a, x)),
            num(1));

        // Even with a valid index k=0, parametric coefficients prevent evaluation
        auto rootof_k0 = SymbolicExpr::root_of(poly_expr, "x", 0);
        auto result = lamina::rootof_evaluate(rootof_k0);
        EXPECT_TRUE(!result.has_value(),
            "rootof_evaluate with parametric coeff 'a' returns nullopt (k=0)");

        auto rootof_k1 = SymbolicExpr::root_of(poly_expr, "x", 1);
        result = lamina::rootof_evaluate(rootof_k1);
        EXPECT_TRUE(!result.has_value(),
            "rootof_evaluate with parametric coeff 'a' returns nullopt (k=1)");

        auto rootof_k2 = SymbolicExpr::root_of(poly_expr, "x", 2);
        result = lamina::rootof_evaluate(rootof_k2);
        EXPECT_TRUE(!result.has_value(),
            "rootof_evaluate with parametric coeff 'a' returns nullopt (k=2)");
    }

    // =========================================================================
    // Test 4: Parametric coefficients with multiple free variables
    // Validates: Requirements 7.3
    // =========================================================================
    TEST_CASE("RootOf - Multiple parametric coefficients returns nullopt");
    {
        // x^2 + b*x + c (two free parameters)
        auto x = SymbolicExpr::variable("x");
        auto b = SymbolicExpr::variable("b");
        auto c = SymbolicExpr::variable("c");
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(x, num(2)),
                SymbolicExpr::multiply(b, x)),
            c);

        auto rootof_k0 = SymbolicExpr::root_of(poly_expr, "x", 0);
        auto result = lamina::rootof_evaluate(rootof_k0);
        EXPECT_TRUE(!result.has_value(),
            "rootof_evaluate with multiple parametric coeffs returns nullopt");
    }

    // =========================================================================
    // Test 5: Reducible polynomial degree <= 4 - rootof_simplify returns closed-form
    // Validates: Requirements 7.5
    // =========================================================================
    TEST_CASE("RootOf - Simplify degree-2 polynomial to closed-form");
    {
        // x^2 - 4 = 0 has roots -2 and 2
        auto x = SymbolicExpr::variable("x");
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::power(x, num(2)),
            num(-4));

        // RootOf(x^2 - 4, x, 0) should simplify to the smaller root (-2)
        auto rootof_k0 = SymbolicExpr::root_of(poly_expr, "x", 0);
        auto simplified_k0 = lamina::rootof_simplify(rootof_k0);

        // The simplified result should NOT be a RootOf node
        std::string s0 = simplified_k0->to_string();
        EXPECT_TRUE(s0.find("RootOf") == std::string::npos,
            "rootof_simplify(degree-2, k=0) returns non-RootOf expression");

        // Numerically, the k=0 root (smallest real) should be -2
        double val0 = eval_numeric(simplified_k0);
        EXPECT_TRUE(!std::isnan(val0) && std::abs(val0 - (-2.0)) < 1e-10,
            "rootof_simplify(x^2-4, x, 0) = -2");

        // RootOf(x^2 - 4, x, 1) should simplify to the larger root (2)
        auto rootof_k1 = SymbolicExpr::root_of(poly_expr, "x", 1);
        auto simplified_k1 = lamina::rootof_simplify(rootof_k1);

        std::string s1 = simplified_k1->to_string();
        EXPECT_TRUE(s1.find("RootOf") == std::string::npos,
            "rootof_simplify(degree-2, k=1) returns non-RootOf expression");

        double val1 = eval_numeric(simplified_k1);
        EXPECT_TRUE(!std::isnan(val1) && std::abs(val1 - 2.0) < 1e-10,
            "rootof_simplify(x^2-4, x, 1) = 2");
    }

    // =========================================================================
    // Test 6: Reducible polynomial degree 3 - rootof_simplify returns closed-form
    // Validates: Requirements 7.5
    // =========================================================================
    TEST_CASE("RootOf - Simplify degree-3 polynomial to closed-form");
    {
        // x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3) = 0, roots: 1, 2, 3
        auto x = SymbolicExpr::variable("x");
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::add(
                    SymbolicExpr::power(x, num(3)),
                    SymbolicExpr::multiply(num(-6), SymbolicExpr::power(x, num(2)))),
                SymbolicExpr::multiply(num(11), x)),
            num(-6));

        // k=0 should be the smallest root (1)
        auto rootof_k0 = SymbolicExpr::root_of(poly_expr, "x", 0);
        auto simplified_k0 = lamina::rootof_simplify(rootof_k0);
        std::string s0 = simplified_k0->to_string();
        EXPECT_TRUE(s0.find("RootOf") == std::string::npos,
            "rootof_simplify(degree-3, k=0) returns non-RootOf expression");
        double val0 = eval_numeric(simplified_k0);
        EXPECT_TRUE(!std::isnan(val0) && std::abs(val0 - 1.0) < 1e-8,
            "rootof_simplify(cubic, k=0) = 1 (smallest root)");

        // k=1 should be the middle root (2)
        auto rootof_k1 = SymbolicExpr::root_of(poly_expr, "x", 1);
        auto simplified_k1 = lamina::rootof_simplify(rootof_k1);
        std::string s1 = simplified_k1->to_string();
        EXPECT_TRUE(s1.find("RootOf") == std::string::npos,
            "rootof_simplify(degree-3, k=1) returns non-RootOf expression");
        double val1 = eval_numeric(simplified_k1);
        EXPECT_TRUE(!std::isnan(val1) && std::abs(val1 - 2.0) < 1e-8,
            "rootof_simplify(cubic, k=1) = 2 (middle root)");

        // k=2 should be the largest root (3)
        auto rootof_k2 = SymbolicExpr::root_of(poly_expr, "x", 2);
        auto simplified_k2 = lamina::rootof_simplify(rootof_k2);
        std::string s2 = simplified_k2->to_string();
        EXPECT_TRUE(s2.find("RootOf") == std::string::npos,
            "rootof_simplify(degree-3, k=2) returns non-RootOf expression");
        double val2 = eval_numeric(simplified_k2);
        EXPECT_TRUE(!std::isnan(val2) && std::abs(val2 - 3.0) < 1e-8,
            "rootof_simplify(cubic, k=2) = 3 (largest root)");
    }

    // =========================================================================
    // Test 7: Reducible polynomial degree 4 - rootof_simplify returns closed-form
    // Validates: Requirements 7.5
    // =========================================================================
    TEST_CASE("RootOf - Simplify degree-4 polynomial to closed-form");
    {
        // x^4 - 5x^2 + 4 = (x^2-1)(x^2-4) = (x-1)(x+1)(x-2)(x+2)
        // Roots: -2, -1, 1, 2
        auto x = SymbolicExpr::variable("x");
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(x, num(4)),
                SymbolicExpr::multiply(num(-5), SymbolicExpr::power(x, num(2)))),
            num(4));

        // Verify all 4 simplified roots are closed-form (non-RootOf)
        bool all_non_rootof = true;
        for (int k = 0; k < 4; ++k) {
            auto rootof_k = SymbolicExpr::root_of(poly_expr, "x", k);
            auto simplified = lamina::rootof_simplify(rootof_k);
            std::string s = simplified->to_string();
            if (s.find("RootOf") != std::string::npos) {
                all_non_rootof = false;
            }
        }

        EXPECT_TRUE(all_non_rootof,
            "rootof_simplify(degree-4) returns non-RootOf for all indices");

        // Verify via rootof_evaluate that the polynomial has 4 real roots
        // at the expected values {-2, -1, 1, 2}
        std::vector<double> eval_values;
        for (int k = 0; k < 4; ++k) {
            auto rootof_k = SymbolicExpr::root_of(poly_expr, "x", k);
            auto result = lamina::rootof_evaluate(rootof_k);
            if (result.has_value()) {
                eval_values.push_back(result.value());
            }
        }

        // rootof_evaluate should find all 4 real roots
        EXPECT_TRUE(eval_values.size() == 4,
            "rootof_evaluate(quartic) finds 4 real roots");

        if (eval_values.size() == 4) {
            std::sort(eval_values.begin(), eval_values.end());
            double expected[] = {-2.0, -1.0, 1.0, 2.0};
            bool all_match = true;
            for (int i = 0; i < 4; ++i) {
                if (std::abs(eval_values[i] - expected[i]) > 1e-8) {
                    all_match = false;
                }
            }
            EXPECT_TRUE(all_match,
                "rootof_evaluate(quartic) roots are {-2, -1, 1, 2}");
        }
    }

    // =========================================================================
    // Test 8: Valid index on numeric polynomial - rootof_evaluate succeeds
    // Validates: Requirements 7.4 (sanity check that valid cases work)
    // =========================================================================
    TEST_CASE("RootOf - Valid index on numeric polynomial evaluates correctly");
    {
        // x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3), roots: 1, 2, 3
        auto x = SymbolicExpr::variable("x");
        auto poly_expr = SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::add(
                    SymbolicExpr::power(x, num(3)),
                    SymbolicExpr::multiply(num(-6), SymbolicExpr::power(x, num(2)))),
                SymbolicExpr::multiply(num(11), x)),
            num(-6));

        // k=0 should evaluate to the smallest real root (1.0)
        auto rootof_k0 = SymbolicExpr::root_of(poly_expr, "x", 0);
        auto result0 = lamina::rootof_evaluate(rootof_k0);
        EXPECT_TRUE(result0.has_value(),
            "rootof_evaluate with valid k=0 returns a value");
        if (result0.has_value()) {
            EXPECT_TRUE(std::abs(result0.value() - 1.0) < 1e-10,
                "rootof_evaluate(cubic, k=0) = 1.0 (smallest root)");
        }

        // k=2 should evaluate to the largest real root (3.0)
        auto rootof_k2 = SymbolicExpr::root_of(poly_expr, "x", 2);
        auto result2 = lamina::rootof_evaluate(rootof_k2);
        EXPECT_TRUE(result2.has_value(),
            "rootof_evaluate with valid k=2 returns a value");
        if (result2.has_value()) {
            EXPECT_TRUE(std::abs(result2.value() - 3.0) < 1e-10,
                "rootof_evaluate(cubic, k=2) = 3.0 (largest root)");
        }
    }

    return TEST_REPORT();
}
