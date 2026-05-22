// test_sturm_isolation.cpp - Test Sturm sequence real root isolation
// This test directly exercises isolate_real_roots via the header-only Polynomial class

#include "test_common.hpp"
#include "newton_raphson.hpp"

using lamina::Polynomial;

int main() {
    TEST_CASE("Sturm Isolation - Linear polynomial x - 3");
    {
        // x - 3 has one root at x=3
        Polynomial<Rational> p({Rational(-3), Rational(1)}, "x");
        auto intervals = lamina::isolate_real_roots(p);
        EXPECT_TRUE(intervals.size() == 1, "x-3 should have 1 real root");
        if (!intervals.empty()) {
            EXPECT_TRUE(intervals[0].first <= Rational(3) && intervals[0].second >= Rational(3),
                "Root x=3 should be in the interval");
        }
    }

    TEST_CASE("Sturm Isolation - Quadratic x^2 - 4 (roots at -2, 2)");
    {
        // x^2 - 4 = (x-2)(x+2)
        Polynomial<Rational> p({Rational(-4), Rational(0), Rational(1)}, "x");
        auto intervals = lamina::isolate_real_roots(p);
        EXPECT_TRUE(intervals.size() == 2, "x^2-4 should have 2 real roots");
        if (intervals.size() == 2) {
            // First interval should contain -2
            EXPECT_TRUE(intervals[0].first <= Rational(-2) && intervals[0].second >= Rational(-2),
                "First root x=-2 should be in first interval");
            // Second interval should contain 2
            EXPECT_TRUE(intervals[1].first <= Rational(2) && intervals[1].second >= Rational(2),
                "Second root x=2 should be in second interval");
        }
    }

    TEST_CASE("Sturm Isolation - Cubic x^3 - 6x^2 + 11x - 6 (roots at 1, 2, 3)");
    {
        // (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
        Polynomial<Rational> p({Rational(-6), Rational(11), Rational(-6), Rational(1)}, "x");
        auto intervals = lamina::isolate_real_roots(p);
        EXPECT_TRUE(intervals.size() == 3, "x^3-6x^2+11x-6 should have 3 real roots");
        if (intervals.size() == 3) {
            EXPECT_TRUE(intervals[0].first <= Rational(1) && intervals[0].second >= Rational(1),
                "Root x=1 should be in first interval");
            EXPECT_TRUE(intervals[1].first <= Rational(2) && intervals[1].second >= Rational(2),
                "Root x=2 should be in second interval");
            EXPECT_TRUE(intervals[2].first <= Rational(3) && intervals[2].second >= Rational(3),
                "Root x=3 should be in third interval");
        }
    }

    TEST_CASE("Sturm Isolation - No real roots: x^2 + 1");
    {
        // x^2 + 1 has no real roots
        Polynomial<Rational> p({Rational(1), Rational(0), Rational(1)}, "x");
        auto intervals = lamina::isolate_real_roots(p);
        EXPECT_TRUE(intervals.size() == 0, "x^2+1 should have 0 real roots");
    }

    TEST_CASE("Sturm Isolation - Double root: x^2 - 2x + 1 = (x-1)^2");
    {
        // (x-1)^2 = x^2 - 2x + 1, square-free part is (x-1)
        Polynomial<Rational> p({Rational(1), Rational(-2), Rational(1)}, "x");
        auto intervals = lamina::isolate_real_roots(p);
        // square_free_part reduces to (x-1), so 1 distinct root
        EXPECT_TRUE(intervals.size() == 1, "(x-1)^2 should have 1 distinct real root");
        if (!intervals.empty()) {
            EXPECT_TRUE(intervals[0].first <= Rational(1) && intervals[0].second >= Rational(1),
                "Root x=1 should be in the interval");
        }
    }

    TEST_CASE("Sturm Isolation - Quartic (x^2-2)(x^2-3) with 4 real roots");
    {
        // (x^2-2)(x^2-3) = x^4 - 5x^2 + 6
        // Roots: -sqrt(3), -sqrt(2), sqrt(2), sqrt(3)
        Polynomial<Rational> p({Rational(6), Rational(0), Rational(-5), Rational(0), Rational(1)}, "x");
        auto intervals = lamina::isolate_real_roots(p);
        EXPECT_TRUE(intervals.size() == 4, "(x^2-2)(x^2-3) should have 4 real roots");
    }

    TEST_CASE("Sturm Isolation - Constant polynomial");
    {
        Polynomial<Rational> p({Rational(5)}, "x");
        auto intervals = lamina::isolate_real_roots(p);
        EXPECT_TRUE(intervals.size() == 0, "Constant polynomial should have 0 roots");
    }

    TEST_CASE("Sturm Isolation - Zero polynomial");
    {
        Polynomial<Rational> p("x");
        auto intervals = lamina::isolate_real_roots(p);
        EXPECT_TRUE(intervals.size() == 0, "Zero polynomial should return empty");
    }

    return TEST_REPORT();
}
