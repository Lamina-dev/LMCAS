#include "solve_polynomial.hpp"
#include "polynomial.hpp"
#include "rational.hpp"
#include "test_common.hpp"
#include <iostream>
#include <algorithm>
#include <map>

using namespace LMCAS;

Polynomial<Rational> poly_from_roots(const std::vector<Rational>& roots, const std::string& var = "x") {
    Polynomial<Rational> result({Rational(1)}, var);
    for (const auto& r : roots) {

        Polynomial<Rational> factor({-r, Rational(1)}, var);
        result = result * factor;
    }
    return result;
}

int count_occurrences(const std::vector<Rational>& vec, const Rational& val) {
    int count = 0;
    for (const auto& v : vec) {
        if (v == val) count++;
    }
    return count;
}

void test_simple_linear() {
    std::cout << "Test: simple linear (x - 2) ... ";

    Polynomial<Rational> p({Rational(-2), Rational(1)}, "x");
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.size() == 1, "linear polynomial has one rational root");
    EXPECT_TRUE(!roots.empty() && roots[0] == Rational(2), "linear root is 2");
    std::cout << "PASSED\n";
}

void test_quadratic_two_roots() {
    std::cout << "Test: quadratic (x-1)(x-3) ... ";

    auto p = poly_from_roots({Rational(1), Rational(3)});
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.size() == 2, "quadratic has two rational roots");
    EXPECT_TRUE(count_occurrences(roots, Rational(1)) == 1, "quadratic root 1 found once");
    EXPECT_TRUE(count_occurrences(roots, Rational(3)) == 1, "quadratic root 3 found once");
    std::cout << "PASSED\n";
}

void test_repeated_root() {
    std::cout << "Test: repeated root (x-2)^3 ... ";

    auto p = poly_from_roots({Rational(2), Rational(2), Rational(2)});
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.size() == 3, "repeated root preserves multiplicity");
    for (const auto& r : roots) {
        EXPECT_TRUE(r == Rational(2), "repeated root value is 2");
    }
    std::cout << "PASSED\n";
}

void test_rational_roots() {
    std::cout << "Test: rational roots (x - 1/2)(x - 3/4) ... ";

    auto p = poly_from_roots({Rational(1, 2), Rational(3, 4)});
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.size() == 2, "polynomial has two rational fractional roots");
    EXPECT_TRUE(count_occurrences(roots, Rational(1, 2)) == 1, "root 1/2 found once");
    EXPECT_TRUE(count_occurrences(roots, Rational(3, 4)) == 1, "root 3/4 found once");
    std::cout << "PASSED\n";
}

void test_zero_constant_term() {
    std::cout << "Test: zero constant term x^2(x-3) ... ";

    auto p = poly_from_roots({Rational(0), Rational(0), Rational(3)});
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.size() == 3, "zero constant term roots preserve multiplicity");
    EXPECT_TRUE(count_occurrences(roots, Rational(0)) == 2, "zero root found twice");
    EXPECT_TRUE(count_occurrences(roots, Rational(3)) == 1, "root 3 found once");
    std::cout << "PASSED\n";
}

void test_negative_roots() {
    std::cout << "Test: negative roots (x+1)(x+2)(x-3) ... ";
    auto p = poly_from_roots({Rational(-1), Rational(-2), Rational(3)});
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.size() == 3, "negative-root polynomial has three roots");
    EXPECT_TRUE(count_occurrences(roots, Rational(-1)) == 1, "root -1 found once");
    EXPECT_TRUE(count_occurrences(roots, Rational(-2)) == 1, "root -2 found once");
    EXPECT_TRUE(count_occurrences(roots, Rational(3)) == 1, "root 3 found once");
    std::cout << "PASSED\n";
}

void test_no_rational_roots() {
    std::cout << "Test: no rational roots (x^2 + 1) ... ";

    Polynomial<Rational> p({Rational(1), Rational(0), Rational(1)}, "x");
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.empty(), "x^2 + 1 has no rational roots");
    std::cout << "PASSED\n";
}

void test_mixed_multiplicity() {
    std::cout << "Test: mixed multiplicity (x-1)^2(x+1)(x-2) ... ";
    auto p = poly_from_roots({Rational(1), Rational(1), Rational(-1), Rational(2)});
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.size() == 4, "mixed multiplicity roots preserve total multiplicity");
    EXPECT_TRUE(count_occurrences(roots, Rational(1)) == 2, "root 1 found twice");
    EXPECT_TRUE(count_occurrences(roots, Rational(-1)) == 1, "root -1 found once");
    EXPECT_TRUE(count_occurrences(roots, Rational(2)) == 1, "root 2 found once");
    std::cout << "PASSED\n";
}

void test_high_degree_partial() {
    std::cout << "Test: degree 6 with rational + irrational roots ... ";

    Polynomial<Rational> p1({Rational(2), Rational(-3), Rational(1)}, "x");
    Polynomial<Rational> p2({Rational(1), Rational(0), Rational(1)}, "x");
    Polynomial<Rational> p3({Rational(1), Rational(1), Rational(1)}, "x");
    auto p = p1 * p2 * p3;
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.size() == 2, "high-degree mixed polynomial finds rational roots only");
    EXPECT_TRUE(count_occurrences(roots, Rational(1)) == 1, "high-degree root 1 found once");
    EXPECT_TRUE(count_occurrences(roots, Rational(2)) == 1, "high-degree root 2 found once");
    std::cout << "PASSED\n";
}

void test_zero_polynomial() {
    std::cout << "Test: zero polynomial ... ";
    Polynomial<Rational> p("x");
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.empty(), "zero polynomial returns no finite rational roots");
    std::cout << "PASSED\n";
}

void test_constant_polynomial() {
    std::cout << "Test: constant polynomial (5) ... ";
    Polynomial<Rational> p({Rational(5)}, "x");
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.empty(), "constant polynomial returns no rational roots");
    std::cout << "PASSED\n";
}

void test_divisors_beyond_legacy_scan_limit() {
    std::cout << "Test: rational roots with prime divisors above 1000 ... ";
    auto p = poly_from_roots({Rational(1009), Rational(1013)});
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(
        roots.size() == 2,
        "candidate enumeration must not stop before sqrt(constant term)");
    EXPECT_TRUE(
        count_occurrences(roots, Rational(1009)) == 1,
        "prime root 1009 is found");
    EXPECT_TRUE(
        count_occurrences(roots, Rational(1013)) == 1,
        "prime root 1013 is found");
    std::cout << "PASSED\n";
}

void test_large_common_scale() {
    auto p = poly_from_roots({Rational(-1, 2), Rational(-1, 2), Rational(2, 3)});
    const Rational scale(BigInt(std::string("-1000000000000000000000000000000")),
                         BigInt(7));
    for (auto& coefficient : p.coeffs) {
        coefficient = coefficient * scale;
    }
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.size() == 3, "common rational scaling preserves root count");
    EXPECT_TRUE(count_occurrences(roots, Rational(-1, 2)) == 2,
                "common rational scaling preserves repeated fractional roots");
    EXPECT_TRUE(count_occurrences(roots, Rational(2, 3)) == 1,
                "common rational scaling preserves the remaining fractional root");
}

void test_large_linear_root_after_deflation() {
    const Rational root(BigInt(std::string("1000000000000000000000000000001")),
                        BigInt(std::string("1000000000000000000000000000003")));
    auto p = poly_from_roots({Rational(0), root});
    auto roots = find_rational_roots(p);
    EXPECT_TRUE(roots.size() == 2, "linear remainder returns both exact roots");
    EXPECT_TRUE(count_occurrences(roots, Rational(0)) == 1,
                "zero root is retained before linear solution");
    EXPECT_TRUE(count_occurrences(roots, root) == 1,
                "large rational linear root is solved without divisor enumeration");
}

void test_degree_stops_at_4() {
    std::cout << "Test: stops when degree <= 4 ... ";

    Polynomial<Rational> p1({Rational(-1), Rational(1)}, "x");
    Polynomial<Rational> p2({Rational(-2), Rational(1)}, "x");
    Polynomial<Rational> p3({Rational(-3), Rational(1)}, "x");
    Polynomial<Rational> p4({Rational(1), Rational(0), Rational(1)}, "x");
    Polynomial<Rational> p5({Rational(1), Rational(1), Rational(1)}, "x");
    auto p = p1 * p2 * p3 * p4 * p5;
    auto roots = find_rational_roots(p);

    EXPECT_TRUE(roots.size() == 3, "degree-stop case finds three rational roots");
    EXPECT_TRUE(count_occurrences(roots, Rational(1)) == 1, "degree-stop root 1 found once");
    EXPECT_TRUE(count_occurrences(roots, Rational(2)) == 1, "degree-stop root 2 found once");
    EXPECT_TRUE(count_occurrences(roots, Rational(3)) == 1, "degree-stop root 3 found once");
    std::cout << "PASSED\n";
}

#include <random>
#include <set>
#include <sstream>

int count_rational(const std::vector<Rational>& vec, const Rational& val) {
    int c = 0;
    for (const auto& v : vec) {
        if (v == val) c++;
    }
    return c;
}

struct PropertyTestResult {
    int trials_passed;
    int trials_total;
    std::string first_failure_msg;
};

PropertyTestResult test_rational_root_completeness() {
    std::cout << "\n=== Rational root completeness ===\n";

    const int NUM_TRIALS = 40;
    int pass_count = 0;
    std::string first_failure;

    std::mt19937 rng(98765);

    std::uniform_int_distribution<int> num_dist(-5, 5);
    std::uniform_int_distribution<int> den_dist(1, 4);
    std::uniform_int_distribution<int> num_roots_dist(1, 5);
    std::uniform_int_distribution<int> irreducible_choice(0, 4);

    auto make_irreducible = [](int choice) -> Polynomial<Rational> {
        switch (choice) {
            case 0:
                return Polynomial<Rational>(
                    {Rational(1), Rational(1), Rational(0), Rational(0), Rational(0), Rational(1)}, "x");
            case 1:
                return Polynomial<Rational>(
                    {Rational(1), Rational(0), Rational(1), Rational(0), Rational(0), Rational(1)}, "x");
            case 2:
                return Polynomial<Rational>(
                    {Rational(-1), Rational(-1), Rational(0), Rational(0), Rational(0), Rational(1)}, "x");
            case 3:
                return Polynomial<Rational>(
                    {Rational(2), Rational(0), Rational(0), Rational(0), Rational(0), Rational(1)}, "x");
            case 4:
                return Polynomial<Rational>(
                    {Rational(1), Rational(0), Rational(0), Rational(0), Rational(1), Rational(1)}, "x");
            default:
                return Polynomial<Rational>(
                    {Rational(1), Rational(1), Rational(0), Rational(0), Rational(0), Rational(1)}, "x");
        }
    };

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {

        int num_roots = num_roots_dist(rng);
        std::vector<Rational> expected_roots;

        for (int i = 0; i < num_roots; ++i) {
            int p = num_dist(rng);
            int q = den_dist(rng);

            expected_roots.push_back(Rational(p, q));
        }

        Polynomial<Rational> poly({Rational(1)}, "x");
        for (const auto& r : expected_roots) {
            Polynomial<Rational> factor({-r, Rational(1)}, "x");
            poly = poly * factor;
        }

        int choice = irreducible_choice(rng);
        Polynomial<Rational> irr = make_irreducible(choice);
        poly = poly * irr;

        auto found_roots = find_rational_roots(poly);

        bool trial_ok = true;
        std::ostringstream failure_detail;

        std::map<std::string, int> expected_mult;
        for (const auto& r : expected_roots) {
            expected_mult[r.to_string()]++;
        }

        std::map<std::string, int> found_mult;
        for (const auto& r : found_roots) {
            found_mult[r.to_string()]++;
        }

        for (const auto& [root_str, exp_count] : expected_mult) {
            int found_count = 0;
            if (found_mult.find(root_str) != found_mult.end()) {
                found_count = found_mult[root_str];
            }
            if (found_count < exp_count) {
                trial_ok = false;
                failure_detail << "Trial " << trial << ": root " << root_str
                               << " expected multiplicity " << exp_count
                               << " but found " << found_count
                               << " (poly degree=" << poly.degree()
                               << ", roots=[";
                for (size_t i = 0; i < expected_roots.size(); ++i) {
                    if (i > 0) failure_detail << ", ";
                    failure_detail << expected_roots[i].to_string();
                }
                failure_detail << "])";
                break;
            }
        }

        if (trial_ok) {
            pass_count++;
        } else {
            if (first_failure.empty()) {
                first_failure = failure_detail.str();
            }
            std::cout << "[FAIL] " << failure_detail.str() << "\n";
        }
    }

    std::cout << "Rational root completeness: " << pass_count
              << "/" << NUM_TRIALS << " trials passed\n";

    return {pass_count, NUM_TRIALS, first_failure};
}

int main() {
    std::cout << "=== find_rational_roots tests ===\n";
    test_simple_linear();
    test_quadratic_two_roots();
    test_repeated_root();
    test_rational_roots();
    test_zero_constant_term();
    test_negative_roots();
    test_no_rational_roots();
    test_mixed_multiplicity();
    test_high_degree_partial();
    test_zero_polynomial();
    test_constant_polynomial();
    test_divisors_beyond_legacy_scan_limit();
    test_degree_stops_at_4();
    test_large_common_scale();
    test_large_linear_root_after_deflation();

    auto result = test_rational_root_completeness();

    if (result.trials_passed < result.trials_total) {
        std::cerr << "\nRational root completeness failed: "
                  << result.trials_passed << "/" << result.trials_total
                  << " trials passed\n";
        if (!result.first_failure_msg.empty()) {
            std::cerr << "First failure: " << result.first_failure_msg << "\n";
        }
        return 1;
    }

    std::cout << "\nAll tests PASSED!\n";
    return TEST_REPORT();
}
