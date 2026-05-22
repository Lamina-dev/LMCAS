// test_find_rational_roots.cpp - Unit tests for find_rational_roots
#include "solve_polynomial.hpp"
#include "polynomial.hpp"
#include "rational.hpp"
#include <iostream>
#include <cassert>
#include <algorithm>
#include <map>

using namespace lamina;

// Helper: build polynomial from roots (x - r1)(x - r2)...(x - rn)
Polynomial<Rational> poly_from_roots(const std::vector<Rational>& roots, const std::string& var = "x") {
    Polynomial<Rational> result({Rational(1)}, var);
    for (const auto& r : roots) {
        // Multiply by (x - r)
        Polynomial<Rational> factor({-r, Rational(1)}, var);
        result = result * factor;
    }
    return result;
}

// Helper: count occurrences of a value in a vector
int count_occurrences(const std::vector<Rational>& vec, const Rational& val) {
    int count = 0;
    for (const auto& v : vec) {
        if (v == val) count++;
    }
    return count;
}

void test_simple_linear() {
    std::cout << "Test: simple linear (x - 2) ... ";
    // x - 2 = 0 => root at 2
    Polynomial<Rational> p({Rational(-2), Rational(1)}, "x");
    auto roots = find_rational_roots(p);
    assert(roots.size() == 1);
    assert(roots[0] == Rational(2));
    std::cout << "PASSED\n";
}

void test_quadratic_two_roots() {
    std::cout << "Test: quadratic (x-1)(x-3) ... ";
    // (x-1)(x-3) = x^2 - 4x + 3
    auto p = poly_from_roots({Rational(1), Rational(3)});
    auto roots = find_rational_roots(p);
    assert(roots.size() == 2);
    assert(count_occurrences(roots, Rational(1)) == 1);
    assert(count_occurrences(roots, Rational(3)) == 1);
    std::cout << "PASSED\n";
}

void test_repeated_root() {
    std::cout << "Test: repeated root (x-2)^3 ... ";
    // (x-2)^3 = x^3 - 6x^2 + 12x - 8
    auto p = poly_from_roots({Rational(2), Rational(2), Rational(2)});
    auto roots = find_rational_roots(p);
    assert(roots.size() == 3);
    for (const auto& r : roots) {
        assert(r == Rational(2));
    }
    std::cout << "PASSED\n";
}

void test_rational_roots() {
    std::cout << "Test: rational roots (x - 1/2)(x - 3/4) ... ";
    // (x - 1/2)(x - 3/4) = x^2 - 5/4 x + 3/8
    auto p = poly_from_roots({Rational(1, 2), Rational(3, 4)});
    auto roots = find_rational_roots(p);
    assert(roots.size() == 2);
    assert(count_occurrences(roots, Rational(1, 2)) == 1);
    assert(count_occurrences(roots, Rational(3, 4)) == 1);
    std::cout << "PASSED\n";
}

void test_zero_constant_term() {
    std::cout << "Test: zero constant term x^2(x-3) ... ";
    // x^2(x-3) = x^3 - 3x^2
    auto p = poly_from_roots({Rational(0), Rational(0), Rational(3)});
    auto roots = find_rational_roots(p);
    assert(roots.size() == 3);
    assert(count_occurrences(roots, Rational(0)) == 2);
    assert(count_occurrences(roots, Rational(3)) == 1);
    std::cout << "PASSED\n";
}

void test_negative_roots() {
    std::cout << "Test: negative roots (x+1)(x+2)(x-3) ... ";
    auto p = poly_from_roots({Rational(-1), Rational(-2), Rational(3)});
    auto roots = find_rational_roots(p);
    assert(roots.size() == 3);
    assert(count_occurrences(roots, Rational(-1)) == 1);
    assert(count_occurrences(roots, Rational(-2)) == 1);
    assert(count_occurrences(roots, Rational(3)) == 1);
    std::cout << "PASSED\n";
}

void test_no_rational_roots() {
    std::cout << "Test: no rational roots (x^2 + 1) ... ";
    // x^2 + 1 has no rational roots
    Polynomial<Rational> p({Rational(1), Rational(0), Rational(1)}, "x");
    auto roots = find_rational_roots(p);
    assert(roots.empty());
    std::cout << "PASSED\n";
}

void test_mixed_multiplicity() {
    std::cout << "Test: mixed multiplicity (x-1)^2(x+1)(x-2) ... ";
    auto p = poly_from_roots({Rational(1), Rational(1), Rational(-1), Rational(2)});
    auto roots = find_rational_roots(p);
    assert(roots.size() == 4);
    assert(count_occurrences(roots, Rational(1)) == 2);
    assert(count_occurrences(roots, Rational(-1)) == 1);
    assert(count_occurrences(roots, Rational(2)) == 1);
    std::cout << "PASSED\n";
}

void test_high_degree_partial() {
    std::cout << "Test: degree 6 with rational + irrational roots ... ";
    // (x-1)(x-2)(x^2+1)(x^2+x+1) - only x=1 and x=2 are rational
    // Build: (x-1)(x-2) = x^2 - 3x + 2
    // (x^2+1) has no rational roots
    // (x^2+x+1) has no rational roots
    Polynomial<Rational> p1({Rational(2), Rational(-3), Rational(1)}, "x"); // (x-1)(x-2)
    Polynomial<Rational> p2({Rational(1), Rational(0), Rational(1)}, "x"); // x^2+1
    Polynomial<Rational> p3({Rational(1), Rational(1), Rational(1)}, "x"); // x^2+x+1
    auto p = p1 * p2 * p3;
    auto roots = find_rational_roots(p);
    assert(roots.size() == 2);
    assert(count_occurrences(roots, Rational(1)) == 1);
    assert(count_occurrences(roots, Rational(2)) == 1);
    std::cout << "PASSED\n";
}

void test_zero_polynomial() {
    std::cout << "Test: zero polynomial ... ";
    Polynomial<Rational> p("x");
    auto roots = find_rational_roots(p);
    assert(roots.empty());
    std::cout << "PASSED\n";
}

void test_constant_polynomial() {
    std::cout << "Test: constant polynomial (5) ... ";
    Polynomial<Rational> p({Rational(5)}, "x");
    auto roots = find_rational_roots(p);
    assert(roots.empty());
    std::cout << "PASSED\n";
}

void test_degree_stops_at_4() {
    std::cout << "Test: stops when degree <= 4 ... ";
    // (x-1)(x-2)(x^2+1)(x^2+x+1)(x-3) = degree 7
    // Rational roots: 1, 2, 3
    // After finding them, remaining is degree 4 (x^2+1)(x^2+x+1), should stop
    Polynomial<Rational> p1({Rational(-1), Rational(1)}, "x"); // (x-1)
    Polynomial<Rational> p2({Rational(-2), Rational(1)}, "x"); // (x-2)
    Polynomial<Rational> p3({Rational(-3), Rational(1)}, "x"); // (x-3)
    Polynomial<Rational> p4({Rational(1), Rational(0), Rational(1)}, "x"); // x^2+1
    Polynomial<Rational> p5({Rational(1), Rational(1), Rational(1)}, "x"); // x^2+x+1
    auto p = p1 * p2 * p3 * p4 * p5;
    auto roots = find_rational_roots(p);
    // Should find all 3 rational roots
    assert(roots.size() == 3);
    assert(count_occurrences(roots, Rational(1)) == 1);
    assert(count_occurrences(roots, Rational(2)) == 1);
    assert(count_occurrences(roots, Rational(3)) == 1);
    std::cout << "PASSED\n";
}

// =========================================================================
// Property 7: Rational root completeness
// Validates: Requirements 4.1, 4.2
//
// For any polynomial constructed as a product of linear factors with rational
// roots (e.g., (x-p1/q1)(x-p2/q2)...), optionally multiplied by an irreducible
// factor, the Preprocessor SHALL find all rational roots present in the polynomial.
//
// Strategy: Generate 30+ random polynomials from known rational roots with small
// numerators/denominators, optionally multiply by irreducible factors like (x^2+1)
// or (x^2+x+1). Verify all known rational roots are found by find_rational_roots.
// =========================================================================

#include <random>
#include <set>
#include <sstream>

// Helper: check if a rational value appears in a vector (counting multiplicity)
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

PropertyTestResult test_property7_rational_root_completeness() {
    std::cout << "\n=== Property 7: Rational root completeness ===\n";
    std::cout << "**Validates: Requirements 4.1, 4.2**\n";

    const int NUM_TRIALS = 40;
    int pass_count = 0;
    std::string first_failure;

    std::mt19937 rng(98765);
    // Small numerators and denominators to keep polynomials manageable
    std::uniform_int_distribution<int> num_dist(-5, 5);
    std::uniform_int_distribution<int> den_dist(1, 4);
    std::uniform_int_distribution<int> num_roots_dist(1, 5);
    std::uniform_int_distribution<int> irreducible_choice(0, 4);

    // Irreducible factors over Q with degree >= 5.
    // The find_rational_roots function stops when the remaining polynomial degree
    // drops to <= 4 (by design, per Requirement 4.2). To test completeness of
    // rational root finding, we multiply by an irreducible factor of degree >= 5
    // so that after all rational roots are extracted, the remaining degree is still > 4
    // and the function doesn't stop prematurely.
    //
    // Irreducible polynomials over Q of degree 5:
    //   x^5 + x + 1 (irreducible by Eisenstein-like arguments / direct check)
    //   x^5 + x^2 + 1 (irreducible over Q)
    //   x^5 - x - 1 (irreducible over Q)
    //   x^5 + 2 (Eisenstein with p=2)
    //   x^5 + x^4 + 1 (irreducible over Q)
    auto make_irreducible = [](int choice) -> Polynomial<Rational> {
        switch (choice) {
            case 0: // x^5 + x + 1 (irreducible over Q)
                return Polynomial<Rational>(
                    {Rational(1), Rational(1), Rational(0), Rational(0), Rational(0), Rational(1)}, "x");
            case 1: // x^5 + x^2 + 1 (irreducible over Q)
                return Polynomial<Rational>(
                    {Rational(1), Rational(0), Rational(1), Rational(0), Rational(0), Rational(1)}, "x");
            case 2: // x^5 - x - 1 (irreducible over Q)
                return Polynomial<Rational>(
                    {Rational(-1), Rational(-1), Rational(0), Rational(0), Rational(0), Rational(1)}, "x");
            case 3: // x^5 + 2 (Eisenstein with p=2, irreducible)
                return Polynomial<Rational>(
                    {Rational(2), Rational(0), Rational(0), Rational(0), Rational(0), Rational(1)}, "x");
            case 4: // x^5 + x^4 + 1 (irreducible over Q)
                return Polynomial<Rational>(
                    {Rational(1), Rational(0), Rational(0), Rational(0), Rational(1), Rational(1)}, "x");
            default:
                return Polynomial<Rational>(
                    {Rational(1), Rational(1), Rational(0), Rational(0), Rational(0), Rational(1)}, "x");
        }
    };

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        // Generate random rational roots
        int num_roots = num_roots_dist(rng);
        std::vector<Rational> expected_roots;

        for (int i = 0; i < num_roots; ++i) {
            int p = num_dist(rng);
            int q = den_dist(rng);
            // Ensure q != 0 (already guaranteed by den_dist range [1,4])
            expected_roots.push_back(Rational(p, q));
        }

        // Build polynomial as product of (x - p_i/q_i)
        Polynomial<Rational> poly({Rational(1)}, "x");
        for (const auto& r : expected_roots) {
            Polynomial<Rational> factor({-r, Rational(1)}, "x");
            poly = poly * factor;
        }

        // Always multiply by an irreducible factor of degree 5 to ensure
        // find_rational_roots doesn't stop early (it stops when degree <= 4).
        // This guarantees the function searches exhaustively for all rational roots.
        int choice = irreducible_choice(rng);
        Polynomial<Rational> irr = make_irreducible(choice);
        poly = poly * irr;

        // Find rational roots
        auto found_roots = find_rational_roots(poly);

        // Verify: every expected root must appear in found_roots with at least
        // the same multiplicity as in expected_roots
        bool trial_ok = true;
        std::ostringstream failure_detail;

        // Count expected multiplicities
        std::map<std::string, int> expected_mult;
        for (const auto& r : expected_roots) {
            expected_mult[r.to_string()]++;
        }

        // Count found multiplicities
        std::map<std::string, int> found_mult;
        for (const auto& r : found_roots) {
            found_mult[r.to_string()]++;
        }

        // Check each expected root is found with correct multiplicity
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

    std::cout << "Property 7: Rational root completeness: " << pass_count
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
    test_degree_stops_at_4();

    // Property-based test
    auto result = test_property7_rational_root_completeness();

    if (result.trials_passed < result.trials_total) {
        std::cerr << "\nProperty 7 FAILED: " << result.trials_passed << "/"
                  << result.trials_total << " trials passed\n";
        if (!result.first_failure_msg.empty()) {
            std::cerr << "First failure: " << result.first_failure_msg << "\n";
        }
        return 1;
    }

    std::cout << "\nAll tests PASSED!\n";
    return 0;
}
