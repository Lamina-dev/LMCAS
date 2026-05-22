#include "solve_polynomial.hpp"
#include "polynomial.hpp"
#include "rational.hpp"
#include <iostream>
#include <cassert>
#include <algorithm>
#include <map>

using namespace lamina;

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
    assert(roots.size() == 1);
    assert(roots[0] == Rational(2));
    std::cout << "PASSED\n";
}

void test_quadratic_two_roots() {
    std::cout << "Test: quadratic (x-1)(x-3) ... ";

    auto p = poly_from_roots({Rational(1), Rational(3)});
    auto roots = find_rational_roots(p);
    assert(roots.size() == 2);
    assert(count_occurrences(roots, Rational(1)) == 1);
    assert(count_occurrences(roots, Rational(3)) == 1);
    std::cout << "PASSED\n";
}

void test_repeated_root() {
    std::cout << "Test: repeated root (x-2)^3 ... ";

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

    auto p = poly_from_roots({Rational(1, 2), Rational(3, 4)});
    auto roots = find_rational_roots(p);
    assert(roots.size() == 2);
    assert(count_occurrences(roots, Rational(1, 2)) == 1);
    assert(count_occurrences(roots, Rational(3, 4)) == 1);
    std::cout << "PASSED\n";
}

void test_zero_constant_term() {
    std::cout << "Test: zero constant term x^2(x-3) ... ";

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

    Polynomial<Rational> p1({Rational(2), Rational(-3), Rational(1)}, "x");
    Polynomial<Rational> p2({Rational(1), Rational(0), Rational(1)}, "x");
    Polynomial<Rational> p3({Rational(1), Rational(1), Rational(1)}, "x");
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

    Polynomial<Rational> p1({Rational(-1), Rational(1)}, "x");
    Polynomial<Rational> p2({Rational(-2), Rational(1)}, "x");
    Polynomial<Rational> p3({Rational(-3), Rational(1)}, "x");
    Polynomial<Rational> p4({Rational(1), Rational(0), Rational(1)}, "x");
    Polynomial<Rational> p5({Rational(1), Rational(1), Rational(1)}, "x");
    auto p = p1 * p2 * p3 * p4 * p5;
    auto roots = find_rational_roots(p);

    assert(roots.size() == 3);
    assert(count_occurrences(roots, Rational(1)) == 1);
    assert(count_occurrences(roots, Rational(2)) == 1);
    assert(count_occurrences(roots, Rational(3)) == 1);
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

PropertyTestResult test_property7_rational_root_completeness() {
    std::cout << "\n=== Property 7: Rational root completeness ===\n";
    std::cout << "**Validates: Requirements 4.1, 4.2**\n";

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
