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

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

static std::shared_ptr<SymbolicExpr> build_poly_expr(const std::vector<int> &coeffs, const std::string &var) {

    auto x = SymbolicExpr::variable(var);
    std::shared_ptr<SymbolicExpr> result = nullptr;

    for (size_t i = 0; i < coeffs.size(); ++i) {
        if (coeffs[i] == 0) {
            continue;
        }
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
    if (!result) {
        result = num(0);
    }
    return result;
}

static bool contains_any(const std::string &s, const std::vector<std::string> &tokens) {
    for (const auto &t : tokens) {
        if (s.find(t) != std::string::npos) {
            return true;
        }
    }
    return false;
}

TEST(DispatcherRouting, LinearDeg1ClosedformReturnsExactly1Root) {
    std::mt19937 rng(314159);
    std::uniform_int_distribution<int> coeff_dist(-5, 5);

    const int NUM_TRIALS = 30;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int a_val = coeff_dist(rng);
        while (a_val == 0)
            a_val = coeff_dist(rng);
        int b_val = coeff_dist(rng);

        auto expr = build_poly_expr({b_val, a_val}, "x");
        SolveOptions opts;
        auto results = solve_vector_for_test(expr, "x", opts);

        EXPECT_EQ(results.size(), 1u)
            << "Trial " << trial << " (a=" << a_val << ", b=" << b_val
            << "): expected 1 root, got " << results.size();
    }
}

TEST(DispatcherRouting, QuadraticDeg2ClosedformReturnsExactly2Roots) {
    std::mt19937 rng(314159);
    std::uniform_int_distribution<int> coeff_dist(-5, 5);

    const int NUM_TRIALS = 30;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int a_val = coeff_dist(rng);
        while (a_val == 0)
            a_val = coeff_dist(rng);
        int b_val = coeff_dist(rng);
        int c_val = coeff_dist(rng);

        auto expr = build_poly_expr({c_val, b_val, a_val}, "x");
        SolveOptions opts;
        auto results = solve_vector_for_test(expr, "x", opts);

        EXPECT_EQ(results.size(), 2u)
            << "Trial " << trial << " (a=" << a_val << ", b=" << b_val
            << ", c=" << c_val << "): expected 2 roots, got " << results.size();
    }
}

TEST(DispatcherRouting, CubicDeg3ClosedformReturnsExactly3Roots) {
    std::mt19937 rng(314159);
    std::uniform_int_distribution<int> coeff_dist(-5, 5);

    const int NUM_TRIALS = 30;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int a_val = coeff_dist(rng);
        while (a_val == 0)
            a_val = coeff_dist(rng);
        int b_val = coeff_dist(rng);
        int c_val = coeff_dist(rng);
        int d_val = coeff_dist(rng);

        auto expr = build_poly_expr({d_val, c_val, b_val, a_val}, "x");
        SolveOptions opts;
        auto results = solve_vector_for_test(expr, "x", opts);

        EXPECT_EQ(results.size(), 3u)
            << "Trial " << trial << " (a=" << a_val << ", b=" << b_val
            << ", c=" << c_val << ", d=" << d_val
            << "): expected 3 roots, got " << results.size();
    }
}

TEST(DispatcherRouting, QuarticDeg4ClosedformReturnsExactly4Roots) {
    std::mt19937 rng(314159);
    std::uniform_int_distribution<int> coeff_dist(-5, 5);

    const int NUM_TRIALS = 30;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int a_val = coeff_dist(rng);
        while (a_val == 0)
            a_val = coeff_dist(rng);
        int b_val = coeff_dist(rng);
        int c_val = coeff_dist(rng);
        int d_val = coeff_dist(rng);
        int e_val = coeff_dist(rng);

        auto expr = build_poly_expr({e_val, d_val, c_val, b_val, a_val}, "x");
        SolveOptions opts;
        auto results = solve_vector_for_test(expr, "x", opts);

        EXPECT_EQ(results.size(), 4u)
            << "Trial " << trial << " (a=" << a_val << ", b=" << b_val
            << ", c=" << c_val << ", d=" << d_val << ", e=" << e_val
            << "): expected 4 roots, got " << results.size();
    }
}

TEST(DispatcherRouting, QuarticIrreducibleReturnsExactly4Roots) {
    auto expr = build_poly_expr({-3, -4, -5, 0, -3}, "x");
    auto results = solve_vector_for_test(expr, "x");

    ASSERT_EQ(results.size(), 4u);
    for (const auto &root : results) {
        EXPECT_NE(root->to_string().find("rootof"), std::string::npos);
    }
}

TEST(DispatcherRouting, Degree5PreprocessingRootof) {
    const int NUM_TRIALS = 5;

    std::vector<std::vector<int>> test_polys = {
        {1, 0, 0, 0, 0, 1},
        {-1, 1, 0, 0, 0, 1},
        {2, 0, 0, 0, 0, 1},
        {0, 0, -1, 0, 0, 1},
        {-1, 0, 0, 0, 1, 1},
    };

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        auto expr = build_poly_expr(test_polys[trial], "x");
        SolveOptions opts;
        opts.return_rootof = true;
        auto results = solve_vector_for_test(expr, "x", opts);

        int degree = (int)test_polys[trial].size() - 1;

        bool correct_count = ((int)results.size() == degree);

        bool has_rootof = false;
        for (const auto &r : results) {
            std::string s = r->to_string();
            if (s.find("rootof") != std::string::npos) {
                has_rootof = true;
                break;
            }
        }

        EXPECT_TRUE(correct_count || has_rootof)
            << "Trial " << trial << " (deg=" << degree
            << "): expected " << degree << " roots, got " << results.size()
            << ", has_rootof=" << has_rootof;
    }
}

TEST(DispatcherRouting, NumericFallbackWithAllowNumericTrue) {
    {

        auto expr = build_poly_expr({-4, 0, 1}, "x");

        SolveOptions opts;
        opts.allow_numeric = true;
        auto results = solve_vector_for_test(expr, "x", opts);

        EXPECT_EQ(results.size(), 2u)
            << "x^2-4 with allow_numeric=true: expected 2 roots, got " << results.size();
    }

    {

        auto expr = build_poly_expr({-5, -2, 0, 1}, "x");

        SolveOptions opts;
        opts.allow_numeric = false;
        auto results = solve_vector_for_test(expr, "x", opts);

        EXPECT_EQ(results.size(), 3u)
            << "x^3-2x-5 with allow_numeric=false: expected 3 roots, got " << results.size();
    }
}

TEST(DispatcherRouting, PriorityOrderPolynomialBeforeTranscendental) {
    std::mt19937 rng(314159);
    std::uniform_int_distribution<int> coeff_dist(-5, 5);

    const int NUM_TRIALS = 20;

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        SCOPED_TRACE(::testing::Message() << "Trial " << trial);

        int a_val = coeff_dist(rng);
        while (a_val == 0)
            a_val = coeff_dist(rng);
        int b_val = coeff_dist(rng);
        int c_val = coeff_dist(rng);
        SCOPED_TRACE(::testing::Message() << "a=" << a_val << ", b=" << b_val
                                          << ", c=" << c_val);

        auto expr = build_poly_expr({c_val, b_val, a_val}, "x");
        SolveOptions opts;
        auto results = solve_vector_for_test(expr, "x", opts);

        EXPECT_EQ(results.size(), 2u) << "polynomials correctly prioritized over transcendental";
        if (results.size() == 2) {
            bool has_transcendental_token = false;
            for (const auto &r : results) {
                std::string s = r->to_string();
                if (contains_any(s, {"arcsin", "arccos", "arctan", "lambertw"})) {
                    has_transcendental_token = true;
                    break;
                }
            }
            EXPECT_FALSE(has_transcendental_token)
                << "Trial " << trial << ": polynomial routed to transcendental solver";
        }
    }
}
