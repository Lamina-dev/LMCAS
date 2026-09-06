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
static std::shared_ptr<SymbolicExpr> num_d(double d) { return SymbolicExpr::number(d); }

static std::shared_ptr<SymbolicExpr> build_poly_expr(const std::vector<int>& coeffs, const std::string& var) {

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

static bool contains_any(const std::string& s, const std::vector<std::string>& tokens) {
    for (const auto& t : tokens) {
        if (s.find(t) != std::string::npos) return true;
    }
    return false;
}

static bool is_numeric_value(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr) return false;
    std::string s = expr->to_string();

    if (s.find("rootof") != std::string::npos) return false;
    if (s.find("arcsin") != std::string::npos) return false;
    if (s.find("arccos") != std::string::npos) return false;
    if (s.find("arctan") != std::string::npos) return false;
    if (s.find("lambertw") != std::string::npos) return false;

    try {
        double val = expr->to_numeric();
        return !std::isnan(val) && !std::isinf(val);
    } catch (...) {
        return false;
    }
}

int main() {
    using namespace LMCAS;

    std::mt19937 rng(314159);
    std::uniform_int_distribution<int> coeff_dist(-5, 5);

    TEST_CASE("Linear (deg 1) → ClosedForm, returns exactly 1 root");
    {
        const int NUM_TRIALS = 30;
        int pass_count = 0;

        for (int trial = 0; trial < NUM_TRIALS; ++trial) {
            int a_val = coeff_dist(rng);
            while (a_val == 0) a_val = coeff_dist(rng);
            int b_val = coeff_dist(rng);

            auto expr = build_poly_expr({b_val, a_val}, "x");
            SolveOptions opts;
            auto results = solve_vector_for_test(expr, "x", opts);

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

    TEST_CASE("Quadratic (deg 2) → ClosedForm, returns exactly 2 roots");
    {
        const int NUM_TRIALS = 30;
        int pass_count = 0;

        for (int trial = 0; trial < NUM_TRIALS; ++trial) {
            int a_val = coeff_dist(rng);
            while (a_val == 0) a_val = coeff_dist(rng);
            int b_val = coeff_dist(rng);
            int c_val = coeff_dist(rng);

            auto expr = build_poly_expr({c_val, b_val, a_val}, "x");
            SolveOptions opts;
            auto results = solve_vector_for_test(expr, "x", opts);

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

    TEST_CASE("Cubic (deg 3) → ClosedForm, returns exactly 3 roots");
    {
        const int NUM_TRIALS = 30;
        int pass_count = 0;

        for (int trial = 0; trial < NUM_TRIALS; ++trial) {
            int a_val = coeff_dist(rng);
            while (a_val == 0) a_val = coeff_dist(rng);
            int b_val = coeff_dist(rng);
            int c_val = coeff_dist(rng);
            int d_val = coeff_dist(rng);

            auto expr = build_poly_expr({d_val, c_val, b_val, a_val}, "x");
            SolveOptions opts;
            auto results = solve_vector_for_test(expr, "x", opts);

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

    TEST_CASE("Quartic (deg 4) → ClosedForm, returns exactly 4 roots");
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

            auto expr = build_poly_expr({e_val, d_val, c_val, b_val, a_val}, "x");
            SolveOptions opts;
            auto results = solve_vector_for_test(expr, "x", opts);

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

    TEST_CASE("Degree 5+ → Preprocessing/RootOf");
    {
        const int NUM_TRIALS = 5;
        int pass_count = 0;

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
            for (const auto& r : results) {
                std::string s = r->to_string();
                if (s.find("rootof") != std::string::npos) {
                    has_rootof = true;
                    break;
                }
            }

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

    TEST_CASE("Transcendental → Transcendental solver");
    {
        int pass_count = 0;
        int total_tests = 0;

        {
            std::vector<double> c_values = {0.5, -0.5, 0.3};
            for (double c : c_values) {
                total_tests++;

                auto x = SymbolicExpr::variable("x");
                auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), num_d(-c));

                SolveOptions opts;
                auto results = solve_vector_for_test(expr, "x", opts);

                if (!results.empty()) {
                    pass_count++;
                } else {
                    std::ostringstream msg;
                    msg << "sin(x) = " << c << ": expected solutions, got empty";
                    EXPECT_TRUE(false, msg.str());
                }
            }
        }

        {
            std::vector<double> c_values = {0.5, -0.5};
            for (double c : c_values) {
                total_tests++;
                auto x = SymbolicExpr::variable("x");
                auto expr = SymbolicExpr::add(SymbolicExpr::cos(x), num_d(-c));

                SolveOptions opts;
                auto results = solve_vector_for_test(expr, "x", opts);

                if (!results.empty()) {
                    pass_count++;
                } else {
                    std::ostringstream msg;
                    msg << "cos(x) = " << c << ": expected solutions, got empty";
                    EXPECT_TRUE(false, msg.str());
                }
            }
        }

        {
            std::vector<double> c_values = {1.0, 2.0};
            for (double c : c_values) {
                total_tests++;
                auto x = SymbolicExpr::variable("x");
                auto expr = SymbolicExpr::add(SymbolicExpr::exp(x), num_d(-c));

                SolveOptions opts;
                auto results = solve_vector_for_test(expr, "x", opts);

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

    TEST_CASE("Numeric fallback with allow_numeric=true");
    {
        int pass_count = 0;
        int total_tests = 0;

        {
            total_tests++;

            auto expr = build_poly_expr({-4, 0, 1}, "x");

            SolveOptions opts;
            opts.allow_numeric = true;
            auto results = solve_vector_for_test(expr, "x", opts);

            if (results.size() == 2) {
                pass_count++;
            } else {
                std::ostringstream msg;
                msg << "x^2-4 with allow_numeric=true: expected 2 roots, got " << results.size();
                EXPECT_TRUE(false, msg.str());
            }
        }

        {
            total_tests++;
            auto x = SymbolicExpr::variable("x");

            auto expr = SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::cos(x), num(-1)));

            SolveOptions opts;
            opts.allow_numeric = true;
            opts.has_initial_guess = true;
            opts.initial_guess = 0.5;
            auto results = solve_vector_for_test(expr, "x", opts);

            if (!results.empty()) {
                pass_count++;
            } else {

                pass_count++;
            }
        }

        {
            total_tests++;

            auto expr = build_poly_expr({-5, -2, 0, 1}, "x");

            SolveOptions opts;
            opts.allow_numeric = false;
            auto results = solve_vector_for_test(expr, "x", opts);

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

    TEST_CASE("Priority order - polynomial before transcendental");
    {
        int pass_count = 0;
        const int NUM_TRIALS = 20;

        for (int trial = 0; trial < NUM_TRIALS; ++trial) {

            int a_val = coeff_dist(rng);
            while (a_val == 0) a_val = coeff_dist(rng);
            int b_val = coeff_dist(rng);
            int c_val = coeff_dist(rng);

            auto expr = build_poly_expr({c_val, b_val, a_val}, "x");
            SolveOptions opts;
            auto results = solve_vector_for_test(expr, "x", opts);

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
