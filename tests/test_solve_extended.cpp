#include "test_common.hpp"
#include "numeric_evaluation.hpp"
#include "solver.hpp"
#include "solve_strategies.hpp"
#include <cmath>
#include <random>
#include <sstream>
#include <algorithm>
#include <optional>

static std::optional<double> real_numeric_value(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr) return std::nullopt;
    lamina::ComputationContext context;
    auto evaluated = lamina::evaluate_numeric(*expr, lamina::NumericBindings{}, context);
    if (!evaluated || !evaluated.value().is_finite() ||
        !std::isfinite(evaluated.value().value)) {
        return std::nullopt;
    }
    return evaluated.value().value;
}

int main() {

    TEST_CASE("Solve Higher Degree Polynomial (RootOf)");
    {

        auto x = SymbolicExpr::variable("x");
        auto x3 = SymbolicExpr::power(x, SymbolicExpr::number(3));
        auto eq = SymbolicExpr::add(x3, SymbolicExpr::number(-2));

        auto sols = SymbolicExpr::solve(eq, "x");
        EXPECT_TRUE(sols.size() == 3, "cubic x^3-2 should return 3 roots");
    }

    TEST_CASE("Solve Linear+Exp (LambertW)");
    {

        auto x = SymbolicExpr::variable("x");
        auto eq = SymbolicExpr::add(x, SymbolicExpr::exp(x));

        auto sols = SymbolicExpr::solve(eq, "x");

        EXPECT_TRUE(!sols.empty(), "x+exp(x)=0 should have a solution");
    }

    TEST_CASE("Solve Rational System (Denominator Filter)");
    {
        auto x = SymbolicExpr::variable("x");
        auto denom = SymbolicExpr::add(x, SymbolicExpr::number(-1));
        auto frac = SymbolicExpr::divide(x, denom);
        auto eq = SymbolicExpr::add(frac, SymbolicExpr::number(-2));

        std::vector<SymbolicExpr> eqs = {*eq};
        auto sols = lamina::Solver::solve_polynomial_system(eqs, {"x"});
        EXPECT_TRUE(sols.size() == 1, "rational system solutions size");
        if (!sols.empty()) {
            auto x_val = lamina::detail::make_expression_ptr(sols[0].at("x"));
            EXPECT_EQ_EXPR(x_val, SymbolicExpr::number(2), "rational system x=2");
        }
    }

    TEST_CASE("Dispatcher: all strategies return empty -> empty result without exception");
    {

        auto x = SymbolicExpr::variable("x");

        auto x_to_x = SymbolicExpr::power(x, x);

        auto eq = SymbolicExpr::add(SymbolicExpr::sin(x), x_to_x);

        lamina::SolveOptions opts;
        opts.allow_numeric = false;
        opts.return_rootof = true;

        auto sols = lamina::solve_dispatch(eq, "x", opts);
        EXPECT_TRUE(sols.empty(), "all strategies fail -> empty result");

        std::cout << "[PASS] no exception thrown on total fallthrough" << std::endl;
    }

    TEST_CASE("Dispatcher: allow_numeric=false skips Numerical_Solver");
    {

        auto x = SymbolicExpr::variable("x");

        auto x_to_x = SymbolicExpr::power(x, x);
        auto eq = SymbolicExpr::add(x_to_x, SymbolicExpr::number(-2));

        lamina::SolveOptions opts_no_numeric;
        opts_no_numeric.allow_numeric = false;
        auto sols_no_numeric = lamina::solve_dispatch(eq, "x", opts_no_numeric);
        EXPECT_TRUE(sols_no_numeric.empty(), "allow_numeric=false -> no solutions for x^x-2");

        lamina::SolveOptions opts_numeric;
        opts_numeric.allow_numeric = true;
        opts_numeric.has_initial_guess = true;
        opts_numeric.initial_guess = 1.5;
        opts_numeric.tolerance = 1e-10;
        auto sols_numeric = lamina::solve_dispatch(eq, "x", opts_numeric);

        EXPECT_TRUE(!sols_numeric.empty(), "allow_numeric=true -> finds numeric solution for x^x-2");
    }

    TEST_CASE("Dispatcher: return_rootof=false suppresses RootOf emission");
    {

        auto x = SymbolicExpr::variable("x");
        auto x5 = SymbolicExpr::power(x, SymbolicExpr::number(5));

        auto eq = SymbolicExpr::add(
            SymbolicExpr::add(x5, SymbolicExpr::multiply(x, SymbolicExpr::number(-1))),
            SymbolicExpr::number(-1));

        lamina::SolveOptions opts_rootof;
        opts_rootof.return_rootof = true;
        opts_rootof.allow_numeric = false;
        auto sols_rootof = lamina::solve_dispatch(eq, "x", opts_rootof);
        EXPECT_TRUE(sols_rootof.size() == 5, "return_rootof=true -> 5 RootOf solutions for degree-5");
        if (!sols_rootof.empty()) {

            EXPECT_CONTAINS(sols_rootof[0]->to_string(), {"rootof"}, "return_rootof=true produces RootOf expressions");
        }

        lamina::SolveOptions opts_no_rootof;
        opts_no_rootof.return_rootof = false;
        opts_no_rootof.allow_numeric = false;
        auto sols_no_rootof = lamina::solve_dispatch(eq, "x", opts_no_rootof);

        bool has_rootof = false;
        for (const auto& sol : sols_no_rootof) {
            if (sol && sol->to_string().find("rootof") != std::string::npos) {
                has_rootof = true;
            }
        }
        EXPECT_TRUE(!has_rootof, "return_rootof=false suppresses RootOf emission");
    }

    TEST_CASE("Dispatcher: simplification converts f(x)=g(x) to f(x)-g(x)=0");
    {

        auto x = SymbolicExpr::variable("x");
        auto lhs = SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), x), SymbolicExpr::number(3));
        auto rhs = SymbolicExpr::add(x, SymbolicExpr::number(5));
        auto eq = SymbolicExpr::eq(lhs, rhs);

        lamina::SolveOptions opts;
        auto sols = lamina::solve_dispatch(eq, "x", opts);
        EXPECT_TRUE(sols.size() == 1, "f(x)=g(x) form produces one solution");
        if (!sols.empty()) {

            auto val = sols[0]->simplify();
            auto num_val = val->to_numeric();
            bool close_to_2 = std::abs(num_val - 2.0) < 1e-10;
            EXPECT_TRUE(close_to_2, "f(x)=g(x) preprocessing: 2x+3=x+5 gives x=2");
        }
    }

    TEST_CASE("Dispatcher: degree-0 non-zero constant returns empty");
    {

        auto five = SymbolicExpr::number(5);

        lamina::SolveOptions opts;
        auto sols = lamina::solve_dispatch(five, "x", opts);
        EXPECT_TRUE(sols.empty(), "degree-0 non-zero constant -> empty result");
    }

    TEST_CASE("Dispatcher: degree-0 non-zero constant via equation form returns empty");
    {

        auto three = SymbolicExpr::number(3);
        auto zero = SymbolicExpr::number(0);
        auto eq = SymbolicExpr::eq(three, zero);

        lamina::SolveOptions opts;
        auto sols = lamina::solve_dispatch(eq, "x", opts);
        EXPECT_TRUE(sols.empty(), "3=0 equation -> empty (no solution exists)");
    }

    TEST_CASE("Root count invariant across strategies");
    {
        const int NUM_TRIALS = 60;
        int pass_count = 0;

        std::mt19937 rng(7777);
        std::uniform_int_distribution<int> degree_dist(1, 8);
        std::uniform_int_distribution<int> root_dist(-5, 5);

        for (int trial = 0; trial < NUM_TRIALS; ++trial) {
            int degree = degree_dist(rng);

            std::vector<int> roots;
            roots.reserve(degree);
            for (int i = 0; i < degree; ++i) {
                roots.push_back(root_dist(rng));
            }

            auto x = SymbolicExpr::variable("x");

            auto poly_expr = SymbolicExpr::add(x, SymbolicExpr::number(-roots[0]));
            for (int i = 1; i < degree; ++i) {
                auto factor = SymbolicExpr::add(x, SymbolicExpr::number(-roots[i]));
                poly_expr = SymbolicExpr::multiply(poly_expr, factor);
            }

            auto expanded = poly_expr->expand();

            auto solutions = SymbolicExpr::solve(expanded, "x");

            if ((int)solutions.size() == degree) {
                pass_count++;
            } else {
                std::ostringstream msg;
                msg << "Trial " << trial << " degree=" << degree
                    << " roots=[";
                for (int i = 0; i < degree; ++i) {
                    if (i > 0) msg << ",";
                    msg << roots[i];
                }
                msg << "]: expected " << degree << " solutions, got "
                    << solutions.size();
                EXPECT_TRUE(false, msg.str());
            }
        }

        {
            std::ostringstream msg;
            msg << "Root count invariant: " << pass_count
                << "/" << NUM_TRIALS << " trials passed";
            EXPECT_TRUE(pass_count == NUM_TRIALS, msg.str());
        }
    }

    TEST_CASE("(Part A): Linear backward compatibility");
    {
        const int NUM_LINEAR_TRIALS = 40;
        const double TOL = 1e-10;
        int linear_pass_count = 0;

        std::mt19937 rng_lin(7777);
        std::uniform_int_distribution<int> coeff_dist(-20, 20);

        for (int trial = 0; trial < NUM_LINEAR_TRIALS; ++trial) {
            int a_val = coeff_dist(rng_lin);
            while (a_val == 0) a_val = coeff_dist(rng_lin);
            int b_val = coeff_dist(rng_lin);

            auto x = SymbolicExpr::variable("x");
            auto expr = SymbolicExpr::add(
                SymbolicExpr::multiply(SymbolicExpr::number(a_val), x),
                SymbolicExpr::number(b_val)
            );

            auto sols = SymbolicExpr::solve(expr, "x");

            if (sols.size() != 1) {
                std::ostringstream msg;
                msg << "Linear Trial " << trial
                    << " (a=" << a_val << ", b=" << b_val
                    << "): expected 1 root, got " << sols.size();
                EXPECT_TRUE(false, msg.str());
                continue;
            }

            double expected_root = -(double)b_val / (double)a_val;
            double actual_root = sols[0]->to_numeric();

            if (std::abs(actual_root - expected_root) < TOL) {
                linear_pass_count++;
            } else {
                std::ostringstream msg;
                msg << "Linear Trial " << trial
                    << " (a=" << a_val << ", b=" << b_val
                    << "): expected root " << expected_root
                    << ", got " << actual_root;
                EXPECT_TRUE(false, msg.str());
            }
        }

        {
            std::ostringstream msg;
            msg << "Linear: " << linear_pass_count << "/"
                << NUM_LINEAR_TRIALS << " trials passed";
            EXPECT_TRUE(linear_pass_count == NUM_LINEAR_TRIALS, msg.str());
        }
    }

    TEST_CASE("(Part B): Quadratic backward compatibility");
    {
        const int NUM_QUAD_TRIALS = 40;
        const double TOL = 1e-10;
        int quad_pass_count = 0;

        std::mt19937 rng_quad(8888);
        std::uniform_int_distribution<int> coeff_dist(-10, 10);

        for (int trial = 0; trial < NUM_QUAD_TRIALS; ++trial) {
            int a_val = coeff_dist(rng_quad);
            while (a_val == 0) a_val = coeff_dist(rng_quad);
            int b_val = coeff_dist(rng_quad);
            int c_val = coeff_dist(rng_quad);

            long long disc = (long long)b_val * b_val - 4LL * a_val * c_val;

            auto x = SymbolicExpr::variable("x");
            auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
            auto expr = SymbolicExpr::add(
                SymbolicExpr::multiply(SymbolicExpr::number(a_val), x2),
                SymbolicExpr::add(
                    SymbolicExpr::multiply(SymbolicExpr::number(b_val), x),
                    SymbolicExpr::number(c_val)
                )
            );

            auto sols = SymbolicExpr::solve(expr, "x");

            size_t expected_count = (disc == 0) ? 1 : 2;

            if (sols.size() != expected_count) {

                if (disc == 0 && sols.size() == 2) {

                    auto res1 = expr->substitute("x", sols[0])->simplify();
                    auto res2 = expr->substitute("x", sols[1])->simplify();
                    auto r1_val = real_numeric_value(res1);
                    auto r2_val = real_numeric_value(res2);
                    if (r1_val && r2_val && std::abs(*r1_val) < TOL && std::abs(*r2_val) < TOL) {
                        quad_pass_count++;
                        continue;
                    }
                }
                std::ostringstream msg;
                msg << "Quadratic Trial " << trial
                    << " (a=" << a_val << ", b=" << b_val << ", c=" << c_val
                    << ", disc=" << disc
                    << "): expected " << expected_count << " roots, got " << sols.size();
                EXPECT_TRUE(false, msg.str());
                continue;
            }

            bool trial_ok = true;
            for (size_t i = 0; i < sols.size(); ++i) {

                auto residual_expr = expr->substitute("x", sols[i])->simplify();
                auto maybe_residual = real_numeric_value(residual_expr);

                if (disc < 0 && (!maybe_residual || std::abs(*maybe_residual) < 1e-6)) {

                    continue;
                }

                if (!maybe_residual) {
                    std::ostringstream msg;
                    msg << "Quadratic Trial " << trial
                        << " root " << i << ": residual is not real-numerically evaluable"
                        << " (a=" << a_val << ", b=" << b_val << ", c=" << c_val << ")";
                    EXPECT_TRUE(false, msg.str());
                    trial_ok = false;
                    break;
                }

                double residual = *maybe_residual;
                if (std::abs(residual) >= TOL) {
                    std::ostringstream msg;
                    msg << "Quadratic Trial " << trial
                        << " root " << i << ": |f(r)| = " << std::abs(residual)
                        << " >= 1e-10"
                        << " (a=" << a_val << ", b=" << b_val << ", c=" << c_val << ")";
                    EXPECT_TRUE(false, msg.str());
                    trial_ok = false;
                    break;
                }
            }

            if (trial_ok) {
                quad_pass_count++;
            }
        }

        {
            std::ostringstream msg;
            msg << "Quadratic: " << quad_pass_count << "/"
                << NUM_QUAD_TRIALS << " trials passed";
            EXPECT_TRUE(quad_pass_count == NUM_QUAD_TRIALS, msg.str());
        }
    }

    return TEST_REPORT();
}
