#include "test_common.hpp"
#include "numeric_evaluation.hpp"
#include "solver.hpp"
#include "solve_strategies.hpp"
#include "newton_raphson.hpp"
#include <cmath>
#include <random>
#include <sstream>
#include <algorithm>
#include <optional>

using namespace LMCAS;

static std::optional<double> real_numeric_value(const std::shared_ptr<SymbolicExpr> &expr) {
    if (!expr) {
        return std::nullopt;
    }
    LMCAS::ComputationContext context;
    auto evaluated = LMCAS::evaluate_numeric(*expr, LMCAS::NumericBindings{}, context);
    if (!evaluated || !evaluated.value().is_finite() ||
        !std::isfinite(evaluated.value().value)) {
        return std::nullopt;
    }
    return evaluated.value().value;
}

static bool quadratic_trial_residuals_are_zero(
    const std::vector<std::shared_ptr<SymbolicExpr>> &sols,
    const std::shared_ptr<SymbolicExpr> &expr, const BigInt &disc,
    int a_val, int b_val, int c_val, int trial, double TOL) {
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
            ADD_FAILURE() << msg.str();
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
            ADD_FAILURE() << msg.str();
            trial_ok = false;
            break;
        }
    }

    return trial_ok;
}

static bool quadratic_duplicate_roots_are_valid(
    const std::vector<std::shared_ptr<SymbolicExpr>> &sols,
    const std::shared_ptr<SymbolicExpr> &expr, const BigInt &disc, double TOL) {
    if (disc == 0 && sols.size() == 2) {

        auto res1 = expr->substitute("x", sols[0])->simplify();
        auto res2 = expr->substitute("x", sols[1])->simplify();
        auto r1_val = real_numeric_value(res1);
        auto r2_val = real_numeric_value(res2);
        if (r1_val && r2_val && std::abs(*r1_val) < TOL && std::abs(*r2_val) < TOL) {
            return true;
        }
    }
    return false;
}

TEST(SolveExtended, SolveHigherDegreePolynomialRootof) {
    auto x = SymbolicExpr::variable("x");
    auto x3 = SymbolicExpr::power(x, SymbolicExpr::number(3));
    auto eq = SymbolicExpr::add(x3, SymbolicExpr::number(-2));

    auto sols = LMCAS::solve_finite_checked(eq, "x").value();
    EXPECT_TRUE((sols.size() == 3)) << "cubic x^3-2 should return 3 roots";
}

TEST(SolveExtended, SolveRationalSystemDenominatorFilter) {
    auto x = SymbolicExpr::variable("x");
    auto denom = SymbolicExpr::add(x, SymbolicExpr::number(-1));
    auto frac = SymbolicExpr::divide(x, denom);
    auto eq = SymbolicExpr::add(frac, SymbolicExpr::number(-2));

    std::vector<SymbolicExpr> eqs = {*eq};
    auto checked_solutions =
        LMCAS::Solver::solve_polynomial_system_checked(eqs, {"x"});
    ASSERT_TRUE(checked_solutions.has_value()) << "checked rational system solve succeeds";
    const auto &sols = checked_solutions.value();
    ASSERT_EQ(sols.size(), 1u) << "rational system solutions size";
    auto x_val = LMCAS::detail::make_expression_ptr(sols[0].at("x"));
    EXPECT_TRUE(test_same_expression(x_val, SymbolicExpr::number(2))) << "rational system x=2";
}

TEST(SolveExtended, DispatcherUnsupportedEquationIsInconclusive) {
    auto x = SymbolicExpr::variable("x");

    auto x_to_x = SymbolicExpr::power(x, x);

    auto eq = SymbolicExpr::add(SymbolicExpr::sin(x), x_to_x);

    LMCAS::SolveOptions opts;
    opts.allow_numeric = false;
    opts.return_rootof = true;

    auto result = solve_equation(eq, "x", opts);
    EXPECT_TRUE((!result && result.error().code == CasErrc::Inconclusive)) << "unsupported sin(x)+x^x with numeric disabled is Inconclusive";
}

TEST(SolveExtended, DispatcherAllowNumericFalseSkipsNumericalSolver) {
    auto x = SymbolicExpr::variable("x");

    auto x_to_x = SymbolicExpr::power(x, x);
    auto eq = SymbolicExpr::add(x_to_x, SymbolicExpr::number(-2));

    LMCAS::SolveOptions opts_no_numeric;
    opts_no_numeric.allow_numeric = false;
    auto result_no_numeric = solve_equation(eq, "x", opts_no_numeric);
    EXPECT_TRUE((!result_no_numeric && result_no_numeric.error().code == CasErrc::Inconclusive)) << "allow_numeric=false on x^x-2 is Inconclusive";

    LMCAS::SolveOptions opts_numeric;
    opts_numeric.allow_numeric = true;
    opts_numeric.has_initial_guess = true;
    opts_numeric.initial_guess = 1.5;
    opts_numeric.tolerance = 1e-10;
    auto sols_numeric = solve_numeric_checked(eq, "x", opts_numeric);
    ASSERT_TRUE((sols_numeric.has_value())) << (sols_numeric ? "numeric candidate search succeeds" : sols_numeric.error().message);
    if (!sols_numeric) {
        return;
    }
    EXPECT_TRUE((!sols_numeric.value().empty())) << "numeric API retains supported x^x-2 candidate search";
    for (const auto &root : sols_numeric.value()) {
        EXPECT_TRUE((std::isfinite(root.value) && root.value > 0)) << "numeric candidate is in the real expression domain";
        {
            const double actual_value = (std::pow(root.value, root.value));
            const double expected_value = (2.0);
            const double tolerance = (opts_numeric.tolerance);
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_GE(tolerance, 0.0);
            EXPECT_NEAR(actual_value, expected_value, tolerance);
        }
    }
    auto complete = solve_finite_checked(eq, "x", opts_numeric);
    EXPECT_TRUE((!complete && complete.error().code == CasErrc::Inconclusive)) << "an initial-guess search cannot establish a globally complete finite set";
}

TEST(SolveExtended, DispatcherReturnRootofFalseSuppressesRootofEmission) {
    auto x = SymbolicExpr::variable("x");
    auto x5 = SymbolicExpr::power(x, SymbolicExpr::number(5));

    auto eq = SymbolicExpr::add(
        SymbolicExpr::add(x5, SymbolicExpr::multiply(x, SymbolicExpr::number(-1))),
        SymbolicExpr::number(-1));

    LMCAS::SolveOptions opts_rootof;
    opts_rootof.return_rootof = true;
    opts_rootof.allow_numeric = false;
    auto sols_rootof = solve_vector_for_test(eq, "x", opts_rootof);
    EXPECT_TRUE((sols_rootof.size() == 5)) << "return_rootof=true -> 5 RootOf solutions for degree-5";
    if (!sols_rootof.empty()) {

        {
            const std::string actual_text = (sols_rootof[0]->to_string());
            for (const auto &token : std::vector<std::string>{"rootof"}) {
                EXPECT_NE(actual_text.find(token), std::string::npos) << "return_rootof=true produces RootOf expressions" << ": missing " << token << " in " << actual_text;
            }
        }
    }

    LMCAS::SolveOptions opts_no_rootof;
    opts_no_rootof.return_rootof = false;
    opts_no_rootof.allow_numeric = false;
    auto result_no_rootof = solve_equation(eq, "x", opts_no_rootof);
    EXPECT_TRUE((!result_no_rootof && result_no_rootof.error().code == CasErrc::Inconclusive)) << "x^5-x-1 is Inconclusive when numeric and RootOf are disabled";
}

TEST(SolveExtended, DispatcherSimplificationConvertsFXGXToFXGX0) {
    auto x = SymbolicExpr::variable("x");
    auto lhs = SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), x), SymbolicExpr::number(3));
    auto rhs = SymbolicExpr::add(x, SymbolicExpr::number(5));
    auto eq = SymbolicExpr::eq(lhs, rhs);

    LMCAS::SolveOptions opts;
    auto sols = solve_vector_for_test(eq, "x", opts);
    EXPECT_TRUE((sols.size() == 1)) << "f(x)=g(x) form produces one solution";
    if (!sols.empty()) {

        auto val = sols[0]->simplify();
        auto num_val = val->to_numeric();
        bool close_to_2 = std::abs(num_val - 2.0) < 1e-10;
        EXPECT_TRUE((close_to_2)) << "f(x)=g(x) preprocessing: 2x+3=x+5 gives x=2";
    }
}

TEST(SolveExtended, DispatcherDegree0NonZeroConstantReturnsEmpty) {
    auto five = SymbolicExpr::number(5);

    LMCAS::SolveOptions opts;
    auto result = solve_equation(five, "x", opts);
    EXPECT_TRUE((result && std::holds_alternative<EmptySolutions>(result.value()))) << "degree-0 non-zero constant has a successful empty solution set";
}

TEST(SolveExtended, DispatcherDegree0NonZeroConstantViaEquationFormReturnsEmpty) {
    auto three = SymbolicExpr::number(3);
    auto zero = SymbolicExpr::number(0);
    auto eq = SymbolicExpr::eq(three, zero);

    LMCAS::SolveOptions opts;
    auto result = solve_equation(eq, "x", opts);
    EXPECT_TRUE((result && std::holds_alternative<EmptySolutions>(result.value()))) << "3=0 has a successful empty solution set";
}

TEST(SolveExtended, RootCountInvariantAcrossStrategies) {
    const int NUM_TRIALS = 60;

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

        auto solutions = LMCAS::solve_finite_checked(expanded, "x").value();

        EXPECT_EQ((int)solutions.size(), degree)
            << "Trial " << trial << " degree=" << degree
            << " roots=" << ::testing::PrintToString(roots);
    }
}

TEST(SolveExtended, PartALinearBackwardCompatibility) {
    const int NUM_LINEAR_TRIALS = 40;
    const double TOL = 1e-10;

    std::mt19937 rng_lin(7777);
    std::uniform_int_distribution<int> coeff_dist(-20, 20);

    for (int trial = 0; trial < NUM_LINEAR_TRIALS; ++trial) {
        int a_val = coeff_dist(rng_lin);
        while (a_val == 0)
            a_val = coeff_dist(rng_lin);
        int b_val = coeff_dist(rng_lin);

        auto x = SymbolicExpr::variable("x");
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(a_val), x),
            SymbolicExpr::number(b_val));

        auto sols = LMCAS::solve_finite_checked(expr, "x").value();

        if (sols.size() != 1) {
            std::ostringstream msg;
            msg << "Linear Trial " << trial
                << " (a=" << a_val << ", b=" << b_val
                << "): expected 1 root, got " << sols.size();
            ADD_FAILURE() << msg.str();
            continue;
        }

        double expected_root = -(double)b_val / (double)a_val;
        double actual_root = sols[0]->to_numeric();

        EXPECT_TRUE(std::abs(actual_root - expected_root) < TOL)
            << "Linear Trial " << trial
            << " (a=" << a_val << ", b=" << b_val
            << "): expected root " << expected_root << ", got " << actual_root;
    }
}

TEST(SolveExtended, PartBQuadraticBackwardCompatibility) {
    const int NUM_QUAD_TRIALS = 40;
    const double TOL = 1e-10;

    std::mt19937 rng_quad(8888);
    std::uniform_int_distribution<int> coeff_dist(-10, 10);

    for (int trial = 0; trial < NUM_QUAD_TRIALS; ++trial) {
        int a_val = coeff_dist(rng_quad);
        while (a_val == 0)
            a_val = coeff_dist(rng_quad);
        int b_val = coeff_dist(rng_quad);
        int c_val = coeff_dist(rng_quad);

        const BigInt disc =
            BigInt(b_val) * BigInt(b_val) - BigInt(4) * BigInt(a_val) * BigInt(c_val);

        auto x = SymbolicExpr::variable("x");
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(a_val), x2),
            SymbolicExpr::add(
                SymbolicExpr::multiply(SymbolicExpr::number(b_val), x),
                SymbolicExpr::number(c_val)));

        auto sols = LMCAS::solve_finite_checked(expr, "x").value();

        size_t expected_count = (disc == 0) ? 1 : 2;

        if (sols.size() != expected_count) {

            EXPECT_TRUE(quadratic_duplicate_roots_are_valid(sols, expr, disc, TOL))
                << "Quadratic Trial " << trial
                << " (a=" << a_val << ", b=" << b_val << ", c=" << c_val
                << ", disc=" << disc.to_string()
                << "): expected " << expected_count << " roots, got " << sols.size();
            continue;
        }

        bool trial_ok = quadratic_trial_residuals_are_zero(
            sols, expr, disc, a_val, b_val, c_val, trial, TOL);

        EXPECT_TRUE(trial_ok)
            << "Quadratic Trial " << trial << " (a=" << a_val
            << ", b=" << b_val << ", c=" << c_val
            << ", disc=" << disc.to_string() << ")";
    }
}
