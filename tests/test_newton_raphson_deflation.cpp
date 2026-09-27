#include "test_newton_raphson_support.hpp"

static void test_polynomial_numeric_root_endpoints(
    const std::shared_ptr<SymbolicExpr> &x) {
    {

        auto poly_f = SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::number(-4));

        LMCAS::SolveOptions opts;
        opts.allow_numeric = true;
        opts.tolerance = 1e-10;
        opts.max_newton_iterations = 100;

        auto roots = LMCAS::solve_numeric_checked(poly_f, "x", opts).value();

        EXPECT_TRUE((roots.size() == 2)) << "Polynomial x^2-4 should find 2 roots via Sturm path";
        for (const auto &root : roots) {
            EXPECT_TRUE((root.residual <= opts.tolerance * 100.0)) << "Every solve_numeric polynomial candidate is residual-verified";
        }
    }

    {
        LMCAS::SolveOptions opts;
        opts.allow_numeric = true;
        opts.tolerance = 1e-10;
        opts.max_newton_iterations = 100;
        const BigInt largest_finite_integer =
            (BigInt(1) << 1024) - (BigInt(1) << 971);
        LMCAS::Polynomial<Rational> endpoint_linear(
            {Rational(-largest_finite_integer), Rational(1)}, "x");
        auto endpoint_roots = LMCAS::solve_numeric_checked(
            LMCAS::poly_to_symbolic(endpoint_linear), "x", opts);
        EXPECT_TRUE((endpoint_roots && endpoint_roots.value().size() == 1 &&
                     endpoint_roots.value()[0].value ==
                         std::numeric_limits<double>::max()))
            << "numeric polynomial solve preserves a finite endpoint root";
    }
}

TEST(NewtonRaphsonDeflation, NewtonRaphsonDeflationCorrectlyContinuesToRemainingRoots) {
    LMCAS::Polynomial<Rational> poly("x");
    poly.coeffs = {Rational(2), Rational(-3), Rational(1)};

    auto intervals = LMCAS::isolate_real_roots_checked(poly).value();
    EXPECT_TRUE((intervals.size() == 2)) << "Sturm should isolate 2 roots for (x-1)(x-2)";

    auto expr = LMCAS::poly_to_symbolic(poly);
    auto df_expr = expr->differentiate("x");

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;

    int roots_found = 0;
    for (const auto &[lo_rat, hi_rat] : intervals) {
        double lo = lo_rat.to_double();
        double hi = hi_rat.to_double();
        double x0 = (lo + hi) * 0.5;

        auto result = LMCAS::newton_raphson_checked(expr, df_expr, "x", x0, lo, hi, opts).value();
        ASSERT_TRUE((result.has_value())) << "Newton should find root in interval [" + std::to_string(lo) + ", " + std::to_string(hi) + "]";
        if (result.has_value()) {
            roots_found++;

            double r = result->value;
            double residual = std::abs(r * r - 3 * r + 2);
            EXPECT_TRUE((residual < 1e-6)) << "Root " + std::to_string(r) + " should satisfy x^2-3x+2=0";
        }
    }
    EXPECT_TRUE((roots_found == 2)) << "Should find a root in each isolated interval";
}

TEST(NewtonRaphsonDeflation, NewtonRaphsonDeflationWithCubicPolynomial) {
    LMCAS::Polynomial<Rational> poly("x");
    poly.coeffs = {Rational(-6), Rational(11), Rational(-6), Rational(1)};

    auto intervals = LMCAS::isolate_real_roots_checked(poly).value();
    EXPECT_TRUE((intervals.size() == 3)) << "Sturm should isolate 3 roots for (x-1)(x-2)(x-3)";

    auto expr = LMCAS::poly_to_symbolic(poly);
    auto df_expr = expr->differentiate("x");

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;

    int roots_found = 0;
    for (const auto &[lo_rat, hi_rat] : intervals) {
        double lo = lo_rat.to_double();
        double hi = hi_rat.to_double();
        double x0 = (lo + hi) * 0.5;

        auto result = LMCAS::newton_raphson_checked(expr, df_expr, "x", x0, lo, hi, opts).value();
        ASSERT_TRUE((result.has_value())) << "Newton should find root in interval [" + std::to_string(lo) + ", " + std::to_string(hi) + "]";
        if (result.has_value()) {
            roots_found++;

            double r = result->value;
            double residual = std::abs((r - 1.0) * (r - 2.0) * (r - 3.0));
            EXPECT_TRUE((residual < 1e-6)) << "Root " + std::to_string(r) + " should satisfy (x-1)(x-2)(x-3)=0";
        }
    }

    EXPECT_TRUE((roots_found == 3)) << "Should find all 3 roots via Newton on isolated intervals";
}

TEST(NewtonRaphsonDeflation, NewtonRaphsonNonPolynomialInputRequiresX0InitialGuess) {
    auto x = SymbolicExpr::variable("x");

    auto f = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::number(-0.5));

    {
        LMCAS::SolveOptions opts;
        opts.allow_numeric = true;
        opts.tolerance = 1e-10;
        opts.max_newton_iterations = 100;
        opts.has_initial_guess = true;
        opts.initial_guess = 0.5;

        auto roots = LMCAS::solve_numeric_checked(f, "x", opts).value();
        EXPECT_TRUE((roots.size() <= 1)) << "Non-polynomial solve_numeric should return at most 1 root";
    }

    {
        LMCAS::SolveOptions opts;
        opts.allow_numeric = true;
        opts.tolerance = 1e-10;
        opts.max_newton_iterations = 100;
        opts.has_initial_guess = true;
        opts.initial_guess = 2.5;

        auto roots = LMCAS::solve_numeric_checked(f, "x", opts).value();
        EXPECT_TRUE((roots.size() <= 1)) << "Non-polynomial with different x0 should still return at most 1 root";
    }

    {
        LMCAS::SolveOptions opts;
        opts.allow_numeric = true;
        opts.tolerance = 1e-10;
        opts.max_newton_iterations = 100;
        opts.has_initial_guess = false;

        auto roots = LMCAS::solve_numeric_checked(f, "x", opts).value();
        EXPECT_TRUE((roots.size() <= 1)) << "Non-polynomial without explicit x0 should return at most 1 root";
    }

    test_polynomial_numeric_root_endpoints(x);
}

TEST(NewtonRaphsonDeflation, NewtonRaphsonNonConvergenceWithLimitedIterationsReturnsEmpty) {
    auto x = SymbolicExpr::variable("x");
    auto f = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(3)),
            SymbolicExpr::multiply(SymbolicExpr::number(-2), x)),
        SymbolicExpr::number(2));
    auto df = SymbolicExpr::add(
        SymbolicExpr::multiply(
            SymbolicExpr::number(3),
            SymbolicExpr::power(x, SymbolicExpr::number(2))),
        SymbolicExpr::number(-2));

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 1;

    auto result = LMCAS::newton_raphson_checked(f, df, "x", 5.0, opts).value();
    EXPECT_TRUE((!result.has_value())) << "Should not converge in 1 iteration from x=5 for x^3-2x+2";
}
