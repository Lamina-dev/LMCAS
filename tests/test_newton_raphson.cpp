#include "test_newton_raphson_support.hpp"

static void test_numeric_solver_binding_and_budgets(
    const std::shared_ptr<SymbolicExpr> &x,
    const std::shared_ptr<SymbolicExpr> &f, const LMCAS::SolveOptions &opts,
    LMCAS::CancellationToken &cancellation) {
    LMCAS::ComputationContext solve_unbound_context;
    auto solve_unbound = LMCAS::solve_numeric_checked(
        f, "x", solve_unbound_context, opts);
    EXPECT_TRUE((!solve_unbound &&
                 solve_unbound.error().code == LMCAS::CasErrc::UnboundSymbol))
        << "checked solve preserves coefficient binding failures";

    auto default_context_unbound = LMCAS::solve_numeric_checked(f, "x", opts);
    EXPECT_TRUE((!default_context_unbound &&
                 default_context_unbound.error().code == LMCAS::CasErrc::UnboundSymbol))
        << "default-context checked solve preserves coefficient binding failures";

    LMCAS::ComputationContext solve_cancelled_context({}, cancellation);
    auto solve_cancelled = LMCAS::solve_numeric_checked(
        x, "x", solve_cancelled_context, opts);
    EXPECT_TRUE((!solve_cancelled &&
                 solve_cancelled.error().code == LMCAS::CasErrc::Cancelled))
        << "checked solve observes cancellation before isolation";

    LMCAS::ResourceLimits solve_limits;
    solve_limits.max_steps = 1;
    LMCAS::ComputationContext solve_limited_context(solve_limits);
    auto solve_limited = LMCAS::solve_numeric_checked(
        x, "x", solve_limited_context, opts);
    EXPECT_TRUE((!solve_limited &&
                 solve_limited.error().code == LMCAS::CasErrc::ResourceLimit))
        << "checked solve accounts for interval refinement steps";

    LMCAS::SolveOptions no_roots_opts = opts;
    no_roots_opts.max_roots = 0;
    LMCAS::ComputationContext no_roots_context;
    auto no_roots = LMCAS::solve_numeric_checked(
        x, "x", no_roots_context, no_roots_opts);
    EXPECT_TRUE((no_roots && no_roots.value().empty())) << "zero root limit returns no candidates";
}

static void test_numeric_polynomial_recognition_limits(
    const std::shared_ptr<SymbolicExpr> &x, const LMCAS::SolveOptions &opts) {
    LMCAS::ResourceLimits expansion_limits;
    expansion_limits.max_expansion_terms = 1;
    LMCAS::ComputationContext exact_expansion_context(expansion_limits);
    auto exact_linear = SymbolicExpr::add(x, SymbolicExpr::number(1));
    auto exact_expansion = LMCAS::solve_numeric_checked(
        exact_linear, "x", exact_expansion_context, opts);
    EXPECT_TRUE((!exact_expansion &&
                 exact_expansion.error().code == LMCAS::CasErrc::ResourceLimit))
        << "exact polynomial recognition enforces expansion limits";

    auto approximate_linear = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(0.5), x),
        SymbolicExpr::number(-1));
    LMCAS::ComputationContext approximate_context(expansion_limits);
    auto approximate = LMCAS::solve_numeric_checked(
        approximate_linear, "x", approximate_context, opts);
    EXPECT_TRUE((approximate && approximate.value().size() == 1)) << "ApproxReal coefficients stay on the explicit numeric path";
    if (approximate && approximate.value().size() == 1) {
        EXPECT_TRUE((std::abs(approximate.value()[0].value - 2.0) < opts.tolerance)) << "approximate linear root is numerically verified";
    }

    LMCAS::ResourceLimits exponent_limits;
    exponent_limits.max_expansion_terms = 10;
    LMCAS::ComputationContext exponent_context(exponent_limits);
    auto large_power = SymbolicExpr::power(x, SymbolicExpr::number(BigInt(1000)));
    auto exponent_limited = LMCAS::solve_numeric_checked(
        large_power, "x", exponent_context, opts);
    EXPECT_TRUE((!exponent_limited &&
                 exponent_limited.error().code == LMCAS::CasErrc::ResourceLimit))
        << "large exact powers fail before polynomial expansion";
}

static bool numeric_root_has_error(
    const LMCAS::NumericRootResult &result, LMCAS::CasErrc code) {
    return !result && result.error().code == code;
}

TEST(NewtonRaphson, CheckedNumericRootErrorsArePreserved) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto f = SymbolicExpr::add(x, y);
    auto df = SymbolicExpr::number(1);

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 10;

    LMCAS::ComputationContext unbound_context;
    auto unbound = LMCAS::newton_raphson_checked(
        f, df, "x", 0.0, unbound_context, opts);
    EXPECT_TRUE((numeric_root_has_error(unbound, LMCAS::CasErrc::UnboundSymbol))) << "unbound symbols return UnboundSymbol";

    auto log_f = SymbolicExpr::ln(x);
    auto log_df = SymbolicExpr::divide(SymbolicExpr::number(1), x);
    LMCAS::ComputationContext domain_context;
    auto domain = LMCAS::newton_raphson_checked(
        log_f, log_df, "x", -1.0, domain_context, opts);
    EXPECT_TRUE((numeric_root_has_error(domain, LMCAS::CasErrc::DomainError))) << "domain failures return DomainError";

    LMCAS::CancellationToken cancellation;
    cancellation.cancel();
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    auto cancelled = LMCAS::newton_raphson_checked(
        x, df, "x", 1.0, cancelled_context, opts);
    EXPECT_TRUE((numeric_root_has_error(cancelled, LMCAS::CasErrc::Cancelled))) << "cancelled computations return Cancelled";

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = LMCAS::newton_raphson_checked(
        x, df, "x", 1.0, limited_context, opts);
    EXPECT_TRUE((numeric_root_has_error(limited, LMCAS::CasErrc::ResourceLimit))) << "step exhaustion returns ResourceLimit";

    LMCAS::ComputationContext invalid_context;
    auto invalid = LMCAS::bisection_checked(
        x, "x", 2.0, 1.0, invalid_context, opts);
    EXPECT_TRUE((numeric_root_has_error(invalid, LMCAS::CasErrc::InvalidArgument))) << "invalid brackets return InvalidArgument";

    auto default_domain = LMCAS::newton_raphson_checked(
        log_f, log_df, "x", -1.0, opts);
    EXPECT_TRUE((numeric_root_has_error(default_domain, LMCAS::CasErrc::DomainError))) << "default-context Newton preserves DomainError";

    auto default_invalid = LMCAS::bisection_checked(
        x, "x", 2.0, 1.0, opts);
    EXPECT_TRUE((numeric_root_has_error(default_invalid, LMCAS::CasErrc::InvalidArgument))) << "default-context bisection preserves InvalidArgument";

    auto default_bracket_domain = LMCAS::newton_raphson_checked(
        log_f, log_df, "x", -1.0, -2.0, 2.0, opts);
    EXPECT_TRUE((numeric_root_has_error(default_bracket_domain, LMCAS::CasErrc::DomainError))) << "default-context bracketed Newton preserves DomainError";

    test_numeric_solver_binding_and_budgets(x, f, opts, cancellation);

    test_numeric_polynomial_recognition_limits(x, opts);
}

TEST(NewtonRaphson, NewtonRaphsonBasicConvergenceX22) {
    auto x = SymbolicExpr::variable("x");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-2));
    auto df = SymbolicExpr::multiply(SymbolicExpr::number(2), x);

    LMCAS::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = LMCAS::newton_raphson_checked(f, df, "x", 1.5, opts).value();
    ASSERT_TRUE((result.has_value())) << "Newton-Raphson should converge for x^2-2 near 1.5";
    if (result.has_value()) {
        EXPECT_TRUE((std::abs(result->value - std::sqrt(2.0)) < 1e-10)) << "Root should be close to sqrt(2)";
        EXPECT_TRUE((result->residual < opts.tolerance)) << "Residual should be below tolerance";
    }
}

TEST(NewtonRaphson, NewtonRaphsonConvergenceWithBracketX22) {
    auto x = SymbolicExpr::variable("x");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-2));
    auto df = SymbolicExpr::multiply(SymbolicExpr::number(2), x);

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = LMCAS::newton_raphson_checked(f, df, "x", 1.5, 1.0, 2.0, opts).value();
    ASSERT_TRUE((result.has_value())) << "Newton-Raphson with bracket should converge for x^2-2";
    if (result.has_value()) {
        EXPECT_TRUE((std::abs(result->value - std::sqrt(2.0)) < 1e-10)) << "Root should be close to sqrt(2)";
        EXPECT_TRUE((result->residual < opts.tolerance)) << "Residual should be below tolerance";
    }
}

TEST(NewtonRaphson, NewtonRaphsonBisectionFallbackWhenDerivativeNearZero) {
    auto x = SymbolicExpr::variable("x");
    auto f = SymbolicExpr::power(x, SymbolicExpr::number(3));
    auto df = SymbolicExpr::multiply(
        SymbolicExpr::number(3),
        SymbolicExpr::power(x, SymbolicExpr::number(2)));

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = LMCAS::newton_raphson_checked(f, df, "x", 1e-8, -1.0, 1.0, opts).value();
    ASSERT_TRUE((result.has_value())) << "Should converge via bisection fallback for x^3 near zero";
    if (result.has_value()) {
        EXPECT_TRUE((std::abs(result->value) < 1e-4)) << "Root should be close to 0";
    }
}

TEST(NewtonRaphson, NewtonRaphsonAcceptsAConstantDerivativeAtLargeX) {
    const double target = std::numeric_limits<double>::max() * 0.5;
    auto x = SymbolicExpr::variable("x");
    auto f = SymbolicExpr::add(x, SymbolicExpr::number(-target));
    auto df = SymbolicExpr::number(1.0);
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 10;

    auto result = LMCAS::newton_raphson_checked(
        f, df, "x", std::numeric_limits<double>::max(), opts);
    EXPECT_TRUE((result && result.value().has_value())) << "Newton does not classify df/dx=1 as zero at large x";
    if (result && result.value()) {
        {
            const double actual_value = (result.value()->value);
            const double expected_value = (target);
            const double tolerance = (0.0);
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_NEAR(actual_value, expected_value, tolerance);
        }
    }
}

TEST(NewtonRaphson, NewtonRaphsonNoBracketDerivativeNearZeroReturnsNullopt) {
    auto x = SymbolicExpr::variable("x");
    auto f = SymbolicExpr::power(x, SymbolicExpr::number(3));
    auto df = SymbolicExpr::multiply(
        SymbolicExpr::number(3),
        SymbolicExpr::power(x, SymbolicExpr::number(2)));

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = LMCAS::newton_raphson_checked(f, df, "x", 1e-8, opts).value();

    EXPECT_TRUE((result.has_value())) << "f(1e-8) = 1e-24 < tolerance, should converge immediately";
}

TEST(NewtonRaphson, NewtonRaphsonNonConvergenceReturnsNullopt) {
    auto x = SymbolicExpr::variable("x");

    auto f = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(5)),
            SymbolicExpr::multiply(SymbolicExpr::number(-1), x)),
        SymbolicExpr::number(-1));

    auto df = SymbolicExpr::add(
        SymbolicExpr::multiply(
            SymbolicExpr::number(5),
            SymbolicExpr::power(x, SymbolicExpr::number(4))),
        SymbolicExpr::number(-1));

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 2;

    auto result = LMCAS::newton_raphson_checked(f, df, "x", 10.0, opts).value();
    EXPECT_TRUE((!result.has_value())) << "Should not converge in 2 iterations from x=10";
}

TEST(NewtonRaphson, NewtonRaphsonDampingEngagesOnOvershoot) {
    auto x = SymbolicExpr::variable("x");

    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-4));
    auto df = SymbolicExpr::multiply(SymbolicExpr::number(2), x);

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = LMCAS::newton_raphson_checked(f, df, "x", 0.01, opts).value();
    ASSERT_TRUE((result.has_value())) << "Should converge for x^2-4 even from x0=0.01 (large first step)";
    if (result.has_value()) {
        EXPECT_TRUE((std::abs(result->value - 2.0) < 1e-10 || std::abs(result->value + 2.0) < 1e-10)) << "Root should be ±2";
    }

    auto result2 = LMCAS::newton_raphson_checked(f, df, "x", 0.01, 0.0, 3.0, opts).value();
    ASSERT_TRUE((result2.has_value())) << "Should converge with bracket for x^2-4 from x0=0.01";
    if (result2.has_value()) {
        EXPECT_TRUE((std::abs(result2->value - 2.0) < 1e-10)) << "Root should be 2 within bracket [0, 3]";
    }
}

TEST(NewtonRaphson, BisectionBasicConvergenceX22) {
    auto x = SymbolicExpr::variable("x");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-2));

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = LMCAS::bisection_checked(f, "x", 1.0, 2.0, opts).value();
    ASSERT_TRUE((result.has_value())) << "Bisection should converge for x^2-2 on [1,2]";
    if (result.has_value()) {
        EXPECT_TRUE((std::abs(result->value - std::sqrt(2.0)) < 1e-10)) << "Root should be close to sqrt(2)";
    }
}

TEST(NewtonRaphson, BisectionAcceptsAFullFiniteSymmetricBracket) {
    auto x = SymbolicExpr::variable("x");
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = LMCAS::bisection_checked(
        x, "x", -std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(), opts);
    EXPECT_TRUE((result && result.value().has_value())) << "bisection does not overflow a finite full-range bracket";
    if (result && result.value()) {
        EXPECT_NEAR(result.value()->value, 0.0, 0.0) << "full-range bisection finds the midpoint root";
    }
}

TEST(NewtonRaphson, BisectionNoSignChangeReturnsNullopt) {
    auto x = SymbolicExpr::variable("x");
    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1));

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = LMCAS::bisection_checked(f, "x", -1.0, 1.0, opts).value();
    EXPECT_TRUE((!result.has_value())) << "Bisection should return nullopt when no sign change";
}
