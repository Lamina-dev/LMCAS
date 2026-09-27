#include "test_mixed_transcendental_support.hpp"

TEST(MixedTranscendentalIsolation, MaxRoots2SinX) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);

    auto derivative = SymbolicExpr::cos(x);

    LMCAS::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-12;
    opts.max_roots = 2;

    LMCAS::SearchInterval interval{-10.0, 10.0};

    auto result = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE((result.size() <= 2)) << "sin(x) on [-10,10] with max_roots=2: got " + std::to_string(result.size()) +
                                             " intervals, expected at most 2";
}

TEST(MixedTranscendentalIsolation, MaxRoots1SinX) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);

    auto derivative = SymbolicExpr::cos(x);

    LMCAS::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-12;
    opts.max_roots = 1;

    LMCAS::SearchInterval interval{-10.0, 10.0};

    auto result = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE((result.size() <= 1)) << "sin(x) on [-10,10] with max_roots=1: got " + std::to_string(result.size()) +
                                             " intervals, expected at most 1";
}

TEST(MixedTranscendentalIsolation, MaxRootsUnlimitedSinX) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);

    auto derivative = SymbolicExpr::cos(x);

    LMCAS::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-12;
    opts.max_roots = -1; /**< 根数量不限。 */

    LMCAS::SearchInterval interval{-10.0, 10.0};

    auto result = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE((result.size() == 7)) << "the uncapped sine search retains all seven roots in this interval";
    for (int k = -3; k <= 3; ++k) {
        const double root = k * std::acos(-1.0);
        bool covered = false;
        for (const auto &bracket : result)
            covered = covered || (bracket.lo <= root && root <= bracket.hi);
        EXPECT_TRUE((covered)) << "each independently known sine root belongs to a returned bracket";
    }
}

TEST(MixedTranscendentalIsolation, MaxRoots3CosXMinusHalf) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::cos(x),
        SymbolicExpr::number(-0.5));

    auto derivative = SymbolicExpr::multiply(
        SymbolicExpr::number(-1),
        SymbolicExpr::sin(x));

    LMCAS::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-12;
    opts.max_roots = 3;

    LMCAS::SearchInterval interval{-10.0, 10.0};

    auto result = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE((result.size() <= 3)) << "cos(x)-0.5 on [-10,10] with max_roots=3: got " + std::to_string(result.size()) +
                                             " intervals, expected at most 3";
}

TEST(MixedTranscendentalIsolation, ExactEndpointRootIsRetained) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(SymbolicExpr::sin(x),
                                        SymbolicExpr::multiply(SymbolicExpr::number(0.1), x));
    const LMCAS::SearchInterval search{0.0, 0.001};
    LMCAS::SolveOptions options;
    options.tolerance = 1e-12;
    auto intervals = LMCAS::isolate_roots(expression, nullptr, "x", search, options);
    EXPECT_TRUE((intervals.size() == 1)) << "the strictly increasing fixture has its one endpoint root";
    if (intervals.size() == 1) {
        EXPECT_TRUE((intervals[0].lo == 0.0 && intervals[0].hi >= 0.0 && intervals[0].hi <= search.hi)) << "the returned bracket contains the exact zero inside the search interval";
        auto residual = eval_at(expression, "x", intervals[0].lo);
        EXPECT_TRUE((residual && *residual == 0.0)) << "the retained endpoint satisfies the original expression";
    }
}

static bool brackets_zero(double lower_value, double upper_value) {
    return lower_value == 0.0 || upper_value == 0.0 ||
           std::signbit(lower_value) != std::signbit(upper_value);
}

TEST(MixedTranscendentalIsolation, SignChangeSinXPlusXDiv10) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);

    auto derivative = compute_derivative(expr, "x");

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    LMCAS::SearchInterval interval{-10.0, 10.0};
    auto intervals = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE((intervals.size() == 7)) << "all seven simple sine roots are isolated in this fixture";

    for (size_t i = 0; i < intervals.size(); ++i) {
        auto f_lo = eval_at(expr, "x", intervals[i].lo);
        auto f_hi = eval_at(expr, "x", intervals[i].hi);

        EXPECT_TRUE((f_lo.has_value() && f_hi.has_value())) << "both bracket endpoints are finite evaluations";
        if (f_lo.has_value() && f_hi.has_value()) {
            const bool sign_change = brackets_zero(*f_lo, *f_hi);
            EXPECT_TRUE((sign_change)) << "Interval [" + std::to_string(intervals[i].lo) + ", " +
                                              std::to_string(intervals[i].hi) + "]: f(lo)=" +
                                              std::to_string(*f_lo) + ", f(hi)=" + std::to_string(*f_hi) +
                                              " must bracket a zero of the continuous function";
        }
    }
}

TEST(MixedTranscendentalIsolation, SignChangeExpXMinusXMinus2) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
            SymbolicExpr::number(-2)));

    auto derivative = compute_derivative(expr, "x");

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    LMCAS::SearchInterval interval{-5.0, 5.0};
    auto intervals = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE((intervals.size() == 2)) << "both simple exponential roots are isolated in this fixture";

    for (size_t i = 0; i < intervals.size(); ++i) {
        auto f_lo = eval_at(expr, "x", intervals[i].lo);
        auto f_hi = eval_at(expr, "x", intervals[i].hi);

        EXPECT_TRUE((f_lo.has_value() && f_hi.has_value())) << "both bracket endpoints are finite evaluations";
        if (f_lo.has_value() && f_hi.has_value()) {
            const bool sign_change = brackets_zero(*f_lo, *f_hi);
            EXPECT_TRUE((sign_change)) << "Interval [" + std::to_string(intervals[i].lo) + ", " +
                                              std::to_string(intervals[i].hi) + "]: f(lo)=" +
                                              std::to_string(*f_lo) + ", f(hi)=" + std::to_string(*f_hi) +
                                              " must bracket a zero of the continuous function";
        }
    }
}

TEST(MixedTranscendentalIsolation, SignChangeXCosXMinus1) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(x, SymbolicExpr::cos(x)),
        SymbolicExpr::number(-1));

    auto derivative = compute_derivative(expr, "x");

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    LMCAS::SearchInterval interval{-10.0, 10.0};
    auto intervals = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    const double pi = std::acos(-1.0);
    bool found_reference_root = false;

    for (size_t i = 0; i < intervals.size(); ++i) {
        auto f_lo = eval_at(expr, "x", intervals[i].lo);
        auto f_hi = eval_at(expr, "x", intervals[i].hi);
        found_reference_root = found_reference_root ||
                               (pi < intervals[i].lo && intervals[i].hi < 2.0 * pi);

        EXPECT_TRUE((f_lo.has_value() && f_hi.has_value())) << "both bracket endpoints are finite evaluations";
        if (f_lo.has_value() && f_hi.has_value()) {
            const bool sign_change = brackets_zero(*f_lo, *f_hi);
            EXPECT_TRUE((sign_change)) << "Interval [" + std::to_string(intervals[i].lo) + ", " +
                                              std::to_string(intervals[i].hi) + "]: f(lo)=" +
                                              std::to_string(*f_lo) + ", f(hi)=" + std::to_string(*f_hi) +
                                              " must bracket a zero of the continuous function";
        }
    }
    EXPECT_TRUE((found_reference_root)) << "the independently known sign change between pi and 2*pi is represented";
}

TEST(MixedTranscendentalIsolation, XSquaredPlusOneAlwaysPositive) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    auto derivative = expr->differentiate("x");

    LMCAS::SearchInterval interval{-10.0, 10.0};
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    auto isolated = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE((isolated.empty())) << "x^2 + 1 is always positive on [-10, 10], isolate_roots should return empty";
}

TEST(MixedTranscendentalIsolation, ExpXPlusOneAlwaysPositive) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::number(1));
    auto derivative = expr->differentiate("x");

    LMCAS::SearchInterval interval{-10.0, 10.0};
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    auto isolated = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE((isolated.empty())) << "exp(x) + 1 is always positive on [-10, 10], isolate_roots should return empty";
}

TEST(MixedTranscendentalIsolation, NegXSquaredMinusOneAlwaysNegative) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::multiply(
        SymbolicExpr::number(-1),
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::number(1)));
    auto derivative = expr->differentiate("x");

    LMCAS::SearchInterval interval{-10.0, 10.0};
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    auto isolated = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE((isolated.empty())) << "-(x^2 + 1) is always negative on [-10, 10], isolate_roots should return empty";
}

TEST(MixedTranscendentalIsolation, SinXPlusFiveAlwaysPositive) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::number(5));
    auto derivative = expr->differentiate("x");

    LMCAS::SearchInterval interval{-10.0, 10.0};
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_roots = -1;

    auto isolated = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    EXPECT_TRUE((isolated.empty())) << "sin(x) + 5 is always positive on [-10, 10], isolate_roots should return empty";
}
