#include "test_mixed_transcendental_support.hpp"

TEST(MixedTranscendentalRefinement, RefineRootSinXNewtonRaphson) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    auto derivative = compute_derivative(expr, "x"); /**< 导数为 cos(x)。 */

    LMCAS::IsolatedInterval interval{3.0, 3.5, true};
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = LMCAS::refine_root(expr, derivative, "x", interval, opts);

    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(std::isfinite(result->value));
    EXPECT_NEAR(result->value, M_PI, 1e-10);
    EXPECT_LT(result->residual, opts.tolerance);
}

TEST(MixedTranscendentalRefinement, RefineRootSinXBisectionFallback) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);

    LMCAS::IsolatedInterval interval{3.0, 3.5, false};
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = LMCAS::refine_root(expr, nullptr, "x", interval, opts);

    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(std::isfinite(result->value));
    EXPECT_NEAR(result->value, M_PI, 1e-10);
    EXPECT_LT(result->residual, opts.tolerance);
}

TEST(MixedTranscendentalRefinement, RefineRootExpXMinusXMinus2) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
            SymbolicExpr::number(-2)));
    auto derivative = compute_derivative(expr, "x"); /**< 导数为 exp(x) - 1。 */

    LMCAS::IsolatedInterval interval{1.0, 1.5, true};
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = LMCAS::refine_root(expr, derivative, "x", interval, opts);

    ASSERT_TRUE(result.has_value());
    EXPECT_LT(result->residual, opts.tolerance);
    EXPECT_GE(result->value, interval.lo);
    EXPECT_LE(result->value, interval.hi);
}

TEST(MixedTranscendentalRefinement, RefineRootXCosXMinus1) {
    auto x = SymbolicExpr::variable("x");
    /**
     * @brief x*cos(x) - 1 在 x ≈ 4.917 处有根。
     * f(4.5) ≈ -1.949 < 0，f(5.0) ≈ 0.418 > 0，区间端点异号。
     */
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(x, SymbolicExpr::cos(x)),
        SymbolicExpr::number(-1));
    auto derivative = compute_derivative(expr, "x");

    LMCAS::IsolatedInterval interval{4.5, 5.0, true};
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;

    auto result = LMCAS::refine_root(expr, derivative, "x", interval, opts);

    ASSERT_TRUE(result.has_value());
    EXPECT_LT(result->residual, opts.tolerance);
    EXPECT_GE(result->value, interval.lo);
    EXPECT_LE(result->value, interval.hi);
}

TEST(MixedTranscendentalRefinement, RefineRootDiscardsWhenNoConvergence) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    auto derivative = compute_derivative(expr, "x");

    LMCAS::IsolatedInterval interval{-1.0, 1.0, false};
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 50;

    auto result = LMCAS::refine_root(expr, derivative, "x", interval, opts);

    EXPECT_FALSE((result.has_value())) << "refine_root should return nullopt for x^2+1 (no real root in [-1, 1])";
}

TEST(MixedTranscendentalRefinement, RefineRootDerivativeZeroBisectionStep) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::power(x, SymbolicExpr::number(3));
    auto derivative = compute_derivative(expr, "x");

    LMCAS::IsolatedInterval interval{-1.0, 1.0, false};
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 200;

    auto result = LMCAS::refine_root(expr, derivative, "x", interval, opts);

    ASSERT_TRUE(result.has_value());
    EXPECT_NEAR(result->value, 0.0, 1e-4);
    EXPECT_LT(result->residual, opts.tolerance);
}

TEST(MixedTranscendentalRefinement, RootsWithinIntervalSinX) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    auto derivative = compute_derivative(expr, "x"); /**< 导数为 cos(x)。 */

    LMCAS::SearchInterval interval{-3.0, 3.0};
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;
    opts.has_search_interval = true;
    opts.search_lo = -3.0;
    opts.search_hi = 3.0;

    auto isolated =
        LMCAS::isolate_roots(expr, derivative, "x", interval, opts);
    ASSERT_EQ(isolated.size(), 1U);
    auto central =
        LMCAS::refine_root(expr, derivative, "x", isolated[0], opts);
    ASSERT_TRUE(central.has_value());
    EXPECT_NEAR(central->value, 0.0, 1e-10);
    EXPECT_LT(central->residual, opts.tolerance);

    LMCAS::SearchInterval interval2{-4.0, 4.0};
    opts.search_lo = -4.0;
    opts.search_hi = 4.0;
    auto isolated2 =
        LMCAS::isolate_roots(expr, derivative, "x", interval2, opts);
    ASSERT_EQ(isolated2.size(), 3U);
    std::vector<double> roots;
    for (const auto &candidate : isolated2) {
        auto refined =
            LMCAS::refine_root(expr, derivative, "x", candidate, opts);
        ASSERT_TRUE(refined.has_value());
        EXPECT_LT(refined->residual, opts.tolerance);
        roots.push_back(refined->value);
    }
    std::sort(roots.begin(), roots.end());
    EXPECT_NEAR(roots[0], -M_PI, 1e-10);
    EXPECT_NEAR(roots[1], 0.0, 1e-10);
    EXPECT_NEAR(roots[2], M_PI, 1e-10);
}

TEST(MixedTranscendentalRefinement, RootsWithinIntervalCosXMinusHalf) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::cos(x),
        SymbolicExpr::number(-0.5));
    auto derivative = SymbolicExpr::multiply(
        SymbolicExpr::number(-1),
        SymbolicExpr::sin(x));

    LMCAS::SearchInterval interval{0.0, 5.0};
    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;
    opts.has_search_interval = true;
    opts.search_lo = 0.0;
    opts.search_hi = 5.0;

    auto isolated =
        LMCAS::isolate_roots(expr, derivative, "x", interval, opts);
    ASSERT_EQ(isolated.size(), 1U);
    auto refined =
        LMCAS::refine_root(expr, derivative, "x", isolated[0], opts);
    ASSERT_TRUE(refined.has_value());
    EXPECT_NEAR(refined->value, M_PI / 3.0, 1e-10);
    EXPECT_LT(refined->residual, opts.tolerance);
}

TEST(MixedTranscendentalRefinement, OutputValiditySinXRefinedRoots) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    auto derivative = compute_derivative(expr, "x"); /**< 导数为 cos(x)。 */

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    LMCAS::SearchInterval interval{-10.0, 10.0};
    auto intervals =
        LMCAS::isolate_roots(expr, derivative, "x", interval, opts);
    ASSERT_EQ(intervals.size(), 7U);

    std::vector<double> roots;
    for (const auto &candidate : intervals) {
        auto root =
            LMCAS::refine_root(expr, derivative, "x", candidate, opts);
        ASSERT_TRUE(root.has_value());
        ASSERT_TRUE(std::isfinite(root->value));
        EXPECT_LT(root->residual, opts.tolerance);
        roots.push_back(root->value);
    }
    std::sort(roots.begin(), roots.end());
    for (std::size_t i = 0; i < roots.size(); ++i) {
        EXPECT_NEAR(
            roots[i], (static_cast<double>(i) - 3.0) * M_PI, 1e-10);
    }
}

TEST(MixedTranscendentalRefinement, OutputValidityExpXMinusXMinus2) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
            SymbolicExpr::number(-2)));
    auto derivative = compute_derivative(expr, "x");

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    LMCAS::SearchInterval interval{-5.0, 5.0};
    auto intervals =
        LMCAS::isolate_roots(expr, derivative, "x", interval, opts);
    ASSERT_EQ(intervals.size(), 2U);

    for (const auto &candidate : intervals) {
        auto root =
            LMCAS::refine_root(expr, derivative, "x", candidate, opts);
        ASSERT_TRUE(root.has_value());
        ASSERT_TRUE(std::isfinite(root->value));
        auto substituted =
            expr->substitute("x", SymbolicExpr::number(root->value));
        ASSERT_NE(substituted, nullptr);
        EXPECT_LT(std::abs(substituted->to_numeric()), opts.tolerance);
    }
}

TEST(MixedTranscendentalRefinement, OutputValidityXCosXMinus1) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(x, SymbolicExpr::cos(x)),
        SymbolicExpr::number(-1));
    auto derivative = compute_derivative(expr, "x");

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    LMCAS::SearchInterval interval{-10.0, 10.0};
    auto intervals =
        LMCAS::isolate_roots(expr, derivative, "x", interval, opts);
    ASSERT_FALSE(intervals.empty());

    std::vector<NumericRoot> all_roots;
    for (const auto &candidate : intervals) {
        auto root =
            LMCAS::refine_root(expr, derivative, "x", candidate, opts);
        ASSERT_TRUE(root.has_value());
        ASSERT_TRUE(std::isfinite(root->value));
        EXPECT_LT(root->residual, opts.tolerance);
        all_roots.push_back(*root);
    }

    auto deduped =
        LMCAS::deduplicate_roots(all_roots, opts.tolerance, -1);
    ASSERT_FALSE(deduped.empty());
    for (double root : deduped) {
        EXPECT_TRUE(std::isfinite(root));
    }
}

TEST(MixedTranscendentalRefinement, HugeSameSignIntervalUsesFiniteMidpoint) {
    auto x = SymbolicExpr::variable("x");
    constexpr double expected = 1.3e308;
    auto expression = SymbolicExpr::add(
        x, SymbolicExpr::number(-expected));
    auto derivative = SymbolicExpr::number(1);
    IsolatedInterval interval{1.0e308, 1.6e308, true};
    SolveOptions options;
    options.tolerance = 1.0e292;
    options.max_newton_iterations = 20;

    auto result = refine_root(
        expression, derivative, "x", interval, options);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(std::isfinite(result->value));
    EXPECT_GE(result->value, interval.lo);
    EXPECT_LE(result->value, interval.hi);
    EXPECT_NEAR(result->value / expected, 1.0, 1e-15);
    EXPECT_LE(result->residual, options.tolerance);
}
