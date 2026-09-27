#include "test_mixed_transcendental_support.hpp"

TEST(MixedTranscendentalSearchInterval, SearchIntervalUserSpecified) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = true;
    opts.search_lo = 2.0;
    opts.search_hi = 5.0;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "User-specified valid interval should return a value";
    if (result) {
        EXPECT_NEAR(result->lo, 2.0, 1e-15) << "lo should be 2.0";
        EXPECT_NEAR(result->hi, 5.0, 1e-15) << "hi should be 5.0";
    }
}

TEST(MixedTranscendentalSearchInterval, SearchIntervalUserInvalidLoGeHi) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = true;
    opts.search_lo = 5.0;
    opts.search_hi = 2.0;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    EXPECT_FALSE((result.has_value())) << "Invalid interval (lo >= hi) should return nullopt";
}

TEST(MixedTranscendentalSearchInterval, SearchIntervalUserInvalidEqual) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = true;
    opts.search_lo = 3.0;
    opts.search_hi = 3.0;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    EXPECT_FALSE((result.has_value())) << "Equal bounds (lo == hi) should return nullopt";
}

TEST(MixedTranscendentalSearchInterval, SearchIntervalUserWidthLeTolerance) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = true;
    opts.search_lo = 1.0;
    opts.search_hi = 1.0 + 1e-13; /**< 区间宽度 1e-13 小于容差 1e-12。 */
    opts.tolerance = 1e-12;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    EXPECT_FALSE((result.has_value())) << "Interval width <= tolerance should return nullopt";
}

TEST(MixedTranscendentalSearchInterval, SearchIntervalDefaultNoPeriodic) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::exp(x),
                                  SymbolicExpr::multiply(SymbolicExpr::number(-1), x));

    LMCAS::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "Default interval should be valid";
    if (result) {
        EXPECT_NEAR(result->lo, -10.0, 1e-15) << "Default lo should be -10";
        EXPECT_NEAR(result->hi, 10.0, 1e-15) << "Default hi should be 10";
    }
}

TEST(MixedTranscendentalSearchInterval, SearchIntervalSinXPeriodicExtension) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "Periodic extension should produce valid interval";
    if (result) {
        EXPECT_NEAR(result->lo, -10.0, 1e-10) << "sin(x): period 2π < 10, so default [-10,10] used";
        EXPECT_NEAR(result->hi, 10.0, 1e-10) << "sin(x): period 2π < 10, so default [-10,10] used";
    }
}

TEST(MixedTranscendentalSearchInterval, SearchIntervalSinSmallKPeriodicExtension) {
    auto x = SymbolicExpr::variable("x");
    auto sin_arg = SymbolicExpr::multiply(SymbolicExpr::number(0.2), x);
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(sin_arg), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "Periodic extension should produce valid interval";
    if (result) {
        double expected_period = 2.0 * M_PI / 0.2; /**< 半区间宽度为周期 10pi ≈ 31.4。 */
        {
            const double actual_value = (result->lo);
            const double expected_value = (-expected_period);
            const double tolerance = (1e-10);
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_NEAR(actual_value, expected_value, tolerance);
        }
        {
            const double actual_value = (result->hi);
            const double expected_value = (expected_period);
            const double tolerance = (1e-10);
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_NEAR(actual_value, expected_value, tolerance);
        }
    }
}

TEST(MixedTranscendentalSearchInterval, SearchIntervalTanPeriodicExtension) {
    auto x = SymbolicExpr::variable("x");
    auto tan_arg = SymbolicExpr::multiply(SymbolicExpr::number(1.0 / 3.0), x);
    auto expr = SymbolicExpr::add(SymbolicExpr::tan(tan_arg), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "Periodic extension should produce valid interval";
    if (result) {
        EXPECT_NEAR(result->lo, -10.0, 1e-10) << "tan(x/3): period 3π < 10, default used";
        EXPECT_NEAR(result->hi, 10.0, 1e-10) << "tan(x/3): period 3π < 10, default used";
    }
}

TEST(MixedTranscendentalSearchInterval, SearchIntervalTanSmallKExtension) {
    auto x = SymbolicExpr::variable("x");
    auto tan_arg = SymbolicExpr::multiply(SymbolicExpr::number(0.1), x);
    auto expr = SymbolicExpr::add(SymbolicExpr::tan(tan_arg), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "Periodic extension should produce valid interval";
    if (result) {
        double expected_period = M_PI / 0.1; /**< 周期为 10pi ≈ 31.4。 */
        {
            const double actual_value = (result->lo);
            const double expected_value = (-expected_period);
            const double tolerance = (1e-10);
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_NEAR(actual_value, expected_value, tolerance);
        }
        {
            const double actual_value = (result->hi);
            const double expected_value = (expected_period);
            const double tolerance = (1e-10);
            EXPECT_TRUE(std::isfinite(actual_value));
            EXPECT_NEAR(actual_value, expected_value, tolerance);
        }
    }
}

TEST(MixedTranscendentalSearchInterval, SearchIntervalClampTo100) {
    auto x = SymbolicExpr::variable("x");
    auto sin_arg = SymbolicExpr::multiply(SymbolicExpr::number(0.01), x);
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(sin_arg), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "Clamped interval should be valid";
    if (result) {
        EXPECT_NEAR(result->lo, -100.0, 1e-15) << "Should be clamped to -100";
        EXPECT_NEAR(result->hi, 100.0, 1e-15) << "Should be clamped to +100";
    }
}

TEST(MixedTranscendentalSearchInterval, SearchIntervalNonlinearArgDefault) {
    auto x = SymbolicExpr::variable("x");
    auto sin_arg = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(sin_arg), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "Non-linear periodic arg should fall back to default";
    if (result) {
        EXPECT_NEAR(result->lo, -10.0, 1e-15) << "Non-linear arg: default lo = -10";
        EXPECT_NEAR(result->hi, 10.0, 1e-15) << "Non-linear arg: default hi = 10";
    }
}

TEST(MixedTranscendentalSearchInterval, SearchIntervalCos2xPlus1) {
    auto x = SymbolicExpr::variable("x");
    auto cos_arg = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::number(1));
    auto expr = SymbolicExpr::add(SymbolicExpr::cos(cos_arg), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "Linear periodic arg should produce valid interval";
    if (result) {
        EXPECT_NEAR(result->lo, -10.0, 1e-10) << "cos(2x+1): period π < 10, default used";
        EXPECT_NEAR(result->hi, 10.0, 1e-10) << "cos(2x+1): period π < 10, default used";
    }
}

TEST(MixedTranscendentalSearchInterval, SearchIntervalUserOverridesPeriodic) {
    auto x = SymbolicExpr::variable("x");
    auto sin_arg = SymbolicExpr::multiply(SymbolicExpr::number(0.01), x);
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(sin_arg), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = true;
    opts.search_lo = -1.0;
    opts.search_hi = 1.0;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "User-specified interval should override periodic extension";
    if (result) {
        EXPECT_NEAR(result->lo, -1.0, 1e-15) << "User override: lo = -1";
        EXPECT_NEAR(result->hi, 1.0, 1e-15) << "User override: hi = 1";
    }
}

TEST(MixedTranscendentalSearchInterval, SinXPlusXTwoPeriods) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "sin(x)+x should produce valid interval";
    if (result) {
        double width = result->hi - result->lo;
        double min_required = 2.0 * (2.0 * M_PI / 1.0); /**< 两个周期宽度为 4pi ≈ 12.57。 */
        EXPECT_TRUE((width >= min_required - 1e-10)) << "sin(x)+x: interval width " + std::to_string(width) +
                                                            " should be >= 4π ≈ " + std::to_string(min_required);
    }
}

TEST(MixedTranscendentalSearchInterval, CosHalfXMinusXTwoPeriods) {
    auto x = SymbolicExpr::variable("x");
    auto cos_arg = SymbolicExpr::multiply(SymbolicExpr::number(0.5), x);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::cos(cos_arg),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x));

    LMCAS::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "cos(0.5*x)-x should produce valid interval";
    if (result) {
        double width = result->hi - result->lo;
        double min_required = 2.0 * (2.0 * M_PI / 0.5); /**< 两个周期宽度为 8pi ≈ 25.13。 */
        EXPECT_TRUE((width >= min_required - 1e-10)) << "cos(0.5*x)-x: interval width " + std::to_string(width) +
                                                            " should be >= 8π ≈ " + std::to_string(min_required);
    }
}

TEST(MixedTranscendentalSearchInterval, TanXPlusXTwoPeriods) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::tan(x), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "tan(x)+x should produce valid interval";
    if (result) {
        double width = result->hi - result->lo;
        double min_required = 2.0 * (M_PI / 1.0); /**< 两个周期宽度为 2pi ≈ 6.28。 */
        EXPECT_TRUE((width >= min_required - 1e-10)) << "tan(x)+x: interval width " + std::to_string(width) +
                                                            " should be >= 2π ≈ " + std::to_string(min_required);
    }
}

TEST(MixedTranscendentalSearchInterval, SinSmallKClamped) {
    auto x = SymbolicExpr::variable("x");
    /**
     * @brief sin(0.1*x) 的周期为 20pi ≈ 62.83，两个周期宽度为 40pi ≈ 125.66。
     * 半区间宽度为一个周期且小于 100，使用 [-20pi, 20pi]。
     */
    auto sin_arg = SymbolicExpr::multiply(SymbolicExpr::number(0.1), x);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(sin_arg),
        SymbolicExpr::multiply(SymbolicExpr::number(0.01), x));

    LMCAS::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "sin(0.1*x)+x/100 should produce valid interval";
    if (result) {
        double width = result->hi - result->lo;
        double min_required = 2.0 * (2.0 * M_PI / 0.1); /**< 两个周期宽度约为 125.66。 */
        EXPECT_TRUE((width >= min_required - 1e-10)) << "sin(0.1*x)+x/100: interval width " + std::to_string(width) +
                                                            " should be >= 2*(2π/0.1) ≈ " + std::to_string(min_required);
        EXPECT_TRUE((result->lo >= -100.0 - 1e-10)) << "sin(0.1*x): lo should be >= -100 (clamping)";
        EXPECT_TRUE((result->hi <= 100.0 + 1e-10)) << "sin(0.1*x): hi should be <= 100 (clamping)";
    }
}

TEST(MixedTranscendentalSearchInterval, SinXSquaredNonlinearDefault) {
    auto x = SymbolicExpr::variable("x");
    auto sin_arg = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(sin_arg), x);

    LMCAS::SolveOptions opts;
    opts.has_search_interval = false;

    auto result = LMCAS::determine_search_interval(expr, "x", opts);
    ASSERT_TRUE(result.has_value()) << "sin(x^2)+x should produce valid interval";
    if (result) {
        EXPECT_NEAR(result->lo, -10.0, 1e-15) << "sin(x^2)+x: non-linear arg, default lo=-10";
        EXPECT_NEAR(result->hi, 10.0, 1e-15) << "sin(x^2)+x: non-linear arg, default hi=10";
    }
}
