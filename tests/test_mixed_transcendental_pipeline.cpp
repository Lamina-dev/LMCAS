#include "test_mixed_transcendental_support.hpp"

void expect_factor_roots_sin_x_minus_half_times_x_minus_3(
    const std::vector<std::shared_ptr<SymbolicExpr>> &results) {
    double expected_pi_6 = M_PI / 6.0;        /**< 约为 0.5236。 */
    double expected_5pi_6 = 5.0 * M_PI / 6.0; /**< 约为 2.6180。 */

    bool found_pi_6 = false;
    bool found_5pi_6 = false;
    bool found_3 = false;

    for (const auto &root : results) {
        auto value = test_numeric_eval(root);
        ASSERT_TRUE(value.has_value()) << "every returned root is numerically evaluable";
        found_pi_6 = found_pi_6 || std::abs(*value - expected_pi_6) < 1e-4;
        found_5pi_6 = found_5pi_6 || std::abs(*value - expected_5pi_6) < 1e-4;
        found_3 = found_3 || std::abs(*value - 3.0) < 1e-4;
    }

    EXPECT_TRUE((found_pi_6)) << "Should find root from factor sin(x)-0.5 at π/6 ≈ " + std::to_string(expected_pi_6);
    EXPECT_TRUE((found_5pi_6)) << "Should find root from factor sin(x)-0.5 at 5π/6 ≈ " + std::to_string(expected_5pi_6);
    EXPECT_TRUE((found_3)) << "Should find root from factor x-3 at x = 3";
}

TEST(MixedTranscendentalPipeline, FactoredSinXMinusHalfTimesXMinus3) {
    auto x = SymbolicExpr::variable("x");
    /**
     * @brief (sin(x) - 0.5)*(x - 3) 在 [0, 4] 内的结果须包含两个因子的根。
     * sin(x) - 0.5 的根为 pi/6 ≈ 0.524 和 5pi/6 ≈ 2.618；x - 3 的根为 3。
     */
    auto factor1 = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::number(-0.5));
    auto factor2 = SymbolicExpr::add(x, SymbolicExpr::number(-3));
    auto expr = SymbolicExpr::multiply(factor1, factor2);

    LMCAS::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;
    opts.has_search_interval = true;
    opts.search_lo = 0.0;
    opts.search_hi = 4.0;

    auto checked = solve_mixed_transcendental_checked(expr, "x", opts);
    ASSERT_TRUE(checked) << checked.error().message;
    const auto &results = checked.value().value;
    ASSERT_GE(results.size(), 3U);

    expect_factor_roots_sin_x_minus_half_times_x_minus_3(results);
}

void expect_factor_roots_sin_x_times_cos_x(
    const std::vector<std::shared_ptr<SymbolicExpr>> &results) {
    bool found_sin_zero = false;
    bool found_sin_pi = false;

    bool found_cos_pi_2 = false;
    bool found_cos_3pi_2 = false;

    for (const auto &root : results) {
        auto value = test_numeric_eval(root);
        ASSERT_TRUE(value.has_value()) << "every returned root is numerically evaluable";
        found_sin_zero = found_sin_zero || std::abs(*value) < 1e-4;
        found_sin_pi = found_sin_pi || std::abs(*value - M_PI) < 1e-4;
        found_cos_pi_2 = found_cos_pi_2 || std::abs(*value - M_PI / 2.0) < 1e-4;
        found_cos_3pi_2 = found_cos_3pi_2 ||
                         std::abs(*value - 3.0 * M_PI / 2.0) < 1e-4;
    }

    EXPECT_TRUE((found_sin_zero)) << "Should find root from sin(x) factor at x = 0";
    EXPECT_TRUE((found_sin_pi)) << "Should find root from sin(x) factor at x = π ≈ " + std::to_string(M_PI);
    EXPECT_TRUE((found_cos_pi_2)) << "Should find root from cos(x) factor at x = π/2 ≈ " + std::to_string(M_PI / 2.0);
    EXPECT_TRUE((found_cos_3pi_2)) << "Should find root from cos(x) factor at x = 3π/2 ≈ " + std::to_string(3.0 * M_PI / 2.0);
}

TEST(MixedTranscendentalPipeline, FactoredSinXTimesCosX) {
    auto x = SymbolicExpr::variable("x");
    /**
     * @brief sin(x)*cos(x) 在 [-10, 10] 内的结果须包含两个因子的根。
     * sin(x) 的根为 0、±pi、±2pi、±3pi；cos(x) 的根为 ±pi/2、±3pi/2、±5pi/2。
     */
    auto expr = SymbolicExpr::multiply(
        SymbolicExpr::sin(x),
        SymbolicExpr::cos(x));

    LMCAS::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;
    opts.has_search_interval = true;
    opts.search_lo = -10.0;
    opts.search_hi = 10.0;

    auto checked = solve_mixed_transcendental_checked(expr, "x", opts);
    ASSERT_TRUE(checked) << checked.error().message;
    const auto &results = checked.value().value;
    ASSERT_GE(results.size(), 6U);

    expect_factor_roots_sin_x_times_cos_x(results);
}

TEST(MixedTranscendentalPipeline, E2eSinXPlusXRootAtZero) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::multiply(SymbolicExpr::number(-0.5), x));

    LMCAS::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;
    opts.has_search_interval = true;
    opts.search_lo = -3.0;
    opts.search_hi = 3.0;

    auto checked = solve_mixed_transcendental_checked(expr, "x", opts);
    ASSERT_TRUE(checked) << checked.error().message;
    const auto &results = checked.value().value;
    ASSERT_FALSE(results.empty());

    for (const auto &root : results) {
        auto value = test_numeric_eval(root);
        ASSERT_TRUE(value.has_value()) << "every returned root is numerically evaluable";
        const double residual = std::sin(*value) - 0.5 * *value;
        EXPECT_LT(std::abs(residual), opts.tolerance);
    }
}

TEST(MixedTranscendentalPipeline, E2eExpXMinusXMinus2) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), x),
            SymbolicExpr::number(-2)));

    LMCAS::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    auto checked = solve_mixed_transcendental_checked(expr, "x", opts);
    ASSERT_TRUE(checked) << checked.error().message;
    const auto &results = checked.value().value;
    ASSERT_EQ(results.size(), 2U);

    bool found_neg = false;
    bool found_pos = false;
    for (const auto &root : results) {
        auto value = test_numeric_eval(root);
        ASSERT_TRUE(value.has_value()) << "every returned root is numerically evaluable";
        found_neg = found_neg || std::abs(*value + 1.8414) < 0.01;
        found_pos = found_pos || std::abs(*value - 1.1462) < 0.01;
    }
    EXPECT_TRUE((found_neg)) << "exp(x)-x-2=0 should have root near x ≈ -1.84";
    EXPECT_TRUE((found_pos)) << "exp(x)-x-2=0 should have root near x ≈ 1.15";

    for (const auto &root : results) {
        auto value = test_numeric_eval(root);
        ASSERT_TRUE(value.has_value()) << "every returned root is numerically evaluable";
        const double residual = std::exp(*value) - *value - 2.0;
        EXPECT_LT(std::abs(residual), opts.tolerance);
    }
}

TEST(MixedTranscendentalPipeline, E2eXCosXMinus1) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(x, SymbolicExpr::cos(x)),
        SymbolicExpr::number(-1));

    LMCAS::SolveOptions opts;
    opts.allow_numeric = true;
    opts.tolerance = 1e-10;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    auto checked = solve_mixed_transcendental_checked(expr, "x", opts);
    ASSERT_TRUE(checked) << checked.error().message;
    const auto &results = checked.value().value;
    ASSERT_GE(results.size(), 2U);

    for (const auto &root : results) {
        auto value = test_numeric_eval(root);
        ASSERT_TRUE(value.has_value()) << "every returned root is numerically evaluable";
        const double residual = *value * std::cos(*value) - 1.0;
        EXPECT_LT(std::abs(residual), opts.tolerance);
    }

    for (size_t i = 0; i + 1 < results.size(); ++i) {
        auto left = test_numeric_eval(results[i]);
        auto right = test_numeric_eval(results[i + 1]);
        ASSERT_TRUE(left.has_value());
        ASSERT_TRUE(right.has_value());
        EXPECT_LT(*left, *right);
    }
}
