#include "test_mixed_transcendental_support.hpp"

TEST(MixedTranscendentalDeduplication, DeduplicateRootsAscendingOrder) {
    std::vector<NumericRoot> roots = {
        {3.0, 1e-14, 5},
        {1.0, 1e-14, 3},
        {5.0, 1e-14, 7},
        {2.0, 1e-14, 4},
        {4.0, 1e-14, 6}};

    lmmc_real_t tolerance = 1e-12;
    int max_roots = -1;

    auto sorted = deduplicate_roots(roots, tolerance, max_roots);

    ASSERT_EQ(sorted.size(), 5u) << "All 5 distinct roots should survive deduplication";

    for (size_t i = 0; i + 1 < sorted.size(); ++i) {
        EXPECT_TRUE((sorted[i] < sorted[i + 1])) << "sorted[" + std::to_string(i) + "]=" + std::to_string(sorted[i]) +
                                                        " should be < sorted[" + std::to_string(i + 1) + "]=" +
                                                        std::to_string(sorted[i + 1]);
    }

    for (size_t i = 0; i < sorted.size(); ++i) {
        const double expected = static_cast<double>(i + 1);
        SCOPED_TRACE(i);
        EXPECT_TRUE(std::isfinite(sorted[i]));
        EXPECT_TRUE(std::isfinite(expected));
        EXPECT_NEAR(sorted[i], expected, 1e-10);
    }
}

TEST(MixedTranscendentalDeduplication, PipelineSinXAscending) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    auto derivative = compute_derivative(expr, "x");

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    LMCAS::SearchInterval interval{-10.0, 10.0};

    auto intervals = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    std::vector<NumericRoot> refined_roots;
    for (const auto &candidate : intervals) {
        auto root =
            LMCAS::refine_root(expr, derivative, "x", candidate, opts);
        ASSERT_TRUE(root.has_value());
        refined_roots.push_back(*root);
    }

    auto sorted =
        deduplicate_roots(refined_roots, opts.tolerance, opts.max_roots);
    ASSERT_EQ(sorted.size(), 7U);
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        EXPECT_NEAR(
            sorted[i], (static_cast<double>(i) - 3.0) * M_PI, 1e-10);
    }
}

TEST(MixedTranscendentalDeduplication, PipelineCosXAscending) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::cos(x);
    auto derivative = compute_derivative(expr, "x");

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    opts.max_roots = -1;

    LMCAS::SearchInterval interval{-10.0, 10.0};

    auto intervals = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    std::vector<NumericRoot> refined_roots;
    for (const auto &candidate : intervals) {
        auto root =
            LMCAS::refine_root(expr, derivative, "x", candidate, opts);
        ASSERT_TRUE(root.has_value());
        refined_roots.push_back(*root);
    }

    auto sorted =
        deduplicate_roots(refined_roots, opts.tolerance, opts.max_roots);
    ASSERT_EQ(sorted.size(), 6U);
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        EXPECT_NEAR(
            sorted[i],
            (static_cast<double>(i) - 2.5) * M_PI,
            1e-10);
    }
}

TEST(MixedTranscendentalDeduplication, DeduplicateWithDuplicatesAscending) {
    lmmc_real_t tolerance = 1e-10;
    std::vector<NumericRoot> roots = {
        {5.0, 1e-13, 4},
        {1.0, 1e-13, 2},
        {1.0 + 5e-10, 1e-14, 3},
        {3.0, 1e-13, 5},
        {3.0 + 2e-10, 1e-12, 6},
        {-2.0, 1e-13, 1}};

    int max_roots = -1;
    auto sorted = deduplicate_roots(roots, tolerance, max_roots);

    EXPECT_TRUE((sorted.size() == 4)) << "Should have 4 distinct roots after deduplication, got " +
                                             std::to_string(sorted.size());

    for (size_t i = 0; i + 1 < sorted.size(); ++i) {
        EXPECT_TRUE((sorted[i] < sorted[i + 1])) << "After dedup: sorted[" + std::to_string(i) + "]=" +
                                                        std::to_string(sorted[i]) + " should be < sorted[" +
                                                        std::to_string(i + 1) + "]=" + std::to_string(sorted[i + 1]);
    }
}

TEST(MixedTranscendentalDeduplication, MaxRootsLimitAfterDeduplication) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    auto derivative = compute_derivative(expr, "x");

    LMCAS::SolveOptions opts;
    opts.tolerance = 1e-12;
    opts.max_newton_iterations = 100;
    opts.max_roots = 2;

    LMCAS::SearchInterval interval{-10.0, 10.0};
    auto intervals = LMCAS::isolate_roots(expr, derivative, "x", interval, opts);

    std::vector<NumericRoot> refined_roots;
    for (const auto &candidate : intervals) {
        auto root_result =
            LMCAS::refine_root(expr, derivative, "x", candidate, opts);
        ASSERT_TRUE(root_result.has_value());
        refined_roots.push_back(*root_result);
    }

    auto deduped = LMCAS::deduplicate_roots(
        refined_roots, opts.tolerance, opts.max_roots);
    ASSERT_EQ(deduped.size(), 2U);
    EXPECT_LT(deduped[0], deduped[1]);
    EXPECT_TRUE(std::isfinite(deduped[0]));
    EXPECT_TRUE(std::isfinite(deduped[1]));
}

TEST(MixedTranscendentalDeduplication, DeduplicateRootsDirect) {
    lmmc_real_t tolerance = 1e-12;
    lmmc_real_t dedup_threshold = 10.0 * tolerance;

    std::vector<NumericRoot> roots;
    roots.push_back(NumericRoot{1.0, 1e-13, 5});
    roots.push_back(NumericRoot{1.0 + 5e-12, 2e-13, 7});
    roots.push_back(NumericRoot{2.0, 3e-13, 4});
    roots.push_back(NumericRoot{2.0 + 3e-12, 1e-14, 6});
    roots.push_back(NumericRoot{3.0, 5e-14, 3});

    auto result = LMCAS::deduplicate_roots(roots, tolerance, -1);

    EXPECT_TRUE((result.size() == 3)) << "Expected 3 roots after deduplication, got " + std::to_string(result.size());

    for (size_t i = 0; i < result.size(); ++i) {
        for (size_t j = i + 1; j < result.size(); ++j) {
            double diff = std::abs(result[i] - result[j]);
            EXPECT_TRUE((diff >= dedup_threshold)) << "Pair (" + std::to_string(i) + "," + std::to_string(j) +
                                                          "): |" + std::to_string(result[i]) + " - " + std::to_string(result[j]) +
                                                          "| = " + std::to_string(diff) + " should be >= " + std::to_string(dedup_threshold);
        }
    }
}

TEST(MixedTranscendentalDeduplication, DeduplicateRootsNoDuplicates) {
    lmmc_real_t tolerance = 1e-12;
    lmmc_real_t dedup_threshold = 10.0 * tolerance;

    std::vector<NumericRoot> roots;
    roots.push_back(NumericRoot{-5.0, 1e-13, 3});
    roots.push_back(NumericRoot{0.0, 2e-14, 4});
    roots.push_back(NumericRoot{5.0, 1e-13, 5});

    auto result = LMCAS::deduplicate_roots(roots, tolerance, -1);

    EXPECT_TRUE((result.size() == 3)) << "All well-separated roots should be kept, got " + std::to_string(result.size());

    for (size_t i = 0; i < result.size(); ++i) {
        for (size_t j = i + 1; j < result.size(); ++j) {
            double diff = std::abs(result[i] - result[j]);
            EXPECT_TRUE((diff >= dedup_threshold)) << "Pair (" + std::to_string(i) + "," + std::to_string(j) +
                                                          "): |" + std::to_string(result[i]) + " - " + std::to_string(result[j]) +
                                                          "| = " + std::to_string(diff) + " should be >= " + std::to_string(dedup_threshold);
        }
    }
}

TEST(MixedTranscendentalDeduplication, DeduplicateRootsAllDuplicates) {
    lmmc_real_t tolerance = 1e-12;

    std::vector<NumericRoot> roots;
    roots.push_back(NumericRoot{1.0, 5e-13, 10});
    roots.push_back(NumericRoot{1.0 + 2e-12, 1e-13, 8});
    roots.push_back(NumericRoot{1.0 + 4e-12, 3e-13, 9});
    roots.push_back(NumericRoot{1.0 + 7e-12, 2e-13, 7});

    auto result = LMCAS::deduplicate_roots(roots, tolerance, -1);

    EXPECT_TRUE((result.size() == 1)) << "All roots within threshold should collapse to 1, got " + std::to_string(result.size());
}

TEST(MixedTranscendentalDeduplication, DeduplicateUnsortedValues) {
    std::vector<LMCAS::NumericRoot> roots = {
        {3.14159, 1e-13, 5},
        {-1.5, 2e-14, 3},
        {0.0, 0.0, 1},
        {2.71828, 5e-14, 7}};

    lmmc_real_t tolerance = 1e-12;
    int max_roots = -1;

    auto deduped = LMCAS::deduplicate_roots(roots, tolerance, max_roots);

    ASSERT_EQ(deduped.size(), 4U);
    EXPECT_NEAR(deduped[0], -1.5, 1e-15);
    EXPECT_NEAR(deduped[1], 0.0, 1e-15);
    EXPECT_NEAR(deduped[2], 2.71828, 1e-10);
    EXPECT_NEAR(deduped[3], 3.14159, 1e-10);
}

TEST(MixedTranscendentalDeduplication, DeduplicateEmptyValues) {
    std::vector<LMCAS::NumericRoot> roots;
    lmmc_real_t tolerance = 1e-12;
    int max_roots = -1;

    auto deduped = LMCAS::deduplicate_roots(roots, tolerance, max_roots);
    EXPECT_TRUE((deduped.empty())) << "Deduplication of empty input should return empty";
}

TEST(MixedTranscendentalDeduplication, DeduplicateClusterValues) {
    std::vector<LMCAS::NumericRoot> roots = {
        {1.0, 1e-13, 5},
        {1.0 + 5e-12, 2e-13, 6},
        {2.0, 3e-14, 3},
        {2.0 + 1e-13, 1e-14, 4},
        {-3.0, 0.0, 1}};

    lmmc_real_t tolerance = 1e-12;
    int max_roots = -1;

    auto deduped = LMCAS::deduplicate_roots(roots, tolerance, max_roots);

    ASSERT_EQ(deduped.size(), 3U);
    EXPECT_NEAR(deduped[0], -3.0, 1e-10);
    EXPECT_NEAR(deduped[1], 1.0, 1e-10);
    EXPECT_NEAR(deduped[2], 2.0, 1e-10);
}
