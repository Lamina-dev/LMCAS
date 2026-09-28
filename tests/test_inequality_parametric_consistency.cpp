#include "test_inequality_support.hpp"

TEST(InequalityParametricConsistency, PositiveQuadraticHasEntireLine) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1));

    auto result = InequalitySolver::solve_parametric_inequality_checked(
        expr, InequalityType::GreaterThan, "x", {"a"});

    ASSERT_TRUE(result) << result.error().message;
    ASSERT_EQ(result.value().cases.size(), 1U);
    EXPECT_EQ(result.value().cases.front().condition, nullptr);
    EXPECT_TRUE(result.value().cases.front().solution.is_entire_line());
}

TEST(InequalityParametricConsistency, DistinctRootsGiveOpenExteriors) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-1));

    auto result = InequalitySolver::solve_parametric_inequality_checked(
        expr, InequalityType::GreaterThan, "x", {"a"});

    ASSERT_TRUE(result) << result.error().message;
    ASSERT_EQ(result.value().cases.size(), 1U);
    EXPECT_EQ(result.value().cases.front().condition, nullptr);
    const auto& intervals = result.value().cases.front().solution.intervals();
    ASSERT_EQ(intervals.size(), 2U);
    EXPECT_TRUE(intervals[0].lower.is_neg_infinity);
    EXPECT_TRUE(intervals[0].upper.is_open);
    ASSERT_TRUE(intervals[0].upper.value);
    EXPECT_TRUE(test_proved_equivalent(intervals[0].upper.value,
                                       SymbolicExpr::number(-1)));
    EXPECT_TRUE(intervals[1].lower.is_open);
    ASSERT_TRUE(intervals[1].lower.value);
    EXPECT_TRUE(test_proved_equivalent(intervals[1].lower.value,
                                       SymbolicExpr::number(1)));
    EXPECT_TRUE(intervals[1].upper.is_pos_infinity);
    EXPECT_FALSE(result.value().cases.front().solution.contains(-1));
    EXPECT_FALSE(result.value().cases.front().solution.contains(1));
    EXPECT_TRUE(result.value().cases.front().solution.contains(-2));
    EXPECT_FALSE(result.value().cases.front().solution.contains(0));
    EXPECT_TRUE(result.value().cases.front().solution.contains(2));
}

TEST(InequalityParametricConsistency, RepeatedRootIsExcluded) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::power(x, SymbolicExpr::number(2));

    auto result = InequalitySolver::solve_parametric_inequality_checked(
        expr, InequalityType::GreaterThan, "x", {});

    ASSERT_TRUE(result) << result.error().message;
    ASSERT_EQ(result.value().cases.size(), 1U);
    const auto& intervals = result.value().cases.front().solution.intervals();
    ASSERT_EQ(intervals.size(), 2U);
    EXPECT_TRUE(intervals[0].lower.is_neg_infinity);
    EXPECT_TRUE(intervals[0].upper.is_open);
    ASSERT_TRUE(intervals[0].upper.value);
    EXPECT_TRUE(test_proved_equivalent(intervals[0].upper.value,
                                       SymbolicExpr::number(0)));
    EXPECT_TRUE(intervals[1].lower.is_open);
    ASSERT_TRUE(intervals[1].lower.value);
    EXPECT_TRUE(test_proved_equivalent(intervals[1].lower.value,
                                       SymbolicExpr::number(0)));
    EXPECT_TRUE(intervals[1].upper.is_pos_infinity);
    EXPECT_TRUE(result.value().cases.front().solution.contains(-1));
    EXPECT_FALSE(result.value().cases.front().solution.contains(0));
    EXPECT_TRUE(result.value().cases.front().solution.contains(1));
}
