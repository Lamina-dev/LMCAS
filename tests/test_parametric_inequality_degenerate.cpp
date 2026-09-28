#include "test_common.hpp"
#include "assumption_context.hpp"
#include "inequality_solver.hpp"
#include "interval.hpp"
#include "symbolic.hpp"

using namespace LMCAS;

TEST(ParametricInequalityDegenerate, SignedLeadingCoefficientIncludesZeroBranch) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto expr = SymbolicExpr::multiply(
        a, SymbolicExpr::power(x, SymbolicExpr::number(2)));

    auto result = InequalitySolver::solve_parametric_inequality_checked(
        expr, InequalityType::GreaterThan, "x", {"a"});

    ASSERT_TRUE(result) << result.error().message;
    const auto& cases = result.value().cases;
    ASSERT_EQ(cases.size(), 3U);
    AssumptionContext facts;
    for (int sign : {-1, 0, 1}) {
        unsigned active = 0;
        for (const auto& branch : cases) {
            ASSERT_TRUE(branch.condition);
            auto condition = branch.condition->substitute(
                "a", SymbolicExpr::number(sign));
            auto truth = facts.evaluate_condition(*condition);
            ASSERT_NE(truth, Tribool::Unknown) << condition->to_string();
            if (truth != Tribool::True) continue;
            ++active;
            if (sign > 0) {
                const auto& intervals = branch.solution.intervals();
                ASSERT_EQ(intervals.size(), 2U);
                EXPECT_TRUE(intervals[0].lower.is_neg_infinity);
                EXPECT_TRUE(intervals[0].upper.is_open);
                ASSERT_TRUE(intervals[0].upper.value);
                EXPECT_TRUE(test_proved_equivalent(
                    intervals[0].upper.value, SymbolicExpr::number(0)));
                EXPECT_TRUE(intervals[1].lower.is_open);
                ASSERT_TRUE(intervals[1].lower.value);
                EXPECT_TRUE(test_proved_equivalent(
                    intervals[1].lower.value, SymbolicExpr::number(0)));
                EXPECT_TRUE(intervals[1].upper.is_pos_infinity);
                EXPECT_TRUE(branch.solution.contains(-1));
                EXPECT_FALSE(branch.solution.contains(0));
                EXPECT_TRUE(branch.solution.contains(1));
            } else {
                EXPECT_TRUE(branch.solution.is_empty());
            }
        }
        EXPECT_EQ(active, 1U) << "a=" << sign;
    }
}

TEST(ParametricInequalityDegenerate, InvalidAndNonPolynomialInputsKeepErrorCodes) {
    auto x = SymbolicExpr::variable("x");
    auto null_result = InequalitySolver::solve_parametric_inequality_checked(
        nullptr, InequalityType::GreaterThan, "x", {"a"});
    ASSERT_FALSE(null_result);
    EXPECT_EQ(null_result.error().code, CasErrc::InvalidArgument);

    auto invalid_variable = InequalitySolver::solve_parametric_inequality_checked(
        x, InequalityType::GreaterThan, "", {"a"});
    ASSERT_FALSE(invalid_variable);
    EXPECT_EQ(invalid_variable.error().code, CasErrc::InvalidArgument);

    auto nonpolynomial = InequalitySolver::solve_parametric_inequality_checked(
        SymbolicExpr::sin(x), InequalityType::GreaterThan, "x", {"a"});
    ASSERT_FALSE(nonpolynomial);
    EXPECT_EQ(nonpolynomial.error().code, CasErrc::UnsupportedExpression);
}

TEST(ParametricInequalityDegenerate, ContextErrorsAreNotReclassified) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto expr = SymbolicExpr::add(x, a);

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled({}, cancellation);
    auto stopped = InequalitySolver::solve_parametric_inequality_checked(
        expr, InequalityType::GreaterThan, "x", {"a"}, cancelled);
    ASSERT_FALSE(stopped);
    EXPECT_EQ(stopped.error().code, CasErrc::Cancelled);

    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext exhausted(limits);
    auto limited = InequalitySolver::solve_parametric_inequality_checked(
        expr, InequalityType::GreaterThan, "x", {"a"}, exhausted);
    ASSERT_FALSE(limited);
    EXPECT_EQ(limited.error().code, CasErrc::ResourceLimit);
}
