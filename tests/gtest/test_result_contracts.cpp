#include <gtest/gtest.h>

#include "solve_strategies.hpp"
#include "symbolic.hpp"

namespace lamina {
namespace {

TEST(ResultContractsTest, SolutionSetKindsAreExplicit) {
    EXPECT_EQ(SolutionSet::empty().kind(), SolutionSet::Kind::Empty);

    auto x = SymbolicExpr::variable("x");
    auto finite = SolutionSet::finite({FiniteSolution{x, 2, {}}});
    ASSERT_EQ(finite.kind(), SolutionSet::Kind::Finite);
    ASSERT_EQ(finite.finite_solutions().size(), 1U);
    EXPECT_EQ(finite.finite_solutions().front().multiplicity, 2);

    auto inconclusive = SolutionSet::inconclusive("unsupported equation class");
    EXPECT_EQ(inconclusive.kind(), SolutionSet::Kind::Inconclusive);
    EXPECT_EQ(inconclusive.reason(), "unsupported equation class");
}

TEST(ResultContractsTest, SolveSeparatesMathematicalOutcomesFromErrors) {
    auto x = SymbolicExpr::variable("x");
    ComputationContext context;

    auto empty = solve_equation(
        SymbolicExpr::number(1), "x", context, SolveOptions{});
    ASSERT_TRUE(empty);
    EXPECT_EQ(empty.value().kind(), SolutionSet::Kind::Empty);

    ComputationContext universal_context;
    auto universal = solve_equation(
        SymbolicExpr::number(0), "x", universal_context, SolveOptions{});
    ASSERT_TRUE(universal);
    EXPECT_EQ(universal.value().kind(), SolutionSet::Kind::Universal);

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled({}, cancellation);
    auto error = solve_equation(x, "x", cancelled, SolveOptions{});
    ASSERT_FALSE(error);
    EXPECT_EQ(error.error().code, CasErrc::Cancelled);
}

TEST(ResultContractsTest, ResourceLimitsAreErrors) {
    auto x = SymbolicExpr::variable("x");
    auto equation = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1));
    ResourceLimits limits;
    limits.max_steps = 1;
    ComputationContext context(limits);

    auto result = solve_equation(equation, "x", context, SolveOptions{});
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::ResourceLimit);
}

} // namespace
} // namespace lamina
