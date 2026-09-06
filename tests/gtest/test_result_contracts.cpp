#include <gtest/gtest.h>

#include "solve_strategies.hpp"
#include "symbolic.hpp"

using namespace LMCAS;

namespace LMCAS {
namespace {

TEST(ResultContractsTest, SolutionSetAlternativesAreExplicit) {
    SolutionSet empty = EmptySolutions{};
    EXPECT_TRUE(std::holds_alternative<EmptySolutions>(empty));

    auto x = SymbolicExpr::variable("x");
    SolutionSet finite = FiniteSolutions{{FiniteSolution{x, 2, {}}}};
    const auto* values = std::get_if<FiniteSolutions>(&finite);
    ASSERT_NE(values, nullptr);
    ASSERT_EQ(values->values.size(), 1U);
    EXPECT_EQ(values->values.front().multiplicity, 2U);
}

TEST(ResultContractsTest, SolveSeparatesMathematicalOutcomesFromErrors) {
    auto x = SymbolicExpr::variable("x");
    ComputationContext context;

    auto empty = solve_equation(
        SymbolicExpr::number(1), "x", context, SolveOptions{});
    ASSERT_TRUE(empty);
    EXPECT_TRUE(std::holds_alternative<EmptySolutions>(empty.value()));

    ComputationContext universal_context;
    auto universal = solve_equation(
        SymbolicExpr::number(0), "x", universal_context, SolveOptions{});
    ASSERT_TRUE(universal);
    EXPECT_TRUE(std::holds_alternative<UniversalSolutions>(universal.value()));

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
} // namespace LMCAS
