#include "test_common.hpp"
#include "inequality_solver.hpp"
#include "interval.hpp"
#include "symbolic.hpp"

using namespace LMCAS;

TEST(ParametricInequality, UnitSlopeKeepsSymbolicRoot) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto expr = SymbolicExpr::add(x, a);

    auto result = InequalitySolver::solve_parametric_inequality_checked(
        expr, InequalityType::GreaterThan, "x", {"a"});

    ASSERT_TRUE(result) << result.error().message;
    ASSERT_EQ(result.value().cases.size(), 1U);
    EXPECT_EQ(result.value().cases.front().condition, nullptr);
    const auto& intervals = result.value().cases.front().solution.intervals();
    ASSERT_EQ(intervals.size(), 1U);
    EXPECT_TRUE(intervals.front().lower.is_open);
    EXPECT_FALSE(intervals.front().lower.is_neg_infinity);
    EXPECT_TRUE(intervals.front().upper.is_pos_infinity);
    ASSERT_TRUE(intervals.front().lower.value);
    EXPECT_TRUE(test_proved_equivalent(
        intervals.front().lower.value,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), a)));
    auto root = intervals.front().lower.value->substitute(
        "a", SymbolicExpr::number(3));
    EXPECT_TRUE(test_proved_equivalent(root, SymbolicExpr::number(-3)));
}

TEST(ParametricInequality, UndeterminedDiscriminantsAndDegreeAreInconclusive) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");
    auto c = SymbolicExpr::variable("c");
    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                std::vector<std::string>>> unsupported = {
        {SymbolicExpr::add(x2, SymbolicExpr::multiply(SymbolicExpr::number(-1), a)), {"a"}},
        {SymbolicExpr::add(SymbolicExpr::add(x2, SymbolicExpr::multiply(b, x)), c), {"b", "c"}},
        {SymbolicExpr::add(SymbolicExpr::multiply(a, x), b), {"a", "b"}},
        {SymbolicExpr::add(
             SymbolicExpr::add(SymbolicExpr::multiply(a, x2),
                               SymbolicExpr::multiply(SymbolicExpr::number(3), x)),
             SymbolicExpr::number(-2)), {"a"}},
        {SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(3)), a), {"a"}},
    };

    for (const auto& [expr, parameters] : unsupported) {
        SCOPED_TRACE(expr->to_string());
        auto result = InequalitySolver::solve_parametric_inequality_checked(
            expr, InequalityType::GreaterThan, "x", parameters);
        ASSERT_FALSE(result);
        EXPECT_EQ(result.error().code, CasErrc::Inconclusive);
    }
}
