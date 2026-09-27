#include "symbolic.hpp"
#include "test_common.hpp"

using namespace LMCAS;

TEST(LmcasSolver, NegativeLinearRoot) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::number(6));
    auto result = solve_finite_checked(expression, "x");
    ASSERT_TRUE(result.has_value());
    ASSERT_EQ(result.value().size(), 1u);
    EXPECT_TRUE(test_expression_text(result.value()[0]->simplify(), "-3"));
}

TEST(LmcasSolver, PositiveLinearRoot) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(3), x),
        SymbolicExpr::number(-9));
    auto result = solve_finite_checked(expression, "x");
    ASSERT_TRUE(result.has_value());
    ASSERT_EQ(result.value().size(), 1u);
    EXPECT_TRUE(test_expression_text(result.value()[0]->simplify(), "3"));
}

TEST(LmcasSolver, MonicLinearRoot) {
    auto expression = SymbolicExpr::add(
        SymbolicExpr::variable("x"), SymbolicExpr::number(-5));
    auto result = solve_finite_checked(expression, "x");
    ASSERT_TRUE(result.has_value());
    ASSERT_EQ(result.value().size(), 1u);
    EXPECT_TRUE(test_expression_text(result.value()[0]->simplify(), "5"));
}

TEST(LmcasSolver, TwoQuadraticRoots) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-4));
    auto result = solve_finite_checked(expression, "x");
    ASSERT_TRUE(result.has_value());
    ASSERT_EQ(result.value().size(), 2u);
    const bool first_order =
        test_expression_text(result.value()[0]->simplify(), "-2") &&
        test_expression_text(result.value()[1]->simplify(), "2");
    const bool second_order =
        test_expression_text(result.value()[0]->simplify(), "2") &&
        test_expression_text(result.value()[1]->simplify(), "-2");
    EXPECT_TRUE(first_order || second_order);
}

TEST(LmcasSolver, RepeatedQuadraticRoot) {
    auto x = SymbolicExpr::variable("x");
    auto n2 = SymbolicExpr::number(2);
    auto n1 = SymbolicExpr::number(1);
    auto x2 = SymbolicExpr::power(x, n2);
    auto term2 = SymbolicExpr::multiply(n2, x);
    auto part1 = SymbolicExpr::add(x2, term2);
    auto expr5 = SymbolicExpr::add(part1, n1);
    auto result = solve_finite_checked(expr5, "x");
    ASSERT_TRUE(result.has_value());
    ASSERT_EQ(result.value().size(), 2u)
        << "x^2 + 2x + 1 preserves the double root";
    EXPECT_TRUE(test_expression_text(result.value()[0]->simplify(), "-1"));
    EXPECT_TRUE(test_expression_text(result.value()[1]->simplify(), "-1"));
}
