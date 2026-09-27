#include "expr.hpp"
#include "test_common.hpp"

using namespace LMCAS;

TEST(LmcasComplexArithmetic, DifferenceOfSquares) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto two = SymbolicExpr::number(2);

    auto A = SymbolicExpr::power(SymbolicExpr::add(x, y), two);
    auto B = SymbolicExpr::power(SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), y)), two);
    auto expr = SymbolicExpr::add(A, SymbolicExpr::multiply(SymbolicExpr::number(-1), B));

    auto expanded = expr->expand();
    auto parsed_expected = parse_expr("4*x*y");
    ASSERT_TRUE(parsed_expected.has_value());
    auto expected = parsed_expected.value()->simplify();
    EXPECT_TRUE(test_same_expression(expanded, expected));
    EXPECT_TRUE(test_proved_equivalent(expanded, expected)) << "difference of squares is 4xy";
}

TEST(LmcasComplexArithmetic, RationalCoefficients) {
    auto x = SymbolicExpr::variable("x");

    auto half_x = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)), x);
    auto third_x = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 3)), x);
    auto sum = SymbolicExpr::add(half_x, third_x);

    auto result = sum->simplify();
    auto parsed_expected = parse_expr("(5/6)*x");
    ASSERT_TRUE(parsed_expected.has_value());
    auto expected = parsed_expected.value()->simplify();
    EXPECT_TRUE(test_same_expression(result, expected));
    EXPECT_TRUE(test_proved_equivalent(result, expected)) << "rational coefficients add exactly";
}

TEST(LmcasComplexArithmetic, CubicCancellation) {
    auto x = SymbolicExpr::variable("x");
    auto one = SymbolicExpr::number(1);

    auto term1 = SymbolicExpr::power(SymbolicExpr::add(x, one), SymbolicExpr::number(3));
    auto term2 = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(x, SymbolicExpr::number(3)));
    auto term3 = SymbolicExpr::number(-1);

    auto poly = SymbolicExpr::add(SymbolicExpr::add(term1, term2), term3);
    auto res = poly->expand();

    auto parsed_expected = parse_expr("3*x^2+3*x");
    ASSERT_TRUE(parsed_expected.has_value());
    auto expected = parsed_expected.value()->simplify();
    EXPECT_TRUE(test_same_expression(res, expected));
    EXPECT_TRUE(test_proved_equivalent(res, expected)) << "cubic leading and constant terms cancel";
}

TEST(LmcasComplexArithmetic, AdditiveInverses) {
    auto x = SymbolicExpr::variable("x");
    auto one = SymbolicExpr::number(1);

    auto p = SymbolicExpr::add(x, one);
    auto neg_p = SymbolicExpr::multiply(SymbolicExpr::number(-1), p);
    auto zero_expr = SymbolicExpr::add(p, neg_p);

    auto res = zero_expr->expand()->simplify();
    auto parsed_expected = parse_expr("0");
    ASSERT_TRUE(parsed_expected.has_value());
    auto expected = parsed_expected.value()->simplify();
    EXPECT_TRUE(test_same_expression(res, expected));
    EXPECT_TRUE(test_proved_equivalent(res, expected)) << "additive inverses cancel";
}

TEST(LmcasComplexArithmetic, NestedIntegerPowers) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);

    auto p2 = SymbolicExpr::power(x, two);
    auto p6 = SymbolicExpr::power(p2, SymbolicExpr::number(3));
    auto res = p6->simplify();
    auto parsed_expected = parse_expr("x^6");
    ASSERT_TRUE(parsed_expected.has_value());
    auto expected = parsed_expected.value()->simplify();
    EXPECT_TRUE(test_same_expression(res, expected));
    EXPECT_TRUE(test_proved_equivalent(res, expected)) << "nested integer powers multiply their exponents";
}

TEST(LmcasComplexArithmetic, AdditiveIdentity) {
    auto z = SymbolicExpr::number(0);
    auto term = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(z, term);
    auto sim = expr->simplify();

    auto parsed_expected = parse_expr("x");
    ASSERT_TRUE(parsed_expected.has_value());
    auto expected = parsed_expected.value()->simplify();
    EXPECT_TRUE(test_same_expression(sim, expected));
    EXPECT_TRUE(test_proved_equivalent(sim, expected)) << "zero is the additive identity";
}

TEST(LmcasComplexArithmetic, SymbolicSummands) {
    auto sum = SymbolicExpr::number(0);
    sum = SymbolicExpr::add(sum, SymbolicExpr::variable("a"));
    sum = SymbolicExpr::add(sum, SymbolicExpr::variable("b"));
    auto sim = sum->simplify();

    auto parsed_expected = parse_expr("a+b");
    ASSERT_TRUE(parsed_expected.has_value());
    auto expected = parsed_expected.value()->simplify();
    EXPECT_TRUE(test_same_expression(sim, expected));
    EXPECT_TRUE(test_proved_equivalent(sim, expected)) << "normalization preserves both symbolic summands";
}
