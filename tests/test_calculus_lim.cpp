#include "limit_result.hpp"
#include "symbolic.hpp"
#include "test_common.hpp"

using namespace LMCAS;

TEST(LmcasCalculusLim, PolynomialPrimitive) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto integrated = expr->integrate("x");
    ASSERT_TRUE((integrated != nullptr)) << ("integral of x^2 exists");
}

TEST(LmcasCalculusLim, PolynomialLimit) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto point = SymbolicExpr::number(2);
    auto lim = LMCAS::limit_expression_checked(expr, "x", point).value();
    EXPECT_TRUE((lim != nullptr && std::abs(lim->to_numeric() - 4.0) < 1e-12)) << ("limit x^2 at 2 equals 4");
}

TEST(LmcasCalculusLim, ConstantSumSimplification) {
    auto x = SymbolicExpr::variable("x");
    auto five = SymbolicExpr::number(5);
    auto sum = SymbolicExpr::add(five, x);
    auto simple_sum = sum->simplify();
    EXPECT_TRUE(test_proved_equivalent(simple_sum, sum));
}

TEST(LmcasCalculusLim, RemovableSingularity) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto four = SymbolicExpr::number(4);
    auto two = SymbolicExpr::number(2);
    auto neg_four = SymbolicExpr::multiply(four, SymbolicExpr::number(-1));
    auto neg_two = SymbolicExpr::multiply(two, SymbolicExpr::number(-1));
    auto num = SymbolicExpr::add(expr, neg_four);
    auto den = SymbolicExpr::add(x, neg_two);
    auto minus_one = SymbolicExpr::number(-1);
    auto rational = SymbolicExpr::multiply(num, SymbolicExpr::power(den, minus_one));
    auto lim_rational = LMCAS::limit_expression_checked(rational, "x", two).value();
    EXPECT_TRUE((lim_rational != nullptr && std::abs(lim_rational->to_numeric() - 4.0) < 1e-12)) << ("removable singularity limit equals 4");
}
