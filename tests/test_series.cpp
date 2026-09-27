#include <iostream>
#include <vector>
#include "expr.hpp"
#include "symbolic.hpp"
#include "test_common.hpp"

using namespace LMCAS;

TEST(LmcasSeries, MaclaurinSin) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);

    auto series_2 = expr->series("x", SymbolicExpr::number(0), 2);

    EXPECT_TRUE(test_proved_equivalent(series_2, x));

    auto series_5 = expr->series("x", SymbolicExpr::number(0), 5);

    ASSERT_NE(series_5, nullptr);
    EXPECT_TRUE(test_proved_equivalent(series_5, parse_expr("x-x^3/6+x^5/120").value()));

    auto expanded = series_5->expand();

    EXPECT_TRUE(test_proved_equivalent(expanded, parse_expr("x-x^3/6+x^5/120").value()));
}

TEST(LmcasSeries, MaclaurinExp) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::exp(x);

    auto series = expr->series("x", SymbolicExpr::number(0), 4);

    EXPECT_TRUE(test_proved_equivalent(series, parse_expr("1+x+x^2/2+x^3/6+x^4/24").value()));
}

TEST(LmcasSeries, TaylorLn) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::ln(x);

    auto series = expr->series("x", SymbolicExpr::number(1), 3);

    EXPECT_TRUE(test_proved_equivalent(series, parse_expr("(x-1)-(x-1)^2/2+(x-1)^3/3").value()));
}

TEST(LmcasSeries, PolySeries) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), x), SymbolicExpr::number(1)));

    auto series = expr->series("x", SymbolicExpr::number(0), 3);

    ASSERT_NE(series, nullptr);
    EXPECT_TRUE(test_proved_equivalent(series, expr));
    auto expanded = series->expand();

    EXPECT_TRUE(test_proved_equivalent(expanded, expr));
}
