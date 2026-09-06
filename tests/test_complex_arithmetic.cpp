#include "test_common.hpp"

using namespace LMCAS;

int main() {
    TEST_CASE("Complex Arithmetic");

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);

    {
        auto A = SymbolicExpr::power(SymbolicExpr::add(x, y), two);
        auto B = SymbolicExpr::power(SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), y)), two);
        auto expr = SymbolicExpr::add(A, SymbolicExpr::multiply(SymbolicExpr::number(-1), B));

        auto expanded = expr->expand();
        EXPECT_CONTAINS(expanded->to_string(), {"4", "x", "y"}, "Diff Squares (4xy)");
    }

    {
        auto half_x = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)), x);
        auto third_x = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 3)), x);
        auto sum = SymbolicExpr::add(half_x, third_x);

        auto result = sum->simplify();
        EXPECT_EQ_EXPR_STR(result, "(5/6)*x", "Rational Coeffs");
    }

    {
        auto term1 = SymbolicExpr::power(SymbolicExpr::add(x, one), SymbolicExpr::number(3));
        auto term2 = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(x, SymbolicExpr::number(3)));
        auto term3 = SymbolicExpr::number(-1);

        auto poly = SymbolicExpr::add(SymbolicExpr::add(term1, term2), term3);
        auto res = poly->expand();

        EXPECT_CONTAINS(res->simplify()->to_string(), {"3*x", "3*(x^2)"}, "Cubic Cancel");
    }

    {
        auto p = SymbolicExpr::add(x, one);
        auto neg_p = SymbolicExpr::multiply(SymbolicExpr::number(-1), p);
        auto zero_expr = SymbolicExpr::add(p, neg_p);

        auto res = zero_expr->expand()->simplify();
        EXPECT_EQ_EXPR_STR(res, "0", "Zero Cancel");
    }

    {
        auto p2 = SymbolicExpr::power(x, two);
        auto p6 = SymbolicExpr::power(p2, SymbolicExpr::number(3));
        auto res = p6->simplify();
        EXPECT_EQ_EXPR_STR(res, "x^6", "Nested Powers");
    }

    {
        auto z = SymbolicExpr::number(0);
        auto term = SymbolicExpr::variable("x");
        auto expr = SymbolicExpr::add(z, term);
        auto sim = expr->simplify();

        EXPECT_EQ_EXPR_STR(sim, "x", "Identity Add (0+x)");
    }

    {
        auto sum = SymbolicExpr::number(0);
        sum = SymbolicExpr::add(sum, SymbolicExpr::variable("a"));
        sum = SymbolicExpr::add(sum, SymbolicExpr::variable("b"));
        auto sim = sum->simplify();

        EXPECT_CONTAINS(sim->to_string(), {"a", "b"}, "Identity Chain");
    }

    return TEST_REPORT();
}
