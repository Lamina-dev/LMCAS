#include "test_common.hpp"
#include "symbolic_implicit_diff.hpp"

using namespace LMCAS;

TEST(SymbolicImplicitDiff, ImplicitDiffCircle) {
    {
        // F(x,y) = x^2 + y^2 - r^2
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto r = SymbolicExpr::variable("r");
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto y2 = SymbolicExpr::power(y, SymbolicExpr::number(2));
        auto r2 = SymbolicExpr::power(r, SymbolicExpr::number(2));
        auto neg_r2 = SymbolicExpr::multiply(SymbolicExpr::number(-1), r2);
        auto F = SymbolicExpr::add(SymbolicExpr::add(x2, y2), neg_r2);

        auto result = LMCAS::implicit_diff(F, "x", "y");
        auto expected = SymbolicExpr::divide(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), x), y);
        EXPECT_TRUE(test_proved_equivalent(result, expected));
    }
}

TEST(SymbolicImplicitDiff, ImplicitDiffEllipse) {
    {
        // F(x,y) = x^2/a^2 + y^2/b^2 - 1
        // dy/dx = -F_x / F_y = -(2x/a^2) / (2y/b^2) = -(b^2 * x) / (a^2 * y)
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto a = SymbolicExpr::variable("a");
        auto b = SymbolicExpr::variable("b");
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto y2 = SymbolicExpr::power(y, SymbolicExpr::number(2));
        auto a2 = SymbolicExpr::power(a, SymbolicExpr::number(2));
        auto b2 = SymbolicExpr::power(b, SymbolicExpr::number(2));
        auto term1 = SymbolicExpr::divide(x2, a2);
        auto term2 = SymbolicExpr::divide(y2, b2);
        auto neg_one = SymbolicExpr::number(-1);
        auto F = SymbolicExpr::add(SymbolicExpr::add(term1, term2), neg_one);

        auto result = LMCAS::implicit_diff(F, "x", "y");
        auto expected = SymbolicExpr::divide(
            SymbolicExpr::multiply(
                SymbolicExpr::number(-1), SymbolicExpr::multiply(b2, x)),
            SymbolicExpr::multiply(a2, y));
        EXPECT_TRUE(test_proved_equivalent(result, expected));
    }
}

TEST(SymbolicImplicitDiff, ImplicitDiffPolynomial) {
    {
        // F(x,y) = x^3 + y^3 - 3xy
        // F_x = 3x^2 - 3y
        // F_y = 3y^2 - 3x
        // dy/dx = -(3x^2 - 3y) / (3y^2 - 3x) = -(x^2 - y) / (y^2 - x)
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto x3 = SymbolicExpr::power(x, SymbolicExpr::number(3));
        auto y3 = SymbolicExpr::power(y, SymbolicExpr::number(3));
        auto three = SymbolicExpr::number(3);
        auto xy = SymbolicExpr::multiply(x, y);
        auto three_xy = SymbolicExpr::multiply(three, xy);
        auto neg_three_xy = SymbolicExpr::multiply(SymbolicExpr::number(-1), three_xy);
        auto F = SymbolicExpr::add(SymbolicExpr::add(x3, y3), neg_three_xy);

        auto result = LMCAS::implicit_diff(F, "x", "y");
        auto numerator = SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::multiply(SymbolicExpr::number(-1), y));
        auto denominator = SymbolicExpr::add(
            SymbolicExpr::power(y, SymbolicExpr::number(2)),
            SymbolicExpr::multiply(SymbolicExpr::number(-1), x));
        auto expected = SymbolicExpr::divide(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), numerator),
            denominator);
        EXPECT_TRUE(test_proved_equivalent(result, expected));
    }
}

TEST(SymbolicImplicitDiff, ImplicitDiffTranscendental) {
    {
        // F(x,y) = sin(x) + y^2
        // F_x = cos(x)
        // F_y = 2y
        // dy/dx = -cos(x) / (2y)
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto sin_x = SymbolicExpr::sin(x);
        auto y2 = SymbolicExpr::power(y, SymbolicExpr::number(2));
        auto F = SymbolicExpr::add(sin_x, y2);

        auto result = LMCAS::implicit_diff(F, "x", "y");
        auto expected = SymbolicExpr::divide(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::cos(x)),
            SymbolicExpr::multiply(SymbolicExpr::number(2), y));
        EXPECT_TRUE(test_proved_equivalent(result, expected));
    }
}
