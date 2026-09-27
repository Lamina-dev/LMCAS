#include "test_common.hpp"

using namespace LMCAS;

bool symbolic_equal_check(const std::shared_ptr<SymbolicExpr> &a, const std::shared_ptr<SymbolicExpr> &b) {
    auto diff = SymbolicExpr::add(a, SymbolicExpr::multiply(SymbolicExpr::number(-1), b))->expand();

    auto simp = diff->simplify();
    return simp->is_number() && simp->convert_rational() == Rational(0);
}

TEST(LmcasNorm, CanonicalOrdering) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto e1 = SymbolicExpr::add(x, y)->simplify();
    auto e2 = SymbolicExpr::add(y, x)->simplify();

    EXPECT_EQ(e1->to_string(), e2->to_string()) << "x+y == y+x (string check)";
}

TEST(LmcasNorm, AssociativeFlattening) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto z = SymbolicExpr::variable("z");
    auto e1 = SymbolicExpr::add(SymbolicExpr::add(x, y), z)->simplify();
    auto e2 = SymbolicExpr::add(x, SymbolicExpr::add(y, z))->simplify();
    EXPECT_EQ(e1->to_string(), e2->to_string()) << "(x+y)+z == x+(y+z) (string check)";
}

TEST(LmcasNorm, SquareExpansion) {
    auto x = SymbolicExpr::variable("x");
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);

    auto sum = SymbolicExpr::add(x, one);
    auto pow = SymbolicExpr::power(sum, two);
    auto expanded = pow->expand();

    auto x2 = SymbolicExpr::power(x, two);
    auto twox = SymbolicExpr::multiply(two, x);
    auto expected = SymbolicExpr::add(SymbolicExpr::add(x2, twox), one)->simplify();

    EXPECT_TRUE(symbolic_equal_check(expanded, expected)) << "(x+1)^2 == x^2+2x+1";
}

TEST(LmcasNorm, CombineLikeTerms) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);

    auto t1 = SymbolicExpr::multiply(two, x);
    auto t3 = SymbolicExpr::multiply(SymbolicExpr::number(3), x);
    auto sum = SymbolicExpr::add(t1, t3)->simplify();

    auto expected = SymbolicExpr::multiply(SymbolicExpr::number(5), x)->simplify();
    EXPECT_TRUE(symbolic_equal_check(sum, expected)) << "2x + 3x == 5x";
}
