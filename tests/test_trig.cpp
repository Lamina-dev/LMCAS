#include "test_common.hpp"
#include "symbolic.hpp"
#include "expr.hpp"

using namespace LMCAS;

TEST(Trig, BasicValues) {
    {

        EXPECT_EQ(LMCAS::simplify(SymbolicExpr::sin(SymbolicExpr::number(0))).value()->to_string(), "0") << "sin(0)";

        EXPECT_EQ(LMCAS::simplify(SymbolicExpr::cos(SymbolicExpr::number(0))).value()->to_string(), "1") << "cos(0)";

        EXPECT_EQ(LMCAS::simplify(SymbolicExpr::tan(SymbolicExpr::number(0))).value()->to_string(), "0") << "tan(0)";
    }
}

TEST(Trig, Parity) {
    auto x = SymbolicExpr::variable("x");
    {

        auto sin_neg = LMCAS::simplify(SymbolicExpr::sin(SymbolicExpr::multiply(SymbolicExpr::number(-1), x))).value();

        std::string s = sin_neg->to_string();
        bool ok = (s == "-1*sin(x)" || s == "-sin(x)" || s == "-1*(sin(x))");
        EXPECT_TRUE(ok) << "sin(-x) -> -sin(x), got: " << s;

        auto cos_neg = LMCAS::simplify(SymbolicExpr::cos(SymbolicExpr::multiply(SymbolicExpr::number(-1), x))).value();
        EXPECT_EQ(cos_neg->to_string(), "cos(x)") << "cos(-x)";
    }
}

TEST(Trig, PiValues) {
    auto pi_constant = LMCAS::pi();
    ASSERT_TRUE(pi_constant);
    auto pi = pi_constant.value();
    {

        EXPECT_EQ(LMCAS::simplify(SymbolicExpr::sin(pi)).value()->to_string(), "0") << "sin(pi)";

        EXPECT_EQ(LMCAS::simplify(SymbolicExpr::cos(pi)).value()->to_string(), "-1") << "cos(pi)";

        auto pi_6 = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 6)), pi);
        EXPECT_EQ(LMCAS::simplify(SymbolicExpr::sin(pi_6)).value()->to_string(), "1/2") << "sin(pi/6)";

        auto pi_3 = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 3)), pi);
        EXPECT_EQ(LMCAS::simplify(SymbolicExpr::cos(pi_3)).value()->to_string(), "1/2") << "cos(pi/3)";

        auto pi_4 = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 4)), pi);
        EXPECT_EQ(LMCAS::simplify(SymbolicExpr::tan(pi_4)).value()->to_string(), "1") << "tan(pi/4)";
    }
}
