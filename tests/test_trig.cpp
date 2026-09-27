#include "test_common.hpp"
#include "symbolic.hpp"

using namespace LMCAS;

TEST(Trig, BasicValues) {
    {

        EXPECT_EQ(SymbolicExpr::sin(SymbolicExpr::number(0))->simplify()->to_string(), "0") << "sin(0)";

        EXPECT_EQ(SymbolicExpr::cos(SymbolicExpr::number(0))->simplify()->to_string(), "1") << "cos(0)";

        EXPECT_EQ(SymbolicExpr::tan(SymbolicExpr::number(0))->simplify()->to_string(), "0") << "tan(0)";
    }
}

TEST(Trig, Parity) {
    auto x = SymbolicExpr::variable("x");
    {

        auto sin_neg = SymbolicExpr::sin(SymbolicExpr::multiply(SymbolicExpr::number(-1), x))->simplify();

        std::string s = sin_neg->to_string();
        bool ok = (s == "-1*sin(x)" || s == "-sin(x)" || s == "-1*(sin(x))");
        EXPECT_TRUE(ok) << "sin(-x) -> -sin(x), got: " << s;

        auto cos_neg = SymbolicExpr::cos(SymbolicExpr::multiply(SymbolicExpr::number(-1), x))->simplify();
        EXPECT_EQ(cos_neg->to_string(), "cos(x)") << "cos(-x)";
    }
}

TEST(Trig, PiValues) {
    auto pi = SymbolicExpr::variable("pi");
    {

        EXPECT_EQ(SymbolicExpr::sin(pi)->simplify()->to_string(), "0") << "sin(pi)";

        EXPECT_EQ(SymbolicExpr::cos(pi)->simplify()->to_string(), "-1") << "cos(pi)";

        auto pi_6 = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 6)), pi);
        EXPECT_EQ(SymbolicExpr::sin(pi_6)->simplify()->to_string(), "1/2") << "sin(pi/6)";

        auto pi_3 = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 3)), pi);
        EXPECT_EQ(SymbolicExpr::cos(pi_3)->simplify()->to_string(), "1/2") << "cos(pi/3)";

        auto pi_4 = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 4)), pi);
        EXPECT_EQ(SymbolicExpr::tan(pi_4)->simplify()->to_string(), "1") << "tan(pi/4)";
    }
}
