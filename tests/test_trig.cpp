#include "test_common.hpp"
#include "symbolic.hpp"

using namespace LMCAS;

int main() {
    auto x = SymbolicExpr::variable("x");
    auto pi = SymbolicExpr::variable("pi");

    TEST_CASE("Trig Basic Values");
    {

        EXPECT_EQ_STR(SymbolicExpr::sin(SymbolicExpr::number(0))->simplify()->to_string(), "0", "sin(0)");

        EXPECT_EQ_STR(SymbolicExpr::cos(SymbolicExpr::number(0))->simplify()->to_string(), "1", "cos(0)");

        EXPECT_EQ_STR(SymbolicExpr::tan(SymbolicExpr::number(0))->simplify()->to_string(), "0", "tan(0)");
    }

    TEST_CASE("Trig Parity");
    {

        auto sin_neg = SymbolicExpr::sin(SymbolicExpr::multiply(SymbolicExpr::number(-1), x))->simplify();

        std::string s = sin_neg->to_string();
        bool ok = (s == "-1*sin(x)" || s == "-sin(x)" || s == "-1*(sin(x))");
        if(!ok) std::cout << "  Got: " << s << std::endl;
        EXPECT_TRUE(ok, "sin(-x) -> -sin(x)");

        auto cos_neg = SymbolicExpr::cos(SymbolicExpr::multiply(SymbolicExpr::number(-1), x))->simplify();
        EXPECT_EQ_STR(cos_neg->to_string(), "cos(x)", "cos(-x)");
    }

    TEST_CASE("Trig Pi Values");
    {

        EXPECT_EQ_STR(SymbolicExpr::sin(pi)->simplify()->to_string(), "0", "sin(pi)");

        EXPECT_EQ_STR(SymbolicExpr::cos(pi)->simplify()->to_string(), "-1", "cos(pi)");

        auto pi_6 = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1,6)), pi);
        EXPECT_EQ_STR(SymbolicExpr::sin(pi_6)->simplify()->to_string(), "1/2", "sin(pi/6)");

        auto pi_3 = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1,3)), pi);
        EXPECT_EQ_STR(SymbolicExpr::cos(pi_3)->simplify()->to_string(), "1/2", "cos(pi/3)");

        auto pi_4 = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1,4)), pi);
        EXPECT_EQ_STR(SymbolicExpr::tan(pi_4)->simplify()->to_string(), "1", "tan(pi/4)");
    }

    return TEST_REPORT();
}
