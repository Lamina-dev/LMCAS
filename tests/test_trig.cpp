#include "test_common.hpp"
#include "symbolic.hpp"

int main() {
    auto x = SymbolicExpr::variable("x");
    auto pi = SymbolicExpr::variable("pi"); // Assuming "pi" or "π" is recognized

    TEST_CASE("Trig Basic Values");
    {
        // sin(0) -> 0
        EXPECT_EQ_STR(SymbolicExpr::sin(SymbolicExpr::number(0))->simplify()->to_string(), "0", "sin(0)");
        // cos(0) -> 1
        EXPECT_EQ_STR(SymbolicExpr::cos(SymbolicExpr::number(0))->simplify()->to_string(), "1", "cos(0)");
        // tan(0) -> 0
        EXPECT_EQ_STR(SymbolicExpr::tan(SymbolicExpr::number(0))->simplify()->to_string(), "0", "tan(0)");
    }

    TEST_CASE("Trig Parity");
    {
        // sin(-x) -> -1*sin(x) or -sin(x) depending on output format
        // The implementation specifically returns -1 * sin(x)
        auto sin_neg = SymbolicExpr::sin(SymbolicExpr::multiply(SymbolicExpr::number(-1), x))->simplify();
        // to_string for -1*sin(x) typically is "-1*sin(x)" or "-sin(x)"?
        // Let's check logic: SymbolicExpr::multiply(SymbolicExpr::number(-1), ... )
        // The toString might be "-1*sin(x)"
        // Let's try to verify if it contains sin(x) and starts with -
        std::string s = sin_neg->to_string();
        bool ok = (s == "-1*sin(x)" || s == "-sin(x)" || s == "-1*(sin(x))");
        if(!ok) std::cout << "  Got: " << s << std::endl;
        EXPECT_TRUE(ok, "sin(-x) -> -sin(x)");

        // cos(-x) -> cos(x)
        auto cos_neg = SymbolicExpr::cos(SymbolicExpr::multiply(SymbolicExpr::number(-1), x))->simplify();
        EXPECT_EQ_STR(cos_neg->to_string(), "cos(x)", "cos(-x)");
    }

    TEST_CASE("Trig Pi Values");
    {
        // sin(pi) -> 0
        EXPECT_EQ_STR(SymbolicExpr::sin(pi)->simplify()->to_string(), "0", "sin(pi)");
        // cos(pi) -> -1
        EXPECT_EQ_STR(SymbolicExpr::cos(pi)->simplify()->to_string(), "-1", "cos(pi)");
        
        // sin(pi/6) -> 1/2
        // pi/6 is pi * (1/6)
        auto pi_6 = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1,6)), pi);
        EXPECT_EQ_STR(SymbolicExpr::sin(pi_6)->simplify()->to_string(), "1/2", "sin(pi/6)");
        
        // cos(pi/3) -> 1/2
        auto pi_3 = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1,3)), pi);
        EXPECT_EQ_STR(SymbolicExpr::cos(pi_3)->simplify()->to_string(), "1/2", "cos(pi/3)");
        
        // tan(pi/4) -> 1
        auto pi_4 = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1,4)), pi);
        EXPECT_EQ_STR(SymbolicExpr::tan(pi_4)->simplify()->to_string(), "1", "tan(pi/4)");
    }

    return 0;
}
