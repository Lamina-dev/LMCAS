#include "test_common.hpp"
#include "../symbolic.hpp"

// 隐式依赖简化版的 symbolic.hpp
// 测试新微积分功能的基本能力

int main() {
    auto x = SymbolicExpr::variable("x");
    
    // 极限
    TEST_CASE("Limit Basic");
    {
        // limit(x+1, x, 2) -> 3
        auto f = SymbolicExpr::add(x, SymbolicExpr::number(1));
        auto lim = f->limit("x", SymbolicExpr::number(2));
        EXPECT_EQ_STR(lim->to_string(), "3", "limit(x+1, x->2)");
    }
    
    TEST_CASE("Limit Indeterminate (L'Hopital Heuristic)");
    {
        // limit((x^2-1)/(x-1), x, 1) -> 2
        // f = (x^2-1) * (x-1)^-1
        auto num = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-1));
        auto den = SymbolicExpr::add(x, SymbolicExpr::number(-1));
        auto f = SymbolicExpr::divide(num, den);
        
        auto lim = f->limit("x", SymbolicExpr::number(1));
        EXPECT_EQ_STR(lim->to_string(), "2", "limit((x^2-1)/(x-1), x->1)");
    }

    // 积分
    TEST_CASE("Integral Polynomial");
    {
        // int(x, x) -> x^2/2
        auto integ = x->integrate("x");
        EXPECT_EQ_STR(integ->to_string(), "(1/2)*(x^2)", "int(x) dx"); // 化简格式可能不同
        
        // int(x^2, x) -> x^3/3
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto integ2 = x2->integrate("x");
        EXPECT_EQ_STR(integ2->to_string(), "(1/3)*(x^3)", "int(x^2) dx");
    }
    
    TEST_CASE("Integral Log Rule");
    {
        // int(x^-1, x) -> ln(x)
        auto x_inv = SymbolicExpr::power(x, SymbolicExpr::number(-1));
        auto integ = x_inv->integrate("x");
        EXPECT_EQ_STR(integ->to_string(), "ln(x)", "int(1/x) dx");
    }
    
    TEST_CASE("Integral Sum");
    {
        // int(x + 1, x) -> x^2/2 + x
        auto f = SymbolicExpr::add(x, SymbolicExpr::number(1));
        auto integ = f->integrate("x");
        
        // 预期: ((1/2)*(x^2))+x
        // 直接检查字符串表示，因为 simplify() 的减法检查不稳定
        EXPECT_EQ_STR(integ->to_string(), "((1/2)*(x^2))+x", "int(x+1) dx");
    }

    return TEST_REPORT();
}
