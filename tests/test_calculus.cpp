#include "test_common.hpp"
#include "symbolic.hpp"




int main() {
    auto x = SymbolicExpr::variable("x");
    
    
    TEST_CASE("Limit Basic");
    {
        
        auto f = SymbolicExpr::add(x, SymbolicExpr::number(1));
        auto lim = f->limit("x", SymbolicExpr::number(2));
        EXPECT_EQ_STR(lim->to_string(), "3", "limit(x+1, x->2)");
    }
    
    TEST_CASE("Limit Indeterminate (L'Hopital Heuristic)");
    {
        
        
        auto num = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-1));
        auto den = SymbolicExpr::add(x, SymbolicExpr::number(-1));
        auto f = SymbolicExpr::divide(num, den);
        
        auto lim = f->limit("x", SymbolicExpr::number(1));
        EXPECT_EQ_STR(lim->to_string(), "2", "limit((x^2-1)/(x-1), x->1)");
    }

    
    TEST_CASE("Integral Polynomial");
    {
        
        auto integ = x->integrate("x");
        EXPECT_EQ_STR(integ->to_string(), "(1/2)*(x^2)", "int(x) dx"); 
        
        
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto integ2 = x2->integrate("x");
        EXPECT_EQ_STR(integ2->to_string(), "(1/3)*(x^3)", "int(x^2) dx");
    }
    
    TEST_CASE("Integral Log Rule");
    {
        
        auto x_inv = SymbolicExpr::power(x, SymbolicExpr::number(-1));
        auto integ = x_inv->integrate("x");
        EXPECT_EQ_STR(integ->to_string(), "ln(x)", "int(1/x) dx");
    }
    
    TEST_CASE("Integral Sum");
    {
        
        auto f = SymbolicExpr::add(x, SymbolicExpr::number(1));
        auto integ = f->integrate("x");
        
        
        
        EXPECT_EQ_STR(integ->to_string(), "((1/2)*(x^2))+x", "int(x+1) dx");
    }

    return TEST_REPORT();
}
