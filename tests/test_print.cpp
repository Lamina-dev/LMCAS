#include <iostream>
#include <cassert>
#include "../include/symbolic.hpp"


int main() {
    auto x = SymbolicExpr::variable("x");
    auto n1 = SymbolicExpr::number(1);
    auto n2 = SymbolicExpr::number(2);
    
    auto sum = SymbolicExpr::add(x, n1);
    auto prod = SymbolicExpr::multiply(sum, n2);
    
    std::cout << "x + 1 = " << sum->to_string() << std::endl;
    std::cout << "(x + 1) * 2 = " << prod->to_string() << std::endl;
    
    
    assert(sum->to_string() == "(x + 1)" || sum->to_string() == "x + 1"); 
    
    std::cout << "Print test passed!" << std::endl;
    return 0;
}
