#include <iostream>
#include <cassert>
#include "../include/symbolic.hpp"

// Simple test for PrintVisitor
int main() {
    auto x = SymbolicExpr::variable("x");
    auto n1 = SymbolicExpr::number(1);
    auto n2 = SymbolicExpr::number(2);
    
    auto sum = SymbolicExpr::add(x, n1);
    auto prod = SymbolicExpr::multiply(sum, n2);
    
    std::cout << "x + 1 = " << sum->to_string() << std::endl;
    std::cout << "(x + 1) * 2 = " << prod->to_string() << std::endl;
    
    // Check specific output format
    assert(sum->to_string() == "(x + 1)" || sum->to_string() == "x + 1"); // Parentheses might depend on implementation details
    
    std::cout << "Print test passed!" << std::endl;
    return 0;
}
