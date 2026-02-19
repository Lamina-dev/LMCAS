#include "../include/symbolic.hpp"
#include <iostream>
#include <cassert>

int main() {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::power(x, SymbolicExpr::number(2));
    
    // 1. Integration
    auto integrated = expr->integrate("x");
    std::cout << "int(x^2) = " << integrated->to_string() << "\n";
    
    // 2. Limit Substitution
    auto point = SymbolicExpr::number(2);
    auto lim = expr->limit("x", point);
    std::cout << "limit(x^2, x->2) = " << lim->to_string() << "\n";

    // 3. Normalization
    auto five = SymbolicExpr::number(5);
    auto sum = SymbolicExpr::add(five, x); 
    auto simple_sum = sum->simplify();
    std::cout << "5 + x simplified -> " << simple_sum->to_string() << "\n";

    // 4. LHopital Components check
    std::cout << "--- LHopital Check ---\n";
    auto four = SymbolicExpr::number(4);
    auto two = SymbolicExpr::number(2);
    auto neg_four = SymbolicExpr::multiply(four, SymbolicExpr::number(-1));
    auto neg_two = SymbolicExpr::multiply(two, SymbolicExpr::number(-1));
    
    auto num = SymbolicExpr::add(expr, neg_four); // x^2 - 4
    auto den = SymbolicExpr::add(x, neg_two); // x - 2
    
    std::cout << "Num: " << num->to_string() << "\n";
    std::cout << "Den: " << den->to_string() << "\n";
    
    auto dN = num->differentiate("x");
    auto dD = den->differentiate("x");
    std::cout << "d(Num)/dx: " << dN->to_string() << "\n";
    std::cout << "d(Den)/dx: " << dD->to_string() << "\n";
    
    auto minus_one = SymbolicExpr::number(-1);
    auto dD_inv = SymbolicExpr::power(dD, minus_one); 
    std::cout << "dD^(-1): " << dD_inv->to_string() << "\n";
    
    auto ratio = SymbolicExpr::multiply(dN, dD_inv); 
    std::cout << "Ratio dN/dD: " << ratio->to_string() << "\n";
    
    auto lim_ratio = ratio->limit("x", two);
    std::cout << "Limit of Ratio: " << lim_ratio->to_string() << "\n";

    auto rational = SymbolicExpr::multiply(num, SymbolicExpr::power(den, minus_one));
    std::cout << "Rational: " << rational->to_string() << "\n";
    auto lim_rational = rational->limit("x", two);
    std::cout << "Result: " << lim_rational->to_string() << "\n";
    
    return 0;
}
