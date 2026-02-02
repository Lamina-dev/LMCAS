#include "../symbolic.hpp"
#include <iostream>
#include <cassert>
#include <string>
#include <algorithm>

// Simple helper to check if string contains substring
bool contains(const std::string& s, const std::string& sub) {
    return s.find(sub) != std::string::npos;
}

int main() {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto two = SymbolicExpr::number(2);
    
    // --- Test 1: Common factor 2x + 2y ---
    {
        auto t1 = SymbolicExpr::multiply(two, x);
        auto t2 = SymbolicExpr::multiply(two, y);
        auto expr1 = SymbolicExpr::add(t1, t2); // 2x + 2y
        
        std::cout << "Test 1 Input: " << expr1->to_string() << std::endl;
        auto factored1 = expr1->factor();
        std::cout << "Test 1 Factored: " << factored1->to_string() << std::endl;
        
        bool passed1 = false;
        // Expect 2 * (x+y) or similar
        std::string s = factored1->to_string();
        if (contains(s, "2*(") || contains(s, "2 (")) {
             passed1 = true;
        }
        
        if (passed1) std::cout << "[PASS] Test 1" << std::endl;
        else std::cout << "[FAIL] Test 1" << std::endl;
    }

    // --- Test 2: Quadratic x^2 + 5x + 6 ---
    {
        // x^2
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        // 5x
        auto x5 = SymbolicExpr::multiply(SymbolicExpr::number(5), x);
        // 6
        auto n6 = SymbolicExpr::number(6);
        
        std::vector<std::shared_ptr<SymbolicExpr>> terms = {x2, x5, n6};
        auto expr2 = SymbolicExpr::add(terms[0], SymbolicExpr::add(terms[1], terms[2]));
        
        std::cout << "Test 2 Input: " << expr2->to_string() << std::endl;
        auto factored2 = expr2->factor();
        std::cout << "Test 2 Factored: " << factored2->to_string() << std::endl;
        
        // Expect (x+2)*(x+3) or (x+3)*(x+2)
        std::string s = factored2->to_string();
        // Check contains (x+2) and (x+3)
        bool passed2 = contains(s, "x+2") && contains(s, "x+3");
        
        if (passed2) std::cout << "[PASS] Test 2" << std::endl;
        else std::cout << "[FAIL] Test 2 (Expected (x+2)(x+3))" << std::endl;
    }
    
    // --- Test 3: Quadratic Difference of Squares x^2 - 4 ---
    {
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto n4 = SymbolicExpr::number(-4);
        auto expr3 = SymbolicExpr::add(x2, n4);
        
        std::cout << "Test 3 Input: " << expr3->to_string() << std::endl;
        auto factored3 = expr3->factor();
        std::cout << "Test 3 Factored: " << factored3->to_string() << std::endl;
        
        // Expect (x+2)*(x-2) or (x+2)*(x+-2)
        std::string s = factored3->to_string();
        bool passed3 = contains(s, "x+2") && (contains(s, "x-2") || contains(s, "x+-2"));
        if (passed3) std::cout << "[PASS] Test 3" << std::endl;
        else std::cout << "[FAIL] Test 3 (Expected (x-2)(x+2))" << std::endl;
    }

    return 0;
}
