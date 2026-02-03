#include "symbolic.hpp"
#include <iostream>
#include <cassert>
#include <vector>

int g_failures = 0;

// Helper to print PASS/FAIL and track failure count
void CHECK_EQ(const std::string& name, const std::string& actual, const std::string& expected) {
    if (actual == expected) {
        std::cout << "[PASS] " << name << std::endl;
    } else {
        std::cout << "[FAIL] " << name << " | Expected: " << expected << ", Got: " << actual << std::endl;
        g_failures++;
    }
}

void CHECK_CONTAINS(const std::string& name, const std::string& result, const std::vector<std::string>& tokens) {
    bool ok = true;
    for(const auto& t : tokens) {
        if (result.find(t) == std::string::npos) {
            ok = false; break;
        }
    }
    if (ok) {
        std::cout << "[PASS] " << name << std::endl;
    } else {
        std::cout << "[FAIL] " << name << " | Result: " << result << std::endl;
        g_failures++;
    }
}

int main() {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    
    // 1. (x+y)^2 - (x-y)^2 = 4xy
    {
        auto A = SymbolicExpr::power(SymbolicExpr::add(x, y), two); 
        auto B = SymbolicExpr::power(SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), y)), two); 
        auto expr = SymbolicExpr::add(A, SymbolicExpr::multiply(SymbolicExpr::number(-1), B)); 
        
        auto expanded = expr->expand();
        CHECK_CONTAINS("Diff Squares (4xy)", expanded->to_string(), {"4", "x", "y"});
    }

    // 2. Rational Coefficients: x/2 + x/3 = 5x/6
    {
        auto half_x = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)), x);
        auto third_x = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 3)), x);
        auto sum = SymbolicExpr::add(half_x, third_x);
        
        auto result = sum->simplify();
        CHECK_EQ("Rational Coeffs", result->to_string(), "(5/6)*x");
    }

    // 3. Complex Polynomial: (x+1)^3 - x^3 - 1 = 3x^2 + 3x
    {
        auto term1 = SymbolicExpr::power(SymbolicExpr::add(x, one), SymbolicExpr::number(3));
        auto term2 = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(x, SymbolicExpr::number(3)));
        auto term3 = SymbolicExpr::number(-1);
        
        auto poly = SymbolicExpr::add(SymbolicExpr::add(term1, term2), term3);
        auto res = poly->expand();
        
        CHECK_CONTAINS("Cubic Cancel", res->simplify()->to_string(), {"3*x", "3*(x^2)"});
    }

    // 4. Zero Cancellation: (x + 1) - (x + 1)
    {
        auto p = SymbolicExpr::add(x, one);
        auto neg_p = SymbolicExpr::multiply(SymbolicExpr::number(-1), p);
        auto zero_expr = SymbolicExpr::add(p, neg_p);
        
        auto res = zero_expr->expand()->simplify();
        CHECK_EQ("Zero Cancel", res->to_string(), "0");
    }
    
    // 5. Nested Powers: (x^2)^3 -> x^6
    {
        auto p2 = SymbolicExpr::power(x, two);
        auto p6 = SymbolicExpr::power(p2, SymbolicExpr::number(3));
        auto res = p6->simplify();
        CHECK_EQ("Nested Powers", res->to_string(), "x^6");
    }

    // 6. Debug 0+x
    {
        auto z = SymbolicExpr::number(0);
        auto term = SymbolicExpr::variable("x");
        auto expr = SymbolicExpr::add(z, term);
        auto sim = expr->simplify();
        
        CHECK_EQ("Identity Add (0+x)", sim->to_string(), "x");
    }

    // 7. Debug 0 + x^2 + ...
    {
        auto sum = SymbolicExpr::number(0);
        sum = SymbolicExpr::add(sum, SymbolicExpr::variable("a"));
        sum = SymbolicExpr::add(sum, SymbolicExpr::variable("b"));
        auto sim = sum->simplify();
        // expect a+b
        CHECK_CONTAINS("Identity Chain", sim->to_string(), {"a", "b"});
    }

    return g_failures > 0 ? 1 : 0;
}
