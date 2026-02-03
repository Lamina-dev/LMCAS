#include "symbolic.hpp"
#include <iostream>
#include <cassert>
#include <vector>

int g_failures = 0;

void EXPECT_EQ(const std::string& actual, const std::string& expected, const std::string& test_name) {
    if (actual == expected) {
        std::cout << "[PASS] " << test_name << std::endl;
    } else {
        std::cout << "[FAIL] " << test_name << " | Expected: " << expected << ", Got: " << actual << std::endl;
        g_failures++;
    }
}

// Loose check for commutative adds
void EXPECT_WILD(const std::string& actual, const std::vector<std::string>& substrings, const std::string& test_name) {
    bool fail = false;
    for(const auto& s : substrings) {
        if (actual.find(s) == std::string::npos) {
            fail = true;
            break;
        }
    }
    if (!fail) {
        std::cout << "[PASS] " << test_name << std::endl;
    } else {
        std::cout << "[FAIL] " << test_name << " | Got: " << actual << std::endl;
        g_failures++;
    }
}

int main() {
    auto x = SymbolicExpr::variable("x");
    
    // 1. Expand (x+1)^3
    {
        auto expr1 = SymbolicExpr::power(SymbolicExpr::add(x, SymbolicExpr::number(1)), SymbolicExpr::number(3));
        auto expanded = expr1->expand();
        // 1 + 3x^2 + x^3 + 3x
        EXPECT_WILD(expanded->to_string(), {"x^3", "3*(x^2)", "3*x", "1"}, "Expand (x+1)^3");
    }
    
    // 3. GCD
    {
        // P = x^2 + 2x + 1 = (x+1)^2
        auto P = SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), x), SymbolicExpr::number(1))
        );
        // Q = x^2 - 1 = (x-1)(x+1)
        auto Q = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-1));
        
        auto gcd = SymbolicExpr::poly_gcd(P, Q);
        // Should be x+1 (or 1+x)
        EXPECT_WILD(gcd->to_string(), {"x", "1"}, "GCD(x^2+2x+1, x^2-1)");
    }
    
    // 4. Solve
    {
        // x^2 - 4 = 0
        auto eq = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-4));
        auto sol = SymbolicExpr::solve(eq, "x");
        if(sol.empty()) {
             std::cout << "[FAIL] Solve x^2-4=0 returned empty" << std::endl;
             g_failures++;
        } else {
             // 2 or -2
             std::string s = sol[0]->to_string();
             if (s == "2" || s == "-2") std::cout << "[PASS] Solve x^2-4=0" << std::endl;
             else {
                 std::cout << "[FAIL] Solve x^2-4=0 got " << s << std::endl;
                 g_failures++;
             }
        }
    }

    // 5. System
    {
        auto y = SymbolicExpr::variable("y");
        // x + y - 3 = 0
        auto eq1 = SymbolicExpr::add(SymbolicExpr::add(x, y), SymbolicExpr::number(-3));
        // x - y - 1 = 0
        auto eq2 = SymbolicExpr::add(SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), y)), SymbolicExpr::number(-1));
        
        auto sys_sol = SymbolicExpr::solve_system({eq1, eq2}, {"x", "y"});
        if (!sys_sol.empty()) {
            std::string sx = (sys_sol[0].count("x") ? sys_sol[0]["x"]->to_string() : "?");
            std::string sy = (sys_sol[0].count("y") ? sys_sol[0]["y"]->to_string() : "?");
            if (sx == "2" && sy == "1") {
                std::cout << "[PASS] System Solver" << std::endl;
            } else {
                std::cout << "[FAIL] System Solver | Got x=" << sx << ", y=" << sy << std::endl;
                g_failures++;
            }
        } else {
            std::cout << "[FAIL] System Solver failed" << std::endl;
            g_failures++;
        }
    }

    // 6. Resultant
    {
        // Res(x^2-1, x+1) should be 0
        auto res_poly1 = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-1));
        auto res_poly2 = SymbolicExpr::add(x, SymbolicExpr::number(1));
        
        auto resultant = SymbolicExpr::poly_resultant(res_poly1, res_poly2, "x");
        EXPECT_EQ(resultant->to_string(), "0", "Res(x^2-1, x+1)");

        // Res(x^2-1, x-2) should be 3
        auto res_poly3 = SymbolicExpr::add(x, SymbolicExpr::number(-2));
        auto resultant2 = SymbolicExpr::poly_resultant(res_poly1, res_poly3, "x");
        EXPECT_EQ(resultant2->to_string(), "3", "Res(x^2-1, x-2)");
    }

    return g_failures > 0 ? 1 : 0;
}
