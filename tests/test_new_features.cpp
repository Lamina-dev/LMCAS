#include "test_common.hpp"

int main() {
    TEST_CASE("New Features (Expand, GCD, Solve, Resultant)");
    auto x = SymbolicExpr::variable("x");
    
    // 1. Expand (x+1)^3
    {
        auto expr1 = SymbolicExpr::power(SymbolicExpr::add(x, SymbolicExpr::number(1)), SymbolicExpr::number(3));
        auto expanded = expr1->expand();
        // 1 + 3x^2 + x^3 + 3x
        EXPECT_CONTAINS(expanded->to_string(), {"x^3", "3*(x^2)", "3*x", "1"}, "Expand (x+1)^3");
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
        EXPECT_CONTAINS(gcd->to_string(), {"x", "1"}, "GCD(x^2+2x+1, x^2-1)");
    }
    
    // 4. Solve
    {
        // x^2 - 4 = 0
        auto eq = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-4));
        auto sol = SymbolicExpr::solve(eq, "x");

        EXPECT_TRUE(!sol.empty(), "Solve x^2-4=0 returned empty");
        
        // 2 and -2
        // We now return two roots, order may vary
        bool has2 = false, hasNeg2 = false;
        for(auto& s : sol) {
            if (s->to_string() == "2") has2 = true;
            if (s->to_string() == "-2") hasNeg2 = true;
        }
        
        EXPECT_TRUE(has2, "Root 2 found");
        EXPECT_TRUE(hasNeg2, "Root -2 found");
    }

    // 5. System
    {
        auto y = SymbolicExpr::variable("y");
        // x + y - 3 = 0
        auto eq1 = SymbolicExpr::add(SymbolicExpr::add(x, y), SymbolicExpr::number(-3));
        // x - y - 1 = 0
        auto eq2 = SymbolicExpr::add(SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), y)), SymbolicExpr::number(-1));
        
        auto sys_sol = SymbolicExpr::solve_system({eq1, eq2}, {"x", "y"});
        EXPECT_TRUE(!sys_sol.empty(), "System Solver returned solution");
        
        if (!sys_sol.empty()) {
             EXPECT_TRUE(sys_sol[0].count("x") > 0, "x in solution");
             EXPECT_TRUE(sys_sol[0].count("y") > 0, "y in solution");
             if (sys_sol[0].count("x") && sys_sol[0].count("y")) {
                 EXPECT_EQ_EXPR(sys_sol[0]["x"], SymbolicExpr::number(2), "x=2");
                 EXPECT_EQ_EXPR(sys_sol[0]["y"], SymbolicExpr::number(1), "y=1");
             }
        }
    }

    // 6. Resultant
    {
        // Res(x^2-1, x+1) should be 0
        auto res_poly1 = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-1));
        auto res_poly2 = SymbolicExpr::add(x, SymbolicExpr::number(1));
        
        auto resultant = SymbolicExpr::poly_resultant(res_poly1, res_poly2, "x");
        EXPECT_EQ_EXPR_STR(resultant, "0", "Res(x^2-1, x+1)");

        // Res(x^2-1, x-2) should be 3
        auto res_poly3 = SymbolicExpr::add(x, SymbolicExpr::number(-2));
        auto resultant2 = SymbolicExpr::poly_resultant(res_poly1, res_poly3, "x");
        EXPECT_EQ_EXPR_STR(resultant2, "3", "Res(x^2-1, x-2)");
    }

    return TEST_REPORT();
}
