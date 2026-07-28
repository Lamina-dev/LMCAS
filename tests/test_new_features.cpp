#include "test_common.hpp"

int main() {
    TEST_CASE("New Features (Expand, GCD, Solve, Resultant)");
    auto x = SymbolicExpr::variable("x");

    {
        auto expr1 = SymbolicExpr::power(SymbolicExpr::add(x, SymbolicExpr::number(1)), SymbolicExpr::number(3));
        auto expanded = expr1->expand();

        EXPECT_CONTAINS(expanded->to_string(), {"x^3", "3*(x^2)", "3*x", "1"}, "Expand (x+1)^3");
    }

    {

        auto P = SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), x), SymbolicExpr::number(1))
        );

        auto Q = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-1));

        auto gcd = SymbolicExpr::poly_gcd(P, Q);

        EXPECT_CONTAINS(gcd->to_string(), {"x", "1"}, "GCD(x^2+2x+1, x^2-1)");
    }

    {

        auto eq = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-4));
        auto sol = SymbolicExpr::solve(eq, "x");

        EXPECT_TRUE(!sol.empty(), "Solve x^2-4=0 returned empty");

        bool has2 = false, hasNeg2 = false;
        for(auto& s : sol) {
            if (s->to_string() == "2") has2 = true;
            if (s->to_string() == "-2") hasNeg2 = true;
        }

        EXPECT_TRUE(has2, "Root 2 found");
        EXPECT_TRUE(hasNeg2, "Root -2 found");
    }

    {
        auto y = SymbolicExpr::variable("y");

        auto eq1 = SymbolicExpr::add(SymbolicExpr::add(x, y), SymbolicExpr::number(-3));

        auto eq2 = SymbolicExpr::add(SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), y)), SymbolicExpr::number(-1));

        auto sys_sol = SymbolicExpr::solve_system({eq1, eq2}, {"x", "y"});
        EXPECT_TRUE(!sys_sol.empty(), "System Solver returned solution");

        if (!sys_sol.empty()) {
             EXPECT_TRUE(sys_sol[0].count("x") > 0, "x in solution");
             EXPECT_TRUE(sys_sol[0].count("y") > 0, "y in solution");
             if (sys_sol[0].count("x") && sys_sol[0].count("y")) {
                 EXPECT_EQ_EXPR(sys_sol[0].at("x"), SymbolicExpr::number(2), "x=2");
                 EXPECT_EQ_EXPR(sys_sol[0].at("y"), SymbolicExpr::number(1), "y=1");
             }
        }
    }

    {

        auto res_poly1 = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-1));
        auto res_poly2 = SymbolicExpr::add(x, SymbolicExpr::number(1));

        auto resultant = SymbolicExpr::poly_resultant(res_poly1, res_poly2, "x");
        EXPECT_EQ_EXPR_STR(resultant, "0", "Res(x^2-1, x+1)");

        auto res_poly3 = SymbolicExpr::add(x, SymbolicExpr::number(-2));
        auto resultant2 = SymbolicExpr::poly_resultant(res_poly1, res_poly3, "x");
        EXPECT_EQ_EXPR_STR(resultant2, "3", "Res(x^2-1, x-2)");

        // 一般情形：两个二次多项式（之前只支持线性第二多项式）
        // Res(x^2-1, x^2-4) = (1-4)(1-4)... 标准值 = ((-1)-(-4))^2 * ... 用根验证：
        // x^2-1 根 ±1；Res = lead(b)^deg(a) ∏ b(根) = 1 * (1-4)(1-4) = 9
        auto res_q1 = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-1));
        auto res_q2 = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(-4));
        auto resultant3 = SymbolicExpr::poly_resultant(res_q1, res_q2, "x");
        EXPECT_EQ_EXPR_STR(resultant3, "9", "Res(x^2-1, x^2-4) = 9");

        // 共享公因子 → 结式为 0：Res(x^2-1, x^2+x-2)，二者都有根 x=1
        auto res_q3 = SymbolicExpr::add(SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), x), SymbolicExpr::number(-2));
        auto resultant4 = SymbolicExpr::poly_resultant(res_q1, res_q3, "x");
        EXPECT_EQ_EXPR_STR(resultant4, "0", "Res(x^2-1, x^2+x-2) = 0 (shared root x=1)");

        // 三次 × 二次
        // Res(x^3-1, x^2+1)：x^3-1 根含 1；x^2+1 在这些根上的积。标准值 = 2
        auto res_cubic = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(3)), SymbolicExpr::number(-1));
        auto res_quad = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(1));
        auto resultant5 = SymbolicExpr::poly_resultant(res_cubic, res_quad, "x");
        EXPECT_EQ_EXPR_STR(resultant5, "2", "Res(x^3-1, x^2+1) = 2");
    }

    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto int_x2 = x2->integrate("x");

        EXPECT_CONTAINS(int_x2->to_string(), {"x^3"}, "Integral(x^2)");

        auto sinx = SymbolicExpr::sin(x);
        auto int_sinx = sinx->integrate("x");

        EXPECT_CONTAINS(int_sinx->to_string(), {"cos(x)"}, "Integral(sin(x))");

        auto static_int = SymbolicExpr::integral(x2, "x");
        if (static_int) {
             EXPECT_CONTAINS(static_int->to_string(), {"x^3"}, "Static Integral(x^2)");
        } else {

        }

        auto sinx_x = SymbolicExpr::divide(sinx, x);

        auto zero = SymbolicExpr::number(0);
        auto lim_res = SymbolicExpr::limit_func(sinx_x, "x", zero);
        if (lim_res) {
             EXPECT_EQ_EXPR_STR(lim_res, "1", "Limit(sin(x)/x, x->0) = 1");
        }
    }

    return TEST_REPORT();
}
