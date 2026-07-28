#include "test_common.hpp"
#include "solve_transcendental.hpp"

int main() {
    TEST_CASE("Solve Transcendental - Basic callable");

    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);

    auto roots = lamina::solve_transcendental(expr, "x");
    EXPECT_TRUE(!roots.empty(), "solve_transcendental returns solutions for sin(x)=0");

    TEST_CASE("Solve Transcendental - numeric expression RHS");
    auto two_as_expr = SymbolicExpr::add(SymbolicExpr::number(1), SymbolicExpr::number(1));
    auto exp_eq = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), two_as_expr));
    auto exp_roots = lamina::solve_transcendental(exp_eq, "x");
    EXPECT_TRUE(!exp_roots.empty(),
                "solve_transcendental accepts finite numeric expression RHS");
    if (!exp_roots.empty()) {
        auto expected = SymbolicExpr::ln(SymbolicExpr::number(2))->simplify();
        EXPECT_EQ_EXPR(exp_roots.front()->simplify(), expected, "exp(x)=1+1 gives ln(2)");
    }

    return TEST_REPORT();
}
