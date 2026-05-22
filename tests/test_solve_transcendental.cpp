#include "test_common.hpp"
#include "solve_transcendental.hpp"

int main() {
    TEST_CASE("Solve Transcendental - Basic callable");

    // Verify the function is callable with sin(x) = 0
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x); // sin(x) = 0

    auto roots = lamina::solve_transcendental(expr, "x");
    EXPECT_TRUE(!roots.empty(), "solve_transcendental returns solutions for sin(x)=0");

    return TEST_REPORT();
}
