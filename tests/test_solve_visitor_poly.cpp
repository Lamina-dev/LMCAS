#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "../include/symbolic.hpp"
#include "test_common.hpp"

using namespace LMCAS;

void test_solve_numeric() {
    std::cout << "Testing Solve Numeric..." << std::endl;

    auto x = SymbolicExpr::variable("x");
    auto eq1 = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::number(-4)
    );

    auto solutions = LMCAS::solve_finite_checked(eq1, "x").value();
    EXPECT_TRUE(solutions.size() == 1, "linear equation has one solution");
    std::cout << "2x - 4 = 0 => x = " << solutions[0]->to_string() << std::endl;

    auto eq2 = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-3), x),
            SymbolicExpr::number(2)
        )
    );

    auto sol2 = LMCAS::solve_finite_checked(eq2, "x").value();
    EXPECT_TRUE(sol2.size() == 2, "quadratic equation has two solutions");
    std::cout << "x^2 - 3x + 2 = 0 => x1=" << sol2[0]->to_string() << ", x2=" << sol2[1]->to_string() << std::endl;
}

void test_solve_symbolic() {
    std::cout << "Testing Solve Symbolic..." << std::endl;
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");

    auto eq = SymbolicExpr::add(
        SymbolicExpr::multiply(a, x),
        b
    );

    auto sols = LMCAS::solve_finite_checked(eq, "x").value();
    EXPECT_TRUE(sols.size() == 1, "symbolic linear equation has one solution");
    std::cout << "ax + b = 0 => x = " << sols[0]->to_string() << std::endl;

}

int main() {
    try {
        test_solve_numeric();
        test_solve_symbolic();
        std::cout << "All solve visitor tests passed!" << std::endl;
    } catch(const std::exception& e) {
        EXPECT_TRUE(false, std::string("unexpected exception: ") + e.what());
    }
    return TEST_REPORT();
}
