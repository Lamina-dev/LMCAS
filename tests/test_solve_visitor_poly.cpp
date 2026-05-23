#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>
#include "../include/symbolic.hpp"

void test_solve_numeric() {
    std::cout << "Testing Solve Numeric..." << std::endl;

    auto x = SymbolicExpr::variable("x");
    auto eq1 = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::number(-4)
    );

    auto solutions = SymbolicExpr::solve(eq1, "x");
    assert(solutions.size() == 1);
    std::cout << "2x - 4 = 0 => x = " << solutions[0]->to_string() << std::endl;

    auto eq2 = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-3), x),
            SymbolicExpr::number(2)
        )
    );

    auto sol2 = SymbolicExpr::solve(eq2, "x");
    assert(sol2.size() == 2);
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

    auto sols = SymbolicExpr::solve(eq, "x");
    assert(sols.size() == 1);
    std::cout << "ax + b = 0 => x = " << sols[0]->to_string() << std::endl;

}

int main() {
    try {
        test_solve_numeric();
        test_solve_symbolic();
        std::cout << "All solve visitor tests passed!" << std::endl;
    } catch(const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
