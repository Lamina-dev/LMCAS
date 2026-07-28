#include <iostream>
#include "../include/symbolic.hpp"
#include "test_common.hpp"

int main() {
    auto x = SymbolicExpr::variable("x");
    auto n1 = SymbolicExpr::number(1);
    auto n2 = SymbolicExpr::number(2);

    auto sum = SymbolicExpr::add(x, n1);
    auto prod = SymbolicExpr::multiply(sum, n2);

    std::cout << "x + 1 = " << sum->to_string() << std::endl;
    std::cout << "(x + 1) * 2 = " << prod->to_string() << std::endl;

    EXPECT_TRUE(sum->to_string() == "(x + 1)" || sum->to_string() == "x + 1",
                "sum prints as x + 1");
    EXPECT_TRUE(prod && !prod->to_string().empty(), "product prints a non-empty expression");

    std::cout << "Print test passed!" << std::endl;
    return TEST_REPORT();
}
