#include "test_common.hpp"
#include <iostream>

int main() {
    auto expr = SymbolicExpr::variable("x");
    for (int i = 0; i < 510; ++i) {
        expr = SymbolicExpr::sin(expr);
    }

    bool caught = false;
    try {
        expr->simplify();
    } catch (const std::runtime_error& e) {
        std::string msg = e.what();
        if (msg == "AST traversal depth limit exceeded") {
            caught = true;
        }
    }

    EXPECT_TRUE(caught, "AST depth limit should trigger an exception on deep trees");

    return TEST_REPORT();
}
