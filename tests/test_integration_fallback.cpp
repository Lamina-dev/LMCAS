#include "test_common.hpp"
#include <iostream>

int main() {
    auto y = SymbolicExpr::variable("y");
    auto res = y->integrate("x");
    std::cout << "Integration result: " << res->to_string() << std::endl;
    
    // y is constant w.r.t. x, so integral(y, x) = y*x
    // The result should contain both y and x (as a product)
    EXPECT_CONTAINS(res->to_string(), {"y", "x"}, "Integration of y w.r.t x should be y*x");
    
    return TEST_REPORT();
}