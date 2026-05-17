#include "test_common.hpp"
#include <iostream>

int main() {
    auto y = SymbolicExpr::variable("y");
    auto res = y->integrate("x");
    std::cout << "Integration result: " << res->to_string() << std::endl;
    
    // Should contain Integral or similar representation based on print_visitor
    EXPECT_CONTAINS(res->to_string(), {"integral", "y", "x"}, "Integration of y w.r.t x should be unevaluated Integral");
    
    return TEST_REPORT();
}