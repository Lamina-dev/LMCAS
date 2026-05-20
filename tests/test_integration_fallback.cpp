#include "test_common.hpp"
#include <iostream>
#include <string>

int main() {
    auto y = SymbolicExpr::variable("y");
    auto res = y->integrate("x");
    std::string s = res->to_string();
    std::cout << "Integration result: " << s << std::endl;
    
    // y is constant w.r.t. x, so integral(y, x) must be y*x (or x*y).
    // Verify the expected product appears AND that the result is not an
    // unevaluated integral (no "integral(" or "Integral(" substring).
    bool has_xy = (s.find("y*x") != std::string::npos)
                || (s.find("x*y") != std::string::npos);
    bool has_unevaluated = (s.find("integral(") != std::string::npos)
                        || (s.find("Integral(") != std::string::npos);
    
    EXPECT_TRUE(has_xy, "Integration of y w.r.t x should produce a product of x and y");
    EXPECT_TRUE(!has_unevaluated, "Result should not be an unevaluated integral");
    
    return TEST_REPORT();
}