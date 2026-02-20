#include "test_common.hpp"





#include "value.hpp"

int main() {
    TEST_CASE("Prove: log_b(x) = ln(x)/ln(b)");

    
    auto x = SymbolicExpr::variable("x");
    auto b = SymbolicExpr::variable("b");

    
    auto lhs = SymbolicExpr::log(x, b);
    std::cout << "LHS: " << lhs->to_string() << std::endl;

    
    auto ln_x = SymbolicExpr::ln(x);
    auto ln_b = SymbolicExpr::ln(b);
    auto rhs = SymbolicExpr::multiply(ln_x, SymbolicExpr::power(ln_b, SymbolicExpr::number(-1)));
    std::cout << "RHS: " << rhs->to_string() << std::endl;

    
    auto lhs_simplified = lhs->simplify();
    std::cout << "LHS Simplified: " << lhs_simplified->to_string() << std::endl; 

    auto rhs_simplified = rhs->simplify();
    std::cout << "RHS Simplified: " << rhs_simplified->to_string() << std::endl;

    
    bool string_match = (lhs_simplified->to_string() == rhs_simplified->to_string());
    
    
    auto diff = SymbolicExpr::add(lhs_simplified, SymbolicExpr::multiply(rhs_simplified, SymbolicExpr::number(-1)));
    auto diff_simplified = diff->simplify();
    bool diff_is_zero = (diff_simplified->is_number() && diff_simplified->convert_rational() == Rational(0));

    EXPECT_TRUE(string_match || diff_is_zero, "log_b(x) == ln(x)/ln(b)");
    
    return TEST_REPORT();
}
