#include "test_common.hpp"

// Implicit dependency on proper symbolic.hpp, value.hpp, etc. from test_common inclusion or relative path.
// test_common.hpp includes ../symbolic.hpp, which should suffice if headers are self-contained.
// The original file included "cas.hpp" and "value.hpp" explicitly.
#include "../cas.hpp"
#include "../value.hpp"

int main() {
    TEST_CASE("Prove: log_b(x) = ln(x)/ln(b)");

    // 1. Define x and b
    auto x = SymbolicExpr::variable("x");
    auto b = SymbolicExpr::variable("b");

    // 2. LHS: log_b(x)
    auto lhs = SymbolicExpr::log(x, b);
    std::cout << "LHS: " << lhs->to_string() << std::endl;

    // 3. RHS: ln(x) / ln(b)
    auto ln_x = SymbolicExpr::ln(x);
    auto ln_b = SymbolicExpr::ln(b);
    auto rhs = SymbolicExpr::multiply(ln_x, SymbolicExpr::power(ln_b, SymbolicExpr::number(-1)));
    std::cout << "RHS: " << rhs->to_string() << std::endl;

    // 4. Simplify
    auto lhs_simplified = lhs->simplify();
    std::cout << "LHS Simplified: " << lhs_simplified->to_string() << std::endl; // Should be ln(x)/ln(b)

    auto rhs_simplified = rhs->simplify();
    std::cout << "RHS Simplified: " << rhs_simplified->to_string() << std::endl;

    // 5. Verify equality
    bool string_match = (lhs_simplified->to_string() == rhs_simplified->to_string());
    
    // Difference check: LHS - RHS = 0
    auto diff = SymbolicExpr::add(lhs_simplified, SymbolicExpr::multiply(rhs_simplified, SymbolicExpr::number(-1)));
    auto diff_simplified = diff->simplify();
    bool diff_is_zero = (diff_simplified->is_number() && diff_simplified->convert_rational() == Rational(0));

    EXPECT_TRUE(string_match || diff_is_zero, "log_b(x) == ln(x)/ln(b)");
    
    return TEST_REPORT();
}
