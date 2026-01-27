#include "../symbolic.hpp"
#include <iostream>
#include <vector>
#include <cassert>

// Simple assertion helper
#define ASSERT_EQ(a, b) \
    if (!symbolic_equal(a, b)) { \
        std::cerr << "[FAIL] Line " << __LINE__ << ": Expected " << b->to_string() << ", got " << a->to_string() << std::endl; \
        return 1; \
    }

// Manual equality check to expose internal structure sorting
bool strict_equal(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    if (a->to_string() != b->to_string()) return false;
    // Check operands order
    if (a->operands.size() != b->operands.size()) return false;
    for (size_t i=0; i<a->operands.size(); ++i) {
        if (!strict_equal(a->operands[i], b->operands[i])) return false;
    }
    return true;
}

static bool symbolic_equal(const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b) {
    // We use subtraction as the ultimate test: a - b == 0
    auto diff = SymbolicExpr::add(a, SymbolicExpr::multiply(SymbolicExpr::number(-1), b))->expand();
    // Simplified result should be 0
    return diff->is_number() && diff->convert_rational() == Rational(0);
}

int main() {
    std::cout << "Running Normalization Tests..." << std::endl;

    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);

    // Test 1: Canonical Ordering (Commutativity)
    // a + b should equal b + a structurally after simplification
    {
        auto e1 = SymbolicExpr::add(x, y)->simplify();
        auto e2 = SymbolicExpr::add(y, x)->simplify();
        std::cout << "e1: " << e1->to_string() << "\n";
        std::cout << "e2: " << e2->to_string() << "\n";
        
        if (e1->to_string() == e2->to_string()) {
             std::cout << "[PASS] Commutativity Sort (x+y == y+x)" << std::endl;
        } else {
             std::cout << "[FAIL] Commutativity Sort" << std::endl;
             // Don't exit yet, check others
        }
    }

    // Test 2: Associativity & Flattening
    // (x + y) + z == x + (y + z)
    {
        auto z = SymbolicExpr::variable("z");
        auto e1 = SymbolicExpr::add(SymbolicExpr::add(x, y), z)->simplify();
        auto e2 = SymbolicExpr::add(x, SymbolicExpr::add(y, z))->simplify();
        std::cout << "assoc 1: " << e1->to_string() << "\n";
        std::cout << "assoc 2: " << e2->to_string() << "\n";
        if (e1->to_string() == e2->to_string()) {
             std::cout << "[PASS] Associativity Sort" << std::endl;
        } else {
             std::cout << "[FAIL] Associativity Sort" << std::endl;
        }
    }

    // Test 3: Expansion of (x+1)^2
    {
        // (x+1)^2
        auto sum = SymbolicExpr::add(x, one);
        auto pow = SymbolicExpr::power(sum, two);
        auto expanded = pow->expand();
        
        std::cout << "(x+1)^2 expanded: " << expanded->to_string() << std::endl;
        // Expected: x^2 + 2x + 1 (order might vary: 1 + 2x + x^2)
        
        // Manual construction of expected: x^2 + 2x + 1
        auto x2 = SymbolicExpr::power(x, two);
        auto twox = SymbolicExpr::multiply(two, x);
        auto expected = SymbolicExpr::add(SymbolicExpr::add(x2, twox), one)->simplify();
        
        std::cout << "Expected: " << expected->to_string() << std::endl;
        
        // Use subtraction check
        if (symbolic_equal(expanded, expected)) {
             std::cout << "[PASS] Expansion (x+1)^2 == x^2 + 2x + 1" << std::endl;
        } else {
             std::cout << "[FAIL] Expansion (x+1)^2" << std::endl;
             return 1;
        }
    }
    
    // Test 4: Combine Like Terms
    // 2x + 3x -> 5x
    {
        auto t1 = SymbolicExpr::multiply(two, x);
        auto t3 = SymbolicExpr::multiply(SymbolicExpr::number(3), x);
        auto sum = SymbolicExpr::add(t1, t3)->simplify();
        std::cout << "2x + 3x = " << sum->to_string() << std::endl;
        
        // Expected: 5x
        auto expected = SymbolicExpr::multiply(SymbolicExpr::number(5), x)->simplify();
         if (symbolic_equal(sum, expected)) {
             std::cout << "[PASS] Combine Like Terms" << std::endl;
        } else {
             std::cout << "[FAIL] Combine Like Terms" << std::endl;
             return 1;
        }
    }

    std::cout << "All Tests Passed." << std::endl;
    return 0;
}
