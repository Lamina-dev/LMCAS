#include "cas.hpp"
#include "symbolic.hpp"
#include "value.hpp"
#include <iostream>
#include <vector>

void print_value(const Value& v) {
    std::cout << "Result: " << v.to_string() << std::endl;
}

int main() {
    std::cout << "=== Testing LMCAS ===" << std::endl;

    // Test 1: Simplify 1 + 1
    // Manually construct expression because parser is TODO
    // 1 + 1
    auto one = SymbolicExpr::number(1);
    auto expr1 = SymbolicExpr::add(one, one);
    
    std::vector<Value> args;
    args.push_back(Value(expr1));
    
    std::cout << "[Test 1] Simplifying 1 + 1..." << std::endl;
    Value res1 = cas_simplify(args);
    print_value(res1); // Should be 2

    // Test 2: Sqrt(4)
    auto four = SymbolicExpr::number(4);
    auto expr2 = SymbolicExpr::sqrt(four);
    args[0] = Value(expr2);
    
    std::cout << "[Test 2] Simplifying sqrt(4)..." << std::endl;
    Value res2 = cas_simplify(args);
    print_value(res2); // Should be 2

    // Test 3: Sqrt(8) -> 2*sqrt(2)
    auto eight = SymbolicExpr::number(8);
    auto expr3 = SymbolicExpr::sqrt(eight);
    args[0] = Value(expr3);

    std::cout << "[Test 3] Simplifying sqrt(8)..." << std::endl;
    Value res3 = cas_simplify(args);
    print_value(res3); // Should be 2*(2)  or similar string repr

    // Test 4: BigInt operations
    // 100! or just large number multiplication
    // Since we don't have factorial exposed easily, let's do 123456789 * 123456789
    auto num = SymbolicExpr::number(BigInt("123456789123456789"));
    auto expr4 = SymbolicExpr::multiply(num, num);
    args[0] = Value(expr4);
    
    std::cout << "[Test 4] BigInt multiplication..." << std::endl;
    Value res4 = cas_simplify(args);
    print_value(res4);

    return 0;
}
