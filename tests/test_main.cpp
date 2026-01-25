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

    // ==========================================
    // New Calculus Tests
    // ==========================================
    std::cout << "\n=== Calculus & Solver Tests ===" << std::endl;

    auto x = SymbolicExpr::variable("x");
    
    // Test 5: Differentiation d/dx (x^2)
    // x^2
    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    std::cout << "[Test 5] Diff (x^2) wrt x..." << std::endl;
    std::cout << "Expr: " << x2->to_string() << std::endl;
    
    // cas_differentiate takes [expr, var_str]
    std::vector<Value> diff_args;
    diff_args.push_back(Value(x2));
    diff_args.push_back(Value("x"));
    
    Value res5 = cas_differentiate(diff_args);
    print_value(res5); // Expected: 2*x

    // Test 6: Differentiation d/dx (sin(x))
    auto sin_x = SymbolicExpr::sin(x);
    std::cout << "[Test 6] Diff (sin(x)) wrt x..." << std::endl;
    std::cout << "Expr: " << sin_x->to_string() << std::endl;
    
    diff_args[0] = Value(sin_x);
    Value res6 = cas_differentiate(diff_args);
    print_value(res6); // Expected: cos(x)

    // Test 7: Integration Representation int(x^2, x)
    std::cout << "[Test 7] Integrate (x^2) wrt x (Symbolic Rep)..." << std::endl;
    std::vector<Value> int_args;
    int_args.push_back(Value(x2));
    int_args.push_back(Value("x"));
    
    Value res7 = cas_integrate(int_args);
    print_value(res7); // Expected: Integral representation (not implemented execution yet)

    // Test 8: Solver linear equation 2*x + 4 = 10 -> 2*x - 6 = 0
    // Construct 2*x - 6
    auto two = SymbolicExpr::number(2);
    auto six = SymbolicExpr::number(6);
    auto eq = SymbolicExpr::add(
        SymbolicExpr::multiply(two, x),
        SymbolicExpr::number(-6) // represented as + (-6)
    );
    
    std::cout << "[Test 8] Solve 2*x - 6 = 0 for x..." << std::endl;
    std::cout << "Expr: " << eq->to_string() << " = 0" << std::endl;
    std::vector<Value> solve_args;
    solve_args.push_back(Value(eq));
    solve_args.push_back(Value("x"));
    
    // Debug output to see if it even returns
    Value res8 = cas_solve(solve_args);
    std::cout << "Solve returned." << std::endl;
    print_value(res8); // Expected: 3

    return 0;
}
