#include "cas.hpp"
#include "symbolic.hpp"
#include "value.hpp"
#include <iostream>
#include <vector>
#include <string>

int g_failures = 0;

void CHECK_EQ(const std::string& name, const std::string& actual, const std::string& expected) {
    if (actual == expected) {
        std::cout << "[PASS] " << name << std::endl;
    } else {
        std::cout << "[FAIL] " << name << " | Expected: " << expected << ", Got: " << actual << std::endl;
        g_failures++;
    }
}

void CHECK_CONTAINS(const std::string& name, const std::string& result, const std::vector<std::string>& tokens) {
    bool ok = true;
    for(const auto& t : tokens) {
        if (result.find(t) == std::string::npos) {
            ok = false; break;
        }
    }
    if (ok) {
        std::cout << "[PASS] " << name << std::endl;
    } else {
        std::cout << "[FAIL] " << name << " | Result: " << result << std::endl;
        g_failures++;
    }
}

int main() {
    // Test 1: Simplify 1 + 1
    {
        auto one = SymbolicExpr::number(1);
        auto expr1 = SymbolicExpr::add(one, one);
        std::vector<Value> args = { Value(expr1) };
        Value res = cas_simplify(args);
        CHECK_EQ("Simplify 1+1", res.to_string(), "2");
    }

    // Test 2: Sqrt(4)
    {
        auto four = SymbolicExpr::number(4);
        auto expr = SymbolicExpr::sqrt(four);
        std::vector<Value> args = { Value(expr) };
        Value res = cas_simplify(args);
        CHECK_EQ("Sqrt(4)", res.to_string(), "2");
    }

    // Test 3: Sqrt(8) -> 2*sqrt(2) or 2√2
    {
        auto eight = SymbolicExpr::number(8);
        auto expr = SymbolicExpr::sqrt(eight);
        std::vector<Value> args = { Value(expr) };
        Value res = cas_simplify(args);
        // Result is typically 2√2 or 2*sqrt(2)
        // Check for '2' and ensuring it expanded (no '8')
        std::string s = res.to_string();
        bool has_2 = s.find("2") != std::string::npos;
        bool no_8 = s.find("8") == std::string::npos;
        if (has_2 && no_8) {
             std::cout << "[PASS] Sqrt(8)" << std::endl;
        } else {
             std::cout << "[FAIL] Sqrt(8) | Result: " << s << std::endl;
             g_failures++;
        }
    }

    // Test 4: BigInt
    {
        auto num = SymbolicExpr::number(BigInt("123456789123456789"));
        auto expr = SymbolicExpr::multiply(num, num);
        std::vector<Value> args = { Value(expr) };
        Value res = cas_simplify(args);
        // Check prefix 152415787
        CHECK_CONTAINS("BigInt Mul", res.to_string(), {"152415787"}); 
    }

    auto x = SymbolicExpr::variable("x");

    // Test 5: Diff x^2
    {
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        std::vector<Value> args = { Value(x2), Value("x") };
        Value res = cas_differentiate(args);
        CHECK_EQ("Diff x^2", res.to_string(), "2*x");
    }

    // Test 6: Diff sin(x)
    {
        auto sin_x = SymbolicExpr::sin(x);
        std::vector<Value> args = { Value(sin_x), Value("x") };
        Value res = cas_differentiate(args);
        CHECK_EQ("Diff sin(x)", res.to_string(), "cos(x)");
    }

    // Test 8: Solve 2x - 6 = 0
    {
        auto eq = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::number(-6)
        );
        std::vector<Value> args = { Value(eq), Value("x") };
        Value res = cas_solve(args);
        // Expect list { 3 }
        CHECK_CONTAINS("Solve 2x-6=0", res.to_string(), {"3"});
    }

    // Test 9: Matrix Addition
    {
        auto m1 = SymbolicExpr::matrix({
            {SymbolicExpr::number(1), SymbolicExpr::number(2)},
            {SymbolicExpr::number(3), SymbolicExpr::number(4)}
        });
        auto m2 = SymbolicExpr::matrix({
            {SymbolicExpr::number(3), SymbolicExpr::number(4)},
            {SymbolicExpr::number(5), SymbolicExpr::number(6)}
        });
        auto mat_add = SymbolicExpr::add(m1, m2);
        std::vector<Value> args = { Value(mat_add) };
        Value res = cas_simplify(args);
        // [4, 6]
        // [8, 10]
        CHECK_CONTAINS("Matrix Add", res.to_string(), {"4", "6", "8", "10"});
    }

    // Test 13: Logarithm
    {
        auto x_squared = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto log_expr = SymbolicExpr::log(x_squared, x);
        std::vector<Value> args = { Value(log_expr) };
        Value res = cas_simplify(args);
        CHECK_EQ("Log(x^2, x)", res.to_string(), "2");
    }

    return g_failures > 0 ? 1 : 0;
}
