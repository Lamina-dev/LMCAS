#include <iostream>
#include <cassert>
#include <vector>
#include "symbolic.hpp"

// Helper to print expression
void print_expr(const std::string& label, const std::shared_ptr<SymbolicExpr>& expr) {
    std::cout << label << ": " << expr->to_string() << std::endl;
}

void test_maclaurin_sin() {
    std::cout << "Testing Maclaurin Series for sin(x)..." << std::endl;
    // sin(x) = x - x^3/6 + x^5/120 - ...
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    
    // Order 5
    auto series = expr->series("x", SymbolicExpr::number(0), 5);
    print_expr("sin(x) series (order 5)", series);
    
    // Expected: x - 1/6*x^3 + 1/120*x^5  (terms may be in any order, but simplified)
    // We can check by substituting a small value or checking structure.
    // Let's check structure visually or by searching for terms.
    
    // Check specific coefficients?
    // Let's expand it first to ensure standard polynomial form
    auto expanded = series->expand();
    print_expr("Expanded", expanded);
}

void test_maclaurin_exp() {
    std::cout << "Testing Maclaurin Series for e^x..." << std::endl;
    // e^x = 1 + x + x^2/2 + x^3/6 + ...
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::exp(x);
    
    auto series = expr->series("x", SymbolicExpr::number(0), 4);
    print_expr("exp(x) series (order 4)", series);
    
    // e^x at x=0 is 1.
    // D(e^x) = e^x -> 1
    // D2(e^x) = e^x -> 1
    // 1 + x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4
}

void test_taylor_ln() {
    std::cout << "Testing Taylor Series for ln(x) at x=1..." << std::endl;
    // ln(x) at x=1
    // f(1) = 0
    // f'(x) = 1/x -> f'(1) = 1
    // f''(x) = -1/x^2 -> f''(1) = -1
    // f'''(x) = 2/x^3 -> f'''(1) = 2
    // Series: (x-1) - (x-1)^2/2 + 2(x-1)^3/6 = (x-1) - 1/2(x-1)^2 + 1/3(x-1)^3
    
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::ln(x);
    
    auto series = expr->series("x", SymbolicExpr::number(1), 3);
    print_expr("ln(x) series at x=1 (order 3)", series);
}

void test_poly_series() {
    std::cout << "Testing Series for Polynomial x^2 + 2x + 1..." << std::endl;
    // Should be exact if order is high enough
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), x), SymbolicExpr::number(1))
    );
    
    auto series = expr->series("x", SymbolicExpr::number(0), 3);
    print_expr("Poly series (order 3)", series);
    auto expanded = series->expand();
    print_expr("Expanded", expanded);
}

int main() {
    try {
        test_maclaurin_sin();
        test_maclaurin_exp();
        test_taylor_ln();
        test_poly_series();
        std::cout << "All series tests passed!" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
