#include <iostream>
#include <cassert>
#include <vector>
#include "symbolic.hpp"


void print_expr(const std::string& label, const std::shared_ptr<SymbolicExpr>& expr) {
    std::cout << label << ": " << expr->to_string() << std::endl;
}

void test_maclaurin_sin() {
    std::cout << "Testing Maclaurin Series for sin(x)..." << std::endl;
    
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(x);
    
    
    auto series = expr->series("x", SymbolicExpr::number(0), 5);
    print_expr("sin(x) series (order 5)", series);
    
    
    
    
    
    
    
    auto expanded = series->expand();
    print_expr("Expanded", expanded);
}

void test_maclaurin_exp() {
    std::cout << "Testing Maclaurin Series for e^x..." << std::endl;
    
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::exp(x);
    
    auto series = expr->series("x", SymbolicExpr::number(0), 4);
    print_expr("exp(x) series (order 4)", series);
    
    
    
    
    
}

void test_taylor_ln() {
    std::cout << "Testing Taylor Series for ln(x) at x=1..." << std::endl;
    
    
    
    
    
    
    
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::ln(x);
    
    auto series = expr->series("x", SymbolicExpr::number(1), 3);
    print_expr("ln(x) series at x=1 (order 3)", series);
}

void test_poly_series() {
    std::cout << "Testing Series for Polynomial x^2 + 2x + 1..." << std::endl;
    
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
