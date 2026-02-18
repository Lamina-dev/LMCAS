#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include "../include/symbolic_builder.hpp"
#include "../include/legacy_converter.hpp"
#include "../include/visitors/print_visitor.hpp"

using namespace symbolic::builder;

void test_builder() {
    std::cout << "Testing SymbolicBuilder..." << std::endl;
    auto x = variable("x");
    auto expr = add(power(x, number(2)), add(multiply(number(2), x), number(1)));
    
    PrintVisitor printer;
    expr->accept(printer);
    std::cout << "Builder created: " << printer.get_string() << std::endl;
    // (x^2 + (2*x + 1)) or similar
}

void test_node_to_legacy_and_back() {
    std::cout << "Testing Node -> Legacy -> Node..." << std::endl;
    auto x = variable("x");
    auto original = add(x, number(1));

    // Convert to legacy
    auto legacy = ToLegacyConverter::convert(original);
    assert(legacy != nullptr);
    assert(legacy->type == SymbolicExpr::Type::Add);
    std::cout << "Converted to Legacy type: " << (int)legacy->type << std::endl;

    // Convert back
    auto converted_back = ToNodeConverter::convert(legacy);
    assert(converted_back != nullptr);
    
    PrintVisitor printer;
    converted_back->accept(printer);
    std::cout << "Converted back: " << printer.get_string() << std::endl;
}

void test_legacy_to_node() {
    std::cout << "Testing Legacy -> Node..." << std::endl;
    auto leg_x = SymbolicExpr::variable("y");
    auto leg_num = SymbolicExpr::number(5);
    auto leg_expr = SymbolicExpr::multiply(leg_x, leg_num);
    
    auto node = ToNodeConverter::convert(leg_expr);
    assert(node != nullptr);
    
    PrintVisitor printer;
    node->accept(printer);
    std::cout << "Legacy(y*5) -> Node: " << printer.get_string() << std::endl;
    // Expected: y * 5 or 5 * y depending on implementation
}

void test_complex_conversion() {
    std::cout << "Testing Complex Conversion (Sin, Add, Mul)..." << std::endl;
    // sin(x) + cos(y)
    auto expr = add(sin(variable("x")), cos(variable("y")));
    
    auto legacy = ToLegacyConverter::convert(expr);
    assert(legacy->type == SymbolicExpr::Type::Add);
    
    auto back = ToNodeConverter::convert(legacy);
    assert(back != nullptr);

    PrintVisitor printer;
    back->accept(printer);
    std::cout << "Complex conversion result: " << printer.get_string() << std::endl;
}

int main() {
    try {
        test_builder();
        test_node_to_legacy_and_back();
        test_legacy_to_node();
        test_complex_conversion();
        std::cout << "ALL MIGRATION TESTS PASSED" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Test Failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
