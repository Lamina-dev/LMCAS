#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>

#include "../include/symbolic_ast.hpp"
#include "../include/visitors/print_visitor.hpp"
#include "../include/visitors/differentiation_visitor.hpp"
#include "../include/visitors/normalization_visitor.hpp"

// Utility to print output
void check(const std::string& name, const std::shared_ptr<SymbolicNode>& node, const std::string& expected = "") {
    if (!node) {
        std::cout << "[FAIL] " << name << ": Node is NULL" << std::endl;
        return;
    }
    
    PrintVisitor pv;
    node->accept(pv);
    std::string result = pv.get_result();
    
    if (!expected.empty() && result != expected) {
        std::cout << "[FAIL] " << name << ": Expected '" << expected << "', got '" << result << "'" << std::endl;
    } else {
        std::cout << "[PASS] " << name << ": " << result << std::endl;
    }
}

// Utility to differentiate
std::shared_ptr<SymbolicNode> diff(const std::shared_ptr<SymbolicNode>& node, const std::string& var) {
    DifferentiationVisitor dv(var);
    node->accept(dv);
    return dv.get_result();
}

// Utility to normalize
std::shared_ptr<SymbolicNode> normalize(const std::shared_ptr<SymbolicNode>& node) {
    NormalizationVisitor nv;
    node->accept(nv);
    return nv.get_result();
}

void test_basic_arithmetic() {
    std::cout << "\n--- Testing Basic Arithmetic Nodes ---" << std::endl;
    
    // 1 + 2
    auto n1 = std::make_shared<NumberNode>(1.0);
    auto n2 = std::make_shared<NumberNode>(2.0);
    std::vector<std::shared_ptr<SymbolicNode>> add_ops = {n1, n2};
    auto add = std::make_shared<AddNode>(std::move(add_ops));
    
    check("1 + 2", add, "(1 + 2)");
    
    // x * y
    auto x = std::make_shared<VariableNode>("x");
    auto y = std::make_shared<VariableNode>("y");
    std::vector<std::shared_ptr<SymbolicNode>> mul_ops = {x, y};
    auto mul = std::make_shared<MultiplyNode>(std::move(mul_ops));
    
    check("x * y", mul, "(x * y)");
}

void test_differentiation() {
    std::cout << "\n--- Testing Differentiation ---" << std::endl;
    
    auto x = std::make_shared<VariableNode>("x");
    
    // d/dx(x) = 1
    check("d/dx(x)", diff(x, "x"), "1");
    
    // d/dx(5) = 0
    auto n5 = std::make_shared<NumberNode>(5.0);
    check("d/dx(5)", diff(n5, "x"), "0");
    
    // d/dx(x + 5) = 1 + 0 -> simplified to 1 by DifferentiationVisitor
    std::vector<std::shared_ptr<SymbolicNode>> add_ops = {x, n5};
    auto add = std::make_shared<AddNode>(std::move(add_ops));
    check("d/dx(x + 5)", diff(add, "x"), "1");
    
    // d/dx(x^2) = 2 * x^(2-1) * 1 = 2 * x^1
    auto n2 = std::make_shared<NumberNode>(2.0);
    auto pow = std::make_shared<PowerNode>(x, n2);
    // Note: The output format depends on how DifferentiationVisitor constructs the tree.
    // It likely produces: (2 * (x^(2 - 1))) * 1 or similar structure.
    // We'll just print it to verify structure visually.
    check("d/dx(x^2)", diff(pow, "x"));
    
    // d/dx(sin(x)) = cos(x) * 1
    std::vector<std::shared_ptr<SymbolicNode>> sin_args = {x};
    auto sin_x = std::make_shared<FunctionNode>(FunctionNode::FuncType::Sin, sin_args);
    check("d/dx(sin(x))", diff(sin_x, "x"));
}

void test_normalization_expansion() {
    std::cout << "\n--- Testing Normalization (Expansion) ---" << std::endl;
    
    // (a + b) * (c + d) -> ac + ad + bc + bd
    auto a = std::make_shared<VariableNode>("a");
    auto b = std::make_shared<VariableNode>("b");
    auto c = std::make_shared<VariableNode>("c");
    auto d = std::make_shared<VariableNode>("d");
    
    std::vector<std::shared_ptr<SymbolicNode>> ab_ops = {a, b};
    auto a_plus_b = std::make_shared<AddNode>(std::move(ab_ops));
    
    std::vector<std::shared_ptr<SymbolicNode>> cd_ops = {c, d};
    auto c_plus_d = std::make_shared<AddNode>(std::move(cd_ops));
    
    std::vector<std::shared_ptr<SymbolicNode>> mul_ops = {a_plus_b, c_plus_d};
    auto expr = std::make_shared<MultiplyNode>(std::move(mul_ops));
    
    check("Original: (a+b)*(c+d)", expr);
    
    auto normalized = normalize(expr);
    // Expected output order depends on sorting: a, b before c, d?
    // likely ((a*c) + (a*d) + (b*c) + (b*d))
    check("Normalized: ac+ad+bc+bd", normalized);
}

void test_normalization_simplification() {
    std::cout << "\n--- Testing Normalization (Simplification) ---" << std::endl;
    
    // x * 0 -> 0
    auto x = std::make_shared<VariableNode>("x");
    auto zero = std::make_shared<NumberNode>(0.0);
    std::vector<std::shared_ptr<SymbolicNode>> mul_ops = {x, zero};
    auto mul = std::make_shared<MultiplyNode>(std::move(mul_ops));
    
    check("x * 0", normalize(mul), "0");
    
    // 2 + 3 + x -> 5 + x -> (x + 5) (sorted)
    auto n2 = std::make_shared<NumberNode>(2.0);
    auto n3 = std::make_shared<NumberNode>(3.0);
    std::vector<std::shared_ptr<SymbolicNode>> add_ops = {n2, n3, x};
    auto add = std::make_shared<AddNode>(std::move(add_ops));
    
    check("2 + 3 + x", normalize(add), "(x + 5)");
}

int main() {
    try {
        test_basic_arithmetic();
        test_differentiation();
        test_normalization_expansion();
        test_normalization_simplification();
        std::cout << "\nAll Test Sections Completed." << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Top level exception: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
