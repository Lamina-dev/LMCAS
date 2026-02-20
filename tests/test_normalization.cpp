#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cassert>
#include "../include/symbolic_ast.hpp"
#include "../include/visitors/print_visitor.hpp"
#include "../include/visitors/normalization_visitor.hpp"


void print_expr(const std::string& label, SymbolicNode& node) {
    try {
        PrintVisitor pv;
        node.accept(pv);
        std::cout << label << ": " << pv.get_result() << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error printing " << label << ": " << e.what() << std::endl;
    }
}

std::shared_ptr<SymbolicNode> normalize(std::shared_ptr<SymbolicNode> node) {
    NormalizationVisitor v;
    node->accept(v);
    return v.get_result();
}

int main() {
    try {
        std::cout << "Testing Normalization..." << std::endl;

        
        
        std::vector<std::shared_ptr<SymbolicNode>> ops1;
        ops1.push_back(std::make_shared<NumberNode>(2.0));
        ops1.push_back(std::make_shared<VariableNode>("x"));
        ops1.push_back(std::make_shared<NumberNode>(3.0));
        auto expr1 = std::make_shared<AddNode>(std::move(ops1));

        print_expr("Expr 1 Original", *expr1);
        auto norm1 = normalize(expr1);
        print_expr("Expr 1 Normalized", *norm1);

        
        
        auto x = std::make_shared<VariableNode>("x");
        auto y = std::make_shared<VariableNode>("y");
        auto one = std::make_shared<NumberNode>(1.0);
        
        std::vector<std::shared_ptr<SymbolicNode>> add_ops;
        add_ops.push_back(y);
        add_ops.push_back(one);
        auto y_plus_1 = std::make_shared<AddNode>(std::move(add_ops));
        
        std::vector<std::shared_ptr<SymbolicNode>> mul_ops;
        mul_ops.push_back(x);
        mul_ops.push_back(y_plus_1);
        auto expr2 = std::make_shared<MultiplyNode>(std::move(mul_ops));
        
        print_expr("Expr 2 Original", *expr2);
        auto norm2 = normalize(expr2);
        print_expr("Expr 2 Normalized", *norm2);

        
        std::vector<std::shared_ptr<SymbolicNode>> zero_ops;
        zero_ops.push_back(x);
        zero_ops.push_back(std::make_shared<NumberNode>(0.0));
        auto expr3 = std::make_shared<MultiplyNode>(std::move(zero_ops));
        
        print_expr("Expr 3 Original", *expr3);
        auto norm3 = normalize(expr3);
        print_expr("Expr 3 Normalized", *norm3);
        
        
        auto a = std::make_shared<VariableNode>("a");
        auto b = std::make_shared<VariableNode>("b");
        auto c = std::make_shared<VariableNode>("c");
        auto d = std::make_shared<VariableNode>("d");

        std::vector<std::shared_ptr<SymbolicNode>> ab_ops = {a, b};
        auto a_plus_b = std::make_shared<AddNode>(std::move(ab_ops));

        std::vector<std::shared_ptr<SymbolicNode>> cd_ops = {c, d};
        auto c_plus_d = std::make_shared<AddNode>(std::move(cd_ops));

        std::vector<std::shared_ptr<SymbolicNode>> poly_mul_ops = {a_plus_b, c_plus_d};
        auto expr4 = std::make_shared<MultiplyNode>(std::move(poly_mul_ops));

        print_expr("Expr 4 Original", *expr4);
        
        std::cout << "Normalizing Expr 4..." << std::endl;
        auto norm4 = normalize(expr4);
        if (!norm4) std::cout << "Normalization returned NULL!" << std::endl;
        else std::cout << "Normalization Done. Printing..." << std::endl;
        
        print_expr("Expr 4 Normalized", *norm4);

        
        std::cout << "Testing Inverse Cancellation..." << std::endl;
        auto var_a = std::make_shared<VariableNode>("a");
        auto num_neg2 = std::make_shared<NumberNode>(BigInt(-2));
        
        
        std::vector<std::shared_ptr<SymbolicNode>> inner_ops;
        inner_ops.push_back(var_a);
        inner_ops.push_back(num_neg2);
        auto inner_mul = std::make_shared<MultiplyNode>(inner_ops);
        
        
        auto inv_inner = std::make_shared<PowerNode>(inner_mul, std::make_shared<NumberNode>(BigInt(-1)));
        
        
        std::vector<std::shared_ptr<SymbolicNode>> full_ops;
        full_ops.push_back(var_a);
        full_ops.push_back(inv_inner);
        full_ops.push_back(num_neg2);
        
        auto expr5 = std::make_shared<MultiplyNode>(full_ops);
        print_expr("Expr 5 Original", *expr5);
        auto norm5 = normalize(expr5);
        print_expr("Expr 5 Normalized", *norm5);
        if (norm5->is_one()) std::cout << "PASS: Cancellation successful" << std::endl;
        else std::cout << "FAIL: Cancellation failed" << std::endl;

    } catch (const std::exception& e) {
        std::cout << "Exception: " << e.what() << std::endl;
    } catch (...) {
        std::cout << "Unknown Exception!" << std::endl;
    }

    return 0;
}
