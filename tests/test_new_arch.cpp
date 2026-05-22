#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <stdexcept>
#include "../include/symbolic_ast.hpp"
#include "../include/visitors/print_visitor.hpp"
#include "../include/visitors/differentiation_visitor.hpp"

void print_expr(const std::string& label, SymbolicNode& node) {
    try {
        PrintVisitor pv;
        node.accept(pv);
        std::cout << label << ": " << pv.get_string() << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error printing " << label << ": " << e.what() << std::endl;
    }
}

int main() {
    try {
        std::cout << "Testing New Architecture..." << std::endl;

        auto x = std::make_shared<VariableNode>("x");

        auto two = std::make_shared<NumberNode>(2.0);

        auto one = std::make_shared<NumberNode>(1.0);

        auto x_sq = std::make_shared<PowerNode>(x, two);

        std::vector<std::shared_ptr<SymbolicNode>> mul_ops;
        mul_ops.push_back(two);
        mul_ops.push_back(x);
        auto two_x = std::make_shared<MultiplyNode>(std::move(mul_ops));

        std::vector<std::shared_ptr<SymbolicNode>> add_ops;
        add_ops.push_back(x_sq);
        add_ops.push_back(two_x);
        add_ops.push_back(one);
        auto poly = std::make_shared<AddNode>(std::move(add_ops));

        print_expr("Original Expression", *poly);

        DifferentiationVisitor diff_v("x");
        poly->accept(diff_v);
        auto derivative = diff_v.get_result();

        if (derivative) {
            print_expr("Derivative (Raw)", *derivative);
        } else {
            std::cout << "Derivative failed (nullptr)" << std::endl;
        }

        std::vector<std::shared_ptr<SymbolicNode>> sin_args;
        sin_args.push_back(x);
        auto sin_x = std::make_shared<FunctionNode>(FunctionNode::FuncType::Sin, std::move(sin_args));

        print_expr("Trig Expression", *sin_x);

        DifferentiationVisitor diff_v2("x");
        std::cout << "Diffing trig... (Variable: x)" << std::endl;
        sin_x->accept(diff_v2);
        std::cout << "Diff done." << std::endl;
        auto sin_deriv = diff_v2.get_result();

        if (sin_deriv) {
            std::cout << "Result acquired. Printing..." << std::endl;
            print_expr("Sin Derivative", *sin_deriv);
            std::cout << "Printed." << std::endl;
        } else {
             std::cout << "Sin deriv null" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "Exception: " << e.what() << std::endl;
    } catch (...) {
        std::cout << "Unknown exception" << std::endl;
    }

    return 0;
}
