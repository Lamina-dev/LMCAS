#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <string>

#include "../include/symbolic_ast.hpp"
#include "../include/visitors/print_visitor.hpp"
#include "../include/visitors/normalization_visitor.hpp"
#include "test_common.hpp"

int main() {
    try {
        std::cout << "Starting Debug Expansion Test..." << std::endl;

        auto a = lamina::detail::make_node<VariableNode>("a");
        auto b = lamina::detail::make_node<VariableNode>("b");
        std::vector<std::shared_ptr<const SymbolicNode>> ab_ops = {a, b};
        auto a_plus_b = lamina::detail::make_node<AddNode>(std::move(ab_ops));

        auto c = lamina::detail::make_node<VariableNode>("c");
        auto d = lamina::detail::make_node<VariableNode>("d");
        std::vector<std::shared_ptr<const SymbolicNode>> cd_ops = {c, d};
        auto c_plus_d = lamina::detail::make_node<AddNode>(std::move(cd_ops));

        std::cout << "Created generic operands." << std::endl;

        NormalizationVisitor v;

        std::cout << "Test 1: (a+b)*c" << std::endl;
        auto res1 = v.expand_product(a_plus_b, c);

        PrintVisitor pv;
        if(res1) {
            res1->accept(pv);
            std::cout << "Result: " << pv.get_result() << std::endl;
        } else {
            std::cout << "Result: NULL" << std::endl;
        }
        EXPECT_TRUE(res1 != nullptr, "expand_product((a+b), c) returns a result");

        std::cout << "Test 2: (a+b)*(c+d)" << std::endl;

        auto res2 = v.expand_product(a_plus_b, c_plus_d);
        std::cout << "Expansion returned." << std::endl;

        PrintVisitor pv2;
        if(res2) {
            std::cout << "Printing result..." << std::endl;
            if (auto add = dynamic_cast<const AddNode*>(res2.get())) {
                std::cout << "Result is AddNode with " << add->operands().size() << " operands." << std::endl;
            }
            res2->accept(pv2);
            std::cout << "Result: " << pv2.get_result() << std::endl;
        } else {
            std::cout << "Result: NULL" << std::endl;
        }
        EXPECT_TRUE(res2 != nullptr, "expand_product((a+b), (c+d)) returns a result");

    } catch (const std::exception& e) {
        EXPECT_TRUE(false, std::string("unexpected exception: ") + e.what());
    }

    return TEST_REPORT();
}
