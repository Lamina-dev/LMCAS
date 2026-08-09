#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cassert>
#include "../include/symbolic_ast.hpp"
#include "../include/visitors/print_visitor.hpp"
#include "../include/visitors/normalization_visitor.hpp"
#include "test_common.hpp"

void print_expr(const std::string& label, const SymbolicNode& node) {
    try {
        PrintVisitor pv;
        node.accept(pv);
        std::cout << label << ": " << pv.get_result() << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error printing " << label << ": " << e.what() << std::endl;
    }
}

std::shared_ptr<const SymbolicNode> normalize(std::shared_ptr<const SymbolicNode> node) {
    NormalizationVisitor v;
    node->accept(v);
    return v.get_result();
}

int main() {
    try {
        std::cout << "Testing Normalization..." << std::endl;

        std::vector<std::shared_ptr<const SymbolicNode>> ops1;
        ops1.push_back(lamina::detail::make_node<NumberNode>(2.0));
        ops1.push_back(lamina::detail::make_node<VariableNode>("x"));
        ops1.push_back(lamina::detail::make_node<NumberNode>(3.0));
        auto expr1 = lamina::detail::make_node<AddNode>(std::move(ops1));

        print_expr("Expr 1 Original", *expr1);
        auto norm1 = normalize(expr1);
        EXPECT_TRUE(norm1 != nullptr, "numeric/add normalization returns a result");
        if (norm1) {
            print_expr("Expr 1 Normalized", *norm1);
        }

        auto x = lamina::detail::make_node<VariableNode>("x");
        auto y = lamina::detail::make_node<VariableNode>("y");
        auto one = lamina::detail::make_node<NumberNode>(1.0);

        std::vector<std::shared_ptr<const SymbolicNode>> add_ops;
        add_ops.push_back(y);
        add_ops.push_back(one);
        auto y_plus_1 = lamina::detail::make_node<AddNode>(std::move(add_ops));

        std::vector<std::shared_ptr<const SymbolicNode>> mul_ops;
        mul_ops.push_back(x);
        mul_ops.push_back(y_plus_1);
        auto expr2 = lamina::detail::make_node<MultiplyNode>(std::move(mul_ops));

        print_expr("Expr 2 Original", *expr2);
        auto norm2 = normalize(expr2);
        EXPECT_TRUE(norm2 != nullptr, "product normalization returns a result");
        if (norm2) {
            print_expr("Expr 2 Normalized", *norm2);
        }

        std::vector<std::shared_ptr<const SymbolicNode>> zero_ops;
        zero_ops.push_back(x);
        zero_ops.push_back(lamina::detail::make_node<NumberNode>(0.0));
        auto expr3 = lamina::detail::make_node<MultiplyNode>(std::move(zero_ops));

        print_expr("Expr 3 Original", *expr3);
        auto norm3 = normalize(expr3);
        EXPECT_TRUE(norm3 != nullptr && norm3->is_zero(),
                    "multiplication by zero normalizes to zero");
        if (norm3) {
            print_expr("Expr 3 Normalized", *norm3);
        }

        auto a = lamina::detail::make_node<VariableNode>("a");
        auto b = lamina::detail::make_node<VariableNode>("b");
        auto c = lamina::detail::make_node<VariableNode>("c");
        auto d = lamina::detail::make_node<VariableNode>("d");

        std::vector<std::shared_ptr<const SymbolicNode>> ab_ops = {a, b};
        auto a_plus_b = lamina::detail::make_node<AddNode>(std::move(ab_ops));

        std::vector<std::shared_ptr<const SymbolicNode>> cd_ops = {c, d};
        auto c_plus_d = lamina::detail::make_node<AddNode>(std::move(cd_ops));

        std::vector<std::shared_ptr<const SymbolicNode>> poly_mul_ops = {a_plus_b, c_plus_d};
        auto expr4 = lamina::detail::make_node<MultiplyNode>(std::move(poly_mul_ops));

        print_expr("Expr 4 Original", *expr4);

        std::cout << "Normalizing Expr 4..." << std::endl;
        auto norm4 = normalize(expr4);
        EXPECT_TRUE(norm4 != nullptr, "polynomial product normalization returns a result");
        if (norm4) std::cout << "Normalization Done. Printing..." << std::endl;

        if (norm4) print_expr("Expr 4 Normalized", *norm4);

        std::cout << "Testing Inverse Cancellation..." << std::endl;
        auto var_a = lamina::detail::make_node<VariableNode>("a");
        auto num_neg2 = lamina::detail::make_node<NumberNode>(BigInt(-2));

        std::vector<std::shared_ptr<const SymbolicNode>> inner_ops;
        inner_ops.push_back(var_a);
        inner_ops.push_back(num_neg2);
        auto inner_mul = lamina::detail::make_node<MultiplyNode>(inner_ops);

        auto inv_inner = lamina::detail::make_node<PowerNode>(inner_mul, lamina::detail::make_node<NumberNode>(BigInt(-1)));

        std::vector<std::shared_ptr<const SymbolicNode>> full_ops;
        full_ops.push_back(var_a);
        full_ops.push_back(inv_inner);
        full_ops.push_back(num_neg2);

        auto expr5 = lamina::detail::make_node<MultiplyNode>(full_ops);
        print_expr("Expr 5 Original", *expr5);
        auto norm5 = normalize(expr5);
        EXPECT_TRUE(norm5 != nullptr, "inverse normalization returns a result");
        if (norm5) {
            print_expr("Expr 5 Normalized", *norm5);
        }

    } catch (const std::exception& e) {
        EXPECT_TRUE(false, std::string("unexpected exception: ") + e.what());
    } catch (...) {
        EXPECT_TRUE(false, "unexpected unknown exception");
    }

    return TEST_REPORT();
}
