#include "expr.hpp"
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>

#include "internal/symbolic_ast.hpp"
#include "internal/visitors/differentiation_visitor.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "test_common.hpp"

using namespace LMCAS;

void check(const std::string &name, const std::shared_ptr<const SymbolicNode> &node,
           const std::string &expected) {
    EXPECT_TRUE(test_proved_equivalent(node ? detail::make_expression_ptr(node) : nullptr, parse_expr(expected).value())) << (name);
}

std::shared_ptr<const SymbolicNode> diff(const std::shared_ptr<const SymbolicNode> &node, const std::string &var) {
    DifferentiationVisitor dv(var);
    node->accept(dv);
    return dv.get_result();
}

std::shared_ptr<const SymbolicNode> normalize(const std::shared_ptr<const SymbolicNode> &node) {
    NormalizationVisitor nv;
    node->accept(nv);
    return nv.get_result();
}

TEST(LmcasSuite, BasicArithmetic) {
    auto n1 = LMCAS::detail::make_node<NumberNode>(1.0);
    auto n2 = LMCAS::detail::make_node<NumberNode>(2.0);
    std::vector<std::shared_ptr<const SymbolicNode>> add_ops = {n1, n2};
    auto add = LMCAS::detail::make_node<AddNode>(std::move(add_ops));

    check("1 + 2", add, "1 + 2");

    auto x = LMCAS::detail::make_node<VariableNode>("x");
    auto y = LMCAS::detail::make_node<VariableNode>("y");
    std::vector<std::shared_ptr<const SymbolicNode>> mul_ops = {x, y};
    auto mul = LMCAS::detail::make_node<MultiplyNode>(std::move(mul_ops));

    check("x * y", mul, "x*y");
}

TEST(LmcasSuite, Differentiation) {
    auto x = LMCAS::detail::make_node<VariableNode>("x");

    check("d/dx(x)", diff(x, "x"), "1");

    auto n5 = LMCAS::detail::make_node<NumberNode>(5.0);
    check("d/dx(5)", diff(n5, "x"), "0");

    std::vector<std::shared_ptr<const SymbolicNode>> add_ops = {x, n5};
    auto add = LMCAS::detail::make_node<AddNode>(std::move(add_ops));
    check("d/dx(x + 5)", diff(add, "x"), "1");

    auto n2 = LMCAS::detail::make_node<NumberNode>(2.0);
    auto pow = LMCAS::detail::make_node<PowerNode>(x, n2);

    check("d/dx(x^2)", diff(pow, "x"), "2*x");

    std::vector<std::shared_ptr<const SymbolicNode>> sin_args = {x};
    auto sin_x = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Sin, sin_args);
    check("d/dx(sin(x))", diff(sin_x, "x"), "cos(x)");
}

TEST(LmcasSuite, NormalizationExpansion) {
    auto a = LMCAS::detail::make_node<VariableNode>("a");
    auto b = LMCAS::detail::make_node<VariableNode>("b");
    auto c = LMCAS::detail::make_node<VariableNode>("c");
    auto d = LMCAS::detail::make_node<VariableNode>("d");

    std::vector<std::shared_ptr<const SymbolicNode>> ab_ops = {a, b};
    auto a_plus_b = LMCAS::detail::make_node<AddNode>(std::move(ab_ops));

    std::vector<std::shared_ptr<const SymbolicNode>> cd_ops = {c, d};
    auto c_plus_d = LMCAS::detail::make_node<AddNode>(std::move(cd_ops));

    std::vector<std::shared_ptr<const SymbolicNode>> mul_ops = {a_plus_b, c_plus_d};
    auto expr = LMCAS::detail::make_node<MultiplyNode>(std::move(mul_ops));

    check("Original: (a+b)*(c+d)", expr, "(a+b)*(c+d)");

    auto normalized = normalize(expr);

    check("Normalization preserves the product value", normalized, "a*c+a*d+b*c+b*d");
}

TEST(LmcasSuite, NormalizationSimplification) {
    auto x = LMCAS::detail::make_node<VariableNode>("x");
    auto zero = LMCAS::detail::make_node<NumberNode>(0.0);
    std::vector<std::shared_ptr<const SymbolicNode>> mul_ops = {x, zero};
    auto mul = LMCAS::detail::make_node<MultiplyNode>(std::move(mul_ops));

    check("x * 0", normalize(mul), "0");

    auto n2 = LMCAS::detail::make_node<NumberNode>(2.0);
    auto n3 = LMCAS::detail::make_node<NumberNode>(3.0);
    std::vector<std::shared_ptr<const SymbolicNode>> add_ops = {n2, n3, x};
    auto add = LMCAS::detail::make_node<AddNode>(std::move(add_ops));

    check("2 + 3 + x", normalize(add), "x + 5");
}
