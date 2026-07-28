#include "integration.hpp"
#include "symbolic_ast.hpp"
#include "test_common.hpp"
#include <iostream>
#include <cmath>
#include <string>

using namespace lamina;

std::shared_ptr<SymbolicExpr> MakeSymbolicExprPtr(const SymbolicExpr& e) {
    return lamina::detail::make_expression_ptr(e);
}

double evaluate_symbolic(const SymbolicExpr& expr) {
    if (!lamina::detail::node(expr)) return 0.0;
    if (auto n = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(expr))) {
        if (std::holds_alternative<double>(n->value())) return std::get<double>(n->value());
        if (std::holds_alternative<BigInt>(n->value())) return std::get<BigInt>(n->value()).to_double();
        if (std::holds_alternative<Rational>(n->value())) return std::get<Rational>(n->value()).to_double();
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(expr))) {
        double sum = 0;
        for (auto& op : add->operands()) sum += evaluate_symbolic(lamina::detail::expression_from_node(op));
        return sum;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
        double prod = 1;
        for (auto& op : mul->operands()) prod *= evaluate_symbolic(lamina::detail::expression_from_node(op));
        return prod;
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr))) {
        return std::pow(evaluate_symbolic(lamina::detail::expression_from_node(pow->base())), evaluate_symbolic(lamina::detail::expression_from_node(pow->exponent())));
    }
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
        double arg = evaluate_symbolic(lamina::detail::expression_from_node(func->arguments()[0]));
        if (func->type() == FunctionNode::FuncType::Sin) return std::sin(arg);
        if (func->type() == FunctionNode::FuncType::Cos) return std::cos(arg);
        if (func->type() == FunctionNode::FuncType::Tan) return std::tan(arg);
        if (func->type() == FunctionNode::FuncType::Exp) return std::exp(arg);
        if (func->type() == FunctionNode::FuncType::Ln) return std::log(arg);
    }
    return 0.0;
}

void test_polynomial() {
    std::cout << "Test Case 1: Polynomial x^2 from 0 to 1" << std::endl;
    Integrator integrator;
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::power(MakeSymbolicExprPtr(*x), SymbolicExpr::number(2));

    auto lower = SymbolicExpr::number(0);
    auto upper = SymbolicExpr::number(1);

    auto res = integrator.integrate_def(*expr, "x", *lower, *upper);
    std::cout << "Definite Integral result: " << res.to_string() << std::endl;

    EXPECT_NEAR(evaluate_symbolic(res), 1.0/3.0, 1e-9,
                "integral of x^2 from 0 to 1 is 1/3");

    std::cout << "[PASS]" << std::endl;
}

void test_trig() {
    std::cout << "Test Case 2: Sin(x) from 0 to Pi" << std::endl;
    Integrator integrator;
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::sin(MakeSymbolicExprPtr(*x));

    auto lower = SymbolicExpr::number(0);
    auto upper = SymbolicExpr::number(3.14159265358979323846);

    auto res = integrator.integrate_def(*expr, "x", *lower, *upper);
    std::cout << "Definite Integral result: " << res.to_string() << std::endl;

    EXPECT_NEAR(evaluate_symbolic(res), 2.0, 1e-5,
                "integral of sin(x) from 0 to pi is 2");

    std::cout << "[PASS]" << std::endl;
}

void test_symbolic_limits() {
    std::cout << "Test Case 3: x from a to b" << std::endl;
    Integrator integrator;
    auto x = SymbolicExpr::variable("x");

    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");

    auto res = integrator.integrate_def(*x, "x", *a, *b);
    std::cout << "Definite Integral result: " << res.to_string() << std::endl;

    auto check = res.substitute("a", SymbolicExpr::number(2))->substitute("b", SymbolicExpr::number(4))->simplify();

    EXPECT_NEAR(evaluate_symbolic(*check), 6.0, 1e-9,
                "symbolic integral of x from a to b evaluates correctly at a=2,b=4");

    std::cout << "[PASS]" << std::endl;
}

void test_improper_split_ignores_nonrepresentable_exact_bounds() {
    std::cout << "Test Case 4: 1/x with nonrepresentable exact bounds does not throw" << std::endl;
    Integrator integrator;
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::power(MakeSymbolicExprPtr(*x), SymbolicExpr::number(-1));

    std::string huge_digits = "1" + std::string(400, '0');
    auto lower = SymbolicExpr::number(BigInt("-" + huge_digits));
    auto upper = SymbolicExpr::number(BigInt(huge_digits));

    bool threw = false;
    std::shared_ptr<SymbolicExpr> res;
    try {
        res = lamina::detail::make_expression_ptr(
            integrator.integrate_def(*expr, "x", *lower, *upper));
    } catch (...) {
        threw = true;
    }

    EXPECT_TRUE(!threw && res,
                "integrate_def preserves symbolic result for nonrepresentable exact bounds");

    if (res) {
        std::cout << "Definite Integral result: " << res->to_string() << std::endl;
    }
    std::cout << "[PASS]" << std::endl;
}

int main() {
    test_polynomial();
    test_trig();
    test_symbolic_limits();
    test_improper_split_ignores_nonrepresentable_exact_bounds();
    return TEST_REPORT();
}
