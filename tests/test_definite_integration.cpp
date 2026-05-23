#include "integration.hpp"
#include "symbolic_ast.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

#define ASSERT_EQ(a, b) \
    if ((a) != (b)) { \
        std::cerr << "Assertion failed: " << (a) << " != " << (b) << std::endl; \
        std::exit(1); \
    }

#define ASSERT_NEAR(a, b, eps) \
    if (std::abs((a) - (b)) > (eps)) { \
        std::cerr << "Assertion failed: " << (a) << " != " << (b) << std::endl; \
        std::exit(1); \
    }

using namespace lamina;

std::shared_ptr<SymbolicExpr> MakeSymbolicExprPtr(const SymbolicExpr& e) {
    return std::make_shared<SymbolicExpr>(e);
}

double evaluate_symbolic(const SymbolicExpr& expr) {
    if (!expr.root) return 0.0;
    if (auto n = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (std::holds_alternative<double>(n->value)) return std::get<double>(n->value);
        if (std::holds_alternative<BigInt>(n->value)) return std::get<BigInt>(n->value).to_double();
        if (std::holds_alternative<Rational>(n->value)) return std::get<Rational>(n->value).to_double();
    }
    if (auto add = std::dynamic_pointer_cast<AddNode>(expr.root)) {
        double sum = 0;
        for (auto& op : add->operands) sum += evaluate_symbolic(SymbolicExpr(op));
        return sum;
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        double prod = 1;
        for (auto& op : mul->operands) prod *= evaluate_symbolic(SymbolicExpr(op));
        return prod;
    }
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        return std::pow(evaluate_symbolic(SymbolicExpr(pow->base)), evaluate_symbolic(SymbolicExpr(pow->exponent)));
    }
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(expr.root)) {
        double arg = evaluate_symbolic(SymbolicExpr(func->arguments[0]));
        if (func->type == FunctionNode::FuncType::Sin) return std::sin(arg);
        if (func->type == FunctionNode::FuncType::Cos) return std::cos(arg);
        if (func->type == FunctionNode::FuncType::Tan) return std::tan(arg);
        if (func->type == FunctionNode::FuncType::Exp) return std::exp(arg);
        if (func->type == FunctionNode::FuncType::Ln) return std::log(arg);
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

    ASSERT_NEAR(evaluate_symbolic(res), 1.0/3.0, 1e-9);

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

    ASSERT_NEAR(evaluate_symbolic(res), 2.0, 1e-5);

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

    ASSERT_NEAR(evaluate_symbolic(*check), 6.0, 1e-9);

    std::cout << "[PASS]" << std::endl;
}

int main() {
    test_polynomial();
    test_trig();
    test_symbolic_limits();
    return 0;
}
