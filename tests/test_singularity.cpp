#include "integration.hpp"
#include "symbolic_ast.hpp"
#include "test_common.hpp"
#include <iostream>
#include <cmath>

using namespace lamina;

std::shared_ptr<SymbolicExpr> MakeSymbolicExprPtr(const SymbolicExpr& e) {
    return lamina::detail::make_expression_ptr(e);
}

void test_one_sided_limit() {
    std::cout << "Test Case 1: One-sided limit 1/x at 0" << std::endl;

    auto x = SymbolicExpr::variable("x");

    auto expr = SymbolicExpr::power(MakeSymbolicExprPtr(*x), SymbolicExpr::number(-1));
    auto zero = SymbolicExpr::number(0);

    auto limit_right = expr->limit("x", MakeSymbolicExprPtr(*zero), "+");
    std::cout << "Limit x->0+ 1/x = " << limit_right->to_string() << std::endl;

    auto limit_left = expr->limit("x", MakeSymbolicExprPtr(*zero), "-");
    std::cout << "Limit x->0- 1/x = " << limit_left->to_string() << std::endl;

    bool has_inf = false;
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(limit_right))) {
        if(func->type() == FunctionNode::FuncType::Infinity) has_inf = true;
    }
    EXPECT_TRUE(has_inf, "right-hand limit of 1/x at 0 is infinity");
}

void test_improper_integral_singularity() {
    std::cout << "Test Case 2: Improper Integral 1/x from -1 to 1" << std::endl;
    Integrator integrator;
    auto x = SymbolicExpr::variable("x");

    auto expr = SymbolicExpr::power(MakeSymbolicExprPtr(*x), SymbolicExpr::number(-1));

    auto lower = SymbolicExpr::number(-1);
    auto upper = SymbolicExpr::number(1);

    auto res = integrator.integrate_def(*expr, "x", *lower, *upper);
    std::cout << "Result: " << res.to_string() << std::endl;

    EXPECT_TRUE(lamina::detail::node(res) != nullptr, "improper integral returns a non-null expression");
}

int main() {
    test_one_sided_limit();
    test_improper_integral_singularity();
    return TEST_REPORT();
}
