#include "integration.hpp"
#include "symbolic_ast.hpp"
#include "test_common.hpp"
#include <iostream>
#include <cmath>

using namespace LMCAS;

static std::shared_ptr<SymbolicExpr> MakeSymbolicExprPtr(const SymbolicExpr& e) {
    return LMCAS::detail::make_expression_ptr(e);
}

void test_one_sided_limit() {
    std::cout << "Test Case 1: One-sided limit 1/x at 0" << std::endl;

    auto x = SymbolicExpr::variable("x");

    auto expr = SymbolicExpr::power(MakeSymbolicExprPtr(*x), SymbolicExpr::number(-1));
    auto zero = SymbolicExpr::number(0);

    auto limit_right = LMCAS::limit_expression_checked(expr, "x", MakeSymbolicExprPtr(*zero), LimitDirection::FromAbove).value();
    std::cout << "Limit x->0+ 1/x = " << limit_right->to_string() << std::endl;

    auto limit_left = LMCAS::limit_expression_checked(expr, "x", MakeSymbolicExprPtr(*zero), LimitDirection::FromBelow).value();
    std::cout << "Limit x->0- 1/x = " << limit_left->to_string() << std::endl;

    bool has_inf = false;
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(limit_right))) {
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
    EXPECT_TRUE(res.has_value(), "improper integration returns a Result");
    if (!res) return;
    std::cout << "Result: " << res.value().to_string() << std::endl;

    EXPECT_TRUE(LMCAS::detail::node(res.value()) != nullptr,
                "improper integral returns a non-null expression");
}

int main() {
    test_one_sided_limit();
    test_improper_integral_singularity();
    return TEST_REPORT();
}
