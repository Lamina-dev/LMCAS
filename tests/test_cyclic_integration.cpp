#include <iostream>
#include <vector>
#include <string>
#include "integration.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "test_common.hpp"

using namespace lamina;

bool test_cyclic_ibp() {
    std::cout << "Test Case: Cyclic IBP (e^x * sin(x))" << std::endl;

    auto x = *(SymbolicExpr::variable("x"));
    auto ex = *(SymbolicExpr::exp(lamina::detail::make_expression_ptr(x)));
    auto sinx = *(SymbolicExpr::sin(lamina::detail::make_expression_ptr(x)));
    auto expr = *(SymbolicExpr::multiply(lamina::detail::make_expression_ptr(ex), lamina::detail::make_expression_ptr(sinx)));

    std::cout << "Integration of: " << expr.to_string() << std::endl;
    Integrator integrator;
    auto result = integrator.integrate(expr, "x");
    if (!result) {
        std::cout << "[FAIL] Integration failed: " << result.error().message << std::endl;
        return false;
    }
    std::cout << "Result: " << result.value().to_string() << std::endl;

    auto diff = result.value().differentiate("x")->simplify();
    std::cout << "Derivative of result: " << diff->to_string() << std::endl;

    auto diff_check = test_normalized_delta(diff, lamina::detail::make_expression_ptr(expr));

    if (diff_check && diff_check->is_zero()) {
        std::cout << "[PASS]" << std::endl;
        return true;
    } else {
        std::cout << "[FAIL] Derivative check failed: "
                  << (diff_check ? diff_check->to_string() : "null") << std::endl;
        return false;
    }
}

int main() {
    EXPECT_TRUE(test_cyclic_ibp(), "cyclic IBP result differentiates back to the integrand");

    return TEST_REPORT();
}
