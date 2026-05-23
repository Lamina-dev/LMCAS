#include <iostream>
#include <vector>
#include <string>
#include "integration.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"

using namespace lamina;

bool test_cyclic_ibp() {
    std::cout << "Test Case: Cyclic IBP (e^x * sin(x))" << std::endl;

    auto x = *(SymbolicExpr::variable("x"));
    auto ex = *(SymbolicExpr::exp(std::make_shared<SymbolicExpr>(x)));
    auto sinx = *(SymbolicExpr::sin(std::make_shared<SymbolicExpr>(x)));
    auto expr = *(SymbolicExpr::multiply(std::make_shared<SymbolicExpr>(ex), std::make_shared<SymbolicExpr>(sinx)));

    std::cout << "Integration of: " << expr.to_string() << std::endl;
    Integrator integrator;
    auto result = integrator.integrate(expr, "x");
    std::cout << "Result: " << result.to_string() << std::endl;

    auto diff = result.differentiate("x")->simplify();
    std::cout << "Derivative of result: " << diff->to_string() << std::endl;

    auto diff_check = SymbolicExpr::add(diff, SymbolicExpr::multiply(std::make_shared<SymbolicExpr>(expr), SymbolicExpr::number(-1)))->simplify();

    if (diff_check->is_zero()) {
        std::cout << "[PASS]" << std::endl;
        return true;
    } else {
        std::cout << "[FAIL] Derivative check failed." << std::endl;
        return false;
    }
}

int main() {
    bool all_passed = true;
    if (!test_cyclic_ibp()) all_passed = false;

    return all_passed ? 0 : 1;
}
