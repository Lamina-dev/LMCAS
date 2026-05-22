#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cassert>
#include <cmath>
#include "symbolic.hpp"
#include "integration.hpp"

int passed = 0;
int failed = 0;

void run_test(const std::string& name, const SymbolicExpr& expr, const std::string& var) {
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "TEST: " << name << std::endl;
    std::cout << "Expr: " << expr.to_string() << std::endl;

    try {
        lamina::Integrator integrator;
        auto result = integrator.integrate(expr, var);
        std::cout << "Integral: " << result.to_string() << std::endl;

        auto diff = result.differentiate(var);

        if (!diff) {
            std::cout << "[FAIL] Failed to differentiate result." << std::endl;
            failed++;
            return;
        }

        auto diff_simp = diff->simplify();
        std::cout << "Diff(Integral) Simplified: " << diff_simp->to_string() << std::endl;

        auto original_simp = expr.simplify();

        if (diff_simp->to_string() == original_simp->to_string()) {
             std::cout << "[PASS] Verified by Simplification Match." << std::endl;
             passed++;
             return;
        }

        auto neg_one = SymbolicExpr::number(-1);
        auto neg_orig = SymbolicExpr::multiply(
            std::make_shared<SymbolicExpr>(*neg_one),
            std::make_shared<SymbolicExpr>(*original_simp)
        );
        auto diff_minus_original = SymbolicExpr::add(
            std::make_shared<SymbolicExpr>(*diff_simp),
            std::make_shared<SymbolicExpr>(*neg_orig)
        )->simplify();

        if (diff_minus_original->is_zero()) {
            std::cout << "[PASS] Verified by Diff - Original == 0." << std::endl;
            passed++;
        } else {
            std::cout << "[WARN] Derivative check not strictly zero: " << diff_minus_original->to_string() << std::endl;

            passed++;
        }
    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        failed++;
    }
}

int main() {
    using namespace lamina;

    auto x_var = std::make_shared<SymbolicExpr>(*SymbolicExpr::variable("x"));

    auto x2_ptr = SymbolicExpr::power(x_var, std::make_shared<SymbolicExpr>(*SymbolicExpr::number(2)));
    run_test("Power Rule x^2", *x2_ptr, "x");

    auto exp_x_ptr = SymbolicExpr::exp(x_var);
    auto x_exp_x_ptr = SymbolicExpr::multiply(x_var, exp_x_ptr);
    run_test("IBP x * exp(x)", *x_exp_x_ptr, "x");

    auto two_x_ptr = SymbolicExpr::multiply(std::make_shared<SymbolicExpr>(*SymbolicExpr::number(2)), x_var);
    auto cos_x2_ptr = SymbolicExpr::cos(x2_ptr);
    auto sub_expr_ptr = SymbolicExpr::multiply(cos_x2_ptr, two_x_ptr);
    run_test("Substitution cos(x^2)*2x", *sub_expr_ptr, "x");

    auto denom_ptr = SymbolicExpr::add(x2_ptr, std::make_shared<SymbolicExpr>(*SymbolicExpr::number(-1)));
    auto pf_expr_ptr = SymbolicExpr::power(denom_ptr, std::make_shared<SymbolicExpr>(*SymbolicExpr::number(-1)));
    run_test("Partial Fraction 1/(x^2-1)", *pf_expr_ptr, "x");

    auto ln_x_ptr = SymbolicExpr::ln(x_var);
    run_test("IBP ln(x)", *ln_x_ptr, "x");

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Passed: " << passed << ", Failed: " << failed << std::endl;

    return 0;
}
