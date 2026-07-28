#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cassert>
#include <cmath>
#include "symbolic.hpp"
#include "integration.hpp"
#include "test_common.hpp"

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
            EXPECT_TRUE(false, name + ": failed to differentiate result");
            return;
        }

        auto diff_simp = diff->simplify();
        std::cout << "Diff(Integral) Simplified: " << diff_simp->to_string() << std::endl;

        auto original_simp = expr.simplify();
        auto diff_verified = test_expr_equivalent(diff_simp, original_simp);

        if (diff_verified) {
            EXPECT_TRUE(diff_verified, name + ": integral differentiates back by simplification match");
            return;
        }

        auto diff_minus_original = test_normalized_delta(diff_simp, original_simp);
        if (diff_minus_original && diff_minus_original->is_zero()) {
            EXPECT_TRUE(diff_minus_original && diff_minus_original->is_zero(),
                        name + ": integral differentiates back by normalized delta");
        } else {
            const std::string delta = diff_minus_original ? diff_minus_original->to_string() : "null";
            EXPECT_TRUE(false, name + ": derivative check not zero: " + delta);
        }
    } catch (const std::exception& e) {
        EXPECT_TRUE(false, name + ": unexpected exception: " + std::string(e.what()));
    }
}

int main() {
    using namespace lamina;

    auto x_var = lamina::detail::make_expression_ptr(*SymbolicExpr::variable("x"));

    auto x2_ptr = SymbolicExpr::power(x_var, lamina::detail::make_expression_ptr(*SymbolicExpr::number(2)));
    run_test("Power Rule x^2", *x2_ptr, "x");

    auto exp_x_ptr = SymbolicExpr::exp(x_var);
    auto x_exp_x_ptr = SymbolicExpr::multiply(x_var, exp_x_ptr);
    run_test("IBP x * exp(x)", *x_exp_x_ptr, "x");

    auto two_x_ptr = SymbolicExpr::multiply(lamina::detail::make_expression_ptr(*SymbolicExpr::number(2)), x_var);
    auto cos_x2_ptr = SymbolicExpr::cos(x2_ptr);
    auto sub_expr_ptr = SymbolicExpr::multiply(cos_x2_ptr, two_x_ptr);
    run_test("Substitution cos(x^2)*2x", *sub_expr_ptr, "x");

    auto denom_ptr = SymbolicExpr::add(x2_ptr, lamina::detail::make_expression_ptr(*SymbolicExpr::number(-1)));
    auto pf_expr_ptr = SymbolicExpr::power(denom_ptr, lamina::detail::make_expression_ptr(*SymbolicExpr::number(-1)));
    run_test("Partial Fraction 1/(x^2-1)", *pf_expr_ptr, "x");

    auto ln_x_ptr = SymbolicExpr::ln(x_var);
    run_test("IBP ln(x)", *ln_x_ptr, "x");

    std::cout << "--------------------------------------------------------" << std::endl;
    return TEST_REPORT();
}
