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
    std::cout << "case: " << name << std::endl;
    std::cout << "Expr: " << expr.to_string() << std::endl;

    LMCAS::Integrator integrator;
    auto result = integrator.integrate(expr, var);
    if (!result) {
        EXPECT_TRUE(false, name + ": integration failed: " + result.error().message);
        return;
    }
    std::cout << "Integral: " << result.value().to_string() << std::endl;

    auto diff = result.value().differentiate(var);

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

}

class WrongIntegrationStrategy final : public LMCAS::IntegrationStrategy {
public:
    LMCAS::Result<std::shared_ptr<SymbolicExpr>> try_integrate_raw(
        const SymbolicExpr&,
        const std::string& variable,
        LMCAS::Integrator&,
        LMCAS::ComputationContext&,
        int) override {
        return SymbolicExpr::variable(variable);
    }

    std::string name() const override { return "DeliberatelyWrong"; }
};

int main() {
    using namespace LMCAS;

    auto x_var = LMCAS::detail::make_expression_ptr(*SymbolicExpr::variable("x"));
    TEST_CASE("Generated integration candidates require exact residual proof");
    Integrator gated_integrator;
    auto added = gated_integrator.add_strategy(
        std::make_unique<WrongIntegrationStrategy>(), 0);
    EXPECT_TRUE(added.has_value(), "wrong strategy is installed for the gate test");
    ComputationContext gate_context;
    auto gated = gated_integrator.integrate_checked(
        *x_var, "x", gate_context);
    EXPECT_TRUE(gated && gated.value().to_string() != "x",
                "wrong integration candidate is rejected");


    auto x2_ptr = SymbolicExpr::power(x_var, LMCAS::detail::make_expression_ptr(*SymbolicExpr::number(2)));
    run_test("Power Rule x^2", *x2_ptr, "x");

    auto exp_x_ptr = SymbolicExpr::exp(x_var);
    auto x_exp_x_ptr = SymbolicExpr::multiply(x_var, exp_x_ptr);
    run_test("IBP x * exp(x)", *x_exp_x_ptr, "x");

    auto two_x_ptr = SymbolicExpr::multiply(LMCAS::detail::make_expression_ptr(*SymbolicExpr::number(2)), x_var);
    auto cos_x2_ptr = SymbolicExpr::cos(x2_ptr);
    auto sub_expr_ptr = SymbolicExpr::multiply(cos_x2_ptr, two_x_ptr);
    run_test("Substitution cos(x^2)*2x", *sub_expr_ptr, "x");

    auto denom_ptr = SymbolicExpr::add(x2_ptr, LMCAS::detail::make_expression_ptr(*SymbolicExpr::number(-1)));
    auto pf_expr_ptr = SymbolicExpr::power(denom_ptr, LMCAS::detail::make_expression_ptr(*SymbolicExpr::number(-1)));
    run_test("Partial Fraction 1/(x^2-1)", *pf_expr_ptr, "x");

    auto ln_x_ptr = SymbolicExpr::ln(x_var);
    run_test("IBP ln(x)", *ln_x_ptr, "x");

    std::cout << "--------------------------------------------------------" << std::endl;
    return TEST_REPORT();
}
