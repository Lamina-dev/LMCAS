#include <iostream>
#include <vector>
#include <string>
#include "integration.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "test_common.hpp"
#include "test_expression_equivalence.hpp"

using namespace LMCAS;

TEST(CyclicIntegration, DifferentiatesBackToIntegrand) {
    auto x = *(SymbolicExpr::variable("x"));
    auto ex = *(SymbolicExpr::exp(LMCAS::detail::make_expression_ptr(x)));
    auto sinx = *(SymbolicExpr::sin(LMCAS::detail::make_expression_ptr(x)));
    auto expr = *(SymbolicExpr::multiply(LMCAS::detail::make_expression_ptr(ex), LMCAS::detail::make_expression_ptr(sinx)));

    SCOPED_TRACE(expr.to_string());
    Integrator integrator;
    auto result = integrator.integrate(expr, "x");
    ASSERT_TRUE(result) << result.error().message;

    auto diff = result.value().differentiate("x")->simplify();
    ASSERT_NE(diff, nullptr);

    auto diff_check = test_normalized_delta(diff, LMCAS::detail::make_expression_ptr(expr));

    ASSERT_NE(diff_check, nullptr);
    EXPECT_TRUE(diff_check->is_zero()) << diff_check->to_string();
}

TEST(CyclicIntegration, IntegerExponentialProduct) {
    auto x = SymbolicExpr::variable("x");
    auto reciprocal = SymbolicExpr::power(
        SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(-2), x)),
        SymbolicExpr::number(-1));
    auto integrand = SymbolicExpr::multiply(reciprocal,
                                            SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(-1), x)));
    Integrator integrator;
    auto primitive = integrator.integrate(*integrand, "x");
    ASSERT_TRUE((primitive.has_value())) << "exponential collection terminates without a coefficient";
    if (primitive) {
        EXPECT_TRUE(test_proved_equivalent(primitive.value().differentiate("x"), SymbolicExpr::exp(x))) << "integer exponential product primitive differentiates to exp(x)";
    }
}
