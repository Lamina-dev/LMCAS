#include "expr.hpp"
#include <iostream>
#include <vector>
#include <string>
#include "integration.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "test_common.hpp"
#include "test_expression_equivalence.hpp"
#include "residual_verification.hpp"
#include "internal/expression_analysis.hpp"
#include "assumption_context.hpp"

using namespace LMCAS;

TEST(CyclicIntegration, DifferentiatesBackToIntegrand) {
    auto x = *(SymbolicExpr::variable("x"));
    auto ex = *(SymbolicExpr::exp(LMCAS::detail::make_expression_ptr(x)));
    auto sinx = *(SymbolicExpr::sin(LMCAS::detail::make_expression_ptr(x)));
    auto expr = *(SymbolicExpr::multiply(LMCAS::detail::make_expression_ptr(ex), LMCAS::detail::make_expression_ptr(sinx)));

    SCOPED_TRACE(expr.to_string());
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(expr, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    ASSERT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));

    auto diff = LMCAS::simplify(LMCAS::differentiate(std::make_shared<SymbolicExpr>(result.value()), "x").value()).value();
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
    ComputationContext context;
    auto primitive = integrator.integrate_checked(*integrand, "x", context);
    ASSERT_TRUE((primitive.has_value())) << "exponential collection terminates without a coefficient";
    ASSERT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(primitive.value())));
    if (primitive) {
        EXPECT_TRUE(test_proved_equivalent(LMCAS::differentiate(std::make_shared<SymbolicExpr>(primitive.value()), "x").value(), SymbolicExpr::exp(x))) << "integer exponential product primitive differentiates to exp(x)";
    }
}

TEST(CyclicIntegration, PrivatePlaceholdersAreEliminated) {
    auto x = SymbolicExpr::variable("x");
    auto occupied = SymbolicExpr::variable("INT_CYCLE_0");
    auto exponential = SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(2), x));
    auto sin3x = SymbolicExpr::sin(SymbolicExpr::multiply(SymbolicExpr::number(3), x));
    auto cos3x = SymbolicExpr::cos(SymbolicExpr::multiply(SymbolicExpr::number(3), x));
    auto integrand = SymbolicExpr::add(SymbolicExpr::multiply(exponential, sin3x), occupied);
    auto expected = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 13)),
            SymbolicExpr::multiply(exponential,
                SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), sin3x),
                    SymbolicExpr::multiply(SymbolicExpr::number(-3), cos3x)))),
        SymbolicExpr::multiply(occupied, x));
    Integrator integrator;
    ComputationContext context;
    auto assumptions = std::make_shared<AssumptionContext>();
    ASSERT_TRUE(assumptions->assume_sign("INT_CYCLE_1", Sign::Zero));
    assumptions->push();
    ASSERT_TRUE(assumptions->assume_sign("INT_CYCLE_2", Sign::Positive));
    ASSERT_TRUE(context.set_assumptions(assumptions));
    auto result = integrator.integrate_checked(*integrand, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    ASSERT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));
    auto names = detail::all_variable_names(detail::node(result.value()));
    EXPECT_EQ(names, (std::set<std::string>{"x", "INT_CYCLE_0"}));
    EXPECT_TRUE(test_proved_equivalent(detail::make_expression_ptr(result.value()), expected));
    auto derivative = LMCAS::differentiate(detail::make_expression_ptr(result.value()), "x", context);
    ASSERT_TRUE(derivative);
    auto proof = check_equivalent(derivative.value(), integrand, context);
    ASSERT_TRUE(proof);
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(proof.value()));
}
