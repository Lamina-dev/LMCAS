#include "test_common.hpp"
#include "transform_engine.hpp"
#include "assumption_context.hpp"
#include "residual_verification.hpp"

using namespace LMCAS;

static std::shared_ptr<SymbolicExpr> decay_input(
    const std::shared_ptr<SymbolicExpr> &rate) {
    auto absolute = detail::make_expression_ptr(detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<const SymbolicNode>>{SymbolicFactory::create_variable("t")}));
    auto exponent = SymbolicExpr::multiply(SymbolicExpr::number(-2),
                                           SymbolicExpr::multiply(rate, absolute));
    return SymbolicExpr::multiply(SymbolicExpr::number(3), SymbolicExpr::exp(exponent));
}

TEST(TransformExactDecay, ExactLaplaceFactorial) {
    ComputationContext context;
    auto input = SymbolicExpr::power(SymbolicExpr::variable("t"), SymbolicExpr::number(21));
    auto result = laplace_transform_checked(input, "t", "s", context);
    ASSERT_TRUE((result.has_value())) << "exact polynomial Laplace succeeds";
    if (!result)
        return;
    auto value = result.value().value.expression->substitute("s", SymbolicExpr::number(1))->simplify();
    auto number = std::dynamic_pointer_cast<const NumberNode>(detail::node(value));
    EXPECT_TRUE((number != nullptr)) << "Laplace value is an exact number";
    if (!number)
        return;
    const auto *integer = std::get_if<BigInt>(&number->value());
    EXPECT_TRUE((integer != nullptr)) << "factorial value is BigInt, not rounded floating point";
    if (integer) {
        EXPECT_TRUE((*integer == BigInt("51090942171709440000"))) << "21! is exact";
    }
}

TEST(TransformExactDecay, FactorialContextBudget) {
    ResourceLimits limits;
    limits.max_integer_bits = 32;
    ComputationContext context(limits);
    auto input = SymbolicExpr::power(SymbolicExpr::variable("t"), SymbolicExpr::number(21));
    auto result = laplace_transform_checked(input, "t", "s", context);
    EXPECT_TRUE((!result && result.error().code == CasErrc::ResourceLimit)) << "factorial cannot escape the transform context's integer budget";
    if (!result) {
        EXPECT_TRUE((result.error().operation == "BigInt::factorial_checked")) << "checked factorial error metadata survives the transform";
    }
}

TEST(TransformExactDecay, LargerExactFactorial) {
    ComputationContext context;
    auto result = BigInt::factorial_checked(BigInt(30), context);
    ASSERT_TRUE((result.has_value())) << "30! succeeds";
    if (result) {
        EXPECT_TRUE((result.value() == BigInt("265252859812191058636308480000000"))) << "30! retains every integer digit";
    }
}

TEST(TransformExactDecay, SymbolicDecayCoefficientAndCondition) {
    auto assumptions = std::make_shared<AssumptionContext>();
    EXPECT_TRUE((assumptions->assume_domain_checked("a", Domain::Real).has_value())) << "a is real";
    ComputationContext context;
    EXPECT_TRUE((context.set_assumptions(assumptions).has_value())) << "attach real facts";
    auto a = SymbolicExpr::variable("a");
    auto result = fourier_transform_checked(decay_input(a), "t", "omega", context);
    ASSERT_TRUE((result.has_value())) << "unknown positive real decay is a conditional result";
    if (!result)
        return;
    auto expression = result.value().value.expression;
    auto denominator = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(4), SymbolicExpr::power(a, SymbolicExpr::number(2))),
        SymbolicExpr::power(SymbolicExpr::variable("omega"), SymbolicExpr::number(2)));
    auto expected = SymbolicExpr::divide(SymbolicExpr::multiply(SymbolicExpr::number(12), a), denominator);
    auto equivalent = check_equivalent(expression, expected, context);
    EXPECT_TRUE((equivalent && std::holds_alternative<ProvedZeroResidual>(equivalent.value()))) << "Fourier result equals 12*a/(4*a^2+omega^2)";
    {
        EXPECT_NEAR(expression->substitute("a", SymbolicExpr::number(1))
                                         ->substitute("omega", SymbolicExpr::number(0))
                                         ->simplify()
                                         ->to_numeric(), 3.0, 1e-12) << "a=1 and omega=0 gives 3, retaining both independent coefficients";
    }
    const auto &conditions = result.value().value.conditions;
    EXPECT_TRUE((conditions.size() == 1)) << "one unknown convergence constraint survives scaling";
    if (conditions.size() != 1)
        return;
    auto relation = std::dynamic_pointer_cast<const RelationalNode>(detail::node(conditions[0]));
    EXPECT_TRUE((relation != nullptr)) << "convergence condition is relational";
    if (!relation)
        return;
    EXPECT_TRUE((relation->op() == RelationalNode::Op::GT)) << "strict positive convergence boundary";
    EXPECT_TRUE((relation->left()->equals(*detail::node(a)) && relation->right()->is_zero())) << "convergence constraint is a > 0";
}

TEST(TransformExactDecay, NonpositiveDecay) {
    auto assumptions = std::make_shared<AssumptionContext>();
    EXPECT_TRUE((assumptions->assume_domain_checked("a", Domain::Real).has_value())) << "a is real";
    EXPECT_TRUE((assumptions->assume_sign_checked("a", Sign::NonPositive).has_value())) << "a <= 0";
    ComputationContext context;
    EXPECT_TRUE((context.set_assumptions(assumptions).has_value())) << "attach nonpositive facts";
    auto result = fourier_transform_checked(decay_input(SymbolicExpr::variable("a")), "t", "omega", context);
    EXPECT_TRUE((!result && result.error().code == CasErrc::Inconclusive)) << "nonpositive parameter cannot use the convergent Fourier rule";
    auto growing = fourier_transform_checked(decay_input(SymbolicExpr::number(-1)), "t", "omega");
    EXPECT_TRUE((!growing && growing.error().code == CasErrc::Inconclusive)) << "numeric exponential growth is not a convergent transform";
    auto zero = fourier_transform_checked(decay_input(SymbolicExpr::number(0)), "t", "omega");
    EXPECT_TRUE((!zero && zero.error().code == CasErrc::Inconclusive)) << "zero decay is not integrable";
}

TEST(TransformExactDecay, UnprovedAndComplexDecay) {
    auto unknown = fourier_transform_checked(decay_input(SymbolicExpr::variable("a")), "t", "omega");
    EXPECT_TRUE((!unknown && unknown.error().code == CasErrc::Inconclusive)) << "unconstrained symbol is not assumed real";
    auto complex_rate = detail::make_expression_ptr(SymbolicFactory::create_complex(
        SymbolicFactory::create_number(BigInt(1)), SymbolicFactory::create_number(BigInt(1))));
    auto complex_result = fourier_transform_checked(decay_input(complex_rate), "t", "omega");
    EXPECT_TRUE((!complex_result && complex_result.error().code == CasErrc::Inconclusive)) << "nonreal decay is explicitly unsupported";
}

TEST(TransformExactDecay, DecaySumConditions) {
    auto assumptions = std::make_shared<AssumptionContext>();
    EXPECT_TRUE((assumptions->assume_domain_checked("a", Domain::Real).has_value())) << "a is real";
    EXPECT_TRUE((assumptions->assume_domain_checked("b", Domain::Real).has_value())) << "b is real";
    ComputationContext context;
    EXPECT_TRUE((context.set_assumptions(assumptions).has_value())) << "attach both real parameters";
    auto sum = SymbolicExpr::add(decay_input(SymbolicExpr::variable("a")),
                                 decay_input(SymbolicExpr::variable("b")));
    auto result = fourier_transform_checked(sum, "t", "omega", context);
    ASSERT_TRUE((result.has_value())) << "sum has a conditional transform";
    if (!result)
        return;
    const auto &conditions = result.value().value.conditions;
    EXPECT_TRUE((conditions.size() == 2)) << "both convergence constraints survive addition";
    AssumptionContext values;
    EXPECT_TRUE((values.assume_sign_checked("a", Sign::Positive).has_value())) << "positive a";
    EXPECT_TRUE((values.assume_sign_checked("b", Sign::Negative).has_value())) << "negative b";
    std::size_t rejected = 0;
    for (const auto &condition : conditions) {
        if (values.evaluate_condition(*condition) == Tribool::False)
            ++rejected;
    }
    EXPECT_TRUE((rejected == 1)) << "the negative b summand remains excluded by its condition";
}
