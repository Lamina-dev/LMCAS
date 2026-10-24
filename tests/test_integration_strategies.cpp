#include "expr.hpp"
#include <cmath>
#include <variant>
#include <vector>
#include <memory>
#include <string>
#include "symbolic.hpp"
#include "integration.hpp"
#include "test_common.hpp"
#include "residual_verification.hpp"
#include "numeric_evaluation.hpp"
#include "assumption_context.hpp"
#include <limits>
#include "internal/normalization_utils.hpp"
#include <utility>
#include <set>

using namespace LMCAS;

static void expect_primitive_identity(const SymbolicExpr &expression, const std::string &variable) {
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(expression, variable, context);
    ASSERT_TRUE(result.has_value()) << result.error().message;
    ASSERT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));

    auto derivative = LMCAS::differentiate(std::make_shared<SymbolicExpr>(result.value()), variable).value();
    ASSERT_NE(derivative, nullptr) << "failed to differentiate integral";

    auto proof = check_equivalent(derivative, detail::make_expression_ptr(expression), context);
    ASSERT_TRUE(proof.has_value()) << "primitive identity proof completes";
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(proof.value()))
        << "derivative equals integrand on their defined domain";
}

TEST(IntegrationStrategies, DecimalCancellationIntegrates) {
    auto x = SymbolicExpr::variable("x");
    auto positive = SymbolicExpr::number(1e-13);
    auto negative = SymbolicExpr::number(-1e-13);
    auto integrand = detail::make_expression_ptr(detail::make_node<AddNode>(
        std::vector<std::shared_ptr<const SymbolicNode>>{
            detail::node(x), detail::node(positive), detail::node(negative)}));
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(*integrand, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    EXPECT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));
    auto value = evaluate_numeric(result.value(), {{"x", 2.0}}, context);
    ASSERT_TRUE(value) << value.error().message;
    EXPECT_DOUBLE_EQ(value.value().value, 2.0);
}


TEST(IntegrationStrategies, PowerRule) {
    auto x = detail::make_expression_ptr(*SymbolicExpr::variable("x"));
    auto expression = SymbolicExpr::power(x, detail::make_expression_ptr(*SymbolicExpr::number(2)));
    expect_primitive_identity(*expression, "x");
}

TEST(IntegrationStrategies, ExponentialIntegrationByParts) {
    auto x = detail::make_expression_ptr(*SymbolicExpr::variable("x"));
    auto expression = SymbolicExpr::multiply(x, SymbolicExpr::exp(x));
    expect_primitive_identity(*expression, "x");
}

TEST(IntegrationStrategies, CosineSubstitution) {
    auto x = detail::make_expression_ptr(*SymbolicExpr::variable("x"));
    auto squared = SymbolicExpr::power(x, detail::make_expression_ptr(*SymbolicExpr::number(2)));
    auto twice_x = SymbolicExpr::multiply(detail::make_expression_ptr(*SymbolicExpr::number(2)), x);
    auto expression = SymbolicExpr::multiply(SymbolicExpr::cos(squared), twice_x);
    expect_primitive_identity(*expression, "x");
}

TEST(IntegrationStrategies, PartialFractions) {
    auto x = detail::make_expression_ptr(*SymbolicExpr::variable("x"));
    auto squared = SymbolicExpr::power(x, detail::make_expression_ptr(*SymbolicExpr::number(2)));
    auto denominator = SymbolicExpr::add(squared, detail::make_expression_ptr(*SymbolicExpr::number(-1)));
    auto expression = SymbolicExpr::power(denominator, detail::make_expression_ptr(*SymbolicExpr::number(-1)));
    expect_primitive_identity(*expression, "x");
}

TEST(IntegrationStrategies, PositiveDiscriminantQuadraticUsesRealLogBranch) {
    auto x = SymbolicExpr::variable("x");
    auto denominator = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-2));
    auto integrand = SymbolicExpr::power(
        denominator, SymbolicExpr::number(-1));

    Integrator integrator;
    ComputationContext context;
    auto primitive = integrator.integrate_checked(
        *integrand, "x", context);
    ASSERT_TRUE(primitive.has_value()) << primitive.error().message;
    ASSERT_FALSE(detail::contains_node_type<IntegralNode>(
        detail::node(primitive.value())));
    auto derivative = LMCAS::differentiate(std::make_shared<SymbolicExpr>(primitive.value()), "x").value();
    ASSERT_NE(derivative, nullptr);

    for (double sample : {-3.0, -2.0, 2.0, 3.0}) {
        SCOPED_TRACE("x=" + std::to_string(sample));
        auto value = SymbolicExpr::number(sample);
        auto expected = test_numeric_eval(
            LMCAS::simplify(LMCAS::substitute(integrand, "x", value).value()).value());
        auto actual = test_numeric_eval(
            LMCAS::simplify(LMCAS::substitute(derivative, "x", value).value()).value());
        ASSERT_TRUE(expected.has_value());
        ASSERT_TRUE(actual.has_value());
        ASSERT_TRUE(std::isfinite(*expected));
        ASSERT_TRUE(std::isfinite(*actual));
        EXPECT_NEAR(*actual, *expected, 1e-10);
    }
}

TEST(IntegrationStrategies, LogarithmicIntegrationByParts) {
    auto x = detail::make_expression_ptr(*SymbolicExpr::variable("x"));
    auto expression = SymbolicExpr::ln(x);
    expect_primitive_identity(*expression, "x");
}

static void expect_power_exponent(const SymbolicExpr &expression, const Rational &exponent) {
    class FindPower final : public detail::RecursiveSymbolicVisitor {
      public:
        detail::SymbolicNodePtr expected;
        bool found = false;
        explicit FindPower(const Rational &value)
            : expected(detail::node(SymbolicExpr::number(value))) {}
        void visit(const PowerNode &node) override {
            auto base = std::dynamic_pointer_cast<const VariableNode>(node.base());
            if (base && base->name() == "x" && node.exponent()->equals(*expected)) {
                found = true;
            }
            detail::RecursiveSymbolicVisitor::visit(node);
        }
    } visitor(exponent);
    detail::node(expression)->accept(visitor);
    EXPECT_TRUE(visitor.found) << expression.to_string();
}

TEST(IntegrationStrategies, VariableExponentRetainsObligation) {
    auto x = SymbolicExpr::variable("x");
    auto source = SymbolicExpr::power(x, x);
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(*source, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    EXPECT_TRUE(structurally_equal(result.value(), *SymbolicExpr::make_integral(source, "x")));
}

TEST(IntegrationStrategies, ParameterConditions) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto exponent = SymbolicExpr::add(a, SymbolicExpr::number(-1));
    auto power = SymbolicExpr::power(x, exponent);
    auto affine = SymbolicExpr::sin(SymbolicExpr::add(
        SymbolicExpr::multiply(a, x), SymbolicExpr::number(1)));
    for (const auto &source : {power, affine}) {
        SCOPED_TRACE(source->to_string());
        Integrator integrator;
        ComputationContext unknown;
        auto unresolved = integrator.integrate_checked(*source, "x", unknown);
        ASSERT_TRUE(unresolved) << unresolved.error().message;
        EXPECT_TRUE(structurally_equal(unresolved.value(), *SymbolicExpr::make_integral(source, "x")));

        auto assumptions = std::make_shared<AssumptionContext>();
        ASSERT_TRUE(assumptions->assume_sign("a", Sign::Positive));
        ComputationContext known;
        ASSERT_TRUE(known.set_assumptions(assumptions));
        auto result = integrator.integrate_checked(*source, "x", known);
        ASSERT_TRUE(result) << result.error().message;
        ASSERT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));
        auto derivative = differentiate(detail::make_expression_ptr(result.value()), "x", known);
        ASSERT_TRUE(derivative) << derivative.error().message;
        auto proof = check_equivalent(derivative.value(), source, known);
        ASSERT_TRUE(proof) << proof.error().message;
        EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(proof.value()));
    }
}

TEST(IntegrationStrategies, DegenerateParametersAreConstants) {
    auto x = SymbolicExpr::variable("x");
    auto argument = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(0), x), SymbolicExpr::number(1));
    expect_primitive_identity(*SymbolicExpr::sin(argument), "x");
    expect_primitive_identity(*SymbolicExpr::power(SymbolicExpr::number(1), x), "x");
}

TEST(IntegrationStrategies, NegativeRealRadicalScale) {
    auto x = SymbolicExpr::variable("x");
    auto scale = SymbolicExpr::number(-2);
    auto source = SymbolicExpr::power(SymbolicExpr::add(
        SymbolicExpr::power(scale, SymbolicExpr::number(2)),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::power(x, SymbolicExpr::number(2)))),
        SymbolicExpr::number(Rational(-1, 2)));
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(*source, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    ASSERT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));
    auto derivative = differentiate(detail::make_expression_ptr(result.value()), "x", context);
    ASSERT_TRUE(derivative) << derivative.error().message;
    for (double point : {-1.5, -0.5, 0.0, 0.5, 1.5}) {
        auto actual = evaluate_numeric(*derivative.value(), {{"x", point}}, context);
        auto expected = evaluate_numeric(*source, {{"x", point}}, context);
        ASSERT_TRUE(actual) << actual.error().message;
        ASSERT_TRUE(expected) << expected.error().message;
        EXPECT_GT(actual.value().value, 0.0);
        EXPECT_NEAR(actual.value().value, expected.value().value, 1e-12);
    }
}

TEST(IntegrationStrategies, DecimalCoefficients) {
    auto x = SymbolicExpr::variable("x");
    for (const auto &entry : std::vector<std::pair<double, Rational>>{
             {0.1, Rational(1, 10)}, {0.5, Rational(1, 2)},
             {1.0000000000000002, Rational("1.0000000000000002")},
             {0.3333333333333333, Rational("0.3333333333333333")},
             {1.25e-20, Rational("1.25e-20")}, {1.25e20, Rational("1.25e20")}}) {
        const double coefficient = entry.first;
        SCOPED_TRACE(coefficient);
        auto source = SymbolicExpr::multiply(SymbolicExpr::number(coefficient), x);
        Integrator integrator;
        ComputationContext context;
        auto result = integrator.integrate_checked(*source, "x", context);
        ASSERT_TRUE(result) << result.error().message;
        ASSERT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));
        expect_power_exponent(result.value(), Rational(2));
        auto expected = SymbolicExpr::multiply(
            SymbolicExpr::number(entry.second / Rational(2)),
            SymbolicExpr::power(x, SymbolicExpr::number(2)));
        EXPECT_TRUE(test_proved_equivalent(detail::make_expression_ptr(result.value()), expected));
        auto value = evaluate_numeric(result.value(), {{"x", 2.0}}, context);
        ASSERT_TRUE(value) << value.error().message;
        EXPECT_DOUBLE_EQ(value.value().value, coefficient * 2.0);
    }
}

TEST(IntegrationStrategies, DecimalTenthsFitSmallIntegerBudget) {
    auto source = SymbolicExpr::multiply(SymbolicExpr::number(0.1), SymbolicExpr::variable("x"));
    ResourceLimits limits;
    limits.max_integer_bits = 8;
    ComputationContext context(limits);
    Integrator integrator;
    auto result = integrator.integrate_checked(*source, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    expect_power_exponent(result.value(), Rational(2));
    auto value = evaluate_numeric(result.value(), {{"x", 2.0}}, context);
    ASSERT_TRUE(value) << value.error().message;
    EXPECT_DOUBLE_EQ(value.value().value, 0.2);
}

TEST(IntegrationStrategies, DecimalExponentDoesNotSnap) {
    auto x = SymbolicExpr::variable("x");
    for (const auto &entry : std::vector<std::pair<double, Rational>>{
             {0.5, Rational(3, 2)},
             {1.0000000000000002, Rational("2.0000000000000002")},
             {0.3333333333333333, Rational("1.3333333333333333")}}) {
        SCOPED_TRACE(entry.first);
        Integrator integrator;
        ComputationContext context;
        auto result = integrator.integrate_checked(
            *SymbolicExpr::power(x, SymbolicExpr::number(entry.first)), "x", context);
        ASSERT_TRUE(result) << result.error().message;
        ASSERT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));
        expect_power_exponent(result.value(), entry.second);
    }
}

TEST(IntegrationStrategies, ExactThirdStaysExact) {
    auto x = SymbolicExpr::variable("x");
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(
        *SymbolicExpr::power(x, SymbolicExpr::number(2)), "x", context);
    ASSERT_TRUE(result) << result.error().message;
    auto expected = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 3)),
        SymbolicExpr::power(x, SymbolicExpr::number(3)));
    EXPECT_TRUE(test_proved_equivalent(detail::make_expression_ptr(result.value()), expected));
}

TEST(IntegrationStrategies, DecimalConversionErrors) {
    auto x = SymbolicExpr::variable("x");
    auto source = SymbolicExpr::multiply(SymbolicExpr::number(0.3333333333333333), x);
    Integrator integrator;
    ResourceLimits limits;
    limits.max_integer_bits = 8;
    ComputationContext small(limits);
    auto limited = integrator.integrate_checked(*source, "x", small);
    ASSERT_FALSE(limited);
    EXPECT_EQ(limited.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(small.recursion_depth(), 0u);

    CancellationToken token;
    token.cancel();
    ComputationContext cancelled({}, token);
    auto stopped = integrator.integrate_checked(*source, "x", cancelled);
    ASSERT_FALSE(stopped);
    EXPECT_EQ(stopped.error().code, CasErrc::Cancelled);
    EXPECT_EQ(cancelled.recursion_depth(), 0u);
}

TEST(IntegrationStrategies, TinyRationalIntegratesBeforeNumericUnderflow) {
    auto source = SymbolicExpr::multiply(
        SymbolicExpr::number(std::numeric_limits<double>::denorm_min()),
        SymbolicExpr::power(SymbolicExpr::variable("x"), SymbolicExpr::number(2)));
    auto source_value = evaluate_numeric(*source, {{"x", 1.0}});
    ASSERT_TRUE(source_value) << source_value.error().message;
    EXPECT_DOUBLE_EQ(source_value.value().value, std::numeric_limits<double>::denorm_min());
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(*source, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    expect_power_exponent(result.value(), Rational(3));
    auto expected = SymbolicExpr::multiply(
        SymbolicExpr::number(Rational("5e-324") / Rational(3)),
        SymbolicExpr::power(SymbolicExpr::variable("x"), SymbolicExpr::number(3)));
    EXPECT_TRUE(test_proved_equivalent(detail::make_expression_ptr(result.value()), expected));
    auto value = evaluate_numeric(result.value(), {{"x", 1.0}}, context);
    ASSERT_FALSE(value);
    EXPECT_EQ(value.error().code, CasErrc::NumericFailure);
}

TEST(IntegrationStrategies, DecimalCoefficientPreservesHugeOddExponent) {
    const BigInt degree(std::string("9007199254740992"));
    auto x = SymbolicExpr::variable("x");
    auto source = SymbolicExpr::multiply(SymbolicExpr::number(0.5),
        SymbolicExpr::power(x, SymbolicExpr::number(degree)));
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(*source, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    expect_power_exponent(result.value(), Rational(degree + BigInt(1)));
    auto value = evaluate_numeric(result.value(), {{"x", -1.0}}, context);
    ASSERT_TRUE(value) << value.error().message;
    EXPECT_LT(value.value().value, 0.0);
}

TEST(IntegrationStrategies, AbsoluteValueRequiresRealVariable) {
    auto x = SymbolicExpr::variable("x");
    auto source = LMCAS::abs(x);
    ASSERT_TRUE(source) << source.error().message;
    Integrator integrator;
    ComputationContext unknown;
    auto unresolved = integrator.integrate_checked(*source.value(), "x", unknown);
    ASSERT_TRUE(unresolved) << unresolved.error().message;
    EXPECT_TRUE(structurally_equal(unresolved.value(),
        *SymbolicExpr::make_integral(source.value(), "x")));
    auto assumptions = std::make_shared<AssumptionContext>();
    ASSERT_TRUE(assumptions->assume_domain_checked("x", Domain::Real));
    ComputationContext real;
    ASSERT_TRUE(real.set_assumptions(assumptions));
    auto result = integrator.integrate_checked(*source.value(), "x", real);
    ASSERT_TRUE(result) << result.error().message;
    ASSERT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));
    auto derivative = differentiate(detail::make_expression_ptr(result.value()), "x", real);
    ASSERT_TRUE(derivative) << derivative.error().message;
    for (double point : {-2.0, 2.0}) {
        auto value = evaluate_numeric(*derivative.value(), {{"x", point}}, real);
        ASSERT_TRUE(value) << value.error().message;
        EXPECT_DOUBLE_EQ(value.value().value, 2.0);
    }
}

TEST(IntegrationStrategies, LinearSubstitutionPreservesOccupiedNames) {
    auto x = SymbolicExpr::variable("x");
    auto occupied = SymbolicExpr::variable("__lin_sub_u__");
    auto source = SymbolicExpr::sin(SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x), occupied));
    auto assumptions = std::make_shared<AssumptionContext>();
    ASSERT_TRUE(assumptions->assume_sign("__lin_sub_u__1", Sign::Zero));
    ComputationContext context;
    ASSERT_TRUE(context.set_assumptions(assumptions));
    Integrator integrator;
    auto result = integrator.integrate_checked(*source, "x", context);
    ASSERT_TRUE(result) << result.error().message;
    ASSERT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));
    EXPECT_EQ(detail::all_variable_names(detail::node(result.value())),
        (std::set<std::string>{"x", "__lin_sub_u__"}));
    auto derivative = differentiate(detail::make_expression_ptr(result.value()), "x", context);
    ASSERT_TRUE(derivative) << derivative.error().message;
    auto proof = check_equivalent(derivative.value(), source, context);
    ASSERT_TRUE(proof) << proof.error().message;
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(proof.value()));
}
