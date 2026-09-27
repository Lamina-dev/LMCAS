#include <cmath>
#include <stdexcept>
#include <variant>
#include <memory>
#include <string>
#include "symbolic.hpp"
#include "integration.hpp"
#include "test_common.hpp"
#include "residual_verification.hpp"

using namespace LMCAS;

static void expect_primitive_identity(const SymbolicExpr &expression, const std::string &variable) {
    Integrator integrator;
    auto result = integrator.integrate(expression, variable);
    ASSERT_TRUE(result.has_value()) << result.error().message;

    auto derivative = result.value().differentiate(variable);
    ASSERT_NE(derivative, nullptr) << "failed to differentiate integral";

    ComputationContext context;
    auto proof = check_equivalent(derivative, detail::make_expression_ptr(expression), context);
    ASSERT_TRUE(proof.has_value()) << "primitive identity proof completes";
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(proof.value()))
        << "derivative equals integrand on their defined domain";
}

class WrongIntegrationStrategy final : public IntegrationStrategy {
  public:
    Result<std::shared_ptr<SymbolicExpr>> try_integrate_raw(
        const SymbolicExpr &,
        const std::string &variable,
        Integrator &,
        ComputationContext &,
        int) override {
        return SymbolicExpr::variable(variable);
    }

    std::string name() const override { return "DeliberatelyWrong"; }
};

class ThrowingIntegrationStrategy final : public IntegrationStrategy {
  public:
    enum class Failure { Runtime, Allocation };

    explicit ThrowingIntegrationStrategy(Failure failure)
        : failure_(failure) {}

    Result<std::shared_ptr<SymbolicExpr>> try_integrate_raw(
        const SymbolicExpr &,
        const std::string &,
        Integrator &,
        ComputationContext &,
        int) override {
        if (failure_ == Failure::Allocation) {
            throw std::bad_alloc{};
        }
        throw std::runtime_error("injected integration failure");
    }

    std::string name() const override { return "Throwing"; }

  private:
    Failure failure_;
};

TEST(IntegrationStrategies, GeneratedCandidatesRequireExactResidualProof) {
    auto x = detail::make_expression_ptr(*SymbolicExpr::variable("x"));
    Integrator integrator;
    auto added = integrator.add_strategy(std::make_unique<WrongIntegrationStrategy>(), 0);
    ASSERT_TRUE(added.has_value()) << "wrong strategy is installed for the gate test";
    ComputationContext context;
    auto result = integrator.integrate_checked(*x, "x", context);
    ASSERT_TRUE(result.has_value());
    EXPECT_NE(result.value().to_string(), "x") << "wrong integration candidate is rejected";
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
    auto derivative = primitive.value().differentiate("x");
    ASSERT_NE(derivative, nullptr);

    for (double sample : {-3.0, -2.0, 2.0, 3.0}) {
        SCOPED_TRACE("x=" + std::to_string(sample));
        auto value = SymbolicExpr::number(sample);
        auto expected = test_numeric_eval(
            integrand->substitute("x", value)->simplify());
        auto actual = test_numeric_eval(
            derivative->substitute("x", value)->simplify());
        ASSERT_TRUE(expected.has_value());
        ASSERT_TRUE(actual.has_value());
        ASSERT_TRUE(std::isfinite(*expected));
        ASSERT_TRUE(std::isfinite(*actual));
        EXPECT_NEAR(*actual, *expected, 1e-10);
    }
}

TEST(IntegrationStrategies, StrategyBoundaryPreservesFallbackAndTypedErrors) {
    auto x = SymbolicExpr::variable("x");
    Integrator integrator;
    ComputationContext context;

    PartialFractionStrategy partial_fraction;
    auto not_applicable = partial_fraction.try_integrate(
        *SymbolicExpr::sin(x), "x", integrator, context);
    ASSERT_TRUE(not_applicable.has_value());
    EXPECT_TRUE(std::holds_alternative<IntegrationNotApplicable>(
        not_applicable.value()));

    ThrowingIntegrationStrategy runtime(
        ThrowingIntegrationStrategy::Failure::Runtime);
    auto invariant = runtime.try_integrate(
        *x, "x", integrator, context);
    ASSERT_FALSE(invariant.has_value());
    EXPECT_EQ(invariant.error().code, CasErrc::InternalInvariant);
    EXPECT_EQ(invariant.error().message, "injected integration failure");
    EXPECT_EQ(invariant.error().operation, "integrate.strategy.Throwing");

    ThrowingIntegrationStrategy allocation(
        ThrowingIntegrationStrategy::Failure::Allocation);
    auto resource = allocation.try_integrate(
        *x, "x", integrator, context);
    ASSERT_FALSE(resource.has_value());
    EXPECT_EQ(resource.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(resource.error().operation, "integrate.strategy.Throwing");
}

TEST(IntegrationStrategies, LogarithmicIntegrationByParts) {
    auto x = detail::make_expression_ptr(*SymbolicExpr::variable("x"));
    auto expression = SymbolicExpr::ln(x);
    expect_primitive_identity(*expression, "x");
}
