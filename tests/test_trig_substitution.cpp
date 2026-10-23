#include "expr.hpp"
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include "symbolic.hpp"
#include "integration.hpp"
#include "test_common.hpp"
#include "internal/normalization_utils.hpp"
#include "residual_verification.hpp"

using namespace LMCAS;

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

static void check_roundtrip(const std::string &name,
                            const std::shared_ptr<SymbolicExpr> &integrand,
                            const std::string &var,
                            const std::vector<double> &points) {
    LMCAS::Integrator integrator;
    ComputationContext context;
    auto primitive = integrator.integrate_checked(*integrand, var, context);
    ASSERT_TRUE(primitive.has_value())
        << name << ": integration failed: " << primitive.error().message;
    EXPECT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(primitive.value())));
    auto differentiated = LMCAS::differentiate(
        std::make_shared<SymbolicExpr>(std::move(primitive.value())), var, context);
    ASSERT_TRUE(differentiated) << name << ": " << differentiated.error().message;
    auto derivative = std::move(differentiated.value());
    auto residual = SymbolicExpr::add(derivative, SymbolicExpr::multiply(num(-1), integrand));
    auto proof = check_zero_residual(residual, context);
    ASSERT_TRUE(proof) << name << ": " << proof.error().message;
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(proof.value())) << name;

    for (double point : points) {
        SCOPED_TRACE(name + " at " + std::to_string(point));
        auto value = SymbolicExpr::number(point);
        auto derivative_at = LMCAS::substitute(derivative, var, value).value();
        auto integrand_at = LMCAS::substitute(integrand, var, value).value();
        ASSERT_NE(derivative_at, nullptr);
        ASSERT_NE(integrand_at, nullptr);
        auto derivative_value = test_numeric_eval(LMCAS::simplify(derivative_at).value());
        auto integrand_value = test_numeric_eval(LMCAS::simplify(integrand_at).value());
        ASSERT_TRUE(derivative_value.has_value());
        ASSERT_TRUE(integrand_value.has_value());
        ASSERT_TRUE(std::isfinite(*derivative_value));
        ASSERT_TRUE(std::isfinite(*integrand_value));
        EXPECT_NEAR(*derivative_value, *integrand_value, 1e-8);
    }
}


TEST(LmcasTrigSubstitution, InverseCircularRadical) {
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, num(2));
    // ∫ 1/√(1-x²) dx = arcsin(x)  (domain |x|<1)

    auto rad = SymbolicExpr::add(num(1), SymbolicExpr::multiply(num(-1), x_sq));
    auto f = SymbolicExpr::power(rad, SymbolicExpr::number(Rational(-1, 2)));
    check_roundtrip("∫1/√(1-x²) = arcsin(x)", f, "x", {0.2, 0.4, 0.6, 0.8});
}

TEST(LmcasTrigSubstitution, InverseHyperbolicRadical) {
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, num(2));
    // ∫ 1/√(1+x²) dx = arcsinh(x) = ln(x+√(x²+1))

    auto rad = SymbolicExpr::add(num(1), x_sq);
    auto f = SymbolicExpr::power(rad, SymbolicExpr::number(Rational(-1, 2)));
    check_roundtrip("∫1/√(1+x²) = arcsinh(x)", f, "x", {0.3, 0.7, 1.5, 2.5});
}

TEST(LmcasTrigSubstitution, ExteriorHyperbolicRadical) {
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, num(2));
    // ∫ 1/√(x²-1) dx = ln(x+√(x²-1))  (domain x>1)

    auto rad = SymbolicExpr::add(x_sq, num(-1));
    auto f = SymbolicExpr::power(rad, SymbolicExpr::number(Rational(-1, 2)));
    check_roundtrip("∫1/√(x²-1)", f, "x", {1.5, 2.0, 3.0, 4.0});
}

TEST(LmcasTrigSubstitution, CircularRadical) {
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, num(2));
    auto rad = SymbolicExpr::add(num(4), SymbolicExpr::multiply(num(-1), x_sq));
    auto f = SymbolicExpr::power(rad, SymbolicExpr::number(Rational(1, 2)));
    check_roundtrip("∫√(4-x²)", f, "x", {-1.5, -0.5, 0.5, 1.5});
}

TEST(LmcasTrigSubstitution, CosineHalfAngle) {
    auto x = SymbolicExpr::variable("x");
    // ∫ 1/(1+cos(x)) dx = tan(x/2)

    auto cosx = SymbolicExpr::cos(x);
    auto denom = SymbolicExpr::add(num(1), cosx);
    auto f = SymbolicExpr::divide(num(1), denom);
    check_roundtrip("∫1/(1+cos x)", f, "x", {0.3, 0.7, 1.0, 1.5});
}

TEST(LmcasTrigSubstitution, CotangentHalfAngle) {
    auto x = SymbolicExpr::variable("x");
    // ∫ 1/(1-cos(x)) dx = -cot(x/2)  (another fully-reducible rational-trig case)

    auto cosx = SymbolicExpr::cos(x);
    auto denom = SymbolicExpr::add(num(1), SymbolicExpr::multiply(num(-1), cosx));
    auto f = SymbolicExpr::divide(num(1), denom);
    check_roundtrip("∫1/(1-cos x)", f, "x", {0.5, 1.0, 1.5, 2.0});
}

TEST(LmcasTrigSubstitution, HyperbolicRadical) {
    auto x = SymbolicExpr::variable("x");
    auto radicand = SymbolicExpr::add(num(4), SymbolicExpr::power(x, num(2)));
    check_roundtrip("sqrt(x^2+4)",
        SymbolicExpr::power(radicand, SymbolicExpr::number(Rational(1, 2))),
        "x", {-2.0, -0.5, 0.5, 2.0});
}

TEST(LmcasTrigSubstitution, ExteriorRadical) {
    auto x = SymbolicExpr::variable("x");
    auto radicand = SymbolicExpr::add(SymbolicExpr::power(x, num(2)), num(-4));
    check_roundtrip("sqrt(x^2-4)",
        SymbolicExpr::power(radicand, SymbolicExpr::number(Rational(1, 2))),
        "x", {2.5, 3.0, 4.0, 5.0});
}

static void expect_original_integral(const ExprPtr& integrand) {
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(*integrand, "x", context);
    ASSERT_TRUE(result);
    auto integral = std::dynamic_pointer_cast<const IntegralNode>(detail::node(result.value()));
    ASSERT_NE(integral, nullptr);
    EXPECT_TRUE(integral->body()->equals(*detail::node(integrand)));
}

TEST(LmcasTrigSubstitution, ApproximateHalfPowerRetainsObligation) {
    auto x = SymbolicExpr::variable("x");
    auto radicand = SymbolicExpr::add(num(1), SymbolicExpr::power(x, num(2)));
    expect_original_integral(SymbolicExpr::power(
        radicand, SymbolicExpr::number(-0.5 + 1e-12)));
}

TEST(LmcasTrigSubstitution, ApproximateQuadraticCoefficientRetainsObligation) {
    auto x = SymbolicExpr::variable("x");
    auto radicand = SymbolicExpr::add(num(1), SymbolicExpr::multiply(
        SymbolicExpr::number(1.0 + 1e-12), SymbolicExpr::power(x, num(2))));
    expect_original_integral(SymbolicExpr::power(
        radicand, SymbolicExpr::number(Rational(-1, 2))));
}

TEST(LmcasTrigSubstitution, RationalThirdStaysExact) {
    auto x = SymbolicExpr::variable("x");
    auto radicand = SymbolicExpr::add(SymbolicExpr::number(Rational(1, 3)),
        SymbolicExpr::power(x, num(2)));
    auto integrand = SymbolicExpr::power(radicand, SymbolicExpr::number(Rational(-1, 2)));
    Integrator integrator;
    ComputationContext context;
    auto result = integrator.integrate_checked(*integrand, "x", context);
    ASSERT_TRUE(result);
    EXPECT_FALSE(detail::contains_node_type<IntegralNode>(detail::node(result.value())));
    EXPECT_FALSE(contains_inexact_number(detail::node(result.value())));
    auto derivative = differentiate(std::make_shared<SymbolicExpr>(result.value()), "x", context);
    ASSERT_TRUE(derivative);
    auto residual = SymbolicExpr::add(derivative.value(),
        SymbolicExpr::multiply(num(-1), integrand));
    auto proof = check_zero_residual(residual, context);
    ASSERT_TRUE(proof);
    EXPECT_TRUE(std::holds_alternative<ProvedZeroResidual>(proof.value()));
}
