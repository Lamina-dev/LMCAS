#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include "symbolic.hpp"
#include "integration.hpp"
#include "test_common.hpp"

using namespace LMCAS;

static std::shared_ptr<SymbolicExpr> num(int n) { return SymbolicExpr::number(n); }

static void check_roundtrip(const std::string &name,
                            const std::shared_ptr<SymbolicExpr> &integrand,
                            const std::string &var,
                            const std::vector<double> &points) {
    LMCAS::Integrator integrator;
    auto primitive = integrator.integrate(*integrand, var);
    ASSERT_TRUE(primitive.has_value())
        << name << ": integration failed: " << primitive.error().message;
    auto derivative = primitive.value().differentiate(var);
    ASSERT_NE(derivative, nullptr) << name << ": cannot differentiate result";

    for (double point : points) {
        SCOPED_TRACE(name + " at " + std::to_string(point));
        auto value = SymbolicExpr::number(point);
        auto derivative_at = derivative->substitute(var, value);
        auto integrand_at = integrand->substitute(var, value);
        ASSERT_NE(derivative_at, nullptr);
        ASSERT_NE(integrand_at, nullptr);
        auto derivative_value = test_numeric_eval(derivative_at->simplify());
        auto integrand_value = test_numeric_eval(integrand_at->simplify());
        ASSERT_TRUE(derivative_value.has_value());
        ASSERT_TRUE(integrand_value.has_value());
        ASSERT_TRUE(std::isfinite(*derivative_value));
        ASSERT_TRUE(std::isfinite(*integrand_value));
        EXPECT_NEAR(*derivative_value, *integrand_value, 1e-8);
    }
}

static void check_evaluated(const std::string &name,
                            const std::shared_ptr<SymbolicExpr> &integrand,
                            const std::string &var) {
    LMCAS::Integrator integrator;
    auto primitive = integrator.integrate(*integrand, var);
    ASSERT_TRUE(primitive.has_value())
        << name << ": integration failed: " << primitive.error().message;
    EXPECT_FALSE(LMCAS::detail::contains_node_type<IntegralNode>(
        LMCAS::detail::node(primitive.value())))
        << name << " left an unevaluated integral";
}

TEST(LmcasTrigSubstitution, InverseCircularRadical) {
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, num(2));
    // ∫ 1/√(1-x²) dx = arcsin(x)  (domain |x|<1)

    auto rad = SymbolicExpr::add(num(1), SymbolicExpr::multiply(num(-1), x_sq));
    auto f = SymbolicExpr::power(rad, SymbolicExpr::number(Rational(-1, 2)));
    check_evaluated("∫1/√(1-x²)", f, "x");
    check_roundtrip("∫1/√(1-x²) = arcsin(x)", f, "x", {0.2, 0.4, 0.6, 0.8});
}

TEST(LmcasTrigSubstitution, InverseHyperbolicRadical) {
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, num(2));
    // ∫ 1/√(1+x²) dx = arcsinh(x) = ln(x+√(x²+1))

    auto rad = SymbolicExpr::add(num(1), x_sq);
    auto f = SymbolicExpr::power(rad, SymbolicExpr::number(Rational(-1, 2)));
    check_evaluated("∫1/√(1+x²)", f, "x");
    check_roundtrip("∫1/√(1+x²) = arcsinh(x)", f, "x", {0.3, 0.7, 1.5, 2.5});
}

TEST(LmcasTrigSubstitution, ExteriorHyperbolicRadical) {
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, num(2));
    // ∫ 1/√(x²-1) dx = ln(x+√(x²-1))  (domain x>1)

    auto rad = SymbolicExpr::add(x_sq, num(-1));
    auto f = SymbolicExpr::power(rad, SymbolicExpr::number(Rational(-1, 2)));
    check_evaluated("∫1/√(x²-1)", f, "x");
    check_roundtrip("∫1/√(x²-1)", f, "x", {1.5, 2.0, 3.0, 4.0});
}

TEST(LmcasTrigSubstitution, CircularRadical) {
    auto x = SymbolicExpr::variable("x");
    auto x_sq = SymbolicExpr::power(x, num(2));
    auto rad = SymbolicExpr::add(num(4), SymbolicExpr::multiply(num(-1), x_sq));
    auto f = SymbolicExpr::power(rad, SymbolicExpr::number(Rational(1, 2)));
    check_evaluated("∫√(4-x²)", f, "x");
    check_roundtrip("∫√(4-x²)", f, "x", {-1.5, -0.5, 0.5, 1.5});
}

TEST(LmcasTrigSubstitution, CosineHalfAngle) {
    auto x = SymbolicExpr::variable("x");
    // ∫ 1/(1+cos(x)) dx = tan(x/2)

    auto cosx = SymbolicExpr::cos(x);
    auto denom = SymbolicExpr::add(num(1), cosx);
    auto f = SymbolicExpr::divide(num(1), denom);
    check_evaluated("∫1/(1+cos x)", f, "x");
    check_roundtrip("∫1/(1+cos x)", f, "x", {0.3, 0.7, 1.0, 1.5});
}

TEST(LmcasTrigSubstitution, CotangentHalfAngle) {
    auto x = SymbolicExpr::variable("x");
    // ∫ 1/(1-cos(x)) dx = -cot(x/2)  (another fully-reducible rational-trig case)

    auto cosx = SymbolicExpr::cos(x);
    auto denom = SymbolicExpr::add(num(1), SymbolicExpr::multiply(num(-1), cosx));
    auto f = SymbolicExpr::divide(num(1), denom);
    check_evaluated("∫1/(1-cos x)", f, "x");
    check_roundtrip("∫1/(1-cos x)", f, "x", {0.5, 1.0, 1.5, 2.0});
}
