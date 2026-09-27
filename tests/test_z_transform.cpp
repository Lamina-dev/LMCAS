#include "test_common.hpp"
#include "transform_engine.hpp"
#include <limits>

using namespace LMCAS;

static double z_value(const std::shared_ptr<SymbolicExpr> &sequence, double z) {
    auto result = z_transform_checked(sequence, "n", "z");
    EXPECT_TRUE((result.has_value())) << "supported sequence has an evaluated Z transform";
    if (!result)
        return std::numeric_limits<double>::quiet_NaN();
    return result.value().value.expression->substitute("z", SymbolicExpr::number(z))->simplify()->to_numeric();
}

TEST(ZTransform, ConstantSequences) {
    EXPECT_NEAR(z_value(SymbolicExpr::number(1), 2), 2.0, 1e-12) << "Z{1}(2)=2";
    EXPECT_NEAR(z_value(SymbolicExpr::number(5), 3), 7.5, 1e-12) << "Z{5}(3)=15/2";
}

TEST(ZTransform, ExponentialSequence) {
    auto sequence = SymbolicExpr::power(SymbolicExpr::number(0.5), SymbolicExpr::variable("n"));
    {
        const double actual_value = (z_value(sequence, 2));
        const double expected_value = (4.0 / 3.0);
        const double tolerance = (1e-12);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
}

TEST(ZTransform, PolynomialSequences) {
    auto n = SymbolicExpr::variable("n");
    EXPECT_NEAR(z_value(n, 3), 0.75, 1e-12) << "Z{n}(3)=3/4";
    EXPECT_NEAR(z_value(SymbolicExpr::power(n, SymbolicExpr::number(2)), 3), 1.5, 1e-12) << "Z{n^2}(3)=3*(3+1)/(3-1)^3";
}

TEST(ZTransform, PolynomialExponentMustBeExactIntegerTwo) {
    auto near_two = SymbolicExpr::number(Rational(
        BigInt("20000000000001"), BigInt("10000000000000")));
    auto sequence = SymbolicExpr::power(
        SymbolicExpr::variable("n"), near_two);

    auto result = z_transform_checked(sequence, "n", "z");
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::Inconclusive)
        << "n^(2+10^-13) cannot use the quadratic-sequence formula";
}

TEST(ZTransform, TrigonometricSequences) {
    auto angle = SymbolicExpr::multiply(SymbolicExpr::number(0.7), SymbolicExpr::variable("n"));
    const double denominator = 5.0 - 4.0 * std::cos(0.7);
    {
        const double actual_value = (z_value(SymbolicExpr::sin(angle), 2));
        const double expected_value = (2.0 * std::sin(0.7) / denominator);
        const double tolerance = (1e-12);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
    {
        const double actual_value = (z_value(SymbolicExpr::cos(angle), 2));
        const double expected_value = (2.0 * (2.0 - std::cos(0.7)) / denominator);
        const double tolerance = (1e-12);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
}

TEST(ZTransform, Linearity) {
    auto n = SymbolicExpr::variable("n");
    auto sum = SymbolicExpr::add(SymbolicExpr::number(1), n);
    EXPECT_NEAR(z_value(sum, 3), 2.25, 1e-12) << "Z{1+n}(3)=3/2+3/4";
    auto scaled = SymbolicExpr::multiply(SymbolicExpr::number(3),
                                         SymbolicExpr::power(SymbolicExpr::number(2), n));
    EXPECT_NEAR(z_value(scaled, 4), 6.0, 1e-12) << "Z{3*2^n}(4)=6";
}

TEST(ZTransform, WeightedExponential) {
    auto n = SymbolicExpr::variable("n");
    auto sequence = SymbolicExpr::multiply(n, SymbolicExpr::power(SymbolicExpr::number(2), n));
    EXPECT_NEAR(z_value(sequence, 4), 2.0, 1e-12) << "Z{n*2^n}(4)=2*4/(4-2)^2";
}

TEST(ZTransform, UnsupportedSequence) {
    auto result = z_transform_checked(SymbolicExpr::ln(SymbolicExpr::variable("n")), "n", "z");
    EXPECT_TRUE((!result && result.error().code == CasErrc::Inconclusive)) << "Z{ln(n)} cannot fabricate an evaluated result";
}
