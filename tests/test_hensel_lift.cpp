#include "test_hensel_lift_support.hpp"

namespace {

std::vector<Polynomial<ModInt>> factors_of_x_squared_minus_one(
    std::int64_t modulus) {
    Polynomial<ModInt> x_plus_one("x");
    x_plus_one.coeffs = {
        ModInt(1, modulus), ModInt(1, modulus)};
    Polynomial<ModInt> x_minus_one("x");
    x_minus_one.coeffs = {
        ModInt(modulus - 1, modulus), ModInt(1, modulus)};
    return {x_plus_one, x_minus_one};
}

} // namespace

TEST(HenselLift, RejectsEmptyFactorList) {
    Polynomial<BigInt> polynomial(
        {BigInt(-1), BigInt(0), BigInt(1)}, "x");

    auto result = hensel_lift_checked(polynomial, {}, 3, 0);

    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::InvalidArgument);
}

TEST(HenselLift, RejectsZeroPolynomial) {
    Polynomial<BigInt> zero("x");

    auto result = hensel_lift_checked(
        zero, factors_of_x_squared_minus_one(3), 3, 0);

    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::InvalidArgument);
}

TEST(HenselLift, RejectsNonPrimeModulus) {
    Polynomial<BigInt> polynomial(
        {BigInt(-1), BigInt(0), BigInt(1)}, "x");

    auto result = hensel_lift_checked(
        polynomial, factors_of_x_squared_minus_one(4), 4, 0);

    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::InvalidArgument);
}

TEST(HenselLift, AutomaticLiftBoundExactlyReconstructsPolynomial) {
    Polynomial<BigInt> polynomial(
        {BigInt(-1), BigInt(0), BigInt(1)}, "x");

    auto result = hensel_lift_checked(
        polynomial, factors_of_x_squared_minus_one(3), 3, 0);

    ASSERT_TRUE(result) << result.error().message;
    ASSERT_EQ(result.value().size(), 2U);
    auto reconstructed = test_poly_mul(
        result.value()[0].coeffs, result.value()[1].coeffs);
    while (!reconstructed.empty() && reconstructed.back().is_zero()) {
        reconstructed.pop_back();
    }
    EXPECT_EQ(reconstructed, polynomial.coeffs);
}
