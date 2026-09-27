#include "poly_utils.hpp"
#include "polynomial.hpp"
#include "symbolic.hpp"
#include "test_common.hpp"

using namespace LMCAS;

TEST(LmcasPolyBridge, Bridge) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::number(1)));

    auto converted = symbolic_to_poly<BigInt>(expression, "x");

    ASSERT_TRUE(converted) << converted.error().message;
    const auto &polynomial = converted.value();
    EXPECT_EQ(
        polynomial.coeffs,
        (std::vector<BigInt>{BigInt(1), BigInt(2), BigInt(1)}));
    auto roundtrip = poly_to_symbolic(polynomial);
    EXPECT_TRUE(test_proved_equivalent(roundtrip, expression));

    Polynomial<BigInt> x_plus_one(
        {BigInt(1), BigInt(1)}, "x");
    auto gcd = Polynomial<BigInt>::gcd(polynomial, x_plus_one);
    EXPECT_TRUE(gcd == x_plus_one);

    auto polynomial_division = polynomial.div_mod(gcd);
    auto factor_division = x_plus_one.div_mod(gcd);
    EXPECT_TRUE(polynomial_division.second.is_zero());
    EXPECT_TRUE(factor_division.second.is_zero());
    EXPECT_TRUE(polynomial_division.first * gcd == polynomial);
    EXPECT_TRUE(factor_division.first * gcd == x_plus_one);

    auto gcd_expression = poly_to_symbolic(gcd);
    EXPECT_TRUE(test_proved_equivalent(
        gcd_expression, SymbolicExpr::add(
                            x, SymbolicExpr::number(1))));

    ComputationContext context;
    auto symbolic_gcd = symbolic_polynomial_gcd(
        *expression, *gcd_expression, context);
    ASSERT_TRUE(symbolic_gcd) << symbolic_gcd.error().message;
    EXPECT_TRUE(test_proved_equivalent(
        symbolic_gcd.value(), gcd_expression));
}
