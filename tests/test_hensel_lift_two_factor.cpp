#include "test_hensel_lift_support.hpp"

/**
 * @brief 测试 x^2 - 1 = (x+1)(x-1) mod 3 的 Hensel 提升到 mod 9.
 *
 * Bezout: s=2, t=1 满足 2*(x+1) + 1*(x-1) = 3x+1 = 1 (mod 3)
 * 提升后应满足 f = g'*h' (mod 9) 且 s'*g' + t'*h' = 1 (mod 9)
 */
TEST(HenselLiftTwoFactor, HenselLiftX2Minus1Mod3) {
    std::vector<BigInt> f = {BigInt(-1), BigInt(0), BigInt(1)};
    std::vector<BigInt> g = {BigInt(1), BigInt(1)};
    std::vector<BigInt> h = {BigInt(-1), BigInt(1)};
    std::vector<BigInt> s = {BigInt(2)};
    std::vector<BigInt> t = {BigInt(1)};
    BigInt m(3);

    EXPECT_TRUE((test_verify_factorization(f, g, h, m))) << "initial: g*h = f (mod 3)";
    EXPECT_TRUE((test_verify_bezout(s, g, t, h, m))) << "initial: s*g + t*h = 1 (mod 3)";

    HenselLiftPair initial{g, h, s, t, m};
    HenselLiftPair lifted = hl_two_factor_lift(f, initial);

    BigInt m2(9);
    EXPECT_TRUE((lifted.modulus == m2)) << "lifted modulus = 9";

    EXPECT_TRUE((test_verify_factorization(f, lifted.g, lifted.h, m2))) << "lifted: g'*h' = f (mod 9)";
    EXPECT_TRUE((test_verify_bezout(lifted.s, lifted.g, lifted.t, lifted.h, m2))) << "lifted: s'*g' + t'*h' = 1 (mod 9)";

    EXPECT_TRUE((lifted.g.size() == 2)) << "g' has degree 1";
    EXPECT_TRUE((lifted.h.size() == 2)) << "h' has degree 1";
}

/**
 * @brief 测试 x^2+3x+2 = (x+1)(x+2) mod 5 的 Hensel 提升到 mod 25.
 *
 * Bezout: s=4, t=1 满足 4*(x+1) + 1*(x+2) = 5x+6 = 1 (mod 5)
 */
TEST(HenselLiftTwoFactor, HenselLiftX2Plus3xPlus2Mod5) {
    std::vector<BigInt> f = {BigInt(2), BigInt(3), BigInt(1)};
    std::vector<BigInt> g = {BigInt(1), BigInt(1)};
    std::vector<BigInt> h = {BigInt(2), BigInt(1)};
    std::vector<BigInt> s = {BigInt(4)};
    std::vector<BigInt> t_coeff = {BigInt(1)};
    BigInt m(5);
    BigInt m2(25);

    EXPECT_TRUE((test_verify_factorization(f, g, h, m))) << "initial: g*h = f (mod 5)";
    EXPECT_TRUE((test_verify_bezout(s, g, t_coeff, h, m))) << "initial: s*g + t*h = 1 (mod 5)";

    HenselLiftPair initial{g, h, s, t_coeff, m};
    HenselLiftPair lifted = hl_two_factor_lift(f, initial);

    EXPECT_TRUE((lifted.modulus == m2)) << "lifted modulus = 25";

    EXPECT_TRUE((test_verify_factorization(f, lifted.g, lifted.h, m2))) << "lifted: g'*h' = f (mod 25)";
    EXPECT_TRUE((test_verify_bezout(lifted.s, lifted.g, lifted.t, lifted.h, m2))) << "lifted: s'*g' + t'*h' = 1 (mod 25)";
}

/**
 * @brief 将 x^2+6x+5 的模 3 因子 (x+1)(x+2) 提升为模 9 因子 (x+1)(x+5)。
 */
TEST(HenselLiftTwoFactor, HenselLiftNonExactMod3) {
    std::vector<BigInt> f = {BigInt(5), BigInt(6), BigInt(1)};
    std::vector<BigInt> g = {BigInt(1), BigInt(1)};
    std::vector<BigInt> h = {BigInt(2), BigInt(1)};
    std::vector<BigInt> s = {BigInt(2)};
    std::vector<BigInt> t_coeff = {BigInt(1)};
    BigInt m(3);
    BigInt m2(9);

    EXPECT_TRUE((test_verify_factorization(f, g, h, m))) << "initial: g*h = f (mod 3)";
    EXPECT_TRUE((test_verify_bezout(s, g, t_coeff, h, m))) << "initial: s*g + t*h = 1 (mod 3)";

    HenselLiftPair initial{g, h, s, t_coeff, m};
    HenselLiftPair lifted = hl_two_factor_lift(f, initial);

    EXPECT_TRUE((lifted.modulus == m2)) << "lifted modulus = 9";

    EXPECT_TRUE((test_verify_factorization(f, lifted.g, lifted.h, m2))) << "lifted: g'*h' = f (mod 9)";
    EXPECT_TRUE((test_verify_bezout(lifted.s, lifted.g, lifted.t, lifted.h, m2))) << "lifted: s'*g' + t'*h' = 1 (mod 9)";

    EXPECT_TRUE((test_verify_factorization(f, {BigInt(1), BigInt(1)}, {BigInt(5), BigInt(1)}, m2))) << "exact: (x+1)*(x+5) = f (mod 9)";
}

TEST(HenselLiftTwoFactor, HenselLiftIterated) {
    /**
     * @brief 取 f = x^2+10x+21 = (x+3)(x+7)，模 5 下为 (x+3)(x+2)。
     * 模 5 根为 2、3；Bezout 系数 s=1、t=4 满足
     * 1*(x+3) + 4*(x+2) = 5x+11 = 1 (mod 5)。
     */

    std::vector<BigInt> f = {BigInt(21), BigInt(10), BigInt(1)};
    std::vector<BigInt> g = {BigInt(3), BigInt(1)};
    std::vector<BigInt> h = {BigInt(2), BigInt(1)};
    std::vector<BigInt> s = {BigInt(1)};
    std::vector<BigInt> t_coeff = {BigInt(4)};
    BigInt m(5);

    EXPECT_TRUE((test_verify_factorization(f, g, h, m))) << "initial: g*h = f (mod 5)";
    EXPECT_TRUE((test_verify_bezout(s, g, t_coeff, h, m))) << "initial: s*g + t*h = 1 (mod 5)";

    HenselLiftPair state{g, h, s, t_coeff, m};
    HenselLiftPair lifted1 = hl_two_factor_lift(f, state);

    BigInt m2(25);
    EXPECT_TRUE((lifted1.modulus == m2)) << "first lift modulus = 25";
    EXPECT_TRUE((test_verify_factorization(f, lifted1.g, lifted1.h, m2))) << "first lift: g'*h' = f (mod 25)";
    EXPECT_TRUE((test_verify_bezout(lifted1.s, lifted1.g, lifted1.t, lifted1.h, m2))) << "first lift: s'*g' + t'*h' = 1 (mod 25)";

    HenselLiftPair lifted2 = hl_two_factor_lift(f, lifted1);

    BigInt m4(625);
    EXPECT_TRUE((lifted2.modulus == m4)) << "second lift modulus = 625";
    EXPECT_TRUE((test_verify_factorization(f, lifted2.g, lifted2.h, m4))) << "second lift: g'*h' = f (mod 625)";
    EXPECT_TRUE((test_verify_bezout(lifted2.s, lifted2.g, lifted2.t, lifted2.h, m4))) << "second lift: s'*g' + t'*h' = 1 (mod 625)";
}

/**
 * @brief 测试三次多项式的 Hensel 提升.
 *
 * f = x^3 - x = x(x+1)(x-1),取两因子 g=x, h=x^2-1 mod 5.
 */
TEST(HenselLiftTwoFactor, HenselLiftCubic) {
    std::vector<BigInt> f = {BigInt(0), BigInt(-1), BigInt(0), BigInt(1)};
    std::vector<BigInt> g = {BigInt(0), BigInt(1)};
    std::vector<BigInt> h = {BigInt(-1), BigInt(0), BigInt(1)};

    /**
     * @brief 取 s=x、t=4，满足 x*x + 4*(x^2-1) = 5x^2-4 = 1 (mod 5)。
     * 次数满足 deg(s) < deg(h)=2、deg(t) < deg(g)=1。
     */
    std::vector<BigInt> s = {BigInt(0), BigInt(1)};
    std::vector<BigInt> t_coeff = {BigInt(4)};
    BigInt m(5);

    EXPECT_TRUE((test_verify_factorization(f, g, h, m))) << "initial: g*h = f (mod 5)";
    EXPECT_TRUE((test_verify_bezout(s, g, t_coeff, h, m))) << "initial: s*g + t*h = 1 (mod 5)";

    HenselLiftPair state{g, h, s, t_coeff, m};
    HenselLiftPair lifted = hl_two_factor_lift(f, state);

    BigInt m2(25);
    EXPECT_TRUE((lifted.modulus == m2)) << "lifted modulus = 25";
    EXPECT_TRUE((test_verify_factorization(f, lifted.g, lifted.h, m2))) << "lifted: g'*h' = f (mod 25)";
    EXPECT_TRUE((test_verify_bezout(lifted.s, lifted.g, lifted.t, lifted.h, m2))) << "lifted: s'*g' + t'*h' = 1 (mod 25)";
}

TEST(HenselLiftTwoFactor, ValidFactorsAndMismatchedProduct) {
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(-1), BigInt(0), BigInt(1)}, "x");

    int64_t p = 3;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(1, p), ModInt(1, p)};
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(2, p), ModInt(1, p)};
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2};
    auto checked = hensel_lift_checked(
        poly, mod_factors, p, 2);
    ASSERT_TRUE(checked) << checked.error().message;
    EXPECT_EQ(checked.value().size(), 2U);

    Polynomial<ModInt> bad_factor("x");
    bad_factor.coeffs = {ModInt(0, p), ModInt(1, p)};
    auto mismatch = hensel_lift_checked(
        poly, {f1, bad_factor}, p, 2);
    ASSERT_FALSE(mismatch);
    EXPECT_EQ(
        mismatch.error().code, CasErrc::InvalidArgument);
}
