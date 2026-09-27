#include "test_common.hpp"
#include "exact_factorization.hpp"
#include "transcendental_factor.hpp"
#include <cstdint>

using namespace LMCAS;

static std::vector<Polynomial<Rational>> checked_zassenhaus(
    const Polynomial<Rational> &poly,
    const std::vector<Polynomial<BigInt>> &lifted,
    std::int64_t modulus) {
    auto result = zassenhaus_combine_checked(
        poly, lifted, BigInt(modulus));
    if (!result) {
        throw std::runtime_error(result.error().message);
    }
    return std::move(result.value().value);
}

TEST(TranscendentalFactorZassenhaus, ZassenhausSingleFactor) {
    Polynomial<Rational> poly({Rational(1), Rational(1)}, "x");
    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")};

    auto result = checked_zassenhaus(poly, lifted, 125); /**< 模数 p^k = 5^3。 */

    EXPECT_TRUE((result.size() == 1)) << "single factor should return 1 true factor";
    EXPECT_TRUE((result[0].degree() == 1)) << "factor should be degree 1";
}

TEST(TranscendentalFactorZassenhaus, ZassenhausTwoLinearFactors) {
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")};

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE((result.size() == 2)) << "x^2 - 1 should have 2 true factors";

    if (result.size() == 2) {
        auto product = result[0] * result[1];
        auto monic_poly = poly.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE((monic_poly == monic_product)) << "product of factors should equal original polynomial";
    }
}

TEST(TranscendentalFactorZassenhaus, ZassenhausIrreducibleQuadratic) {
    Polynomial<Rational> poly({Rational(1), Rational(1), Rational(1)}, "x");

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-18), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(-30), BigInt(1)}, "x")};

    auto result = checked_zassenhaus(poly, lifted, 49);

    EXPECT_TRUE((result.size() == 1)) << "irreducible polynomial should return 1 factor";
    EXPECT_TRUE((result[0].degree() == 2)) << "factor should be degree 2";
}

TEST(TranscendentalFactorZassenhaus, ZassenhausThreeFactors) {
    Polynomial<Rational> poly({Rational(0), Rational(-1), Rational(0), Rational(1)}, "x");

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(0), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")};

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE((result.size() == 3)) << "x^3 - x should have 3 true factors";

    int linear_count = 0;
    for (const auto &f : result) {
        if (f.degree() == 1)
            linear_count++;
    }
    EXPECT_TRUE((linear_count == 3)) << "all 3 factors should be linear";
}

TEST(TranscendentalFactorZassenhaus, ZassenhausEmptyFactors) {
    Polynomial<Rational> poly({Rational(1), Rational(0), Rational(1)}, "x");
    std::vector<Polynomial<BigInt>> lifted;

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE((result.size() == 1)) << "empty factors should return original poly";
}

TEST(TranscendentalFactorZassenhaus, ZassenhausZeroPolynomial) {
    Polynomial<Rational> poly("x");
    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")};

    auto result =
        zassenhaus_combine_checked(poly, lifted, BigInt(125));
    EXPECT_FALSE((result.has_value())) << "零多项式与非空提升因子不一致";
    EXPECT_TRUE((result.error().code == CasErrc::InvalidArgument)) << "不一致的提升因子应报告 InvalidArgument";
}

TEST(TranscendentalFactorZassenhaus, ZassenhausQuadraticTimesLinear) {
    Polynomial<Rational> x2_plus_1({Rational(1), Rational(0), Rational(1)}, "x");
    Polynomial<Rational> x_minus_1({Rational(-1), Rational(1)}, "x");
    Polynomial<Rational> poly = x2_plus_1 * x_minus_1;

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-57), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(-68), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x")};

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE((result.size() >= 1)) << "should find at least 1 factor";

    if (result.size() >= 1) {
        Polynomial<Rational> product = result[0];
        for (size_t i = 1; i < result.size(); ++i) {
            product = product * result[i];
        }
        auto monic_poly = poly.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE((monic_poly == monic_product)) << "product of all returned factors should equal original polynomial";
    }
}

TEST(TranscendentalFactorZassenhaus, ZassenhausEarlyTerminationIrreducible) {
    Polynomial<Rational> poly({Rational(1), Rational(1), Rational(1)}, "x");

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-18), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(-30), BigInt(1)}, "x")};

    auto result = checked_zassenhaus(poly, lifted, 49);

    EXPECT_TRUE((result.size() == 1)) << "irreducible polynomial should return 1 factor";
    EXPECT_TRUE((result[0].degree() == 2)) << "factor should be the original degree-2 polynomial";
    EXPECT_TRUE((result[0].coeffs[2] == Rational(1))) << "factor should be monic";
}

TEST(TranscendentalFactorZassenhaus, ZassenhausEarlyTerminationRemainingIrreducible) {
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(0), Rational(1)}, "x");

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(-18), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(-30), BigInt(1)}, "x")};

    auto result = checked_zassenhaus(poly, lifted, 49);

    EXPECT_TRUE((result.size() == 2)) << "should find 2 true factors";

    if (result.size() == 2) {
        auto product = result[0] * result[1];
        auto monic_poly = poly.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE((monic_poly == monic_product)) << "product of factors should equal original polynomial";
    }
}

TEST(TranscendentalFactorZassenhaus, ZassenhausEarlyTerminationLinearRemaining) {
    Polynomial<Rational> x2_minus_1({Rational(-1), Rational(0), Rational(1)}, "x");
    Polynomial<Rational> x_plus_2({Rational(2), Rational(1)}, "x");
    Polynomial<Rational> poly = x2_minus_1 * x_plus_2;

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(2), BigInt(1)}, "x")};

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE((result.size() >= 2)) << "should find at least 2 factors";

    Polynomial<Rational> product = result[0];
    for (size_t i = 1; i < result.size(); ++i) {
        product = product * result[i];
    }
    auto monic_poly = poly.make_monic();
    auto monic_product = product.make_monic();
    EXPECT_TRUE((monic_poly == monic_product)) << "product of all factors should equal original polynomial";
}

TEST(TranscendentalFactorZassenhaus, ZassenhausEarlyTerminationSingleActiveFactor) {
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")};

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE((result.size() == 2)) << "should find 2 factors";

    if (result.size() == 2) {
        auto product = result[0] * result[1];
        auto monic_poly = poly.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE((monic_poly == monic_product)) << "product of factors should equal original polynomial";
    }
}

TEST(TranscendentalFactorZassenhaus, ZassenhausEarlyTerminationDegreeOneInput) {
    Polynomial<Rational> poly({Rational(5), Rational(1)}, "x");

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(5), BigInt(1)}, "x")};

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE((result.size() == 1)) << "linear polynomial should return 1 factor";
    EXPECT_TRUE((result[0].degree() == 1)) << "factor should be linear";
}

TEST(TranscendentalFactorZassenhaus, ZassenhausRationalReconstructionIntegerCoeffs) {
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")};

    auto result = checked_zassenhaus(poly, lifted, 125);

    EXPECT_TRUE((result.size() == 2)) << "x^2 - 1 should still have 2 factors with rational reconstruction";

    if (result.size() == 2) {
        auto product = result[0] * result[1];
        auto monic_poly = poly.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE((monic_poly == monic_product)) << "product of factors should equal original polynomial";
    }
}

TEST(TranscendentalFactorZassenhaus, ZassenhausLargePrimeModulus) {
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");
    std::vector<Polynomial<BigInt>> lifted = {
        Polynomial<BigInt>({BigInt(-1), BigInt(1)}, "x"),
        Polynomial<BigInt>({BigInt(1), BigInt(1)}, "x")};

    auto result = checked_zassenhaus(poly, lifted, 1000000007);

    EXPECT_TRUE((result.size() == 2)) << "should find 2 factors with large prime power";

    if (result.size() == 2) {
        auto product = result[0] * result[1];
        auto monic_poly = poly.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE((monic_poly == monic_product)) << "product of factors should equal original with large prime power";
    }
}

TEST(TranscendentalFactorZassenhaus, ZassenhausBoundedEnumerationManyFactors) {
    Polynomial<Rational> poly({Rational(1)}, "x");
    for (int i = 1; i <= 16; ++i) {
        Polynomial<Rational> factor({Rational(-i), Rational(1)}, "x");
        poly = poly * factor;
    }

    std::vector<Polynomial<BigInt>> lifted;
    for (int i = 1; i <= 16; ++i) {
        lifted.push_back(Polynomial<BigInt>({BigInt(-i), BigInt(1)}, "x"));
    }

    const std::int64_t prime_power = 1000000007;

    auto result = checked_zassenhaus(poly, lifted, prime_power);

    EXPECT_TRUE((result.size() == 16)) << "bounded enumeration finds 16 factors";

    int linear_count = 0;
    for (const auto &f : result) {
        if (f.degree() == 1)
            linear_count++;
    }
    EXPECT_TRUE((linear_count == 16)) << "all 16 factors should be linear";
}

TEST(TranscendentalFactorZassenhaus, ZassenhausBoundedEnumerationProduct) {
    Polynomial<Rational> poly({Rational(1)}, "x");
    for (int i = 1; i <= 16; ++i) {
        Polynomial<Rational> factor({Rational(-i), Rational(1)}, "x");
        poly = poly * factor;
    }

    std::vector<Polynomial<BigInt>> lifted;
    for (int i = 1; i <= 16; ++i) {
        lifted.push_back(Polynomial<BigInt>({BigInt(-i), BigInt(1)}, "x"));
    }

    const std::int64_t prime_power = 1000000007;

    auto result = checked_zassenhaus(poly, lifted, prime_power);

    if (!result.empty()) {
        Polynomial<Rational> product = result[0];
        for (size_t i = 1; i < result.size(); ++i) {
            product = product * result[i];
        }
        auto monic_poly = poly.make_monic();
        auto monic_product = product.make_monic();
        EXPECT_TRUE((monic_poly == monic_product)) << "bounded-enumeration factors reconstruct original polynomial";
    }
}

TEST(TranscendentalFactorZassenhaus, ZassenhausTinyBudgetIsInconclusive) {
    Polynomial<Rational> poly({Rational(1)}, "x");
    std::vector<Polynomial<BigInt>> lifted;
    for (int i = 1; i <= 4; ++i) {
        poly = poly *
               Polynomial<Rational>({Rational(-i), Rational(1)}, "x");
        lifted.emplace_back(
            std::vector<BigInt>{BigInt(-i), BigInt(1)}, "x");
    }
    ResourceLimits limits;
    limits.max_steps = lifted.size();
    ComputationContext context(limits);
    auto result = zassenhaus_combine_checked(
        poly, lifted, BigInt(1000000007), context);

    ASSERT_TRUE((result.has_value())) << "预算耗尽应返回精确部分结果";
    if (result) {
        EXPECT_TRUE((result.value().completeness == Completeness::Inconclusive)) << "预算耗尽应标记为 Inconclusive";
        EXPECT_TRUE((result.value().value.size() == 1 &&
                     result.value().value[0] == poly.make_monic()))
            << "未决结果仍应精确重构输入";
    }
}
