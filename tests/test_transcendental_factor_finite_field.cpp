#include "test_common.hpp"
#include "exact_factorization.hpp"
#include "transcendental_factor.hpp"
#include "internal/multivariate_factor_support.hpp"

using namespace LMCAS;

TEST(TranscendentalFactorFiniteField, SquareFreeAlreadySquareFree) {
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = tf_square_free(poly);

    EXPECT_FALSE((result.had_repeated_factors)) << "x^2 - 1 should be square-free";
    EXPECT_TRUE((result.square_free.degree() == 2)) << "square-free part degree should be 2";
    EXPECT_TRUE((result.repeated_factor.degree() == 0)) << "repeated factor should be constant";
}

TEST(TranscendentalFactorFiniteField, SquareFreeRepeatedRoot) {
    Polynomial<Rational> poly({Rational(1), Rational(-2), Rational(1)}, "x");

    auto result = tf_square_free(poly);

    EXPECT_TRUE((result.had_repeated_factors)) << "should detect repeated factors";
    EXPECT_TRUE((result.square_free.degree() == 1)) << "square-free part should be linear (x-1)";
    EXPECT_TRUE((result.square_free.coeffs[1] == Rational(1))) << "leading coeff should be 1";
    EXPECT_TRUE((result.square_free.coeffs[0] == Rational(-1))) << "constant term should be -1";
}

TEST(TranscendentalFactorFiniteField, SquareFreeHigherMultiplicity) {
    Polynomial<Rational> x_minus_1({Rational(-1), Rational(1)}, "x");
    Polynomial<Rational> x_plus_1({Rational(1), Rational(1)}, "x");
    auto x_minus_1_sq = x_minus_1 * x_minus_1;
    auto x_minus_1_cubed = x_minus_1_sq * x_minus_1;
    auto poly = x_minus_1_cubed * x_plus_1;

    auto result = tf_square_free(poly);

    EXPECT_TRUE((result.had_repeated_factors)) << "should detect repeated factors";
    EXPECT_TRUE((result.square_free.degree() == 2)) << "square-free part should be degree 2";
}

TEST(TranscendentalFactorFiniteField, SquareFreeConstantPolynomial) {
    Polynomial<Rational> poly({Rational(5)}, "x");

    auto result = tf_square_free(poly);

    EXPECT_FALSE((result.had_repeated_factors)) << "constant should have no repeated factors";
    EXPECT_TRUE((result.square_free.degree() == 0)) << "square-free part should be constant";
}

TEST(TranscendentalFactorFiniteField, SquareFreeLinearPolynomial) {
    Polynomial<Rational> poly({Rational(3), Rational(1)}, "x");

    auto result = tf_square_free(poly);

    EXPECT_FALSE((result.had_repeated_factors)) << "linear polynomial should be square-free";
    EXPECT_TRUE((result.square_free.degree() == 1)) << "square-free part should be linear";
    EXPECT_TRUE((result.square_free.coeffs[0] == Rational(3))) << "constant term preserved";
    EXPECT_TRUE((result.square_free.coeffs[1] == Rational(1))) << "linear term preserved";
}

TEST(TranscendentalFactorFiniteField, SquareFreeZeroPolynomial) {
    Polynomial<Rational> poly("x");

    auto result = tf_square_free(poly);

    EXPECT_FALSE((result.had_repeated_factors)) << "zero polynomial has no repeated factors";
    EXPECT_TRUE((result.square_free.is_zero())) << "square-free part of zero is zero";
}

TEST(TranscendentalFactorFiniteField, BerlekampLinearPoly) {
    Polynomial<Rational> poly({Rational(1), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 0);

    EXPECT_TRUE((result.prime > 0)) << "should select a valid prime";
    EXPECT_TRUE((result.factors.size() == 1)) << "linear polynomial should have 1 factor";
    EXPECT_TRUE((result.factors[0].degree() == 1)) << "factor should be degree 1";
}

TEST(TranscendentalFactorFiniteField, ExplicitPrimeRejectsRepeatedReduction) {
    Polynomial<Rational> polynomial(
        {Rational(4), Rational(5), Rational(1)}, "x");

    auto reduction = berlekamp_factor(polynomial, 3);

    EXPECT_TRUE(reduction.factors.empty())
        << "x^2+5x+4 becomes (x+1)^2 modulo three and is unsuitable";
}

TEST(TranscendentalFactorFiniteField, ExplicitPrimeRejectsDegreeDrop) {
    Polynomial<Rational> polynomial(
        {Rational(1), Rational(4), Rational(3)}, "x");

    auto reduction = berlekamp_factor(polynomial, 3);

    EXPECT_TRUE(reduction.factors.empty())
        << "a prime that erases the leading coefficient is unsuitable";
}

TEST(TranscendentalFactorFiniteField, BridgeSkipsUnsuitablePrimeAndReconstructs) {
    for (const auto &polynomial : {
             Polynomial<Rational>(
                 {Rational(4), Rational(5), Rational(1)}, "x"),
             Polynomial<Rational>(
                 {Rational(1), Rational(4), Rational(3)}, "x")}) {
        ComputationContext context;
        auto factored =
            factor_univariate_bridge_checked(polynomial, context);
        ASSERT_TRUE(factored) << factored.error().message;
        ASSERT_EQ(
            factored.value().completeness, Completeness::Complete);
        ASSERT_EQ(factored.value().value.size(), 2U);

        Polynomial<Rational> reconstructed({Rational(1)}, "x");
        for (const auto &factor : factored.value().value) {
            reconstructed = reconstructed * factor;
        }
        EXPECT_TRUE(reconstructed.make_monic() == polynomial.make_monic());
    }
}

TEST(TranscendentalFactorFiniteField, BerlekampQuadraticIrreducible) {
    Polynomial<Rational> poly({Rational(1), Rational(1), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 2);

    EXPECT_TRUE((result.prime == 2)) << "should use prime 2";
    EXPECT_TRUE((result.factors.size() >= 1)) << "should return at least 1 factor";
    int total_deg = 0;
    for (const auto &f : result.factors) {
        total_deg += f.degree();
    }
    EXPECT_TRUE((total_deg == 2)) << "total degree of factors should equal original degree";
}

TEST(TranscendentalFactorFiniteField, BerlekampX2Minus1Mod3) {
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 3);

    EXPECT_TRUE((result.prime == 3)) << "should use prime 3";
    EXPECT_TRUE((result.factors.size() == 2)) << "should split into 2 factors";
    int total_deg = 0;
    for (const auto &f : result.factors) {
        total_deg += f.degree();
    }
    EXPECT_TRUE((total_deg == 2)) << "total degree of factors should equal 2";
    for (const auto &f : result.factors) {
        EXPECT_TRUE((f.degree() == 1)) << "each factor should be linear";
    }
}

TEST(TranscendentalFactorFiniteField, BerlekampCubicPoly) {
    Polynomial<Rational> poly({Rational(0), Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 5);

    EXPECT_TRUE((result.prime == 5)) << "should use prime 5";
    EXPECT_TRUE((result.factors.size() == 3)) << "should split into 3 factors";
    int total_deg = 0;
    for (const auto &f : result.factors) {
        total_deg += f.degree();
    }
    EXPECT_TRUE((total_deg == 3)) << "total degree of factors should equal 3";
    for (const auto &f : result.factors) {
        EXPECT_TRUE((f.degree() == 1)) << "each factor should be linear";
    }
}

TEST(TranscendentalFactorFiniteField, BerlekampConstantPoly) {
    Polynomial<Rational> poly({Rational(5)}, "x");

    auto result = berlekamp_factor(poly, 0);

    EXPECT_TRUE((result.factors.empty())) << "constant polynomial should have no factors";
}

TEST(TranscendentalFactorFiniteField, NullSpaceX2Minus1Mod3) {
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 3);

    EXPECT_TRUE((result.prime == 3)) << "should use prime 3";
    EXPECT_TRUE((result.null_space_dim == 2)) << "null space dimension should be 2 (two factors)";
    EXPECT_TRUE((result.null_space_basis.size() == 2)) << "should have 2 basis vectors";

    for (const auto &v : result.null_space_basis) {
        EXPECT_TRUE((v.size() == 2)) << "basis vector length should equal polynomial degree";
    }
}

TEST(TranscendentalFactorFiniteField, NullSpaceX2PlusXPlus1Mod2) {
    Polynomial<Rational> poly({Rational(1), Rational(1), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 2);

    EXPECT_TRUE((result.prime == 2)) << "should use prime 2";
    EXPECT_TRUE((result.null_space_dim == 1)) << "null space dimension should be 1 (irreducible)";
    EXPECT_TRUE((result.null_space_basis.size() == 1)) << "should have 1 basis vector";

    if (!result.null_space_basis.empty()) {
        EXPECT_TRUE((result.null_space_basis[0][0] == 1)) << "first basis vector starts with 1";
        EXPECT_TRUE((result.null_space_basis[0][1] == 0)) << "first basis vector second component is 0";
    }
}

TEST(TranscendentalFactorFiniteField, NullSpaceX3MinusXMod5) {
    Polynomial<Rational> poly({Rational(0), Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 5);

    EXPECT_TRUE((result.prime == 5)) << "should use prime 5";
    EXPECT_TRUE((result.null_space_dim == 3)) << "null space dimension should be 3 (three factors)";
    EXPECT_TRUE((result.null_space_basis.size() == 3)) << "should have 3 basis vectors";

    for (const auto &v : result.null_space_basis) {
        EXPECT_TRUE((v.size() == 3)) << "basis vector length should equal polynomial degree";
    }
}

TEST(TranscendentalFactorFiniteField, NullSpaceBasisVectorsInKernel) {
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 3);

    /**
     * @brief 通过零空间维度间接检查 Q*v = v，即 (Q-I)*v = 0。
     * F_3 上的系数为 [2, 0, 1]，其中 -1 mod 3 = 2。
     */
    EXPECT_TRUE((result.null_space_dim == 2)) << "x^2 - 1 mod 3 has null space dim 2";

    for (const auto &v : result.null_space_basis) {
        bool all_zero = true;
        for (int64_t c : v) {
            if (c != 0) {
                all_zero = false;
                break;
            }
        }
        EXPECT_FALSE((all_zero)) << "basis vector should be non-zero";
    }
}

TEST(TranscendentalFactorFiniteField, SplitX2Minus1Mod3) {
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 3);

    EXPECT_TRUE((result.prime == 3)) << "should use prime 3";
    EXPECT_TRUE((result.factors.size() == 2)) << "should produce 2 irreducible factors";

    int total_deg = 0;
    for (const auto &f : result.factors) {
        total_deg += f.degree();
        EXPECT_TRUE((f.degree() == 1)) << "each factor should be linear";
    }
    EXPECT_TRUE((total_deg == 2)) << "total degree should be 2";
}

TEST(TranscendentalFactorFiniteField, SplitX3MinusXMod5) {
    Polynomial<Rational> poly({Rational(0), Rational(-1), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 5);

    EXPECT_TRUE((result.prime == 5)) << "should use prime 5";
    EXPECT_TRUE((result.factors.size() == 3)) << "should produce 3 irreducible factors";

    for (const auto &f : result.factors) {
        EXPECT_TRUE((f.degree() == 1)) << "each factor should be linear";
    }
}

TEST(TranscendentalFactorFiniteField, SplitIrreducibleX2PlusXPlus1Mod2) {
    Polynomial<Rational> poly({Rational(1), Rational(1), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 2);

    EXPECT_TRUE((result.prime == 2)) << "should use prime 2";
    EXPECT_TRUE((result.factors.size() == 1)) << "irreducible polynomial should have 1 factor";
    EXPECT_TRUE((result.factors[0].degree() == 2)) << "factor should have degree 2";
}

TEST(TranscendentalFactorFiniteField, SplitProductVerificationMod3) {
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(1)}, "x");
    int64_t p = 3;

    auto result = berlekamp_factor(poly, p);

    EXPECT_TRUE((result.factors.size() == 2)) << "should have 2 factors";

    /**
     * @brief 按系数验证 (x+a)(x+b) = x^2 + (a+b)x + ab。
     * 原多项式为 x^2 + 2，系数为 [2, 0, 1]；
     * 因此 a+b == 0、a*b == 2 (mod 3)。
     */
    if (result.factors.size() == 2) {
        int64_t a = result.factors[0].coeffs[0].value();
        int64_t b = result.factors[1].coeffs[0].value();
        EXPECT_TRUE(((a + b) % p == 0)) << "sum of constant terms should be 0 mod 3";
        EXPECT_TRUE(((a * b) % p == 2)) << "product of constant terms should be 2 mod 3";
    }
}

TEST(TranscendentalFactorFiniteField, SplitProductVerificationMod5) {
    Polynomial<Rational> poly({Rational(0), Rational(-1), Rational(0), Rational(1)}, "x");
    int64_t p = 5;

    auto result = berlekamp_factor(poly, p);

    EXPECT_TRUE((result.factors.size() == 3)) << "should have 3 factors";

    /**
     * @brief F_5 上 x^3 - x = x(x-1)(x+1) = x(x+4)(x+1)。
     * (x+a)(x+b)(x+c) = x^3 + (a+b+c)x^2 + (ab+ac+bc)x + abc。
     * 原多项式为 x^3 + 4x，系数为 [0, 4, 0, 1]；
     * 因此 a+b+c == 0、ab+ac+bc == 4、abc == 0 (mod 5)。
     */
    if (result.factors.size() == 3) {
        int64_t a = result.factors[0].coeffs[0].value();
        int64_t b = result.factors[1].coeffs[0].value();
        int64_t c = result.factors[2].coeffs[0].value();
        EXPECT_TRUE(((a + b + c) % p == 0)) << "sum of constant terms should be 0 mod 5";
        EXPECT_TRUE(((a * b % p + a * c % p + b * c % p) % p == 4)) << "sum of pairwise products should be 4 mod 5";
        EXPECT_TRUE(((a * b % p * c % p) % p == 0)) << "product of constant terms should be 0 mod 5";
    }
}

TEST(TranscendentalFactorFiniteField, SplitX4Minus1Mod5) {
    /**
     * @brief F_5 上 x^4 - 1 = (x-1)(x+1)(x^2+1)。
     * x^2+1 的根为 2 和 3（4+1 == 9+1 == 0 mod 5），故
     * x^4-1 = (x-1)(x+1)(x-2)(x-3) = (x+4)(x+1)(x+3)(x+2)。
     */
    Polynomial<Rational> poly({Rational(-1), Rational(0), Rational(0), Rational(0), Rational(1)}, "x");

    auto result = berlekamp_factor(poly, 5);

    EXPECT_TRUE((result.prime == 5)) << "should use prime 5";
    EXPECT_TRUE((result.factors.size() == 4)) << "x^4 - 1 should split into 4 linear factors mod 5";

    int total_deg = 0;
    for (const auto &f : result.factors) {
        total_deg += f.degree();
    }
    EXPECT_TRUE((total_deg == 4)) << "total degree should be 4";
}
