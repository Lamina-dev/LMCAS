#include "test_common.hpp"
#include "bigint.hpp"
#include "rational.hpp"

using namespace LMCAS;

TEST(BigintCore, BigintStrings) {
    BigInt two(2);
    EXPECT_EQ((two.to_string()), ("2")) << "BigInt(2)";

    BigInt x("123");
    EXPECT_EQ((x.to_string()), ("123")) << "BigInt(\"123\")";
}

TEST(BigintCore, BigintOps) {
    BigInt a(2), b(3);
    BigInt c = a * b;
    EXPECT_TRUE(c == BigInt(6)) << "2 * 3 = 6";

    std::string large = "123456789123456789";
    BigInt lx(large);
    EXPECT_EQ((lx.to_string()), (large)) << "Large BigInt String check";

    BigInt y = lx * lx;

    BigInt one(1);
    BigInt zero = y % one;
    EXPECT_TRUE(zero == BigInt(0)) << "y % 1 == 0";

    BigInt n1("123456789");
    BigInt n2("987654321");
    BigInt rem = n2 % n1;
    EXPECT_TRUE(rem == BigInt(9)) << "987654321 % 123456789 = 9";
}

TEST(BigintCore, RationalLargeMultiplication) {
    BigInt n1("123456789");
    BigInt n2("987654321");

    Rational r1(n1);
    Rational r2(n2);

    Rational r3 = r1 * r2;

    EXPECT_TRUE(r3 == Rational(BigInt("121932631112635269")))
        << "Rational Mult Large";
}

TEST(BigintCore, GcdLogic) {
    BigInt two(2);
    BigInt three(3);
    BigInt twelve(12);
    BigInt eighteen(18);

    EXPECT_TRUE(BigInt::gcd(twelve, eighteen) == BigInt(6))
        << "gcd(12, 18)";

    BigInt c(101);
    BigInt d(103);
    EXPECT_TRUE(BigInt::gcd(c, d) == BigInt(1))
        << "gcd(101, 103) - Coprime";

    BigInt zero(0);
    EXPECT_TRUE(BigInt::gcd(zero, twelve) == BigInt(12));
    EXPECT_TRUE(BigInt::gcd(twelve, zero) == BigInt(12));
    EXPECT_TRUE(BigInt::gcd(zero, zero) == BigInt(0));
}
