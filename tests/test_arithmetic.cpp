#include "test_common.hpp"
#include "bigint.hpp"
#include "rational.hpp"

using namespace LMCAS;

TEST(Arithmetic, Bigint) {
    BigInt a(1);
    BigInt b(5);
    BigInt c(6);

    BigInt b2 = b * b;

    EXPECT_EQ((b2.to_string()), ("25")) << "5^2 = 25";

    BigInt ac4 = BigInt(4) * a * c;
    EXPECT_EQ((ac4.to_string()), ("24")) << "4*1*6 = 24";

    BigInt D = b2 - ac4;
    EXPECT_EQ((D.to_string()), ("1")) << "25 - 24 = 1";

    EXPECT_TRUE((D.is_perfect_square())) << "1 is perfect square";
    if (D.is_perfect_square()) {
        EXPECT_EQ((D.sqrt().to_string()), ("1")) << "sqrt(1) = 1";
    }
}

TEST(Arithmetic, Rational) {
    Rational a(1);
    Rational b(5);
    Rational c(6);

    Rational b2 = b.power(BigInt(2));
    EXPECT_EQ((b2.to_string()), ("25")) << "rat 5^2 = 25";

    Rational ac4 = Rational(4) * a * c;
    EXPECT_EQ((ac4.to_string()), ("24")) << "rat 4*1*6 = 24";

    Rational D = b2 - ac4;
    EXPECT_EQ((D.to_string()), ("1")) << "rat 25 - 24 = 1";
}

TEST(Arithmetic, LargeMul) {
    std::string s = "123456789123456789";
    BigInt a(s);
    BigInt b(s);

    BigInt c = a * b;

    std::string expected = "15241578780673678515622620750190521";
    EXPECT_EQ((c.to_string()), (expected)) << "Large Mul Check";
}
