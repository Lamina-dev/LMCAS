#include "test_common.hpp"
#include "bigint.hpp"
#include "rational.hpp"

// Note: symbolic.hpp (included by test_common.hpp) likely includes bigint/rational or they are in scope.
// If not, we include them above.

void test_bigint() {
    TEST_CASE("BigInt Basic");
    BigInt a(1);
    BigInt b(5);
    BigInt c(6);

    BigInt b2 = b * b;
    // std::cout << "b^2 = " << b2.to_string() << std::endl;
    EXPECT_EQ_STR(b2.to_string(), "25", "5^2 = 25");

    BigInt ac4 = BigInt(4) * a * c;
    EXPECT_EQ_STR(ac4.to_string(), "24", "4*1*6 = 24");

    BigInt D = b2 - ac4;
    EXPECT_EQ_STR(D.to_string(), "1", "25 - 24 = 1");

    EXPECT_TRUE(D.is_perfect_square(), "1 is perfect square");
    if (D.is_perfect_square()) {
        EXPECT_EQ_STR(D.sqrt().to_string(), "1", "sqrt(1) = 1");
    }
}

void test_rational() {
    TEST_CASE("Rational Basic");
    Rational a(1);
    Rational b(5);
    Rational c(6);

    Rational b2 = b.power(BigInt(2));
    EXPECT_EQ_STR(b2.to_string(), "25", "rat 5^2 = 25");

    Rational ac4 = Rational(4) * a * c;
    EXPECT_EQ_STR(ac4.to_string(), "24", "rat 4*1*6 = 24");

    Rational D = b2 - ac4;
    EXPECT_EQ_STR(D.to_string(), "1", "rat 25 - 24 = 1");
}

void test_large_mul() {
    TEST_CASE("Large Multiplication");
    // std::cout << "sizeof(mp_limb_t) = " << sizeof(mp_limb_t) << std::endl;
    std::string s = "123456789123456789";
    BigInt a(s);
    BigInt b(s);
    // std::cout << "a = " << a.to_string() << std::endl;
    BigInt c = a * b;
    // std::cout << "a*a = " << c.to_string() << std::endl;
    
    // Check known result
    // 123456789123456789 ^ 2 = 15241578780673678515622620750190521
    std::string expected = "15241578780673678515622620750190521";
    EXPECT_EQ_STR(c.to_string(), expected, "Large Mul Check");
}

int main() {
    try {
        test_bigint();
        test_rational();
        test_large_mul();
    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        g_failures++;
    } catch (...) {
        std::cout << "[FAIL] Unknown Exception!" << std::endl;
        g_failures++;
    }
    return TEST_REPORT();
}
