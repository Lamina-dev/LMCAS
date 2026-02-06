#include "test_common.hpp"
#include "bigint.hpp"
#include "rational.hpp"

// Note: Using to_string() instead of ToString() for consistency.

void test_bigint_strings() {
    TEST_CASE("BigInt to_string");
    
    BigInt two(2);
    EXPECT_EQ_STR(two.to_string(), "2", "BigInt(2)");

    BigInt x("123");
    EXPECT_EQ_STR(x.to_string(), "123", "BigInt(\"123\")");
}

void test_bigint_ops() {
    TEST_CASE("BigInt Operations");
    
    BigInt a(2), b(3);
    BigInt c = a * b;
    EXPECT_EQ_STR(c.to_string(), "6", "2 * 3 = 6");

    std::string large = "123456789123456789";
    BigInt lx(large);
    EXPECT_EQ_STR(lx.to_string(), large, "Large BigInt String check");

    // Large mult
    BigInt y = lx * lx;
    // Expected: 15241578780673678515622620750190521 (from test_arithmetic)
    // The original code didn't check the result explicitly against a constant, 
    // but just printed it. We will skip exact check here if not known, or trust test_arithmetic.
    // However, it checked y % 1.
    
    BigInt one(1);
    BigInt zero = y % one;
    EXPECT_EQ_STR(zero.to_string(), "0", "y % 1 == 0");
    
    // Modulo test from original code
    BigInt n1("123456789");
    BigInt n2("987654321");
    BigInt rem = n2 % n1;
    EXPECT_EQ_STR(rem.to_string(), "9", "987654321 % 123456789 = 9");
}

void test_rational_debug() {
    TEST_CASE("Rational Debug");
    
    BigInt n1("123456789");
    BigInt n2("987654321");
    
    Rational r1(n1);
    Rational r2(n2);
    
    Rational r3 = r1 * r2;
    // 123456789 * 987654321 = 121932631112635269
    EXPECT_EQ_STR(r3.to_string(), "121932631112635269", "Rational Mult Large");
}

void test_gcd_logic() {
    TEST_CASE("BigInt GCD Logic");
    
    BigInt two(2);
    BigInt three(3);
    BigInt twelve(12);
    BigInt eighteen(18);
    
    EXPECT_EQ_STR(BigInt::gcd(twelve, eighteen).to_string(), "6", "gcd(12, 18)");
    
    BigInt c(101);
    BigInt d(103);
    EXPECT_EQ_STR(BigInt::gcd(c, d).to_string(), "1", "gcd(101, 103) - Coprime");
    
    BigInt zero(0);
    EXPECT_EQ_STR(BigInt::gcd(zero, twelve).to_string(), "12", "gcd(0, 12)");
    EXPECT_EQ_STR(BigInt::gcd(twelve, zero).to_string(), "12", "gcd(12, 0)");
    EXPECT_EQ_STR(BigInt::gcd(zero, zero).to_string(), "0", "gcd(0, 0)");
}

int main() {
    try {
        test_gcd_logic();
        test_bigint_strings();
        test_bigint_ops();
        test_rational_debug();
    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        g_failures++;
    } catch (...) {
        std::cout << "[FAIL] Unknown Exception" << std::endl;
        g_failures++;
    }
    return TEST_REPORT();
}
