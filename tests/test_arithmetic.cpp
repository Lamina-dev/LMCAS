#include "bigint.hpp"
#include "rational.hpp"
#include <iostream>

void test_bigint() {
    std::cout << "Enter test_bigint" << std::endl;
    std::cout << "Testing BigInt..." << std::endl;
    BigInt a(1);
    std::cout << "BigInt a(1) created" << std::endl;
    BigInt b(5);
    BigInt c(6);

    BigInt b2 = b * b;
    std::cout << "b^2 = " << b2.to_string() << std::endl;

    BigInt ac4 = BigInt(4) * a * c;
    std::cout << "4ac = " << ac4.to_string() << std::endl;

    BigInt D = b2 - ac4;
    std::cout << "D = " << D.to_string() << std::endl;

    if (D.is_perfect_square()) {
        std::cout << "D is square" << std::endl;
        std::cout << "sqrt(D) = " << D.sqrt().to_string() << std::endl;
    } else {
        std::cout << "D is not square" << std::endl;
    }
}

void test_rational() {
    std::cout << "Testing Rational..." << std::endl;
    Rational a(1);
    Rational b(5);
    Rational c(6);

    Rational b2 = b.power(BigInt(2));
    std::cout << "rat b^2 = " << b2.to_string() << std::endl;

    Rational ac4 = Rational(4) * a * c;
    std::cout << "rat 4ac = " << ac4.to_string() << std::endl;

    Rational D = b2 - ac4;
    std::cout << "rat D = " << D.to_string() << std::endl;
}

void test_large_mul() {
    std::cout << "Testing Large Multiplication..." << std::endl;
    std::cout << "sizeof(mp_limb_t) = " << sizeof(mp_limb_t) << std::endl;
    // std::cout << "sizeof(lamp_ui) = " << sizeof(lamp_ui) << std::endl;
    std::string s = "123456789123456789";
    BigInt a(s);
    BigInt b(s);
    std::cout << "a = " << a.to_string() << std::endl;
    BigInt c = a * b;
    std::cout << "a*a = " << c.to_string() << std::endl;
    
    // Check known result
    // 123456789123456789 ^ 2 = 15241578780673678515622620750190521
    std::string expected = "15241578780673678515622620750190521";
    if (c.to_string() == expected) {
        std::cout << "[PASS] Large Mul" << std::endl;
    } else {
        std::cout << "[FAIL] Large Mul" << std::endl;
        std::cout << "Expected: " << expected << std::endl;
        std::cout << "Got:      " << c.to_string() << std::endl;
    }
}

int main() {
    std::cout << "Starting test_arithmetic..." << std::endl;
    try {
        test_bigint();
        test_rational();
        test_large_mul();
    } catch (const std::exception& e) {
        std::cout << "Exception: " << e.what() << std::endl;
    } catch (...) {
        std::cout << "Unknown Exception!" << std::endl;
    }
    return 0;
}
