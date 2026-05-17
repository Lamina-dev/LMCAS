#include "test_common.hpp"
#include <iostream>

int main() {
    
    // Test 1: 2^128 - 1 + 1
    BigInt a1 = BigInt(2).power(128) - BigInt(1);
    BigInt b1 = BigInt(1);
    BigInt c1 = a1 + b1;
    EXPECT_TRUE(c1 == BigInt(2).power(128), "2^128 - 1 + 1 should equal 2^128");

    // Test 2: 2^192 - 1 + 1
    BigInt a2 = BigInt(2).power(192) - BigInt(1);
    BigInt b2 = BigInt(1);
    BigInt c2 = a2 + b2;
    EXPECT_TRUE(c2 == BigInt(2).power(192), "2^192 - 1 + 1 should equal 2^192");

    // Test 3: Large carry propagation
    BigInt a3 = BigInt(2).power(256) - BigInt(1);
    BigInt b3 = BigInt(2).power(64) - BigInt(1);
    BigInt c3 = a3 + b3;
    EXPECT_TRUE(c3 == BigInt(2).power(256) + BigInt(2).power(64) - BigInt(2), "Carry propagation across multiple limbs");

    // Test 4: b is larger than a
    BigInt a4 = BigInt(1);
    BigInt b4 = BigInt(2).power(192) - BigInt(1);
    BigInt c4 = a4 + b4;
    EXPECT_TRUE(c4 == BigInt(2).power(192), "1 + 2^192 - 1 should equal 2^192");

    return TEST_REPORT();
}