#include "test_common.hpp"
#include "bigint.hpp"

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

using namespace LMCAS;

int main() {
    TEST_CASE("BigInt exact narrowing");
    auto zero = BigInt(0).try_to_int64();
    EXPECT_TRUE(zero && *zero == 0, "zero converts exactly");

    auto maximum = BigInt("9223372036854775807").try_to_int64();
    EXPECT_TRUE(maximum && *maximum == std::numeric_limits<std::int64_t>::max(),
                "INT64_MAX converts exactly");

    auto minimum = BigInt("-9223372036854775808").try_to_int64();
    EXPECT_TRUE(minimum && *minimum == std::numeric_limits<std::int64_t>::min(),
                "INT64_MIN converts exactly");

    EXPECT_FALSE(BigInt("9223372036854775808").try_to_int64().has_value(),
                 "positive overflow is reported");
    EXPECT_FALSE(BigInt("-9223372036854775809").try_to_int64().has_value(),
                 "negative overflow is reported");

    TEST_CASE("BigInt exponent zero");
    EXPECT_EQ_STR(BigInt(0).power(0).to_string(), "1", "0^0 follows integer power convention");
    EXPECT_EQ_STR(BigInt(123).power(0).to_string(), "1", "n^0 equals one");

    TEST_CASE("BigInt integer square root");
    EXPECT_EQ_STR(BigInt("2000000000000").sqrt().to_string(), "1414213",
                  "floor(sqrt(2e12))");
    EXPECT_EQ_STR(BigInt("15241578750190521").sqrt().to_string(), "123456789",
                  "large perfect square");

    TEST_CASE("BigInt right shift canonicalizes zero");
    BigInt shifted_one(1);
    shifted_one >>= 64;
    const BigInt canonical_zero(0);
    EXPECT_EQ_STR(shifted_one.to_string(), "0",
                  "exact-limb shift stringifies as zero");
    EXPECT_TRUE(shifted_one == canonical_zero,
                "shifted zero equals constructed zero");
    EXPECT_FALSE(shifted_one < canonical_zero,
                 "shifted zero is not less than zero");
    EXPECT_FALSE(canonical_zero < shifted_one,
                 "shifted zero is not greater than zero");
    EXPECT_TRUE(shifted_one.hash() == canonical_zero.hash(),
                "canonical zeros hash identically");

    BigInt multi_limb("18446744073709551616");
    multi_limb >>= 128;
    EXPECT_TRUE(multi_limb == canonical_zero,
                "multi-limb overshift canonicalizes zero");
    BigInt negative(-1);
    negative >>= 64;
    EXPECT_TRUE(negative == canonical_zero,
                "negative overshift canonicalizes zero");


    TEST_CASE("BigInt modular power validates signs and canonicalizes residues");
    EXPECT_EQ_STR(BigInt::pow_mod(BigInt(-2), BigInt(3), BigInt(5)).to_string(),
                  "2", "negative bases produce canonical positive residues");
    EXPECT_EQ_STR(
        BigInt::pow_mod(BigInt("-18446744073709551618"), BigInt(3),
                        BigInt("18446744073709551617")).to_string(),
        "18446744073709551616",
        "multi-limb modular power canonicalizes a negative base");
    EXPECT_EQ_STR(BigInt::pow_mod(BigInt(99), BigInt(0), BigInt(1)).to_string(),
                  "0", "modulus one always yields zero");

    bool negative_exponent_rejected = false;
    try {
        (void)BigInt::pow_mod(BigInt(2), BigInt(-1), BigInt(5));
    } catch (const std::domain_error&) {
        negative_exponent_rejected = true;
    }
    EXPECT_TRUE(negative_exponent_rejected,
                "negative modular exponents are rejected");

    bool nonpositive_modulus_rejected = false;
    try {
        (void)BigInt::pow_mod(BigInt(2), BigInt(3), BigInt(-5));
    } catch (const std::domain_error&) {
        nonpositive_modulus_rejected = true;
    }
    EXPECT_TRUE(nonpositive_modulus_rejected,
                "nonpositive modular moduli are rejected");

    TEST_CASE("BigInt moved-from storage is reusable");
    BigInt source("18446744073709551616");
    BigInt destination;
    destination = std::move(source);
    source = BigInt(7);
    EXPECT_EQ_STR(destination.to_string(), "18446744073709551616",
                  "move assignment preserves destination");
    EXPECT_EQ_STR(source.to_string(), "7",
                  "moved-from value can allocate storage again");
    return TEST_REPORT();
}
