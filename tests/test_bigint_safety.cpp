#include "test_common.hpp"
#include "bigint.hpp"
#include "computation_context.hpp"

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

using namespace LMCAS;

TEST(BigintSafety, WordImportBoundaries) {
    const auto signed_minimum = std::numeric_limits<std::int64_t>::min();
    const auto unsigned_maximum = std::numeric_limits<std::uint64_t>::max();
    EXPECT_TRUE(
        BigInt(signed_minimum) == BigInt("-9223372036854775808"))
        << "signed minimum imports without signed negation";
    const BigInt imported(unsigned_maximum);
    EXPECT_TRUE(imported == BigInt("18446744073709551615"))
        << "unsigned import preserves all 64 bits";
    const auto narrowed = imported.try_to_uint64();
    EXPECT_TRUE((narrowed && *narrowed == unsigned_maximum)) << "unsigned maximum round-trips exactly";
    EXPECT_FALSE(((imported + BigInt(1)).try_to_uint64().has_value())) << "unsigned narrowing rejects the first out-of-range integer";
    EXPECT_FALSE((BigInt(-1).try_to_uint64().has_value())) << "unsigned narrowing rejects negative integers";
    EXPECT_TRUE((imported.bit_length() == 64 &&
                 (imported + BigInt(1)).bit_length() == 65))
        << "bit length crosses the limb boundary exactly";
}

static void test_combinatorial_domain_errors() {
    auto negative = BigInt::factorial_checked(BigInt(-1));
    EXPECT_TRUE((!negative && negative.error().code == CasErrc::DomainError)) << "negative factorial arguments are domain errors";
    auto oversized = BigInt::factorial_checked(BigInt("4294967296"));
    EXPECT_TRUE((!oversized && oversized.error().code == CasErrc::ResourceLimit)) << "kernel arguments are checked before narrowing";
}

static void test_combinatorial_resource_limits() {
    ResourceLimits limits;
    limits.max_integer_bits = 32;
    ComputationContext bounded(limits);
    auto exhausted = BigInt::factorial_checked(BigInt(21), bounded);
    EXPECT_TRUE((!exhausted && exhausted.error().code == CasErrc::ResourceLimit)) << "integer bit budget prevents oversized combinatorial results";
    limits.max_steps = 0;
    ComputationContext no_steps(limits);
    auto no_work = BigInt::ncr_checked(BigInt(10), BigInt(2), no_steps);
    EXPECT_TRUE((!no_work && no_work.error().code == CasErrc::ResourceLimit)) << "checked combinatorics honor the step budget";
}

TEST(BigintSafety, CheckedCombinatorics) {
    auto factorial = BigInt::factorial_checked(BigInt(21));
    EXPECT_TRUE((factorial && factorial.value() == BigInt("51090942171709440000"))) << "21 factorial is exact beyond uint64";
    auto combinations = BigInt::ncr_checked(BigInt(100), BigInt(50));
    EXPECT_TRUE((combinations &&
                 combinations.value() == BigInt("100891344545564193334812497256")))
        << "binomial kernel retains its full arbitrary precision result";
    auto permutations = BigInt::npr_checked(BigInt(25), BigInt(24));
    EXPECT_TRUE((permutations &&
                 permutations.value() == BigInt("15511210043330985984000000")))
        << "permutation kernel retains its full arbitrary precision result";
    test_combinatorial_domain_errors();
    test_combinatorial_resource_limits();
}

TEST(BigintSafety, SignedMachineBoundaries) {
    const long long_min = std::numeric_limits<long>::min();
    const long long_max = std::numeric_limits<long>::max();
    EXPECT_TRUE(BigInt(long_min) == BigInt(std::to_string(long_min)));
    EXPECT_TRUE(BigInt(long_max) == BigInt(std::to_string(long_max)));
}

TEST(BigintSafety, ExactNarrowing) {
    auto zero = BigInt(0).try_to_int64();
    EXPECT_TRUE((zero && *zero == 0)) << "zero converts exactly";

    auto maximum = BigInt("9223372036854775807").try_to_int64();
    EXPECT_TRUE((maximum && *maximum == std::numeric_limits<std::int64_t>::max())) << "INT64_MAX converts exactly";

    auto minimum = BigInt("-9223372036854775808").try_to_int64();
    EXPECT_TRUE((minimum && *minimum == std::numeric_limits<std::int64_t>::min())) << "INT64_MIN converts exactly";

    EXPECT_FALSE((BigInt("9223372036854775808").try_to_int64().has_value())) << "positive overflow is reported";
    EXPECT_FALSE((BigInt("-9223372036854775809").try_to_int64().has_value())) << "negative overflow is reported";
}

TEST(BigintSafety, IntegerPowersAndRoots) {
    EXPECT_TRUE(BigInt(0).power(0) == BigInt(1))
        << "0^0 follows integer power convention";
    EXPECT_TRUE(BigInt(123).power(0) == BigInt(1))
        << "n^0 equals one";
    EXPECT_TRUE(BigInt("2000000000000").sqrt() == BigInt("1414213"))
        << "floor(sqrt(2e12))";
    EXPECT_TRUE(
        BigInt("15241578750190521").sqrt() == BigInt("123456789"))
        << "large perfect square";
}

TEST(BigintSafety, ShiftedZero) {
    BigInt shifted_one(1);
    shifted_one >>= 64;
    const BigInt canonical_zero(0);
    EXPECT_EQ((shifted_one.to_string()), ("0")) << "exact-limb shift stringifies as zero";
    EXPECT_TRUE((shifted_one == canonical_zero)) << "shifted zero equals constructed zero";
    EXPECT_FALSE((shifted_one < canonical_zero)) << "shifted zero is not less than zero";
    EXPECT_FALSE((canonical_zero < shifted_one)) << "shifted zero is not greater than zero";
    EXPECT_TRUE((shifted_one.hash() == canonical_zero.hash())) << "canonical zeros hash identically";

    BigInt multi_limb("18446744073709551616");
    multi_limb >>= 128;
    EXPECT_TRUE((multi_limb == canonical_zero)) << "multi-limb overshift canonicalizes zero";
    BigInt negative(-1);
    negative >>= 64;
    EXPECT_TRUE((negative == canonical_zero)) << "negative overshift canonicalizes zero";
}

TEST(BigintSafety, ModularPowerDomain) {
    EXPECT_TRUE(
        BigInt::pow_mod(BigInt(-2), BigInt(3), BigInt(5)) == BigInt(2))
        << "negative bases produce canonical positive residues";
    EXPECT_TRUE(
        BigInt::pow_mod(
            BigInt("-18446744073709551618"), BigInt(3),
            BigInt("18446744073709551617")) ==
        BigInt("18446744073709551616"))
        << "multi-limb modular power canonicalizes a negative base";
    EXPECT_TRUE(
        BigInt::pow_mod(BigInt(99), BigInt(0), BigInt(1)) == BigInt(0))
        << "modulus one always yields zero";

    bool negative_exponent_rejected = false;
    try {
        (void)BigInt::pow_mod(BigInt(2), BigInt(-1), BigInt(5));
    } catch (const std::domain_error &) {
        negative_exponent_rejected = true;
    }
    EXPECT_TRUE((negative_exponent_rejected)) << "negative modular exponents are rejected";

    bool nonpositive_modulus_rejected = false;
    try {
        (void)BigInt::pow_mod(BigInt(2), BigInt(3), BigInt(-5));
    } catch (const std::domain_error &) {
        nonpositive_modulus_rejected = true;
    }
    EXPECT_TRUE((nonpositive_modulus_rejected)) << "nonpositive modular moduli are rejected";
}

TEST(BigintSafety, MovedValueReuse) {
    BigInt source("18446744073709551616");
    BigInt destination;
    destination = std::move(source);
    source = BigInt(7);
    EXPECT_TRUE(destination == BigInt("18446744073709551616"))
        << "move assignment preserves destination";
    EXPECT_TRUE(source == BigInt(7))
        << "moved-from value can allocate storage again";
}
