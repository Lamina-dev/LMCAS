#include "test_common.hpp"
#include "bigint.hpp"

#include <cstdint>
#include <limits>

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

    return TEST_REPORT();
}
