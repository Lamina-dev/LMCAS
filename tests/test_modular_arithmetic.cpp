#include "test_common.hpp"
#include "modular_arithmetic.hpp"

#include <limits>
#include <stdexcept>
#include <vector>

using namespace lamina;

static void test_crt_checked_contracts() {
    TEST_CASE("modular_arithmetic checked CRT: explicit errors and context");

    {
        auto result = crt_checked(2, 3, 3, 5);
        EXPECT_TRUE(result.has_value(), "checked CRT succeeds for coprime moduli");
        if (result) {
            EXPECT_TRUE(result.value().first == 8, "checked CRT returns x = 8");
            EXPECT_TRUE(result.value().second == 15, "checked CRT returns modulus 15");
        }
    }

    {
        auto result = multi_crt_checked(std::vector<int64_t>{2, 3, 2},
                                        std::vector<int64_t>{3, 5, 7});
        EXPECT_TRUE(result.has_value(), "checked multi_crt succeeds for coprime moduli");
        if (result) {
            EXPECT_TRUE(result.value().first == 23, "checked multi_crt returns x = 23");
            EXPECT_TRUE(result.value().second == 105, "checked multi_crt returns modulus 105");
        }
    }

    {
        auto result = crt_checked(1, 6, 3, 9);
        EXPECT_TRUE(!result.has_value(), "checked CRT rejects non-coprime moduli");
        EXPECT_TRUE(result.error().code == CasErrc::InvalidArgument,
                    "checked CRT reports InvalidArgument for non-coprime moduli");
    }

    {
        auto result = crt_checked(1, 0, 2, 5);
        EXPECT_TRUE(!result.has_value(), "checked CRT rejects zero modulus");
        EXPECT_TRUE(result.error().code == CasErrc::InvalidArgument,
                    "checked CRT reports InvalidArgument for zero modulus");
    }

    {
        const int64_t too_large = std::numeric_limits<int64_t>::max() / 2 + 1;
        auto result = crt_checked(1, too_large, 2, 3);
        EXPECT_TRUE(!result.has_value(), "checked CRT rejects product overflow");
        EXPECT_TRUE(result.error().code == CasErrc::ResourceLimit,
                    "checked CRT reports ResourceLimit for modulus product overflow");
    }

    {
        auto result = multi_crt_checked(std::vector<int64_t>{1},
                                        std::vector<int64_t>{});
        EXPECT_TRUE(!result.has_value(), "checked multi_crt rejects size mismatch");
        EXPECT_TRUE(result.error().code == CasErrc::InvalidArgument,
                    "checked multi_crt reports InvalidArgument for size mismatch");
    }

    {
        CancellationToken token;
        token.cancel();
        ComputationContext cancelled_context({}, token);
        auto result = crt_checked(2, 3, 3, 5, cancelled_context);
        EXPECT_TRUE(!result.has_value(), "checked CRT observes cancellation");
        EXPECT_TRUE(result.error().code == CasErrc::Cancelled,
                    "checked CRT reports Cancelled");
    }

    {
        ResourceLimits limits;
        limits.max_steps = 0;
        ComputationContext limited_context(limits);
        auto result = multi_crt_checked(std::vector<int64_t>{2, 3},
                                        std::vector<int64_t>{3, 5},
                                        limited_context);
        EXPECT_TRUE(!result.has_value(), "checked multi_crt observes exhausted step budget");
        EXPECT_TRUE(result.error().code == CasErrc::ResourceLimit,
                    "checked multi_crt reports ResourceLimit");
    }
}

static void test_rational_reconstruction_checked_contracts() {
    TEST_CASE("modular_arithmetic checked rational reconstruction");

    {
        auto result = rational_reconstruction_checked(68, 101);
        EXPECT_TRUE(result.has_value(), "checked rational reconstruction succeeds");
        if (result) {
            const auto [a, b] = result.value();
            EXPECT_TRUE(b != 0, "checked rational reconstruction returns nonzero denominator");
            EXPECT_TRUE(((68 * b - a) % 101 + 101) % 101 == 0,
                        "checked rational reconstruction satisfies modular image");
        }
    }

    {
        auto result = rational_reconstruction_checked(1, 0);
        EXPECT_TRUE(!result.has_value(), "checked rational reconstruction rejects invalid modulus");
        EXPECT_TRUE(result.error().code == CasErrc::InvalidArgument,
                    "checked rational reconstruction reports InvalidArgument");
    }

    {
        auto result = rational_reconstruction_checked(0, 1);
        EXPECT_TRUE(!result.has_value(), "checked rational reconstruction reports unsupported case");
        EXPECT_TRUE(result.error().code == CasErrc::Inconclusive,
                    "checked rational reconstruction returns Inconclusive instead of (0,0)");
    }
}

static void test_legacy_compatibility() {
    TEST_CASE("modular_arithmetic legacy compatibility");

    bool threw = false;
    try {
        (void)crt(1, 6, 3, 9);
    } catch (const std::domain_error&) {
        threw = true;
    }
    EXPECT_TRUE(threw, "legacy CRT preserves domain_error for non-coprime moduli");
}

int main() {
    test_crt_checked_contracts();
    test_rational_reconstruction_checked_contracts();
    test_legacy_compatibility();
    return TEST_REPORT();
}
