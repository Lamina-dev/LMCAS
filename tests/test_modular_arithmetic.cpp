#include "test_common.hpp"
#include "modular_arithmetic.hpp"
#include "internal/exact_factorization_support.hpp"

#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace LMCAS;

TEST(ModularArithmetic, MachineResidueBoundaries) {
    const auto maximum = std::numeric_limits<int64_t>::max();
    const ModInt minus_one(maximum - 1, maximum);
    EXPECT_TRUE(((minus_one + minus_one).value() == maximum - 2)) << "modular addition reduces before signed overflow";
    EXPECT_TRUE(((ModInt(0, maximum) - minus_one).value() == 1)) << "modular subtraction stays in the residue interval";
    EXPECT_TRUE(((minus_one * minus_one).value() == 1 &&
                 minus_one.inverse() == minus_one))
        << "full-word modular products and inverse retain exact residues";
    const auto minimum = std::numeric_limits<int64_t>::min();
    auto merged = crt_checked(minimum, 3, maximum, 5);
    EXPECT_TRUE((merged && merged.value().first == 7 &&
                 merged.value().second == 15))
        << "CRT subtracts extreme signed inputs without overflow";
    auto reconstructed = rational_reconstruction_checked(maximum - 1, maximum);
    EXPECT_TRUE((reconstructed && reconstructed.value().first == -1 &&
                 reconstructed.value().second == 1))
        << "large residue normalization reconstructs minus one exactly";
    const BigInt coefficient = (BigInt(1) << 100) * BigInt(3);
    EXPECT_FALSE((is_good_prime(3, {coefficient}))) << "prime filtering retains an unbounded divisible coefficient";
    EXPECT_TRUE((is_good_prime(5, {coefficient}))) << "prime filtering retains an unbounded nondivisible coefficient";
}

TEST(ModularArithmetic, ExtendedGcdSupportsSignedInputs) {
    struct Case {
        int64_t a;
        int64_t b;
        int64_t expected_gcd;
    };
    const Case cases[] = {
        {-3, 5, 1},
        {-3, -5, 1},
        {0, -5, 5},
        {0, 0, 0},
        {std::numeric_limits<int64_t>::max(),
         std::numeric_limits<int64_t>::max() - 1, 1},
    };

    for (const auto& test : cases) {
        int64_t s = 0;
        int64_t t = 0;
        const int64_t gcd = extended_gcd(test.a, test.b, s, t);
        EXPECT_EQ(gcd, test.expected_gcd)
            << "extended_gcd returns the nonnegative gcd";
        EXPECT_EQ(BigInt(test.a) * BigInt(s) + BigInt(test.b) * BigInt(t),
                  BigInt(gcd))
            << "signed Bezout coefficients satisfy the identity without a narrowing oracle";
    }
}

TEST(ModularArithmetic, ExtendedGcdRejectsUnrepresentableMagnitude) {
    const auto minimum = std::numeric_limits<int64_t>::min();
    for (const auto& inputs :
         {std::pair<int64_t, int64_t>{minimum, 1},
          std::pair<int64_t, int64_t>{1, minimum},
          std::pair<int64_t, int64_t>{minimum, minimum}}) {
        int64_t s = 17;
        int64_t t = -23;
        EXPECT_THROW((void)extended_gcd(inputs.first, inputs.second, s, t),
                     std::overflow_error);
        EXPECT_EQ(s, 17) << "failed conversion leaves the first output untouched";
        EXPECT_EQ(t, -23) << "failed conversion leaves the second output untouched";
    }
}

TEST(ModularArithmetic, SymmetricRemainderBoundaries) {
    EXPECT_TRUE(detail::symmetric_remainder(BigInt(5), BigInt(10)) ==
                BigInt(5));
    EXPECT_TRUE(detail::symmetric_remainder(BigInt(-5), BigInt(10)) ==
                BigInt(5));
    EXPECT_TRUE(detail::symmetric_remainder(BigInt(6), BigInt(10)) ==
                BigInt(-4));
    EXPECT_TRUE(detail::symmetric_remainder(BigInt(3), BigInt(5)) ==
                BigInt(-2));
    EXPECT_TRUE(detail::symmetric_remainder(BigInt(-7), BigInt(0)) ==
                BigInt(-7));
}

TEST(ModularArithmetic, CrtCheckedContracts) {
    {
        auto result = crt_checked(2, 3, 3, 5);
        ASSERT_TRUE((result.has_value())) << "checked CRT succeeds for coprime moduli";
        if (result) {
            EXPECT_TRUE((result.value().first == 8)) << "checked CRT returns x = 8";
            EXPECT_TRUE((result.value().second == 15)) << "checked CRT returns modulus 15";
        }
    }

    {
        auto result = multi_crt_checked(std::vector<int64_t>{2, 3, 2},
                                        std::vector<int64_t>{3, 5, 7});
        ASSERT_TRUE((result.has_value())) << "checked multi_crt succeeds for coprime moduli";
        if (result) {
            EXPECT_TRUE((result.value().first == 23)) << "checked multi_crt returns x = 23";
            EXPECT_TRUE((result.value().second == 105)) << "checked multi_crt returns modulus 105";
        }
    }

    {
        auto result = crt_checked(1, 6, 3, 9);
        EXPECT_TRUE((!result.has_value())) << "checked CRT rejects non-coprime moduli";
        EXPECT_TRUE((result.error().code == CasErrc::InvalidArgument)) << "checked CRT reports InvalidArgument for non-coprime moduli";
    }

    {
        auto result = crt_checked(1, 0, 2, 5);
        EXPECT_TRUE((!result.has_value())) << "checked CRT rejects zero modulus";
        EXPECT_TRUE((result.error().code == CasErrc::InvalidArgument)) << "checked CRT reports InvalidArgument for zero modulus";
    }

    {
        const int64_t too_large = std::numeric_limits<int64_t>::max() / 2 + 1;
        auto result = crt_checked(1, too_large, 2, 3);
        EXPECT_TRUE((!result.has_value())) << "checked CRT rejects product overflow";
        EXPECT_TRUE((result.error().code == CasErrc::ResourceLimit)) << "checked CRT reports ResourceLimit for modulus product overflow";
    }

    {
        auto result = multi_crt_checked(std::vector<int64_t>{1},
                                        std::vector<int64_t>{});
        EXPECT_TRUE((!result.has_value())) << "checked multi_crt rejects size mismatch";
        EXPECT_TRUE((result.error().code == CasErrc::InvalidArgument)) << "checked multi_crt reports InvalidArgument for size mismatch";
    }
}

TEST(ModularArithmetic, CrtResourceErrors) {
    {
        CancellationToken token;
        token.cancel();
        ComputationContext cancelled_context({}, token);
        auto result = crt_checked(2, 3, 3, 5, cancelled_context);
        EXPECT_TRUE((!result.has_value())) << "checked CRT observes cancellation";
        EXPECT_TRUE((result.error().code == CasErrc::Cancelled)) << "checked CRT reports Cancelled";
    }

    {
        ResourceLimits limits;
        limits.max_steps = 0;
        ComputationContext limited_context(limits);
        auto result = multi_crt_checked(std::vector<int64_t>{2, 3},
                                        std::vector<int64_t>{3, 5},
                                        limited_context);
        EXPECT_TRUE((!result.has_value())) << "checked multi_crt observes exhausted step budget";
        EXPECT_TRUE((result.error().code == CasErrc::ResourceLimit)) << "checked multi_crt reports ResourceLimit";
    }
}

TEST(ModularArithmetic, RationalReconstructionCheckedContracts) {
    {
        auto result = rational_reconstruction_checked(68, 101);
        ASSERT_TRUE((result.has_value())) << "checked rational reconstruction succeeds";
        if (result) {
            const auto [a, b] = result.value();
            EXPECT_TRUE((b != 0)) << "checked rational reconstruction returns nonzero denominator";
            EXPECT_TRUE((((68 * b - a) % 101 + 101) % 101 == 0)) << "checked rational reconstruction satisfies modular image";
        }
    }

    {
        auto result = rational_reconstruction_checked(1, 0);
        EXPECT_TRUE((!result.has_value())) << "checked rational reconstruction rejects invalid modulus";
        EXPECT_TRUE((result.error().code == CasErrc::InvalidArgument)) << "checked rational reconstruction reports InvalidArgument";
    }

    {
        auto result = rational_reconstruction_checked(0, 1);
        EXPECT_TRUE((!result.has_value())) << "checked rational reconstruction reports unsupported case";
        EXPECT_TRUE((result.error().code == CasErrc::Inconclusive)) << "checked rational reconstruction returns Inconclusive instead of (0,0)";
    }
}

TEST(ModularArithmetic, LegacyCompatibility) {
    bool threw = false;
    try {
        (void)crt(1, 6, 3, 9);
    } catch (const std::domain_error &) {
        threw = true;
    }
    EXPECT_TRUE((threw)) << "legacy CRT preserves domain_error for non-coprime moduli";
}
