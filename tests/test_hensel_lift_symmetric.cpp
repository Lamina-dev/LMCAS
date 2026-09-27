#include "test_hensel_lift_support.hpp"

TEST(HenselLiftSymmetric, SymmetricReprBasicRange) {
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(-1), BigInt(0), BigInt(1)}, "x");

    int64_t p = 3;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(1, p), ModInt(1, p)}; /**< 因子 x+1。 */
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(2, p), ModInt(1, p)}; /**< 模 3 因子 x+2 = x-1。 */
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2};

    int lift_bound = 4; /**< 提升至模 3^4 = 81。 */
    auto checked =
        hensel_lift_checked(poly, mod_factors, p, lift_bound);
    ASSERT_TRUE(checked) << checked.error().message;
    auto lifted = std::move(checked.value());

    EXPECT_TRUE((lifted.size() == 2)) << "should produce 2 lifted factors";

    BigInt mod81(81);
    BigInt half(40);

    bool all_in_range = true;
    for (const auto &factor : lifted) {
        for (const auto &c : factor.coeffs) {
            BigInt abs_c = c.abs();
            if (abs_c > half) {
                all_in_range = false;
                break;
            }
        }
        if (!all_in_range) {
            break;
        }
    }
    EXPECT_TRUE((all_in_range)) << "all coefficients in [-40, 40] for mod 81";

    std::vector<BigInt> product = test_poly_mul(lifted[0].coeffs, lifted[1].coeffs);
    auto reduced_product = test_reduce(product, mod81);
    auto reduced_f = test_reduce(poly.coeffs, mod81);

    size_t n = std::max(reduced_f.size(), reduced_product.size());
    reduced_f.resize(n, BigInt(0));
    reduced_product.resize(n, BigInt(0));

    bool match = true;
    for (size_t i = 0; i < n; ++i) {
        if (reduced_f[i] != reduced_product[i]) {
            match = false;
            break;
        }
    }
    EXPECT_TRUE((match)) << "product of lifted factors = f (mod 81)";
}

TEST(HenselLiftSymmetric, SymmetricReprLargeCoefficients) {
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(105), BigInt(71), BigInt(15), BigInt(1)}, "x");

    int64_t p = 7;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(0, p), ModInt(1, p)}; /**< 因子 x。 */
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(3, p), ModInt(1, p)}; /**< 因子 x+3。 */
    Polynomial<ModInt> f3("x");
    f3.coeffs = {ModInt(5, p), ModInt(1, p)}; /**< 因子 x+5。 */
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2, f3};

    int lift_bound = 3; /**< 提升至模 7^3 = 343。 */
    auto checked =
        hensel_lift_checked(poly, mod_factors, p, lift_bound);
    ASSERT_TRUE(checked) << checked.error().message;
    auto lifted = std::move(checked.value());

    EXPECT_TRUE((lifted.size() == 3)) << "should produce 3 lifted factors";

    BigInt mod343(343);
    BigInt half(171);

    bool all_in_range = true;
    for (const auto &factor : lifted) {
        for (const auto &c : factor.coeffs) {
            BigInt abs_c = c.abs();
            if (abs_c > half) {
                all_in_range = false;
                break;
            }
        }
        if (!all_in_range) {
            break;
        }
    }
    EXPECT_TRUE((all_in_range)) << "all coefficients in [-171, 171] for mod 343";

    std::vector<BigInt> product = lifted[0].coeffs;
    for (size_t i = 1; i < lifted.size(); ++i) {
        product = test_poly_mul(product, lifted[i].coeffs);
        for (auto &c : product) {
            c = test_sym_mod(c, mod343);
        }
        while (!product.empty() && product.back().is_zero()) {
            product.pop_back();
        }
    }

    auto reduced_f = test_reduce(poly.coeffs, mod343);
    auto reduced_product = test_reduce(product, mod343);

    size_t n2 = std::max(reduced_f.size(), reduced_product.size());
    reduced_f.resize(n2, BigInt(0));
    reduced_product.resize(n2, BigInt(0));

    bool match = true;
    for (size_t i = 0; i < n2; ++i) {
        if (reduced_f[i] != reduced_product[i]) {
            match = false;
            break;
        }
    }
    EXPECT_TRUE((match)) << "product of lifted factors = f (mod 343)";
}

TEST(HenselLiftSymmetric, SymmetricReprRoundtrip) {
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(21), BigInt(10), BigInt(1)}, "x");

    int64_t p = 5;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(3, p), ModInt(1, p)}; /**< 因子 x+3。 */
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(2, p), ModInt(1, p)}; /**< 模 5 因子 x+2 = x+7。 */
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2};

    int lift_bound = 3; /**< 提升至模 5^3 = 125。 */
    auto checked =
        hensel_lift_checked(poly, mod_factors, p, lift_bound);
    ASSERT_TRUE(checked) << checked.error().message;
    auto lifted = std::move(checked.value());

    EXPECT_TRUE((lifted.size() == 2)) << "should produce 2 lifted factors";

    BigInt mod125(125);
    BigInt half(62);

    bool all_in_range = true;
    for (const auto &factor : lifted) {
        for (const auto &c : factor.coeffs) {
            BigInt abs_c = c.abs();
            if (abs_c > half) {
                all_in_range = false;
                break;
            }
        }
        if (!all_in_range) {
            break;
        }
    }
    EXPECT_TRUE((all_in_range)) << "all coefficients in [-62, 62] for mod 125";

    bool roundtrip_ok = true;
    for (const auto &factor : lifted) {
        for (const auto &c : factor.coeffs) {
            BigInt reduced = test_sym_mod(c, mod125);
            if (reduced != c) {
                roundtrip_ok = false;
                break;
            }
        }
        if (!roundtrip_ok) {
            break;
        }
    }
    EXPECT_TRUE((roundtrip_ok)) << "all coefficients are already in symmetric form (idempotent)";
}

TEST(HenselLiftSymmetric, SymmetricReprBoundary) {
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(6), BigInt(5), BigInt(1)}, "x");

    int64_t p = 2;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(0, p), ModInt(1, p)}; /**< 因子 x。 */
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(1, p), ModInt(1, p)}; /**< 因子 x+1。 */
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2};

    int lift_bound = 3; /**< 提升至模 2^3 = 8。 */
    auto checked =
        hensel_lift_checked(poly, mod_factors, p, lift_bound);
    ASSERT_TRUE(checked) << checked.error().message;
    auto lifted = std::move(checked.value());

    EXPECT_TRUE((lifted.size() == 2)) << "should produce 2 lifted factors";

    BigInt mod8(8);
    BigInt half(4);

    bool all_in_range = true;
    for (const auto &factor : lifted) {
        for (const auto &c : factor.coeffs) {
            BigInt abs_c = c.abs();
            if (abs_c > half) {
                all_in_range = false;
                break;
            }
        }
        if (!all_in_range) {
            break;
        }
    }
    EXPECT_TRUE((all_in_range)) << "all coefficients in [-4, 4] for mod 8";

    std::vector<BigInt> product = test_poly_mul(lifted[0].coeffs, lifted[1].coeffs);
    auto reduced_product = test_reduce(product, mod8);
    auto reduced_f = test_reduce(poly.coeffs, mod8);

    size_t n = std::max(reduced_f.size(), reduced_product.size());
    reduced_f.resize(n, BigInt(0));
    reduced_product.resize(n, BigInt(0));

    bool match = true;
    for (size_t i = 0; i < n; ++i) {
        if (reduced_f[i] != reduced_product[i]) {
            match = false;
            break;
        }
    }
    EXPECT_TRUE((match)) << "product of lifted factors = f (mod 8)";
}

TEST(HenselLiftSymmetric, SymmetricReprHighLift) {
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(2), BigInt(3), BigInt(1)}, "x");

    int64_t p = 5;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(1, p), ModInt(1, p)}; /**< 因子 x+1。 */
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(2, p), ModInt(1, p)}; /**< 因子 x+2。 */
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2};

    int lift_bound = 5; /**< 提升至模 5^5 = 3125。 */
    auto checked =
        hensel_lift_checked(poly, mod_factors, p, lift_bound);
    ASSERT_TRUE(checked) << checked.error().message;
    auto lifted = std::move(checked.value());

    EXPECT_TRUE((lifted.size() == 2)) << "should produce 2 lifted factors";

    BigInt mod3125(3125);
    BigInt half(1562);

    bool all_in_range = true;
    for (const auto &factor : lifted) {
        for (const auto &c : factor.coeffs) {
            BigInt abs_c = c.abs();
            if (abs_c > half) {
                all_in_range = false;
                break;
            }
        }
        if (!all_in_range) {
            break;
        }
    }
    EXPECT_TRUE((all_in_range)) << "all coefficients in [-1562, 1562] for mod 3125";

    std::vector<BigInt> product = test_poly_mul(lifted[0].coeffs, lifted[1].coeffs);
    auto reduced_product = test_reduce(product, mod3125);
    auto reduced_f = test_reduce(poly.coeffs, mod3125);

    size_t n = std::max(reduced_f.size(), reduced_product.size());
    reduced_f.resize(n, BigInt(0));
    reduced_product.resize(n, BigInt(0));

    bool match = true;
    for (size_t i = 0; i < n; ++i) {
        if (reduced_f[i] != reduced_product[i]) {
            match = false;
            break;
        }
    }
    EXPECT_TRUE((match)) << "product of lifted factors = f (mod 3125)";
}
