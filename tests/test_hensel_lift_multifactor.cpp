#include "test_hensel_lift_support.hpp"

TEST(HenselLiftMultifactor, MultiFactorLift3Factors) {
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(0), BigInt(-1), BigInt(0), BigInt(1)}, "x");

    int64_t p = 5;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(0, p), ModInt(1, p)}; /**< 因子 x。 */
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(4, p), ModInt(1, p)}; /**< 模 5 因子 x+4 = x-1。 */
    Polynomial<ModInt> f3("x");
    f3.coeffs = {ModInt(1, p), ModInt(1, p)}; /**< 因子 x+1。 */
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2, f3};

    int lift_bound = 2; /**< 提升至模 5^2 = 25。 */
    auto checked =
        hensel_lift_checked(poly, mod_factors, p, lift_bound);
    ASSERT_TRUE(checked) << checked.error().message;
    auto lifted = std::move(checked.value());

    EXPECT_TRUE((lifted.size() == 3)) << "should produce 3 lifted factors";

    BigInt mod25(25);
    std::vector<BigInt> product = lifted[0].coeffs;
    for (size_t i = 1; i < lifted.size(); ++i) {
        product = test_poly_mul(product, lifted[i].coeffs);
        for (auto &c : product) {
            c = test_sym_mod(c, mod25);
        }
        while (!product.empty() && product.back().is_zero()) {
            product.pop_back();
        }
    }

    auto reduced_f = test_reduce(poly.coeffs, mod25);
    auto reduced_product = test_reduce(product, mod25);

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
    EXPECT_TRUE((match)) << "product of lifted factors = f (mod 25)";
}

TEST(HenselLiftMultifactor, MultiFactorLift2FactorsViaApi) {
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(-1), BigInt(0), BigInt(1)}, "x");

    int64_t p = 3;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(1, p), ModInt(1, p)}; /**< 因子 x+1。 */
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(2, p), ModInt(1, p)}; /**< 模 3 因子 x+2 = x-1。 */
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2};

    int lift_bound = 2; /**< 提升至模 3^2 = 9。 */
    auto checked =
        hensel_lift_checked(poly, mod_factors, p, lift_bound);
    ASSERT_TRUE(checked) << checked.error().message;
    auto lifted = std::move(checked.value());

    EXPECT_TRUE((lifted.size() == 2)) << "should produce 2 lifted factors";

    BigInt mod9(9);
    std::vector<BigInt> product = test_poly_mul(lifted[0].coeffs, lifted[1].coeffs);
    auto reduced_product = test_reduce(product, mod9);
    auto reduced_f = test_reduce(poly.coeffs, mod9);

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
    EXPECT_TRUE((match)) << "product of 2 lifted factors = f (mod 9)";

    bool exact = test_verify_factorization(poly.coeffs, lifted[0].coeffs, lifted[1].coeffs, mod9);
    EXPECT_TRUE((exact)) << "lifted factors are exact for x^2-1";
}

TEST(HenselLiftMultifactor, MultiFactorLiftSymmetricCoeffs) {
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(6), BigInt(11), BigInt(6), BigInt(1)}, "x");

    int64_t p = 5;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(1, p), ModInt(1, p)}; /**< 因子 x+1。 */
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(2, p), ModInt(1, p)}; /**< 因子 x+2。 */
    Polynomial<ModInt> f3("x");
    f3.coeffs = {ModInt(3, p), ModInt(1, p)}; /**< 因子 x+3。 */
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2, f3};

    int lift_bound = 2; /**< 提升至模 25。 */
    auto checked =
        hensel_lift_checked(poly, mod_factors, p, lift_bound);
    ASSERT_TRUE(checked) << checked.error().message;
    auto lifted = std::move(checked.value());

    EXPECT_TRUE((lifted.size() == 3)) << "should produce 3 lifted factors";

    BigInt mod25(25);
    BigInt half(12);
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
    EXPECT_TRUE((all_in_range)) << "all coefficients in symmetric range [-12, 12]";

    std::vector<BigInt> product = lifted[0].coeffs;
    for (size_t i = 1; i < lifted.size(); ++i) {
        product = test_poly_mul(product, lifted[i].coeffs);
        for (auto &c : product) {
            c = test_sym_mod(c, mod25);
        }
        while (!product.empty() && product.back().is_zero()) {
            product.pop_back();
        }
    }

    auto reduced_f = test_reduce(poly.coeffs, mod25);
    auto reduced_product = test_reduce(product, mod25);

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
    EXPECT_TRUE((match)) << "product of 3 lifted factors = f (mod 25)";
}
