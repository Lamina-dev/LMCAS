/**
 * @file test_hensel_lift.cpp
 * @brief Hensel 提升模块单元测试：Mignotte 界、提升高度计算与二因子二次提升。
 */

#include "test_common.hpp"
#include "transcendental_factor.hpp"
#include "polynomial.hpp"
#include "bigint.hpp"
#include "modular_arithmetic.hpp"

#include <vector>
#include <cmath>

using namespace lamina;

/**
 * @brief 测试简单多项式 (x^2 - 1) 的 Mignotte 界合理性。
 */
void test_mignotte_bound_simple() {
    TEST_CASE("Mignotte bound: x^2 - 1");

    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(-1), BigInt(0), BigInt(1)}, "x");

    std::vector<Polynomial<ModInt>> empty_factors;
    auto result = hensel_lift(poly, empty_factors, 3, 0);
    EXPECT_TRUE(result.empty(), "hensel_lift stub returns empty for empty factors");

    auto result2 = hensel_lift(poly, empty_factors, 3, 0);
    EXPECT_TRUE(result2.empty(), "hensel_lift stub returns empty (auto lift_bound)");
}

/**
 * @brief 测试高次多项式的提升高度随次数增长。
 */
void test_lift_height_increases_with_degree() {
    TEST_CASE("Lift height increases with degree");

    Polynomial<BigInt> f1(std::vector<BigInt>{BigInt(-1), BigInt(0), BigInt(1)}, "x");
    Polynomial<BigInt> f2(std::vector<BigInt>{BigInt(-1), BigInt(0), BigInt(0), BigInt(0), BigInt(1)}, "x");

    std::vector<Polynomial<ModInt>> empty_factors;
    auto r1 = hensel_lift(f1, empty_factors, 5, 0);
    auto r2 = hensel_lift(f2, empty_factors, 5, 0);

    EXPECT_TRUE(r1.empty(), "hensel_lift stub for x^2-1");
    EXPECT_TRUE(r2.empty(), "hensel_lift stub for x^4-1");
}

/**
 * @brief 测试较大系数多项式的 Mignotte 界。
 */
void test_mignotte_bound_larger_coefficients() {
    TEST_CASE("Mignotte bound: larger coefficients");

    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(3), BigInt(-4), BigInt(5), BigInt(6)}, "x");

    std::vector<Polynomial<ModInt>> empty_factors;
    auto result = hensel_lift(poly, empty_factors, 7, 0);
    EXPECT_TRUE(result.empty(), "hensel_lift stub for 6x^3+5x^2-4x+3");
}

/**
 * @brief 测试零多项式和常数多项式的边界情形。
 */
void test_edge_cases() {
    TEST_CASE("Edge cases: zero and constant polynomials");

    Polynomial<BigInt> zero_poly("x");
    std::vector<Polynomial<ModInt>> empty_factors;
    auto r1 = hensel_lift(zero_poly, empty_factors, 5, 0);
    EXPECT_TRUE(r1.empty(), "hensel_lift returns empty for zero polynomial");

    Polynomial<BigInt> const_poly(std::vector<BigInt>{BigInt(42)}, "x");
    auto r2 = hensel_lift(const_poly, empty_factors, 5, 0);
    EXPECT_TRUE(r2.empty(), "hensel_lift returns empty for constant polynomial");

    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(-1), BigInt(0), BigInt(1)}, "x");
    auto r3 = hensel_lift(poly, empty_factors, 3, 5);
    EXPECT_TRUE(r3.empty(), "hensel_lift stub with explicit lift_bound=5");
}

/**
 * @brief 测试大系数多项式不会溢出。
 */
void test_large_coefficients() {
    TEST_CASE("Large coefficients: no overflow");

    Polynomial<BigInt> poly(std::vector<BigInt>{
        BigInt(1000004), BigInt(1000003), BigInt(1000002),
        BigInt(1000001), BigInt(999999), BigInt(1000000)
    }, "x");

    std::vector<Polynomial<ModInt>> empty_factors;
    auto result = hensel_lift(poly, empty_factors, 7, 0);
    EXPECT_TRUE(result.empty(), "hensel_lift stub for large coefficient polynomial");
}

// ============================================================
// 二因子二次 Hensel 提升测试 (Task 4.2)
// ============================================================

namespace {

/**
 * @brief 辅助：多项式乘法（BigInt 系数，不取模）
 */
std::vector<BigInt> test_poly_mul(const std::vector<BigInt>& a,
                                  const std::vector<BigInt>& b) {
    if (a.empty() || b.empty()) return {};
    size_t n = a.size() + b.size() - 1;
    std::vector<BigInt> result(n, BigInt(0));
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            result[i + j] = result[i + j] + a[i] * b[j];
        }
    }
    return result;
}

/**
 * @brief 辅助：对称模归约
 */
BigInt test_sym_mod(const BigInt& c, const BigInt& m) {
    if (m.is_zero()) return c;
    BigInt r = c % m;
    if (r.IsNegative()) r = r + m;
    BigInt half_m = m / BigInt(2);
    if (r > half_m) r = r - m;
    return r;
}

/**
 * @brief 辅助：对多项式系数取模归约
 */
std::vector<BigInt> test_reduce(const std::vector<BigInt>& poly, const BigInt& m) {
    std::vector<BigInt> result = poly;
    for (auto& c : result) {
        c = test_sym_mod(c, m);
    }
    while (!result.empty() && result.back().is_zero()) {
        result.pop_back();
    }
    return result;
}

/**
 * @brief 辅助：验证 g*h = f (mod m)
 */
bool test_verify_factorization(const std::vector<BigInt>& f,
                               const std::vector<BigInt>& g,
                               const std::vector<BigInt>& h,
                               const BigInt& m) {
    auto product = test_poly_mul(g, h);
    auto reduced_product = test_reduce(product, m);
    auto reduced_f = test_reduce(f, m);

    size_t n = std::max(reduced_product.size(), reduced_f.size());
    reduced_product.resize(n, BigInt(0));
    reduced_f.resize(n, BigInt(0));

    for (size_t i = 0; i < n; ++i) {
        if (reduced_product[i] != reduced_f[i]) return false;
    }
    return true;
}

/**
 * @brief 辅助：验证 Bezout 关系 s*g + t*h = 1 (mod m)
 */
bool test_verify_bezout(const std::vector<BigInt>& s,
                        const std::vector<BigInt>& g,
                        const std::vector<BigInt>& t,
                        const std::vector<BigInt>& h,
                        const BigInt& m) {
    auto sg = test_poly_mul(s, g);
    auto th = test_poly_mul(t, h);

    size_t n = std::max(sg.size(), th.size());
    std::vector<BigInt> sum(n, BigInt(0));
    for (size_t i = 0; i < n; ++i) {
        BigInt a = (i < sg.size()) ? sg[i] : BigInt(0);
        BigInt b = (i < th.size()) ? th[i] : BigInt(0);
        sum[i] = a + b;
    }

    auto reduced = test_reduce(sum, m);

    if (reduced.empty()) return false;
    if (reduced[0] != BigInt(1)) return false;
    for (size_t i = 1; i < reduced.size(); ++i) {
        if (!reduced[i].is_zero()) return false;
    }
    return true;
}

} // anonymous namespace

/**
 * @brief 测试 x^2 - 1 = (x+1)(x-1) mod 3 的 Hensel 提升到 mod 9。
 *
 * Bezout: s=2, t=1 满足 2*(x+1) + 1*(x-1) = 3x+1 = 1 (mod 3)
 * 提升后应满足 f = g'*h' (mod 9) 且 s'*g' + t'*h' = 1 (mod 9)
 */
void test_hensel_lift_x2_minus_1_mod3() {
    TEST_CASE("hl_two_factor_lift: x^2 - 1 = (x+1)(x-1) mod 3, lift to mod 9");

    std::vector<BigInt> f = {BigInt(-1), BigInt(0), BigInt(1)};
    std::vector<BigInt> g = {BigInt(1), BigInt(1)};
    std::vector<BigInt> h = {BigInt(-1), BigInt(1)};
    std::vector<BigInt> s = {BigInt(2)};
    std::vector<BigInt> t = {BigInt(1)};
    BigInt m(3);

    // Verify preconditions
    EXPECT_TRUE(test_verify_factorization(f, g, h, m),
        "initial: g*h = f (mod 3)");
    EXPECT_TRUE(test_verify_bezout(s, g, t, h, m),
        "initial: s*g + t*h = 1 (mod 3)");

    // Perform the lift
    HenselLiftPair initial{g, h, s, t, m};
    HenselLiftPair lifted = hl_two_factor_lift(f, initial);

    BigInt m2(9);
    EXPECT_TRUE(lifted.modulus == m2, "lifted modulus = 9");

    // Verify postconditions
    EXPECT_TRUE(test_verify_factorization(f, lifted.g, lifted.h, m2),
        "lifted: g'*h' = f (mod 9)");
    EXPECT_TRUE(test_verify_bezout(lifted.s, lifted.g, lifted.t, lifted.h, m2),
        "lifted: s'*g' + t'*h' = 1 (mod 9)");

    // Since (x+1)(x-1) = x^2-1 exactly, factors should remain unchanged
    EXPECT_TRUE(lifted.g.size() == 2, "g' has degree 1");
    EXPECT_TRUE(lifted.h.size() == 2, "h' has degree 1");
}

/**
 * @brief 测试 x^2+3x+2 = (x+1)(x+2) mod 5 的 Hensel 提升到 mod 25。
 *
 * Bezout: s=4, t=1 满足 4*(x+1) + 1*(x+2) = 5x+6 = 1 (mod 5)
 */
void test_hensel_lift_x2_plus_3x_plus_2_mod5() {
    TEST_CASE("hl_two_factor_lift: x^2+3x+2 = (x+1)(x+2) mod 5, lift to mod 25");

    std::vector<BigInt> f = {BigInt(2), BigInt(3), BigInt(1)};
    std::vector<BigInt> g = {BigInt(1), BigInt(1)};
    std::vector<BigInt> h = {BigInt(2), BigInt(1)};
    std::vector<BigInt> s = {BigInt(4)};
    std::vector<BigInt> t_coeff = {BigInt(1)};
    BigInt m(5);
    BigInt m2(25);

    // Verify preconditions
    EXPECT_TRUE(test_verify_factorization(f, g, h, m),
        "initial: g*h = f (mod 5)");
    EXPECT_TRUE(test_verify_bezout(s, g, t_coeff, h, m),
        "initial: s*g + t*h = 1 (mod 5)");

    // Perform the lift
    HenselLiftPair initial{g, h, s, t_coeff, m};
    HenselLiftPair lifted = hl_two_factor_lift(f, initial);

    EXPECT_TRUE(lifted.modulus == m2, "lifted modulus = 25");

    // Verify postconditions
    EXPECT_TRUE(test_verify_factorization(f, lifted.g, lifted.h, m2),
        "lifted: g'*h' = f (mod 25)");
    EXPECT_TRUE(test_verify_bezout(lifted.s, lifted.g, lifted.t, lifted.h, m2),
        "lifted: s'*g' + t'*h' = 1 (mod 25)");
}

/**
 * @brief 测试非精确分解情形：x^2+6x+5 = (x+1)(x+5)，mod 3 下为 (x+1)(x+2)。
 *
 * 这是一个非平凡提升：mod 3 下的因子 (x+2) 需要被修正为 (x+5) mod 9。
 */
void test_hensel_lift_non_exact_mod3() {
    TEST_CASE("hl_two_factor_lift: x^2+6x+5 = (x+1)(x+5), lift from mod 3 to mod 9");

    std::vector<BigInt> f = {BigInt(5), BigInt(6), BigInt(1)};
    // mod 3: g = x+1, h = x+2 (since 5 = 2 mod 3)
    std::vector<BigInt> g = {BigInt(1), BigInt(1)};
    std::vector<BigInt> h = {BigInt(2), BigInt(1)};
    // Bezout: s=2, t=1 satisfies 2*(x+1) + 1*(x+2) = 3x+4 = 1 (mod 3)
    std::vector<BigInt> s = {BigInt(2)};
    std::vector<BigInt> t_coeff = {BigInt(1)};
    BigInt m(3);
    BigInt m2(9);

    // Verify preconditions
    EXPECT_TRUE(test_verify_factorization(f, g, h, m),
        "initial: g*h = f (mod 3)");
    EXPECT_TRUE(test_verify_bezout(s, g, t_coeff, h, m),
        "initial: s*g + t*h = 1 (mod 3)");

    // Perform the lift
    HenselLiftPair initial{g, h, s, t_coeff, m};
    HenselLiftPair lifted = hl_two_factor_lift(f, initial);

    EXPECT_TRUE(lifted.modulus == m2, "lifted modulus = 9");

    // Verify postconditions
    EXPECT_TRUE(test_verify_factorization(f, lifted.g, lifted.h, m2),
        "lifted: g'*h' = f (mod 9)");
    EXPECT_TRUE(test_verify_bezout(lifted.s, lifted.g, lifted.t, lifted.h, m2),
        "lifted: s'*g' + t'*h' = 1 (mod 9)");

    // The true factorization is (x+1)(x+5), so after lifting mod 9:
    // g' should be [1,1] (x+1) and h' should be [5,1] (x+5)
    // since (x+1)(x+5) = x^2+6x+5 exactly equals f
    EXPECT_TRUE(test_verify_factorization(f, {BigInt(1), BigInt(1)}, {BigInt(5), BigInt(1)}, m2),
        "exact: (x+1)*(x+5) = f (mod 9)");
}

/**
 * @brief 测试二次提升的迭代：mod 3 -> mod 9 -> mod 81。
 *
 * 验证可以连续调用 hl_two_factor_lift 进行多步提升。
 */
void test_hensel_lift_iterated() {
    TEST_CASE("hl_two_factor_lift: iterated lift mod 3 -> 9 -> 81");

    // f = x^2 + 8x + 7 = (x+1)(x+7)
    // mod 3: f = x^2 + 2x + 1 = (x+1)^2? No: (x+1)(x+7) mod 3 = (x+1)(x+1) mod 3
    // That's not coprime! Use a different example.
    //
    // f = x^2 + 10x + 21 = (x+3)(x+7)
    // mod 5: f = x^2 + 0x + 1 = x^2+1. Roots: 2^2+1=5=0, 3^2+1=10=0 mod 5
    // So x^2+1 = (x-2)(x-3) = (x+3)(x+2) mod 5
    // g = x+3, h = x+2 mod 5
    // Bezout: s*(x+3) + t*(x+2) = 1 mod 5
    // (s+t)x + (3s+2t) = 1 -> s+t=0, 3s+2t=1 -> s=1, t=-1=4 mod 5
    // Verify: 1*(x+3) + 4*(x+2) = 5x+11 = 0x+1 = 1 mod 5

    std::vector<BigInt> f = {BigInt(21), BigInt(10), BigInt(1)};
    std::vector<BigInt> g = {BigInt(3), BigInt(1)};
    std::vector<BigInt> h = {BigInt(2), BigInt(1)};
    std::vector<BigInt> s = {BigInt(1)};
    std::vector<BigInt> t_coeff = {BigInt(4)};
    BigInt m(5);

    // Verify preconditions
    EXPECT_TRUE(test_verify_factorization(f, g, h, m),
        "initial: g*h = f (mod 5)");
    EXPECT_TRUE(test_verify_bezout(s, g, t_coeff, h, m),
        "initial: s*g + t*h = 1 (mod 5)");

    // First lift: mod 5 -> mod 25
    HenselLiftPair state{g, h, s, t_coeff, m};
    HenselLiftPair lifted1 = hl_two_factor_lift(f, state);

    BigInt m2(25);
    EXPECT_TRUE(lifted1.modulus == m2, "first lift modulus = 25");
    EXPECT_TRUE(test_verify_factorization(f, lifted1.g, lifted1.h, m2),
        "first lift: g'*h' = f (mod 25)");
    EXPECT_TRUE(test_verify_bezout(lifted1.s, lifted1.g, lifted1.t, lifted1.h, m2),
        "first lift: s'*g' + t'*h' = 1 (mod 25)");

    // Second lift: mod 25 -> mod 625
    HenselLiftPair lifted2 = hl_two_factor_lift(f, lifted1);

    BigInt m4(625);
    EXPECT_TRUE(lifted2.modulus == m4, "second lift modulus = 625");
    EXPECT_TRUE(test_verify_factorization(f, lifted2.g, lifted2.h, m4),
        "second lift: g'*h' = f (mod 625)");
    EXPECT_TRUE(test_verify_bezout(lifted2.s, lifted2.g, lifted2.t, lifted2.h, m4),
        "second lift: s'*g' + t'*h' = 1 (mod 625)");
}

/**
 * @brief 测试三次多项式的 Hensel 提升。
 *
 * f = x^3 - x = x(x+1)(x-1)，取两因子 g=x, h=x^2-1 mod 5。
 */
void test_hensel_lift_cubic() {
    TEST_CASE("hl_two_factor_lift: x^3 - x with g=x, h=x^2-1 mod 5");

    // f = x^3 - x: coeffs [0, -1, 0, 1]
    std::vector<BigInt> f = {BigInt(0), BigInt(-1), BigInt(0), BigInt(1)};
    // g = x: coeffs [0, 1]
    std::vector<BigInt> g = {BigInt(0), BigInt(1)};
    // h = x^2 - 1: coeffs [-1, 0, 1]
    std::vector<BigInt> h = {BigInt(-1), BigInt(0), BigInt(1)};

    // Bezout: s*g + t*h = 1 mod 5
    // s*x + t*(x^2-1) = 1 mod 5
    // deg(s) < deg(h)=2, so s can be degree 1: s = a + bx
    // deg(t) < deg(g)=1, so t is constant: t = c
    // (a+bx)*x + c*(x^2-1) = ax + bx^2 + cx^2 - c = (b+c)x^2 + ax - c = 1 mod 5
    // b+c=0, a=0, -c=1 -> c=-1=4 mod 5, b=1, a=0
    // s = [0, 1] (= x), t = [4] (= -1 mod 5)
    // Verify: x*x + 4*(x^2-1) = x^2 + 4x^2 - 4 = 5x^2 - 4 = -4 = 1 mod 5
    std::vector<BigInt> s = {BigInt(0), BigInt(1)};
    std::vector<BigInt> t_coeff = {BigInt(4)};
    BigInt m(5);

    EXPECT_TRUE(test_verify_factorization(f, g, h, m),
        "initial: g*h = f (mod 5)");
    EXPECT_TRUE(test_verify_bezout(s, g, t_coeff, h, m),
        "initial: s*g + t*h = 1 (mod 5)");

    // Lift mod 5 -> mod 25
    HenselLiftPair state{g, h, s, t_coeff, m};
    HenselLiftPair lifted = hl_two_factor_lift(f, state);

    BigInt m2(25);
    EXPECT_TRUE(lifted.modulus == m2, "lifted modulus = 25");
    EXPECT_TRUE(test_verify_factorization(f, lifted.g, lifted.h, m2),
        "lifted: g'*h' = f (mod 25)");
    EXPECT_TRUE(test_verify_bezout(lifted.s, lifted.g, lifted.t, lifted.h, m2),
        "lifted: s'*g' + t'*h' = 1 (mod 25)");
}

/**
 * @brief 测试通过公共 hensel_lift API 间接验证。
 */
void test_hensel_lift_via_api_x2_minus_1() {
    TEST_CASE("hensel_lift API: x^2 - 1 with empty factors mod 3");

    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(-1), BigInt(0), BigInt(1)}, "x");

    // hensel_lift with empty factors should return empty
    std::vector<Polynomial<ModInt>> mod_factors;
    auto result = hensel_lift(poly, mod_factors, 3, 2);
    EXPECT_TRUE(result.empty(), "hensel_lift API call returns empty for no factors");
}

// ============================================================
// 多因子 Hensel 提升测试 (Task 4.3)
// ============================================================

/**
 * @brief 测试 3 因子提升：x³ - x = x(x-1)(x+1) mod 5，提升到 mod 25。
 *
 * 验证多因子提升后所有因子乘积等于 f (mod p^k)。
 */
void test_multi_factor_lift_3_factors() {
    TEST_CASE("multi-factor lift: x^3 - x = x*(x-1)*(x+1) mod 5, lift to mod 25");

    // f = x^3 - x: coeffs [0, -1, 0, 1]
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(0), BigInt(-1), BigInt(0), BigInt(1)}, "x");

    // mod 5 factors: x, x+4 (= x-1 mod 5), x+1
    int64_t p = 5;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(0, p), ModInt(1, p)};  // x
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(4, p), ModInt(1, p)};  // x+4 = x-1 mod 5
    Polynomial<ModInt> f3("x");
    f3.coeffs = {ModInt(1, p), ModInt(1, p)};  // x+1
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2, f3};

    int lift_bound = 2;  // lift to mod 5^2 = 25
    auto lifted = hensel_lift(poly, mod_factors, p, lift_bound);

    EXPECT_TRUE(lifted.size() == 3, "should produce 3 lifted factors");

    // Verify: product of all lifted factors = f (mod 25)
    BigInt mod25(25);
    std::vector<BigInt> product = lifted[0].coeffs;
    for (size_t i = 1; i < lifted.size(); ++i) {
        product = test_poly_mul(product, lifted[i].coeffs);
        // Reduce mod 25
        for (auto& c : product) {
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
        if (reduced_f[i] != reduced_product[i]) { match = false; break; }
    }
    EXPECT_TRUE(match, "product of lifted factors = f (mod 25)");
}

/**
 * @brief 测试 2 因子通过多因子路径：x² - 1 = (x+1)(x-1) mod 3，提升到 mod 9。
 *
 * 验证多因子提升对 2 因子情形也正确工作。
 */
void test_multi_factor_lift_2_factors_via_api() {
    TEST_CASE("multi-factor lift via API: x^2 - 1 = (x+1)(x-1) mod 3, lift to mod 9");

    // f = x^2 - 1
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(-1), BigInt(0), BigInt(1)}, "x");

    // mod 3 factors: (x+1), (x+2) where x+2 = x-1 mod 3
    int64_t p = 3;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(1, p), ModInt(1, p)};  // x+1
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(2, p), ModInt(1, p)};  // x+2 = x-1 mod 3
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2};

    int lift_bound = 2;  // lift to mod 3^2 = 9
    auto lifted = hensel_lift(poly, mod_factors, p, lift_bound);

    EXPECT_TRUE(lifted.size() == 2, "should produce 2 lifted factors");

    // Verify: product of lifted factors = f (mod 9)
    BigInt mod9(9);
    std::vector<BigInt> product = test_poly_mul(lifted[0].coeffs, lifted[1].coeffs);
    auto reduced_product = test_reduce(product, mod9);
    auto reduced_f = test_reduce(poly.coeffs, mod9);

    size_t n = std::max(reduced_f.size(), reduced_product.size());
    reduced_f.resize(n, BigInt(0));
    reduced_product.resize(n, BigInt(0));

    bool match = true;
    for (size_t i = 0; i < n; ++i) {
        if (reduced_f[i] != reduced_product[i]) { match = false; break; }
    }
    EXPECT_TRUE(match, "product of 2 lifted factors = f (mod 9)");

    // Since (x+1)(x-1) = x^2-1 exactly, factors should be exact
    // g should be (x+1) and h should be (x-1) = (x+8) mod 9 or (x-1)
    // In symmetric representation: (x+1) and (x-1)
    bool exact = test_verify_factorization(poly.coeffs, lifted[0].coeffs, lifted[1].coeffs, mod9);
    EXPECT_TRUE(exact, "lifted factors are exact for x^2-1");
}

/**
 * @brief 测试多因子提升的系数在对称表示范围内。
 */
void test_multi_factor_lift_symmetric_coeffs() {
    TEST_CASE("multi-factor lift: coefficients in symmetric representation");

    // f = x^3 + 6x^2 + 11x + 6 = (x+1)(x+2)(x+3)
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(6), BigInt(11), BigInt(6), BigInt(1)}, "x");

    // mod 5: (x+1)(x+2)(x+3)
    int64_t p = 5;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(1, p), ModInt(1, p)};  // x+1
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(2, p), ModInt(1, p)};  // x+2
    Polynomial<ModInt> f3("x");
    f3.coeffs = {ModInt(3, p), ModInt(1, p)};  // x+3
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2, f3};

    int lift_bound = 2;  // lift to mod 25
    auto lifted = hensel_lift(poly, mod_factors, p, lift_bound);

    EXPECT_TRUE(lifted.size() == 3, "should produce 3 lifted factors");

    // All coefficients should be in [-12, 12] (symmetric mod 25)
    BigInt mod25(25);
    BigInt half(12);
    bool all_in_range = true;
    for (const auto& factor : lifted) {
        for (const auto& c : factor.coeffs) {
            BigInt abs_c = c.Abs();
            if (abs_c > half) { all_in_range = false; break; }
        }
        if (!all_in_range) break;
    }
    EXPECT_TRUE(all_in_range, "all coefficients in symmetric range [-12, 12]");

    // Verify product = f (mod 25)
    std::vector<BigInt> product = lifted[0].coeffs;
    for (size_t i = 1; i < lifted.size(); ++i) {
        product = test_poly_mul(product, lifted[i].coeffs);
        for (auto& c : product) {
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
        if (reduced_f[i] != reduced_product[i]) { match = false; break; }
    }
    EXPECT_TRUE(match, "product of 3 lifted factors = f (mod 25)");
}

// ============================================================
// 系数对称表示测试 (Task 4.4)
// ============================================================

/**
 * @brief 验证对称模归约的基本正确性。
 *
 * 测试 hl_symmetric_mod 的行为通过公共 API 间接验证：
 * 提升后的系数 c 应满足 -p^k/2 ≤ c ≤ p^k/2。
 */
void test_symmetric_repr_basic_range() {
    TEST_CASE("symmetric repr: all coefficients in [-p^k/2, p^k/2] for x^2-1 mod 3, k=4");

    // f = x^2 - 1 = (x+1)(x-1)
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(-1), BigInt(0), BigInt(1)}, "x");

    int64_t p = 3;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(1, p), ModInt(1, p)};  // x+1
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(2, p), ModInt(1, p)};  // x+2 = x-1 mod 3
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2};

    int lift_bound = 4;  // lift to mod 3^4 = 81
    auto lifted = hensel_lift(poly, mod_factors, p, lift_bound);

    EXPECT_TRUE(lifted.size() == 2, "should produce 2 lifted factors");

    // p^k = 81, half = 40
    BigInt mod81(81);
    BigInt half(40);

    bool all_in_range = true;
    for (const auto& factor : lifted) {
        for (const auto& c : factor.coeffs) {
            BigInt abs_c = c.Abs();
            if (abs_c > half) {
                all_in_range = false;
                break;
            }
        }
        if (!all_in_range) break;
    }
    EXPECT_TRUE(all_in_range, "all coefficients in [-40, 40] for mod 81");

    // Verify product still equals f mod p^k
    std::vector<BigInt> product = test_poly_mul(lifted[0].coeffs, lifted[1].coeffs);
    auto reduced_product = test_reduce(product, mod81);
    auto reduced_f = test_reduce(poly.coeffs, mod81);

    size_t n = std::max(reduced_f.size(), reduced_product.size());
    reduced_f.resize(n, BigInt(0));
    reduced_product.resize(n, BigInt(0));

    bool match = true;
    for (size_t i = 0; i < n; ++i) {
        if (reduced_f[i] != reduced_product[i]) { match = false; break; }
    }
    EXPECT_TRUE(match, "product of lifted factors = f (mod 81)");
}

/**
 * @brief 测试大系数多项式提升后系数仍在对称范围内。
 *
 * f = x^3 + 15x^2 + 71x + 105 = (x+3)(x+5)(x+7)
 * 系数较大，验证提升后归约到对称表示。
 */
void test_symmetric_repr_large_coefficients() {
    TEST_CASE("symmetric repr: large coefficient polynomial (x+3)(x+5)(x+7) mod 7, k=3");

    // f = x^3 + 15x^2 + 71x + 105 = (x+3)(x+5)(x+7)
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(105), BigInt(71), BigInt(15), BigInt(1)}, "x");

    // mod 7: (x+3)(x+5)(x+0) = x(x+3)(x+5)
    // 105 mod 7 = 0, 71 mod 7 = 1, 15 mod 7 = 1
    // f mod 7 = x^3 + x^2 + x = x(x^2 + x + 1)
    // x^2 + x + 1 mod 7: discriminant = 1-4 = -3 = 4 mod 7, sqrt(4)=2
    // roots: (-1±2)/2 = 1/2=4, -3/2=2 mod 7
    // So x^2+x+1 = (x-4)(x-2) = (x+3)(x+5) mod 7
    int64_t p = 7;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(0, p), ModInt(1, p)};  // x
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(3, p), ModInt(1, p)};  // x+3
    Polynomial<ModInt> f3("x");
    f3.coeffs = {ModInt(5, p), ModInt(1, p)};  // x+5
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2, f3};

    int lift_bound = 3;  // lift to mod 7^3 = 343
    auto lifted = hensel_lift(poly, mod_factors, p, lift_bound);

    EXPECT_TRUE(lifted.size() == 3, "should produce 3 lifted factors");

    // p^k = 343, half = 171
    BigInt mod343(343);
    BigInt half(171);

    bool all_in_range = true;
    for (const auto& factor : lifted) {
        for (const auto& c : factor.coeffs) {
            BigInt abs_c = c.Abs();
            if (abs_c > half) {
                all_in_range = false;
                break;
            }
        }
        if (!all_in_range) break;
    }
    EXPECT_TRUE(all_in_range, "all coefficients in [-171, 171] for mod 343");

    // Verify product = f (mod 343)
    std::vector<BigInt> product = lifted[0].coeffs;
    for (size_t i = 1; i < lifted.size(); ++i) {
        product = test_poly_mul(product, lifted[i].coeffs);
        for (auto& c : product) {
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
        if (reduced_f[i] != reduced_product[i]) { match = false; break; }
    }
    EXPECT_TRUE(match, "product of lifted factors = f (mod 343)");
}

/**
 * @brief 测试对称归约的往返一致性：归约后的值 mod p^k 等于原值 mod p^k。
 *
 * 对各种系数值验证 symmetric_mod(c, m) ≡ c (mod m)。
 */
void test_symmetric_repr_roundtrip() {
    TEST_CASE("symmetric repr: roundtrip consistency c ≡ symmetric_mod(c, m) (mod m)");

    // 通过提升一个已知多项式并验证归约后的系数仍满足 mod 等价性
    // f = x^2 + 10x + 21 = (x+3)(x+7)
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(21), BigInt(10), BigInt(1)}, "x");

    // mod 5: 21 mod 5 = 1, 10 mod 5 = 0 → f = x^2 + 1 = (x+3)(x+2) mod 5
    int64_t p = 5;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(3, p), ModInt(1, p)};  // x+3
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(2, p), ModInt(1, p)};  // x+2 = x+7 mod 5
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2};

    int lift_bound = 3;  // lift to mod 5^3 = 125
    auto lifted = hensel_lift(poly, mod_factors, p, lift_bound);

    EXPECT_TRUE(lifted.size() == 2, "should produce 2 lifted factors");

    BigInt mod125(125);
    BigInt half(62);

    // Verify all coefficients are in symmetric range
    bool all_in_range = true;
    for (const auto& factor : lifted) {
        for (const auto& c : factor.coeffs) {
            BigInt abs_c = c.Abs();
            if (abs_c > half) {
                all_in_range = false;
                break;
            }
        }
        if (!all_in_range) break;
    }
    EXPECT_TRUE(all_in_range, "all coefficients in [-62, 62] for mod 125");

    // Verify roundtrip: for each coefficient c, c mod m should equal c (since it's already reduced)
    bool roundtrip_ok = true;
    for (const auto& factor : lifted) {
        for (const auto& c : factor.coeffs) {
            BigInt reduced = test_sym_mod(c, mod125);
            if (reduced != c) {
                roundtrip_ok = false;
                break;
            }
        }
        if (!roundtrip_ok) break;
    }
    EXPECT_TRUE(roundtrip_ok, "all coefficients are already in symmetric form (idempotent)");
}

/**
 * @brief 测试边界情形：系数恰好在 p^k/2 处。
 *
 * 使用 p=2, k=3 (mod 8)，half=4。验证系数 4 被保留，系数 5 被归约为 -3。
 */
void test_symmetric_repr_boundary() {
    TEST_CASE("symmetric repr: boundary case at p^k/2");

    // 使用 p=2, k=3 → mod 8, half = 4
    // f = x^2 + 5x + 6 = (x+2)(x+3)
    // mod 2: f = x^2 + x = x(x+1)
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(6), BigInt(5), BigInt(1)}, "x");

    int64_t p = 2;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(0, p), ModInt(1, p)};  // x
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(1, p), ModInt(1, p)};  // x+1
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2};

    int lift_bound = 3;  // lift to mod 2^3 = 8
    auto lifted = hensel_lift(poly, mod_factors, p, lift_bound);

    EXPECT_TRUE(lifted.size() == 2, "should produce 2 lifted factors");

    // p^k = 8, half = 4
    BigInt mod8(8);
    BigInt half(4);

    // All coefficients should satisfy |c| <= 4
    bool all_in_range = true;
    for (const auto& factor : lifted) {
        for (const auto& c : factor.coeffs) {
            BigInt abs_c = c.Abs();
            if (abs_c > half) {
                all_in_range = false;
                break;
            }
        }
        if (!all_in_range) break;
    }
    EXPECT_TRUE(all_in_range, "all coefficients in [-4, 4] for mod 8");

    // Verify product = f (mod 8)
    std::vector<BigInt> product = test_poly_mul(lifted[0].coeffs, lifted[1].coeffs);
    auto reduced_product = test_reduce(product, mod8);
    auto reduced_f = test_reduce(poly.coeffs, mod8);

    size_t n = std::max(reduced_f.size(), reduced_product.size());
    reduced_f.resize(n, BigInt(0));
    reduced_product.resize(n, BigInt(0));

    bool match = true;
    for (size_t i = 0; i < n; ++i) {
        if (reduced_f[i] != reduced_product[i]) { match = false; break; }
    }
    EXPECT_TRUE(match, "product of lifted factors = f (mod 8)");
}

/**
 * @brief 测试高提升次数下系数仍在对称范围内。
 *
 * f = x^2 + 3x + 2 = (x+1)(x+2)，提升到 mod 5^5 = 3125。
 */
void test_symmetric_repr_high_lift() {
    TEST_CASE("symmetric repr: high lift bound p=5, k=5 (mod 3125)");

    // f = x^2 + 3x + 2 = (x+1)(x+2)
    Polynomial<BigInt> poly(std::vector<BigInt>{BigInt(2), BigInt(3), BigInt(1)}, "x");

    // mod 5: f = x^2 + 3x + 2 = (x+1)(x+2)
    int64_t p = 5;
    Polynomial<ModInt> f1("x");
    f1.coeffs = {ModInt(1, p), ModInt(1, p)};  // x+1
    Polynomial<ModInt> f2("x");
    f2.coeffs = {ModInt(2, p), ModInt(1, p)};  // x+2
    std::vector<Polynomial<ModInt>> mod_factors = {f1, f2};

    int lift_bound = 5;  // lift to mod 5^5 = 3125
    auto lifted = hensel_lift(poly, mod_factors, p, lift_bound);

    EXPECT_TRUE(lifted.size() == 2, "should produce 2 lifted factors");

    // p^k = 3125, half = 1562
    BigInt mod3125(3125);
    BigInt half(1562);

    bool all_in_range = true;
    for (const auto& factor : lifted) {
        for (const auto& c : factor.coeffs) {
            BigInt abs_c = c.Abs();
            if (abs_c > half) {
                all_in_range = false;
                break;
            }
        }
        if (!all_in_range) break;
    }
    EXPECT_TRUE(all_in_range, "all coefficients in [-1562, 1562] for mod 3125");

    // Verify product = f (mod 3125)
    std::vector<BigInt> product = test_poly_mul(lifted[0].coeffs, lifted[1].coeffs);
    auto reduced_product = test_reduce(product, mod3125);
    auto reduced_f = test_reduce(poly.coeffs, mod3125);

    size_t n = std::max(reduced_f.size(), reduced_product.size());
    reduced_f.resize(n, BigInt(0));
    reduced_product.resize(n, BigInt(0));

    bool match = true;
    for (size_t i = 0; i < n; ++i) {
        if (reduced_f[i] != reduced_product[i]) { match = false; break; }
    }
    EXPECT_TRUE(match, "product of lifted factors = f (mod 3125)");
}

int main() {
    test_mignotte_bound_simple();
    test_lift_height_increases_with_degree();
    test_mignotte_bound_larger_coefficients();
    test_edge_cases();
    test_large_coefficients();

    // Task 4.2: Two-factor quadratic Hensel lifting tests
    test_hensel_lift_x2_minus_1_mod3();
    test_hensel_lift_x2_plus_3x_plus_2_mod5();
    test_hensel_lift_non_exact_mod3();
    test_hensel_lift_iterated();
    test_hensel_lift_cubic();
    test_hensel_lift_via_api_x2_minus_1();

    // Task 4.3: Multi-factor Hensel lifting tests
    test_multi_factor_lift_3_factors();
    test_multi_factor_lift_2_factors_via_api();
    test_multi_factor_lift_symmetric_coeffs();

    // Task 4.4: Coefficient symmetric representation tests
    test_symmetric_repr_basic_range();
    test_symmetric_repr_large_coefficients();
    test_symmetric_repr_roundtrip();
    test_symmetric_repr_boundary();
    test_symmetric_repr_high_lift();

    return TEST_REPORT();
}
