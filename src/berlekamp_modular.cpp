#include "internal/berlekamp_support.hpp"
#include <algorithm>
#include <limits>

namespace LMCAS::berlekamp_detail {

/**
 * @brief Berlekamp 模分解的小素数候选表。
 * 按升序选取首个不整除首项系数且使模约化结果无平方因子的素数。
 * @internal
 */
static constexpr int64_t BK_SMALL_PRIMES[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229
};
static constexpr size_t BK_NUM_SMALL_PRIMES =
    sizeof(BK_SMALL_PRIMES) / sizeof(BK_SMALL_PRIMES[0]);

/**
 * @brief 以系数分母的 LCM 清除分母，再逐项取模 p。
 * @param[in] poly 有理系数多项式
 * @param[in] p 素数模数
 * @return 升幂排列的模 p 系数向量，去除高次零系数。
 * @internal
 */
std::vector<int64_t> bk_reduce_to_mod_coeffs(
    const Polynomial<Rational>& poly,
    int64_t p) {

    if (poly.is_zero()) {
        return {};
    }

    BigInt lcm_den(1);
    for (const auto& c : poly.coeffs) {
        BigInt den = c.get_denominator();
        BigInt g = BigInt::gcd(lcm_den, den);
        lcm_den = lcm_den * (den / g);
    }

    BigInt big_p(static_cast<std::int64_t>(p));

    std::vector<int64_t> mod_coeffs;
    mod_coeffs.reserve(poly.coeffs.size());

    for (const auto& c : poly.coeffs) {
        BigInt num = c.get_numerator();
        BigInt den = c.get_denominator();
        BigInt integer_coeff = num * (lcm_den / den);

        BigInt rem = integer_coeff % big_p;
        int64_t rem_val = rem.try_to_int64().value();
        if (rem_val < 0) rem_val += p;

        mod_coeffs.push_back(rem_val);
    }

    while (!mod_coeffs.empty() && mod_coeffs.back() == 0) {
        mod_coeffs.pop_back();
    }

    return mod_coeffs;
}

/**
 * @brief 计算模 p 多项式的形式导数.
 *
 * @param[in] f 模 p 系数向量
 * @param[in] p 素数模数
 * @return f' 的系数向量
 * @internal
 */
static std::vector<int64_t> bk_derivative_mod(
    const std::vector<int64_t>& f,
    int64_t p) {

    if (f.size() <= 1) {
        return {};
    }

    std::vector<int64_t> deriv;
    deriv.reserve(f.size() - 1);

    for (size_t i = 1; i < f.size(); ++i) {
        const auto index_residue = static_cast<std::int64_t>(
            i % static_cast<std::uint64_t>(p));
        int64_t coeff = (f[i] * index_residue) % p;
        deriv.push_back(coeff);
    }

    while (!deriv.empty() && deriv.back() == 0) {
        deriv.pop_back();
    }

    return deriv;
}

/**
 * @brief 模 p 下求乘法逆元.
 *
 * 使用扩展欧几里得算法计算 a^(-1) mod p.
 *
 * @param[in] a 待求逆的值(非零)
 * @param[in] p 素数模数
 * @return a 在模 p 下的逆元
 * @internal
 */
int64_t bk_mod_inverse(int64_t a, int64_t p) {
    int64_t s, t;
    extended_gcd(((a % p) + p) % p, p, s, t);
    return ((s % p) + p) % p;
}

/**
 * @brief 模 p 多项式带余除法.
 *
 * 计算 a / b 的商和余数,使得 a = q * b + r,deg(r) < deg(b).
 *
 * @param[in]  a 被除数系数向量
 * @param[in]  b 除数系数向量(非空)
 * @param[in]  p 素数模数
 * @param[out] q 商系数向量
 * @param[out] r 余数系数向量
 * @internal
 */
void bk_div_mod(
    const std::vector<int64_t>& a,
    const std::vector<int64_t>& b,
    int64_t p,
    std::vector<int64_t>& q,
    std::vector<int64_t>& r) {

    q.clear();
    r = a;

    if (b.empty()) return;

    int deg_b = static_cast<int>(b.size()) - 1;
    int64_t lc_b_inv = bk_mod_inverse(b.back(), p);

    int deg_r = static_cast<int>(r.size()) - 1;

    if (deg_r < deg_b) {
        return;
    }

    q.resize(deg_r - deg_b + 1, 0);

    while (deg_r >= deg_b) {
        int64_t factor = (r[deg_r] * lc_b_inv) % p;
        int shift = deg_r - deg_b;
        q[shift] = factor;

        for (int i = 0; i <= deg_b; ++i) {
            r[shift + i] = (r[shift + i] - factor * b[i] % p + p) % p;
        }

        while (deg_r >= 0 && r[deg_r] == 0) {
            r.pop_back();
            deg_r--;
        }
    }
}

/**
 * @brief 计算两个模 p 多项式的最大公因式.
 *
 * 使用欧几里得算法.结果为首一多项式.
 *
 * @param[in] a 第一个多项式系数向量
 * @param[in] b 第二个多项式系数向量
 * @param[in] p 素数模数
 * @return gcd(a, b) 的首一化系数向量
 * @internal
 */
std::vector<int64_t> bk_gcd_mod(
    std::vector<int64_t> a,
    std::vector<int64_t> b,
    int64_t p) {

    while (!b.empty()) {
        std::vector<int64_t> q, r;
        bk_div_mod(a, b, p, q, r);
        a = std::move(b);
        b = std::move(r);
    }

    if (!a.empty()) {
        int64_t lc_inv = bk_mod_inverse(a.back(), p);
        for (auto& c : a) {
            c = (c * lc_inv) % p;
        }
    }

    return a;
}

/**
 * @brief 检查模 p 多项式是否为 无平方因子.
 *
 * 多项式 f 在 F_p 上 无平方因子 当且仅当 gcd(f, f') 的次数为 0.
 *
 * @param[in] f 模 p 系数向量
 * @param[in] p 素数模数
 * @return f 为 无平方因子 返回 true
 * @internal
 */
static bool bk_is_square_free_mod(
    const std::vector<int64_t>& f,
    int64_t p) {

    if (f.size() <= 1) return true;

    std::vector<int64_t> f_prime = bk_derivative_mod(f, p);
    if (f_prime.empty()) return false;

    std::vector<int64_t> g = bk_gcd_mod(f, f_prime, p);
    return g.size() <= 1;
}

/**
 * @brief 判断素数是否保持多项式次数和无平方性质。
 *
 * Berlekamp 分解要求模约化仍位于同一次数的多项式环中，并且
 * gcd(f mod p, f' mod p) 为常数。显式和自动素数选择共用此检查。
 *
 * @param[in] poly 非零有理系数多项式，次数 >= 1。
 * @param[in] prime 候选模数。
 * @return 候选为内核支持的素数且满足分解前置条件时返回 true。
 * @internal
 */
bool bk_is_suitable_prime(
    const Polynomial<Rational>& poly,
    int64_t prime) {
    if (poly.is_zero() || poly.degree() < 1 ||
        prime < 2 ||
        prime > std::numeric_limits<std::int64_t>::max() / prime ||
        !lmmp_is_prime_ulong_(static_cast<ulong>(prime))) {
        return false;
    }

    auto reduced = bk_reduce_to_mod_coeffs(poly, prime);
    if (static_cast<int>(reduced.size()) - 1 != poly.degree()) {
        return false;
    }
    return bk_is_square_free_mod(reduced, prime);
}

/**
 * @brief 从小素数表中选取适合 Berlekamp 分解的素数。
 * @param[in] poly 非零有理系数多项式，次数 >= 1。
 * @return 首个合适的素数 p；候选均不满足时返回 -1。
 * @internal
 */
int64_t bk_select_prime(const Polynomial<Rational>& poly) {
    if (poly.is_zero() || poly.degree() < 1) {
        return -1;
    }
    for (size_t i = 0; i < BK_NUM_SMALL_PRIMES; ++i) {
        if (bk_is_suitable_prime(poly, BK_SMALL_PRIMES[i])) {
            return BK_SMALL_PRIMES[i];
        }
    }
    return -1;
}

}
