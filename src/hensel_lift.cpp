/**
 * @file hensel_lift.cpp
 * @brief Hensel 提升算法实现:Mignotte 界计算与二次提升.
 *
 * 本文件实现 Phase 4(Hensel 提升)的核心逻辑:
 * - Mignotte 界计算:确定提升高度 k 使得 p^k 足以恢复真因子系数
 * - 二次 Hensel 提升:将模 p 因子逐步提升到 mod p^k
 * - 多因子提升:通过二叉树结构递归配对提升
 *
 * @see Mignotte, M. "An inequality about factors of polynomials."
 *      Mathematics of Computation, 28(128), 1974.
 * @see Zassenhaus, H. "On Hensel factorization, I."
 *      Journal of Number Theory, 1(3), 1969.
 */

#include "transcendental_factor.hpp"

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <exception>

namespace lamina {


/**
 * @brief 计算二项式系数 C(n, k).
 *
 * 使用乘法公式 C(n, k) = n! / (k! * (n-k)!) 逐步计算;
 * BigInt 承载全部中间结果,k = min(k, n-k) 缩短迭代路径.
 *
 * @param[in] n 上标(非负整数)
 * @param[in] k 下标(0 <= k <= n)
 * @return C(n, k) 的精确值
 * @internal
 */
static BigInt hl_binomial(int n, int k) {
    if (k < 0 || k > n) return BigInt(0);
    if (k == 0 || k == n) return BigInt(1);

    /// 利用对称性 C(n, k) = C(n, n-k)
    if (k > n - k) {
        k = n - k;
    }

    BigInt result(1);
    for (int i = 0; i < k; ++i) {
        result = result * BigInt(n - i);
        result = result / BigInt(i + 1);
    }
    return result;
}

/**
 * @brief 计算多项式系数的 L2 范数的平方.
 *
 * 即 ||f||_2^2 = sum(a_i^2),其中 a_i 为多项式各项系数.
 *
 * @param[in] poly 整系数多项式
 * @return 系数平方和(BigInt 精确值)
 * @internal
 */
static BigInt hl_l2_norm_squared(const Polynomial<BigInt>& poly) {
    BigInt sum(0);
    for (const auto& c : poly.coeffs) {
        sum = sum + c * c;
    }
    return sum;
}

/**
 * @brief 计算多项式的 Mignotte 界.
 *
 * 对于次数为 n 的多项式 f(x),其任意因子的系数绝对值上界为:
 *   B = C(n, floor(n/2)) * ||f||_2
 *
 * 其中 C(n, k) 为二项式系数,||f||_2 为系数向量的 L2 范数.
 *
 * 实现使用 B^2 = C(n, floor(n/2))^2 * ||f||_2^2 在整数域内计算,
 * 最后由 BigInt::sqrt() 求向上取整的整数平方根.
 *
 * @param[in] poly 整系数多项式(非零,次数 >= 1)
 * @return Mignotte 界 B
 *
 * @see Mignotte, M. "An inequality about factors of polynomials."
 *      Mathematics of Computation, 28(128), 1974.
 * @internal
 */
static BigInt hl_mignotte_bound(const Polynomial<BigInt>& poly) {
    int n = poly.degree();
    if (n <= 0) return BigInt(1);

    /// C(n, floor(n/2))
    BigInt binom = hl_binomial(n, n / 2);

    /// ||f||_2^2
    BigInt norm_sq = hl_l2_norm_squared(poly);

    /// B^2 = binom^2 * norm_sq
    BigInt b_squared = binom * binom * norm_sq;

    /// B = ceil(sqrt(B^2))
    /// BigInt::sqrt() 返回 floor(sqrt(x)),需要检查是否精确
    BigInt b = b_squared.sqrt();
    if (b * b < b_squared) {
        b = b + BigInt(1);
    }

    return b;
}

/**
 * @brief 计算 Hensel 提升所需的高度 k.
 *
 * 提升高度 k 需满足 p^k > 2 * B * |lc(f)|,其中:
 * - B 为 Mignotte 界
 * - lc(f) 为多项式首项系数
 * - p 为选定的素数
 *
 * 通过连续乘以 p 直至超过阈值确定 k,使界计算保持在整数域内.
 *
 * @param[in] poly  整系数多项式(非零,次数 >= 1)
 * @param[in] prime 选定的素数 p
 * @return 所需的提升高度 k(正整数)
 * @internal
 */
static int hl_compute_lift_height(const Polynomial<BigInt>& poly, int64_t prime) {
    if (poly.is_zero() || poly.degree() <= 0) return 1;

    BigInt B = hl_mignotte_bound(poly);
    BigInt lc = poly.lead_coeff().Abs();

    /// 阈值 threshold = 2 * B * |lc(f)|
    BigInt threshold = BigInt(2) * B * lc;

    /// 计算最小 k 使得 p^k > threshold
    BigInt p_power(1);
    BigInt big_p(static_cast<long long>(prime));
    int k = 0;

    while (p_power <= threshold) {
        p_power = p_power * big_p;
        k++;
    }

    /// 至少提升 1 次
    if (k < 1) k = 1;

    return k;
}

/**
 * @brief 将 BigInt 系数归约到对称表示 [-m/2, m/2).
 *
 * 对给定系数 c 和模数 m,计算 c mod m 并映射到 [-m/2, m/2) 区间.
 * 这确保提升后的系数保持最小绝对值表示.
 *
 * @param[in] c 待归约的系数
 * @param[in] m 模数(正整数)
 * @return 对称表示下的归约值
 * @internal
 */
static BigInt hl_symmetric_mod(const BigInt& c, const BigInt& m) {
    if (m.is_zero()) return c;

    BigInt r = c % m;
    /// 确保 r 在 [0, m) 范围内
    if (r.IsNegative()) {
        r = r + m;
    }

    /// 映射到 [-m/2, m/2)
    BigInt half_m = m / BigInt(2);
    if (r > half_m) {
        r = r - m;
    }
    return r;
}

/**
 * @brief 对多项式系数向量执行模归约(对称表示).
 *
 * 将向量中每个系数归约到 [-m/2, m/2) 区间,并去除高次零系数.
 *
 * @param[in,out] poly 系数向量(就地修改)
 * @param[in] m 模数
 * @internal
 */
static void hl_reduce_coeffs(std::vector<BigInt>& poly, const BigInt& m) {
    for (auto& c : poly) {
        c = hl_symmetric_mod(c, m);
    }
    /// 去除高次零系数
    while (!poly.empty() && poly.back().is_zero()) {
        poly.pop_back();
    }
}

/**
 * @brief 多项式乘法(BigInt 系数),结果模 m 归约.
 *
 * 计算两个系数向量表示的多项式之积,所有系数归约到对称表示 [-m/2, m/2).
 *
 * @param[in] a 第一个多项式的系数向量
 * @param[in] b 第二个多项式的系数向量
 * @param[in] m 模数
 * @return 乘积多项式的系数向量
 * @internal
 */
static std::vector<BigInt> hl_poly_mul_mod(const std::vector<BigInt>& a,
                                           const std::vector<BigInt>& b,
                                           const BigInt& m) {
    if (a.empty() || b.empty()) return {};

    size_t n = a.size() + b.size() - 1;
    std::vector<BigInt> result(n, BigInt(0));

    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].is_zero()) continue;
        for (size_t j = 0; j < b.size(); ++j) {
            if (b[j].is_zero()) continue;
            result[i + j] = result[i + j] + a[i] * b[j];
        }
    }

    hl_reduce_coeffs(result, m);
    return result;
}

/**
 * @brief 多项式减法(BigInt 系数),结果模 m 归约.
 *
 * @param[in] a 被减数多项式系数向量
 * @param[in] b 减数多项式系数向量
 * @param[in] m 模数
 * @return 差多项式的系数向量
 * @internal
 */
static std::vector<BigInt> hl_poly_sub_mod(const std::vector<BigInt>& a,
                                           const std::vector<BigInt>& b,
                                           const BigInt& m) {
    size_t n = std::max(a.size(), b.size());
    std::vector<BigInt> result(n, BigInt(0));

    for (size_t i = 0; i < n; ++i) {
        BigInt ai = (i < a.size()) ? a[i] : BigInt(0);
        BigInt bi = (i < b.size()) ? b[i] : BigInt(0);
        result[i] = ai - bi;
    }

    hl_reduce_coeffs(result, m);
    return result;
}

/**
 * @brief 多项式加法(BigInt 系数),结果模 m 归约.
 *
 * @param[in] a 第一个多项式系数向量
 * @param[in] b 第二个多项式系数向量
 * @param[in] m 模数
 * @return 和多项式的系数向量
 * @internal
 */
static std::vector<BigInt> hl_poly_add_mod(const std::vector<BigInt>& a,
                                           const std::vector<BigInt>& b,
                                           const BigInt& m) {
    size_t n = std::max(a.size(), b.size());
    std::vector<BigInt> result(n, BigInt(0));

    for (size_t i = 0; i < n; ++i) {
        BigInt ai = (i < a.size()) ? a[i] : BigInt(0);
        BigInt bi = (i < b.size()) ? b[i] : BigInt(0);
        result[i] = ai + bi;
    }

    hl_reduce_coeffs(result, m);
    return result;
}

/**
 * @brief 计算 BigInt 在模 m 下的乘法逆元.
 *
 * 使用扩展欧几里得算法求 a 在 Z/mZ 中的逆元.
 * 要求 gcd(a, m) = 1.
 *
 * @param[in] a 待求逆的整数
 * @param[in] m 模数
 * @return a 的模逆元,归约到 [0, m)
 * @throw std::domain_error 若 a 不可逆
 * @internal
 */
static BigInt hl_mod_inverse(const BigInt& a, const BigInt& m) {
    /// 扩展欧几里得:求 s 使得 a*s == 1 (mod m)
    BigInt r0 = m, r1 = a % m;
    if (r1.IsNegative()) r1 = r1 + m;

    BigInt s0(0), s1(1);

    while (!r1.is_zero()) {
        BigInt q = r0 / r1;
        BigInt r_new = r0 - q * r1;
        BigInt s_new = s0 - q * s1;

        r0 = r1;
        r1 = r_new;
        s0 = s1;
        s1 = s_new;
    }

    /// r0 应为 gcd(a, m) = 1
    if (r0 != BigInt(1) && r0 != BigInt(-1)) {
        throw std::domain_error("hl_mod_inverse: element not invertible");
    }

    /// 若 gcd 为 -1,调整符号
    if (r0.IsNegative()) {
        s0 = s0.negate();
    }

    BigInt result = s0 % m;
    if (result.IsNegative()) result = result + m;
    return result;
}

/**
 * @brief 多项式带余除法(BigInt 系数),模 m 下运算.
 *
 * 计算 a(x) = q(x) * b(x) + r(x),其中 deg(r) < deg(b).
 * 要求 b 的首项系数在 Z/mZ 中可逆.
 * 所有系数归约到对称表示.
 *
 * @param[in] a 被除数多项式系数向量
 * @param[in] b 除数多项式系数向量
 * @param[in] m 模数
 * @return pair(商, 余数) 的系数向量
 * @internal
 */
static std::pair<std::vector<BigInt>, std::vector<BigInt>>
hl_poly_divmod(const std::vector<BigInt>& a,
               const std::vector<BigInt>& b,
               const BigInt& m) {
    if (b.empty()) {
        throw std::runtime_error("hl_poly_divmod: division by zero polynomial");
    }

    int deg_a = static_cast<int>(a.size()) - 1;
    int deg_b = static_cast<int>(b.size()) - 1;

    if (deg_a < deg_b) {
        std::vector<BigInt> r = a;
        hl_reduce_coeffs(r, m);
        return {{}, r};
    }

    /// 首项系数的逆元
    BigInt lc_inv = hl_mod_inverse(b.back(), m);

    std::vector<BigInt> remainder = a;
    std::vector<BigInt> quotient(deg_a - deg_b + 1, BigInt(0));

    for (int i = deg_a; i >= deg_b; --i) {
        BigInt coeff = hl_symmetric_mod(remainder[i], m);
        if (coeff.is_zero()) continue;

        BigInt factor = hl_symmetric_mod(coeff * lc_inv, m);
        quotient[i - deg_b] = factor;

        for (int j = 0; j <= deg_b; ++j) {
            remainder[i - deg_b + j] = remainder[i - deg_b + j] - factor * b[j];
            remainder[i - deg_b + j] = hl_symmetric_mod(remainder[i - deg_b + j], m);
        }
    }

    /// 余数截断到 deg < deg_b
    std::vector<BigInt> rem(remainder.begin(), remainder.begin() + deg_b);
    hl_reduce_coeffs(quotient, m);
    hl_reduce_coeffs(rem, m);
    return {quotient, rem};
}

static bool hl_is_prime(int64_t value) {
    if (value < 2) return false;
    if (value == 2) return true;
    if (value % 2 == 0) return false;
    for (int64_t divisor = 3; divisor <= value / divisor; divisor += 2) {
        if (value % divisor == 0) return false;
    }
    return true;
}

static bool hl_coeff_vectors_equal(std::vector<BigInt> left,
                                   std::vector<BigInt> right) {
    while (!left.empty() && left.back().is_zero()) left.pop_back();
    while (!right.empty() && right.back().is_zero()) right.pop_back();
    const std::size_t size = std::max(left.size(), right.size());
    left.resize(size, BigInt(0));
    right.resize(size, BigInt(0));
    return left == right;
}

static Result<void> hl_validate_lift_inputs(
    const Polynomial<BigInt>& poly,
    const std::vector<std::vector<BigInt>>& factor_vecs,
    int64_t prime) {
    constexpr const char* operation = "hensel_lift";
    if (!hl_is_prime(prime)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "Hensel lifting requires a prime modulus",
                                     operation);
    }
    if (poly.is_zero() || poly.degree() <= 0) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "Hensel lifting requires a non-constant polynomial",
                                     operation);
    }
    if (factor_vecs.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "Hensel lifting requires at least one mod-p factor",
                                     operation);
    }
    for (const auto& factor : factor_vecs) {
        if (factor.empty()) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "Hensel lifting factor cannot be zero",
                                         operation);
        }
    }

    const BigInt modulus(static_cast<long long>(prime));
    std::vector<BigInt> product{BigInt(1)};
    for (const auto& factor : factor_vecs) {
        product = hl_poly_mul_mod(product, factor, modulus);
    }
    std::vector<BigInt> reduced_f = poly.coeffs;
    hl_reduce_coeffs(reduced_f, modulus);
    if (!hl_coeff_vectors_equal(product, reduced_f)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "mod-p factors do not multiply to the input polynomial",
                                     operation);
    }
    return Result<void>::success();
}

/// HenselLiftPair 结构体已在 transcendental_factor.hpp 中声明

/**
 * @brief 执行一步二次 Hensel 提升:mod m -> mod m^2.
 *
 * 给定 f == g*h (mod m) 且 s*g + t*h == 1 (mod m),
 * 计算 g', h', s', t' 使得 f == g'*h' (mod m^2) 且 s'*g' + t'*h' == 1 (mod m^2).
 *
 * 算法步骤:
 * 1. 计算误差 e = f - g*h(所有系数应被 m 整除)
 * 2. 计算修正项:(q, r) = divmod(s*e, h) over Z/m^2Z
 *    - h' = h + r (mod m^2)
 *    - g' = g + t*e + q*g (mod m^2)
 * 3. 更新 Bezout 系数:
 *    - b = s*g' + t*h' - 1 (mod m^2)
 *    - (c, d) = divmod(s*b, h') over Z/m^2Z
 *    - s' = s - d (mod m^2)
 *    - t' = t - t*b - c*g' (mod m^2)
 *
 * @param[in] f       原始多项式系数向量
 * @param[in] current 当前提升状态(g, h, s, t, modulus=m)
 * @return 提升后的状态(g', h', s', t', modulus=m^2)
 *
 * @pre f == g*h (mod m)
 * @pre s*g + t*h == 1 (mod m)
 * @pre deg(s) < deg(h), deg(t) < deg(g)
 *
 * @see Zassenhaus, H. "On Hensel factorization, I."
 *      Journal of Number Theory, 1(3), 1969.
 * @internal
 */
HenselLiftPair hl_two_factor_lift(
    const std::vector<BigInt>& f,
    const HenselLiftPair& current) {

    const auto& g = current.g;
    const auto& h = current.h;
    const auto& s = current.s;
    const auto& t = current.t;
    const BigInt& m = current.modulus;

    /// 新模数 m^2 = m * m
    BigInt m2 = m * m;

    /// 步骤 1:计算误差 e = f - g*h (mod m^2)
    std::vector<BigInt> gh = hl_poly_mul_mod(g, h, m2);
    std::vector<BigInt> e = hl_poly_sub_mod(f, gh, m2);

    /// 步骤 2:计算修正项
    /// se = s * e (mod m^2)
    std::vector<BigInt> se = hl_poly_mul_mod(s, e, m2);

    /// (q, r) = divmod(se, h) over Z/m^2Z,其中 deg(r) < deg(h)
    auto [q, r] = hl_poly_divmod(se, h, m2);

    /// g' = g + t*e + q*g (mod m^2)
    std::vector<BigInt> te = hl_poly_mul_mod(t, e, m2);
    std::vector<BigInt> qg = hl_poly_mul_mod(q, g, m2);
    std::vector<BigInt> g_new = hl_poly_add_mod(g, hl_poly_add_mod(te, qg, m2), m2);
    hl_reduce_coeffs(g_new, m2);

    /// h' = h + r (mod m^2)
    std::vector<BigInt> h_new = hl_poly_add_mod(h, r, m2);
    hl_reduce_coeffs(h_new, m2);

    /// 步骤 3:更新 Bezout 系数
    /// b = s*g' + t*h' - 1 (mod m^2)
    std::vector<BigInt> sg_new = hl_poly_mul_mod(s, g_new, m2);
    std::vector<BigInt> th_new = hl_poly_mul_mod(t, h_new, m2);
    std::vector<BigInt> sg_plus_th = hl_poly_add_mod(sg_new, th_new, m2);

    /// 减去常数 1
    std::vector<BigInt> one_poly = {BigInt(1)};
    std::vector<BigInt> b = hl_poly_sub_mod(sg_plus_th, one_poly, m2);

    /// (c, d) = divmod(s*b, h') over Z/m^2Z
    std::vector<BigInt> sb = hl_poly_mul_mod(s, b, m2);
    auto [c, d] = hl_poly_divmod(sb, h_new, m2);

    /// s' = s - d (mod m^2)
    std::vector<BigInt> s_new = hl_poly_sub_mod(s, d, m2);
    hl_reduce_coeffs(s_new, m2);

    /// t' = t - t*b - c*g' (mod m^2)
    std::vector<BigInt> tb = hl_poly_mul_mod(t, b, m2);
    std::vector<BigInt> cg_new = hl_poly_mul_mod(c, g_new, m2);
    std::vector<BigInt> t_new = hl_poly_sub_mod(t, hl_poly_add_mod(tb, cg_new, m2), m2);
    hl_reduce_coeffs(t_new, m2);

    return HenselLiftPair{g_new, h_new, s_new, t_new, m2};
}

/**
 * @brief 计算两个多项式在模 m 下的扩展 GCD 的 Bezout 系数.
 *
 * 给定互素多项式 g, h(mod m),求 s, t 使得 s*g + t*h == 1 (mod m),
 * 且 deg(s) < deg(h), deg(t) < deg(g).
 *
 * 使用扩展欧几里得算法在 Z/mZ[x] 上运算.
 *
 * @param[in] g 第一个多项式系数向量
 * @param[in] h 第二个多项式系数向量
 * @param[in] m 模数
 * @return pair(s, t) 满足 s*g + t*h == 1 (mod m)
 * @internal
 */
static std::pair<std::vector<BigInt>, std::vector<BigInt>>
hl_extended_gcd_poly(const std::vector<BigInt>& g,
                     const std::vector<BigInt>& h,
                     const BigInt& m) {
    /// 扩展欧几里得算法:r0 = g, r1 = h
    /// s0*g + t0*h = r0, s1*g + t1*h = r1
    std::vector<BigInt> r0 = g, r1 = h;
    std::vector<BigInt> s0 = {BigInt(1)}, s1 = {};  // s0=1, s1=0
    std::vector<BigInt> t0 = {}, t1 = {BigInt(1)};  // t0=0, t1=1

    while (!r1.empty()) {
        auto [q, r] = hl_poly_divmod(r0, r1, m);

        /// r_new = r0 - q*r1
        std::vector<BigInt> r_new = r;

        /// s_new = s0 - q*s1
        std::vector<BigInt> qs1 = hl_poly_mul_mod(q, s1, m);
        std::vector<BigInt> s_new = hl_poly_sub_mod(s0, qs1, m);

        /// t_new = t0 - q*t1
        std::vector<BigInt> qt1 = hl_poly_mul_mod(q, t1, m);
        std::vector<BigInt> t_new = hl_poly_sub_mod(t0, qt1, m);

        r0 = r1; r1 = r_new;
        s0 = s1; s1 = s_new;
        t0 = t1; t1 = t_new;
    }

    /// r0 = gcd,应为常数(可逆元)
    /// 归一化使 gcd = 1
    if (!r0.empty() && r0[0] != BigInt(1)) {
        BigInt inv = hl_mod_inverse(r0[0], m);
        for (auto& c : s0) {
            c = hl_symmetric_mod(c * inv, m);
        }
        for (auto& c : t0) {
            c = hl_symmetric_mod(c * inv, m);
        }
    }

    /// 去除高次零系数
    while (!s0.empty() && s0.back().is_zero()) s0.pop_back();
    while (!t0.empty() && t0.back().is_zero()) t0.pop_back();

    return {s0, t0};
}

/**
 * @brief 多因子 Hensel 提升:将 k 个模 p 因子同时提升到 mod p^target_k.
 *
 * 使用顺序剥离策略:
 * 1. 将 f 分解为 g_1 和 h_1 = g_2*g₃*...*gₖ
 * 2. 对 (g_1, h_1) 执行二因子提升
 * 3. 递归对 h_1 继续剥离下一个因子
 *
 * 每步提升通过反复调用 hl_two_factor_lift 实现二次提升(mod p -> p^2 -> p^4 -> ...),
 * 直到模数达到或超过 p^target_k.
 *
 * @param[in] f         原始多项式系数向量(升幂排列)
 * @param[in] factors   模 p 下的因子列表(各因子系数向量)
 * @param[in] prime     素数 p
 * @param[in] target_k  目标提升高度(提升到 mod p^target_k)
 * @return 提升后的因子列表(各因子系数在对称表示下)
 *
 * @pre f == factors[0] * factors[1] * ... * factors[k-1] (mod p)
 * @pre 各因子两两互素 (mod p)
 *
 * @see Zassenhaus, H. "On Hensel factorization, I."
 *      Journal of Number Theory, 1(3), 1969.
 * @internal
 */
static std::vector<std::vector<BigInt>> hl_multi_factor_lift(
    const std::vector<BigInt>& f,
    const std::vector<std::vector<BigInt>>& factors,
    int64_t prime,
    int target_k) {

    if (factors.empty()) return {};
    if (factors.size() == 1) {
        /// 单因子:直接归约到目标模数
        BigInt target_mod(1);
        BigInt big_p(static_cast<long long>(prime));
        for (int i = 0; i < target_k; ++i) {
            target_mod = target_mod * big_p;
        }
        std::vector<BigInt> result = f;
        hl_reduce_coeffs(result, target_mod);
        return {result};
    }

    BigInt big_p(static_cast<long long>(prime));
    BigInt initial_mod(static_cast<long long>(prime));

    /// 计算目标模数 p^target_k
    BigInt target_mod(1);
    for (int i = 0; i < target_k; ++i) {
        target_mod = target_mod * big_p;
    }

    /// 顺序剥离策略:逐个提升因子
    std::vector<std::vector<BigInt>> lifted_factors;
    lifted_factors.reserve(factors.size());

    /// 当前待分解的多项式(初始为 f)
    std::vector<BigInt> remaining = f;
    hl_reduce_coeffs(remaining, target_mod);

    /// 剩余因子列表
    std::vector<std::vector<BigInt>> remaining_factors = factors;

    for (size_t i = 0; i < factors.size() - 1; ++i) {
        /// 取出第一个因子 g
        std::vector<BigInt> g = remaining_factors[0];

        /// 计算 h = 剩余因子的乘积 (mod p)
        std::vector<BigInt> h = remaining_factors[1];
        for (size_t j = 2; j < remaining_factors.size(); ++j) {
            h = hl_poly_mul_mod(h, remaining_factors[j], initial_mod);
        }

        /// 计算 Bezout 系数 s, t 使得 s*g + t*h == 1 (mod p)
        auto [s, t] = hl_extended_gcd_poly(g, h, initial_mod);

        /// 构造初始提升状态
        HenselLiftPair state{g, h, s, t, initial_mod};

        /// 反复二次提升直到模数 >= target_mod
        while (state.modulus < target_mod) {
            state = hl_two_factor_lift(remaining, state);
        }

        /// 提取提升后的 g(已提升到 target_mod)
        std::vector<BigInt> lifted_g = state.g;
        hl_reduce_coeffs(lifted_g, target_mod);
        lifted_factors.push_back(lifted_g);

        /// 更新 remaining 为提升后的 h
        remaining = state.h;
        hl_reduce_coeffs(remaining, target_mod);

        /// 更新 remaining_factors:去掉第一个,后续因子保持不变
        remaining_factors.erase(remaining_factors.begin());
    }

    /// 最后一个因子就是 remaining
    hl_reduce_coeffs(remaining, target_mod);
    lifted_factors.push_back(remaining);

    return lifted_factors;
}


/**
 * @brief Hensel 提升:将模 p 因子提升到 mod p^k.
 *
 * 将 Berlekamp 分解得到的模 p 不可约因子通过二次 Hensel 提升
 * 升至 mod p^k,其中 k 由 Mignotte 界确定.
 *
 * 流程:
 * 1. 将 Polynomial<ModInt> 因子转换为 vector<BigInt> 内部表示
 * 2. 调用 hl_multi_factor_lift 执行多因子提升
 * 3. 将结果转换回 Polynomial<BigInt>
 *
 * @param[in] poly        整系数多项式
 * @param[in] mod_factors 模 p 下的不可约因子
 * @param[in] prime       素数 p
 * @param[in] lift_bound  提升次数上界 k(若为 0 则自动计算)
 * @return 提升后的整系数因子列表(系数在对称表示 [-p^k/2, p^k/2] 下)
 */
HenselLiftResult hensel_lift_checked(
    const Polynomial<BigInt>& poly,
    const std::vector<Polynomial<ModInt>>& mod_factors,
    int64_t prime,
    int lift_bound) {
    constexpr const char* operation = "hensel_lift";

    /// 将原始多项式转换为系数向量
    std::vector<BigInt> f_vec = poly.coeffs;

    /// 将 Polynomial<ModInt> 因子转换为 vector<BigInt> 表示
    std::vector<std::vector<BigInt>> factor_vecs;
    factor_vecs.reserve(mod_factors.size());
    for (const auto& mf : mod_factors) {
        std::vector<BigInt> fv;
        fv.reserve(mf.coeffs.size());
        for (const auto& c : mf.coeffs) {
            fv.push_back(BigInt(static_cast<long long>(c.value())));
        }
        factor_vecs.push_back(std::move(fv));
    }

    auto valid = hl_validate_lift_inputs(poly, factor_vecs, prime);
    if (!valid) return HenselLiftResult::failure(valid.error());

    /// 若调用方未指定提升界,则自动计算
    int k = lift_bound;
    if (k <= 0) {
        k = hl_compute_lift_height(poly, prime);
    }

    /// 执行多因子 Hensel 提升
    std::vector<std::vector<BigInt>> lifted_vecs;
    try {
        lifted_vecs = hl_multi_factor_lift(f_vec, factor_vecs, prime, k);
    } catch (const std::exception& ex) {
        return HenselLiftResult::failure(CasErrc::InternalInvariant,
                                         ex.what(),
                                         operation);
    }

    /// 将结果转换回 Polynomial<BigInt>
    std::vector<Polynomial<BigInt>> result;
    result.reserve(lifted_vecs.size());
    for (auto& lv : lifted_vecs) {
        result.emplace_back(std::move(lv), poly.variable_name);
    }

    return HenselLiftResult::success(std::move(result));
}


} // namespace lamina
