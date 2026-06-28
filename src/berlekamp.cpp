/**
 * @file berlekamp.cpp
 * @brief Berlekamp 模分解算法实现：素数选取与模多项式因式分解。
 *
 * 本文件实现 Phase 3（Berlekamp 模分解）的核心逻辑：
 * - 从小素数表中选取满足条件的素数 p
 * - 将有理系数多项式约化为 F_p 上的多项式
 * - 执行 Berlekamp 算法求模不可约因子（矩阵构造与分裂在后续任务中实现）
 *
 * @see Berlekamp, E.R. "Factoring polynomials over finite fields."
 *      Bell System Technical Journal, 46(8), 1967.
 */

#include "transcendental_factor.hpp"
#include "poly_utils.hpp"

#include <vector>
#include <cstdint>
#include <algorithm>

namespace lamina {

// ============================================================
/// 文件局部常量与辅助函数
// ============================================================

/**
 * @brief Berlekamp 算法使用的小素数表。
 *
 * 用于模分解时的素数选取。从最小素数开始迭代，
 * 选取第一个满足条件（不整除首项系数且模约化后仍 square-free）的素数。
 *
 * @internal
 */
static constexpr int64_t BK_SMALL_PRIMES[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229
};

/** @brief 小素数表的长度 */
static constexpr size_t BK_NUM_SMALL_PRIMES =
    sizeof(BK_SMALL_PRIMES) / sizeof(BK_SMALL_PRIMES[0]);

/**
 * @brief 将有理系数多项式约化为模 p 的整数系数向量。
 *
 * 算法：
 *   1. 清除分母：计算所有系数分母的 LCM，乘以 LCM 得到整系数多项式
 *   2. 对每个整系数取模 p
 *
 * @param[in] poly 有理系数多项式
 * @param[in] p    素数模数
 * @return 模 p 下的系数向量（升幂排列），去除高次零系数
 * @internal
 */
static std::vector<int64_t> bk_reduce_to_mod_coeffs(
    const Polynomial<Rational>& poly,
    int64_t p) {

    if (poly.is_zero()) {
        return {};
    }

    /// 计算所有分母的 LCM 以清除分母
    BigInt lcm_den(1);
    for (const auto& c : poly.coeffs) {
        BigInt den = c.get_denominator();
        BigInt g = BigInt::gcd(lcm_den, den);
        lcm_den = lcm_den * (den / g);
    }

    BigInt big_p(static_cast<long long>(p));

    /// 将每个系数乘以 lcm_den 得到整系数，再取模 p
    std::vector<int64_t> mod_coeffs;
    mod_coeffs.reserve(poly.coeffs.size());

    for (const auto& c : poly.coeffs) {
        /// integer_coeff = c * lcm_den = (num/den) * lcm_den
        BigInt num = c.get_numerator();
        BigInt den = c.get_denominator();
        BigInt integer_coeff = num * (lcm_den / den);

        /// 取模 p
        BigInt rem = integer_coeff % big_p;
        int64_t rem_val = rem.to_int();
        /// 确保非负
        if (rem_val < 0) rem_val += p;

        mod_coeffs.push_back(rem_val);
    }

    /// 去除高次零系数
    while (!mod_coeffs.empty() && mod_coeffs.back() == 0) {
        mod_coeffs.pop_back();
    }

    return mod_coeffs;
}

/**
 * @brief 计算模 p 多项式的形式导数。
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
        int64_t coeff = (f[i] * static_cast<int64_t>(i)) % p;
        deriv.push_back(coeff);
    }

    /// 去除高次零系数
    while (!deriv.empty() && deriv.back() == 0) {
        deriv.pop_back();
    }

    return deriv;
}

/**
 * @brief 模 p 下求乘法逆元。
 *
 * 使用扩展欧几里得算法计算 a^(-1) mod p。
 *
 * @param[in] a 待求逆的值（非零）
 * @param[in] p 素数模数
 * @return a 在模 p 下的逆元
 * @internal
 */
static int64_t bk_mod_inverse(int64_t a, int64_t p) {
    int64_t s, t;
    extended_gcd(((a % p) + p) % p, p, s, t);
    return ((s % p) + p) % p;
}

/**
 * @brief 模 p 多项式带余除法。
 *
 * 计算 a / b 的商和余数，使得 a = q * b + r，deg(r) < deg(b)。
 *
 * @param[in]  a 被除数系数向量
 * @param[in]  b 除数系数向量（非空）
 * @param[in]  p 素数模数
 * @param[out] q 商系数向量
 * @param[out] r 余数系数向量
 * @internal
 */
static void bk_div_mod(
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

        /// 更新 deg_r
        while (deg_r >= 0 && r[deg_r] == 0) {
            r.pop_back();
            deg_r--;
        }
    }
}

/**
 * @brief 计算两个模 p 多项式的最大公因式。
 *
 * 使用欧几里得算法。结果为首一多项式。
 *
 * @param[in] a 第一个多项式系数向量
 * @param[in] b 第二个多项式系数向量
 * @param[in] p 素数模数
 * @return gcd(a, b) 的首一化系数向量
 * @internal
 */
static std::vector<int64_t> bk_gcd_mod(
    std::vector<int64_t> a,
    std::vector<int64_t> b,
    int64_t p) {

    while (!b.empty()) {
        std::vector<int64_t> q, r;
        bk_div_mod(a, b, p, q, r);
        a = std::move(b);
        b = std::move(r);
    }

    /// 首一化
    if (!a.empty()) {
        int64_t lc_inv = bk_mod_inverse(a.back(), p);
        for (auto& c : a) {
            c = (c * lc_inv) % p;
        }
    }

    return a;
}

/**
 * @brief 检查模 p 多项式是否为 square-free。
 *
 * 多项式 f 在 F_p 上 square-free 当且仅当 gcd(f, f') 的次数为 0。
 *
 * @param[in] f 模 p 系数向量
 * @param[in] p 素数模数
 * @return f 为 square-free 返回 true
 * @internal
 */
static bool bk_is_square_free_mod(
    const std::vector<int64_t>& f,
    int64_t p) {

    if (f.size() <= 1) return true;

    std::vector<int64_t> f_prime = bk_derivative_mod(f, p);

    /// 若导数为零多项式，则 f 不是 square-free（特征 p 下的情形）
    if (f_prime.empty()) return false;

    std::vector<int64_t> g = bk_gcd_mod(f, f_prime, p);

    /// gcd 次数为 0（即 g 为非零常数）意味着 square-free
    return g.size() <= 1;
}

/**
 * @brief 从小素数表中选取适合 Berlekamp 分解的素数。
 *
 * 遍历 BK_SMALL_PRIMES 表，对每个素数 p 检查：
 *   1. p 不整除多项式的首项系数
 *   2. 多项式模 p 约化后仍为 square-free（gcd(f mod p, f' mod p) 次数为 0）
 *
 * 返回第一个满足条件的素数；若所有候选素数均不满足，返回 -1。
 *
 * @param[in] poly 有理系数多项式（假设非零且次数 ≥ 1）
 * @return 选定的素数 p；无合适素数时返回 -1
 * @internal
 */
static int64_t bk_select_prime(const Polynomial<Rational>& poly) {
    if (poly.is_zero() || poly.degree() < 1) {
        return -1;
    }

    /// 提取首项系数的分子绝对值
    Rational lc = poly.lead_coeff();
    BigInt lc_num = lc.get_numerator();
    if (lc_num.IsNegative()) {
        lc_num = lc_num.Abs();
    }

    for (size_t i = 0; i < BK_NUM_SMALL_PRIMES; ++i) {
        int64_t p = BK_SMALL_PRIMES[i];
        BigInt big_p(static_cast<long long>(p));

        /// 条件 1：p 不整除首项系数的分子
        BigInt rem = lc_num % big_p;
        if (rem == BigInt(0)) {
            continue;
        }

        /// 条件 2：多项式模 p 约化后仍为 square-free
        std::vector<int64_t> f_mod_p = bk_reduce_to_mod_coeffs(poly, p);

        /// 约化后次数降低说明清除分母后首项系数被 p 整除
        if (static_cast<int>(f_mod_p.size()) - 1 < poly.degree()) {
            continue;
        }

        if (bk_is_square_free_mod(f_mod_p, p)) {
            return p;
        }
    }

    return -1;
}

/**
 * @brief 模 p 多项式乘法并对 f 取余。
 *
 * 计算 a(x) * b(x) mod f(x)，所有系数运算在 F_p 上进行。
 *
 * @param[in] a 第一个多项式系数向量（升幂排列）
 * @param[in] b 第二个多项式系数向量（升幂排列）
 * @param[in] f 模多项式系数向量（首一，升幂排列）
 * @param[in] p 素数模数
 * @return a * b mod f 的系数向量
 * @internal
 */
static std::vector<int64_t> bk_poly_mul_mod(
    const std::vector<int64_t>& a,
    const std::vector<int64_t>& b,
    const std::vector<int64_t>& f,
    int64_t p) {

    if (a.empty() || b.empty()) {
        return {};
    }

    /// 计算 a * b（卷积）
    size_t prod_size = a.size() + b.size() - 1;
    std::vector<int64_t> product(prod_size, 0);

    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] == 0) continue;
        for (size_t j = 0; j < b.size(); ++j) {
            product[i + j] = (product[i + j] + a[i] * b[j]) % p;
            if (product[i + j] < 0) product[i + j] += p;
        }
    }

    /// 对 f 取余
    std::vector<int64_t> q, r;
    bk_div_mod(product, f, p, q, r);

    return r;
}

/**
 * @brief 模 p 多项式幂运算：计算 base^exp mod f(x)。
 *
 * 使用平方-乘法（repeated squaring）算法高效计算多项式幂次模约化。
 *
 * @param[in] base 底多项式系数向量（升幂排列）
 * @param[in] exp  指数（非负整数）
 * @param[in] f    模多项式系数向量（首一，升幂排列）
 * @param[in] p    素数模数
 * @return base^exp mod f 的系数向量
 *
 * @see Berlekamp, E.R. "Factoring polynomials over finite fields."
 * @internal
 */
static std::vector<int64_t> bk_poly_pow_mod(
    const std::vector<int64_t>& base,
    int64_t exp,
    const std::vector<int64_t>& f,
    int64_t p) {

    /// base^0 = 1
    if (exp == 0) {
        return {1};
    }

    /// base^1 mod f
    if (exp == 1) {
        std::vector<int64_t> q, r;
        bk_div_mod(base, f, p, q, r);
        return r;
    }

    /// Repeated squaring: 从高位到低位扫描 exp 的二进制位
    std::vector<int64_t> result = {1};  // 累积结果，初始为 1
    std::vector<int64_t> cur = base;    // 当前底数

    /// 先对 cur 取余确保 deg(cur) < deg(f)
    {
        std::vector<int64_t> q, r;
        bk_div_mod(cur, f, p, q, r);
        cur = std::move(r);
    }

    int64_t e = exp;
    while (e > 0) {
        if (e & 1) {
            result = bk_poly_mul_mod(result, cur, f, p);
        }
        e >>= 1;
        if (e > 0) {
            cur = bk_poly_mul_mod(cur, cur, f, p);
        }
    }

    return result;
}

/**
 * @brief 构造 Berlekamp 矩阵 Q。
 *
 * 设 f(x) 为 F_p 上次数为 n 的首一多项式。Berlekamp 矩阵 Q 为 n×n 矩阵，
 * 其中第 i 行为 x^(i*p) mod f(x) 的系数向量（i = 0, 1, ..., n-1）。
 *
 * 算法：
 *   1. 计算 h(x) = x^p mod f(x)（使用 repeated squaring）
 *   2. 第 0 行为 x^0 mod f(x) = 1（即 [1, 0, ..., 0]）
 *   3. 第 i 行 = 第 (i-1) 行 * h(x) mod f(x)
 *
 * @param[in] f 首一多项式系数向量（升幂排列，deg(f) = n）
 * @param[in] p 素数模数
 * @return n×n Berlekamp 矩阵（行优先存储）
 *
 * @see Berlekamp, E.R. "Factoring polynomials over finite fields."
 *      Bell System Technical Journal, 46(8), 1967.
 * @internal
 */
static std::vector<std::vector<int64_t>> bk_build_q_matrix(
    const std::vector<int64_t>& f,
    int64_t p) {

    int n = static_cast<int>(f.size()) - 1;  // deg(f)

    std::vector<std::vector<int64_t>> Q(n, std::vector<int64_t>(n, 0));

    if (n <= 0) {
        return Q;
    }

    /// 第 0 行：x^0 mod f(x) = 1
    Q[0][0] = 1;

    if (n == 1) {
        return Q;
    }

    /// 计算 h(x) = x^p mod f(x)
    std::vector<int64_t> x_poly = {0, 1};  // x
    std::vector<int64_t> h = bk_poly_pow_mod(x_poly, p, f, p);

    /// 第 1 行：x^p mod f(x) = h(x)
    /// 将 h 的系数填入第 1 行（不足 n 位补零）
    for (size_t j = 0; j < h.size() && j < static_cast<size_t>(n); ++j) {
        Q[1][j] = h[j];
    }

    /// 第 i 行（i >= 2）：x^(ip) mod f(x) = (x^((i-1)*p) mod f) * h mod f
    std::vector<int64_t> current = h;
    for (int i = 2; i < n; ++i) {
        current = bk_poly_mul_mod(current, h, f, p);
        for (size_t j = 0; j < current.size() && j < static_cast<size_t>(n); ++j) {
            Q[i][j] = current[j];
        }
    }

    return Q;
}

/**
 * @brief 计算 (Q - I)^T 在 F_p 上的零空间（核）。
 *
 * 算法：
 *   1. 构造 M = (Q - I)^T（n×n 矩阵）
 *   2. 对 M 执行列主元高斯消元，化为行阶梯形
 *   3. 识别自由变量（无主元的列）
 *   4. 对每个自由变量，通过回代构造零空间基向量
 *
 * 零空间维度 k 等于 f(x) 在 F_p 上的不可约因子数。
 * 第一个基向量始终为 [1, 0, 0, ..., 0]（对应平凡因子）。
 *
 * @param[in] Q Berlekamp 矩阵（n×n，行优先存储）
 * @param[in] p 素数模数
 * @return 零空间基向量集合，每个向量长度为 n
 *
 * @see Berlekamp, E.R. "Factoring polynomials over finite fields."
 *      Bell System Technical Journal, 46(8), 1967.
 * @internal
 */
static std::vector<std::vector<int64_t>> bk_null_space(
    const std::vector<std::vector<int64_t>>& Q,
    int64_t p) {

    int n = static_cast<int>(Q.size());
    if (n == 0) {
        return {};
    }

    /// 构造 M = (Q - I)^T
    /// M[i][j] = Q[j][i] - (i == j ? 1 : 0)
    std::vector<std::vector<int64_t>> M(n, std::vector<int64_t>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int64_t val = Q[j][i];  // 转置：M[i][j] = Q[j][i]
            if (i == j) {
                val = (val - 1 + p) % p;  // 减去单位矩阵
            }
            M[i][j] = val;
        }
    }

    /// 高斯消元（列主元），记录每行的主元列位置
    /// pivot_col[row] = 该行主元所在列，-1 表示该行为零行
    std::vector<int> pivot_col(n, -1);
    /// pivot_row[col] = 该列主元所在行，-1 表示该列为自由变量
    std::vector<int> pivot_row(n, -1);

    int current_row = 0;
    for (int col = 0; col < n && current_row < n; ++col) {
        /// 寻找列主元（当前列中绝对值最大的非零元素）
        int pivot = -1;
        for (int row = current_row; row < n; ++row) {
            if (M[row][col] != 0) {
                pivot = row;
                break;
            }
        }

        if (pivot == -1) {
            /// 该列无主元，为自由变量
            continue;
        }

        /// 交换行
        if (pivot != current_row) {
            std::swap(M[pivot], M[current_row]);
        }

        /// 记录主元位置
        pivot_col[current_row] = col;
        pivot_row[col] = current_row;

        /// 将主元归一化为 1
        int64_t lc_inv = bk_mod_inverse(M[current_row][col], p);
        for (int j = 0; j < n; ++j) {
            M[current_row][j] = (M[current_row][j] * lc_inv) % p;
        }

        /// 消去其他行中该列的元素
        for (int row = 0; row < n; ++row) {
            if (row == current_row) continue;
            int64_t factor = M[row][col];
            if (factor == 0) continue;
            for (int j = 0; j < n; ++j) {
                M[row][j] = (M[row][j] - factor * M[current_row][j] % p + p) % p;
            }
        }

        current_row++;
    }

    /// 识别自由变量（无主元的列）并构造零空间基向量
    std::vector<std::vector<int64_t>> basis;

    for (int col = 0; col < n; ++col) {
        if (pivot_row[col] != -1) {
            /// 该列有主元，不是自由变量
            continue;
        }

        /// 该列为自由变量，构造对应的零空间基向量
        std::vector<int64_t> vec(n, 0);
        vec[col] = 1;  // 自由变量设为 1

        /// 回代：对每个有主元的行，计算对应分量
        for (int row = 0; row < n; ++row) {
            int pc = pivot_col[row];
            if (pc == -1) continue;  // 零行
            /// M[row] 已化为行简化阶梯形，M[row][pc] = 1
            /// 方程：x[pc] + sum_{free cols j} M[row][j] * x[j] = 0
            /// 所以 x[pc] = -M[row][col] mod p
            vec[pc] = (p - M[row][col]) % p;
        }

        basis.push_back(std::move(vec));
    }

    return basis;
}

/**
 * @brief 将模 p 系数向量转换为 Polynomial<ModInt>。
 *
 * @param[in] coeffs 系数向量（升幂排列）
 * @param[in] p      素数模数
 * @param[in] var    变量名
 * @return 对应的 Polynomial<ModInt>
 * @internal
 */
static Polynomial<ModInt> bk_to_poly_modint(
    const std::vector<int64_t>& coeffs,
    int64_t p,
    const std::string& var) {

    std::vector<ModInt> mod_coeffs;
    mod_coeffs.reserve(coeffs.size());
    for (int64_t c : coeffs) {
        mod_coeffs.emplace_back(c, p);
    }
    /// 直接构造，不调用 trim()（已保证无高次零系数）
    Polynomial<ModInt> result(var);
    result.coeffs = std::move(mod_coeffs);
    return result;
}

/**
 * @brief 使用零空间基向量将多项式分裂为不可约因子。
 *
 * 对每个非平凡零空间基向量 v，将其解释为多项式 v(x)，
 * 然后对当前因子列表中每个次数 > 1 的因子 g，
 * 计算 gcd(g, v(x) - c) 对 c = 0, 1, ..., p-1，
 * 若得到非平凡因子则分裂 g。
 *
 * @param[in] f_coeffs       首一多项式系数向量（升幂排列）
 * @param[in] null_basis     零空间基向量集合
 * @param[in] null_dim       零空间维度（= 期望的不可约因子数）
 * @param[in] p              素数模数
 * @return 不可约因子的系数向量列表
 *
 * @see Berlekamp, E.R. "Factoring polynomials over finite fields."
 *      Bell System Technical Journal, 46(8), 1967.
 * @internal
 */
static std::vector<std::vector<int64_t>> bk_split_factors(
    const std::vector<int64_t>& f_coeffs,
    const std::vector<std::vector<int64_t>>& null_basis,
    int null_dim,
    int64_t p) {

    /// 初始因子列表：仅包含 f 本身
    std::vector<std::vector<int64_t>> factor_list;
    factor_list.push_back(f_coeffs);

    /// 遍历每个零空间基向量
    for (const auto& basis_vec : null_basis) {
        /// 跳过平凡基向量 [1, 0, 0, ..., 0]
        bool is_trivial = true;
        if (!basis_vec.empty() && basis_vec[0] == 1) {
            for (size_t i = 1; i < basis_vec.size(); ++i) {
                if (basis_vec[i] != 0) {
                    is_trivial = false;
                    break;
                }
            }
        } else {
            is_trivial = false;
        }
        if (is_trivial) continue;

        /// 已找到足够因子则提前终止
        if (static_cast<int>(factor_list.size()) >= null_dim) break;

        /// 将基向量解释为多项式 v(x) 的系数（去除高次零系数）
        std::vector<int64_t> v_poly = basis_vec;
        while (!v_poly.empty() && v_poly.back() == 0) {
            v_poly.pop_back();
        }

        /// 对当前因子列表中每个次数 > 1 的因子尝试分裂
        for (size_t fi = 0; fi < factor_list.size(); ++fi) {
            /// 已找到足够因子则提前终止
            if (static_cast<int>(factor_list.size()) >= null_dim) break;

            /// 只对次数 > 1 的因子尝试分裂
            if (factor_list[fi].size() <= 2) continue;

            /// 对 c = 0, 1, ..., p-1 计算 gcd(g, v(x) - c)
            for (int64_t c = 0; c < p; ++c) {
                if (static_cast<int>(factor_list.size()) >= null_dim) break;

                /// 构造 v(x) - c：从常数项减去 c
                std::vector<int64_t> v_minus_c = v_poly;
                if (v_minus_c.empty()) {
                    /// v(x) = 0，则 v(x) - c = -c
                    v_minus_c.push_back((p - c) % p);
                } else {
                    v_minus_c[0] = (v_minus_c[0] - c % p + p) % p;
                }
                /// 去除高次零系数
                while (!v_minus_c.empty() && v_minus_c.back() == 0) {
                    v_minus_c.pop_back();
                }

                /// 若 v(x) - c 为零多项式则跳过
                if (v_minus_c.empty()) continue;

                /// 计算 h = gcd(g, v(x) - c)
                std::vector<int64_t> h = bk_gcd_mod(factor_list[fi], v_minus_c, p);

                /// 检查是否为非平凡因子：0 < deg(h) < deg(g)
                int deg_h = static_cast<int>(h.size()) - 1;
                int deg_g = static_cast<int>(factor_list[fi].size()) - 1;

                if (deg_h > 0 && deg_h < deg_g) {
                    /// 计算商 g / h
                    std::vector<int64_t> quotient, remainder;
                    bk_div_mod(factor_list[fi], h, p, quotient, remainder);

                    /// 替换 g 为 h，并将商加入因子列表
                    factor_list[fi] = h;
                    factor_list.push_back(quotient);

                    /// 当前因子 g 已被分裂，跳出 c 循环
                    /// （新的 factor_list[fi] = h 可能还能继续分裂，
                    ///   但由外层 fi 循环在下一轮处理）
                    break;
                }
            }
        }
    }

    return factor_list;
}

// ============================================================
/// 公共 API 实现
// ============================================================

/**
 * @brief Berlekamp 模素数分解。
 *
 * 在有限域 F_p 上对有理系数多项式执行 Berlekamp 算法，
 * 返回模 p 下的不可约因子列表。
 *
 * 当前为部分实现：完成素数选取与模约化，矩阵构造与因子分裂
 * 将在后续任务（3.2–3.4）中补充。
 *
 * @param[in] poly  有理系数多项式
 * @param[in] prime 指定的素数 p（若为 0 则自动选取）
 * @return Berlekamp 分解结果
 */
BerlekampResult berlekamp_factor(
    const Polynomial<Rational>& poly,
    int64_t prime) {

    BerlekampResult result;

    /// 零多项式或常数多项式：无因子
    if (poly.is_zero() || poly.degree() <= 0) {
        return result;
    }

    /// 素数选取：若调用方指定了 prime > 0 则使用之，否则自动选取
    int64_t p = prime;
    if (p <= 0) {
        p = bk_select_prime(poly);
        if (p < 0) {
            return result;
        }
    }

    result.prime = p;

    /// 将多项式约化到 F_p
    std::vector<int64_t> f_coeffs = bk_reduce_to_mod_coeffs(poly, p);

    /// 首一化
    if (!f_coeffs.empty()) {
        int64_t lc = f_coeffs.back();
        if (lc != 1) {
            int64_t lc_inv = bk_mod_inverse(lc, p);
            for (auto& c : f_coeffs) {
                c = (c * lc_inv) % p;
            }
        }
    }

    /// 线性多项式：本身不可约
    if (f_coeffs.size() <= 2) {
        result.factors.push_back(bk_to_poly_modint(f_coeffs, p, poly.variable_name));
        return result;
    }

    /// 构造 Berlekamp 矩阵 Q
    auto Q = bk_build_q_matrix(f_coeffs, p);

    /// 计算 (Q - I)^T 的零空间
    auto null_basis = bk_null_space(Q, p);
    int null_dim = static_cast<int>(null_basis.size());

    /// 存储零空间信息
    result.null_space_dim = null_dim;
    result.null_space_basis = null_basis;

    /// 零空间维度 = 1：多项式在 F_p 上不可约
    if (null_dim <= 1) {
        result.factors.push_back(bk_to_poly_modint(f_coeffs, p, poly.variable_name));
        return result;
    }

    /// 零空间维度 > 1：多项式可分解，使用基向量分裂
    auto split = bk_split_factors(f_coeffs, null_basis, null_dim, p);
    for (const auto& factor_coeffs : split) {
        result.factors.push_back(bk_to_poly_modint(factor_coeffs, p, poly.variable_name));
    }
    return result;
}

} // namespace lamina
