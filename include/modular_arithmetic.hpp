/**
 * @file modular_arithmetic.hpp
 * @brief 模运算工具：ModInt 类、扩展欧几里得、CRT、有理重构。
 */
#pragma once

#include <cstdint>
#include <vector>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <string>
#include "computation_context.hpp"
#include "lammp/lmmp.h"
#include "lammp/numth.h"
#include "result.hpp"

namespace lamina {

using CrtResult = Result<std::pair<int64_t, int64_t>>;
using RationalReconstructionResult = Result<std::pair<int64_t, int64_t>>;

namespace modular_detail {

inline uint64_t abs_to_u64(int64_t value) {
    if (value >= 0) return static_cast<uint64_t>(value);
    return static_cast<uint64_t>(-(value + 1)) + 1;
}

inline bool checked_positive_product(int64_t a, int64_t b, int64_t& out) {
    if (a <= 0 || b <= 0) return false;
    if (static_cast<uint64_t>(a) >
        static_cast<uint64_t>(std::numeric_limits<int64_t>::max()) / static_cast<uint64_t>(b)) {
        return false;
    }
    out = a * b;
    return true;
}

} // namespace modular_detail

/**
 * @brief 扩展欧几里得算法，求 gcd(a, b) 及 Bezout 系数
 * @param a 第一个整数
 * @param b 第二个整数
 * @param s 输出 Bezout 系数 s，满足 a*s + b*t = gcd(a,b)
 * @param t 输出 Bezout 系数 t
 * @return gcd(a, b)
 */
inline int64_t extended_gcd(int64_t a, int64_t b, int64_t& s, int64_t& t) {
    int64_t old_r = a, r = b;
    int64_t old_s = 1, ss = 0;
    int64_t old_t = 0, tt = 1;

    while (r != 0) {
        int64_t q = old_r / r;
        int64_t tmp = r;
        r = old_r - q * r;
        old_r = tmp;

        tmp = ss;
        ss = old_s - q * ss;
        old_s = tmp;

        tmp = tt;
        tt = old_t - q * tt;
        old_t = tmp;
    }

    s = old_s;
    t = old_t;
    return old_r;
}

/**
 * @brief 计算两个无符号整数的最大公约数（调用 LAMMP）
 * @param a 第一个无符号整数
 * @param b 第二个无符号整数
 * @return gcd(a, b)
 */
inline uint64_t lammp_gcd(uint64_t a, uint64_t b) {
    if (a == 0) return b;
    if (b == 0) return a;
    return lmmp_gcd_11_(a, b);
}

/**
 * @brief 判断无符号整数是否为素数
 * @param n 待判断的整数
 * @return 是素数返回 true
 */
inline bool lammp_is_prime(uint64_t n) {
    return lmmp_is_prime_ulong_(n);
}

/**
 * @brief 求大于 n 的下一个素数
 * @param n 起始值
 * @return 大于 n 的最小素数
 */
inline uint64_t lammp_next_prime(uint64_t n) {
    return lmmp_next_prime_ulong_(n);
}

/** @brief 模整数类，封装模 p 下的算术运算 */
class ModInt {
    int64_t val_;
    int64_t mod_;

public:
    /** @brief 默认构造，值为 0，模为 2 */
    ModInt() : val_(0), mod_(2) {}

    /**
     * @brief 构造模整数
     * @param v 整数值（自动取模归约到 [0, p)）
     * @param p 模数
     */
    ModInt(int64_t v, int64_t p) : mod_(p) {
        val_ = v % p;
        if (val_ < 0) val_ += p;
    }

    /** @brief 获取当前值 */
    int64_t value() const { return val_; }

    /** @brief 获取模数 */
    int64_t modulus() const { return mod_; }

    /** @brief 模加法 */
    ModInt operator+(const ModInt& other) const {
        return ModInt(val_ + other.val_, mod_);
    }

    /** @brief 模减法 */
    ModInt operator-(const ModInt& other) const {
        return ModInt(val_ - other.val_ + mod_, mod_);
    }

    /** @brief 模乘法 */
    ModInt operator*(const ModInt& other) const {

        uint64_t q;
        uint64_t result = lmmp_mulmod_ulong_(
            static_cast<uint64_t>(val_),
            static_cast<uint64_t>(other.val_),
            static_cast<uint64_t>(mod_), &q);
        return ModInt(static_cast<int64_t>(result), mod_);
    }

    /** @brief 模除法（乘以逆元） */
    ModInt operator/(const ModInt& other) const {
        return *this * other.inverse();
    }

    /** @brief 取负 */
    ModInt operator-() const {
        return ModInt(val_ == 0 ? 0 : mod_ - val_, mod_);
    }

    /** @brief 判等 */
    bool operator==(const ModInt& other) const {
        return val_ == other.val_ && mod_ == other.mod_;
    }

    /** @brief 判不等 */
    bool operator!=(const ModInt& other) const {
        return !(*this == other);
    }

    /** @brief 小于比较（先比模数，再比值） */
    bool operator<(const ModInt& other) const {
        if (mod_ != other.mod_) return mod_ < other.mod_;
        return val_ < other.val_;
    }

    /** @brief 大于比较 */
    bool operator>(const ModInt& other) const {
        return other < *this;
    }

    /** @brief 判断是否为零 */
    bool is_zero() const { return val_ == 0; }

    /**
     * @brief 求模逆元
     * @return 当前值在模 mod_ 下的乘法逆元
     * @throw std::domain_error 若值为 0 或不可逆
     */
    ModInt inverse() const {
        if (val_ == 0) {
            throw std::domain_error("ModInt::inverse(): zero is not invertible");
        }
        int64_t s, t;
        int64_t g = extended_gcd(val_, mod_, s, t);
        if (g != 1) {
            throw std::domain_error(
                "ModInt::inverse(): element is not invertible (gcd != 1)");
        }
        return ModInt(s, mod_);
    }

    /**
     * @brief 模幂运算
     * @param base 底数
     * @param exp 指数（可为负，负指数先求逆）
     * @return base^exp mod p
     */
    static ModInt pow(ModInt base, int64_t exp) {
        if (exp < 0) {
            base = base.inverse();
            exp = -exp;
        }
        if (exp == 0) return ModInt(1, base.modulus());

        uint64_t result = lmmp_powmod_ulong_(
            static_cast<uint64_t>(base.val_),
            static_cast<uint64_t>(exp),
            static_cast<uint64_t>(base.mod_));
        return ModInt(static_cast<int64_t>(result), base.mod_);
    }
};

/**
 * @brief 中国剩余定理（两模数）
 * @param r1 第一个余数
 * @param m1 第一个模数
 * @param r2 第二个余数
 * @param m2 第二个模数
 * @return pair(合并余数, 合并模数 m1*m2)
 * @throw std::domain_error 若模数不互素
 */
inline std::pair<int64_t, int64_t> crt(int64_t r1, int64_t m1,
                                        int64_t r2, int64_t m2) {
    int64_t s, t;
    int64_t g = extended_gcd(m1, m2, s, t);
    if (g != 1) {
        throw std::domain_error("crt(): moduli must be coprime");
    }

    int64_t diff = ((r2 - r1) % m2 + m2) % m2;

    int64_t s_mod = ((s % m2) + m2) % m2;

    uint64_t q_unused;
    uint64_t factor = lmmp_mulmod_ulong_(
        static_cast<uint64_t>(s_mod),
        static_cast<uint64_t>(diff),
        static_cast<uint64_t>(m2), &q_unused);

#if defined(__GNUC__) || defined(__clang__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wpedantic"
    __int128 x = static_cast<__int128>(r1) + static_cast<__int128>(m1) * factor;
    __int128 bigM = static_cast<__int128>(m1) * m2;
    x = ((x % bigM) + bigM) % bigM;
    #pragma GCC diagnostic pop
    return {static_cast<int64_t>(x), static_cast<int64_t>(bigM)};
#else
    int64_t M = m1 * m2;
    int64_t x = r1 + m1 * static_cast<int64_t>(factor);
    x = ((x % M) + M) % M;
    return {x, M};
#endif
}

/**
 * @brief Checked Chinese remainder theorem for two positive coprime moduli.
 */
inline CrtResult crt_checked(int64_t r1, int64_t m1,
                             int64_t r2, int64_t m2,
                             ComputationContext& context) {
    constexpr const char* operation = "crt";
    auto step = context.consume_steps(1, operation);
    if (!step) return CrtResult::failure(step.error());
    if (m1 <= 0 || m2 <= 0) {
        return CrtResult::failure(CasErrc::InvalidArgument,
                                  "CRT moduli must be positive", operation);
    }

    int64_t product = 0;
    if (!modular_detail::checked_positive_product(m1, m2, product)) {
        return CrtResult::failure(CasErrc::ResourceLimit,
                                  "CRT modulus product exceeds int64 range", operation);
    }

    int64_t s = 0;
    int64_t t = 0;
    int64_t g = extended_gcd(m1, m2, s, t);
    if (g < 0) g = -g;
    if (g != 1) {
        return CrtResult::failure(CasErrc::InvalidArgument,
                                  "CRT moduli must be coprime", operation);
    }

#if defined(__GNUC__) || defined(__clang__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wpedantic"
    const __int128 big_m2 = static_cast<__int128>(m2);
    const __int128 diff = ((static_cast<__int128>(r2) - r1) % big_m2 + big_m2) % big_m2;
    const __int128 s_mod = ((static_cast<__int128>(s) % big_m2) + big_m2) % big_m2;
    __int128 factor = (s_mod * diff) % big_m2;
    __int128 x = static_cast<__int128>(r1) + static_cast<__int128>(m1) * factor;
    const __int128 big_product = static_cast<__int128>(product);
    x = ((x % big_product) + big_product) % big_product;
    #pragma GCC diagnostic pop
    return CrtResult::success({static_cast<int64_t>(x), product});
#else
    try {
        return CrtResult::success(crt(r1, m1, r2, m2));
    } catch (const std::bad_alloc&) {
        return CrtResult::failure(CasErrc::ResourceLimit,
                                  "CRT allocation failed", operation);
    } catch (const std::exception& ex) {
        return CrtResult::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
#endif
}

inline CrtResult crt_checked(int64_t r1, int64_t m1, int64_t r2, int64_t m2) {
    ComputationContext context;
    return crt_checked(r1, m1, r2, m2, context);
}

/**
 * @brief 多模数中国剩余定理
 * @param residues 余数列表
 * @param primes 对应的模数列表（两两互素）
 * @return pair(合并余数, 合并模数)
 * @throw std::invalid_argument 若输入为空或大小不匹配
 */
inline std::pair<int64_t, int64_t> multi_crt(const std::vector<int64_t>& residues,
                                              const std::vector<int64_t>& primes) {
    if (residues.empty() || residues.size() != primes.size()) {
        throw std::invalid_argument("multi_crt(): residues and primes must be non-empty and same size");
    }

    int64_t combined = residues[0];
    int64_t modulus = primes[0];

    for (size_t i = 1; i < residues.size(); ++i) {
        auto [x, m] = crt(combined, modulus, residues[i], primes[i]);
        combined = x;
        modulus = m;
    }

    return {combined, modulus};
}

/**
 * @brief Checked CRT for a non-empty list of pairwise coprime positive moduli.
 */
inline CrtResult multi_crt_checked(const std::vector<int64_t>& residues,
                                   const std::vector<int64_t>& primes,
                                   ComputationContext& context) {
    constexpr const char* operation = "multi_crt";
    auto step = context.consume_steps(residues.size() + 1, operation);
    if (!step) return CrtResult::failure(step.error());
    if (residues.empty() || residues.size() != primes.size()) {
        return CrtResult::failure(
            CasErrc::InvalidArgument,
            "multi_crt residues and moduli must be non-empty and the same size",
            operation);
    }
    for (int64_t modulus : primes) {
        if (modulus <= 0) {
            return CrtResult::failure(CasErrc::InvalidArgument,
                                      "multi_crt moduli must be positive", operation);
        }
    }

    int64_t combined = residues[0];
    int64_t modulus = primes[0];
    for (size_t i = 1; i < residues.size(); ++i) {
        auto merged = crt_checked(combined, modulus, residues[i], primes[i], context);
        if (!merged) return merged;
        combined = merged.value().first;
        modulus = merged.value().second;
    }

    return CrtResult::success({combined, modulus});
}

inline CrtResult multi_crt_checked(const std::vector<int64_t>& residues,
                                   const std::vector<int64_t>& primes) {
    ComputationContext context;
    return multi_crt_checked(residues, primes, context);
}

/**
 * @brief 有理数重构：从模像 x mod m 恢复有理数 a/b
 * @param x 模像值
 * @param m 模数
 * @return pair(分子 a, 分母 b)；失败时返回 (0, 0)
 */
inline std::pair<int64_t, int64_t> rational_reconstruction(int64_t x, int64_t m) {
    if (m <= 0) {
        return {0, 0};
    }

    x = ((x % m) + m) % m;

    int64_t bound = static_cast<int64_t>(std::floor(std::sqrt(static_cast<double>(m) / 2.0)));
    if (bound == 0) bound = 1;

    int64_t r0 = m, r1 = x;
    int64_t s0 = 0, s1 = 1;

    while (r1 > bound) {
        int64_t q = r0 / r1;
        int64_t r_new = r0 - q * r1;
        int64_t s_new = s0 - q * s1;

        r0 = r1;
        r1 = r_new;
        s0 = s1;
        s1 = s_new;
    }

    int64_t a = r1;
    int64_t b = s1;

    if (b < 0) {
        a = -a;
        b = -b;
    }

    if (std::abs(a) >= bound || b == 0 || b >= bound) {

        return {0, 0};
    }

    int64_t s_unused, t_unused;
    int64_t g = extended_gcd(std::abs(a), b, s_unused, t_unused);
    if (g != 1) {
        return {0, 0};
    }

    return {a, b};
}

/**
 * @brief Checked rational reconstruction. Returns Inconclusive when no
 * supported coprime numerator/denominator can be reconstructed.
 */
inline RationalReconstructionResult rational_reconstruction_checked(
    int64_t x,
    int64_t m,
    ComputationContext& context) {
    constexpr const char* operation = "rational_reconstruction";
    auto step = context.consume_steps(1, operation);
    if (!step) return RationalReconstructionResult::failure(step.error());
    if (m <= 0) {
        return RationalReconstructionResult::failure(
            CasErrc::InvalidArgument, "rational reconstruction modulus must be positive", operation);
    }

    auto reconstructed = rational_reconstruction(x, m);
    if (reconstructed.second == 0) {
        return RationalReconstructionResult::failure(
            CasErrc::Inconclusive,
            "no rational reconstruction exists in the supported bound",
            operation);
    }
    return RationalReconstructionResult::success(reconstructed);
}

inline RationalReconstructionResult rational_reconstruction_checked(int64_t x, int64_t m) {
    ComputationContext context;
    return rational_reconstruction_checked(x, m, context);
}

/** @brief 预定义的大素数表，用于模运算多素数方案 */
constexpr int64_t MODULAR_PRIMES[] = {
    1000000007LL, 1000000009LL, 1000000021LL, 1000000033LL,
    1000000087LL, 1000000093LL, 1000000097LL, 1000000103LL,
    1000000123LL, 1000000181LL, 1000000207LL, 1000000223LL,
    1000000231LL, 1000000271LL, 1000000289LL, 1000000297LL
};
/** @brief 预定义素数表的长度 */
constexpr size_t NUM_MODULAR_PRIMES = sizeof(MODULAR_PRIMES) / sizeof(MODULAR_PRIMES[0]);

/**
 * @brief 判断素数 p 是否为"好素数"（不整除任何首项系数）
 * @param p 待检测素数
 * @param leading_coeffs 首项系数列表
 * @return 若 p 不整除任何系数则返回 true
 */
inline bool is_good_prime(int64_t p, const std::vector<int64_t>& leading_coeffs) {
    return std::none_of(leading_coeffs.begin(), leading_coeffs.end(),
        [p](int64_t c) { return (c % p) == 0; });
}

/**
 * @brief 生成指定数量的好素数
 * @param count 需要的素数个数
 * @param leading_coeffs 首项系数列表（用于过滤）
 * @param start_from 搜索起始值
 * @return 好素数列表
 */
inline std::vector<int64_t> generate_good_primes(int count,
                                                  const std::vector<int64_t>& leading_coeffs,
                                                  uint64_t start_from = 1000000000ULL) {
    std::vector<int64_t> result;
    result.reserve(count);
    uint64_t candidate = start_from;
    while ((int)result.size() < count) {
        candidate = lmmp_next_prime_ulong_(candidate);
        if (is_good_prime(static_cast<int64_t>(candidate), leading_coeffs)) {
            result.push_back(static_cast<int64_t>(candidate));
        }
    }
    return result;
}

}
