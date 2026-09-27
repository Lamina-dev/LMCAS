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
#include <lmmp.h>
#include <numth.h>
#include "result.hpp"

namespace LMCAS {

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
 * @brief 用扩展欧几里得算法求 gcd(a,b) 及 Bezout 系数。
 * @param a 第一个有符号机器整数
 * @param b 第二个有符号机器整数
 * @param s 输出 Bezout 系数 s，满足 a*s + b*t = gcd(a,b)
 * @param t 输出 Bezout 系数 t
 * @return 非负的 gcd(a, b)
 * @throws std::overflow_error 任一输入为 INT64_MIN，其绝对值无法由返回类型表示
 *
 * 输出参数仅在计算成功后写入。
 */
inline int64_t extended_gcd(int64_t a, int64_t b, int64_t& s, int64_t& t) {
    const uint64_t magnitude_a = modular_detail::abs_to_u64(a);
    const uint64_t magnitude_b = modular_detail::abs_to_u64(b);
    const auto maximum =
        static_cast<uint64_t>(std::numeric_limits<int64_t>::max());
    if (magnitude_a > maximum || magnitude_b > maximum) {
        throw std::overflow_error(
            "modular extended_gcd does not support INT64_MIN");
    }

    /** @brief 余数及交替 Bezout 系数以 max(|a|,|b|) 为界，包括余数归零的最后一步。 */
    int64_t old_r = static_cast<int64_t>(magnitude_a);
    int64_t r = static_cast<int64_t>(magnitude_b);
    int64_t old_s = 1;
    int64_t next_s = 0;
    int64_t old_t = 0;
    int64_t next_t = 1;

    while (r != 0) {
        const int64_t quotient = old_r / r;
        int64_t temporary = r;
        r = old_r - quotient * r;
        old_r = temporary;

        temporary = next_s;
        next_s = old_s - quotient * next_s;
        old_s = temporary;

        temporary = next_t;
        next_t = old_t - quotient * next_t;
        old_t = temporary;
    }

    const int64_t result_s = a < 0 ? -old_s : old_s;
    const int64_t result_t = b < 0 ? -old_t : old_t;
    s = result_s;
    t = result_t;
    return old_r;
}

/**
 * @brief 计算两个无符号整数的最大公约数（调用 LMMP）
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

    void require_same_modulus(const ModInt& other) const {
        if (mod_ != other.mod_) {
            throw std::domain_error("modular operands must have the same modulus");
        }
    }

public:
    /** @brief 默认构造，值为 0，模为 2 */
    ModInt() : val_(0), mod_(2) {}

    /**
     * @brief 构造模整数
     * @param v 整数值（自动取模归约到 [0, p)）
     * @param p 模数
     */
    ModInt(int64_t v, int64_t p) : mod_(p) {
        if (p <= 0) {
            throw std::domain_error("modulus must be positive");
        }
        val_ = v % p;
        if (val_ < 0) val_ += p;
    }

    /** @brief 获取当前值 */
    int64_t value() const { return val_; }

    /** @brief 获取模数 */
    int64_t modulus() const { return mod_; }

    /** @brief 模加法 */
    ModInt operator+(const ModInt& other) const {
        require_same_modulus(other);
        const int64_t distance = mod_ - other.val_;
        return ModInt(val_ >= distance ? val_ - distance : val_ + other.val_, mod_);
    }

    /** @brief 模减法 */
    ModInt operator-(const ModInt& other) const {
        require_same_modulus(other);
        return ModInt(val_ >= other.val_ ? val_ - other.val_
                                        : mod_ - (other.val_ - val_), mod_);
    }

    /** @brief 模乘法 */
    ModInt operator*(const ModInt& other) const {
        require_same_modulus(other);

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
        uint64_t power = 0;
        if (exp < 0) {
            base = base.inverse();
            power = static_cast<uint64_t>(-(exp + 1)) + 1;
        } else {
            power = static_cast<uint64_t>(exp);
        }
        ModInt result(1, base.mod_);
        while (power != 0) {
            if ((power & 1U) != 0) result = result * base;
            power >>= 1U;
            if (power != 0) base = base * base;
        }
        return result;
    }
};


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
    if (g != 1) {
        return CrtResult::failure(CasErrc::InvalidArgument,
                                  "CRT moduli must be coprime", operation);
    }

    const auto factor = ((ModInt(r2, m2) - ModInt(r1, m2)) * ModInt(s, m2)).value();
    /** @brief factor < m2 且 m1*m2 <= INT64_MAX，保证 m1*factor 在机器整数范围内。 */
    const auto merged = ModInt(r1, product) + ModInt(m1 * factor, product);
    return CrtResult::success({merged.value(), product});
}

inline CrtResult crt_checked(int64_t r1, int64_t m1, int64_t r2, int64_t m2) {
    ComputationContext context;
    return crt_checked(r1, m1, r2, m2, context);
}
/**
 * @brief 中国剩余定理（两模数），返回规范余数与模数乘积。
 * @throw std::domain_error 模数非正或不互素。
 * @throw std::overflow_error 合并模数超出机器内核范围。
 */
inline std::pair<int64_t, int64_t> crt(
    int64_t r1, int64_t m1, int64_t r2, int64_t m2) {
    auto merged = crt_checked(r1, m1, r2, m2);
    if (merged) return merged.value();
    if (merged.error().code == CasErrc::ResourceLimit) {
        throw std::overflow_error(merged.error().message);
    }
    throw std::domain_error(merged.error().message);
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

    int64_t modulus = primes[0];
    int64_t combined = ModInt(residues[0], modulus).value();

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

    int64_t modulus = primes[0];
    int64_t combined = ModInt(residues[0], modulus).value();
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
 * @brief Checked rational reconstruction. Returns Inconclusive when no
 * supported coprime numerator/denominator can be reconstructed.
 */
inline RationalReconstructionResult rational_reconstruction_checked(
    int64_t x,
    int64_t m,
    ComputationContext& context) {
    constexpr const char* operation = "rational_reconstruction";
    auto step = context.consume_steps(1, operation);
    if (!step) {
        return RationalReconstructionResult::failure(step.error());
    }
    if (m <= 0) {
        return RationalReconstructionResult::failure(
            CasErrc::InvalidArgument,
            "rational reconstruction modulus must be positive",
            operation);
    }

    x %= m;
    if (x < 0) {
        x += m;
    }
    int64_t bound = static_cast<int64_t>(
        std::floor(std::sqrt(static_cast<double>(m) / 2.0)));
    if (bound == 0) {
        bound = 1;
    }
    int64_t r0 = m;
    int64_t r1 = x;
    int64_t s0 = 0;
    int64_t s1 = 1;
    while (r1 > bound) {
        step = context.consume_steps(1, operation);
        if (!step) {
            return RationalReconstructionResult::failure(step.error());
        }
        const int64_t quotient = r0 / r1;
        const int64_t next_r = r0 - quotient * r1;
        const int64_t next_s = s0 - quotient * s1;
        r0 = r1;
        r1 = next_r;
        s0 = s1;
        s1 = next_s;
    }
    int64_t numerator = r1;
    int64_t denominator = s1;
    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }
    int64_t unused_s = 0;
    int64_t unused_t = 0;
    const bool outside_bound =
        std::abs(numerator) >= bound || denominator == 0 ||
        denominator >= bound;
    const bool not_coprime = !outside_bound &&
        extended_gcd(std::abs(numerator), denominator,
                     unused_s, unused_t) != 1;
    if (outside_bound || not_coprime) {
        return RationalReconstructionResult::failure(
            CasErrc::Inconclusive,
            "no rational reconstruction exists in the supported bound",
            operation);
    }
    return RationalReconstructionResult::success(
        {numerator, denominator});
}

inline RationalReconstructionResult rational_reconstruction_checked(int64_t x, int64_t m) {
    ComputationContext context;
    return rational_reconstruction_checked(x, m, context);
}

/** @brief 预定义的大素数表，用于模运算多素数方案 */
constexpr int64_t MODULAR_PRIMES[] = {
    1000000007, 1000000009, 1000000021, 1000000033,
    1000000087, 1000000093, 1000000097, 1000000103,
    1000000123, 1000000181, 1000000207, 1000000223,
    1000000231, 1000000271, 1000000289, 1000000297
};
/** @brief 预定义素数表的长度 */
constexpr size_t NUM_MODULAR_PRIMES = sizeof(MODULAR_PRIMES) / sizeof(MODULAR_PRIMES[0]);

/**
 * @brief 判断素数 p 是否为"好素数"（不整除任何首项系数）
 * @param p 待检测素数
 * @param leading_coeffs 首项系数列表
 * @return 若 p 不整除任何系数则返回 true
 */
inline bool is_good_prime(int64_t p, const std::vector<BigInt>& leading_coeffs) {
    if (p < 2) {
        throw std::domain_error("prime modulus must be at least two");
    }
    const BigInt modulus(p);
    return std::none_of(leading_coeffs.begin(), leading_coeffs.end(),
        [&modulus](const BigInt& coefficient) {
            return (coefficient % modulus).is_zero();
        });
}

/**
 * @brief 生成指定数量的好素数
 * @param count 需要的素数个数
 * @param leading_coeffs 首项系数列表（用于过滤）
 * @param start_from 搜索起始值
 * @return 好素数列表
 */
inline std::vector<int64_t> generate_good_primes(
    int count, const std::vector<BigInt>& leading_coeffs,
    uint64_t start_from = UINT64_C(1000000000)) {
    if (count < 0) {
        throw std::invalid_argument("prime count must be nonnegative");
    }
    const auto requested = static_cast<std::size_t>(count);
    const auto maximum = static_cast<uint64_t>(std::numeric_limits<int64_t>::max());
    std::vector<int64_t> result;
    result.reserve(requested);
    uint64_t candidate = start_from;
    while (result.size() < requested) {
        if (candidate >= maximum) {
            throw std::overflow_error("prime search exceeds the signed residue range");
        }
        const uint64_t next = lmmp_next_prime_ulong_(candidate);
        if (next <= candidate || next > maximum) {
            throw std::overflow_error("prime search exhausted the signed residue range");
        }
        candidate = next;
        const auto prime = static_cast<int64_t>(candidate);
        if (is_good_prime(prime, leading_coeffs)) {
            result.push_back(prime);
        }
    }
    return result;
}

}
