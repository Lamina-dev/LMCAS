#pragma once

#include <cstdint>
#include <vector>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "lammp/lmmp.h"
#include "lammp/numth.h"

namespace lamina {

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

inline uint64_t lammp_gcd(uint64_t a, uint64_t b) {
    if (a == 0) return b;
    if (b == 0) return a;
    return lmmp_gcd_11_(a, b);
}

inline bool lammp_is_prime(uint64_t n) {
    return lmmp_is_prime_ulong_(n);
}

inline uint64_t lammp_next_prime(uint64_t n) {
    return lmmp_next_prime_ulong_(n);
}

class ModInt {
    int64_t val_;
    int64_t mod_;

public:
    ModInt() : val_(0), mod_(2) {}

    ModInt(int64_t v, int64_t p) : mod_(p) {
        val_ = v % p;
        if (val_ < 0) val_ += p;
    }

    int64_t value() const { return val_; }
    int64_t modulus() const { return mod_; }

    ModInt operator+(const ModInt& other) const {
        return ModInt(val_ + other.val_, mod_);
    }

    ModInt operator-(const ModInt& other) const {
        return ModInt(val_ - other.val_ + mod_, mod_);
    }

    ModInt operator*(const ModInt& other) const {

        uint64_t q;
        uint64_t result = lmmp_mulmod_ulong_(
            static_cast<uint64_t>(val_),
            static_cast<uint64_t>(other.val_),
            static_cast<uint64_t>(mod_), &q);
        return ModInt(static_cast<int64_t>(result), mod_);
    }

    ModInt operator/(const ModInt& other) const {
        return *this * other.inverse();
    }

    ModInt operator-() const {
        return ModInt(val_ == 0 ? 0 : mod_ - val_, mod_);
    }

    bool operator==(const ModInt& other) const {
        return val_ == other.val_ && mod_ == other.mod_;
    }

    bool operator!=(const ModInt& other) const {
        return !(*this == other);
    }

    bool operator<(const ModInt& other) const {
        if (mod_ != other.mod_) return mod_ < other.mod_;
        return val_ < other.val_;
    }

    bool operator>(const ModInt& other) const {
        return other < *this;
    }

    bool is_zero() const { return val_ == 0; }

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

inline std::pair<int64_t, int64_t> crt(int64_t r1, int64_t m1,
                                        int64_t r2, int64_t m2) {
    int64_t s, t;
    int64_t g = extended_gcd(m1, m2, s, t);
    if (g != 1) {
        throw std::domain_error("crt(): moduli must be coprime");
    }

    int64_t M = m1 * m2;
    int64_t diff = ((r2 - r1) % m2 + m2) % m2;

    int64_t s_mod = ((s % m2) + m2) % m2;

    uint64_t q_unused;
    uint64_t factor = lmmp_mulmod_ulong_(
        static_cast<uint64_t>(s_mod),
        static_cast<uint64_t>(diff),
        static_cast<uint64_t>(m2), &q_unused);

#if defined(__GNUC__) || defined(__clang__)
    __int128 x = static_cast<__int128>(r1) + static_cast<__int128>(m1) * factor;
    __int128 bigM = static_cast<__int128>(m1) * m2;
    x = ((x % bigM) + bigM) % bigM;
    return {static_cast<int64_t>(x), static_cast<int64_t>(bigM)};
#else
    int64_t x = r1 + m1 * static_cast<int64_t>(factor);
    x = ((x % M) + M) % M;
    return {x, M};
#endif
}

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

constexpr int64_t MODULAR_PRIMES[] = {
    1000000007LL, 1000000009LL, 1000000021LL, 1000000033LL,
    1000000087LL, 1000000093LL, 1000000097LL, 1000000103LL,
    1000000123LL, 1000000181LL, 1000000207LL, 1000000223LL,
    1000000231LL, 1000000271LL, 1000000289LL, 1000000297LL
};
constexpr size_t NUM_MODULAR_PRIMES = sizeof(MODULAR_PRIMES) / sizeof(MODULAR_PRIMES[0]);

inline bool is_good_prime(int64_t p, const std::vector<int64_t>& leading_coeffs) {
    return std::none_of(leading_coeffs.begin(), leading_coeffs.end(),
        [p](int64_t c) { return (c % p) == 0; });
}

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
