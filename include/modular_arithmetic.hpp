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

// ---------------------------------------------------------------------------
// Extended GCD
// ---------------------------------------------------------------------------

/// Extended Euclidean algorithm: returns gcd(a, b) and sets s, t such that
/// a*s + b*t = gcd(a, b).
/// Uses LAMMP's lmmp_gcd_11_ for the GCD computation internally.
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
    return old_r; // gcd
}

/// Compute gcd of two positive integers using LAMMP's optimized implementation.
inline uint64_t lammp_gcd(uint64_t a, uint64_t b) {
    if (a == 0) return b;
    if (b == 0) return a;
    return lmmp_gcd_11_(a, b);
}

/// Check primality using LAMMP's Miller-Rabin implementation.
inline bool lammp_is_prime(uint64_t n) {
    return lmmp_is_prime_ulong_(n);
}

/// Get next prime after n using LAMMP.
inline uint64_t lammp_next_prime(uint64_t n) {
    return lmmp_next_prime_ulong_(n);
}

// ---------------------------------------------------------------------------
// ModInt — arithmetic in Z/pZ
// ---------------------------------------------------------------------------

/// Represents an element of Z/pZ with full arithmetic support.
class ModInt {
    int64_t val_;  // value in [0, p)
    int64_t mod_;  // the prime modulus

public:
    ModInt() : val_(0), mod_(2) {}

    /// Construct from value v and prime modulus p. Normalizes v into [0, p).
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
        // Use LAMMP's optimized modular multiplication
        uint64_t q;
        uint64_t result = lmmp_mulmod_ulong_(
            static_cast<uint64_t>(val_), 
            static_cast<uint64_t>(other.val_), 
            static_cast<uint64_t>(mod_), &q);
        return ModInt(static_cast<int64_t>(result), mod_);
    }

    /// Division via modular inverse. Throws if other is not invertible.
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

    /// Comparison for use in ordered containers (lexicographic on (mod, val)).
    bool operator<(const ModInt& other) const {
        if (mod_ != other.mod_) return mod_ < other.mod_;
        return val_ < other.val_;
    }

    bool operator>(const ModInt& other) const {
        return other < *this;
    }

    bool is_zero() const { return val_ == 0; }

    /// Modular inverse using extended Euclidean algorithm.
    /// Throws std::domain_error if the element is not invertible (gcd != 1).
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

    /// Modular exponentiation using LAMMP's optimized implementation.
    static ModInt pow(ModInt base, int64_t exp) {
        if (exp < 0) {
            base = base.inverse();
            exp = -exp;
        }
        if (exp == 0) return ModInt(1, base.modulus());
        // Use LAMMP's lmmp_powmod_ulong_ for the heavy lifting
        uint64_t result = lmmp_powmod_ulong_(
            static_cast<uint64_t>(base.val_),
            static_cast<uint64_t>(exp),
            static_cast<uint64_t>(base.mod_));
        return ModInt(static_cast<int64_t>(result), base.mod_);
    }
};

// ---------------------------------------------------------------------------
// Chinese Remainder Theorem
// ---------------------------------------------------------------------------

/// Two-prime CRT: given r1 mod m1 and r2 mod m2 (m1, m2 coprime),
/// find x such that x ≡ r1 (mod m1) and x ≡ r2 (mod m2).
/// Returns (x, m1*m2).
/// Uses LAMMP's lmmp_mulmod_ulong_ for overflow-safe modular multiplication.
inline std::pair<int64_t, int64_t> crt(int64_t r1, int64_t m1,
                                        int64_t r2, int64_t m2) {
    int64_t s, t;
    int64_t g = extended_gcd(m1, m2, s, t);
    if (g != 1) {
        throw std::domain_error("crt(): moduli must be coprime");
    }

    // x = r1 + m1 * s * (r2 - r1) mod (m1 * m2)
    // Use LAMMP's mulmod for the intermediate product to avoid overflow
    int64_t M = m1 * m2; // This may overflow for very large primes; use __int128 as fallback
    int64_t diff = ((r2 - r1) % m2 + m2) % m2;
    
    // Compute s mod m2 (ensure positive)
    int64_t s_mod = ((s % m2) + m2) % m2;
    
    // Use LAMMP mulmod: s_mod * diff mod m2
    uint64_t q_unused;
    uint64_t factor = lmmp_mulmod_ulong_(
        static_cast<uint64_t>(s_mod),
        static_cast<uint64_t>(diff),
        static_cast<uint64_t>(m2), &q_unused);
    
    // x = r1 + m1 * factor
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

/// Multi-prime CRT: combine multiple residues.
/// residues[i] is the value mod primes[i].
/// Returns (combined_value, product_of_all_primes).
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

// ---------------------------------------------------------------------------
// Rational Reconstruction
// ---------------------------------------------------------------------------

/// Given x mod m, find a/b in Q such that a/b ≡ x (mod m),
/// with |a| < bound and 0 < b < bound, where bound = floor(sqrt(m/2)).
/// Returns (numerator, denominator) or (0, 0) if reconstruction fails.
///
/// Uses the half-GCD / lattice reduction approach via the Euclidean algorithm.
inline std::pair<int64_t, int64_t> rational_reconstruction(int64_t x, int64_t m) {
    if (m <= 0) {
        return {0, 0};
    }

    // Normalize x into [0, m)
    x = ((x % m) + m) % m;

    // bound = floor(sqrt(m / 2))
    int64_t bound = static_cast<int64_t>(std::floor(std::sqrt(static_cast<double>(m) / 2.0)));
    if (bound == 0) bound = 1;

    // Run extended Euclidean algorithm on (m, x)
    // r0 = m, r1 = x, s0 = 0, s1 = 1
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

    // a = r1, b = s1
    int64_t a = r1;
    int64_t b = s1;

    // Ensure b > 0
    if (b < 0) {
        a = -a;
        b = -b;
    }

    // Check |a| < bound and 0 < b < bound
    if (std::abs(a) >= bound || b == 0 || b >= bound) {
        // Relaxed check: still valid if gcd(a, b) == 1
        // but bounds are not met — reconstruction failed
        return {0, 0};
    }

    // Verify gcd(|a|, b) == 1
    int64_t s_unused, t_unused;
    int64_t g = extended_gcd(std::abs(a), b, s_unused, t_unused);
    if (g != 1) {
        return {0, 0};
    }

    return {a, b};
}

// ---------------------------------------------------------------------------
// Prime selection utilities
// ---------------------------------------------------------------------------

/// Good primes for modular computation (close to 10^9, fit in int64_t).
/// Products of two such primes fit in int64_t; products of three require __int128.
/// Verified prime using LAMMP's lmmp_is_prime_ulong_.
constexpr int64_t MODULAR_PRIMES[] = {
    1000000007LL, 1000000009LL, 1000000021LL, 1000000033LL,
    1000000087LL, 1000000093LL, 1000000097LL, 1000000103LL,
    1000000123LL, 1000000181LL, 1000000207LL, 1000000223LL,
    1000000231LL, 1000000271LL, 1000000289LL, 1000000297LL
};
constexpr size_t NUM_MODULAR_PRIMES = sizeof(MODULAR_PRIMES) / sizeof(MODULAR_PRIMES[0]);

/// Check if a prime p is "good" for a given polynomial system.
/// A prime is good if none of the leading coefficients vanish mod p.
inline bool is_good_prime(int64_t p, const std::vector<int64_t>& leading_coeffs) {
    return std::none_of(leading_coeffs.begin(), leading_coeffs.end(),
        [p](int64_t c) { return (c % p) == 0; });
}

/// Generate a sequence of good primes for modular computation using LAMMP.
/// Starts from a given base and finds primes that don't divide any leading coefficient.
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

} // namespace lamina
