#include "bigint.hpp"
#include "internal/lmmc_lifecycle.hpp"

namespace LMCAS {
namespace {

BigInt binary_modular_power(const BigInt& base, const BigInt& exp, const BigInt& mod) {
    BigInt result = 1;
    BigInt residue = base % mod;
    if (residue.is_negative()) {
        residue += mod;
    }
    BigInt exponent = exp;
    while (!exponent.is_zero()) {
        if (exponent.is_odd()) {
            result = (result * residue) % mod;
        }
        residue = (residue * residue) % mod;
        exponent >>= 1;
    }
    return result;
}

}

BigInt BigInt::multinomial(unsigned int n, const std::vector<unsigned int>& r) {
        LMCAS::detail::ensure_lmmc_lifecycle();
        if (r.empty()) return BigInt(1);

        std::vector<uint> r_uints;
        r_uints.reserve(r.size());
        ulong sum = 0;
        for(auto val : r) {
            r_uints.push_back((uint)val);
            sum += val;
        }
        if (sum != n) throw std::invalid_argument("multinomial: sum of ranks must equal n");

        BigInt res;

        ulong n_calc = 0;
        mp_size_t needed = lmmp_multinomial_size_(r_uints.data(), (uint)r_uints.size(), &n_calc);

        res.realloc_to(needed);
        res.size_ = lmmp_multinomial_(res.data_, needed, (uint)sum, r_uints.data(), (uint)r_uints.size());
        res.sign_ = POSITIVE;
        res.normalize();
        return res;
    }

BigInt BigInt::gcd(const BigInt& a, const BigInt& b) {
        LMCAS::detail::ensure_lmmc_lifecycle();
        if (a.is_zero()) return b.abs();
        if (b.is_zero()) return a.abs();

        BigInt abs_a = a.abs();
        BigInt abs_b = b.abs();

        BigInt res;
        mp_size_t na = abs_a.size_;
        mp_size_t nb = abs_b.size_;

        mp_size_t min_n = (na < nb) ? na : nb;
        res.realloc_to(min_n);

        res.size_ = lmmp_gcd_lehmer_(res.data_, abs_a.data_, na, abs_b.data_, nb);

        res.sign_ = POSITIVE;
        res.normalize();
        return res;
    }

BigInt BigInt::lcm(const BigInt& a, const BigInt& b) {
        if (a.is_zero() || b.is_zero()) return BigInt(0);
        return (a.abs() / gcd(a, b)) * b.abs();
    }

BigInt BigInt::pow_mod(const BigInt& base, const BigInt& exp, const BigInt& mod) {
        LMCAS::detail::ensure_lmmc_lifecycle();
        if (exp.is_negative()) {
            throw std::domain_error(
                "Negative exponent in modular integer power");
        }
        if (mod.is_zero() || mod.is_negative()) {
            throw std::domain_error(
                "Modulus must be positive in modular integer power");
        }
        if (mod == BigInt(1)) return BigInt(0);

        if (base.size_ <= 1 && exp.size_ <= 1 && mod.size_ <= 1) {
            const ulong modulus = mod.data_[0];
            ulong residue = base.is_zero() ? 0 : base.data_[0] % modulus;
            const ulong exponent = exp.is_zero() ? 0 : exp.data_[0];
            if (base.is_negative() && residue != 0) {
                residue = modulus - residue;
            }
            if ((modulus & 1U) != 0) {
                return BigInt(
                    lmmp_powmod_ulong_odd_(residue, exponent, modulus));
            }
        }

        return binary_modular_power(base, exp, mod);
    }

Result<bool> BigInt::is_prime_checked() const {
        LMCAS::detail::ensure_lmmc_lifecycle();
        if (sign_ == NEGATIVE) return false;

        if (size_ <= 1) {
            return lmmp_is_prime_ulong_(size_ == 0 ? 0 : data_[0]);
        }
        if (is_even()) return false;

        static const uint64_t small_primes[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
        for (auto p : small_primes) {
            if ((*this % p).is_zero()) return false;
        }

        const BigInt& n = *this;
        const BigInt n_minus_one = n - 1;
        BigInt d = n_minus_one;
        BigInt two(2);
        int s = 0;
        while (d.is_even()) {
            d = d / two;
            s++;
        }

        static const uint64_t mr_bases[] = {2, 3, 5, 7, 11, 13, 17};
        for (auto base : mr_bases) {
            BigInt x = BigInt::pow_mod(BigInt(base), d, n);
            if (x == 1 || x == n_minus_one) continue;

            bool composite = true;
            for (int r = 1; r < s; ++r) {
                x = BigInt::pow_mod(x, two, n);
                if (x == n_minus_one) {
                    composite = false;
                    break;
                }
            }
            if (composite) return false;
        }
        return Result<bool>::failure(
            CasErrc::Inconclusive, "primality is unproved beyond the LMMP word range",
            "BigInt::is_prime_checked");
    }

bool BigInt::is_perfect_square() const {
        if (sign_ == NEGATIVE) return false;
        if (size_ == 0) return true;
        BigInt s = this->sqrt();
        return (s * s) == *this;
    }

}
