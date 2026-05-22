#pragma once

#ifndef _STATIC_ASSERT
#define _STATIC_ASSERT(x) static_assert(x, #x)
#endif

#include "lammp/lmmp.h"
#include "lammp/lmmpn.h"
#include "lammp/numth.h"
#include "lmmc/init.h"
#include <vector>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include "lmmc/config.h"

inline void bigint_free(void* p) {
    lmmp_free(p);
}

inline void* bigint_alloc(size_t size) {
    return lmmp_alloc(size);
}

inline mp_size_t lmmp_rlz(mp_srcptr p, mp_size_t n) {
    while (n > 0 && p[n-1] == 0) n--;
    return n;
}

inline void bigint_mul_(mp_ptr dst, mp_srcptr numa, mp_size_t na, mp_srcptr numb, mp_size_t nb) {
    ::lmmp_mul_(dst, numa, na, numb, nb);
}

inline void bigint_div_(mp_ptr q, mp_ptr r, mp_srcptr n, mp_size_t nn, mp_srcptr d, mp_size_t dn) {
    ::lmmp_div_(q, r, n, nn, d, dn);
}

class BigInt {
public:
    enum Sign {
        POSITIVE = 1,
        ZERO = 0,
        NEGATIVE = -1
    };

    mp_ptr _data = nullptr;
    mp_size_t _size = 0;
    mp_size_t _alloc = 0;
    int _sign = ZERO;
    bool negative = false;

    void realloc_to(mp_size_t new_alloc) {
        if (new_alloc <= _alloc) return;
        new_alloc = (new_alloc + 3) & ~3;

        mp_ptr new_data = (mp_ptr)bigint_alloc(new_alloc * sizeof(mp_limb_t));

        if (_size > 0 && _data) {
             std::memcpy(new_data, _data, _size * sizeof(mp_limb_t));
        }
        if (_data) bigint_free(_data);
        _data = new_data;
        _alloc = new_alloc;
    }

    void normalize() {
        while (_size > 0 && _data[_size - 1] == 0) {
            _size--;
        }
        if (_size == 0) {
            _sign = ZERO;
            negative = false;
        } else {
             if (_sign == ZERO) _sign = POSITIVE;
             negative = (_sign == NEGATIVE);
        }
    }

    void zero() {
        _size = 0;
        _sign = ZERO;
        negative = false;
    }

public:
    BigInt() : _data(nullptr), _size(0), _alloc(0), _sign(ZERO), negative(false) {}

    ~BigInt() {
        if (_data) bigint_free(_data);
    }

    BigInt(const BigInt& other) {
        if (other._size > 0) {
            realloc_to(other._size);
            std::memcpy(_data, other._data, other._size * sizeof(mp_limb_t));
            _size = other._size;
            _sign = other._sign;
            negative = other.negative;
        } else {
            zero();
        }
    }

    BigInt(BigInt&& other) noexcept
        : _data(other._data), _size(other._size), _alloc(other._alloc), _sign(other._sign), negative(other.negative) {
        other._data = nullptr;
        other._size = 0;
        other._alloc = 0;
        other._sign = ZERO;
        other.negative = false;
    }

    BigInt& operator=(const BigInt& other) {
        if (this != &other) {
            if (other._size > _alloc) {
                realloc_to(other._size);
            }
            if (other._size > 0)
                std::memcpy(_data, other._data, other._size * sizeof(mp_limb_t));
            _size = other._size;
            _sign = other._sign;
            negative = other.negative;
        }
        return *this;
    }

    BigInt& operator=(BigInt&& other) noexcept {
        if (this != &other) {
            if (_data) bigint_free(_data);
            _data = other._data;
            _size = other._size;
            _alloc = other._alloc;
            _sign = other._sign;
            negative = other.negative;

            other._data = nullptr;
            other.zero();
        }
        return *this;
    }

    BigInt(long long val) {
        if (val == 0) {
            zero();
            return;
        }
        realloc_to(1);
        if (val < 0) {
            _sign = NEGATIVE;
            negative = true;

            if (val == std::numeric_limits<long long>::min()) {
                 _data[0] = (uint64_t)(-(val + 1)) + 1;
            } else {
                _data[0] = -val;
            }
        } else {
            _sign = POSITIVE;
            negative = false;
            _data[0] = val;
        }
        _size = 1;
    }
    BigInt(int val) : BigInt((long long)val) {}

    BigInt(unsigned long long val) {
        if (val == 0) {
            zero();
            return;
        }
        realloc_to(1);
        _sign = POSITIVE;
        negative = false;
        _data[0] = val;
        _size = 1;
    }
    BigInt(unsigned int val) : BigInt((unsigned long long)val) {}
    BigInt(unsigned long val) : BigInt((unsigned long long)val) {}

    BigInt(const std::string& str) {
        if (str.empty()) { zero(); return; }
        size_t start = 0;
        int sign = POSITIVE;
        if (str[0] == '-') {
            sign = NEGATIVE;
            start = 1;
        } else if (str[0] == '+') {
            start = 1;
        }

        if (start == str.length()) { zero(); return; }

        size_t len = str.length() - start;

        mp_size_t needed = len / 19 + 2;
        realloc_to(needed);

        std::vector<mp_byte_t> digit_buf(len);
        for (size_t i = 0; i < len; ++i) {
            char c = str[start + i];
            uint8_t d = 0;
            if (c >= '0' && c <= '9') {
                d = c - '0';
            }

            digit_buf[len - 1 - i] = d;
        }

        _size = lmmp_from_str_(_data, digit_buf.data(), len, 10);

        if (_size == 0) {
            zero();
        } else {
            _sign = sign;
            negative = (_sign == NEGATIVE);
        }
        normalize();
    }

    std::string ToString() const {
        if (_size == 0) return "0";

        size_t len_needed = _size * 20 + 5;
        std::vector<mp_byte_t> buf(len_needed);

        mp_size_t str_len = lmmp_to_str_((mp_byte_t*)buf.data(), _data, _size, 10);

        if (str_len == 0) return "0";

        std::string res;
        res.reserve(str_len + 2);
        if (_sign == NEGATIVE) res += '-';

        for(mp_size_t i = str_len; i > 0; --i) {
            res += (char)(buf[i - 1] + '0');
        }
        return res;
    }
    std::string to_string() const { return ToString(); }

    int to_int() const {
        if (_size == 0) return 0;

        long long val = _data[0];
        if (_size > 1) {

             return _sign == POSITIVE ? 2147483647 : -2147483648;
        }
        if (_sign == NEGATIVE) val = -val;
        return (int)val;
    }

    lmmc_real_t to_double() const {
        if (_size == 0) return 0.0;
        lmmc_real_t res = 0.0;
        lmmc_real_t base_mul = 18446744073709551616.0;

        for (mp_size_t i = _size; i > 0; --i) {
            LMMC_REAL_MUL(&res, &res, &base_mul);
            lmmc_real_t val = (lmmc_real_t)_data[i-1];
            LMMC_REAL_ADD(&res, &res, &val);
        }
        if (_sign == NEGATIVE) {
            lmmc_real_t zero = 0.0;
            LMMC_REAL_SUB(&res, &zero, &res);
        }
        return res;
    }

    static int cmp_abs(const BigInt& a, const BigInt& b) {
        if (a._size != b._size) return a._size > b._size ? 1 : -1;
        if (a._size == 0) return 0;
        for (mp_size_t i = a._size; i > 0; --i) {
             if (a._data[i-1] != b._data[i-1])
                 return a._data[i-1] > b._data[i-1] ? 1 : -1;
        }
        return 0;
    }

    bool operator==(const BigInt& other) const {
        return _sign == other._sign && cmp_abs(*this, other) == 0;
    }
    bool operator!=(const BigInt& other) const { return !(*this == other); }
    bool operator<(const BigInt& other) const {
        if (_sign != other._sign) return _sign < other._sign;
        if (_sign == ZERO) return false;
        int cmp = cmp_abs(*this, other);
        return _sign == POSITIVE ? cmp < 0 : cmp > 0;
    }
    bool operator>(const BigInt& other) const { return other < *this; }
    bool operator<=(const BigInt& other) const { return !(*this > other); }
    bool operator>=(const BigInt& other) const { return !(*this < other); }
    bool operator!() const { return _size == 0; }
    explicit operator bool() const { return _size != 0; }
    bool is_zero() const { return _size == 0; }

    std::vector<uint64_t> get_digits() const {
        std::vector<uint64_t> d;
        d.reserve(_size);
        for(mp_size_t i = 0; i < _size; ++i) {
            d.push_back(_data[i]);
        }
        return d;
    }

    std::size_t hash() const {
        std::size_t seed = 0;
        for (mp_size_t i = 0; i < _size; ++i) {
            seed ^= std::hash<uint64_t>{}(_data[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        if (negative) {
            seed ^= std::hash<int>{}(-1) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    BigInt Abs() const {
        BigInt ret = *this;
        if (ret._size > 0) {
            ret._sign = POSITIVE;
            ret.negative = false;
        }
        return ret;
    }

    bool IsNegative() const { return _sign == NEGATIVE; }

    BigInt negate() const {
        BigInt ret = *this;
        if (ret._size > 0) {
            ret._sign = (ret._sign == POSITIVE) ? NEGATIVE : POSITIVE;
            ret.negative = (ret._sign == NEGATIVE);
        }
        return ret;
    }

    BigInt operator-() const { return negate(); }

    static void add_abs(BigInt& dst, const BigInt& a, const BigInt& b) {
        mp_size_t n = std::max(a._size, b._size);
        dst.realloc_to(n + 1);

        mp_limb_t cy = 0;

        mp_srcptr ap = a._data;
        mp_srcptr bp = b._data;
        mp_size_t na = a._size;
        mp_size_t nb = b._size;

        if (na < nb) { std::swap(ap, bp); std::swap(na, nb); }

        if (nb > 0) {
             cy = lmmp_add_n_(dst._data, ap, bp, nb);
        }

        if (na > nb) {
             if (dst._data != ap) std::memcpy(dst._data + nb, ap + nb, (na - nb) * sizeof(mp_limb_t));

             mp_size_t k = nb;
             while (cy && k < na) {
                 dst._data[k]++;
                 cy = (dst._data[k] == 0);
                 k++;
             }
        }
        dst._data[na] = cy;
        dst._size = na + cy;
        dst.normalize();
    }

    static void sub_abs(BigInt& dst, const BigInt& a, const BigInt& b) {
        mp_size_t na = a._size;
        mp_size_t nb = b._size;
        dst.realloc_to(na);

        mp_limb_t bw = 0;
        if (nb > 0)
            bw = lmmp_sub_n_(dst._data, a._data, b._data, nb);

        if (na > nb) {
             if (dst._data != a._data) std::memcpy(dst._data + nb, a._data + nb, (na - nb) * sizeof(mp_limb_t));
             mp_size_t k = nb;
             while (bw && k < na) {
                 dst._data[k]--;
                 bw = (dst._data[k] == (mp_limb_t)-1);
                 k++;
             }
        }
        dst._size = na;
        dst.normalize();
    }

    BigInt operator+(const BigInt& other) const {
        if (_sign == ZERO) return other;
        if (other._sign == ZERO) return *this;

        BigInt res;
        if (_sign == other._sign) {
            add_abs(res, *this, other);
            res._sign = _sign;
            res.negative = (_sign == NEGATIVE);
        } else {

            int cmp = cmp_abs(*this, other);
            if (cmp == 0) {
                return BigInt(0);
            }
            if (cmp > 0) {
                sub_abs(res, *this, other);
                res._sign = _sign;
                res.negative = (_sign == NEGATIVE);
            } else {
                sub_abs(res, other, *this);
                res._sign = other._sign;
                res.negative = (res._sign == NEGATIVE);
            }
        }
        return res;
    }

    BigInt operator-(const BigInt& other) const {
        return *this + (-other);
    }

    BigInt operator*(const BigInt& other) const {
        if (_size == 0 || other._size == 0) return BigInt(0);

        BigInt res;
        mp_size_t na = _size;
        mp_size_t nb = other._size;
        res.realloc_to(na + nb);

        if (na >= nb) {
            bigint_mul_(res._data, _data, na, other._data, nb);
        } else {
            bigint_mul_(res._data, other._data, nb, _data, na);
        }

        res._size = na + nb;
        res._sign = (_sign == other._sign) ? POSITIVE : NEGATIVE;
        res.negative = (res._sign == NEGATIVE);
        res.normalize();
        return res;
    }

    BigInt operator/(const BigInt& other) const {
        if (other._size == 0) throw std::domain_error("Division by zero");
        if (_size < other._size) return BigInt(0);

        BigInt q, r;
        mp_size_t na = _size;
        mp_size_t nb = other._size;

        q.realloc_to(na - nb + 1);
        r.realloc_to(nb);

        bigint_div_(q._data, r._data, _data, na, other._data, nb);

        q._size = na - nb + 1;
        q._sign = (_sign == other._sign) ? POSITIVE : NEGATIVE;
        q.negative = (q._sign == NEGATIVE);
        q.normalize();
        return q;
    }

    BigInt operator%(const BigInt& other) const {
        if (other._size == 0) throw std::domain_error("Division by zero");
        if (_size < other._size) return *this;

        BigInt q;
        mp_size_t na = _size;
        mp_size_t nb = other._size;
        q.realloc_to(na - nb + 1);

        BigInt r;
        r.realloc_to(nb);

        bigint_div_(q._data, r._data, _data, na, other._data, nb);

        r._size = nb;
        r._sign = _sign;
        r.negative = (r._sign == NEGATIVE);
        r.normalize();
        return r;
    }

    BigInt& operator+=(const BigInt& other) { *this = *this + other; return *this; }
    BigInt& operator-=(const BigInt& other) { *this = *this - other; return *this; }
    BigInt& operator*=(const BigInt& other) { *this = *this * other; return *this; }
    BigInt& operator/=(const BigInt& other) { *this = *this / other; return *this; }
    BigInt& operator%=(const BigInt& other) { *this = *this % other; return *this; }

    BigInt power(unsigned long exp) const {
        if (_size == 0) return exp == 0 ? BigInt(1) : BigInt(0);

        BigInt res;
        mp_size_t needed = lmmp_pow_size_(_data, _size, exp);
        res.realloc_to(needed);

        res._size = lmmp_pow_(res._data, needed, _data, _size, exp);

        if (_sign == NEGATIVE && (exp & 1)) {
            res._sign = NEGATIVE;
            res.negative = true;
        } else {
            res._sign = POSITIVE;
            res.negative = false;
        }
        res.normalize();
        return res;
    }

    BigInt power(BigInt exp) const {
        if (exp._sign == NEGATIVE) throw std::domain_error("Negative exponent in integer power");
        if (exp._size == 0) return BigInt(1);

        if (exp._size <= 1) {
            return power((unsigned long)exp._data[0]);
        }

        BigInt base = *this;
        BigInt res(1);

        while (!exp.is_zero()) {
            if (exp._data[0] & 1) {
                res = res * base;
            }
            base = base * base;

             mp_limb_t carry = 0;
             for (mp_size_t i = exp._size; i > 0; --i) {
                 mp_limb_t cur = exp._data[i-1];
                 mp_limb_t next_carry = cur & 1;
                 exp._data[i-1] = (cur >> 1) | (carry << 63);
                 carry = next_carry;
             }
             exp.normalize();
        }
        return res;
    }

    BigInt sqrt() const {
        if (_sign == NEGATIVE) throw std::domain_error("Sqrt of negative number");
        if (_size == 0) return BigInt(0);
        if (*this == BigInt(1)) return BigInt(1);

        BigInt res;

        mp_size_t res_alloc = (_size / 2) + 2;
        res.realloc_to(res_alloc);

        lmmp_sqrt_(res._data, nullptr, _data, _size, 0);

        res._size = res_alloc;
        res._sign = POSITIVE;
        res.negative = false;
        res.normalize();

        return res;
    }

    bool is_odd() const {
        if (_size == 0) return false;
        return (_data[0] & 1);
    }

    bool is_even() const {
        return !is_odd();
    }

    mp_size_t trailing_zeros() const {
        if (is_zero()) return 0;
        mp_size_t count = 0;
        for (mp_size_t i = 0; i < _size; ++i) {
            if (_data[i] == 0) {
                count += LIMB_BITS;
            } else {
                count += lmmp_tailing_zeros_(_data[i]);
                break;
            }
        }
        return count;
    }

    BigInt& operator>>=(mp_size_t shift) {
        if (shift == 0) return *this;
        if (is_zero()) return *this;

        mp_size_t limb_shift = shift / LIMB_BITS;
        mp_size_t bit_shift = shift % LIMB_BITS;

        if (limb_shift >= _size) {
            _size = 0;
            _sign = POSITIVE;
            negative = false;
            return *this;
        }

        if (limb_shift > 0) {
            std::memmove(_data, _data + limb_shift, (_size - limb_shift) * sizeof(mp_limb_t));
            _size -= limb_shift;
        }

        if (bit_shift > 0) {
             lmmp_shr_(_data, _data, _size, bit_shift);
             if (_size > 0 && _data[_size-1] == 0) _size--;
        }
        normalize();
        return *this;
    }

    BigInt& operator<<=(mp_size_t shift) {
        if (shift == 0) return *this;
        if (is_zero()) return *this;

        mp_size_t limb_shift = shift / LIMB_BITS;
        mp_size_t bit_shift = shift % LIMB_BITS;

        mp_size_t old_size = _size;
        mp_size_t needed = old_size + limb_shift + (bit_shift > 0 ? 1 : 0);
        realloc_to(needed);

        if (bit_shift > 0) {
            mp_limb_t carry = lmmp_shl_(_data + limb_shift, _data, old_size, bit_shift);
            if (carry) {
                _data[old_size + limb_shift] = carry;
                _size = old_size + limb_shift + 1;
            } else {
                 _size = old_size + limb_shift;
            }
        } else {
             std::memmove(_data + limb_shift, _data, old_size * sizeof(mp_limb_t));
             _size += limb_shift;
        }

        if (limb_shift > 0) {
            std::memset(_data, 0, limb_shift * sizeof(mp_limb_t));
        }

        normalize();
        return *this;
    }

    BigInt operator>>(mp_size_t shift) const {
        BigInt res = *this;
        res >>= shift;
        return res;
    }

    BigInt operator<<(mp_size_t shift) const {
        BigInt res = *this;
        res <<= shift;
        return res;
    }

    static BigInt factorial(unsigned int n) {
        BigInt res;
        mp_bitcnt_t bits = 0;
        mp_size_t needed = lmmp_factorial_size_(n, &bits);
        res.realloc_to(needed);
        res._size = lmmp_factorial_(res._data, bits, needed, n);
        res._sign = POSITIVE;
        res.negative = false;
        res.normalize();
        return res;
    }

    static BigInt nPr(unsigned int n, unsigned int r) {
        if (r > n) return BigInt(0);
        BigInt res;
        mp_bitcnt_t bits = 0;
        mp_size_t needed = lmmp_nPr_size_(n, r, &bits);
        res.realloc_to(needed);
        res._size = lmmp_nPr_(res._data, bits, needed, n, r);
        res._sign = POSITIVE;
        res.negative = false;
        res.normalize();
        return res;
    }

    static BigInt nCr(unsigned int n, unsigned int r) {
        if (r > n) return BigInt(0);
        BigInt res;
        mp_bitcnt_t bits = 0;
        mp_size_t needed = lmmp_nCr_size_(n, r, &bits);
        res.realloc_to(needed);
        res._size = lmmp_nCr_(res._data, bits, needed, n, r);
        res._sign = POSITIVE;
        res.negative = false;
        res.normalize();
        return res;
    }

    static BigInt multinomial(unsigned int n, const std::vector<unsigned int>& r) {
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
        res._size = lmmp_multinomial_(res._data, needed, (uint)sum, r_uints.data(), (uint)r_uints.size());
        res._sign = POSITIVE;
        res.negative = false;
        res.normalize();
        return res;
    }

    static BigInt gcd(const BigInt& a, const BigInt& b) {
        if (a.is_zero()) return b.Abs();
        if (b.is_zero()) return a.Abs();

        BigInt abs_a = a.Abs();
        BigInt abs_b = b.Abs();

        BigInt res;
        mp_size_t na = abs_a._size;
        mp_size_t nb = abs_b._size;

        mp_size_t min_n = (na < nb) ? na : nb;
        res.realloc_to(min_n);

        res._size = lmmp_gcd_lehmer_(res._data, abs_a._data, na, abs_b._data, nb);

        res._sign = POSITIVE;
        res.negative = false;
        res.normalize();
        return res;
    }

    static BigInt lcm(const BigInt& a, const BigInt& b) {
        if (a.is_zero() || b.is_zero()) return BigInt(0);
        return (a.Abs() / gcd(a, b)) * b.Abs();
    }

    static BigInt pow_mod(const BigInt& base, const BigInt& exp, const BigInt& mod) {
         if (mod.is_zero()) throw std::runtime_error("Modulo by zero");

         if (base._size <= 1 && exp._size <= 1 && mod._size <= 1) {
             ulong b = base.is_zero() ? 0 : base._data[0];
             ulong e = exp.is_zero() ? 0 : exp._data[0];
             ulong m = mod._data[0];

             return BigInt(lmmp_powmod_ulong_(b % m, e, m));
         }

         BigInt res = 1;
         BigInt b = base % mod;
         BigInt e = exp;

         while (!e.is_zero()) {
             if (e.is_odd()) res = (res * b) % mod;
             b = (b * b) % mod;
             e >>= 1;
         }
         return res;
    }

    bool is_prime() const {
        if (_sign == NEGATIVE) return false;

        if (_size <= 1) {
             return lmmp_is_prime_ulong_(_size == 0 ? 0 : _data[0]);
        }

        static const uint64_t mr_bases[] = {2, 3, 5, 7, 11, 13, 17};

        static const uint64_t small_primes[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
        for (auto p : small_primes) {
            if ((*this % p).is_zero()) return false;
        }

        BigInt n = *this;
        BigInt d = n - 1;
        BigInt two(2);
        int s = 0;
        while (d.is_even()) {
            d = d / two;
            s++;
        }

        const int num_bases = sizeof(mr_bases) / sizeof(mr_bases[0]);
        for (int i = 0; i < num_bases; ++i) {
            BigInt a(mr_bases[i]);
            if (a >= n) continue;

            BigInt x = BigInt::pow_mod(a, d, n);
            if (x == 1 || x == n - 1) continue;

            bool composite = true;
            for (int r = 1; r < s; ++r) {
                x = BigInt::pow_mod(x, two, n);
                if (x == n - 1) {
                    composite = false;
                    break;
                }
            }
            if (composite) return false;
        }
        return true;
    }

    bool is_perfect_square() const {
        if (_sign == NEGATIVE) return false;
        if (_size == 0) return true;
        BigInt s = this->sqrt();
        return (s * s) == *this;
    }
};
