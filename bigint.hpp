#pragma once

#include "LAMMP/include/lammp/lmmp.h"
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <limits> // for std::numeric_limits

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
    bool negative = false; // Compatibility with legacy code accessing .negative member

    // Helper: Allocate memory
    void realloc_to(mp_size_t new_alloc) {
        if (new_alloc <= _alloc) return;
        new_alloc = (new_alloc + 3) & ~3; // Align to 4
        mp_ptr new_data = (mp_ptr)lmmp_alloc(new_alloc * sizeof(mp_limb_t));
        if (!new_data) throw std::bad_alloc();
        
        // Initialize new memory to zero
        // std::memset(new_data, 0, new_alloc * sizeof(mp_limb_t));

        if (_size > 0 && _data) {
             std::memcpy(new_data, _data, _size * sizeof(mp_limb_t));
        }
        if (_data) lmmp_free(_data);
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
        if (_data) lmmp_free(_data);
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
            if (_data) lmmp_free(_data);
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
            // Handle INT64_MIN carefully
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
        // Estimate limbs needed (safe upper bound)
        mp_size_t needed = len / 19 + 2; 
        realloc_to(needed);

        // Create a temporary buffer with digit values (0-9), reversed (Little Endian)
        std::vector<mp_byte_t> digit_buf(len);
        for (size_t i = 0; i < len; ++i) {
            char c = str[start + i];
            uint8_t d = 0;
            if (c >= '0' && c <= '9') {
                d = c - '0';
            }
            // Store in reverse order (LSD at index 0)
            digit_buf[len - 1 - i] = d; 
        }

        // lmmp_from_str_ handles base conversion (expects LE digit array)
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
        // 1 limb ~ 19.3 digits
        size_t len_needed = _size * 20 + 5; 
        std::vector<mp_byte_t> buf(len_needed);
        
        mp_size_t str_len = lmmp_to_str_((mp_byte_t*)buf.data(), _data, _size, 10);
        
        std::string res;
        res.reserve(str_len + 2);
        if (_sign == NEGATIVE) res += '-';
        // lmmp_to_str_ returns digits in Little Endian (LSD at index 0).
        // We need to print MSD first, so reverse iteration.
        for(mp_size_t i = str_len; i > 0; --i) {
            res += (char)(buf[i - 1] + '0');
        }
        return res;
    }
    std::string to_string() const { return ToString(); }

    int to_int() const {
        if (_size == 0) return 0;
        // Check overflow?
        long long val = _data[0];
        if (_size > 1) { // Oversimplified check
             // saturated
             return _sign == POSITIVE ? 2147483647 : -2147483648;
        }
        if (_sign == NEGATIVE) val = -val;
        return (int)val;
    }

    double to_double() const {
        if (_size == 0) return 0.0;
        double res = 0.0;
        // Convert limbs to double
        for (mp_size_t i = _size; i > 0; --i) {
            res = res * 18446744073709551616.0; // 2^64
            res += (double)_data[i-1];
        }
        return _sign == NEGATIVE ? -res : res;
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

    // Comparison
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
    
    // Compatibility with legacy code expecting vector of digits
    std::vector<uint64_t> get_digits() const {
        std::vector<uint64_t> d;
        d.reserve(_size);
        for(mp_size_t i = 0; i < _size; ++i) {
            d.push_back(_data[i]);
        }
        return d;
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

    // Arithmetic
    // Note: These implementations assume basic understanding of lmmp functions
    // lmmp_add_n_ requires limbs, not BigInts.
    
    // Unsigned addition: dst = |a| + |b|
    static void add_abs(BigInt& dst, const BigInt& a, const BigInt& b) {
        mp_size_t n = std::max(a._size, b._size);
        dst.realloc_to(n + 1);
        
        // Temporarily clear upper part of dst if realloc didn't zero it
        // and copy data if needed. simpler strategy: use buffer
        
        mp_limb_t cy = 0;
        // Standard add loop since lmmp might not support uneven sizes directly with one call
        
        mp_srcptr ap = a._data;
        mp_srcptr bp = b._data;
        mp_size_t na = a._size;
        mp_size_t nb = b._size;
        
        // Ensure na >= nb for convenience
        if (na < nb) { std::swap(ap, bp); std::swap(na, nb); }
        
        // Now adding [ap, na] + [bp, nb]
        if (nb > 0) {
             cy = lmmp_add_n_(dst._data, ap, bp, nb);
        }
        // Propagate carry through the rest of A
        if (na > nb) {
             // add carry to a[nb...na]
             // lmmp_add_1 or similar could be useful, or manual
             // lmmp_add_nc_ is [dst, n] = [src1, n] + [src2, n] + c
             // We can treat it as adding 0 + c ?? No.
             
             // Just Manual copy + carry prop
             mp_ptr dp = dst._data + nb;
             mp_srcptr asp = ap + nb;
             mp_size_t count = na - nb;
             // Add carry to first limb
             mp_limb_t r = asp[0] + cy;
             dp[0] = r;
             cy = (r < asp[0]); // Carry Occurred? 
             if (r < asp[0] || (cy && r == asp[0])) {} // Check correctness of carry check
             // Actually: cy_out = (sum < op1)
             if (cy) {
                 // Propagate
                 for (mp_size_t i = 1; i < count; ++i) {
                     r = asp[i] + 1;
                     dp[i] = r;
                     if (r != 0) { cy = 0; break; } // No more carry
                 }
                 if (cy) {
                      // Remaining copy
                      if (count > 0) {
                           // If we broke loop, we need to copy rest? 
                           // wait, loop logic is tricky.
                      }
                 }
             } 
             // Better: Use lmmp logic. 
             // We can use a zero array? No.
             // Let's implement manually for safety or look closer at lmmp.h
             // Actually, simplest is to use a scratch space of size na for B, padded with 0.
             // But valid for performance? 
             
             // Re-reading lmmp.h... "lmmp_add_nc_"
             // Can we just use logic: 
             // Copy A to Dest. Add B to Dest.
             if (dst._data != ap) std::memcpy(dst._data + nb, ap + nb, (na - nb) * sizeof(mp_limb_t));
             
             // Propagate carry manually on dst
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
    
    // Unsigned subtraction: dst = |a| - |b|. Assumes |a| >= |b|
    static void sub_abs(BigInt& dst, const BigInt& a, const BigInt& b) {
        // ... similar to add_abs ...
        mp_size_t na = a._size;
        mp_size_t nb = b._size;
        dst.realloc_to(na);
        
        mp_limb_t bw = 0;
        if (nb > 0)
            bw = lmmp_sub_n_(dst._data, a._data, b._data, nb);
        
        // Propagate borrow
        if (na > nb) {
             if (dst._data != a._data) std::memcpy(dst._data + nb, a._data + nb, (na - nb) * sizeof(mp_limb_t));
             mp_size_t k = nb;
             while (bw && k < na) {
                 dst._data[k]--;
                 bw = (dst._data[k] == (mp_limb_t)-1); // Borrowed wrapped around
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
            // Different signs: abs(A) - abs(B)
            int cmp = cmp_abs(*this, other);
            if (cmp == 0) {
                return BigInt(0);
            }
            if (cmp > 0) {
                sub_abs(res, *this, other);
                res._sign = _sign; // Sign of larger
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
        
        // lmmp_mul_ handles unequal sizes: [dst, na+nb] = [a, na] * [b, nb]
        // Note: lmmp_mul_ expects 0 < nb <= na. Swap if needed.
        if (na >= nb) {
            lmmp_mul_(res._data, _data, na, other._data, nb);
        } else {
            lmmp_mul_(res._data, other._data, nb, _data, na);
        }
        
        res._size = na + nb;
        res._sign = (_sign == other._sign) ? POSITIVE : NEGATIVE;
        res.negative = (res._sign == NEGATIVE);
        res.normalize();
        return res;
    }

    BigInt operator/(const BigInt& other) const {
        if (other._size == 0) throw std::runtime_error("Division by zero");
        if (_size < other._size) return BigInt(0); // Integer division
        
        BigInt q, r;
        mp_size_t na = _size;
        mp_size_t nb = other._size;
        
        // q size: na - nb + 1
        q.realloc_to(na - nb + 1);
        r.realloc_to(nb); // remainder size is nb
        
        // Assume q._data, r._data, a._data, b._data do not overlap in bad ways
        // Allocate temp buffer for A because lmmp_div_ might modify it? 
        // Doc says: [dstq] = [numa]/[numb]. numa is srcptr.
        
        // IMPORTANT: lmmp_div_ usually requires numa to be copied if it's destructively modified 
        // BUT here params are (mp_ptr dstq, mp_ptr dstr, mp_srcptr numa, ...) 
        // so numa is const.
        
        lmmp_div_(q._data, r._data, _data, na, other._data, nb);
        
        q._size = na - nb + 1;
        q._sign = (_sign == other._sign) ? POSITIVE : NEGATIVE;
        q.negative = (q._sign == NEGATIVE);
        q.normalize();
        return q;
    }
    
    BigInt operator%(const BigInt& other) const {
        if (other._size == 0) throw std::runtime_error("Division by zero");
        if (_size < other._size) return *this; 

        BigInt r;
        mp_size_t na = _size;
        mp_size_t nb = other._size;
        r.realloc_to(nb);
        
        // Pass NULL for quotient
        lmmp_div_(nullptr, r._data, _data, na, other._data, nb);
        
        r._size = nb;
        r._sign = _sign; // Remainder sign follows dividend
        r.negative = (r._sign == NEGATIVE);
        r.normalize();
        return r;
    }

    BigInt& operator+=(const BigInt& other) { *this = *this + other; return *this; }
    BigInt& operator-=(const BigInt& other) { *this = *this - other; return *this; }
    BigInt& operator*=(const BigInt& other) { *this = *this * other; return *this; }
    BigInt& operator/=(const BigInt& other) { *this = *this / other; return *this; }
    BigInt& operator%=(const BigInt& other) { *this = *this % other; return *this; }

    // Power
    BigInt power(BigInt exp) const {
        if (exp._sign == NEGATIVE) throw std::domain_error("Negative exponent in integer power");
        if (exp._size == 0) return BigInt(1);
        
        BigInt base = *this;
        BigInt res(1);
        
        // Binary exponentiation
        // We can inspect bits of exp
        // But exp is BigInt. 
        // Simple loop
        while (!exp.is_zero()) {
            if (exp._data[0] & 1) {
                res = res * base;
            }
            base = base * base;
            // exp >>= 1
             // Right shift BigInt is not implemented yet. 
             // Implement simple div 2 or shift.
             // lmmp_shr1... functions exist!
             // lmmp_shr1sub_n_ is complex. 
             // Just use div 2 for now, optimizing later
             // Or implement operator>>=
             
             // Optimized shift right 1
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

    // Sqrt (Integer sqrt)
    BigInt sqrt() const {
        if (_sign == NEGATIVE) throw std::domain_error("Sqrt of negative number");
        if (_size == 0) return BigInt(0);
        
        BigInt res;
        // Output size: (na + 1) / 2
        mp_size_t res_size = (_size + 1) / 2 + 1; // +1 buffer
        res.realloc_to(res_size);
        
        // lmmp_sqrt_(mp_ptr dsts, mp_ptr dstr, mp_srcptr numa, mp_size_t na, mp_size_t nf)
        // nf is precision factor, 0 for int sqrt?
        // Doc says: [dsts] = floor(sqrt([numa] * B^(2*nf)))
        // So use nf=0 for integer sqrt
        
        lmmp_sqrt_(res._data, nullptr, _data, _size, 0);
        
        res._size = res_size; 
        res.normalize();
        res._sign = POSITIVE;
        res.negative = false;
        return res;
    }

    bool is_perfect_square() const {
        if (_sign == NEGATIVE) return false;
        if (_size == 0) return true;
        
        // Calc sqrt, then square it to check ?? Expensive.
        // Or look at dstr from lmmp_sqrt_
        // lmmp_sqrt_ returns remainder if dstr is not null!
        BigInt s, r;
        mp_size_t res_size = (_size + 1) / 2 + 1;
        s.realloc_to(res_size);
        r.realloc_to(res_size + 2); // Remainder can be larger? 
        // Remainder bounded by 2*sqrt(n)
        
        lmmp_sqrt_(s._data, r._data, _data, _size, 0);
        
        // If remainder r is zero, it's perfect square
        // Check if r is zero from data
        // We need to set r._size properly to check?
        // lmmp functions don't set separate size return usually for array like this?
        // Wait, lmmp_sqrt declares: [dstr, nf+na/2+1]
        // We need to check all limbs of dstr
        mp_size_t r_len = 0 + _size/2 + 1;
        for(mp_size_t i=0; i<r_len; ++i) {
            if (r._data[i] != 0) return false;
        }
        return true;
    }
};
