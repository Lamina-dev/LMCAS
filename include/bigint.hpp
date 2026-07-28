/**
 * @file bigint.hpp
 * @brief 任意精度整数 BigInt，底层调用 LAMMP 实现高性能大数运算。
 */
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
#include <cstdint>
#include <cstring>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <optional>
#include "lmmc/config.h"

/** @brief 释放 BigInt 内部分配的内存 */
inline void bigint_free(void* p) {
    lmmp_free(p);
}

/** @brief 为 BigInt 分配内存 */
inline void* bigint_alloc(size_t size) {
    return lmmp_alloc(size);
}

/**
 * @brief 去除 limb 数组的前导零，返回有效长度
 * @param p limb 数组指针
 * @param n 原始长度
 * @return 去除前导零后的有效长度
 */
inline mp_size_t lmmp_rlz(mp_srcptr p, mp_size_t n) {
    while (n > 0 && p[n-1] == 0) n--;
    return n;
}

/**
 * @brief 大整数乘法（底层 LAMMP 调用）
 * @param dst 结果缓冲区
 * @param numa 乘数 a 的 limb 数组
 * @param na 乘数 a 的 limb 数
 * @param numb 乘数 b 的 limb 数组
 * @param nb 乘数 b 的 limb 数
 */
inline void bigint_mul_(mp_ptr dst, mp_srcptr numa, mp_size_t na, mp_srcptr numb, mp_size_t nb) {
    ::lmmp_mul_(dst, numa, na, numb, nb);
}

/**
 * @brief 大整数除法（底层 LAMMP 调用）
 * @param q 商缓冲区
 * @param r 余数缓冲区
 * @param n 被除数 limb 数组
 * @param nn 被除数 limb 数
 * @param d 除数 limb 数组
 * @param dn 除数 limb 数
 */
inline void bigint_div_(mp_ptr q, mp_ptr r, mp_srcptr n, mp_size_t nn, mp_srcptr d, mp_size_t dn) {
    ::lmmp_div_(q, r, n, nn, d, dn);
}

/** @brief 任意精度整数类，底层使用 LAMMP 多精度库 */
class BigInt {
public:
    /** @brief 符号枚举 */
    enum Sign {
        POSITIVE = 1,
        ZERO = 0,
        NEGATIVE = -1
    };

private:
    mp_ptr _data = nullptr;
    mp_size_t _size = 0;
    mp_size_t _alloc = 0;
    int _sign = ZERO;
    bool negative = false;

    void realloc_to(mp_size_t new_alloc) {
        if (new_alloc <= _alloc) return;
        const std::size_t requested = static_cast<std::size_t>(new_alloc);
        if (requested > std::numeric_limits<std::size_t>::max() - 3) {
            throw std::length_error("BigInt allocation size overflow");
        }
        const std::size_t rounded = (requested + 3) & ~std::size_t(3);
        if (rounded > std::numeric_limits<std::size_t>::max() / sizeof(mp_limb_t) ||
            rounded > static_cast<std::size_t>(std::numeric_limits<mp_size_t>::max())) {
            throw std::length_error("BigInt allocation byte size overflow");
        }

        mp_ptr new_data = static_cast<mp_ptr>(bigint_alloc(rounded * sizeof(mp_limb_t)));
        if (!new_data) throw std::bad_alloc();

        if (_size > 0 && _data) {
             std::memcpy(new_data, _data, _size * sizeof(mp_limb_t));
        }
        if (_data) bigint_free(_data);
        _data = new_data;
        _alloc = static_cast<mp_size_t>(rounded);
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
    /** @brief 默认构造，值为 0 */
    BigInt() : _data(nullptr), _size(0), _alloc(0), _sign(ZERO), negative(false) {}

    ~BigInt() {
        if (_data) bigint_free(_data);
    }

    /** @brief 拷贝构造 */
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

    /** @brief 移动构造 */
    BigInt(BigInt&& other) noexcept
        : _data(other._data), _size(other._size), _alloc(other._alloc), _sign(other._sign), negative(other.negative) {
        other._data = nullptr;
        other._size = 0;
        other._alloc = 0;
        other._sign = ZERO;
        other.negative = false;
    }

    /** @brief 拷贝赋值 */
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

    /** @brief 移动赋值 */
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

    /**
     * @brief 从 long long 构造
     * @param val 整数值
     */
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
    /** @brief 从 int 构造 */
    BigInt(int val) : BigInt((long long)val) {}

    /**
     * @brief 从 unsigned long long 构造
     * @param val 无符号整数值
     */
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
    /** @brief 从 unsigned int 构造 */
    BigInt(unsigned int val) : BigInt((unsigned long long)val) {}
    /** @brief 从 unsigned long 构造 */
    BigInt(unsigned long val) : BigInt((unsigned long long)val) {}

    /**
     * @brief 从十进制字符串构造
     * @param str 十进制数字字符串，可带正负号
     */
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
            if (c < '0' || c > '9') {
                throw std::invalid_argument(
                    "BigInt: invalid character '" + std::string(1, c) +
                    "' at position " + std::to_string(start + i) +
                    " in input \"" + str + "\"");
            }
            uint8_t d = static_cast<uint8_t>(c - '0');

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

    /**
     * @brief 转换为十进制字符串
     * @return 十进制表示的字符串
     */
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

    /** @brief 转换为十进制字符串（同 ToString） */
    std::string to_string() const { return ToString(); }

    /**
     * @brief 转换为 int（溢出时饱和到 int 范围）
     * @return 对应的 int 值
     */
    int to_int() const {
        if (_size == 0) return 0;

        constexpr long long INT_MAX_LL = static_cast<long long>(std::numeric_limits<int>::max());
        if (_size > 1) {
            return _sign == POSITIVE
                       ? std::numeric_limits<int>::max()
                       : std::numeric_limits<int>::min();
        }

        // _size == 1: a single limb. Saturate when the magnitude does not fit.
        unsigned long long mag = _data[0];
        if (_sign == NEGATIVE) {
            // Negative: representable range is [INT_MIN, 0].
            // |INT_MIN| = static_cast<unsigned long long>(INT_MAX) + 1.
            unsigned long long min_mag =
                static_cast<unsigned long long>(INT_MAX_LL) + 1ULL;
            if (mag > min_mag) return std::numeric_limits<int>::min();
            if (mag == min_mag) return std::numeric_limits<int>::min();
            return static_cast<int>(-static_cast<long long>(mag));
        } else {
            if (mag > static_cast<unsigned long long>(INT_MAX_LL)) {
                return std::numeric_limits<int>::max();
            }
            return static_cast<int>(mag);
        }
    }

    /** @brief 精确转换为 int64_t；超出范围时返回空。 */
    std::optional<std::int64_t> try_to_int64() const noexcept {
        if (_size == 0) return std::int64_t{0};
        if (_size > 1) return std::nullopt;

        const unsigned long long magnitude =
            static_cast<unsigned long long>(_data[0]);
        const auto max_value =
            static_cast<unsigned long long>(std::numeric_limits<std::int64_t>::max());
        if (_sign == NEGATIVE) {
            const unsigned long long min_magnitude = max_value + 1ULL;
            if (magnitude > min_magnitude) return std::nullopt;
            if (magnitude == min_magnitude) return std::numeric_limits<std::int64_t>::min();
            return -static_cast<std::int64_t>(magnitude);
        }
        if (magnitude > max_value) return std::nullopt;
        return static_cast<std::int64_t>(magnitude);
    }

    /**
     * @brief 转换为浮点数
     * @return 对应的 double 值
     */
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

    /**
     * @brief 比较两个 BigInt 的绝对值
     * @param a 第一个操作数
     * @param b 第二个操作数
     * @return 1 表示 |a|>|b|，-1 表示 |a|<|b|，0 表示相等
     */
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

    /** @brief 判断是否为零 */
    bool is_zero() const { return _size == 0; }

    /**
     * @brief 获取 limb 数组的副本
     * @return 各 limb 值组成的 vector
     */
    std::vector<uint64_t> get_digits() const {
        std::vector<uint64_t> d;
        d.reserve(_size);
        for(mp_size_t i = 0; i < _size; ++i) {
            d.push_back(_data[i]);
        }
        return d;
    }

    /**
     * @brief 计算哈希值
     * @return 哈希值
     */
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

    /**
     * @brief 取绝对值
     * @return |*this|
     */
    BigInt Abs() const {
        BigInt ret = *this;
        if (ret._size > 0) {
            ret._sign = POSITIVE;
            ret.negative = false;
        }
        return ret;
    }

    /** @brief 判断是否为负数 */
    bool IsNegative() const { return _sign == NEGATIVE; }

    /**
     * @brief 取相反数
     * @return -*this
     */
    BigInt negate() const {
        BigInt ret = *this;
        if (ret._size > 0) {
            ret._sign = (ret._sign == POSITIVE) ? NEGATIVE : POSITIVE;
            ret.negative = (ret._sign == NEGATIVE);
        }
        return ret;
    }

    /** @brief 一元负号运算符 */
    BigInt operator-() const { return negate(); }

    /**
     * @brief 绝对值加法（内部使用）
     * @param dst 结果
     * @param a 加数
     * @param b 加数
     */
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

    /**
     * @brief 绝对值减法（内部使用，要求 |a| >= |b|）
     * @param dst 结果
     * @param a 被减数
     * @param b 减数
     */
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

    /** @brief 加法运算符 */
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

    /** @brief 减法运算符 */
    BigInt operator-(const BigInt& other) const {
        return *this + (-other);
    }

    /** @brief 乘法运算符 */
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

    /**
     * @brief 除法运算符（整数除法，截断）
     * @param other 除数
     * @return 商
     * @throw std::domain_error 除数为零时抛出
     */
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

    /**
     * @brief 取模运算符
     * @param other 模数
     * @return 余数
     * @throw std::domain_error 模数为零时抛出
     */
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

    /**
     * @brief 整数幂运算
     * @param exp 非负指数
     * @return this^exp
     */
    BigInt power(unsigned long exp) const {
        if (exp == 0) return BigInt(1);
        if (_size == 0) return BigInt(0);

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

    /**
     * @brief 整数幂运算（BigInt 指数）
     * @param exp 非负 BigInt 指数
     * @return this^exp
     * @throw std::domain_error 指数为负时抛出
     */
    BigInt power(BigInt exp) const {
        if (exp._sign == NEGATIVE) throw std::domain_error("Negative exponent in integer power");
        if (exp._size == 0) return BigInt(1);

        if (exp._size <= 1) {
            // 单 limb 时尽量走窄重载，但 mp_limb_t 的位宽未必等于 unsigned long
            // （例如 LLP64 的 Windows 上 unsigned long 是 32-bit），所以仅在 limb 值
            // 能完整放进 unsigned long 时才走窄路径，否则继续走通用循环避免截断。
            mp_limb_t lo = (exp._size == 0) ? 0 : exp._data[0];
            if (lo <= static_cast<mp_limb_t>(std::numeric_limits<unsigned long>::max())) {
                return power(static_cast<unsigned long>(lo));
            }
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
                 // 不要写死 limb=64：用 LIMB_BITS 计算插入位置，且在 mp_limb_t 类型下移位。
                 exp._data[i-1] = (cur >> 1) | (carry << (LIMB_BITS - 1));
                 carry = next_carry;
             }
             exp.normalize();
        }
        return res;
    }

    /**
     * @brief 整数平方根（向下取整）
     * @return floor(sqrt(*this))
     * @throw std::domain_error 负数时抛出
     */
    BigInt sqrt() const {
        if (_sign == NEGATIVE) throw std::domain_error("Sqrt of negative number");
        if (_size == 0) return BigInt(0);
        if (*this == BigInt(1)) return BigInt(1);

        BigInt res;

        const mp_size_t root_size = (_size + 1) / 2;
        // lmmp_sqrt_ may write one guard limb for an even input limb count.
        res.realloc_to(root_size + 1);

        lmmp_sqrt_(res._data, nullptr, _data, _size, 0);

        res._size = root_size;
        res._sign = POSITIVE;
        res.negative = false;
        res.normalize();

        return res;
    }

    /** @brief 判断是否为奇数 */
    bool is_odd() const {
        if (_size == 0) return false;
        return (_data[0] & 1);
    }

    /** @brief 判断是否为偶数 */
    bool is_even() const {
        return !is_odd();
    }

    /**
     * @brief 计算末尾零比特数
     * @return 二进制表示中末尾连续 0 的个数
     */
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

    /**
     * @brief 右移赋值
     * @param shift 移位比特数
     * @return *this
     */
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

    /**
     * @brief 左移赋值
     * @param shift 移位比特数
     * @return *this
     */
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

    /** @brief 右移运算符 */
    BigInt operator>>(mp_size_t shift) const {
        BigInt res = *this;
        res >>= shift;
        return res;
    }

    /** @brief 左移运算符 */
    BigInt operator<<(mp_size_t shift) const {
        BigInt res = *this;
        res <<= shift;
        return res;
    }

    /**
     * @brief 计算阶乘
     * @param n 非负整数
     * @return n!
     */
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

    /**
     * @brief 计算排列数 P(n, r)
     * @param n 总数
     * @param r 选取数
     * @return n! / (n-r)!
     */
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

    /**
     * @brief 计算组合数 C(n, r)
     * @param n 总数
     * @param r 选取数
     * @return n! / (r! * (n-r)!)
     */
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

    /**
     * @brief 计算多项式系数（多重组合数）
     * @param n 总数（应等于 r 各元素之和）
     * @param r 各组大小
     * @return n! / (r[0]! * r[1]! * ... * r[k]!)
     * @throw std::invalid_argument r 之和不等于 n 时抛出
     */
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

    /**
     * @brief 计算最大公约数
     * @param a 第一个大整数
     * @param b 第二个大整数
     * @return gcd(|a|, |b|)
     */
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

    /**
     * @brief 计算最小公倍数
     * @param a 第一个大整数
     * @param b 第二个大整数
     * @return lcm(|a|, |b|)
     */
    static BigInt lcm(const BigInt& a, const BigInt& b) {
        if (a.is_zero() || b.is_zero()) return BigInt(0);
        return (a.Abs() / gcd(a, b)) * b.Abs();
    }

    /**
     * @brief 模幂运算
     * @param base 底数
     * @param exp 指数
     * @param mod 模数
     * @return base^exp mod mod
     * @throw std::runtime_error 模数为零时抛出
     */
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

    /**
     * @brief 素性测试（Miller-Rabin）
     * @return 是素数返回 true
     */
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

    /**
     * @brief 判断是否为完全平方数
     * @return 是完全平方数返回 true
     */
    bool is_perfect_square() const {
        if (_sign == NEGATIVE) return false;
        if (_size == 0) return true;
        BigInt s = this->sqrt();
        return (s * s) == *this;
    }
};
