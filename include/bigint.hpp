/**
 * @file bigint.hpp
 * @brief 任意精度整数 BigInt，底层调用 LMMP 实现高性能大数运算。
 */
#pragma once

#ifndef _STATIC_ASSERT
#define _STATIC_ASSERT(x) static_assert(x, #x)
#endif

#include <lmmp.h>
#include <lmmpn.h>
#include <numth.h>
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
#include "lamina_export.hpp"

class LAMINA_API BigInt {
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

    void realloc_to(mp_size_t new_alloc);

    void normalize();

    void zero();

public:
    /** @brief 默认构造，值为 0 */
    BigInt();

    ~BigInt();

    /** @brief 拷贝构造 */
    BigInt(const BigInt& other);

    /** @brief 移动构造 */
    BigInt(BigInt&& other) noexcept;

    /** @brief 拷贝赋值 */
    BigInt& operator=(const BigInt& other);

    /** @brief 移动赋值 */
    BigInt& operator=(BigInt&& other) noexcept;

    /**
     * @brief 从 long long 构造
     * @param val 整数值
     */
    BigInt(long long val);
    /** @brief 从 int 构造 */
    BigInt(int val);

    /**
     * @brief 从 unsigned long long 构造
     * @param val 无符号整数值
     */
    BigInt(unsigned long long val);
    /** @brief 从 unsigned int 构造 */
    BigInt(unsigned int val);
    /** @brief 从 unsigned long 构造 */
    BigInt(unsigned long val);

    /**
     * @brief 从十进制字符串构造
     * @param str 十进制数字字符串，可带正负号
     */
    BigInt(const std::string& str);

    /**
     * @brief 转换为十进制字符串
     * @return 十进制表示的字符串
     */
    std::string ToString() const;

    /** @brief 转换为十进制字符串（同 ToString） */
    std::string to_string() const;

    /**
     * @brief 转换为 int（溢出时饱和到 int 范围）
     * @return 对应的 int 值
     */
    int to_int() const;

    /** @brief 精确转换为 int64_t；超出范围时返回空。 */
    std::optional<std::int64_t> try_to_int64() const noexcept;

    /**
     * @brief 转换为浮点数
     * @return 对应的 double 值
     */
    lmmc_real_t to_double() const;

    /**
     * @brief 比较两个 BigInt 的绝对值
     * @param a 第一个操作数
     * @param b 第二个操作数
     * @return 1 表示 |a|>|b|，-1 表示 |a|<|b|，0 表示相等
     */
    static int cmp_abs(const BigInt& a, const BigInt& b);

    bool operator==(const BigInt& other) const;
    bool operator!=(const BigInt& other) const;
    bool operator<(const BigInt& other) const;
    bool operator>(const BigInt& other) const;
    bool operator<=(const BigInt& other) const;
    bool operator>=(const BigInt& other) const;
    bool operator!() const;
    explicit operator bool() const;

    /** @brief 判断是否为零 */
    bool is_zero() const;

    /**
     * @brief 获取 limb 数组的副本
     * @return 各 limb 值组成的 vector
     */
    std::vector<uint64_t> get_digits() const;

    /**
     * @brief 计算哈希值
     * @return 哈希值
     */
    std::size_t hash() const;

    /**
     * @brief 取绝对值
     * @return |*this|
     */
    BigInt Abs() const;

    /** @brief 判断是否为负数 */
    bool IsNegative() const;

    /**
     * @brief 取相反数
     * @return -*this
     */
    BigInt negate() const;

    /** @brief 一元负号运算符 */
    BigInt operator-() const;

    /**
     * @brief 绝对值加法（内部使用）
     * @param dst 结果
     * @param a 加数
     * @param b 加数
     */
    static void add_abs(BigInt& dst, const BigInt& a, const BigInt& b);

    /**
     * @brief 绝对值减法（内部使用，要求 |a| >= |b|）
     * @param dst 结果
     * @param a 被减数
     * @param b 减数
     */
    static void sub_abs(BigInt& dst, const BigInt& a, const BigInt& b);

    /** @brief 加法运算符 */
    BigInt operator+(const BigInt& other) const;

    /** @brief 减法运算符 */
    BigInt operator-(const BigInt& other) const;

    /** @brief 乘法运算符 */
    BigInt operator*(const BigInt& other) const;

    /**
     * @brief 除法运算符（整数除法，截断）
     * @param other 除数
     * @return 商
     * @throw std::domain_error 除数为零时抛出
     */
    BigInt operator/(const BigInt& other) const;

    /**
     * @brief 取模运算符
     * @param other 模数
     * @return 余数
     * @throw std::domain_error 模数为零时抛出
     */
    BigInt operator%(const BigInt& other) const;

    BigInt& operator+=(const BigInt& other);
    BigInt& operator-=(const BigInt& other);
    BigInt& operator*=(const BigInt& other);
    BigInt& operator/=(const BigInt& other);
    BigInt& operator%=(const BigInt& other);

    /**
     * @brief 整数幂运算
     * @param exp 非负指数
     * @return this^exp
     */
    BigInt power(unsigned long exp) const;

    /**
     * @brief 整数幂运算（BigInt 指数）
     * @param exp 非负 BigInt 指数
     * @return this^exp
     * @throw std::domain_error 指数为负时抛出
     */
    BigInt power(BigInt exp) const;

    /**
     * @brief 整数平方根（向下取整）
     * @return floor(sqrt(*this))
     * @throw std::domain_error 负数时抛出
     */
    BigInt sqrt() const;

    /** @brief 判断是否为奇数 */
    bool is_odd() const;

    /** @brief 判断是否为偶数 */
    bool is_even() const;

    /**
     * @brief 计算末尾零比特数
     * @return 二进制表示中末尾连续 0 的个数
     */
    mp_size_t trailing_zeros() const;

    /**
     * @brief 右移赋值
     * @param shift 移位比特数
     * @return *this
     */
    BigInt& operator>>=(mp_size_t shift);

    /**
     * @brief 左移赋值
     * @param shift 移位比特数
     * @return *this
     */
    BigInt& operator<<=(mp_size_t shift);

    /** @brief 右移运算符 */
    BigInt operator>>(mp_size_t shift) const;

    /** @brief 左移运算符 */
    BigInt operator<<(mp_size_t shift) const;

    /**
     * @brief 计算阶乘
     * @param n 非负整数
     * @return n!
     */
    static BigInt factorial(unsigned int n);

    /**
     * @brief 计算排列数 P(n, r)
     * @param n 总数
     * @param r 选取数
     * @return n! / (n-r)!
     */
    static BigInt nPr(unsigned int n, unsigned int r);

    /**
     * @brief 计算组合数 C(n, r)
     * @param n 总数
     * @param r 选取数
     * @return n! / (r! * (n-r)!)
     */
    static BigInt nCr(unsigned int n, unsigned int r);

    /**
     * @brief 计算多项式系数（多重组合数）
     * @param n 总数（应等于 r 各元素之和）
     * @param r 各组大小
     * @return n! / (r[0]! * r[1]! * ... * r[k]!)
     * @throw std::invalid_argument r 之和不等于 n 时抛出
     */
    static BigInt multinomial(unsigned int n, const std::vector<unsigned int>& r);

    /**
     * @brief 计算最大公约数
     * @param a 第一个大整数
     * @param b 第二个大整数
     * @return gcd(|a|, |b|)
     */
    static BigInt gcd(const BigInt& a, const BigInt& b);

    /**
     * @brief 计算最小公倍数
     * @param a 第一个大整数
     * @param b 第二个大整数
     * @return lcm(|a|, |b|)
     */
    static BigInt lcm(const BigInt& a, const BigInt& b);

    /**
     * @brief 模幂运算
     * @param base 底数
     * @param exp 指数
     * @param mod 模数
     * @return base^exp mod mod
     * @throw std::runtime_error 模数为零时抛出
     */
    static BigInt pow_mod(const BigInt& base, const BigInt& exp, const BigInt& mod);

    /**
     * @brief 素性测试（Miller-Rabin）
     * @return 是素数返回 true
     */
    bool is_prime() const;

    /**
     * @brief 判断是否为完全平方数
     * @return 是完全平方数返回 true
     */
    bool is_perfect_square() const;
};
