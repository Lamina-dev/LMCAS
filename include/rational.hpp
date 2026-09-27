/**
 * @file rational.hpp
 * @brief 精确有理数 Rational，基于 BigInt 分子/分母表示，自动约分。
 */
#pragma once
#include "bigint.hpp"
#include "lmcas_export.hpp"
#include "lmmc/config.h"
#include <cstddef>
#include <cstdint>
#include <string>
namespace LMCAS {

/** @brief 精确有理数类，以 BigInt 分子/分母表示，构造时自动约分 */
class LMCAS_API Rational {
private:
    static constexpr std::size_t max_literal_bytes = 1'048'576;
    static constexpr std::size_t max_decimal_scale = 1'000'000;

    BigInt numerator;
    BigInt denominator;

    static BigInt pow10(std::size_t exponent);

    static std::size_t decimal_scale(std::int64_t exponent);

    static bool power_leq(BigInt base, unsigned long exponent,
                          const BigInt& limit);

    static BigInt integer_nth_root(const BigInt& value, unsigned long degree);

    void simplify();

public:

    /** @brief 默认构造，值为 0 */
    Rational() : numerator(0), denominator(1) {}

    /**
     * @brief 从 BigInt 构造整数有理数
     * @param num 分子（分母为 1）
     */
    Rational(const BigInt& num) : numerator(num), denominator(1) {}

    /**
     * @brief 从分子分母构造有理数
     * @param num 分子
     * @param den 分母（不可为零）
     */
    Rational(const BigInt& num, const BigInt& den);

    /**
     * @brief 从 int 构造整数有理数
     * @param num 整数值
     */
    Rational(int num) : numerator(num), denominator(1) {}

    /**
     * @brief 从 int 分子分母构造有理数
     * @param num 分子
     * @param den 分母
     */
    Rational(int num, int den);

    /**
     * @brief 从字符串构造有理数，支持小数、科学计数法、循环小数
     * @param num 数字字符串
     */
    explicit Rational(const std::string& num);

public:

    /**
     * @brief 从 double 构造精确有理数。
     * @param value 有限的 binary64 浮点数值
     * @return binary64 存储值对应的精确有理数，而非十进制显示文本的值
     * @throws std::invalid_argument 若 value 为 NaN 或无穷大
     */
    static Rational from_double(double value);

    /** @brief 获取分子 */
    BigInt get_numerator() const { return numerator; }

    /** @brief 获取分母 */
    BigInt get_denominator() const { return denominator; }

    /** @brief 判断是否为整数（分母为 1） */
    bool is_integer() const {
        return denominator == BigInt(1);
    }

    /** @brief 判断是否为零 */
    bool is_zero() const {
        return !numerator;
    }

    /**
     * @brief 计算哈希值
     * @return 哈希值
     */
    std::size_t hash() const;

    /**
     * @brief 转换为分数字符串（如 "3/4"）
     * @return 字符串表示
     */
    std::string to_string() const;

    /**
     * @brief 转换为十进制小数字符串（自动检测循环节）
     * @return 小数字符串，循环部分用括号标记
     */
    std::string to_float_string();

    /**
     * @brief 转换为指定精度的十进制小数字符串
     * @param n 小数位数（正数为小数点后位数，负数为整数部分截断位数）
     * @return 定精度小数字符串
     */
    std::string to_float_string(std::int64_t n) const;

    /**
     * @brief 就地向零截断到指定十进制精度，区别于数学下取整。
     * @param n 正数为小数位，零为整数位，负数为十进制整数位
     * @throws std::length_error |n| 超过 1,000,000 时抛出，原值不变
     */
    void floor(const std::int64_t& n);

private:

    void floor_without_sim(const std::int64_t& n);

public:

    /**
     * @brief 牛顿迭代法计算平方根（就地修改）
     * @param n 精度位数
     */
    void sqrt_self(std::int64_t n);

    /**
     * @brief 计算平方根（返回新值）
     * @param n 精度位数
     * @return 平方根的有理近似
     */
    Rational sqrt(std::int64_t n) const;

    /**
     * @brief 牛顿迭代法计算 n 次方根（就地修改）
     * @param radical 根次数
     * @param n 精度位数
     */
    void radicand_self(const BigInt& radical, std::int64_t n);

    /**
     * @brief 计算 n 次方根（返回新值）
     * @param radical 根次数
     * @param n 精度位数
     * @return n 次方根的有理近似
     */
    Rational radicand(const BigInt& radical, std::int64_t n);

    /**
     * @brief 转换为 BigInt（仅当为整数时有效）
     * @return 对应的 BigInt
     * @throw std::runtime_error 非整数时抛出
     */
    BigInt to_bigint() const;

    /**
     * @brief 转换为浮点数。
     *
     * 精确值可按 binary64 round-to-nearest 转为有限值时返回有限结果，
     * 包括正负最大有限端点及其半 ULP 内邻域；
     * 超出有限舍入范围或非零值下溢为零时抛出范围异常。
     * @return 对应的 double 值
     */
    lmmc_real_t to_double() const;

    /** @brief 有理数加法 */
    Rational operator+(const Rational& other) const;

    /** @brief 有理数减法 */
    Rational operator-(const Rational& other) const;

    /** @brief 有理数乘法 */
    Rational operator*(const Rational& other) const;

    /**
     * @brief 有理数除法
     * @param other 除数
     * @return 商
     * @throw std::runtime_error 除数为零时抛出
     */
    Rational operator/(const Rational& other) const;

    /**
     * @brief 有理数整数幂
     * @param exponent 指数（可为负）
     * @return this^exponent
     * @throw std::runtime_error 零的负幂或 0^0 时抛出
     */
    Rational power(const BigInt& exponent) const;

    /**
     * @brief 有理数的有理数幂（就地修改）
     * @param exponent 有理数指数
     * @param n 精度位数
     */
    void power_self(const Rational& exponent, std::size_t n);

    /**
     * @brief 有理数的有理数幂（返回新值）
     * @param exponent 有理数指数
     * @param n 精度位数
     * @return 结果的有理近似
     */
    Rational power(const Rational& exponent, std::size_t n);

    bool operator==(const Rational& other) const;

    bool operator!=(const Rational& other) const;

    bool operator<(const Rational& other) const;

    bool operator<=(const Rational& other) const;

    bool operator>(const Rational& other) const;

    bool operator>=(const Rational& other) const;

    /**
     * @brief 取倒数
     * @return 1 / *this
     * @throw std::runtime_error 零的倒数时抛出
     */
    Rational reciprocal() const;

    /**
     * @brief 取绝对值
     * @return |*this|
     */
    Rational abs() const;

    /** @brief 一元负号运算符 */
    Rational operator-() const;
};

} // namespace LMCAS
