/**
 * @file rational.hpp
 * @brief 精确有理数 Rational，基于 BigInt 分子/分母表示，自动约分。
 */
#pragma once
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "bigint.hpp"
#include <iomanip>
#include <cctype>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <map>
#include <stdint.h>

namespace LMCAS {

/** @brief 精确有理数类，以 BigInt 分子/分母表示，构造时自动约分 */
class Rational {
private:
    static constexpr std::size_t max_literal_bytes = 1'048'576;
    static constexpr std::size_t max_decimal_scale = 1'000'000;

    BigInt numerator;
    BigInt denominator;

    static BigInt pow10(std::size_t exponent) {
        if (exponent > max_decimal_scale) {
            throw std::length_error("Rational decimal scale exceeds safety limit");
        }
        if (exponent == 0) return BigInt(1);
        return BigInt(10).power(static_cast<unsigned long>(exponent));
    }

    static bool power_leq(BigInt base, unsigned long exponent,
                          const BigInt& limit) {
        BigInt result(1);
        while (exponent != 0) {
            if ((exponent & 1UL) != 0) {
                result *= base;
                if (result > limit) return false;
            }
            exponent >>= 1UL;
            if (exponent != 0) {
                base *= base;
                if (base > limit) base = limit + BigInt(1);
            }
        }
        return result <= limit;
    }

    static BigInt integer_nth_root(const BigInt& value, unsigned long degree) {
        if (degree == 0) throw std::domain_error("Zeroth root is undefined");
        if (value < BigInt(0)) throw std::domain_error("Root input must be non-negative");
        if (value <= BigInt(1) || degree == 1) return value;

        BigInt low(0);
        BigInt high(1);
        while (power_leq(high, degree, value)) high *= BigInt(2);

        while (high - low > BigInt(1)) {
            BigInt mid = (low + high) / BigInt(2);
            if (power_leq(mid, degree, value)) low = std::move(mid);
            else high = std::move(mid);
        }
        return low;
    }

    void simplify() {
        if (!denominator) {
            throw std::runtime_error("Denominator cannot be zero");
        }

        if (denominator.IsNegative()) {
            numerator = -numerator;
            denominator = denominator.Abs();
        }

        BigInt g = BigInt::gcd(numerator, denominator);
        if (g != 0 && g != 1) {
            numerator = numerator / g;
            denominator = denominator / g;
        }
    }

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
    Rational(const BigInt& num, const BigInt& den) : numerator(num), denominator(den) {
        simplify();
    }

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
    Rational(int num, int den) : numerator(num), denominator(den) {
        simplify();
    }

    /**
     * @brief 从字符串构造有理数，支持小数、科学计数法、循环小数
     * @param num 数字字符串
     */
    explicit Rational(const std::string& num) : Rational() {
        if (num.empty()) throw std::invalid_argument("Rational literal is empty");
        if (num.size() > max_literal_bytes) {
            throw std::length_error("Rational literal exceeds safety limit");
        }

        std::size_t pos = 0;
        bool negative = false;
        if (num[pos] == '+' || num[pos] == '-') {
            negative = num[pos] == '-';
            if (++pos == num.size()) {
                throw std::invalid_argument("Rational literal has no digits");
            }
        }

        auto read_digits = [&](std::string& output) {
            while (pos < num.size() &&
                   std::isdigit(static_cast<unsigned char>(num[pos]))) {
                output.push_back(num[pos++]);
            }
        };

        std::string integer_digits;
        std::string fractional_digits;
        std::string repeating_digits;
        read_digits(integer_digits);

        bool has_decimal_point = false;
        if (pos < num.size() && num[pos] == '.') {
            has_decimal_point = true;
            ++pos;
            read_digits(fractional_digits);
        }

        if (pos < num.size() && num[pos] == '(') {
            if (!has_decimal_point) {
                throw std::invalid_argument("Repeating decimal section requires a decimal point");
            }
            ++pos;
            read_digits(repeating_digits);
            if (repeating_digits.empty() || pos >= num.size() || num[pos] != ')') {
                throw std::invalid_argument("Invalid repeating decimal section");
            }
            ++pos;
        }

        if (integer_digits.empty() && fractional_digits.empty()) {
            throw std::invalid_argument("Rational literal has no digits");
        }

        long long exponent = 0;
        if (pos < num.size() && (num[pos] == 'e' || num[pos] == 'E')) {
            ++pos;
            bool exponent_negative = false;
            if (pos < num.size() && (num[pos] == '+' || num[pos] == '-')) {
                exponent_negative = num[pos] == '-';
                ++pos;
            }
            if (pos == num.size() ||
                !std::isdigit(static_cast<unsigned char>(num[pos]))) {
                throw std::invalid_argument("Rational exponent has no digits");
            }
            while (pos < num.size() &&
                   std::isdigit(static_cast<unsigned char>(num[pos]))) {
                const int digit = num[pos++] - '0';
                if (exponent > static_cast<long long>(max_decimal_scale / 10) ||
                    exponent * 10 + digit > static_cast<long long>(max_decimal_scale)) {
                    throw std::length_error("Rational exponent exceeds safety limit");
                }
                exponent = exponent * 10 + digit;
            }
            if (exponent_negative) exponent = -exponent;
        }

        if (pos != num.size()) {
            throw std::invalid_argument("Invalid character in Rational literal");
        }

        if (integer_digits.empty()) integer_digits = "0";
        const std::string non_repeating = integer_digits + fractional_digits;
        numerator = BigInt(non_repeating);
        denominator = pow10(fractional_digits.size());

        if (!repeating_digits.empty()) {
            const BigInt repeat_scale = pow10(repeating_digits.size());
            const BigInt repeat_factor = repeat_scale - BigInt(1);
            numerator = numerator * repeat_factor + BigInt(repeating_digits);
            denominator *= repeat_factor;
        }

        if (exponent > 0) {
            numerator *= pow10(static_cast<std::size_t>(exponent));
        } else if (exponent < 0) {
            denominator *= pow10(static_cast<std::size_t>(-exponent));
        }
        if (negative) numerator = -numerator;
        simplify();
    }

public:

    /**
     * @brief 从 double 构造有理数
     * @param value 浮点数值
     * @return 对应的精确有理数
     */
    static Rational from_double(double value) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Cannot construct Rational from NaN or infinity");
        }
        if (value == 0.0) {
            return Rational();
        }
        if (std::floor(value) == value) {
            // value 是整数级浮点数；std::to_string(double) 会输出 "123.000000"，
            // 直接喂给 BigInt 会因含小数点被拒绝。先把整数部分提取成字符串。
            std::ostringstream iss;
            iss << std::fixed << std::setprecision(0) << value;
            return Rational(BigInt(iss.str()));
        }
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(15) << value;
        std::string str = oss.str();

        size_t e_pos = str.find('e');
        std::string base = str.substr(0, e_pos);
        bool negative = false;
        if (base[0] == '-') {
            negative = true;
            base.erase(0, 1);
        }
        base.erase(1, 1);
        while (base.back() == '0') {
            base.pop_back();
        }
        int exponent = std::stoi(str.substr(e_pos + 1));

        size_t decimal_places = std::max(0, static_cast<int>(base.length()) - exponent - 1);

        BigInt num(base);
        if (negative) num = -num;
        BigInt den("1" + std::string(decimal_places, '0'));

        Rational r(num, den);
        r.simplify();
        return r;
    }

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
    std::size_t hash() const {
        std::size_t seed = numerator.hash();
        std::size_t d_hash = denominator.hash();
        seed ^= d_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }

    /**
     * @brief 转换为分数字符串（如 "3/4"）
     * @return 字符串表示
     */
    std::string to_string() const {
        if (is_integer()) {
            return numerator.ToString();
        }
        return numerator.ToString() + "/" + denominator.ToString();
    }

    /**
     * @brief 转换为十进制小数字符串（自动检测循环节）
     * @return 小数字符串，循环部分用括号标记
     */
    inline std::string to_float_string() {
        if (numerator % denominator == BigInt(0)) return (numerator / denominator).ToString();
        std::string re;
        const BigInt ten(10);
        const BigInt zero(0);
        if ((numerator < zero) ^ (denominator < zero)) re += '-';
        BigInt num = numerator.Abs(), den = denominator.Abs();
        re += (num / den).ToString();
        num %= den;
        re += '.';

        std::map<BigInt, size_t> num_index;

        while (num != zero && num_index.find(num) == num_index.end()) {
            num_index[num] = re.size();
            re += '0';
            num *= ten;
            while (num >= den) {
                num -= den;
                re.back()++;
            }
        }
        if (num_index.find(num) != num_index.end()) {
            re.insert(re.begin() + num_index[num], '(');
            re += ")...";
        }
        return re;

    }

    /**
     * @brief 转换为指定精度的十进制小数字符串
     * @param n 小数位数（正数为小数点后位数，负数为整数部分截断位数）
     * @return 定精度小数字符串
     */
    inline std::string to_float_string(int64_t n) const{
        if (is_zero()) return "0";
        if (n == 0) return (numerator / denominator).ToString();
        if (n > 0) {
            BigInt pow10n(10);
            pow10n = pow10n.power(BigInt(static_cast<long long>(n)));
            BigInt scaled = numerator * pow10n / denominator;
            bool negative = scaled < BigInt(0);
            if (negative) scaled = scaled * BigInt(-1);
            std::string re = scaled.ToString();
            if ((int64_t)re.size() <= n) re.insert(re.begin(), n - re.size() + 1, '0');
            re.insert(re.end() - n,'.');
            while (re.back() == '0') re.pop_back();
            if (re.back() == '.') re.pop_back();
            if (negative) re.insert(re.begin(), '-');
            return re;
        } else {
            BigInt pow10n(10);
            pow10n = pow10n.power(BigInt(static_cast<long long>(n)).Abs());
            std::string re = (numerator / (denominator * pow10n)).ToString();
            if (re.size() == 1 && re[0] == '0') return re;
            for (; n < 0; n++) re.push_back('0');
            return re;
        }
    }

    /**
     * @brief 截断到指定精度（就地修改）
     * @param n 精度位数
     */
    inline void floor(const int64_t& n) {
        floor_without_sim(n);
        simplify();
    }

private:

    inline void floor_without_sim(const int64_t& n) {
        if (n == 0) return;
        if (n > 0) {
            BigInt pow10n(10);
            pow10n = pow10n.power(BigInt(static_cast<long long>(n)));
            numerator *= pow10n;
            numerator /= denominator;
            denominator = pow10n;
        } else {
            BigInt pow10n(10);
            pow10n = pow10n.power(BigInt(static_cast<long long>(n)).Abs());
            denominator *= pow10n;
            numerator /= denominator;
            numerator *= pow10n;
            denominator = BigInt(1);
        }
    }

public:

    /**
     * @brief 牛顿迭代法计算平方根（就地修改）
     * @param n 精度位数
     */
    inline void sqrt_self(int64_t n) {
        if (n < 0) throw std::invalid_argument("Square-root precision must be non-negative");
        if (numerator < BigInt(0)) throw std::domain_error("Square root of a negative Rational");
        if (numerator == BigInt(0)) {
            denominator = BigInt(1);
            return;
        }
        if (static_cast<std::size_t>(n) > max_decimal_scale / 2) {
            throw std::length_error("Square-root precision exceeds safety limit");
        }

        BigInt numerator_root = numerator.sqrt();
        BigInt denominator_root = denominator.sqrt();
        if (numerator_root * numerator_root == numerator &&
            denominator_root * denominator_root == denominator) {
            numerator = std::move(numerator_root);
            denominator = std::move(denominator_root);
            simplify();
            return;
        }

        const BigInt scale = pow10(static_cast<std::size_t>(n));
        const BigInt scaled = numerator * scale * scale / denominator;
        numerator = scaled.sqrt();
        denominator = scale;
        simplify();
    }

    /**
     * @brief 计算平方根（返回新值）
     * @param n 精度位数
     * @return 平方根的有理近似
     */
    inline Rational sqrt(int64_t n) const{
        Rational re = *this;
        re.sqrt_self(n);
        return re;
    }

    /**
     * @brief 牛顿迭代法计算 n 次方根（就地修改）
     * @param radical 根次数
     * @param n 精度位数
     */
    inline void radicand_self(const BigInt& radical, int64_t n) {
        if (radical <= BigInt(0)) throw std::domain_error("Root degree must be positive");
        if (radical > BigInt(1024)) throw std::length_error("Root degree exceeds safety limit");
        if (n < 0) throw std::invalid_argument("Root precision must be non-negative");

        const unsigned long degree = static_cast<unsigned long>(radical.to_int());
        const bool negative = numerator < BigInt(0);
        if (negative && (degree % 2UL) == 0) {
            throw std::domain_error("Even root of a negative Rational");
        }
        if (static_cast<std::size_t>(n) > max_decimal_scale / degree) {
            throw std::length_error("Root precision exceeds safety limit");
        }

        const BigInt scale = pow10(static_cast<std::size_t>(n));
        const BigInt scaled = numerator.Abs() * scale.power(degree) / denominator;
        numerator = integer_nth_root(scaled, degree);
        if (negative) numerator = -numerator;
        denominator = scale;
        simplify();
    }

    /**
     * @brief 计算 n 次方根（返回新值）
     * @param radical 根次数
     * @param n 精度位数
     * @return n 次方根的有理近似
     */
    inline Rational radicand(const BigInt& radical, int64_t n) {
        Rational re = *this;
        re.radicand_self(radical, n);
        return re;
    }

    /**
     * @brief 转换为 BigInt（仅当为整数时有效）
     * @return 对应的 BigInt
     * @throw std::runtime_error 非整数时抛出
     */
    BigInt to_BigInt() const {
        if (!is_integer()) {
            throw std::runtime_error("Cannot convert non-integer fraction to BigInt");
        }
        return numerator;
    }

    /**
     * @brief 转换为浮点数。
     *
     * 精确值可按 binary64 round-to-nearest 转为有限值时返回有限结果，
     * 包括正负最大有限端点及其半 ULP 内邻域；
     * 超出有限舍入范围或非零值下溢为零时抛出范围异常。
     * @return 对应的 double 值
     */
    lmmc_real_t to_double() const {
        if (is_zero()) return 0.0;
        const std::string num_text = numerator.Abs().ToString();
        const std::string den_text = denominator.ToString();
        const std::size_t num_digits = std::min<std::size_t>(18, num_text.size());
        const std::size_t den_digits = std::min<std::size_t>(18, den_text.size());
        const long double num_head = std::stold(num_text.substr(0, num_digits));
        const long double den_head = std::stold(den_text.substr(0, den_digits));
        const long long decimal_shift =
            static_cast<long long>(num_text.size()) - static_cast<long long>(num_digits) -
            static_cast<long long>(den_text.size()) + static_cast<long long>(den_digits);
        long double value = (num_head / den_head) *
                            std::pow(10.0L, static_cast<long double>(decimal_shift));
        if (numerator.IsNegative()) value = -value;
        const lmmc_real_t converted = static_cast<lmmc_real_t>(value);
        if (!std::isfinite(converted)) {
            /*
             * A decimal-head approximation can round slightly above the
             * binary64 endpoint even when the exact rational still rounds to
             * DBL_MAX. Resolve that boundary case with exact integer
             * arithmetic. Values below max + 0.5 ULP round to a finite value.
             */
            const BigInt max_finite_integer =
                (BigInt(1) << 1024) - (BigInt(1) << 971);
            const BigInt finite_rounding_limit =
                max_finite_integer + (BigInt(1) << 970);
            if (numerator.Abs() <
                denominator * finite_rounding_limit) {
                return std::copysign(
                    std::numeric_limits<lmmc_real_t>::max(),
                    numerator.IsNegative() ? -1.0 : 1.0);
            }
            throw std::overflow_error(
                "Rational cannot be represented as a finite double");
        }
        if (converted == 0.0) {
            throw std::underflow_error("Rational underflow during double conversion");
        }
        return converted;
    }

    /** @brief 有理数加法 */
    Rational operator+(const Rational& other) const {
        BigInt new_num = numerator * other.denominator + other.numerator * denominator;
        BigInt new_den = denominator * other.denominator;
        return Rational(new_num, new_den);
    }

    /** @brief 有理数减法 */
    Rational operator-(const Rational& other) const {
        BigInt new_num = numerator * other.denominator - other.numerator * denominator;
        BigInt new_den = denominator * other.denominator;
        return Rational(new_num, new_den);
    }

    /** @brief 有理数乘法 */
    Rational operator*(const Rational& other) const {
        BigInt new_num = numerator * other.numerator;
        BigInt new_den = denominator * other.denominator;
        return Rational(new_num, new_den);
    }

    /**
     * @brief 有理数除法
     * @param other 除数
     * @return 商
     * @throw std::runtime_error 除数为零时抛出
     */
    Rational operator/(const Rational& other) const {
        if (other.is_zero()) {
            throw std::runtime_error("Division by zero");
        }
        BigInt new_num = numerator * other.denominator;
        BigInt new_den = denominator * other.numerator;
        return Rational(new_num, new_den);
    }

    /**
     * @brief 有理数整数幂
     * @param exponent 指数（可为负）
     * @return this^exponent
     * @throw std::runtime_error 零的负幂或 0^0 时抛出
     */
    Rational power(const BigInt& exponent) const {
        if (exponent < BigInt(0)) {

            if (is_zero()) {
                throw std::runtime_error("Cannot raise zero to negative power");
            }
            BigInt pos_exp = -exponent;
            return Rational(denominator.power(pos_exp), numerator.power(pos_exp));
        }

        if (!exponent) {
            if (is_zero()) {
                throw std::runtime_error("0^0 is undefined");
            }
            return Rational(1);
        }

        return Rational(numerator.power(exponent), denominator.power(exponent));
    }

    /**
     * @brief 有理数的有理数幂（就地修改）
     * @param exponent 有理数指数
     * @param n 精度位数
     */
    inline void power_self(const Rational& exponent,size_t n) {
        *this = power(exponent.numerator);
        radicand_self(exponent.denominator,n);
    }

    /**
     * @brief 有理数的有理数幂（返回新值）
     * @param exponent 有理数指数
     * @param n 精度位数
     * @return 结果的有理近似
     */
    inline Rational power(const Rational& exponent, size_t n) {
        Rational re = *this;
        re.power_self(exponent, n);
        return re;
    }

    bool operator==(const Rational& other) const {
        return numerator * other.denominator == other.numerator * denominator;
    }

    bool operator!=(const Rational& other) const {
        return !(*this == other);
    }

    bool operator<(const Rational& other) const {
        BigInt left = numerator * other.denominator;
        BigInt right = other.numerator * denominator;
        return left < right;
    }

    bool operator<=(const Rational& other) const {
        return *this < other || *this == other;
    }

    bool operator>(const Rational& other) const {
        return !(*this <= other);
    }

    bool operator>=(const Rational& other) const {
        return !(*this < other);
    }

    /**
     * @brief 取倒数
     * @return 1 / *this
     * @throw std::runtime_error 零的倒数时抛出
     */
    Rational reciprocal() const {
        if (is_zero()) {
            throw std::runtime_error("Cannot take reciprocal of zero");
        }
        return Rational(denominator, numerator);
    }

    /**
     * @brief 取绝对值
     * @return |*this|
     */
    Rational abs() const {
        Rational result = *this;
        result.numerator = result.numerator.Abs();
        return result;
    }

    /** @brief 一元负号运算符 */
    Rational operator-() const {
        Rational result = *this;
        result.numerator = -result.numerator;
        return result;
    }
};

} // namespace LMCAS
