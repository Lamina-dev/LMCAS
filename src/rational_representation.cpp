#include "rational.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace LMCAS {
namespace {

std::string format_scaled_decimal(BigInt scaled, std::size_t digits) {
    const bool negative = scaled < BigInt(0);
    if (negative) {
        scaled = scaled * BigInt(-1);
    }
    std::string text = scaled.to_string();
    if (text.size() <= digits) {
        text.insert(0, digits - text.size() + 1, '0');
    }
    text.insert(text.size() - digits, 1, '.');
    while (text.back() == '0') {
        text.pop_back();
    }
    if (text.back() == '.') {
        text.pop_back();
    }
    if (negative) {
        text.insert(text.begin(), '-');
    }
    return text;
}

}

BigInt Rational::pow10(std::size_t exponent) {
    if (exponent > max_decimal_scale) {
        throw std::length_error("Rational decimal scale exceeds safety limit");
    }
    if (exponent == 0) return BigInt(1);
    return BigInt(10).power(static_cast<unsigned long>(exponent));
}

std::size_t Rational::decimal_scale(std::int64_t exponent) {
    constexpr auto limit = static_cast<std::int64_t>(max_decimal_scale);
    if (exponent < -limit || exponent > limit) {
        throw std::length_error("Rational decimal scale exceeds safety limit");
    }
    return static_cast<std::size_t>(exponent < 0 ? -exponent : exponent);
}

Rational Rational::from_double(double value) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("Cannot construct Rational from NaN or infinity");
    }
    if (value == 0.0) {
        return Rational();
    }
    int exponent = 0;
    const double fraction = std::frexp(value, &exponent);
    const auto mantissa = static_cast<std::int64_t>(std::ldexp(fraction, 53));
    const int binary_exponent = exponent - 53;

    BigInt num(mantissa);
    BigInt den(1);
    if (binary_exponent >= 0) {
        num <<= static_cast<mp_size_t>(binary_exponent);
    } else {
        den <<= static_cast<mp_size_t>(-binary_exponent);
    }
    return Rational(num, den);
}

std::size_t Rational::hash() const {
    std::size_t seed = numerator.hash();
    std::size_t d_hash = denominator.hash();
    seed ^= d_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

std::string Rational::to_string() const {
    if (is_integer()) {
        return numerator.to_string();
    }
    return numerator.to_string() + "/" + denominator.to_string();
}

std::string Rational::to_float_string() {
    if (numerator % denominator == BigInt(0)) return (numerator / denominator).to_string();
    std::string re;
    const BigInt ten(10);
    const BigInt zero(0);
    if ((numerator < zero) ^ (denominator < zero)) re += '-';
    BigInt num = numerator.abs(), den = denominator.abs();
    re += (num / den).to_string();
    num %= den;
    re += '.';

    std::map<BigInt, std::size_t> num_index;

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
        re.insert(num_index[num], 1, '(');
        re += ")...";
    }
    return re;

}

std::string Rational::to_float_string(std::int64_t n) const{
    if (is_zero()) {
        return "0";
    }
    if (n == 0) {
        return (numerator / denominator).to_string();
    }
    const std::size_t digits = decimal_scale(n);
    const BigInt pow10n = pow10(digits);
    if (n > 0) {
        return format_scaled_decimal(numerator * pow10n / denominator, digits);
    }
    std::string re = (numerator / (denominator * pow10n)).to_string();
    if (re.size() == 1 && re[0] == '0') {
        return re;
    }
    re.append(digits, '0');
    return re;
}

void Rational::floor(const std::int64_t& n) {
    floor_without_sim(n);
    simplify();
}

void Rational::floor_without_sim(const std::int64_t& n) {
    if (n == 0) {
        numerator /= denominator;
        denominator = BigInt(1);
        return;
    }
    const BigInt pow10n = pow10(decimal_scale(n));
    if (n > 0) {
        numerator *= pow10n;
        numerator /= denominator;
        denominator = pow10n;
    } else {
        denominator *= pow10n;
        numerator /= denominator;
        numerator *= pow10n;
        denominator = BigInt(1);
    }
}

BigInt Rational::to_bigint() const {
    if (!is_integer()) {
        throw std::runtime_error("Cannot convert non-integer fraction to BigInt");
    }
    return numerator;
}

lmmc_real_t Rational::to_double() const {
    if (is_zero()) return 0.0;
    const std::string num_text = numerator.abs().to_string();
    const std::string den_text = denominator.to_string();
    const std::size_t num_digits = std::min<std::size_t>(18, num_text.size());
    const std::size_t den_digits = std::min<std::size_t>(18, den_text.size());
    const long double num_head = std::stold(num_text.substr(0, num_digits));
    const long double den_head = std::stold(den_text.substr(0, den_digits));
    const std::size_t num_tail = num_text.size() - num_digits;
    const std::size_t den_tail = den_text.size() - den_digits;
    const bool negative_shift = num_tail < den_tail;
    const std::size_t shift = negative_shift ? den_tail - num_tail : num_tail - den_tail;
    if (shift > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
        if (negative_shift) {
            throw std::underflow_error("Rational underflow during double conversion");
        }
        throw std::overflow_error("Rational cannot be represented as a finite double");
    }
    const std::int64_t decimal_shift = negative_shift ?
        -static_cast<std::int64_t>(shift) : static_cast<std::int64_t>(shift);
    long double value = (num_head / den_head) *
                        std::pow(10.0L, static_cast<long double>(decimal_shift));
    if (numerator.is_negative()) value = -value;
    const lmmc_real_t converted = static_cast<lmmc_real_t>(value);
    if (!std::isfinite(converted)) {
        /**
         * @brief 用精确整数运算判定 binary64 的有限值舍入边界。
         * 十进制前缀近似可能越过上界，而精确有理数仍舍入为 DBL_MAX；
         * 绝对值小于 DBL_MAX + 0.5 ULP 时仍舍入为有限值。
         */
        const BigInt max_finite_integer =
            (BigInt(1) << 1024) - (BigInt(1) << 971);
        const BigInt finite_rounding_limit =
            max_finite_integer + (BigInt(1) << 970);
        if (numerator.abs() <
            denominator * finite_rounding_limit) {
            return std::copysign(
                std::numeric_limits<lmmc_real_t>::max(),
                numerator.is_negative() ? -1.0 : 1.0);
        }
        throw std::overflow_error(
            "Rational cannot be represented as a finite double");
    }
    if (converted == 0.0) {
        throw std::underflow_error("Rational underflow during double conversion");
    }
    return converted;
}

}
