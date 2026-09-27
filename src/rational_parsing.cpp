#include "rational.hpp"
#include <cctype>
#include <stdexcept>

namespace LMCAS {
namespace {

class RationalLiteral {
public:
    std::string integer_digits;
    std::string fractional_digits;
    std::string repeating_digits;
    bool negative = false;
    std::int64_t exponent = 0;

    RationalLiteral(const std::string& text, std::size_t max_bytes,
                    std::size_t max_scale) : text_(text) {
        if (text.empty()) throw std::invalid_argument("Rational literal is empty");
        if (text.size() > max_bytes) {
            throw std::length_error("Rational literal exceeds safety limit");
        }
        negative = read_sign();
        if (pos_ == text.size()) {
            throw std::invalid_argument("Rational literal has no digits");
        }
        read_digits(integer_digits);
        const bool has_decimal_point = pos_ < text.size() && text[pos_] == '.';
        if (has_decimal_point) {
            ++pos_;
            read_digits(fractional_digits);
        }
        read_repeating_digits(has_decimal_point);
        if (integer_digits.empty() && fractional_digits.empty()) {
            throw std::invalid_argument("Rational literal has no digits");
        }
        exponent = read_exponent(max_scale);
        if (pos_ != text.size()) {
            throw std::invalid_argument("Invalid character in Rational literal");
        }
        if (integer_digits.empty()) integer_digits = "0";
    }

private:
    const std::string& text_;
    std::size_t pos_ = 0;

    bool read_sign() {
        if (pos_ == text_.size() || (text_[pos_] != '+' && text_[pos_] != '-')) {
            return false;
        }
        return text_[pos_++] == '-';
    }

    void read_digits(std::string& output) {
        while (pos_ < text_.size() &&
               std::isdigit(static_cast<unsigned char>(text_[pos_]))) {
            output.push_back(text_[pos_++]);
        }
    }

    void read_repeating_digits(bool has_decimal_point) {
        if (pos_ == text_.size() || text_[pos_] != '(') return;
        if (!has_decimal_point) {
            throw std::invalid_argument("Repeating decimal section requires a decimal point");
        }
        ++pos_;
        read_digits(repeating_digits);
        if (repeating_digits.empty() || pos_ >= text_.size() || text_[pos_] != ')') {
            throw std::invalid_argument("Invalid repeating decimal section");
        }
        ++pos_;
    }

    std::int64_t read_exponent(std::size_t max_scale) {
        if (pos_ == text_.size() || (text_[pos_] != 'e' && text_[pos_] != 'E')) {
            return 0;
        }
        ++pos_;
        const bool exponent_negative = read_sign();
        if (pos_ == text_.size() ||
            !std::isdigit(static_cast<unsigned char>(text_[pos_]))) {
            throw std::invalid_argument("Rational exponent has no digits");
        }
        const std::int64_t value = read_exponent_magnitude(max_scale);
        return exponent_negative ? -value : value;
    }

    std::int64_t read_exponent_magnitude(std::size_t max_scale) {
        std::int64_t value = 0;
        while (pos_ < text_.size() &&
               std::isdigit(static_cast<unsigned char>(text_[pos_]))) {
            const int digit = text_[pos_++] - '0';
            if (value > static_cast<std::int64_t>(max_scale / 10) ||
                value * 10 + digit > static_cast<std::int64_t>(max_scale)) {
                throw std::length_error("Rational exponent exceeds safety limit");
            }
            value = value * 10 + digit;
        }
        return value;
    }
};

}

Rational::Rational(const std::string& num) : Rational() {
    const RationalLiteral literal(num, max_literal_bytes, max_decimal_scale);
    const std::string non_repeating = literal.integer_digits + literal.fractional_digits;
    numerator = BigInt(non_repeating);
    denominator = pow10(literal.fractional_digits.size());

    if (!literal.repeating_digits.empty()) {
        const BigInt repeat_scale = pow10(literal.repeating_digits.size());
        const BigInt repeat_factor = repeat_scale - BigInt(1);
        numerator = numerator * repeat_factor + BigInt(literal.repeating_digits);
        denominator *= repeat_factor;
    }

    if (literal.exponent > 0) {
        numerator *= pow10(static_cast<std::size_t>(literal.exponent));
    } else if (literal.exponent < 0) {
        denominator *= pow10(static_cast<std::size_t>(-literal.exponent));
    }
    if (literal.negative) numerator = -numerator;
    simplify();
}

}
