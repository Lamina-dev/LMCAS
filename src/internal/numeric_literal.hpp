#pragma once

#include "result.hpp"
#include <charconv>
#include <cmath>
#include <string_view>
#include <system_error>

namespace LMCAS::detail {

inline bool consume_decimal_digits(std::string_view token, std::size_t& pos) noexcept {
    const auto start = pos;
    while (pos < token.size() && token[pos] >= '0' && token[pos] <= '9') {
        ++pos;
    }
    return pos != start;
}

inline bool consume_decimal_mantissa(std::string_view token, std::size_t& pos) noexcept {
    bool has_digits = consume_decimal_digits(token, pos);
    if (pos < token.size() && token[pos] == '.') {
        ++pos;
        has_digits = consume_decimal_digits(token, pos) || has_digits;
    }
    return has_digits;
}

inline void consume_decimal_sign(std::string_view token, std::size_t& pos) noexcept {
    if (pos < token.size() && (token[pos] == '+' || token[pos] == '-')) {
        ++pos;
    }
}

inline bool consume_decimal_exponent(std::string_view token, std::size_t& pos) noexcept {
    if (pos < token.size() && (token[pos] == 'e' || token[pos] == 'E')) {
        ++pos;
        consume_decimal_sign(token, pos);
        return consume_decimal_digits(token, pos);
    }
    return true;
}

inline bool is_signed_decimal_token(std::string_view token) noexcept {
    std::size_t pos = 0;
    consume_decimal_sign(token, pos);
    if (!consume_decimal_mantissa(token, pos)) { return false; }
    if (!consume_decimal_exponent(token, pos)) { return false; }
    return pos == token.size();
}

inline Result<double> parse_finite_decimal(std::string_view token) {
    const auto invalid = [] {
        return Result<double>::failure(CasErrc::ParseError,
            "expected a representable finite decimal", "LMCAS.numeric_literal");
    };
    if (!is_signed_decimal_token(token)) { return invalid(); }
    if (token.front() == '+') { token.remove_prefix(1); }
    double value = 0;
    const auto parsed = std::from_chars(token.data(), token.data() + token.size(),
                                        value, std::chars_format::general);
    if (parsed.ec != std::errc{} || parsed.ptr != token.data() + token.size() ||
        !std::isfinite(value)) { return invalid(); }
    if (value == 0) {
        for (char c : token) {
            if (c == 'e' || c == 'E') { break; }
            if (c >= '1' && c <= '9') { return invalid(); }
        }
        if (token.front() == '-') { value = -0.0; }
    }
    return value;
}

}
