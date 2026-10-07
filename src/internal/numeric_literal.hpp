#pragma once

#include "result.hpp"
#include <cmath>
#include <string_view>
#ifdef __APPLE__
#include <cstdlib>
#include <locale.h>
#include <string>
#else
#include <charconv>
#include <system_error>
#endif

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
#ifdef __APPLE__
    // Xcode 16 uses strtod_l with a fixed C locale for decimal conversion.
    static const locale_t c_locale = newlocale(LC_NUMERIC_MASK, "C", nullptr);
    if (!c_locale) { return invalid(); }
    std::string text(token);
    char* end = nullptr;
    value = strtod_l(text.c_str(), &end, c_locale);
    if (end != text.c_str() + text.size() || !std::isfinite(value)) {
        return invalid();
    }
#else
    const auto parsed = std::from_chars(token.data(), token.data() + token.size(),
                                        value, std::chars_format::general);
    if (parsed.ec != std::errc{} || parsed.ptr != token.data() + token.size() ||
        !std::isfinite(value)) { return invalid(); }
#endif
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
