#include "internal/assumption_parser.hpp"
#include "expr.hpp"
#include "internal/numeric_literal.hpp"
#include <cctype>
#include <vector>
#include <cmath>

namespace LMCAS::assumption_detail {

std::string_view trim_ascii_whitespace(std::string_view text) noexcept {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos) return {};
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

bool is_serialized_numeric_atom(std::string_view expression) noexcept {
    bool saw_digit = false;
    for (std::size_t i = 0; i < expression.size(); ++i) {
        const char c = expression[i];
        if (std::isdigit(static_cast<unsigned char>(c))) {
            saw_digit = true;
            continue;
        }
        if (c == '.' || c == '/' || c == 'e' || c == 'E') {
            continue;
        }
        if ((c == '+' || c == '-') &&
            (i == 0 || expression[i - 1] == 'e' ||
             expression[i - 1] == 'E')) {
            continue;
        }
        return false;
    }
    return saw_digit;
}

char matching_closer(char opener) noexcept {
    switch (opener) {
    case '(': return ')';
    case '[': return ']';
    case '{': return '}';
    default: return '\0';
    }
}

bool is_closing_delimiter(char value) noexcept {
    return value == ')' || value == ']' || value == '}';
}

void require_deserialization_update(const Result<void>& result,
                                    int line,
                                    const std::string& keyword) {
    if (result) {
        return;
    }
    if (result.error().code == CasErrc::InvalidArgument) {
        throw std::invalid_argument(
            "Line " + std::to_string(line) + ": " + keyword + ": " +
            result.error().message);
    }
    throw result.error();
}
void require_line_end(std::istringstream& input,
                      int line_num,
                      const std::string& keyword) {
    std::string trailing;
    if (input >> trailing) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": " + keyword +
            " has unexpected trailing field '" + trailing + "'");
    }
}

std::shared_ptr<const SymbolicNode> parse_finite_number_node(
    const std::string& token,
    int line_num) {
    const auto invalid = [&]() {
        return std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": invalid finite number '" + token + "'");
    };
    if (token.empty()) {
        throw invalid();
    }
    if (token.compare(0, 7, "approx(") == 0 && token.back() == ')') {
        std::string_view decimal(token.data() + 7, token.size() - 8);
        decimal = trim_ascii_whitespace(decimal);
        auto value = LMCAS::detail::parse_finite_decimal(decimal);
        if (!value) { throw invalid(); }
        return LMCAS::detail::make_node<NumberNode>(value.value());
    }

    const auto slash = token.find('/');
    if (slash != std::string::npos) {
        if (slash == 0 || slash + 1 == token.size() ||
            token.find('/', slash + 1) != std::string::npos) {
            throw invalid();
        }
        return LMCAS::detail::make_node<NumberNode>(
            Rational(
                BigInt(token.substr(0, slash)),
                BigInt(token.substr(slash + 1))));
    }
    if (token.find_first_of(".eE") == std::string::npos) {
        return LMCAS::detail::make_node<NumberNode>(BigInt(token));
    }

    auto value = LMCAS::detail::parse_finite_decimal(token);
    if (!value) { throw invalid(); }
    return LMCAS::detail::make_node<NumberNode>(value.value());
}

namespace {


SymbolicExpr parse_serialized_numeric_atom(
    const std::string& expression, int line_num) {
    try {
        return LMCAS::detail::expression_from_node(
            parse_finite_number_node(expression, line_num));
    } catch (const std::bad_alloc&) {
        throw;
    } catch (const std::invalid_argument&) {
        throw;
    } catch (const std::exception& error) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": invalid numeric expression '" + expression +
            "': " + error.what());
    }
}

SymbolicExpr parse_general_expression(
    const std::string& expression, int line_num,
    ComputationContext& context) {
    auto parsed = parse_expr(expression, context);
    if (parsed) {
        auto normalized = LMCAS::simplify(parsed.value(), context);
        if (normalized) {
            return *normalized.value();
        }
        throw normalized.error();
    }
    if (parsed.error().code == CasErrc::ParseError ||
        parsed.error().code == CasErrc::InvalidArgument) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": " +
            parsed.error().message);
    }
    throw parsed.error();
}

}

SymbolicExpr parse_serialized_expression(
    const std::string& text, int line_num, ComputationContext& context) {
    const auto trimmed = trim_ascii_whitespace(text);
    if (trimmed.empty()) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": empty expression");
    }
    const std::string expression(trimmed);
    if (is_serialized_numeric_atom(trimmed)) {
        return parse_serialized_numeric_atom(expression, line_num);
    }
    return parse_general_expression(expression, line_num, context);
}

namespace {

std::shared_ptr<SymbolicExpr> parse_interval_expression(
    const std::string& text, int line_num, const std::string& interval,
    ComputationContext& context) {
    try {
        return std::make_shared<SymbolicExpr>(
            parse_serialized_expression(text, line_num, context));
    } catch (const std::invalid_argument&) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": malformed interval '" + interval + "'");
    }
}

Endpoint parse_interval_endpoint(
    const std::string& text, bool open, int line_num,
    const std::string& interval, ComputationContext& context) {
    auto value = parse_interval_expression(
        text, line_num, interval, context);
    return open ? Endpoint::open(std::move(value))
                : Endpoint::closed(std::move(value));
}

bool interval_delimiters_valid(const std::string& s) {
    if (s.size() < 5) return false;
    const bool lower = s.front() == '[' || s.front() == '(';
    const bool upper = s.back() == ']' || s.back() == ')';
    return lower && upper;
}

std::invalid_argument malformed_interval(
    const std::string& text, int line_num) {
    return std::invalid_argument(
        "Line " + std::to_string(line_num) +
        ": malformed interval '" + text + "'");
}


std::size_t find_interval_comma(
    const std::string& inner, const std::string& interval, int line_num) {
    std::vector<char> expected_closers;
    std::size_t comma_pos = std::string::npos;
    for (std::size_t i = 0; i < inner.size(); ++i) {
        const char c = inner[i];
        if (const char closer = matching_closer(c)) {
            expected_closers.push_back(closer);
            continue;
        }
        if (is_closing_delimiter(c)) {
            if (expected_closers.empty() || expected_closers.back() != c) {
                throw malformed_interval(interval, line_num);
            }
            expected_closers.pop_back();
            continue;
        }
        if (c != ',' || !expected_closers.empty()) {
            continue;
        }
        if (comma_pos != std::string::npos) {
            throw malformed_interval(interval, line_num);
        }
        comma_pos = i;
    }
    if (!expected_closers.empty() || comma_pos == std::string::npos) {
        throw malformed_interval(interval, line_num);
    }
    return comma_pos;
}

std::string trim_interval_endpoint(
    std::string value, const std::string& interval, int line_num) {
    const auto trimmed = trim_ascii_whitespace(value);
    if (trimmed.empty()) {
        throw malformed_interval(interval, line_num);
    }
    return std::string(trimmed);
}

bool is_infinity_token(const std::string& endpoint) {
    return endpoint == "nan" || endpoint == "inf" ||
        endpoint == "+inf" || endpoint == "-inf";
}

void validate_interval_endpoints(
    const std::string& lower, const std::string& upper,
    bool lower_open, bool upper_open,
    const std::string& interval, int line_num) {
    if ((lower == "-inf" && !lower_open) ||
        (upper == "+inf" && !upper_open) ||
        (is_infinity_token(lower) && lower != "-inf") ||
        (is_infinity_token(upper) && upper != "+inf")) {
        throw malformed_interval(interval, line_num);
    }
}

}

Interval parse_interval(
    const std::string& text, int line_num, ComputationContext& context) {
    if (!interval_delimiters_valid(text)) {
        throw malformed_interval(text, line_num);
    }

    const bool lower_open = text.front() == '(';
    const bool upper_open = text.back() == ')';
    const std::string inner = text.substr(1, text.size() - 2);
    const std::size_t comma_pos =
        find_interval_comma(inner, text, line_num);
    const std::string lower_text = trim_interval_endpoint(
        inner.substr(0, comma_pos), text, line_num);
    const std::string upper_text = trim_interval_endpoint(
        inner.substr(comma_pos + 1), text, line_num);
    validate_interval_endpoints(
        lower_text, upper_text, lower_open, upper_open, text, line_num);

    Interval interval;
    interval.lower = lower_text == "-inf"
        ? Endpoint::neg_inf()
        : parse_interval_endpoint(
            lower_text, lower_open, line_num, text, context);
    interval.upper = upper_text == "+inf"
        ? Endpoint::pos_inf()
        : parse_interval_endpoint(
            upper_text, upper_open, line_num, text, context);
    return interval;
}

}
