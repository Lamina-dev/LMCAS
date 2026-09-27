#include "interval.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/numeric_literal.hpp"
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace LMCAS {

std::string IntervalUnion::to_string() const {

    if (intervals_.empty()) {
        return "\xe2\x88\x85";
    }

    if (intervals_.size() == 1 && intervals_[0].is_entire_line()) {
        return "(-\xe2\x88\x9e, +\xe2\x88\x9e)";
    }

    std::string result;
    for (size_t i = 0; i < intervals_.size(); ++i) {
        if (i > 0) {
            result += " \xe2\x88\xaa ";
        }

        const Interval& iv = intervals_[i];

        if (iv.lower.is_open) {
            result += "(";
        } else {
            result += "[";
        }

        if (iv.lower.is_neg_infinity) {
            result += "-\xe2\x88\x9e";
        } else if (iv.lower.value) {
            result += iv.lower.value->to_string();
        } else {
            result += "0";
        }

        result += ", ";

        if (iv.upper.is_pos_infinity) {
            result += "+\xe2\x88\x9e";
        } else if (iv.upper.value) {
            result += iv.upper.value->to_string();
        } else {
            result += "0";
        }

        if (iv.upper.is_open) {
            result += ")";
        } else {
            result += "]";
        }
    }

    return result;
}

static bool starts_with_at(const std::string& str, size_t pos, const std::string& prefix) {
    if (pos + prefix.size() > str.size()) {
        return false;
    }
    return str.compare(pos, prefix.size(), prefix) == 0;
}

static size_t skip_ws(const std::string& str, size_t pos) {
    while (pos < str.size() && std::isspace(static_cast<unsigned char>(str[pos]))) {
        ++pos;
    }
    return pos;
}

static bool consume_digits(const std::string& str, size_t& pos) {
    const size_t start = pos;
    while (pos < str.size() &&
           std::isdigit(static_cast<unsigned char>(str[pos]))) {
        ++pos;
    }
    return pos != start;
}

static void consume_sign(const std::string& str, size_t& pos) {
    if (pos < str.size() && (str[pos] == '-' || str[pos] == '+')) ++pos;
}

static std::shared_ptr<SymbolicExpr> parse_fraction(
    const std::string& str, size_t start, size_t& pos) {
    const size_t slash = pos++;
    const size_t denominator_start = pos;
    if (!consume_digits(str, pos)) {
        pos = start;
        return nullptr;
    }
    const BigInt numerator(str.substr(start, slash - start));
    const BigInt denominator(str.substr(denominator_start, pos - denominator_start));
    if (denominator.is_zero()) {
        pos = start;
        return nullptr;
    }
    return SymbolicExpr::number(Rational(numerator, denominator));
}

static std::shared_ptr<SymbolicExpr> numeric_token_value(
    const std::string& str, size_t start, size_t& pos, bool approximate) {
    const std::string token = str.substr(start, pos - start);
    if (!approximate) {
        return SymbolicExpr::number(BigInt(token));
    }
    auto value = detail::parse_finite_decimal(token);
    if (!value) {
        pos = start;
        return nullptr;
    }
    return SymbolicExpr::number(value.value());
}

static std::shared_ptr<SymbolicExpr> parse_approximate_atom(
    const std::string& str, size_t& pos) {
    const size_t start = pos;
    pos = skip_ws(str, pos + 7);
    const auto token_start = pos;
    while (pos < str.size() && str[pos] != ')' &&
           !std::isspace(static_cast<unsigned char>(str[pos]))) { ++pos; }
    const std::string_view token(str.data() + token_start, pos - token_start);
    pos = skip_ws(str, pos);
    if (pos == str.size() || str[pos] != ')') {
        pos = start;
        return nullptr;
    }
    ++pos; /**< 消耗原子的结束分隔符，保留区间的分隔符。 */
    auto value = detail::parse_finite_decimal(token);
    if (!value) {
        pos = start;
        return nullptr;
    }
    return SymbolicExpr::number(value.value());
}

static std::shared_ptr<SymbolicExpr> parse_numeric_value(
    const std::string& str, size_t& pos) {
    const size_t start = pos;
    if (starts_with_at(str, pos, "approx(")) {
        return parse_approximate_atom(str, pos);
    }
    consume_sign(str, pos);
    bool has_digit = consume_digits(str, pos);
    if (has_digit && pos < str.size() && str[pos] == '/') {
        return parse_fraction(str, start, pos);
    }
    bool approximate = false;
    if (pos < str.size() && str[pos] == '.') {
        approximate = true;
        ++pos;
        has_digit = consume_digits(str, pos) || has_digit;
    }
    if (!has_digit) {
        pos = start;
        return nullptr;
    }
    if (pos < str.size() && (str[pos] == 'e' || str[pos] == 'E')) {
        approximate = true;
        ++pos;
        consume_sign(str, pos);
        if (!consume_digits(str, pos)) {
            pos = start;
            return nullptr;
        }
    }
    return numeric_token_value(str, start, pos, approximate);
}

static bool parse_endpoint_value(const std::string& str, size_t& pos, Endpoint& ep) {
    pos = skip_ws(str, pos);

    if (starts_with_at(str, pos, "-\xe2\x88\x9e")) {
        ep.is_neg_infinity = true;
        ep.is_pos_infinity = false;
        ep.is_open = true;
        ep.value = nullptr;
        pos += 4;
        return true;
    }

    if (starts_with_at(str, pos, "+\xe2\x88\x9e")) {
        ep.is_pos_infinity = true;
        ep.is_neg_infinity = false;
        ep.is_open = true;
        ep.value = nullptr;
        pos += 4;
        return true;
    }

    auto val = parse_numeric_value(str, pos);
    if (!val) {
        return false;
    }

    ep.value = val;
    ep.is_neg_infinity = false;
    ep.is_pos_infinity = false;

    return true;
}

static bool parse_interval_lower(
    const std::string& str, size_t& pos, Endpoint& lower) {
    pos = skip_ws(str, pos);
    if (pos >= str.size()) {
        return false;
    }
    const char bracket = str[pos];
    if (bracket != '(' && bracket != '[') {
        return false;
    }
    const bool open = bracket == '(';
    ++pos;
    if (!parse_endpoint_value(str, pos, lower)) {
        return false;
    }
    if (lower.is_pos_infinity || (lower.is_neg_infinity && !open)) {
        return false;
    }
    if (!lower.is_neg_infinity) lower.is_open = open;
    return true;
}

static bool parse_interval_upper(
    const std::string& str, size_t& pos, Endpoint& upper) {
    if (!parse_endpoint_value(str, pos, upper)) {
        return false;
    }
    pos = skip_ws(str, pos);
    if (pos >= str.size()) {
        return false;
    }
    const char bracket = str[pos];
    if (bracket != ')' && bracket != ']') {
        return false;
    }
    const bool open = bracket == ')';
    ++pos;
    if (upper.is_neg_infinity || (upper.is_pos_infinity && !open)) {
        return false;
    }
    if (!upper.is_pos_infinity) upper.is_open = open;
    return true;
}

static bool parse_single_interval(
    const std::string& str, size_t& pos, Interval& iv) {
    Endpoint lower;
    if (!parse_interval_lower(str, pos, lower)) {
        return false;
    }
    pos = skip_ws(str, pos);
    if (pos >= str.size() || str[pos] != ',') {
        return false;
    }
    ++pos;
    Endpoint upper;
    if (!parse_interval_upper(str, pos, upper)) {
        return false;
    }
    iv.lower = std::move(lower);
    iv.upper = std::move(upper);
    return true;
}

std::optional<IntervalUnion> IntervalUnion::parse(const std::string& str) {
    if (str.empty()) {
        return std::nullopt;
    }

    if (str == "\xe2\x88\x85") {
        return IntervalUnion::empty();
    }

    try {
        std::vector<Interval> intervals;
        size_t pos = 0;

        Interval iv;
        if (!parse_single_interval(str, pos, iv)) {
            return std::nullopt;
        }
        intervals.push_back(std::move(iv));

        while (pos < str.size()) {
            pos = skip_ws(str, pos);
            if (pos >= str.size()) break;

            if (!starts_with_at(str, pos, "\xe2\x88\xaa")) {
                return std::nullopt;
            }
            pos += 3;
            pos = skip_ws(str, pos);

            Interval next_iv;
            if (!parse_single_interval(str, pos, next_iv)) {
                return std::nullopt;
            }
            intervals.push_back(std::move(next_iv));
        }

        auto normalized =
            IntervalUnion::from_intervals_checked(std::move(intervals));
        if (!normalized) {
            return std::nullopt;
        }
        return std::move(normalized.value());
    } catch (const std::invalid_argument&) {
        return std::nullopt;
    } catch (const std::out_of_range&) {
        return std::nullopt;
    }
}

static std::shared_ptr<const SymbolicNode> interval_to_expr(
    const Interval& interval, const std::shared_ptr<const SymbolicNode>& variable) {
    std::shared_ptr<const SymbolicNode> lower;
    std::shared_ptr<const SymbolicNode> upper;
    if (!interval.lower.is_neg_infinity) {
        auto bound = interval.lower.value ? LMCAS::detail::node(interval.lower.value)
            : LMCAS::detail::make_node<NumberNode>(0.0);
        const auto op = interval.lower.is_open ? RelationalNode::Op::GT
                                               : RelationalNode::Op::GEQ;
        lower = LMCAS::detail::make_node<RelationalNode>(variable, bound, op);
    }
    if (!interval.upper.is_pos_infinity) {
        auto bound = interval.upper.value ? LMCAS::detail::node(interval.upper.value)
            : LMCAS::detail::make_node<NumberNode>(0.0);
        const auto op = interval.upper.is_open ? RelationalNode::Op::LT
                                               : RelationalNode::Op::LEQ;
        upper = LMCAS::detail::make_node<RelationalNode>(variable, bound, op);
    }
    if (lower && upper) {
        return LMCAS::detail::make_node<LogicalNode>(lower, upper, LogicalNode::Op::And);
    }
    return lower ? lower : upper;
}

std::shared_ptr<SymbolicExpr> IntervalUnion::to_expr(const std::string& var) const {
    if (is_empty()) {
        return LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<FiniteSetNode>(
                std::vector<std::shared_ptr<const SymbolicNode>>{}));
    }

    auto var_node = LMCAS::detail::node(SymbolicExpr::variable(var));
    if (is_entire_line()) {
        auto zero = LMCAS::detail::node(SymbolicExpr::number(0));
        return LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<RelationalNode>(
                zero, zero, RelationalNode::Op::EQ));
    }

    if (intervals_.size() == 1) {
        auto node = interval_to_expr(intervals_[0], var_node);
        if (!node) {
            return nullptr;
        }
        return LMCAS::detail::make_expression_ptr(node);
    }

    auto result = interval_to_expr(intervals_[0], var_node);
    for (size_t i = 1; i < intervals_.size(); ++i) {
        auto next = interval_to_expr(intervals_[i], var_node);
        if (result && next) {
            result = LMCAS::detail::make_node<LogicalNode>(result, next, LogicalNode::Op::Or);
        } else if (next) {
            result = next;
        }

    }

    if (!result) {
        return nullptr;
    }
    return LMCAS::detail::make_expression_ptr(result);
}

}
