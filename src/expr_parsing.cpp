#include "internal/expr_parser.hpp"
#include "computation_context.hpp"
#include "internal/expr_common.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/numeric_literal.hpp"

#include <cctype>
#include <exception>
#include <stdexcept>
#include <utility>

namespace LMCAS::expr_detail {

ExprResult ExprParser::parse() {
    skip_space();
    if (eof()) {
        return fail("empty expression");
    }
    auto result = parse_logical_or();
    if (!result) { return result; }
    skip_space();
    if (!eof()) {
        return fail("unexpected token '" + std::string(1, peek()) + "'");
    }
    return result;
}

bool ExprParser::eof() const noexcept { return pos_ >= source_.size(); }

char ExprParser::peek() const noexcept {
    return eof() ? '\0' : source_[pos_];
}

bool ExprParser::match(char c) {
    skip_space();
    if (peek() != c) { return false; }
    ++pos_;
    return true;
}

bool ExprParser::match_text(const char* text) {
    skip_space();
    const std::size_t start = pos_;
    for (const char* p = text; *p; ++p) {
        if (pos_ >= source_.size() || source_[pos_] != *p) {
            pos_ = start;
            return false;
        }
        ++pos_;
    }
    return true;
}

bool ExprParser::match_keyword(const char* text) {
    skip_space();
    const std::size_t start = pos_;
    if (!match_text(text)) { return false; }
    if (is_ident_continue(static_cast<unsigned char>(peek()))) {
        pos_ = start;
        return false;
    }
    return true;
}

void ExprParser::skip_space() noexcept {
    while (!eof() &&
           std::isspace(static_cast<unsigned char>(source_[pos_]))) {
        ++pos_;
    }
}

ExprResult ExprParser::fail(std::string message) const {
    return expr_detail::expr_common::expression_failure(CasErrc::ParseError,
    std::move(message) + " at byte " + std::to_string(pos_),
    kParseOperation);
}

bool ExprParser::is_ident_start(unsigned char c) noexcept {
    return std::isalpha(c) || c == '_' || c >= 0x80;
}

bool ExprParser::is_ident_continue(unsigned char c) noexcept {
    return std::isalnum(c) || c == '_' || c >= 0x80;
}

ExprResult ExprParser::parse_logical_or() {
    auto lhs = parse_logical_and();
    if (!lhs) { return lhs; }
    while (match_keyword("or")) {
        auto rhs = parse_logical_and();
        if (!rhs) { return rhs; }
        lhs = logical(lhs.value(), rhs.value(), LogicalNode::Op::Or);
        if (!lhs) { return lhs; }
    }
    return lhs;
}

ExprResult ExprParser::parse_logical_and() {
    auto lhs = parse_logical_not();
    if (!lhs) { return lhs; }
    while (match_keyword("and")) {
        auto rhs = parse_logical_not();
        if (!rhs) { return rhs; }
        lhs = logical(lhs.value(), rhs.value(), LogicalNode::Op::And);
        if (!lhs) { return lhs; }
    }
    return lhs;
}

ExprResult ExprParser::parse_logical_not() {
    if (match_keyword("not") || match('!')) {
        auto operand = parse_logical_not();
        if (!operand) { return operand; }
        return logical(operand.value(), nullptr, LogicalNode::Op::Not);
    }
    return parse_membership();
}

ExprResult ExprParser::parse_membership() {
    auto lhs = parse_equality();
    if (!lhs) { return lhs; }
    if (match_keyword("in")) {
        auto rhs = parse_equality();
        if (!rhs) { return rhs; }
        return membership(lhs.value(), rhs.value(), false);
    }
    const auto save = pos_;
    if (match_keyword("not")) {
        if (match_keyword("in")) {
            auto rhs = parse_equality();
            if (!rhs) { return rhs; }
            return membership(lhs.value(), rhs.value(), true);
        }
        pos_ = save;
    }
    return lhs;
}

ExprResult ExprParser::parse_equality() {
    auto lhs = parse_relational();
    if (!lhs) { return lhs; }
    while (true) {
        RelationOp op;
        if (match_text("==")) {
            op = RelationOp::EQ;
        } else if (match_text("!=")) {
            op = RelationOp::NEQ;
        } else {
            return lhs;
        }
        auto rhs = parse_relational();
        if (!rhs) { return rhs; }
        lhs = relational(lhs.value(), rhs.value(), op);
        if (!lhs) { return lhs; }
    }
}

ExprResult ExprParser::parse_relational() {
    auto lhs = parse_additive();
    if (!lhs) { return lhs; }
    while (true) {
        RelationOp op;
        if (match_text("<=")) {
            op = RelationOp::LEQ;
        } else if (match_text(">=")) {
            op = RelationOp::GEQ;
        } else if (match('<')) {
            op = RelationOp::LT;
        } else if (match('>')) {
            op = RelationOp::GT;
        } else {
            return lhs;
        }
        auto rhs = parse_additive();
        if (!rhs) { return rhs; }
        lhs = relational(lhs.value(), rhs.value(), op);
        if (!lhs) { return lhs; }
    }
}

ExprResult ExprParser::parse_additive() {
    auto lhs = parse_multiplicative();
    if (!lhs) { return lhs; }
    while (true) {
        if (match('+')) {
            auto rhs = parse_multiplicative();
            if (!rhs) { return rhs; }
            lhs = add(lhs.value(), rhs.value());
        } else if (match('-')) {
            auto rhs = parse_multiplicative();
            if (!rhs) { return rhs; }
            lhs = sub(lhs.value(), rhs.value());
        } else {
            return lhs;
        }
        if (!lhs) { return lhs; }
    }
}

ExprResult ExprParser::parse_multiplicative() {
    auto lhs = parse_power();
    if (!lhs) { return lhs; }
    while (true) {
        if (match_text("**")) {
            return fail("operator '**' is not supported; use '^'");
        } else if (match('*')) {
            auto rhs = parse_power();
            if (!rhs) { return rhs; }
            lhs = mul(lhs.value(), rhs.value());
        } else if (match('/')) {
            auto rhs = parse_power();
            if (!rhs) { return rhs; }
            lhs = div(lhs.value(), rhs.value());
        } else {
            return lhs;
        }
        if (!lhs) { return lhs; }
    }
}

ExprResult ExprParser::parse_power() {
    auto base = parse_unary();
    if (!base) { return base; }
    if (match_text("**")) {
        return fail("operator '**' is not supported; use '^'");
    }
    if (match('^')) {
        auto exponent = parse_power();
        if (!exponent) { return exponent; }
        return pow(base.value(), exponent.value());
    }
    return base;
}

ExprResult ExprParser::parse_unary() {
    if (match('+')) { return parse_unary(); }
    if (match('-')) {
        auto operand = parse_unary();
        if (!operand) { return operand; }
        return neg(operand.value());
    }
    return parse_primary();
}

ExprResult ExprParser::parse_primary() {
    skip_space();
    if (match('(')) {
        auto inner = parse_logical_or();
        if (!inner) { return inner; }
        if (match(',')) {
            auto rhs = parse_logical_or();
            if (!rhs) { return rhs; }
            if (!match(')')) { return fail("expected ')' after interval"); }
            return interval(inner.value(), rhs.value(), false, false);
        }
        if (!match(')')) { return fail("expected ')'"); }
        return inner;
    }
    if (match('{')) {
        return parse_set_literal();
    }
    if (match('[')) {
        return parse_interval_literal(true);
    }
    if (std::isdigit(static_cast<unsigned char>(peek())) || peek() == '.') {
        return parse_number();
    }
    if (is_ident_start(static_cast<unsigned char>(peek()))) {
        return parse_identifier_or_call();
    }
    return fail("expected expression");
}

ExprResult ExprParser::parse_set_literal() {
    std::vector<ExprPtr> elements;
    if (!match('}')) {
        while (true) {
            auto element = parse_logical_or();
            if (!element) { return element; }
            elements.push_back(element.value());
            if (match('}')) break;
            if (!match(',')) { return fail("expected ',' or '}' in set literal"); }
        }
    }
    return finite_set(elements);
}

ExprResult ExprParser::parse_interval_literal(bool lower_closed) {
    auto lower = parse_logical_or();
    if (!lower) { return lower; }
    if (!match(',')) { return fail("expected ',' in interval literal"); }
    auto upper = parse_logical_or();
    if (!upper) { return upper; }
    bool upper_closed;
    if (match(']')) {
        upper_closed = true;
    } else if (match(')')) {
        upper_closed = false;
    } else {
        return fail("expected ']' or ')' after interval literal");
    }
    return interval(lower.value(), upper.value(), lower_closed, upper_closed);
}

ExprResult ExprParser::parse_number() {
    skip_space();
    const std::size_t start = pos_;
    bool saw_digit = false;
    while (std::isdigit(static_cast<unsigned char>(peek()))) {
        saw_digit = true;
        ++pos_;
    }
    if (peek() == '.') {
        ++pos_;
        while (std::isdigit(static_cast<unsigned char>(peek()))) {
            saw_digit = true;
            ++pos_;
        }
    }
    if (!saw_digit) { return fail("invalid number literal"); }
    if (peek() == 'e' || peek() == 'E') {
        ++pos_;
        if (peek() == '+' || peek() == '-') ++pos_;
        bool saw_exp_digit = false;
        while (std::isdigit(static_cast<unsigned char>(peek()))) {
            saw_exp_digit = true;
            ++pos_;
        }
        if (!saw_exp_digit) { return fail("invalid number exponent"); }
    }
    try {
        return rational(Rational(source_.substr(start, pos_ - start)));
    } catch (const std::bad_alloc&) {
        return expr_detail::expr_common::expression_failure(CasErrc::ResourceLimit,
                                  "number literal allocation failed",
                                  kParseOperation);
    } catch (const std::invalid_argument& error) {
        return expr_detail::expr_common::expression_failure(
            CasErrc::ParseError, error.what(), kParseOperation);
    } catch (const std::out_of_range& error) {
        return expr_detail::expr_common::expression_failure(
            CasErrc::ParseError, error.what(), kParseOperation);
    } catch (const std::exception& error) {
        return expr_detail::expr_common::expression_failure(
            CasErrc::InternalInvariant, error.what(), kParseOperation);
    }
}

std::string ExprParser::parse_identifier() {
    skip_space();
    const std::size_t start = pos_;
    if (!is_ident_start(static_cast<unsigned char>(peek()))) { return {}; }
    ++pos_;
    while (is_ident_continue(static_cast<unsigned char>(peek()))) ++pos_;
    return source_.substr(start, pos_ - start);
}

ExprResult ExprParser::parse_approx_literal() {
    skip_space();
    const auto start = pos_;
    while (!eof() && peek() != ')' &&
           !std::isspace(static_cast<unsigned char>(peek()))) ++pos_;
    const std::string_view token(source_.data() + start, pos_ - start);
    if (!match(')')) return fail("expected ')' after approximate decimal");
    auto value = detail::parse_finite_decimal(token);
    if (!value) return fail(value.error().message);
    return approx_real(value.value());
}

ExprResult ExprParser::parse_call(const std::string& name) {
    if (name == "approx") return parse_approx_literal();
    std::vector<ExprPtr> arguments;
    bool closed = match(')');
    if (!closed) {
        do {
            auto argument = parse_logical_or();
            if (!argument) {
                return argument;
            }
            arguments.push_back(argument.value());
            if (match(')')) {
                closed = true;
                break;
            }
        } while (match(','));
    }
    if (!closed) {
        return fail("expected ')'");
    }
    return apply_function(name, arguments);
}

ExprResult ExprParser::parse_identifier_or_call() {
    auto name = parse_identifier();
    if (name.empty()) {
        return fail("expected identifier");
    }
    skip_space();
    if (match('(')) {
        return parse_call(name);
    }
    if (name == "pi" || name == "π") {
        return pi();
    }
    if (name == "e") {
        return e();
    }
    if (name == "phi") {
        return phi();
    }
    if (detail::is_imaginary_unit_name(name)) {
        return imaginary_unit();
    }
    return sym(name);
}

}

namespace LMCAS {

using expr_detail::ExprParser;
using expr_detail::kParseOperation;

ExprResult parse_expr(const std::string& source) {
    ComputationContext context;
    return parse_expr(source, context);
}

ExprResult parse_expr(const std::string& source,
                      ComputationContext& context) {
    auto input_budget =
        context.require_input_bytes(source.size(), kParseOperation);
    if (!input_budget) {
        return ExprResult::failure(input_budget.error());
    }
    try {
        return ExprParser(source).parse();
    } catch (const std::bad_alloc&) {
        return expr_detail::expr_common::expression_failure(CasErrc::ResourceLimit,
                                  "expression parse allocation failed",
                                  kParseOperation);
    } catch (const std::exception& error) {
        return expr_detail::expr_common::expression_failure(
            CasErrc::InternalInvariant, error.what(), kParseOperation);
    }
}

}
