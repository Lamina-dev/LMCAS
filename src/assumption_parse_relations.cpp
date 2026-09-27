#include "internal/assumption_parser.hpp"
#include <cctype>
#include <string_view>
#include <vector>

namespace LMCAS::assumption_detail {
namespace {

struct TopLevelToken {
    std::size_t position;
    std::string_view token;
};


bool consume_delimiter(
    char value, std::vector<char>& expected_closers, int line_num) {
    if (const char closer = matching_closer(value)) {
        expected_closers.push_back(closer);
        return true;
    }
    if (!is_closing_delimiter(value)) return false;
    if (expected_closers.empty() || expected_closers.back() != value) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": mismatched expression delimiter");
    }
    expected_closers.pop_back();
    return true;
}

bool token_is_separated(
    const std::string& text, std::size_t position,
    std::string_view token) {
    if (position == 0 || position + token.size() >= text.size()) {
        return false;
    }
    return std::isspace(
               static_cast<unsigned char>(text[position - 1])) != 0 &&
        std::isspace(static_cast<unsigned char>(
            text[position + token.size()])) != 0;
}

std::string_view top_level_token_at(
    const std::string& text, std::size_t position,
    const std::vector<std::string_view>& tokens) {
    for (const auto token : tokens) {
        if (position + token.size() <= text.size() &&
            text.compare(
                position, token.size(), token.data(), token.size()) == 0 &&
            token_is_separated(text, position, token)) {
            return token;
        }
    }
    return {};
}

std::vector<TopLevelToken> find_top_level_tokens(
    const std::string& text,
    const std::vector<std::string_view>& tokens,
    int line_num) {
    std::vector<char> expected_closers;
    std::vector<TopLevelToken> matches;
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (consume_delimiter(text[i], expected_closers, line_num) ||
            !expected_closers.empty()) {
            continue;
        }
        const auto token = top_level_token_at(text, i, tokens);
        if (token.empty()) {
            continue;
        }
        matches.push_back({i, token});
        i += token.size() - 1;
    }
    if (!expected_closers.empty()) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": unclosed expression delimiter");
    }
    return matches;
}

std::string trim_copy(const std::string& text, int line_num) {
    const auto trimmed = trim_ascii_whitespace(text);
    if (trimmed.empty()) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": empty expression");
    }
    return std::string(trimmed);
}

std::string unwrap_parenthesized(
    const std::string& text, int line_num, const char* role) {
    const std::string trimmed = trim_copy(text, line_num);
    if (trimmed.size() < 2 || trimmed.front() != '(' ||
        trimmed.back() != ')') {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": CONDITIONAL " + role + " must be parenthesized");
    }
    int depth = 0;
    for (std::size_t i = 0; i < trimmed.size(); ++i) {
        if (trimmed[i] == '(') {
            ++depth;
        } else if (trimmed[i] == ')') {
            --depth;
            if (depth == 0 && i + 1 != trimmed.size()) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) +
                    ": CONDITIONAL has text outside " + role);
            }
        }
    }
    return trimmed.substr(1, trimmed.size() - 2);
}

SymbolicExpr parse_conditional_relation(
    const std::string& text, int line_num,
    ComputationContext& context) {
    auto relation = parse_relation(text, line_num, true, context);
    return LMCAS::detail::expression_from_node(
        LMCAS::detail::make_node<RelationalNode>(
            LMCAS::detail::node(relation.lhs),
            LMCAS::detail::node(relation.rhs), relation.op));
}

}

Relation parse_relation(
    const std::string& text, int line_num, bool conditional,
    ComputationContext& context) {
    static const std::vector<std::string_view> operators = {
        "GEQ", "LEQ", "NEQ", "GT", "LT", "EQ",
        ">=", "<=", "!=", "==", ">", "<"};
    const auto matches =
        find_top_level_tokens(text, operators, line_num);
    if (matches.size() != 1) {
        const std::string message = conditional
            ? "CONDITIONAL cannot parse expression '" + text + "'"
            : "RELATION requires exactly one top-level operator";
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": " + message);
    }

    const auto match = matches.front();
    const std::string lhs_text =
        trim_copy(text.substr(0, match.position), line_num);
    const std::string rhs_text = trim_copy(
        text.substr(match.position + match.token.size()), line_num);
    return {
        parse_serialized_expression(lhs_text, line_num, context),
        parse_serialized_expression(rhs_text, line_num, context),
        parse_relop(std::string(match.token), line_num)};
}

void AssumptionParser::parse_relation(
    std::istringstream& input, const std::string& keyword) {
    std::string text;
    std::getline(input, text);
    const auto trimmed = trim_ascii_whitespace(text);
    if (trimmed.empty()) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": RELATION requires lhs, op, and rhs");
    }
    auto relation = assumption_detail::parse_relation(
        std::string(trimmed), line_num, false, computation);
    require_deserialization_update(
        ctx.current_relations().add_relation_checked(
            relation.lhs, relation.rhs, relation.op,
            ctx.current_properties()),
        line_num, keyword);
}

void AssumptionParser::parse_conditional(
    std::istringstream& input, const std::string& keyword) {
    std::string text;
    std::getline(input, text);
    const auto trimmed = trim_ascii_whitespace(text);
    if (trimmed.empty()) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": CONDITIONAL requires condition and conclusion");
    }
    text = trimmed;
    const auto arrows =
        find_top_level_tokens(text, {"=>"}, line_num);
    if (arrows.size() != 1) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": CONDITIONAL requires one top-level '=>'");
    }

    const auto arrow = arrows.front();
    const std::string condition_text = unwrap_parenthesized(
        text.substr(0, arrow.position), line_num, "condition");
    const std::string conclusion_text = unwrap_parenthesized(
        text.substr(arrow.position + arrow.token.size()),
        line_num, "conclusion");
    SymbolicExpr condition = parse_conditional_relation(
        condition_text, line_num, computation);
    SymbolicExpr conclusion = parse_conditional_relation(
        conclusion_text, line_num, computation);
    require_deserialization_update(
        ctx.assume_conditional_checked(condition, conclusion),
        line_num, keyword);
}

}
