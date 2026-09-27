#include "internal/assumption_parser.hpp"

namespace LMCAS::assumption_detail {

void AssumptionParser::parse_periodic(
    std::istringstream& ls, const std::string& keyword) {
    std::string sym, variable;
    if (!(ls >> sym >> variable)) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": PERIODIC requires symbol and variable");
    }
    std::string period_str;
    std::getline(ls, period_str);
    const auto period_text = trim_ascii_whitespace(period_str);
    if (period_text.empty()) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": PERIODIC requires period value");
    }
    auto period = parse_serialized_expression(
        std::string(period_text), line_num, computation);
    require_deserialization_update(
        ctx.current_properties().declare_periodic_checked(
            sym, variable, period),
        line_num, keyword);
}

void AssumptionParser::parse_continuity(
    std::istringstream& ls, const std::string& keyword) {
    std::string sym;
    if (!(ls >> sym)) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": " + keyword + " requires symbol");
    }
    std::string iv_str;
    std::getline(ls, iv_str);
    const auto interval_text = trim_ascii_whitespace(iv_str);
    if (interval_text.empty()) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": " + keyword + " requires interval");
    }
    Interval iv = parse_interval(
        std::string(interval_text), line_num, computation);
    auto declared = keyword == "DIFFERENTIABLE"
        ? ctx.current_properties().declare_differentiable_checked(
            sym, iv, computation)
        : ctx.current_properties().declare_continuous_checked(
            sym, iv, computation);
    require_deserialization_update(declared, line_num, keyword);
}


void AssumptionParser::parse_monotonicity(
    std::istringstream& ls, const std::string& keyword) {

    std::string sym, var;
    if (!(ls >> sym >> var)) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": MONOTONICITY requires symbol and variable");
    }
    std::string rest_line;
    std::getline(ls, rest_line);
    const auto monotonicity_text = trim_ascii_whitespace(rest_line);
    if (monotonicity_text.empty()) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": MONOTONICITY requires interval and type");
    }
    rest_line = monotonicity_text;
    size_t bracket_end = std::string::npos;
    if (rest_line.front() == '[' || rest_line.front() == '(') {
        size_t depth = 0;
        for (size_t i = 1; i < rest_line.size(); ++i) {
            const char c = rest_line[i];
            if (c == '(' || c == '[') {
                ++depth;
            } else if (c == ')' || c == ']') {
                if (depth == 0) {
                    bracket_end = i;
                    break;
                }
                --depth;
            }
        }
    }
    if (bracket_end == std::string::npos) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": MONOTONICITY malformed interval");
    }
    std::string mono_iv_str = rest_line.substr(0, bracket_end + 1);
    Interval iv = parse_interval(mono_iv_str, line_num, computation);

    std::string mono_str = rest_line.substr(bracket_end + 1);
    const auto type_text = trim_ascii_whitespace(mono_str);
    if (type_text.empty()) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": MONOTONICITY requires type after interval");
    }
    mono_str = type_text;

    Monotonicity m = assumption_detail::parse_monotonicity(mono_str, line_num);
    require_deserialization_update(
        ctx.current_properties().declare_monotonicity_checked(
            sym, var, iv, m, computation),
        line_num, keyword);
}


}
