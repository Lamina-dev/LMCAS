#include "internal/assumption_parser.hpp"

namespace LMCAS::assumption_detail {

AssumptionContext AssumptionParser::parse(const std::string& data) {
    std::istringstream input(data);
    std::string line;
    while (std::getline(input, line)) {
        ++line_num;
        if (line.empty()) continue;
        if (line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (ended) {
            throw std::invalid_argument(
                "Line " + std::to_string(line_num) + ": data after END");
        }
        if (line == "END") {
            ended = true;
            continue;
        }
        parse_record(line);
    }
    if (!ended) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": unexpected end of input (missing END)");
    }
    return std::move(ctx);
}

void AssumptionParser::parse_record(const std::string& line) {
    std::istringstream ls(line);
    std::string keyword;
    ls >> keyword;
    using Handler = void (AssumptionParser::*)(std::istringstream&, const std::string&);
    struct Command {
        const char* name;
        Handler handler;
    };
    static const Command commands[] = {
        {"SCOPE", &AssumptionParser::parse_scope},
        {"DOMAIN", &AssumptionParser::parse_domain},
        {"SIGN", &AssumptionParser::parse_sign},
        {"PARITY", &AssumptionParser::parse_parity},
        {"BOUNDED", &AssumptionParser::parse_bounded},
        {"TRANSCENDENTAL", &AssumptionParser::parse_transcendental},
        {"FINITENESS", &AssumptionParser::parse_finiteness},
        {"DEFINITENESS", &AssumptionParser::parse_definiteness},
        {"PERIODIC", &AssumptionParser::parse_periodic},
        {"RELATION", &AssumptionParser::parse_relation},
        {"CONDITIONAL", &AssumptionParser::parse_conditional},
        {"CONTINUOUS", &AssumptionParser::parse_continuity},
        {"DIFFERENTIABLE", &AssumptionParser::parse_continuity},
        {"MONOTONICITY", &AssumptionParser::parse_monotonicity},
    };
    for (const auto& command : commands) {
        if (keyword == command.name) {
            if (command.handler != &AssumptionParser::parse_scope && current_scope < 0) {
                throw std::invalid_argument(
                    "Line " + std::to_string(line_num) + ": " + keyword + " before SCOPE");
            }
            (this->*command.handler)(ls, keyword);
            return;
        }
    }
    if (keyword == "END") {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": malformed END record");
    }
    throw std::invalid_argument(
        "Line " + std::to_string(line_num) + ": unknown keyword '" + keyword + "'");
}

void AssumptionParser::parse_scope(
    std::istringstream& ls, const std::string& keyword) {
    int idx;
    if (!(ls >> idx)) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": SCOPE missing index");
    }
    require_line_end(ls, line_num, keyword);
    if (idx != current_scope + 1) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": SCOPE indices must be sequential from zero");
    }
    if (idx > 0) {
        ctx.push();
    }
    current_scope = idx;
}


}

#include "computation_context.hpp"

namespace LMCAS {

AssumptionContext AssumptionContext::deserialize_impl(
    const std::string& data, ComputationContext& context) {
    assumption_detail::AssumptionParser parser(context);
    return parser.parse(data);
}

Result<AssumptionContext> AssumptionContext::deserialize(const std::string& data) {
    return deserialize_checked(data);
}

Result<AssumptionContext> AssumptionContext::deserialize_checked(
    const std::string& data) {
    ComputationContext context;
    return deserialize_checked(data, context);
}

Result<AssumptionContext> AssumptionContext::deserialize_checked(
    const std::string& data, ComputationContext& context) {
    constexpr const char* operation = "deserialize";
    auto input_budget = context.require_input_bytes(data.size(), operation);
    if (!input_budget) {
        return Result<AssumptionContext>::failure(input_budget.error());
    }
    try {
        return Result<AssumptionContext>::success(deserialize_impl(data, context));
    } catch (const CasError& error) {
        return Result<AssumptionContext>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<AssumptionContext>::failure(
            CasErrc::ResourceLimit,
            "assumption deserialization allocation failed", operation);
    } catch (const std::invalid_argument& ex) {
        return Result<AssumptionContext>::failure(
            CasErrc::ParseError, ex.what(), operation);
    } catch (const std::exception& ex) {
        return Result<AssumptionContext>::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

}
