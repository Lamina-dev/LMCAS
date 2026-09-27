#include "internal/assumption_parser.hpp"

namespace LMCAS::assumption_detail {

namespace {

template <typename T, std::size_t N>
T parse_named_value(const std::string& text, int line_num, const char* kind,
                    const std::pair<const char*, T> (&values)[N]) {
    for (const auto& value : values) {
        if (text == value.first) return value.second;
    }
    throw std::invalid_argument(
        "Line " + std::to_string(line_num) + ": unknown " + kind + " '" + text + "'");
}

}

Domain parse_domain(const std::string& s, int line_num) {
    static const std::pair<const char*, Domain> values[] = {
        {"Complex", Domain::Complex},
        {"Real", Domain::Real},
        {"Algebraic", Domain::Algebraic},
        {"Rational", Domain::Rational},
        {"Integer", Domain::Integer},
        {"Natural", Domain::Natural},
        {"PositiveInt", Domain::PositiveInt},
    };
    return parse_named_value(s, line_num, "domain", values);
}

Sign parse_sign(const std::string& s, int line_num) {
    static const std::pair<const char*, Sign> values[] = {
        {"Positive", Sign::Positive},
        {"Negative", Sign::Negative},
        {"NonNegative", Sign::NonNegative},
        {"NonPositive", Sign::NonPositive},
        {"Zero", Sign::Zero},
        {"NonZero", Sign::NonZero},
    };
    return parse_named_value(s, line_num, "sign", values);
}

Parity parse_parity(const std::string& s, int line_num) {
    static const std::pair<const char*, Parity> values[] = {
        {"Even", Parity::Even},
        {"Odd", Parity::Odd},
    };
    return parse_named_value(s, line_num, "parity", values);
}

Boundedness parse_boundedness(const std::string& s, int line_num) {
    static const std::pair<const char*, Boundedness> values[] = {
        {"Bounded", Boundedness::Bounded},
        {"Unbounded", Boundedness::Unbounded},
    };
    return parse_named_value(s, line_num, "boundedness", values);
}

Finiteness parse_finiteness(const std::string& s, int line_num) {
    static const std::pair<const char*, Finiteness> values[] = {
        {"Finite", Finiteness::Finite},
        {"Divergent", Finiteness::Divergent},
    };
    return parse_named_value(s, line_num, "finiteness", values);
}

Definiteness parse_definiteness(const std::string& s, int line_num) {
    static const std::pair<const char*, Definiteness> values[] = {
        {"PositiveDefinite", Definiteness::PositiveDefinite},
        {"PositiveSemiDefinite", Definiteness::PositiveSemiDefinite},
        {"NegativeDefinite", Definiteness::NegativeDefinite},
        {"NegativeSemiDefinite", Definiteness::NegativeSemiDefinite},
        {"Indefinite", Definiteness::Indefinite},
    };
    return parse_named_value(s, line_num, "definiteness", values);
}

Monotonicity parse_monotonicity(const std::string& s, int line_num) {
    static const std::pair<const char*, Monotonicity> values[] = {
        {"Increasing", Monotonicity::Increasing},
        {"Decreasing", Monotonicity::Decreasing},
        {"NonDecreasing", Monotonicity::NonDecreasing},
        {"NonIncreasing", Monotonicity::NonIncreasing},
    };
    return parse_named_value(s, line_num, "monotonicity", values);
}

RelationOp parse_relop(const std::string& s, int line_num) {
    static const std::pair<const char*, RelationOp> values[] = {
        {"GT", RelationOp::GT},
        {"LT", RelationOp::LT},
        {"GEQ", RelationOp::GEQ},
        {"LEQ", RelationOp::LEQ},
        {"NEQ", RelationOp::NEQ},
        {"EQ", RelationOp::EQ},
        {">", RelationOp::GT},
        {"<", RelationOp::LT},
        {">=", RelationOp::GEQ},
        {"<=", RelationOp::LEQ},
        {"!=", RelationOp::NEQ},
        {"==", RelationOp::EQ},
    };
    return parse_named_value(s, line_num, "relational operator", values);
}


}
