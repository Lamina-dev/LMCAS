#include "internal/assumption_parser.hpp"

namespace LMCAS::assumption_detail {

void AssumptionParser::parse_domain(
    std::istringstream& ls, const std::string& keyword) {

    std::string sym, dom_str;
    if (!(ls >> sym >> dom_str)) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": DOMAIN requires symbol and domain");
    }
    require_line_end(ls, line_num, keyword);
    Domain dom = assumption_detail::parse_domain(dom_str, line_num);
    require_deserialization_update(
        ctx.current_properties().declare_domain_checked(sym, dom), line_num, keyword);
}

void AssumptionParser::parse_sign(
    std::istringstream& ls, const std::string& keyword) {

    std::string sym, sign_str;
    if (!(ls >> sym >> sign_str)) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": SIGN requires symbol and sign");
    }
    require_line_end(ls, line_num, keyword);
    Sign s = assumption_detail::parse_sign(sign_str, line_num);
    require_deserialization_update(
        ctx.current_properties().declare_sign_checked(sym, s), line_num, keyword);
}

void AssumptionParser::parse_parity(
    std::istringstream& ls, const std::string& keyword) {

    std::string sym, par_str;
    if (!(ls >> sym >> par_str)) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": PARITY requires symbol and parity");
    }
    require_line_end(ls, line_num, keyword);
    Parity p = assumption_detail::parse_parity(par_str, line_num);
    require_deserialization_update(
        ctx.current_properties().declare_parity_checked(sym, p), line_num, keyword);
}

void AssumptionParser::parse_bounded(
    std::istringstream& ls, const std::string& keyword) {
    std::string symbol;
    std::string boundedness_text;
    if (!(ls >> symbol >> boundedness_text)) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) +
            ": BOUNDED requires symbol and boundedness");
    }
    const Boundedness bounded =
        parse_boundedness(boundedness_text, line_num);

    std::string payload;
    std::getline(ls, payload);
    const auto bounds_text = trim_ascii_whitespace(payload);
    std::optional<Interval> bounds;
    if (!bounds_text.empty()) {
        if (bounded != Boundedness::Bounded) {
            throw std::invalid_argument(
                "Line " + std::to_string(line_num) +
                ": only Bounded declarations may carry an interval");
        }
        bounds = parse_interval(
            std::string(bounds_text), line_num, computation);
    }
    require_deserialization_update(
        ctx.current_properties().declare_bounded_checked(
            symbol, bounded, std::move(bounds), computation),
        line_num, keyword);
}

void AssumptionParser::parse_transcendental(
    std::istringstream& ls, const std::string& keyword) {

    std::string sym;
    if (!(ls >> sym)) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": TRANSCENDENTAL requires symbol");
    }
    require_line_end(ls, line_num, keyword);
    require_deserialization_update(
        ctx.current_properties().declare_transcendental_checked(sym), line_num, keyword);
}

void AssumptionParser::parse_finiteness(
    std::istringstream& ls, const std::string& keyword) {

    std::string sym, fin_str;
    if (!(ls >> sym >> fin_str)) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": FINITENESS requires symbol and value");
    }
    Finiteness f = assumption_detail::parse_finiteness(fin_str, line_num);
    require_line_end(ls, line_num, keyword);
    require_deserialization_update(
        ctx.current_properties().declare_finiteness_checked(sym, f), line_num, keyword);
}

void AssumptionParser::parse_definiteness(
    std::istringstream& ls, const std::string& keyword) {

    std::string sym, def_str;
    if (!(ls >> sym >> def_str)) {
        throw std::invalid_argument(
            "Line " + std::to_string(line_num) + ": DEFINITENESS requires symbol and value");
    }
    Definiteness d = assumption_detail::parse_definiteness(def_str, line_num);
    require_line_end(ls, line_num, keyword);
    require_deserialization_update(
        ctx.current_properties().declare_definiteness_checked(sym, d), line_num, keyword);
}


}
