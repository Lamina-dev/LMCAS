#pragma once

#include "assumption_context.hpp"
#include "internal/symbolic_ast.hpp"
#include <sstream>
#include <string_view>

namespace LMCAS::assumption_detail {

std::string_view trim_ascii_whitespace(std::string_view text) noexcept;
bool is_serialized_numeric_atom(std::string_view expression) noexcept;
char matching_closer(char opener) noexcept;
bool is_closing_delimiter(char value) noexcept;

void require_deserialization_update(const Result<void>& result, int line,
                                    const std::string& keyword);
void require_line_end(std::istringstream& input, int line_num,
                      const std::string& keyword);
Domain parse_domain(const std::string& s, int line_num);
Sign parse_sign(const std::string& s, int line_num);
Parity parse_parity(const std::string& s, int line_num);
Boundedness parse_boundedness(const std::string& s, int line_num);
Finiteness parse_finiteness(const std::string& s, int line_num);
Definiteness parse_definiteness(const std::string& s, int line_num);
Monotonicity parse_monotonicity(const std::string& s, int line_num);
RelationOp parse_relop(const std::string& s, int line_num);
SymbolicExpr parse_serialized_expression(
    const std::string& text, int line_num, ComputationContext& context);
Interval parse_interval(
    const std::string& text, int line_num, ComputationContext& context);
Relation parse_relation(
    const std::string& text, int line_num, bool conditional,
    ComputationContext& context);
class AssumptionParser {
public:
    explicit AssumptionParser(ComputationContext& context) : computation(context) {}
    AssumptionContext parse(const std::string& data);

private:
    ComputationContext& computation;
    AssumptionContext ctx;
    int line_num = 0;
    int current_scope = -1;
    bool ended = false;

    void parse_record(const std::string& line);
    void parse_scope(std::istringstream& ls, const std::string& keyword);
    void parse_domain(std::istringstream& ls, const std::string& keyword);
    void parse_sign(std::istringstream& ls, const std::string& keyword);
    void parse_parity(std::istringstream& ls, const std::string& keyword);
    void parse_bounded(std::istringstream& ls, const std::string& keyword);
    void parse_transcendental(std::istringstream& ls, const std::string& keyword);
    void parse_finiteness(std::istringstream& ls, const std::string& keyword);
    void parse_definiteness(std::istringstream& ls, const std::string& keyword);
    void parse_periodic(std::istringstream& ls, const std::string& keyword);
    void parse_relation(std::istringstream& ls, const std::string& keyword);
    void parse_conditional(std::istringstream& ls, const std::string& keyword);
    void parse_continuity(std::istringstream& ls, const std::string& keyword);
    void parse_monotonicity(std::istringstream& ls, const std::string& keyword);
};

}
