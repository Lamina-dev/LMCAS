#pragma once

#include "expr.hpp"
#include "internal/symbolic_ast.hpp"

#include <string>
#include <utility>
#include <vector>

namespace LMCAS::expr_detail {

inline constexpr const char* kParseOperation = "LMCAS.parse_expr";

class ExprParser {
public:
    explicit ExprParser(std::string source) : source_(std::move(source)) {}
    ExprResult parse();

private:
    std::string source_;
    std::size_t pos_ = 0;

    bool eof() const noexcept;
    char peek() const noexcept;
    bool match(char c);
    bool match_text(const char* text);
    bool match_keyword(const char* text);
    void skip_space() noexcept;
    ExprResult fail(std::string message) const;
    static bool is_ident_start(unsigned char c) noexcept;
    static bool is_ident_continue(unsigned char c) noexcept;
    ExprResult parse_logical_or();
    ExprResult parse_logical_and();
    ExprResult parse_logical_not();
    ExprResult parse_membership();
    ExprResult parse_equality();
    ExprResult parse_relational();
    ExprResult parse_additive();
    ExprResult parse_multiplicative();
    ExprResult parse_power();
    ExprResult parse_unary();
    ExprResult parse_primary();
    ExprResult parse_set_literal();
    ExprResult parse_interval_literal(bool lower_closed);
    ExprResult parse_number();
    std::string parse_identifier();
    ExprResult parse_identifier_or_call();
    ExprResult parse_call(const std::string& name);
    ExprResult parse_approx_literal();

    ExprResult relational(const ExprPtr& lhs, const ExprPtr& rhs, RelationOp op);
    ExprResult logical(const ExprPtr& lhs, const ExprPtr& rhs, LogicalNode::Op op);
    ExprResult finite_set(const std::vector<ExprPtr>& elements);
    ExprResult interval(const ExprPtr& lower, const ExprPtr& upper,
                        bool lower_closed, bool upper_closed);
    ExprResult membership(const ExprPtr& element, const ExprPtr& set, bool negated);
    ExprResult function_node(const std::string& name,
                             const std::vector<ExprPtr>& arguments,
                             FunctionNode::FuncType type);
    ExprResult apply_function(const std::string& name,
                              const std::vector<ExprPtr>& arguments);
    ExprResult apply_log(const std::vector<ExprPtr>& arguments);
    ExprResult apply_atan2(const std::vector<ExprPtr>& arguments);
    ExprResult apply_integral(const std::vector<ExprPtr>& arguments);
    ExprResult apply_limit(const std::vector<ExprPtr>& arguments);
    ExprResult apply_rootof(const std::vector<ExprPtr>& arguments);
    ExprResult apply_clamp(const std::vector<ExprPtr>& arguments);
    ExprResult uninterpreted_function(const std::string& name,
                                      const std::vector<ExprPtr>& arguments);
};

}
