#include "lsr_expr.hpp"

#include <cctype>
#include <exception>
#include <utility>

#include "lsr_expr_internal.hpp"
#include "symbolic_ast.hpp"

namespace lamina::lsr {
namespace {

constexpr const char* kParseOperation = "lsr.parse_expr";

class ExprParser {
public:
    explicit ExprParser(std::string source) : source_(std::move(source)) {}

    ExprResult parse() {
        skip_space();
        if (eof()) {
            return fail("empty expression");
        }
        auto result = parse_logical_or();
        if (!result) return result;
        skip_space();
        if (!eof()) {
            return fail("unexpected token '" + std::string(1, peek()) + "'");
        }
        return result;
    }

private:
    std::string source_;
    std::size_t pos_ = 0;

    bool eof() const noexcept { return pos_ >= source_.size(); }

    char peek() const noexcept {
        return eof() ? '\0' : source_[pos_];
    }

    bool match(char c) {
        skip_space();
        if (peek() != c) return false;
        ++pos_;
        return true;
    }

    bool match_text(const char* text) {
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

    bool match_keyword(const char* text) {
        skip_space();
        const std::size_t start = pos_;
        if (!match_text(text)) return false;
        if (is_ident_continue(static_cast<unsigned char>(peek()))) {
            pos_ = start;
            return false;
        }
        return true;
    }

    void skip_space() noexcept {
        while (!eof() &&
               std::isspace(static_cast<unsigned char>(source_[pos_]))) {
            ++pos_;
        }
    }

    ExprResult fail(std::string message) const {
        return expression_failure(
            CasErrc::ParseError,
            std::move(message) + " at byte " + std::to_string(pos_),
            kParseOperation);
    }

    ExprResult unsupported(std::string message) const {
        return expression_failure(
            CasErrc::UnsupportedExpression,
            std::move(message) + " at byte " + std::to_string(pos_),
            kParseOperation);
    }

    static bool is_ident_start(unsigned char c) noexcept {
        return std::isalpha(c) || c == '_' || c >= 0x80;
    }

    static bool is_ident_continue(unsigned char c) noexcept {
        return std::isalnum(c) || c == '_' || c >= 0x80;
    }

    ExprResult parse_logical_or() {
        auto lhs = parse_logical_and();
        if (!lhs) return lhs;
        while (match_keyword("or")) {
            auto rhs = parse_logical_and();
            if (!rhs) return rhs;
            lhs = logical(lhs.value(), rhs.value(), LogicalNode::Op::Or);
            if (!lhs) return lhs;
        }
        return lhs;
    }

    ExprResult parse_logical_and() {
        auto lhs = parse_logical_not();
        if (!lhs) return lhs;
        while (match_keyword("and")) {
            auto rhs = parse_logical_not();
            if (!rhs) return rhs;
            lhs = logical(lhs.value(), rhs.value(), LogicalNode::Op::And);
            if (!lhs) return lhs;
        }
        return lhs;
    }

    ExprResult parse_logical_not() {
        if (match_keyword("not") || match('!')) {
            auto operand = parse_logical_not();
            if (!operand) return operand;
            return logical(operand.value(), nullptr, LogicalNode::Op::Not);
        }
        return parse_membership();
    }

    ExprResult parse_membership() {
        auto lhs = parse_equality();
        if (!lhs) return lhs;
        if (match_keyword("in")) {
            auto rhs = parse_equality();
            if (!rhs) return rhs;
            return membership(lhs.value(), rhs.value(), false);
        }
        const auto save = pos_;
        if (match_keyword("not")) {
            if (match_keyword("in")) {
                auto rhs = parse_equality();
                if (!rhs) return rhs;
                return membership(lhs.value(), rhs.value(), true);
            }
            pos_ = save;
        }
        return lhs;
    }

    ExprResult parse_equality() {
        auto lhs = parse_relational();
        if (!lhs) return lhs;
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
            if (!rhs) return rhs;
            lhs = relational(lhs.value(), rhs.value(), op);
            if (!lhs) return lhs;
        }
    }

    ExprResult parse_relational() {
        auto lhs = parse_additive();
        if (!lhs) return lhs;
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
            if (!rhs) return rhs;
            lhs = relational(lhs.value(), rhs.value(), op);
            if (!lhs) return lhs;
        }
    }

    ExprResult parse_additive() {
        auto lhs = parse_multiplicative();
        if (!lhs) return lhs;
        while (true) {
            if (match('+')) {
                auto rhs = parse_multiplicative();
                if (!rhs) return rhs;
                lhs = add(lhs.value(), rhs.value());
            } else if (match('-')) {
                auto rhs = parse_multiplicative();
                if (!rhs) return rhs;
                lhs = sub(lhs.value(), rhs.value());
            } else {
                return lhs;
            }
            if (!lhs) return lhs;
        }
    }

    ExprResult parse_multiplicative() {
        auto lhs = parse_power();
        if (!lhs) return lhs;
        while (true) {
            if (match_text("**")) {
                return fail("operator '**' is not supported; use '^'");
            } else if (match('*')) {
                auto rhs = parse_power();
                if (!rhs) return rhs;
                lhs = mul(lhs.value(), rhs.value());
            } else if (match('/')) {
                auto rhs = parse_power();
                if (!rhs) return rhs;
                lhs = div(lhs.value(), rhs.value());
            } else {
                return lhs;
            }
            if (!lhs) return lhs;
        }
    }

    ExprResult parse_power() {
        auto base = parse_unary();
        if (!base) return base;
        if (match_text("**")) {
            return fail("operator '**' is not supported; use '^'");
        }
        if (match('^')) {
            auto exponent = parse_power();
            if (!exponent) return exponent;
            return pow(base.value(), exponent.value());
        }
        return base;
    }

    ExprResult parse_unary() {
        if (match('+')) return parse_unary();
        if (match('-')) {
            auto operand = parse_unary();
            if (!operand) return operand;
            return neg(operand.value());
        }
        return parse_primary();
    }

    ExprResult parse_primary() {
        skip_space();
        if (match('(')) {
            auto inner = parse_logical_or();
            if (!inner) return inner;
            if (match(',')) {
                auto rhs = parse_logical_or();
                if (!rhs) return rhs;
                if (!match(')')) return fail("expected ')' after interval");
                return interval(inner.value(), rhs.value(), false, false);
            }
            if (!match(')')) return fail("expected ')'");
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

    ExprResult parse_set_literal() {
        std::vector<ExprPtr> elements;
        if (!match('}')) {
            while (true) {
                auto element = parse_logical_or();
                if (!element) return element;
                elements.push_back(element.value());
                if (match('}')) break;
                if (!match(',')) return fail("expected ',' or '}' in set literal");
            }
        }
        return finite_set(elements);
    }

    ExprResult parse_interval_literal(bool lower_closed) {
        auto lower = parse_logical_or();
        if (!lower) return lower;
        if (!match(',')) return fail("expected ',' in interval literal");
        auto upper = parse_logical_or();
        if (!upper) return upper;
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

    ExprResult parse_number() {
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
        if (!saw_digit) return fail("invalid number literal");
        if (peek() == 'e' || peek() == 'E') {
            ++pos_;
            if (peek() == '+' || peek() == '-') ++pos_;
            bool saw_exp_digit = false;
            while (std::isdigit(static_cast<unsigned char>(peek()))) {
                saw_exp_digit = true;
                ++pos_;
            }
            if (!saw_exp_digit) return fail("invalid number exponent");
        }
        try {
            return rational(Rational(source_.substr(start, pos_ - start)));
        } catch (const std::bad_alloc&) {
            return expression_failure(CasErrc::ResourceLimit,
                                      "number literal allocation failed",
                                      kParseOperation);
        } catch (const std::exception& error) {
            return expression_failure(CasErrc::ParseError, error.what(),
                                      kParseOperation);
        }
    }

    std::string parse_identifier() {
        skip_space();
        const std::size_t start = pos_;
        if (!is_ident_start(static_cast<unsigned char>(peek()))) return {};
        ++pos_;
        while (is_ident_continue(static_cast<unsigned char>(peek()))) ++pos_;
        return source_.substr(start, pos_ - start);
    }

    ExprResult parse_identifier_or_call() {
        auto name = parse_identifier();
        if (name.empty()) return fail("expected identifier");
        skip_space();
        if (match('(')) {
            std::vector<ExprPtr> arguments;
            bool closed = match(')');
            if (!closed) {
                do {
                    auto argument = parse_logical_or();
                    if (!argument) return argument;
                    arguments.push_back(argument.value());
                    if (match(')')) {
                        closed = true;
                        break;
                    }
                } while (match(','));
            }
            if (!closed) return fail("expected ')'");
            return apply_function(name, arguments);
        }
        if (name == "pi" || name == "π") return pi();
        if (name == "e") return e();
        if (name == "phi") return phi();
        if (is_imaginary_unit_name(name)) return imaginary_unit();
        return sym(name);
    }

    ExprResult relational(const ExprPtr& lhs, const ExprPtr& rhs, RelationOp op) {
        if (!lhs || !rhs) return expression_failure(CasErrc::InvalidArgument,
                                                    "relational operands cannot be null",
                                                    kParseOperation);
        try {
            auto node = lamina::detail::make_node<RelationalNode>(
                lamina::detail::node(lhs), lamina::detail::node(rhs), op);
            return ExprResult::success(lamina::detail::make_expression_ptr(std::move(node)));
        } catch (const std::bad_alloc&) {
            return expression_failure(CasErrc::ResourceLimit,
                                      "relational expression allocation failed",
                                      kParseOperation);
        } catch (const std::exception& error) {
            return expression_failure(CasErrc::InvalidArgument, error.what(),
                                      kParseOperation);
        }
    }

    ExprResult logical(const ExprPtr& lhs, const ExprPtr& rhs, LogicalNode::Op op) {
        if (!lhs || (op != LogicalNode::Op::Not && !rhs)) {
            return expression_failure(CasErrc::InvalidArgument,
                                      "logical operands cannot be null",
                                      kParseOperation);
        }
        try {
            auto node = lamina::detail::make_node<LogicalNode>(
                lamina::detail::node(lhs),
                rhs ? lamina::detail::node(rhs) : nullptr,
                op);
            return ExprResult::success(lamina::detail::make_expression_ptr(std::move(node)));
        } catch (const std::bad_alloc&) {
            return expression_failure(CasErrc::ResourceLimit,
                                      "logical expression allocation failed",
                                      kParseOperation);
        } catch (const std::exception& error) {
            return expression_failure(CasErrc::InvalidArgument, error.what(),
                                      kParseOperation);
        }
    }

    ExprResult finite_set(const std::vector<ExprPtr>& elements) {
        try {
            std::vector<std::shared_ptr<const SymbolicNode>> nodes;
            nodes.reserve(elements.size());
            for (const auto& element : elements) {
                if (!element || !lamina::detail::node(element)) {
                    return expression_failure(CasErrc::InvalidArgument,
                                              "set elements cannot be null",
                                              kParseOperation);
                }
                nodes.push_back(lamina::detail::node(element));
            }
            auto node = lamina::detail::make_node<FiniteSetNode>(std::move(nodes));
            return ExprResult::success(lamina::detail::make_expression_ptr(std::move(node)));
        } catch (const std::bad_alloc&) {
            return expression_failure(CasErrc::ResourceLimit,
                                      "set expression allocation failed",
                                      kParseOperation);
        } catch (const std::exception& error) {
            return expression_failure(CasErrc::InvalidArgument, error.what(),
                                      kParseOperation);
        }
    }

    ExprResult interval(const ExprPtr& lower,
                        const ExprPtr& upper,
                        bool lower_closed,
                        bool upper_closed) {
        if (!lower || !upper) {
            return expression_failure(CasErrc::InvalidArgument,
                                      "interval bounds cannot be null",
                                      kParseOperation);
        }
        try {
            auto node = lamina::detail::make_node<IntervalNode>(
                lamina::detail::node(lower), lamina::detail::node(upper),
                lower_closed, upper_closed);
            return ExprResult::success(lamina::detail::make_expression_ptr(std::move(node)));
        } catch (const std::bad_alloc&) {
            return expression_failure(CasErrc::ResourceLimit,
                                      "interval expression allocation failed",
                                      kParseOperation);
        } catch (const std::exception& error) {
            return expression_failure(CasErrc::InvalidArgument, error.what(),
                                      kParseOperation);
        }
    }

    ExprResult membership(const ExprPtr& element,
                          const ExprPtr& set,
                          bool negated) {
        if (!element || !set) {
            return expression_failure(CasErrc::InvalidArgument,
                                      "membership operands cannot be null",
                                      kParseOperation);
        }
        try {
            std::shared_ptr<const SymbolicNode> node =
                lamina::detail::make_node<MembershipNode>(
                    lamina::detail::node(element), lamina::detail::node(set));
            if (negated) {
                node = lamina::detail::make_node<LogicalNode>(
                    std::move(node), nullptr, LogicalNode::Op::Not);
            }
            return ExprResult::success(lamina::detail::make_expression_ptr(std::move(node)));
        } catch (const std::bad_alloc&) {
            return expression_failure(CasErrc::ResourceLimit,
                                      "membership expression allocation failed",
                                      kParseOperation);
        } catch (const std::exception& error) {
            return expression_failure(CasErrc::InvalidArgument, error.what(),
                                      kParseOperation);
        }
    }

    ExprResult function_node(const std::string& name,
                             const std::vector<ExprPtr>& arguments,
                             FunctionNode::FuncType type) {
        try {
            std::vector<std::shared_ptr<const SymbolicNode>> nodes;
            nodes.reserve(arguments.size());
            for (const auto& argument : arguments) {
                if (!argument || !lamina::detail::node(argument)) {
                    return expression_failure(CasErrc::InvalidArgument,
                                              "function arguments cannot be null",
                                              kParseOperation);
                }
                nodes.push_back(lamina::detail::node(argument));
            }
            auto node = lamina::detail::make_node<FunctionNode>(type, std::move(nodes));
            return ExprResult::success(lamina::detail::make_expression_ptr(std::move(node)));
        } catch (const std::bad_alloc&) {
            return expression_failure(CasErrc::ResourceLimit,
                                      name + " expression allocation failed",
                                      kParseOperation);
        } catch (const std::exception& error) {
            return expression_failure(CasErrc::InvalidArgument, error.what(),
                                      kParseOperation);
        }
    }

    ExprResult require_arity(const std::string& name,
                             const std::vector<ExprPtr>& arguments,
                             std::size_t arity,
                             FunctionNode::FuncType type) {
        if (arguments.size() != arity) {
            return fail("function '" + name + "' expects " + std::to_string(arity) + " argument(s)");
        }
        return function_node(name, arguments, type);
    }

    ExprResult apply_function(const std::string& name,
                              const std::vector<ExprPtr>& arguments) {
        if (name == "sin" && arguments.size() == 1) return sin(arguments[0]);
        if (name == "cos" && arguments.size() == 1) return cos(arguments[0]);
        if (name == "tan" && arguments.size() == 1) return tan(arguments[0]);
        if (name == "asin" && arguments.size() == 1) return asin(arguments[0]);
        if (name == "acos" && arguments.size() == 1) return acos(arguments[0]);
        if (name == "atan" && arguments.size() == 1) return atan(arguments[0]);
        if (name == "sqrt" && arguments.size() == 1) return sqrt(arguments[0]);
        if (name == "exp" && arguments.size() == 1) return exp(arguments[0]);
        if ((name == "ln" || name == "log") && arguments.size() == 1) return log(arguments[0]);
        if (name == "log" && arguments.size() == 2) {
            try {
                auto result = SymbolicExpr::log(arguments[0], arguments[1]);
                return ExprResult::success(std::move(result));
            } catch (const std::bad_alloc&) {
                return expression_failure(CasErrc::ResourceLimit,
                                          "log expression allocation failed",
                                          kParseOperation);
            } catch (const std::exception& error) {
                return expression_failure(CasErrc::InvalidArgument, error.what(),
                                          kParseOperation);
            }
        }
        if (name == "atan2" && arguments.size() == 2) {
            try {
                auto result = SymbolicExpr::atan2(arguments[0], arguments[1]);
                return ExprResult::success(std::move(result));
            } catch (const std::bad_alloc&) {
                return expression_failure(CasErrc::ResourceLimit,
                                          "atan2 expression allocation failed",
                                          kParseOperation);
            } catch (const std::exception& error) {
                return expression_failure(CasErrc::InvalidArgument, error.what(),
                                          kParseOperation);
            }
        }
        if (name == "log10" && arguments.size() == 1) return log10(arguments[0]);
        if (name == "floor" && arguments.size() == 1) return floor(arguments[0]);
        if (name == "ceil" && arguments.size() == 1) return ceil(arguments[0]);
        if (name == "round" && arguments.size() == 1) return round(arguments[0]);
        if (name == "abs" && arguments.size() == 1) return abs(arguments[0]);
        if (name == "real" && arguments.size() == 1) return real(arguments[0]);
        if (name == "imag" && arguments.size() == 1) return imag(arguments[0]);
        if (name == "conj" && arguments.size() == 1) return conj(arguments[0]);
        if (name == "max" && !arguments.empty()) {
            return function_node(name, arguments, FunctionNode::FuncType::Max);
        }
        if (name == "min" && !arguments.empty()) {
            return function_node(name, arguments, FunctionNode::FuncType::Min);
        }
        if (name == "clamp" && arguments.size() == 3) {
            return clamp(arguments[0], arguments[1], arguments[2]);
        }
        if (name == "sin" || name == "cos" || name == "tan" ||
            name == "asin" || name == "acos" || name == "atan" ||
            name == "sqrt" || name == "exp" || name == "ln" ||
            name == "log" || name == "log10" || name == "floor" ||
            name == "ceil" || name == "round" || name == "abs" ||
            name == "real" || name == "imag" || name == "conj" ||
            name == "atan2" || name == "clamp") {
            return fail("invalid argument count for function '" + name + "'");
        }
        try {
            std::vector<std::shared_ptr<const SymbolicNode>> nodes;
            nodes.reserve(arguments.size());
            for (const auto& argument : arguments) {
                if (!argument || !lamina::detail::node(argument)) {
                    return expression_failure(CasErrc::InvalidArgument,
                                              "function arguments cannot be null",
                                              kParseOperation);
                }
                nodes.push_back(lamina::detail::node(argument));
            }
            auto node = lamina::detail::make_node<UninterpretedFunctionNode>(
                name, std::move(nodes));
            return ExprResult::success(lamina::detail::make_expression_ptr(std::move(node)));
        } catch (const std::bad_alloc&) {
            return expression_failure(CasErrc::ResourceLimit,
                                      "function expression allocation failed",
                                      kParseOperation);
        } catch (const std::exception& error) {
            return expression_failure(CasErrc::InvalidArgument, error.what(),
                                      kParseOperation);
        }
    }
};

} // namespace

ExprResult parse_expr(const std::string& source) {
    try {
        return ExprParser(source).parse();
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "expression parse allocation failed",
                                  kParseOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::ParseError, error.what(),
                                  kParseOperation);
    }
}

} // namespace lamina::lsr
