#include "lsr_expr.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <unordered_set>
#include <utility>

#include "assumption_context.hpp"
#include "complex_analysis.hpp"
#include "matcher.hpp"
#include "poly_utils.hpp"
#include "symbolic_ast.hpp"

namespace lamina::lsr {
namespace {

constexpr const char* kSymOperation = "lsr.sym";
constexpr const char* kIntegerOperation = "lsr.integer";
constexpr const char* kRationalOperation = "lsr.rational";
constexpr const char* kApproxOperation = "lsr.approx_real";
constexpr const char* kConstantOperation = "lsr.constant";
constexpr const char* kComplexOperation = "lsr.complex";
constexpr const char* kParseOperation = "lsr.parse_expr";
constexpr const char* kExprOperation = "lsr.expr_op";
constexpr const char* kMathOperation = "lsr.math";
constexpr const char* kRealOperation = "lsr.real";
constexpr const char* kImagOperation = "lsr.imag";
constexpr const char* kConjOperation = "lsr.conj";
constexpr const char* kAbsOperation = "lsr.abs";
constexpr const char* kSimplifyOperation = "lsr.simplify";
constexpr const char* kExpandOperation = "lsr.expand";
constexpr const char* kDifferentiateOperation = "lsr.differentiate";
constexpr const char* kSubstituteOperation = "lsr.substitute";
constexpr const char* kExprMatchOperation = "lsr.expr_match";
constexpr const char* kEquivalentOperation = "lsr.equivalent_core";
constexpr const char* kEquivalentProfileOperation = "lsr.equivalent_core.profile";
constexpr const char* kExprSetOperation = "lsr.expr_set";
constexpr const char* kNumberDomainOperation = "lsr.number_domain";
constexpr const char* kSolveExprSetOperation = "lsr.solve_expr_set";
constexpr const char* kEvalfOperation = "lsr.evalf";
constexpr const char* kEvalComplexOperation = "lsr.eval_complex";

ExprResult expression_failure(CasErrc code, std::string message,
                              const char* operation) {
    return ExprResult::failure(code, std::move(message), operation);
}

ExprSetResult expr_set_failure(CasErrc code, std::string message,
                               const char* operation) {
    return ExprSetResult::failure(code, std::move(message), operation);
}

Result<bool> bool_failure(CasErrc code, std::string message,
                          const char* operation) {
    return Result<bool>::failure(code, std::move(message), operation);
}

Result<EqvProfile> eqv_profile_failure(std::string message) {
    return Result<EqvProfile>::failure(
        CasErrc::UnsupportedExpression,
        std::move(message),
        kEquivalentProfileOperation);
}

Result<ApproxComplex> complex_failure(CasErrc code, std::string message,
                                      const char* operation) {
    return Result<ApproxComplex>::failure(code, std::move(message), operation);
}

Result<ApproxComplex> eval_complex_failure(const CasError& error) {
    return complex_failure(error.code, error.message, kEvalComplexOperation);
}

ApproxReal approx_part(double value) {
    ApproxReal part;
    part.value = value;
    part.absolute_error = std::numeric_limits<double>::epsilon() *
                          std::max(1.0, std::abs(value)) * 4.0;
    part.status = NumericStatus::Finite;
    return part;
}

ApproxComplex approx_complex(double real, double imag) {
    return ApproxComplex{approx_part(real), approx_part(imag)};
}

Result<ApproxComplex> checked_complex(double real, double imag,
                                      const char* operation) {
    if (!std::isfinite(real) || !std::isfinite(imag)) {
        return complex_failure(CasErrc::NumericFailure,
                               "complex evaluation produced a non-finite component",
                               operation);
    }
    return Result<ApproxComplex>::success(approx_complex(real, imag));
}

Result<ApproxComplex> real_to_complex(const Result<ApproxReal>& real) {
    if (!real) {
        return eval_complex_failure(real.error());
    }
    if (!real.value().is_finite()) {
        return complex_failure(CasErrc::NumericFailure,
                               "complex evaluation requires finite real components",
                               kEvalComplexOperation);
    }
    return Result<ApproxComplex>::success(
        ApproxComplex{real.value(), approx_part(0.0)});
}

Result<ApproxComplex> add_complex(const ApproxComplex& lhs,
                                  const ApproxComplex& rhs) {
    auto result = checked_complex(lhs.real.value + rhs.real.value,
                                  lhs.imag.value + rhs.imag.value,
                                  kEvalComplexOperation);
    if (!result) return result;
    result.value().real.absolute_error +=
        lhs.real.absolute_error + rhs.real.absolute_error;
    result.value().imag.absolute_error +=
        lhs.imag.absolute_error + rhs.imag.absolute_error;
    return result;
}

Result<ApproxComplex> multiply_complex(const ApproxComplex& lhs,
                                       const ApproxComplex& rhs) {
    const double a = lhs.real.value;
    const double b = lhs.imag.value;
    const double c = rhs.real.value;
    const double d = rhs.imag.value;
    auto result = checked_complex(a * c - b * d, a * d + b * c,
                                  kEvalComplexOperation);
    if (!result) return result;
    result.value().real.absolute_error +=
        std::abs(c) * lhs.real.absolute_error +
        std::abs(a) * rhs.real.absolute_error +
        std::abs(d) * lhs.imag.absolute_error +
        std::abs(b) * rhs.imag.absolute_error;
    result.value().imag.absolute_error +=
        std::abs(d) * lhs.real.absolute_error +
        std::abs(a) * rhs.imag.absolute_error +
        std::abs(c) * lhs.imag.absolute_error +
        std::abs(b) * rhs.real.absolute_error;
    return result;
}

Result<ApproxComplex> divide_complex(const ApproxComplex& lhs,
                                     const ApproxComplex& rhs) {
    const double c = rhs.real.value;
    const double d = rhs.imag.value;
    const double denom = c * c + d * d;
    if (denom == 0.0 || !std::isfinite(denom)) {
        return complex_failure(CasErrc::DomainError,
                               "complex division by zero or overflow",
                               kEvalComplexOperation);
    }
    return checked_complex((lhs.real.value * c + lhs.imag.value * d) / denom,
                           (lhs.imag.value * c - lhs.real.value * d) / denom,
                           kEvalComplexOperation);
}

bool is_integer_double(double value) {
    return std::isfinite(value) && std::floor(value) == value;
}

bool is_imaginary_unit_name(const std::string& name) {
    return name == "I";
}

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

int domain_rank(NumberDomain domain) noexcept {
    switch (domain) {
    case NumberDomain::Integers:
        return 0;
    case NumberDomain::Rationals:
        return 1;
    case NumberDomain::Reals:
        return 2;
    case NumberDomain::Complexes:
        return 3;
    case NumberDomain::Expressions:
        return 4;
    }
    return -1;
}

std::optional<int> exact_small_integer_node(
    const std::shared_ptr<const SymbolicNode>& node,
    int min_value,
    int max_value);

Result<bool> domain_contains_node(
    NumberDomain domain,
    const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) {
        return bool_failure(CasErrc::InvalidArgument,
                            "domain membership element must not be null",
                            kNumberDomainOperation);
    }
    if (domain == NumberDomain::Expressions) {
        return Result<bool>::success(true);
    }

    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (is_imaginary_unit_name(variable->name())) {
            return Result<bool>::success(domain == NumberDomain::Complexes);
        }
    }

    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        const auto& value = number->value();
        if (std::holds_alternative<BigInt>(value)) {
            return Result<bool>::success(true);
        }
        if (std::holds_alternative<Rational>(value)) {
            const auto& rational = std::get<Rational>(value);
            switch (domain) {
            case NumberDomain::Integers:
                return Result<bool>::success(rational.is_integer());
            case NumberDomain::Rationals:
            case NumberDomain::Reals:
            case NumberDomain::Complexes:
            case NumberDomain::Expressions:
                return Result<bool>::success(true);
            }
        }
        const lmmc_real_t real = std::get<lmmc_real_t>(value);
        if (!std::isfinite(static_cast<double>(real))) {
            return bool_failure(CasErrc::NumericFailure,
                                "domain membership requires finite numeric literals",
                                kNumberDomainOperation);
        }
        switch (domain) {
        case NumberDomain::Integers:
        case NumberDomain::Rationals:
            return Result<bool>::success(false);
        case NumberDomain::Reals:
        case NumberDomain::Complexes:
        case NumberDomain::Expressions:
            return Result<bool>::success(true);
        }
    }

    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        if (domain == NumberDomain::Complexes) {
            auto real_part =
                domain_contains_node(NumberDomain::Reals, complex_node->real());
            if (!real_part) return real_part;
            auto imag_part =
                domain_contains_node(NumberDomain::Reals, complex_node->imag());
            if (!imag_part) return imag_part;
            return Result<bool>::success(real_part.value() && imag_part.value());
        }

        if (!complex_node->imag()->is_zero()) {
            return Result<bool>::success(false);
        }
        return domain_contains_node(domain, complex_node->real());
    }

    if (domain == NumberDomain::Complexes) {
        if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
            for (const auto& operand : add->operands()) {
                auto member = domain_contains_node(domain, operand);
                if (!member || !member.value()) return member;
            }
            return Result<bool>::success(true);
        }

        if (auto multiply =
                std::dynamic_pointer_cast<const MultiplyNode>(node)) {
            for (const auto& operand : multiply->operands()) {
                auto member = domain_contains_node(domain, operand);
                if (!member || !member.value()) return member;
            }
            return Result<bool>::success(true);
        }

        if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
            if (!exact_small_integer_node(power->exponent(), 1, 64)) {
                return bool_failure(
                    CasErrc::Inconclusive,
                    "complex domain membership for powers requires a positive exact integer exponent",
                    kNumberDomainOperation);
            }
            auto base_member = domain_contains_node(domain, power->base());
            if (!base_member) return base_member;
            return Result<bool>::success(base_member.value());
        }
    }

    return bool_failure(CasErrc::Inconclusive,
                        "domain membership is only decidable for numeric and explicit complex expressions",
                        kNumberDomainOperation);
}

Result<void> validate_eqv_options(const EqvOptions& options) {
    if (options.budget.max_rewrite_steps == 0 ||
        options.budget.max_rewrite_depth == 0 ||
        options.budget.max_node_growth_factor == 0) {
        return Result<void>::failure(
            CasErrc::ResourceLimit,
            "equivalence rewrite budget exhausted before normalization",
            kEquivalentOperation);
    }
    return Result<void>::success();
}

bool exact_integer_node(const std::shared_ptr<const SymbolicNode>& node,
                        int expected) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return false;
    if (std::holds_alternative<BigInt>(number->value())) {
        return std::get<BigInt>(number->value()) == BigInt(expected);
    }
    if (std::holds_alternative<Rational>(number->value())) {
        return std::get<Rational>(number->value()) == Rational(expected);
    }
    return false;
}

std::optional<int> exact_small_integer_node(
    const std::shared_ptr<const SymbolicNode>& node,
    int min_value,
    int max_value) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return std::nullopt;

    BigInt value;
    if (std::holds_alternative<BigInt>(number->value())) {
        value = std::get<BigInt>(number->value());
    } else if (std::holds_alternative<Rational>(number->value())) {
        const Rational& rational = std::get<Rational>(number->value());
        if (!rational.is_integer()) return std::nullopt;
        value = rational.to_BigInt();
    } else {
        return std::nullopt;
    }

    if (value < BigInt(min_value) || value > BigInt(max_value)) {
        return std::nullopt;
    }
    return value.to_int();
}

bool trig_square_argument(const std::shared_ptr<const SymbolicNode>& node,
                          FunctionNode::FuncType type,
                          std::shared_ptr<const SymbolicNode>& argument) {
    auto power = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!power || !exact_integer_node(power->exponent(), 2)) return false;
    auto function = std::dynamic_pointer_cast<const FunctionNode>(power->base());
    if (!function || function->type() != type ||
        function->arguments().size() != 1) {
        return false;
    }
    argument = function->arguments()[0];
    return true;
}

std::optional<ExprPtr> unwrap_trig_negated_argument(
    const std::shared_ptr<const SymbolicNode>& node) {
    auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!multiply) return std::nullopt;

    bool found_negative_one = false;
    std::vector<std::shared_ptr<const SymbolicNode>> remaining;
    remaining.reserve(multiply->operands().size());
    for (const auto& operand : multiply->operands()) {
        if (!found_negative_one && exact_integer_node(operand, -1)) {
            found_negative_one = true;
            continue;
        }
        remaining.push_back(operand);
    }
    if (!found_negative_one || remaining.empty()) return std::nullopt;
    if (remaining.size() == 1) {
        return lamina::detail::make_expression_ptr(remaining.front());
    }
    return lamina::detail::make_expression_ptr(
        SymbolicFactory::create_multiply(std::move(remaining)))->simplify();
}

ExprPtr rewrite_trig_basic_identity(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return nullptr;

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> rewritten;
        rewritten.reserve(add->operands().size());
        for (const auto& operand : add->operands()) {
            auto child = rewrite_trig_basic_identity(operand);
            rewritten.push_back(lamina::detail::node(child));
        }

        std::vector<bool> used(rewritten.size(), false);
        std::vector<std::shared_ptr<const SymbolicNode>> result_nodes;
        for (std::size_t i = 0; i < rewritten.size(); ++i) {
            if (used[i]) continue;
            std::shared_ptr<const SymbolicNode> sin_arg;
            std::shared_ptr<const SymbolicNode> cos_arg;
            const bool is_sin_square = trig_square_argument(
                rewritten[i], FunctionNode::FuncType::Sin, sin_arg);
            const bool is_cos_square = trig_square_argument(
                rewritten[i], FunctionNode::FuncType::Cos, cos_arg);
            bool matched = false;
            for (std::size_t j = i + 1; j < rewritten.size(); ++j) {
                if (used[j]) continue;
                std::shared_ptr<const SymbolicNode> other_arg;
                if (is_sin_square &&
                    trig_square_argument(rewritten[j],
                                         FunctionNode::FuncType::Cos,
                                         other_arg) &&
                    sin_arg->equals(*other_arg)) {
                    matched = true;
                } else if (is_cos_square &&
                           trig_square_argument(rewritten[j],
                                                FunctionNode::FuncType::Sin,
                                                other_arg) &&
                           cos_arg->equals(*other_arg)) {
                    matched = true;
                }
                if (matched) {
                    used[i] = true;
                    used[j] = true;
                    result_nodes.push_back(
                        lamina::detail::node(SymbolicExpr::number(1)));
                    break;
                }
            }
            if (!matched && !used[i]) {
                result_nodes.push_back(rewritten[i]);
            }
        }

        if (result_nodes.empty()) {
            return SymbolicExpr::number(0);
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_add(std::move(result_nodes)))->simplify();
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        operands.reserve(multiply->operands().size());
        for (const auto& operand : multiply->operands()) {
            auto child = rewrite_trig_basic_identity(operand);
            operands.push_back(lamina::detail::node(child));
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_multiply(std::move(operands)))->simplify();
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base = rewrite_trig_basic_identity(power->base());
        auto exponent = rewrite_trig_basic_identity(power->exponent());
        return SymbolicExpr::power(base, exponent)->simplify();
    }

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> args;
        args.reserve(function->arguments().size());
        for (const auto& argument : function->arguments()) {
            auto child = rewrite_trig_basic_identity(argument);
            args.push_back(lamina::detail::node(child));
        }

        if (args.size() == 1) {
            auto positive_arg = unwrap_trig_negated_argument(args[0]);
            if (positive_arg && *positive_arg) {
                if (function->type() == FunctionNode::FuncType::Sin) {
                    return SymbolicExpr::multiply(
                        SymbolicExpr::number(-1),
                        SymbolicExpr::sin(*positive_arg))->simplify();
                }
                if (function->type() == FunctionNode::FuncType::Cos) {
                    return SymbolicExpr::cos(*positive_arg)->simplify();
                }
            }
        }

        return lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(function->type(),
                                                    std::move(args)))->simplify();
    }

    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real_part = rewrite_trig_basic_identity(complex_node->real());
        auto imag_part = rewrite_trig_basic_identity(complex_node->imag());
        auto value = complex(real_part, imag_part);
        if (!value) throw std::runtime_error(value.error().message);
        return value.value()->simplify();
    }

    return lamina::detail::make_expression_ptr(node);
}

ExprPtr rewrite_exp_log_basic_identity(
    const std::shared_ptr<const SymbolicNode>& node,
    const AssumptionContext* assumptions) {
    if (!node) return nullptr;

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        operands.reserve(add->operands().size());
        for (const auto& operand : add->operands()) {
            auto child = rewrite_exp_log_basic_identity(operand, assumptions);
            operands.push_back(lamina::detail::node(child));
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_add(std::move(operands)))->simplify();
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        operands.reserve(multiply->operands().size());
        for (const auto& operand : multiply->operands()) {
            auto child = rewrite_exp_log_basic_identity(operand, assumptions);
            operands.push_back(lamina::detail::node(child));
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_multiply(std::move(operands)))->simplify();
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base = rewrite_exp_log_basic_identity(power->base(), assumptions);
        auto exponent = rewrite_exp_log_basic_identity(power->exponent(), assumptions);
        return SymbolicExpr::power(base, exponent)->simplify();
    }

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> args;
        args.reserve(function->arguments().size());
        for (const auto& argument : function->arguments()) {
            auto child = rewrite_exp_log_basic_identity(argument, assumptions);
            args.push_back(lamina::detail::node(child));
        }

        if (function->type() == FunctionNode::FuncType::Exp &&
            args.size() == 1 && exact_integer_node(args[0], 0)) {
            return SymbolicExpr::number(1);
        }
        if (function->type() == FunctionNode::FuncType::Ln &&
            args.size() == 1 && exact_integer_node(args[0], 1)) {
            return SymbolicExpr::number(0);
        }
        if (function->type() == FunctionNode::FuncType::Exp &&
            args.size() == 1) {
            auto inner_ln = std::dynamic_pointer_cast<const FunctionNode>(args[0]);
            if (inner_ln && inner_ln->type() == FunctionNode::FuncType::Ln &&
                inner_ln->arguments().size() == 1) {
                auto ln_arg = lamina::detail::make_expression_ptr(
                    inner_ln->arguments()[0]);
                const bool known_positive = inner_ln->arguments()[0]->is_positive() ||
                    (assumptions &&
                     assumptions->is_positive(*ln_arg) == Tribool::True);
                if (known_positive) {
                    return ln_arg->simplify();
                }
            }
        }

        return lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(function->type(),
                                                    std::move(args)))->simplify();
    }

    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real_part = rewrite_exp_log_basic_identity(complex_node->real(),
                                                        assumptions);
        auto imag_part = rewrite_exp_log_basic_identity(complex_node->imag(),
                                                        assumptions);
        auto value = complex(real_part, imag_part);
        if (!value) throw std::runtime_error(value.error().message);
        return value.value()->simplify();
    }

    return lamina::detail::make_expression_ptr(node);
}

void collect_variable_names(
    const std::shared_ptr<const SymbolicNode>& node,
    std::set<std::string>& variables) {
    if (!node) return;
    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        variables.insert(variable->name());
        return;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& operand : add->operands()) {
            collect_variable_names(operand, variables);
        }
        return;
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& operand : multiply->operands()) {
            collect_variable_names(operand, variables);
        }
        return;
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        collect_variable_names(power->base(), variables);
        collect_variable_names(power->exponent(), variables);
        return;
    }
    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        collect_variable_names(complex_node->real(), variables);
        collect_variable_names(complex_node->imag(), variables);
    }
}

Result<std::optional<bool>> prove_rational_polynomial_equivalence(
    const ExprPtr& difference,
    ComputationContext& context,
    const EqvOptions& options) {
    if (!difference) {
        return Result<std::optional<bool>>::failure(
            CasErrc::InternalInvariant,
            "equivalence difference is null",
            kEquivalentOperation);
    }

    if (options.budget.max_rewrite_steps < 4) {
        return Result<std::optional<bool>>::failure(
            CasErrc::ResourceLimit,
            "equivalence rewrite budget exhausted before polynomial normalization",
            kEquivalentOperation);
    }

    auto step = context.consume_steps(4, kEquivalentOperation);
    if (!step) return Result<std::optional<bool>>::failure(step.error());

    auto expanded = difference->expand();
    if (!expanded || !lamina::detail::node(expanded)) {
        return Result<std::optional<bool>>::failure(
            CasErrc::InternalInvariant,
            "equivalence expansion returned null",
            kEquivalentOperation);
    }

    std::set<std::string> variables;
    collect_variable_names(lamina::detail::node(expanded), variables);
    if (variables.size() > 1) {
        return Result<std::optional<bool>>::success(std::nullopt);
    }
    const std::string variable = variables.empty() ? "x" : *variables.begin();

    auto recognized = recognize_rational_polynomial(*expanded, variable, context);
    if (!recognized) {
        return Result<std::optional<bool>>::failure(recognized.error());
    }
    if (!recognized.value()) {
        return Result<std::optional<bool>>::success(std::nullopt);
    }
    return Result<std::optional<bool>>::success(recognized.value()->is_zero());
}

Rational polynomial_coeff_or_zero(const Polynomial<Rational>& polynomial,
                                  std::size_t degree) {
    return degree < polynomial.coeffs.size()
        ? polynomial.coeffs[degree]
        : Rational(0);
}

ExprPtr rational_expression(const Rational& value) {
    return SymbolicExpr::number(value);
}

bool exact_rational_sqrt(const Rational& value, Rational& root) {
    if (value < Rational(0)) return false;
    const BigInt numerator_root = value.get_numerator().sqrt();
    const BigInt denominator_root = value.get_denominator().sqrt();
    if (numerator_root * numerator_root != value.get_numerator() ||
        denominator_root * denominator_root != value.get_denominator()) {
        return false;
    }
    root = Rational(numerator_root, denominator_root);
    return true;
}

ExprPtr sqrt_rational_expression(const Rational& value) {
    Rational root;
    if (exact_rational_sqrt(value, root)) {
        return rational_expression(root);
    }
    return SymbolicExpr::sqrt(rational_expression(value))->simplify();
}

ExprResult verified_lsr_complex(ExprPtr real_part, ExprPtr imag_part) {
    auto value = complex(std::move(real_part), std::move(imag_part));
    if (!value) return value;
    auto simplified = value.value()->simplify();
    if (!simplified || !lamina::detail::node(simplified)) {
        return expression_failure(CasErrc::InternalInvariant,
                                  "complex root simplification returned null",
                                  kSolveExprSetOperation);
    }
    return ExprResult::success(std::move(simplified));
}

Result<std::optional<ExprSet>> try_lsr_closed_form_rational_poly_roots(
    const ExprPtr& equation,
    const std::string& variable,
    ComputationContext& context) {
    if (!equation) {
        return Result<std::optional<ExprSet>>::failure(
            CasErrc::InvalidArgument, "equation cannot be null",
            kSolveExprSetOperation);
    }
    if (variable.empty()) {
        return Result<std::optional<ExprSet>>::failure(
            CasErrc::InvalidArgument, "solve variable cannot be empty",
            kSolveExprSetOperation);
    }

    auto recognized = recognize_rational_polynomial(*equation, variable, context);
    if (!recognized) {
        return Result<std::optional<ExprSet>>::failure(recognized.error());
    }
    if (!recognized.value()) {
        return Result<std::optional<ExprSet>>::success(std::nullopt);
    }

    const Polynomial<Rational>& polynomial = *recognized.value();
    const int degree = polynomial.degree();
    if (degree < 1 || degree > 2) {
        return Result<std::optional<ExprSet>>::success(std::nullopt);
    }

    std::vector<ExprPtr> roots;
    if (degree == 1) {
        const Rational b = polynomial_coeff_or_zero(polynomial, 1);
        if (b.is_zero()) {
            return Result<std::optional<ExprSet>>::success(std::nullopt);
        }
        const Rational c = polynomial_coeff_or_zero(polynomial, 0);
        roots.push_back(rational_expression((-c) / b)->simplify());
    } else {
        const Rational a = polynomial_coeff_or_zero(polynomial, 2);
        const Rational b = polynomial_coeff_or_zero(polynomial, 1);
        const Rational c = polynomial_coeff_or_zero(polynomial, 0);
        if (a.is_zero()) {
            return Result<std::optional<ExprSet>>::success(std::nullopt);
        }

        const Rational two_a = Rational(2) * a;
        const Rational discriminant = b * b - Rational(4) * a * c;
        const Rational real_component = (-b) / two_a;

        if (discriminant < Rational(0)) {
            const Rational positive_discriminant = -discriminant;
            const Rational positive_denominator = two_a.abs();
            auto imag_magnitude = SymbolicExpr::divide(
                sqrt_rational_expression(positive_discriminant),
                rational_expression(positive_denominator))->simplify();
            auto positive = verified_lsr_complex(
                rational_expression(real_component), imag_magnitude);
            if (!positive) {
                return Result<std::optional<ExprSet>>::failure(positive.error());
            }
            auto negative_imag = SymbolicExpr::multiply(
                SymbolicExpr::number(-1), imag_magnitude)->simplify();
            auto negative = verified_lsr_complex(
                rational_expression(real_component), negative_imag);
            if (!negative) {
                return Result<std::optional<ExprSet>>::failure(negative.error());
            }
            roots.push_back(negative.value());
            roots.push_back(positive.value());
        } else {
            auto sqrt_discriminant = sqrt_rational_expression(discriminant);
            auto numerator_left = SymbolicExpr::add(
                rational_expression(-b),
                SymbolicExpr::multiply(SymbolicExpr::number(-1),
                                       sqrt_discriminant))->simplify();
            roots.push_back(SymbolicExpr::divide(
                numerator_left, rational_expression(two_a))->simplify());
            if (!discriminant.is_zero()) {
                auto numerator_right = SymbolicExpr::add(
                    rational_expression(-b), sqrt_discriminant)->simplify();
                roots.push_back(SymbolicExpr::divide(
                    numerator_right, rational_expression(two_a))->simplify());
            }
        }
    }

    auto set = ExprSet::make(std::move(roots));
    if (!set) {
        return Result<std::optional<ExprSet>>::failure(set.error());
    }
    return Result<std::optional<ExprSet>>::success(std::move(set.value()));
}

ExprPtr canonicalize_lsr_complex_product(const SymbolicExpr& expression) {
    const auto& node = lamina::detail::node(expression);

    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (is_imaginary_unit_name(variable->name())) {
            auto unit = imaginary_unit();
            if (!unit) throw std::runtime_error(unit.error().message);
            return unit.value();
        }
        return lamina::detail::make_expression_ptr(node);
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        operands.reserve(add->operands().size());
        for (const auto& operand : add->operands()) {
            auto canonical_operand = canonicalize_lsr_complex_product(
                *lamina::detail::make_expression_ptr(operand));
            operands.push_back(lamina::detail::node(canonical_operand));
        }
        return lamina::detail::make_expression_ptr(
            SymbolicFactory::create_add(std::move(operands)))->simplify();
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent = exact_small_integer_node(power->exponent(), 0, 16);
        if (exponent) {
            auto canonical_base = canonicalize_lsr_complex_product(
                *lamina::detail::make_expression_ptr(power->base()));
            if (std::dynamic_pointer_cast<const ComplexNode>(
                    lamina::detail::node(canonical_base))) {
                auto result = SymbolicExpr::number(1);
                for (int i = 0; i < *exponent; ++i) {
                    result = canonicalize_lsr_complex_product(
                        *SymbolicExpr::multiply(result, canonical_base));
                }
                return result->simplify();
            }
        }
        return lamina::detail::make_expression_ptr(node);
    }

    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real_part = canonicalize_lsr_complex_product(
            *lamina::detail::make_expression_ptr(complex_node->real()));
        auto imag_part = canonicalize_lsr_complex_product(
            *lamina::detail::make_expression_ptr(complex_node->imag()));
        auto value = complex(real_part, imag_part);
        if (!value) throw std::runtime_error(value.error().message);
        return value.value()->simplify();
    }

    auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!multiply) {
        return lamina::detail::make_expression_ptr(node);
    }

    bool saw_complex = false;
    auto real = SymbolicExpr::number(1);
    auto imag = SymbolicExpr::number(0);

    for (const auto& operand : multiply->operands()) {
        auto canonical_operand = canonicalize_lsr_complex_product(
            *lamina::detail::make_expression_ptr(operand));
        const auto& operand_node = lamina::detail::node(canonical_operand);
        ExprPtr factor_real;
        ExprPtr factor_imag;
        if (auto complex_operand = std::dynamic_pointer_cast<const ComplexNode>(
                operand_node)) {
            saw_complex = true;
            factor_real = lamina::detail::make_expression_ptr(complex_operand->real());
            factor_imag = lamina::detail::make_expression_ptr(complex_operand->imag());
        } else {
            factor_real = canonical_operand;
            factor_imag = SymbolicExpr::number(0);
        }

        auto ac = SymbolicExpr::multiply(real, factor_real);
        auto bd = SymbolicExpr::multiply(imag, factor_imag);
        auto ad = SymbolicExpr::multiply(real, factor_imag);
        auto bc = SymbolicExpr::multiply(imag, factor_real);
        auto next_real = SymbolicExpr::add(
            ac, SymbolicExpr::multiply(SymbolicExpr::number(-1), bd))->simplify();
        auto next_imag = SymbolicExpr::add(ad, bc)->simplify();
        real = std::move(next_real);
        imag = std::move(next_imag);
    }

    if (!saw_complex) {
        return lamina::detail::make_expression_ptr(node);
    }
    auto result = complex(real, imag);
    if (!result) {
        throw std::runtime_error(result.error().message);
    }
    return result.value()->simplify();
}

Result<ApproxComplex> evaluate_complex_node(
    const std::shared_ptr<const SymbolicNode>& node,
    const NumericBindings& bindings,
    ComputationContext& context) {
    auto step = context.consume_steps(1, kEvalComplexOperation);
    if (!step) return Result<ApproxComplex>::failure(step.error());
    if (!node) {
        return complex_failure(CasErrc::InvalidArgument,
                               "expression contains a null node",
                               kEvalComplexOperation);
    }

    if (auto complex_node = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real = evaluate_numeric(
            *lamina::detail::make_expression_ptr(complex_node->real()),
            bindings, context);
        if (!real) return eval_complex_failure(real.error());
        auto imag = evaluate_numeric(
            *lamina::detail::make_expression_ptr(complex_node->imag()),
            bindings, context);
        if (!imag) return eval_complex_failure(imag.error());
        if (!real.value().is_finite() || !imag.value().is_finite()) {
            return complex_failure(CasErrc::NumericFailure,
                                   "complex components must be finite",
                                   kEvalComplexOperation);
        }
        return Result<ApproxComplex>::success(
            ApproxComplex{real.value(), imag.value()});
    }

    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (is_imaginary_unit_name(variable->name())) {
            return Result<ApproxComplex>::success(approx_complex(0.0, 1.0));
        }
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        auto sum = Result<ApproxComplex>::success(approx_complex(0.0, 0.0));
        for (const auto& operand : add->operands()) {
            auto term = evaluate_complex_node(operand, bindings, context);
            if (!term) return term;
            sum = add_complex(sum.value(), term.value());
            if (!sum) return sum;
        }
        return sum;
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        auto product = Result<ApproxComplex>::success(approx_complex(1.0, 0.0));
        for (const auto& operand : multiply->operands()) {
            auto factor = evaluate_complex_node(operand, bindings, context);
            if (!factor) return factor;
            product = multiply_complex(product.value(), factor.value());
            if (!product) return product;
        }
        return product;
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base = evaluate_complex_node(power->base(), bindings, context);
        if (!base) return base;
        auto exponent = evaluate_numeric(
            *lamina::detail::make_expression_ptr(power->exponent()),
            bindings, context);
        if (!exponent) return eval_complex_failure(exponent.error());
        if (!exponent.value().is_finite() ||
            !std::isfinite(exponent.value().value)) {
            return complex_failure(CasErrc::NumericFailure,
                                   "complex power exponent must be finite",
                                   kEvalComplexOperation);
        }
        const double exponent_value = exponent.value().value;
        if (!is_integer_double(exponent_value) ||
            std::abs(exponent_value) > 64.0) {
            return complex_failure(CasErrc::UnsupportedExpression,
                                   "complex evaluation only supports integer powers with |n| <= 64",
                                   kEvalComplexOperation);
        }
        int exponent_int = static_cast<int>(exponent_value);
        auto result = Result<ApproxComplex>::success(approx_complex(1.0, 0.0));
        for (int i = 0; i < std::abs(exponent_int); ++i) {
            result = multiply_complex(result.value(), base.value());
            if (!result) return result;
        }
        if (exponent_int < 0) {
            result = divide_complex(approx_complex(1.0, 0.0), result.value());
        }
        return result;
    }

    return real_to_complex(
        evaluate_numeric(*lamina::detail::make_expression_ptr(node),
                         bindings, context));
}

} // namespace

Result<EqvProfile> eqv_profile_from_name(const std::string& name) {
    if (name == "Core") {
        return Result<EqvProfile>::success(EqvProfile::Core);
    }
    if (name == "Trig-Basic") {
        return Result<EqvProfile>::success(EqvProfile::TrigBasic);
    }
    if (name == "ExpLog-Basic") {
        return Result<EqvProfile>::success(EqvProfile::ExpLogBasic);
    }
    return eqv_profile_failure("unsupported equivalence profile: " + name);
}

Result<void> set_eqv_profile(EqvOptions& options,
                             const std::string& name) {
    auto profile = eqv_profile_from_name(name);
    if (!profile) {
        return Result<void>::failure(profile.error());
    }
    options.profile = profile.value();
    return Result<void>::success();
}

Result<void> set_eqv_budget(EqvOptions& options,
                            std::size_t steps,
                            std::size_t depth,
                            std::size_t growth) {
    EqvOptions candidate = options;
    candidate.budget.max_rewrite_steps = steps;
    candidate.budget.max_rewrite_depth = depth;
    candidate.budget.max_node_growth_factor = growth;
    auto valid = validate_eqv_options(candidate);
    if (!valid) {
        return valid;
    }
    options = candidate;
    return Result<void>::success();
}

ExprSet::ExprSet(std::vector<ExprPtr> elements, ExprPtr expression)
    : elements_(std::move(elements)), expression_(std::move(expression)) {}

Result<ExprSet> ExprSet::make(std::vector<ExprPtr> elements) {
    for (const auto& element : elements) {
        if (!element) {
            return expr_set_failure(CasErrc::InvalidArgument,
                                    "set<Expr> elements cannot be null",
                                    kExprSetOperation);
        }
    }
    ComputationContext context;
    auto expression = make_finite_set(elements, context);
    if (!expression) return Result<ExprSet>::failure(expression.error());
    auto node = std::dynamic_pointer_cast<const FiniteSetNode>(
        lamina::detail::node(expression.value()));
    std::vector<ExprPtr> unique;
    unique.reserve(node->elements().size());
    for (const auto& element : node->elements()) {
        unique.push_back(lamina::detail::make_expression_ptr(element));
    }
    return Result<ExprSet>::success(
        ExprSet(std::move(unique), expression.value()));
}

bool ExprSet::contains(const SymbolicExpr& expression) const {
    for (const auto& element : elements_) {
        if (element && structurally_equal(*element, expression)) {
            return true;
        }
    }
    return false;
}

bool ExprSet::subset_of(const ExprSet& other) const {
    for (const auto& element : elements_) {
        if (!element || !other.contains(*element)) {
            return false;
        }
    }
    return true;
}

ExprSet ExprSet::set_union(const ExprSet& other) const {
    std::vector<ExprPtr> result = elements_;
    for (const auto& element : other.elements_) {
        if (element && !contains(*element)) {
            result.push_back(element);
        }
    }
    return ExprSet::make(std::move(result)).value();
}

ExprSet ExprSet::intersection(const ExprSet& other) const {
    std::vector<ExprPtr> result;
    for (const auto& element : elements_) {
        if (element && other.contains(*element)) {
            result.push_back(element);
        }
    }
    return ExprSet::make(std::move(result)).value();
}

ExprSet ExprSet::difference(const ExprSet& other) const {
    std::vector<ExprPtr> result;
    for (const auto& element : elements_) {
        if (element && !other.contains(*element)) {
            result.push_back(element);
        }
    }
    return ExprSet::make(std::move(result)).value();
}

ExprSet ExprSet::symmetric_difference(const ExprSet& other) const {
    auto left_only = difference(other);
    auto right_only = other.difference(*this);
    return left_only.set_union(right_only);
}

const char* NumberDomainSet::name() const noexcept {
    switch (domain_) {
    case NumberDomain::Integers:
        return "Z";
    case NumberDomain::Rationals:
        return "Q";
    case NumberDomain::Reals:
        return "R";
    case NumberDomain::Complexes:
        return "C";
    case NumberDomain::Expressions:
        return "Expr";
    }
    return "?";
}

bool NumberDomainSet::subset_of(const NumberDomainSet& other) const noexcept {
    return domain_rank(domain_) <= domain_rank(other.domain_);
}

Result<bool> NumberDomainSet::contains(const ExprPtr& element) const {
    if (!element || !lamina::detail::node(element)) {
        return bool_failure(CasErrc::InvalidArgument,
                            "domain membership element must not be null",
                            kNumberDomainOperation);
    }
    return domain_contains_node(domain_, lamina::detail::node(element));
}

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

Result<ApproxReal> evalf(const SymbolicExpr& expression,
                         const NumericBindings& bindings,
                         ComputationContext& context) {
    auto evaluated = evaluate_numeric(expression, bindings, context);
    if (!evaluated) return evaluated;
    if (!evaluated.value().is_finite() ||
        !std::isfinite(evaluated.value().value)) {
        return Result<ApproxReal>::failure(
            CasErrc::NumericFailure,
            "LSR evalf produced a non-finite result",
            kEvalfOperation);
    }
    return evaluated;
}

Result<ApproxReal> evalf(const SymbolicExpr& expression,
                         const NumericBindings& bindings) {
    ComputationContext context;
    return evalf(expression, bindings, context);
}

Result<ApproxComplex> eval_complex(const SymbolicExpr& expression,
                                   const NumericBindings& bindings,
                                   ComputationContext& context) {
    try {
        if (!lamina::detail::node(expression)) {
            return complex_failure(CasErrc::InvalidArgument,
                                   "cannot evaluate an empty expression as complex",
                                   kEvalComplexOperation);
        }
        return evaluate_complex_node(lamina::detail::node(expression),
                                     bindings, context);
    } catch (const std::bad_alloc&) {
        return complex_failure(CasErrc::ResourceLimit,
                               "complex evaluation allocation failed",
                               kEvalComplexOperation);
    } catch (const std::exception& error) {
        return complex_failure(CasErrc::UnsupportedExpression, error.what(),
                               kEvalComplexOperation);
    }
}

Result<ApproxComplex> eval_complex(const SymbolicExpr& expression,
                                   const NumericBindings& bindings) {
    ComputationContext context;
    return eval_complex(expression, bindings, context);
}

SolveResult solve_set(const ExprPtr& equation,
                      const std::string& variable,
                      ComputationContext& context,
                      const SolveOptions& options) {
    if (!equation) {
        return SolveResult::failure(CasErrc::InvalidArgument,
                                    "equation cannot be null", "lsr.solve_set");
    }
    if (variable.empty()) {
        return SolveResult::failure(CasErrc::InvalidArgument,
                                    "solve variable cannot be empty",
                                    "lsr.solve_set");
    }
    return solve_equation(equation, variable, context, options);
}

SolveResult solve_set(const ExprPtr& equation,
                      const std::string& variable,
                      const SolveOptions& options) {
    ComputationContext context;
    return solve_set(equation, variable, context, options);
}

ExprSetResult expr_set(std::vector<ExprPtr> elements) {
    try {
        return ExprSet::make(std::move(elements));
    } catch (const std::bad_alloc&) {
        return expr_set_failure(CasErrc::ResourceLimit,
                                "set<Expr> allocation failed",
                                kExprSetOperation);
    } catch (const std::exception& error) {
        return expr_set_failure(CasErrc::InvalidArgument, error.what(),
                                kExprSetOperation);
    }
}

NumberDomainSet integers() {
    return NumberDomainSet(NumberDomain::Integers);
}

NumberDomainSet rationals() {
    return NumberDomainSet(NumberDomain::Rationals);
}

NumberDomainSet reals() {
    return NumberDomainSet(NumberDomain::Reals);
}

NumberDomainSet complexes() {
    return NumberDomainSet(NumberDomain::Complexes);
}

NumberDomainSet expressions() {
    return NumberDomainSet(NumberDomain::Expressions);
}

Result<bool> domain_contains(const NumberDomainSet& domain,
                             const ExprPtr& element) {
    return domain.contains(element);
}

Result<bool> domain_subset(const NumberDomainSet& lhs,
                           const NumberDomainSet& rhs) {
    return Result<bool>::success(lhs.subset_of(rhs));
}

Result<bool> expr_set_contains(const ExprSet& set,
                               const ExprPtr& element) {
    if (!element) {
        return Result<bool>::failure(CasErrc::InvalidArgument,
                                     "set<Expr> membership element cannot be null",
                                     kExprSetOperation);
    }
    return Result<bool>::success(set.contains(*element));
}

Result<bool> expr_set_not_contains(const ExprSet& set,
                                   const ExprPtr& element) {
    auto result = expr_set_contains(set, element);
    if (!result) return result;
    return Result<bool>::success(!result.value());
}

Result<bool> expr_set_subset(const ExprSet& lhs,
                             const ExprSet& rhs) {
    return Result<bool>::success(lhs.subset_of(rhs));
}

Result<bool> expr_set_subset_domain(const ExprSet& set,
                                    const NumberDomainSet& domain) {
    for (const auto& element : set.elements()) {
        auto contains = domain.contains(element);
        if (!contains) {
            return contains;
        }
        if (!contains.value()) {
            return Result<bool>::success(false);
        }
    }
    return Result<bool>::success(true);
}

ExprSetResult expr_set_union(const ExprSet& lhs,
                             const ExprSet& rhs) {
    try {
        return ExprSetResult::success(lhs.set_union(rhs));
    } catch (const std::bad_alloc&) {
        return expr_set_failure(CasErrc::ResourceLimit,
                                "set<Expr> union allocation failed",
                                kExprSetOperation);
    }
}

ExprSetResult expr_set_intersection(const ExprSet& lhs,
                                    const ExprSet& rhs) {
    try {
        return ExprSetResult::success(lhs.intersection(rhs));
    } catch (const std::bad_alloc&) {
        return expr_set_failure(CasErrc::ResourceLimit,
                                "set<Expr> intersection allocation failed",
                                kExprSetOperation);
    }
}

ExprSetResult expr_set_difference(const ExprSet& lhs,
                                  const ExprSet& rhs) {
    try {
        return ExprSetResult::success(lhs.difference(rhs));
    } catch (const std::bad_alloc&) {
        return expr_set_failure(CasErrc::ResourceLimit,
                                "set<Expr> difference allocation failed",
                                kExprSetOperation);
    }
}

ExprSetResult expr_set_symmetric_difference(const ExprSet& lhs,
                                            const ExprSet& rhs) {
    try {
        return ExprSetResult::success(lhs.symmetric_difference(rhs));
    } catch (const std::bad_alloc&) {
        return expr_set_failure(CasErrc::ResourceLimit,
                                "set<Expr> symmetric difference allocation failed",
                                kExprSetOperation);
    }
}

ExprSetResult solve_expr_set(const ExprPtr& equation,
                             const std::string& variable,
                             ComputationContext& context,
                             const SolveOptions& options) {
    auto closed_form = try_lsr_closed_form_rational_poly_roots(
        equation, variable, context);
    if (!closed_form) {
        return ExprSetResult::failure(closed_form.error());
    }
    if (closed_form.value()) {
        return ExprSetResult::success(std::move(*closed_form.value()));
    }

    auto solved = solve_set(equation, variable, context, options);
    if (!solved) {
        return ExprSetResult::failure(solved.error());
    }

    const auto& solution_set = solved.value();
    if (solution_set.kind() == SolutionSet::Kind::Empty) {
        return expr_set({});
    }
    if (solution_set.kind() != SolutionSet::Kind::Finite) {
        std::string reason = solution_set.reason();
        if (reason.empty()) {
            reason = "solution set is not a finite enumerable set<Expr>";
        }
        return expr_set_failure(CasErrc::Inconclusive, std::move(reason),
                                kSolveExprSetOperation);
    }

    std::vector<ExprPtr> elements;
    elements.reserve(solution_set.finite_solutions().size());
    for (const auto& solution : solution_set.finite_solutions()) {
        if (!solution.conditions.empty()) {
            return expr_set_failure(CasErrc::Inconclusive,
                                    "conditional finite solutions cannot be lowered to set<Expr>",
                                    kSolveExprSetOperation);
        }
        elements.push_back(solution.value);
    }
    return expr_set(std::move(elements));
}

ExprSetResult solve_expr_set(const ExprPtr& equation,
                             const std::string& variable,
                             const SolveOptions& options) {
    ComputationContext context;
    return solve_expr_set(equation, variable, context, options);
}

ExprSetResult roots(const ExprPtr& expression,
                    const std::string& variable,
                    ComputationContext& context,
                    const SolveOptions& options) {
    return solve_expr_set(expression, variable, context, options);
}

ExprSetResult roots(const ExprPtr& expression,
                    const std::string& variable,
                    const SolveOptions& options) {
    ComputationContext context;
    return roots(expression, variable, context, options);
}

ExprSetResult solve(const ExprPtr& equation,
                    const std::string& variable,
                    ComputationContext& context,
                    const SolveOptions& options) {
    return solve_expr_set(equation, variable, context, options);
}

ExprSetResult solve(const ExprPtr& equation,
                    const std::string& variable,
                    const SolveOptions& options) {
    ComputationContext context;
    return solve(equation, variable, context, options);
}

const char* error_name(CasErrc code) noexcept {
    switch (code) {
    case CasErrc::InvalidArgument:
        return "InvalidArgument";
    case CasErrc::ParseError:
        return "ParseError";
    case CasErrc::UnboundSymbol:
        return "UnboundSymbol";
    case CasErrc::DomainError:
        return "DomainError";
    case CasErrc::UnsupportedExpression:
        return "UnsupportedExpression";
    case CasErrc::Inconclusive:
        return "Inconclusive";
    case CasErrc::ResourceLimit:
        return "ResourceLimit";
    case CasErrc::Cancelled:
        return "Cancelled";
    case CasErrc::NumericFailure:
        return "NumericFailure";
    case CasErrc::InternalInvariant:
        return "InternalInvariant";
    case CasErrc::DimensionMismatch:
        return "DimensionMismatch";
    case CasErrc::UnitInvalid:
        return "UnitInvalid";
    case CasErrc::UnitStripTypeMismatch:
        return "UnitStripTypeMismatch";
    case CasErrc::SetElementTypeMismatch:
        return "SetElementTypeMismatch";
    case CasErrc::SetOperandTypeMismatch:
        return "SetOperandTypeMismatch";
    case CasErrc::SetElementNotHashable:
        return "SetElementNotHashable";
    }
    return "InternalInvariant";
}

const char* error_name(const CasError& error) noexcept {
    if (error.operation == kSymOperation &&
        error.code == CasErrc::InvalidArgument &&
        error.message.find("imaginary unit") != std::string::npos) {
        return "ImaginaryUnitReserved";
    }
    if (error.operation == kEvalComplexOperation &&
        error.code == CasErrc::UnboundSymbol) {
        return "ComplexEvalUnboundSymbol";
    }
    if (error.operation == kComplexOperation &&
        error.code == CasErrc::InvalidArgument) {
        return "ComplexTypeMismatch";
    }
    if (error.operation == kExprSetOperation &&
        error.code == CasErrc::InvalidArgument) {
        return "SetElementTypeMismatch";
    }
    if (error.operation == kSolveExprSetOperation &&
        error.code == CasErrc::Inconclusive) {
        return "SetResultInconclusive";
    }
    if (error.operation == kEquivalentOperation &&
        error.code == CasErrc::ResourceLimit) {
        return "EqvBudgetExceeded";
    }
    if (error.operation == kEquivalentProfileOperation &&
        error.code == CasErrc::UnsupportedExpression) {
        return "EqvRuleDisabled";
    }
    return error_name(error.code);
}

bool structurally_equal(const SymbolicExpr& lhs, const SymbolicExpr& rhs) {
    const auto& left = lamina::detail::node(lhs);
    const auto& right = lamina::detail::node(rhs);
    if (!left || !right) return left == right;
    return left->equals(*right);
}

Result<bool> equivalent_core(const SymbolicExpr& lhs,
                             const SymbolicExpr& rhs,
                             ComputationContext& context,
                             const EqvOptions& options) {
    auto options_valid = validate_eqv_options(options);
    if (!options_valid) return Result<bool>::failure(options_valid.error());

    auto step = context.consume_steps(1, kEquivalentOperation);
    if (!step) return Result<bool>::failure(step.error());
    try {
        auto lhs_dimension = dimension_of(lhs);
        if (!lhs_dimension) return Result<bool>::failure(lhs_dimension.error());
        auto rhs_dimension = dimension_of(rhs);
        if (!rhs_dimension) return Result<bool>::failure(rhs_dimension.error());
        if (lhs_dimension.value() != rhs_dimension.value()) {
            return Result<bool>::success(false);
        }
        auto lhs_ptr = std::make_shared<SymbolicExpr>(lhs);
        auto rhs_ptr = std::make_shared<SymbolicExpr>(rhs);
        if (std::dynamic_pointer_cast<const QuantityNode>(lamina::detail::node(lhs))) {
            auto stripped = strip_unit(lhs_ptr, UnitStripMode::BaseValue, context);
            if (!stripped) return Result<bool>::failure(stripped.error());
            lhs_ptr = stripped.value();
        }
        if (std::dynamic_pointer_cast<const QuantityNode>(lamina::detail::node(rhs))) {
            auto stripped = strip_unit(rhs_ptr, UnitStripMode::BaseValue, context);
            if (!stripped) return Result<bool>::failure(stripped.error());
            rhs_ptr = stripped.value();
        }
        auto canonical_lhs = canonicalize_lsr_complex_product(*lhs_ptr);
        auto canonical_rhs = canonicalize_lsr_complex_product(*rhs_ptr);
        if (structurally_equal(*canonical_lhs, *canonical_rhs)) {
            return Result<bool>::success(true);
        }
        auto difference = SymbolicExpr::add(
            canonical_lhs,
            SymbolicExpr::multiply(SymbolicExpr::number(-1),
                                   canonical_rhs));
        if (!difference) {
            return Result<bool>::failure(CasErrc::InternalInvariant,
                                         "equivalence difference construction failed",
                                         kEquivalentOperation);
        }
        if (difference->simplify()->is_zero()) {
            return Result<bool>::success(true);
        }

        auto polynomial_proof = prove_rational_polynomial_equivalence(
            difference, context, options);
        if (!polynomial_proof) {
            return Result<bool>::failure(polynomial_proof.error());
        }
        if (polynomial_proof.value()) {
            return Result<bool>::success(*polynomial_proof.value());
        }
        if (options.profile == EqvProfile::TrigBasic) {
            if (options.budget.max_rewrite_steps < 8) {
                return Result<bool>::failure(
                    CasErrc::ResourceLimit,
                    "equivalence rewrite budget exhausted before Trig-Basic normalization",
                    kEquivalentOperation);
            }
            auto trig_step = context.consume_steps(8, kEquivalentOperation);
            if (!trig_step) return Result<bool>::failure(trig_step.error());
            auto trig_lhs = rewrite_trig_basic_identity(lamina::detail::node(lhs));
            auto trig_rhs = rewrite_trig_basic_identity(lamina::detail::node(rhs));
            EqvOptions core_options = options;
            core_options.profile = EqvProfile::Core;
            return equivalent_core(*trig_lhs, *trig_rhs, context,
                                   core_options);
        }
        if (options.profile == EqvProfile::ExpLogBasic) {
            if (options.budget.max_rewrite_steps < 8) {
                return Result<bool>::failure(
                    CasErrc::ResourceLimit,
                    "equivalence rewrite budget exhausted before ExpLog-Basic normalization",
                    kEquivalentOperation);
            }
            auto exp_log_step = context.consume_steps(8, kEquivalentOperation);
            if (!exp_log_step) return Result<bool>::failure(exp_log_step.error());
            auto exp_log_lhs = rewrite_exp_log_basic_identity(
                lamina::detail::node(lhs), context.assumptions().get());
            auto exp_log_rhs = rewrite_exp_log_basic_identity(
                lamina::detail::node(rhs), context.assumptions().get());
            EqvOptions core_options = options;
            core_options.profile = EqvProfile::Core;
            return equivalent_core(*exp_log_lhs, *exp_log_rhs, context,
                                   core_options);
        }
        return Result<bool>::success(false);
    } catch (const std::bad_alloc&) {
        return Result<bool>::failure(CasErrc::ResourceLimit,
                                     "equivalence check allocation failed",
                                     kEquivalentOperation);
    } catch (const std::exception& error) {
        return Result<bool>::failure(CasErrc::Inconclusive, error.what(),
                                     kEquivalentOperation);
    }
}

Result<bool> equivalent_core(const SymbolicExpr& lhs,
                             const SymbolicExpr& rhs,
                             ComputationContext& context) {
    return equivalent_core(lhs, rhs, context, EqvOptions{});
}

Result<bool> equivalent(const SymbolicExpr& lhs,
                        const SymbolicExpr& rhs,
                        ComputationContext& context,
                        const EqvOptions& options) {
    auto checked = equivalent_core(lhs, rhs, context, options);
    if (checked) return checked;

    const auto code = checked.error().code;
    if (code == CasErrc::ResourceLimit ||
        code == CasErrc::Inconclusive ||
        code == CasErrc::UnsupportedExpression) {
        (void)context.add_diagnostic(
            Diagnostic{DiagnosticSeverity::Warning,
                       kEquivalentOperation,
                       checked.error().message});
        return Result<bool>::success(false);
    }
    return Result<bool>::failure(checked.error());
}

Result<bool> equivalent(const SymbolicExpr& lhs,
                        const SymbolicExpr& rhs,
                        ComputationContext& context) {
    return equivalent(lhs, rhs, context, EqvOptions{});
}

} // namespace lamina::lsr
