#include "lmcas_export.hpp"
#include "internal/expr_serialization_common.hpp"

#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace LMCAS {
namespace {
using namespace serialization_detail;

class Reader {
public:
    Reader(const std::string& source, ComputationContext& context)
        : source_(source), context_(context), position_(kHeader.size()) {}

    NodePtr finish() {
        auto result = node();
        if (position_ != source_.size()) invalid("trailing serialized expression data");
        return result;
    }

private:
    [[noreturn]] void invalid(const char* message) const {
        throw CasError{CasErrc::ParseError, message, kDecode};
    }
    std::size_t read_length() {
        std::size_t length = 0;
        while (position_ < source_.size() && source_[position_] >= '0' &&
               source_[position_] <= '9') {
            const auto digit = static_cast<std::size_t>(source_[position_++] - '0');
            if (length > (std::numeric_limits<std::size_t>::max() - digit) / 10)
                invalid("serialized expression length overflow");
            length = length * 10 + digit;
        }
        return length;
    }
    std::string_view field() {
        check(context_.consume_steps(1, kDecode));
        const auto begin = position_;
        const auto length = read_length();
        if (begin == position_) {
            invalid("invalid serialized expression length");
        }
        if (position_ == source_.size()) {
            invalid("invalid serialized expression length");
        }
        if (source_[position_++] != ':') {
            invalid("invalid serialized expression length");
        }
        if (position_ - begin > 2 && source_[begin] == '0') {
            invalid("invalid serialized expression length");
        }
        if (length >= source_.size() - position_) {
            invalid("truncated serialized expression field");
        }
        if (source_[position_ + length] != ',') {
            invalid("truncated serialized expression field");
        }
        auto result = std::string_view(source_).substr(position_, length);
        position_ += length + 1;
        return result;
    }
    std::size_t unsigned_field() {
        const auto text = field();
        std::size_t result = 0;
        auto parsed = std::from_chars(text.data(), text.data() + text.size(), result);
        if (text.empty() || (text.size() > 1 && text[0] == '0') ||
            parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
            invalid("invalid serialized expression integer");
        return result;
    }
    std::size_t count() {
        auto n = unsigned_field();
        if (n > source_.size() - position_)
            invalid("invalid serialized expression count");
        if (n > context_.limits().max_ast_nodes)
            throw CasError{CasErrc::ResourceLimit, "AST node budget exhausted", kDecode};
        return n;
    }
    bool flag() {
        auto text = field();
        if (text == "0") return false;
        if (text == "1") return true;
        invalid("invalid serialized expression flag");
    }
    template <typename Enum, std::size_t N>
    Enum choice(const std::array<std::string_view, N>& names) {
        auto text = field();
        for (std::size_t i = 0; i < names.size(); ++i)
            if (text == names[i]) return static_cast<Enum>(i);
        invalid("invalid serialized expression operator");
    }
    void check_integer_digits(std::string_view text, std::size_t first) const {
        for (std::size_t i = first; i < text.size(); ++i)
            if (text[i] < '0' || text[i] > '9') invalid("invalid serialized expression integer");
    }
    BigInt integer() {
        auto text = field();
        if (text.empty()) {
            invalid("invalid serialized expression integer");
        }
        std::size_t first = text.front() == '-' ? 1 : 0;
        if (first == text.size()) {
            invalid("invalid serialized expression integer");
        }
        if (text.size() - first > 1 && text[first] == '0') {
            invalid("invalid serialized expression integer");
        }
        check_integer_digits(text, first);
        if (first != 0 && text.substr(first) == "0") {
            invalid("invalid serialized expression integer");
        }
        // Decimal length gives an inexpensive lower bound before allocating BigInt.
        if (text.size() - first > context_.limits().max_integer_bits / 3 + 1) {
            throw CasError{CasErrc::ResourceLimit, "integer bit budget exceeded", kDecode};
        }
        BigInt result{std::string(text)};
        check(context_.require_integer_bits(result.bit_length(), kDecode));
        return result;
    }
    Rational rational() {
        auto numerator = integer();
        auto denominator = integer();
        if (denominator <= BigInt(0)) invalid("invalid serialized expression denominator");
        return Rational(numerator, denominator);
    }
    std::vector<NodePtr> nodes() {
        auto n = count();
        std::vector<NodePtr> result;
        result.reserve(n);
        for (std::size_t i = 0; i < n; ++i) result.push_back(node());
        return result;
    }
    int hex_digit(char c) const {
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'a' && c <= 'f') return c - 'a' + 10;
        invalid("invalid binary64 bits");
    }
    NodePtr binary64() {
        auto text = field();
        if (text.size() != 16) invalid("invalid binary64 bits");
        std::uint64_t bits = 0;
        for (char c : text) {
            const auto digit = hex_digit(c);
            bits = bits * 16 + static_cast<unsigned>(digit);
        }
        double value;
        std::memcpy(&value, &bits, sizeof(value));
        if (!std::isfinite(value)) invalid("nonfinite binary64 value");
        return detail::make_node<NumberNode>(value);
    }
    NodePtr matrix() {
        auto rows = unsigned_field(), cols = unsigned_field();
        if (!rows || !cols || rows > std::numeric_limits<std::size_t>::max() / cols ||
            rows * cols > source_.size() - position_)
            invalid("invalid matrix dimensions");
        if (rows * cols > context_.limits().max_ast_nodes)
            throw CasError{CasErrc::ResourceLimit, "AST node budget exhausted", kDecode};
        MatrixNode::DenseStorage elements;
        elements.reserve(rows * cols);
        for (std::size_t i = 0; i < rows * cols; ++i) elements.push_back(node());
        return detail::make_node<MatrixNode>(rows, cols, std::move(elements));
    }
    NodePtr piecewise() {
        auto count_branches = count();
        std::vector<PiecewiseNode::Branch> branches;
        branches.reserve(count_branches);
        for (std::size_t i = 0; i < count_branches; ++i) {
            auto expression = node();
            branches.push_back({expression, node()});
        }
        auto has_default = flag();
        return detail::make_node<PiecewiseNode>(std::move(branches), has_default ? node() : nullptr);
    }
    NodePtr quantity() {
        auto dimensions = count();
        DimensionSignature::Exponents exponents;
        for (std::size_t i = 0; i < dimensions; ++i) {
            auto name = std::string(field());
            auto exponent = rational();
            if (name.empty() || exponent.is_zero() ||
                !exponents.emplace(std::move(name), std::move(exponent)).second)
                invalid("invalid dimension exponent");
        }
        auto dimension = DimensionSignature(std::move(exponents));
        return detail::make_node<QuantityNode>(node(), std::move(dimension), Rational(1), "");
    }
    NodePtr root_of() {
        auto name = std::string(field());
        auto index = unsigned_field(), size = count();
        std::vector<Rational> coefficients;
        coefficients.reserve(size);
        for (std::size_t i = 0; i < size; ++i) coefficients.push_back(rational());
        if (coefficients.size() < 2 || coefficients.back().is_zero() ||
            index >= coefficients.size() - 1) invalid("invalid root identity");
        detail::ExactRootId id{Polynomial<Rational>(coefficients), index};
        return detail::make_node<RootOfNode>(std::move(id), std::move(name));
    }
    NodePtr node() {
        Recursion depth(context_, kDecode);
        check(context_.reserve_nodes(1, kDecode));
        auto kind = field();
        if (kind == "integer") return detail::make_node<NumberNode>(integer());
        if (kind == "rational") return detail::make_node<NumberNode>(rational());
        if (kind == "binary64") return binary64();
        if (kind == "symbol") {
            auto name = field();
            return detail::make_node<VariableNode>(std::string(name));
        }
        if (kind == "constant") {
            auto name = field();
            if (name != "pi" && name != "π" && name != "e" && name != "phi")
                throw CasError{CasErrc::UnsupportedExpression, "unknown mathematical constant", kDecode};
            return detail::make_node<VariableNode>(std::string(name), true);
        }
        if (kind == "add") return detail::make_node<AddNode>(nodes());
        if (kind == "multiply") return detail::make_node<MultiplyNode>(nodes());
        return compound_node(kind);
    }
    NodePtr compound_node(std::string_view kind) {
        if (kind == "power") { auto base = node(); return detail::make_node<PowerNode>(base, node()); }
        if (kind == "complex") { auto real = node(); return detail::make_node<ComplexNode>(real, node()); }
        if (kind == "function") {
            auto type = choice<FunctionNode::FuncType>(kFunctions);
            return detail::make_node<FunctionNode>(type, nodes());
        }
        if (kind == "call") {
            auto name = field();
            if (name.empty()) invalid("empty function name");
            return detail::make_node<UninterpretedFunctionNode>(std::string(name), nodes());
        }
        if (kind == "matrix") return matrix();
        if (kind == "relation") {
            auto op = choice<RelationOp>(kRelations);
            auto lhs = node(); return detail::make_node<RelationalNode>(lhs, node(), op);
        }
        if (kind == "logic") {
            auto op = choice<LogicalNode::Op>(kLogic);
            auto lhs = node();
            return detail::make_node<LogicalNode>(lhs, op == LogicalNode::Op::Not ? nullptr : node(), op);
        }
        if (kind == "piecewise") return piecewise();
        return extended_node(kind);
    }
    NodePtr extended_node(std::string_view kind) {
        if (kind == "sum" || kind == "product") {
            auto name = std::string(field());
            auto lower = node(), upper = node(), body = node();
            if (kind == "sum") return detail::make_node<SummationNode>(body, std::move(name), lower, upper);
            return detail::make_node<ProductNode>(body, std::move(name), lower, upper);
        }
        if (kind == "transform") {
            auto type = choice<TransformNode::TransformType>(kTransforms);
            auto name = std::string(field());
            auto target = node(), body = node();
            return detail::make_node<TransformNode>(type, body, std::move(name), target);
        }
        if (kind == "quantifier") {
            auto type = choice<QuantifierNode::Type>(kQuantifiers);
            auto name = std::string(field());
            auto domain = node(), predicate = node();
            return detail::make_node<QuantifierNode>(type, std::move(name), domain, predicate);
        }
        if (kind == "set_builder") {
            auto name = std::string(field());
            auto domain = node(), predicate = node();
            return detail::make_node<SetBuilderNode>(std::move(name), domain, predicate);
        }
        if (kind == "set") return detail::make_node<FiniteSetNode>(nodes());
        if (kind == "interval") {
            auto lower_closed = flag(), upper_closed = flag();
            auto lower = node();
            return detail::make_node<IntervalNode>(lower, node(), lower_closed, upper_closed);
        }
        if (kind == "membership") {
            auto element = node(); return detail::make_node<MembershipNode>(element, node());
        }
        if (kind == "quantity") return quantity();
        if (kind == "integral") {
            auto name = std::string(field());
            auto bounded = flag();
            NodePtr lower, upper;
            if (bounded) { lower = node(); upper = node(); }
            return detail::make_node<IntegralNode>(node(), std::move(name), lower, upper);
        }
        if (kind == "limit") {
            auto name = std::string(field());
            auto direction = choice<LimitDirection>(kDirections);
            auto point = node();
            return detail::make_node<LimitNode>(node(), std::move(name), point, direction);
        }
        if (kind == "root_of") return root_of();
        throw CasError{CasErrc::UnsupportedExpression, "unknown expression kind", kDecode};
    }

    const std::string& source_;
    ComputationContext& context_;
    std::size_t position_;
};
} // namespace

LMCAS_API ExpressionResult parse_serialized_expr(const std::string& source,
                                                 ComputationContext& context) {
    auto input = context.require_input_bytes(source.size(), kDecode);
    if (!input) return ExpressionResult::failure(input.error());
    if (source.compare(0, kHeader.size(), kHeader.data(), kHeader.size()) != 0)
        return ExpressionResult::failure(CasErrc::ParseError,
                                         "unknown serialized expression version", kDecode);
    try {
        return detail::make_expression_ptr(Reader(source, context).finish());
    } catch (const CasError& error) {
        return ExpressionResult::failure(error);
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                         "expression parse allocation failed", kDecode);
    } catch (const std::exception& error) {
        return ExpressionResult::failure(CasErrc::ParseError, error.what(), kDecode);
    }
}

LMCAS_API ExpressionResult parse_serialized_expr(const std::string& source) {
    ComputationContext context;
    return parse_serialized_expr(source, context);
}
} // namespace LMCAS
