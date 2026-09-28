#include "expr.hpp"
#include "computation_context.hpp"
#include "internal/symbolic_ast.hpp"

#include <array>
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

constexpr std::string_view kHeader = "LMCAS_EXPR/1\n";
constexpr const char* kEncode = "LMCAS.serialize_expr";
constexpr const char* kDecode = "LMCAS.parse_serialized_expr";
using NodePtr = detail::SymbolicNodePtr;

// These names, not the C++ enum ordinals, are the wire values.
constexpr std::array<std::string_view, 36> kFunctions = {
    "sin", "cos", "tan", "cot", "sec", "csc", "asin", "acos", "atan",
    "sinh", "cosh", "tanh", "ln", "log", "abs", "sqrt", "exp",
    "lambert_w", "atan2", "infinity", "erf", "ei", "si", "ci", "li",
    "max", "min", "sgn", "floor", "ceil", "round", "real", "imag",
    "conjugate", "complex_abs", "complex_arg"
};
constexpr std::array<std::string_view, 6> kRelations = {
    "eq", "neq", "lt", "gt", "leq", "geq"
};
constexpr std::array<std::string_view, 4> kLogic = {
    "and", "or", "not", "implies"
};
constexpr std::array<std::string_view, 5> kTransforms = {
    "laplace", "inverse_laplace", "fourier", "inverse_fourier", "z"
};
constexpr std::array<std::string_view, 2> kQuantifiers = {"forall", "exists"};
constexpr std::array<std::string_view, 3> kDirections = {
    "both", "below", "above"
};

template <typename Enum, std::size_t N>
std::string_view enum_name(Enum value, const std::array<std::string_view, N>& names,
                           const char* operation) {
    auto index = static_cast<std::size_t>(value);
    if (index >= names.size()) {
        throw CasError{CasErrc::UnsupportedExpression, "unsupported expression operator", operation};
    }
    return names[index];
}

void check(Result<void> result) {
    if (!result) throw result.error();
}

class Recursion {
public:
    Recursion(ComputationContext& context, const char* operation)
        : context_(context) {
        check(context_.enter_recursion(operation));
    }
    ~Recursion() { context_.leave_recursion(); }
    Recursion(const Recursion&) = delete;
    Recursion& operator=(const Recursion&) = delete;
private:
    ComputationContext& context_;
};

class Writer : public detail::SymbolicVisitor {
public:
    explicit Writer(ComputationContext& context) : context_(context), output_(kHeader) {}

    std::string finish(const NodePtr& root) {
        emit(root);
        return std::move(output_);
    }

    void visit(const NumberNode& n) override {
        const auto& value = n.value();
        if (const auto* integer = std::get_if<BigInt>(&value)) {
            put("integer");
            put_integer(*integer);
        } else if (const auto* fraction = std::get_if<Rational>(&value)) {
            put("rational");
            put_rational(*fraction);
        } else {
            put("binary64");
            static_assert(sizeof(lmmc_real_t) == sizeof(std::uint64_t), "binary64 required");
            std::uint64_t bits;
            auto number = std::get<lmmc_real_t>(value);
            std::memcpy(&bits, &number, sizeof(bits));
            char hex[16];
            constexpr char digits[] = "0123456789abcdef";
            for (int i = 15; i >= 0; --i) {
                hex[i] = digits[bits & 15];
                bits >>= 4;
            }
            put(std::string_view(hex, sizeof(hex)));
        }
    }
    void visit(const VariableNode& n) override {
        put(n.is_constant() ? "constant" : "symbol"); put(n.name());
    }
    void visit(const AddNode& n) override { put("add"); emit_list(n.operands()); }
    void visit(const MultiplyNode& n) override { put("multiply"); emit_list(n.operands()); }
    void visit(const PowerNode& n) override {
        put("power"); emit(n.base()); emit(n.exponent());
    }
    void visit(const FunctionNode& n) override {
        put("function"); put(enum_name(n.type(), kFunctions, kEncode));
        emit_list(n.arguments());
    }
    void visit(const UninterpretedFunctionNode& n) override {
        put("call"); put(n.name()); emit_list(n.arguments());
    }
    void visit(const MatrixNode& n) override {
        put("matrix"); put_size(n.rows()); put_size(n.cols());
        for (std::size_t r = 0; r < n.rows(); ++r)
            for (std::size_t c = 0; c < n.cols(); ++c)
                emit(n.get(r, c));
    }
    void visit(const RelationalNode& n) override {
        put("relation"); put(enum_name(n.op(), kRelations, kEncode));
        emit(n.left()); emit(n.right());
    }
    void visit(const LogicalNode& n) override {
        put("logic"); put(enum_name(n.op(), kLogic, kEncode));
        emit(n.left());
        if (n.op() != LogicalNode::Op::Not) emit(n.right());
    }
    void visit(const PiecewiseNode& n) override {
        put("piecewise"); put_size(n.branches().size());
        for (const auto& branch : n.branches()) {
            emit(branch.expression); emit(branch.condition);
        }
        put(n.default_expr() ? "1" : "0");
        if (n.default_expr()) emit(n.default_expr());
    }
    void visit(const SummationNode& n) override {
        put("sum"); put(n.index_var()); emit(n.lower_bound());
        emit(n.upper_bound()); emit(n.body());
    }
    void visit(const ProductNode& n) override {
        put("product"); put(n.index_var()); emit(n.lower_bound());
        emit(n.upper_bound()); emit(n.body());
    }
    void visit(const TransformNode& n) override {
        put("transform"); put(enum_name(n.transform_type(), kTransforms, kEncode));
        put(n.source_var()); emit(n.target()); emit(n.body());
    }
    void visit(const QuantifierNode& n) override {
        put("quantifier"); put(enum_name(n.quantifier_type(), kQuantifiers, kEncode));
        put(n.bound_var()); emit(n.domain()); emit(n.predicate());
    }
    void visit(const SetBuilderNode& n) override {
        put("set_builder"); put(n.element_var()); emit(n.domain()); emit(n.predicate());
    }
    void visit(const FiniteSetNode& n) override { put("set"); emit_list(n.elements()); }
    void visit(const IntervalNode& n) override {
        put("interval"); put(n.lower_closed() ? "1" : "0");
        put(n.upper_closed() ? "1" : "0"); emit(n.lower()); emit(n.upper());
    }
    void visit(const MembershipNode& n) override {
        put("membership"); emit(n.element()); emit(n.set());
    }
    void visit(const QuantityNode& n) override {
        put("quantity"); put_size(n.dimension().exponents().size());
        for (const auto& [name, exponent] : n.dimension().exponents()) {
            put(name); put_rational(exponent);
        }
        // A quantity's physical value is its displayed value multiplied by its exact scale.
        NodePtr base = n.value();
        const auto* number = dynamic_cast<const NumberNode*>(base.get());
        if (n.scale_to_base() != Rational(1) ||
            dynamic_cast<const MultiplyNode*>(base.get()) ||
            (number && std::holds_alternative<Rational>(number->value()))) {
            Rational coefficient = n.scale_to_base();
            check(context_.require_integer_bits(coefficient.get_numerator().bit_length(), kEncode));
            check(context_.require_integer_bits(coefficient.get_denominator().bit_length(), kEncode));
            std::vector<NodePtr> factors;
            if (const auto* product = dynamic_cast<const MultiplyNode*>(base.get()))
                factors.assign(product->operands().begin(), product->operands().end());
            else
                factors.push_back(base);
            std::vector<NodePtr> remaining;
            remaining.reserve(factors.size() + 1);
            bool exact_factor = false;
            for (const auto& factor : factors) {
                const auto* number = dynamic_cast<const NumberNode*>(factor.get());
                if (!number || std::holds_alternative<lmmc_real_t>(number->value())) {
                    remaining.push_back(factor);
                    continue;
                }
                const Rational value = std::holds_alternative<BigInt>(number->value())
                    ? Rational(std::get<BigInt>(number->value()))
                    : std::get<Rational>(number->value());
                check(context_.require_integer_bits(value.get_numerator().bit_length(), kEncode));
                check(context_.require_integer_bits(value.get_denominator().bit_length(), kEncode));
                coefficient = coefficient * value;
                check(context_.require_integer_bits(coefficient.get_numerator().bit_length(), kEncode));
                check(context_.require_integer_bits(coefficient.get_denominator().bit_length(), kEncode));
                exact_factor = true;
            }
            if (exact_factor || coefficient != Rational(1)) {
                if (coefficient != Rational(1) || remaining.empty()) {
                    remaining.push_back(coefficient.is_integer()
                        ? SymbolicFactory::create_number(coefficient.to_bigint())
                        : SymbolicFactory::create_number(coefficient));
                }
                base = remaining.size() == 1 ? remaining.front()
                    : SymbolicFactory::create_multiply(std::move(remaining));
            }
        }
        emit(base);
    }
    void visit(const ComplexNode& n) override {
        put("complex"); emit(n.real()); emit(n.imag());
    }
    void visit(const IntegralNode& n) override {
        put("integral"); put(n.variable()); put(n.is_definite() ? "1" : "0");
        if (n.is_definite()) { emit(n.lower()); emit(n.upper()); }
        emit(n.body());
    }
    void visit(const LimitNode& n) override {
        put("limit"); put(n.variable()); put(enum_name(n.direction(), kDirections, kEncode));
        emit(n.point()); emit(n.body());
    }
    void visit(const RootOfNode& n) override {
        put("root_of"); put(n.variable()); put_size(n.index());
        put_size(n.exact_id().polynomial.coeffs.size());
        for (const auto& coefficient : n.exact_id().polynomial.coeffs)
            put_rational(coefficient);
    }

private:
    void put(std::string_view field) {
        check(context_.consume_steps(1, kEncode));
        char digits[std::numeric_limits<std::size_t>::digits10 + 2];
        auto end = std::to_chars(digits, digits + sizeof(digits), field.size()).ptr;
        const auto width = static_cast<std::size_t>(end - digits);
        const auto extra = width + 2;
        if (field.size() > context_.limits().max_input_bytes ||
            output_.size() > context_.limits().max_input_bytes - field.size() ||
            extra > context_.limits().max_input_bytes - field.size() - output_.size()) {
            throw CasError{CasErrc::ResourceLimit, "serialized expression byte limit exceeded", kEncode};
        }
        output_.append(digits, width);
        output_ += ':';
        output_.append(field.data(), field.size());
        output_ += ',';
    }
    void put_size(std::size_t value) {
        char digits[std::numeric_limits<std::size_t>::digits10 + 2];
        auto end = std::to_chars(digits, digits + sizeof(digits), value).ptr;
        put(std::string_view(digits, static_cast<std::size_t>(end - digits)));
    }
    void put_integer(const BigInt& value) {
        check(context_.require_integer_bits(value.bit_length(), kEncode));
        put(value.to_string());
    }
    void put_rational(const Rational& value) {
        put_integer(value.get_numerator());
        put_integer(value.get_denominator());
    }
    template <typename Children>
    void emit_list(const Children& children) {
        put_size(children.size());
        for (const auto& child : children) emit(child);
    }
    void emit(const NodePtr& node) {
        Recursion depth(context_, kEncode);
        if (!node) throw CasError{CasErrc::InvalidArgument, "null expression node", kEncode};
        try {
            node->accept(*this);
        } catch (const std::runtime_error& error) {
            if (std::string_view(error.what()) == "AST traversal depth limit exceeded")
                throw CasError{CasErrc::ResourceLimit, error.what(), kEncode};
            throw;
        }
    }

    ComputationContext& context_;
    std::string output_;
};

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
    std::string_view field() {
        check(context_.consume_steps(1, kDecode));
        std::size_t length = 0;
        const auto begin = position_;
        while (position_ < source_.size() && source_[position_] >= '0' &&
               source_[position_] <= '9') {
            const auto digit = static_cast<std::size_t>(source_[position_++] - '0');
            if (length > (std::numeric_limits<std::size_t>::max() - digit) / 10)
                invalid("serialized expression length overflow");
            length = length * 10 + digit;
        }
        if (begin == position_ || position_ == source_.size() ||
            source_[position_++] != ':' ||
            (position_ - begin > 2 && source_[begin] == '0'))
            invalid("invalid serialized expression length");
        if (length >= source_.size() - position_ || source_[position_ + length] != ',')
            invalid("truncated serialized expression field");
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
    BigInt integer() {
        auto text = field();
        if (text.empty()) invalid("invalid serialized expression integer");
        std::size_t first = text.front() == '-' ? 1 : 0;
        if (first == text.size() || (text.size() - first > 1 && text[first] == '0'))
            invalid("invalid serialized expression integer");
        for (std::size_t i = first; i < text.size(); ++i)
            if (text[i] < '0' || text[i] > '9') invalid("invalid serialized expression integer");
        if (first && text.substr(first) == "0") invalid("invalid serialized expression integer");
        // Decimal length gives an inexpensive lower bound before allocating BigInt.
        if (text.size() - first > context_.limits().max_integer_bits / 3 + 1)
            throw CasError{CasErrc::ResourceLimit, "integer bit budget exceeded", kDecode};
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
    NodePtr node() {
        Recursion depth(context_, kDecode);
        check(context_.reserve_nodes(1, kDecode));
        auto kind = field();
        if (kind == "integer") return detail::make_node<NumberNode>(integer());
        if (kind == "rational") return detail::make_node<NumberNode>(rational());
        if (kind == "binary64") {
            auto text = field();
            if (text.size() != 16) invalid("invalid binary64 bits");
            std::uint64_t bits = 0;
            for (char c : text) {
                auto digit = c >= '0' && c <= '9' ? c - '0' :
                             c >= 'a' && c <= 'f' ? c - 'a' + 10 : -1;
                if (digit < 0) invalid("invalid binary64 bits");
                bits = bits * 16 + static_cast<unsigned>(digit);
            }
            double value;
            std::memcpy(&value, &bits, sizeof(value));
            if (!std::isfinite(value)) invalid("nonfinite binary64 value");
            return detail::make_node<NumberNode>(value);
        }
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
        if (kind == "matrix") {
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
        if (kind == "relation") {
            auto op = choice<RelationOp>(kRelations);
            auto lhs = node(); return detail::make_node<RelationalNode>(lhs, node(), op);
        }
        if (kind == "logic") {
            auto op = choice<LogicalNode::Op>(kLogic);
            auto lhs = node();
            return detail::make_node<LogicalNode>(lhs, op == LogicalNode::Op::Not ? nullptr : node(), op);
        }
        if (kind == "piecewise") {
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
        if (kind == "quantity") {
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
        if (kind == "root_of") {
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
        throw CasError{CasErrc::UnsupportedExpression, "unknown expression kind", kDecode};
    }

    const std::string& source_;
    ComputationContext& context_;
    std::size_t position_;
};

} // namespace

Result<std::string> serialize_expr(const ExprPtr& expression,
                                   ComputationContext& context) {
    auto step = context.consume_steps(1, kEncode);
    if (!step) return Result<std::string>::failure(step.error());
    if (!expression || !detail::node(expression))
        return Result<std::string>::failure(CasErrc::InvalidArgument,
                                            "expression cannot be null", kEncode);
    try {
        return Writer(context).finish(detail::node(expression));
    } catch (const CasError& error) {
        return Result<std::string>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<std::string>::failure(CasErrc::ResourceLimit,
                                            "expression serialization allocation failed", kEncode);
    } catch (const std::exception& error) {
        return Result<std::string>::failure(CasErrc::UnsupportedExpression,
                                            error.what(), kEncode);
    }
}

Result<std::string> serialize_expr(const ExprPtr& expression) {
    ComputationContext context;
    return serialize_expr(expression, context);
}

ExprResult parse_serialized_expr(const std::string& source,
                                 ComputationContext& context) {
    auto input = context.require_input_bytes(source.size(), kDecode);
    if (!input) return ExprResult::failure(input.error());
    if (source.compare(0, kHeader.size(), kHeader.data(), kHeader.size()) != 0)
        return ExprResult::failure(CasErrc::ParseError,
                                   "unknown serialized expression version", kDecode);
    try {
        return detail::make_expression_ptr(Reader(source, context).finish());
    } catch (const CasError& error) {
        return ExprResult::failure(error);
    } catch (const std::bad_alloc&) {
        return ExprResult::failure(CasErrc::ResourceLimit,
                                   "expression parse allocation failed", kDecode);
    } catch (const std::exception& error) {
        return ExprResult::failure(CasErrc::ParseError, error.what(), kDecode);
    }
}

ExprResult parse_serialized_expr(const std::string& source) {
    ComputationContext context;
    return parse_serialized_expr(source, context);
}

} // namespace LMCAS
