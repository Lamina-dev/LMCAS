#include "lmcas_export.hpp"
#include "internal/expr_serialization_common.hpp"

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
using namespace serialization_detail;


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
        emit(base_value(n));
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
    bool needs_scale(const QuantityNode& n, const NodePtr& base) {
        if (n.scale_to_base() != Rational(1)) return true;
        if (dynamic_cast<const MultiplyNode*>(base.get())) return true;
        const auto* number = dynamic_cast<const NumberNode*>(base.get());
        return number && std::holds_alternative<Rational>(number->value());
    }
    bool collect_exact_factors(const std::vector<NodePtr>& factors,
                               Rational& coefficient, std::vector<NodePtr>& remaining) {
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
        return exact_factor;
    }
    NodePtr base_value(const QuantityNode& n) {
        NodePtr base = n.value();
        if (!needs_scale(n, base)) return base;
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
        const bool exact_factor = collect_exact_factors(factors, coefficient, remaining);
        return scaled_value(base, coefficient, remaining, exact_factor);
    }
    NodePtr scaled_value(const NodePtr& base, const Rational& coefficient,
                         std::vector<NodePtr>& remaining, bool exact_factor) {
        if (exact_factor || coefficient != Rational(1)) {
            if (coefficient != Rational(1) || remaining.empty()) {
                remaining.push_back(coefficient.is_integer()
                    ? SymbolicFactory::create_number(coefficient.to_bigint())
                    : SymbolicFactory::create_number(coefficient));
            }
            return remaining.size() == 1 ? remaining.front()
                : SymbolicFactory::create_multiply(std::move(remaining));
        }
        return base;
    }
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


} // namespace

LMCAS_API Result<std::string> serialize_expr(const ExprPtr& expression,
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

LMCAS_API Result<std::string> serialize_expr(const ExprPtr& expression) {
    ComputationContext context;
    return serialize_expr(expression, context);
}


} // namespace LMCAS
