#include "polynomial_conversion.hpp"
#include "internal/symbolic_ast.hpp"
#include <cstdint>
#include <utility>

namespace LMCAS {
namespace {

constexpr const char* kOperation = "recognize_rational_polynomial";

class RecursionScope {
public:
    explicit RecursionScope(ComputationContext& context) : context_(context) {}
    ~RecursionScope() { context_.leave_recursion(); }

private:
    ComputationContext& context_;
};

Result<OptionalRationalPolynomial> recognize_node(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& variable, ComputationContext& context);

Result<OptionalRationalPolynomial> recognize_number(
    const NumberNode& number, const std::string& variable) {
        if (std::holds_alternative<BigInt>(number.value())) {
            return Result<OptionalRationalPolynomial>::success(Polynomial<Rational>(
                Rational(std::get<BigInt>(number.value())), variable));
        }
        if (std::holds_alternative<Rational>(number.value())) {
            return Result<OptionalRationalPolynomial>::success(Polynomial<Rational>(
                std::get<Rational>(number.value()), variable));
        }
        return Result<OptionalRationalPolynomial>::success(std::nullopt);
}

Result<OptionalRationalPolynomial> recognize_sum(
    const AddNode& add, const std::string& variable, ComputationContext& context) {
        Polynomial<Rational> result(variable);
        for (const auto& operand : add.operands()) {
            auto child = recognize_node(operand, variable, context);
            if (!child) return child;
            if (!child.value()) {
                return Result<OptionalRationalPolynomial>::success(std::nullopt);
            }
            result = result + *child.value();
            if (result.coeffs.size() > context.limits().max_expansion_terms) {
                return Result<OptionalRationalPolynomial>::failure(
                    CasErrc::ResourceLimit, "polynomial term budget exhausted", kOperation);
            }
        }
        return Result<OptionalRationalPolynomial>::success(std::move(result));
}

Result<OptionalRationalPolynomial> recognize_product(
    const MultiplyNode& multiply, const std::string& variable, ComputationContext& context) {
        Polynomial<Rational> result({Rational(1)}, variable);
        for (const auto& operand : multiply.operands()) {
            auto child = recognize_node(operand, variable, context);
            if (!child) return child;
            if (!child.value()) {
                return Result<OptionalRationalPolynomial>::success(std::nullopt);
            }
            const std::size_t projected_terms = result.is_zero() || child.value()->is_zero()
                ? 0
                : result.coeffs.size() + child.value()->coeffs.size() - 1;
            if (projected_terms > context.limits().max_expansion_terms) {
                return Result<OptionalRationalPolynomial>::failure(
                    CasErrc::ResourceLimit,
                    "polynomial expansion term budget exhausted",
                    kOperation);
            }
            result = result * *child.value();
        }
        return Result<OptionalRationalPolynomial>::success(std::move(result));
}

Result<OptionalRationalPolynomial> recognize_power(
    const PowerNode& power, const std::string& variable, ComputationContext& context) {
        auto exponent_node = std::dynamic_pointer_cast<const NumberNode>(power.exponent());
        if (!exponent_node || std::holds_alternative<lmmc_real_t>(exponent_node->value())) {
            return Result<OptionalRationalPolynomial>::success(std::nullopt);
        }

        BigInt exponent;
        if (std::holds_alternative<BigInt>(exponent_node->value())) {
            exponent = std::get<BigInt>(exponent_node->value());
        } else {
            const Rational& rational = std::get<Rational>(exponent_node->value());
            if (!rational.is_integer()) {
                return Result<OptionalRationalPolynomial>::success(std::nullopt);
            }
            exponent = rational.to_bigint();
        }

        auto exponent_value = exponent.try_to_int64();
        if (!exponent_value || *exponent_value <= 0) {
            return Result<OptionalRationalPolynomial>::success(std::nullopt);
        }
        if (static_cast<std::uint64_t>(*exponent_value) >=
            context.limits().max_expansion_terms) {
            return Result<OptionalRationalPolynomial>::failure(
                CasErrc::ResourceLimit,
                "polynomial exponent exceeds the expansion term budget",
                kOperation);
        }

        auto base = recognize_node(power.base(), variable, context);
        if (!base) return base;
        if (!base.value()) {
            return Result<OptionalRationalPolynomial>::success(std::nullopt);
        }

        Polynomial<Rational> result({Rational(1)}, variable);
        for (std::int64_t i = 0; i < *exponent_value; ++i) {
            auto step = context.consume_steps(1, kOperation);
            if (!step) return Result<OptionalRationalPolynomial>::failure(step.error());
            const std::size_t projected_terms = result.is_zero() || base.value()->is_zero()
                ? 0
                : result.coeffs.size() + base.value()->coeffs.size() - 1;
            if (projected_terms > context.limits().max_expansion_terms) {
                return Result<OptionalRationalPolynomial>::failure(
                    CasErrc::ResourceLimit,
                    "polynomial expansion term budget exhausted",
                    kOperation);
            }
            result = result * *base.value();
        }
        return Result<OptionalRationalPolynomial>::success(std::move(result));
}

Result<OptionalRationalPolynomial> recognize_node(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& variable,
    ComputationContext& context) {
    auto entered = context.enter_recursion(kOperation);
    if (!entered) return Result<OptionalRationalPolynomial>::failure(entered.error());
    RecursionScope scope(context);

    if (!node) {
        return Result<OptionalRationalPolynomial>::failure(
            CasErrc::InternalInvariant,
            "polynomial recognition encountered a null AST node",
            kOperation);
    }

    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return recognize_number(*number, variable);
    }
    if (auto symbol = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (symbol->name() != variable) {
            return Result<OptionalRationalPolynomial>::success(std::nullopt);
        }
        return Result<OptionalRationalPolynomial>::success(
            Polynomial<Rational>({Rational(0), Rational(1)}, variable));
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return recognize_sum(*add, variable, context);
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return recognize_product(*multiply, variable, context);
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return recognize_power(*power, variable, context);
    }
    return Result<OptionalRationalPolynomial>::success(std::nullopt);
}

}

Result<OptionalRationalPolynomial> recognize_rational_polynomial(
    const SymbolicExpr& expression,
    const std::string& variable,
    ComputationContext& context) {
    if (!LMCAS::detail::node(expression)) {
        return Result<OptionalRationalPolynomial>::failure(
            CasErrc::InvalidArgument, "expression cannot be null", kOperation);
    }
    if (variable.empty()) {
        return Result<OptionalRationalPolynomial>::failure(
            CasErrc::InvalidArgument, "polynomial variable cannot be empty", kOperation);
    }
    return recognize_node(LMCAS::detail::node(expression), variable, context);
}

}
