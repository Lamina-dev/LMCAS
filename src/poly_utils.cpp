#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "symbolic_ast.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
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
        if (std::holds_alternative<BigInt>(number->value())) {
            return Result<OptionalRationalPolynomial>::success(Polynomial<Rational>(
                Rational(std::get<BigInt>(number->value())), variable));
        }
        if (std::holds_alternative<Rational>(number->value())) {
            return Result<OptionalRationalPolynomial>::success(Polynomial<Rational>(
                std::get<Rational>(number->value()), variable));
        }
        return Result<OptionalRationalPolynomial>::success(std::nullopt);
    }

    if (auto symbol = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (symbol->name() != variable) {
            return Result<OptionalRationalPolynomial>::success(std::nullopt);
        }
        return Result<OptionalRationalPolynomial>::success(
            Polynomial<Rational>({Rational(0), Rational(1)}, variable));
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        Polynomial<Rational> result(variable);
        for (const auto& operand : add->operands()) {
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

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        Polynomial<Rational> result({Rational(1)}, variable);
        for (const auto& operand : multiply->operands()) {
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

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent_node = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
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
            exponent = rational.to_BigInt();
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

        auto base = recognize_node(power->base(), variable, context);
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

    return Result<OptionalRationalPolynomial>::success(std::nullopt);
}

} // namespace

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

template <>
BigInt extract_coeff_value<BigInt>(
    const std::shared_ptr<SymbolicExpr>& coefficient) {
    auto simplified = coefficient->simplify();
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(
            LMCAS::detail::node(simplified))) {
        if (std::holds_alternative<BigInt>(number->value())) {
            return std::get<BigInt>(number->value());
        }
        if (std::holds_alternative<Rational>(number->value())) {
            const Rational& value = std::get<Rational>(number->value());
            return value.is_integer() ? value.to_BigInt() : BigInt(0);
        }
        const lmmc_real_t value = std::get<lmmc_real_t>(number->value());
        if (std::isfinite(value) && value == std::floor(value) &&
            value >= static_cast<lmmc_real_t>(
                std::numeric_limits<long long>::min()) &&
            value <= static_cast<lmmc_real_t>(
                std::numeric_limits<long long>::max())) {
            return BigInt(static_cast<long long>(value));
        }
    }
    if (simplified->is_one()) return BigInt(1);
    return BigInt(0);
}

template <>
Rational extract_coeff_value<Rational>(
    const std::shared_ptr<SymbolicExpr>& coefficient) {
    auto simplified = coefficient->simplify();
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(
            LMCAS::detail::node(simplified))) {
        if (std::holds_alternative<Rational>(number->value())) {
            return std::get<Rational>(number->value());
        }
        if (std::holds_alternative<BigInt>(number->value())) {
            return Rational(std::get<BigInt>(number->value()));
        }
        return Rational::from_double(std::get<lmmc_real_t>(number->value()));
    }
    if (simplified->is_one()) return Rational(1);
    return Rational(0);
}

bool contains(const SymbolicExpr& expression, const std::string& variable) {
    return expression_depends_on_variable(
        LMCAS::detail::node(expression), variable);
}

template <typename T>
Polynomial<T> symbolic_to_poly_recursive(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& variable) {
    if (!node) return Polynomial<T>(variable);

    if (!expression_depends_on_variable(node, variable)) {
        return Polynomial<T>({extract_coeff_value<T>(
            LMCAS::detail::make_expression_ptr(node))}, variable);
    }

    if (auto symbol = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (symbol->name() == variable) {
            return Polynomial<T>({T(0), T(1)}, variable);
        }
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        Polynomial<T> result(variable);
        for (const auto& operand : add->operands()) {
            result = result + symbolic_to_poly_recursive<T>(operand, variable);
        }
        return result;
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        Polynomial<T> result({T(1)}, variable);
        for (const auto& operand : multiply->operands()) {
            result = result * symbolic_to_poly_recursive<T>(operand, variable);
        }
        return result;
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
        bool is_nonnegative_integer = false;
        int exponent_value = 0;
        if (exponent && std::holds_alternative<BigInt>(exponent->value())) {
            const auto& value = std::get<BigInt>(exponent->value());
            const auto converted = value.try_to_int64();
            if (converted && *converted >= 0 && *converted < 1000) {
                is_nonnegative_integer = true;
                exponent_value = static_cast<int>(*converted);
            }
        } else if (exponent &&
                   std::holds_alternative<Rational>(exponent->value())) {
            const auto& value = std::get<Rational>(exponent->value());
            const auto converted = value.is_integer()
                ? value.to_BigInt().try_to_int64()
                : std::optional<std::int64_t>{};
            if (converted && *converted >= 0 && *converted < 1000) {
                is_nonnegative_integer = true;
                exponent_value = static_cast<int>(*converted);
            }
        } else if (exponent) {
            const lmmc_real_t value = std::get<lmmc_real_t>(exponent->value());
            if (std::isfinite(value) && value >= 0.0 && value < 1000.0 &&
                value == std::floor(value)) {
                is_nonnegative_integer = true;
                exponent_value = static_cast<int>(value);
            }
        }

        if (is_nonnegative_integer) {
            Polynomial<T> result({T(1)}, variable);
            if (exponent_value == 0) return result;
            const auto base = symbolic_to_poly_recursive<T>(power->base(), variable);
            for (int index = 0; index < exponent_value; ++index) {
                result = result * base;
            }
            return result;
        }
    }

    return Polynomial<T>(variable);
}

template <typename T>
LMCAS_API Polynomial<T> symbolic_to_poly(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable) {
    if (!expression || !LMCAS::detail::node(expression)) {
        return Polynomial<T>(variable);
    }
    return symbolic_to_poly_recursive<T>(
        LMCAS::detail::node(expression), variable);
}

template <typename T>
LMCAS_API std::shared_ptr<SymbolicExpr> poly_to_symbolic(
    const Polynomial<T>& polynomial) {
    if (polynomial.is_zero()) return SymbolicExpr::number(0);

    std::vector<std::shared_ptr<SymbolicExpr>> terms;
    for (std::size_t degree = 0; degree < polynomial.coeffs.size(); ++degree) {
        if (polynomial.coeffs[degree] == T(0)) continue;

        std::shared_ptr<SymbolicExpr> coefficient;
        if constexpr (std::is_same_v<T, SymbolicPolyCoeff>) {
            coefficient = polynomial.coeffs[degree].val
                ? polynomial.coeffs[degree].val
                : SymbolicExpr::number(0);
        } else {
            coefficient = SymbolicExpr::number(polynomial.coeffs[degree]);
        }

        if (degree == 0) {
            terms.push_back(coefficient);
            continue;
        }

        auto variable = SymbolicExpr::variable(polynomial.variable_name);
        auto variable_part = degree == 1
            ? variable
            : SymbolicExpr::power(
                variable, SymbolicExpr::number(static_cast<int>(degree)));
        if (polynomial.coeffs[degree] == T(1)) {
            terms.push_back(variable_part);
        } else if (polynomial.coeffs[degree] == T(-1)) {
            terms.push_back(SymbolicExpr::multiply(
                SymbolicExpr::number(-1), variable_part));
        } else {
            terms.push_back(SymbolicExpr::multiply(coefficient, variable_part));
        }
    }

    std::reverse(terms.begin(), terms.end());
    if (terms.empty()) return SymbolicExpr::number(0);
    auto result = terms.front();
    for (std::size_t index = 1; index < terms.size(); ++index) {
        result = SymbolicExpr::add(result, terms[index]);
    }
    return result;
}

template LMCAS_API Polynomial<BigInt> symbolic_to_poly<BigInt>(
    const std::shared_ptr<SymbolicExpr>&, const std::string&);
template LMCAS_API Polynomial<Rational> symbolic_to_poly<Rational>(
    const std::shared_ptr<SymbolicExpr>&, const std::string&);
template LMCAS_API Polynomial<SymbolicPolyCoeff> symbolic_to_poly<SymbolicPolyCoeff>(
    const std::shared_ptr<SymbolicExpr>&, const std::string&);

template LMCAS_API std::shared_ptr<SymbolicExpr> poly_to_symbolic<BigInt>(
    const Polynomial<BigInt>&);
template LMCAS_API std::shared_ptr<SymbolicExpr> poly_to_symbolic<Rational>(
    const Polynomial<Rational>&);
template LMCAS_API std::shared_ptr<SymbolicExpr> poly_to_symbolic<SymbolicPolyCoeff>(
    const Polynomial<SymbolicPolyCoeff>&);

} // namespace LMCAS
