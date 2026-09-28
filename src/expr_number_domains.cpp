#include "expr.hpp"

#include <cmath>
#include <optional>
#include <utility>
#include <vector>

#include "internal/expression_analysis.hpp"
#include "internal/assumption_simplification.hpp"
#include "internal/symbolic_ast.hpp"
#include "root_of_identity.hpp"

namespace LMCAS {
namespace {

constexpr const char* kNumberDomainOperation = "LMCAS.number_domain";

Result<bool> bool_failure(CasErrc code, std::string message,
                          const char* operation) {
    return Result<bool>::failure(code, std::move(message), operation);
}

Result<bool> inconclusive_domain() {
    return bool_failure(CasErrc::Inconclusive,
                        "domain membership is only decidable for numeric and explicit complex expressions",
                        kNumberDomainOperation);
}


int domain_rank(NumberDomain domain) noexcept {
    switch (domain) {
    case NumberDomain::Integers: {
        return 0;
    }
    case NumberDomain::Rationals: {
        return 1;
    }
    case NumberDomain::Reals: {
        return 2;
    }
    case NumberDomain::Complexes: {
        return 3;
    }
    case NumberDomain::Expressions: {
        return 4;
    }
    }
    return -1;
}

Result<bool> domain_contains_node(
    NumberDomain domain, const std::shared_ptr<const SymbolicNode>& node);

Result<bool> domain_contains_number(NumberDomain domain, const NumberNode& number) {
    const auto& value = number.value();
    if (std::holds_alternative<BigInt>(value)) {
        return Result<bool>::success(true);
    }
    if (std::holds_alternative<Rational>(value)) {
        const auto& rational = std::get<Rational>(value);
        if (domain == NumberDomain::Integers) {
            return Result<bool>::success(rational.is_integer());
        }
        if (domain_rank(domain) >= 0) {
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
    case NumberDomain::Rationals: {
        return Result<bool>::success(false);
    }
    case NumberDomain::Reals:
    case NumberDomain::Complexes:
    case NumberDomain::Expressions: {
        return Result<bool>::success(true);
    }
    }
    return inconclusive_domain();
}

Result<bool> domain_contains_complex(NumberDomain domain, const ComplexNode& node) {
    if (domain == NumberDomain::Complexes) {
        auto real = domain_contains_node(NumberDomain::Reals, node.real());
        if (!real) {
            return real;
        }
        auto imag = domain_contains_node(NumberDomain::Reals, node.imag());
        if (!imag) {
            return imag;
        }
        return Result<bool>::success(real.value() && imag.value());
    }
    if (!node.imag()->is_zero()) {
        return Result<bool>::success(false);
    }
    return domain_contains_node(domain, node.real());
}

Result<bool> domain_contains_operands(
    NumberDomain domain,
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands) {
    bool all_members = true;
    for (const auto& operand : operands) {
        auto member = domain_contains_node(domain, operand);
        if (!member) {
            if (member.error().code != CasErrc::Inconclusive) {
                return member;
            }
            all_members = false;
            continue;
        }
        all_members = all_members && member.value();
    }
    return all_members ? Result<bool>::success(true) : inconclusive_domain();
}

bool exact_sqrt_contains(NumberDomain domain, const Rational& value) {
    if (value < Rational(0)) {
        return domain == NumberDomain::Complexes;
    }
    if (domain == NumberDomain::Reals || domain == NumberDomain::Complexes) {
        return true;
    }
    const bool rational_root = value.get_numerator().is_perfect_square() &&
                               value.get_denominator().is_perfect_square();
    return domain == NumberDomain::Rationals
        ? rational_root
        : rational_root && value.get_denominator() == BigInt(1);
}

Result<bool> domain_contains_fractional_power(
    NumberDomain domain, const PowerNode& power, const Rational& exponent) {
    auto base = detail::exact_rational_value(power.base());
    if (exponent == Rational(1, 2) && base) {
        return Result<bool>::success(exact_sqrt_contains(domain, *base));
    }
    return bool_failure(
        CasErrc::Inconclusive,
        "non-integer power membership requires a branch/domain proof",
        kNumberDomainOperation);
}

bool is_undefined_zero_power(
    const Rational& exponent,
    const std::shared_ptr<const NumberNode>& base) {
    return exponent <= Rational(0) && base && base->is_zero();
}

Result<bool> domain_contains_negative_integer_power(
    NumberDomain domain,
    const std::shared_ptr<const NumberNode>& base_number,
    const std::optional<Rational>& exact_base) {
    if (exact_base) {
        if (domain == NumberDomain::Integers) {
            return Result<bool>::success(
                exact_base->get_numerator().abs() == BigInt(1));
        }
        return Result<bool>::success(
            domain == NumberDomain::Rationals ||
            domain == NumberDomain::Reals ||
            domain == NumberDomain::Complexes);
    }
    if (base_number) {
        return domain_contains_number(domain, *base_number);
    }
    return inconclusive_domain();
}

Result<bool> domain_contains_positive_integer_power(
    NumberDomain domain, const PowerNode& power) {
    auto base = domain_contains_node(domain, power.base());
    if (!base) {
        return base.error().code == CasErrc::Inconclusive
            ? inconclusive_domain()
            : Result<bool>::failure(base.error());
    }
    return base.value() ? Result<bool>::success(true) : inconclusive_domain();
}

Result<bool> domain_contains_power(NumberDomain domain, const PowerNode& power) {
    const auto exponent = detail::exact_rational_value(power.exponent());
    if (!exponent) {
        return bool_failure(
            CasErrc::Inconclusive,
            "domain membership for powers requires an exact exponent",
            kNumberDomainOperation);
    }
    if (!exponent->is_integer()) {
        return domain_contains_fractional_power(domain, power, *exponent);
    }

    const auto base_number =
        std::dynamic_pointer_cast<const NumberNode>(power.base());
    const auto exact_base = detail::exact_rational_value(power.base());
    if (is_undefined_zero_power(*exponent, base_number)) {
        return bool_failure(
            CasErrc::DomainError,
            "zero to a nonpositive integer power is undefined",
            kNumberDomainOperation);
    }
    if (exponent->is_zero()) {
        return base_number ? Result<bool>::success(true)
                           : inconclusive_domain();
    }
    if (*exponent < Rational(0)) {
        return domain_contains_negative_integer_power(
            domain, base_number, exact_base);
    }
    return domain_contains_positive_integer_power(domain, power);
}

Result<bool> domain_contains_root(
    NumberDomain domain, const std::shared_ptr<const SymbolicNode>& node) {
    auto expression = detail::make_expression_ptr(node);
    if (domain == NumberDomain::Complexes) {
        auto value = rootof_evaluate_complex_checked(expression);
        return value ? Result<bool>::success(true)
                     : Result<bool>::failure(value.error());
    }
    ComputationContext context;
    auto value = rootof_evaluate_checked(expression, context);
    if (value) {
        return Result<bool>::success(domain == NumberDomain::Reals);
    }
    if (value.error().code == CasErrc::DomainError) {
        return Result<bool>::success(false);
    }
    return Result<bool>::failure(value.error());
}

bool function_preserves_domain(FunctionNode::FuncType type) {
    switch (type) {
        case FunctionNode::FuncType::Sin: {
            return true;
        }
        case FunctionNode::FuncType::Cos: {
            return true;
        }
        case FunctionNode::FuncType::Exp: {
            return true;
        }
        case FunctionNode::FuncType::Sinh: {
            return true;
        }
        case FunctionNode::FuncType::Cosh: {
            return true;
        }
        case FunctionNode::FuncType::Tanh: {
            return true;
        }
        case FunctionNode::FuncType::Abs: {
            return true;
        }
        default:
            return false;
    }
}

Result<bool> domain_contains_function(NumberDomain domain, const FunctionNode& function) {
    if (function.arguments().size() != 1) {
        return inconclusive_domain();
    }
    if (function.type() == FunctionNode::FuncType::Sqrt) {
        auto exact = detail::exact_rational_value(function.arguments()[0]);
        if (exact) {
            if (domain_rank(domain) >= 0) {
                return Result<bool>::success(exact_sqrt_contains(domain, *exact));
            }
            if (*exact < Rational(0)) {
                return Result<bool>::success(false);
            }
        }
    }
    if (function_preserves_domain(function.type())) {
        if (domain == NumberDomain::Integers ||
            domain == NumberDomain::Rationals) {
            return inconclusive_domain();
        }
        return domain_contains_node(
            domain == NumberDomain::Complexes ? NumberDomain::Complexes
                                               : NumberDomain::Reals,
            function.arguments()[0]);
    }
    return inconclusive_domain();
}

Result<bool> domain_contains_node(
    NumberDomain domain, const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) {
        return bool_failure(CasErrc::InvalidArgument,
                            "domain membership element must not be null",
                            kNumberDomainOperation);
    }
    if (domain == NumberDomain::Expressions) {
        return Result<bool>::success(true);
    }
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return domain_contains_number(domain, *number);
    }
    if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        return domain_contains_complex(domain, *complex);
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return domain_contains_operands(domain, add->operands());
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return domain_contains_operands(domain, multiply->operands());
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return domain_contains_power(domain, *power);
    }
    if (std::dynamic_pointer_cast<const RootOfNode>(node)) {
        return domain_contains_root(domain, node);
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return domain_contains_function(domain, *function);
    }
    return inconclusive_domain();
}

}

const char* NumberDomainSet::name() const noexcept {
    switch (domain_) {
    case NumberDomain::Integers: {
        return "Z";
    }
    case NumberDomain::Rationals: {
        return "Q";
    }
    case NumberDomain::Reals: {
        return "R";
    }
    case NumberDomain::Complexes: {
        return "C";
    }
    case NumberDomain::Expressions: {
        return "Expr";
    }
    }
    return "?";
}

bool NumberDomainSet::subset_of(const NumberDomainSet& other) const noexcept {
    return domain_rank(domain_) <= domain_rank(other.domain_);
}

Result<bool> NumberDomainSet::contains(const ConstExprPtr& element) const {
    if (!element || !detail::node(*element)) {
        return bool_failure(CasErrc::InvalidArgument,
                            "domain membership element must not be null",
                            kNumberDomainOperation);
    }
    if (domain_ == NumberDomain::Expressions) {
        return Result<bool>::success(true);
    }

    ComputationContext context;
    auto normalized = detail::simplify_expression(
        detail::make_expression_ptr(detail::node(*element)), context);
    if (!normalized) {
        return Result<bool>::failure(normalized.error());
    }
    return domain_contains_node(domain_, detail::node(normalized.value()));
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

Result<bool> domain_contains(const NumberDomainSet& domain, const ExprPtr& element) {
    return domain.contains(element);
}

Result<bool> domain_subset(const NumberDomainSet& lhs, const NumberDomainSet& rhs) {
    return Result<bool>::success(lhs.subset_of(rhs));
}

}
