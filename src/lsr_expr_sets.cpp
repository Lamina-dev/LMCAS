#include "lsr_expr.hpp"

#include <cmath>
#include <exception>
#include <utility>
#include <vector>

#include "lsr_expr_internal.hpp"
#include "symbolic_ast.hpp"
#include "root_of_utils.hpp"

namespace lamina::lsr {
namespace {

constexpr const char* kNumberDomainOperation = "lsr.number_domain";

Result<bool> bool_failure(CasErrc code, std::string message,
                          const char* operation) {
    return Result<bool>::failure(code, std::move(message), operation);
}

std::optional<Rational> exact_number_value(
    const std::shared_ptr<const SymbolicNode>& node) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number || std::holds_alternative<lmmc_real_t>(number->value())) {
        return std::nullopt;
    }
    if (std::holds_alternative<BigInt>(number->value())) {
        return Rational(std::get<BigInt>(number->value()));
    }
    return std::get<Rational>(number->value());
}

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

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& operand : add->operands()) {
            auto member = domain_contains_node(domain, operand);
            if (!member || !member.value()) return member;
        }
        return Result<bool>::success(true);
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& operand : multiply->operands()) {
            auto member = domain_contains_node(domain, operand);
            if (!member || !member.value()) return member;
        }
        return Result<bool>::success(true);
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent = std::dynamic_pointer_cast<const NumberNode>(
            power->exponent());
        if (!exponent || std::holds_alternative<lmmc_real_t>(exponent->value())) {
            return bool_failure(
                CasErrc::Inconclusive,
                "domain membership for powers requires an exact exponent",
                kNumberDomainOperation);
        }
        Rational exponent_value = std::holds_alternative<BigInt>(exponent->value())
            ? Rational(std::get<BigInt>(exponent->value()))
            : std::get<Rational>(exponent->value());
        if (!exponent_value.is_integer()) {
            auto exact_base = exact_number_value(power->base());
            if (exponent_value == Rational(1, 2) && exact_base) {
                if (*exact_base < Rational(0)) {
                    return Result<bool>::success(
                        domain == NumberDomain::Complexes);
                }
                if (domain == NumberDomain::Reals ||
                    domain == NumberDomain::Complexes) {
                    return Result<bool>::success(true);
                }
                const bool rational_root =
                    exact_base->get_numerator().is_perfect_square() &&
                    exact_base->get_denominator().is_perfect_square();
                return Result<bool>::success(
                    domain == NumberDomain::Rationals
                        ? rational_root
                        : rational_root &&
                              exact_base->get_denominator() == BigInt(1));
            }
            return bool_failure(
                CasErrc::Inconclusive,
                "non-integer power membership requires a branch/domain proof",
                kNumberDomainOperation);
        }
        auto base_member = domain_contains_node(domain, power->base());
        if (!base_member || !base_member.value()) return base_member;
        if (domain == NumberDomain::Integers &&
            exponent_value < Rational(0)) {
            return Result<bool>::success(false);
        }
        return Result<bool>::success(true);
    }

    if (std::dynamic_pointer_cast<const RootOfNode>(node)) {
        auto expression = lamina::detail::make_expression_ptr(node);
        if (domain == NumberDomain::Complexes) {
            auto value = rootof_evaluate_complex_checked(expression);
            return value ? Result<bool>::success(true)
                         : Result<bool>::failure(value.error());
        }
        ComputationContext root_context;
        auto value = rootof_evaluate_checked(expression, root_context);
        if (value) return Result<bool>::success(domain == NumberDomain::Reals);
        if (value.error().code == CasErrc::DomainError) {
            return Result<bool>::success(false);
        }
        return Result<bool>::failure(value.error());
    }

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (function->type() == FunctionNode::FuncType::Sqrt &&
            function->arguments().size() == 1) {
            auto exact = exact_number_value(function->arguments()[0]);
            if (exact) {
                const bool negative = *exact < Rational(0);
                if (domain == NumberDomain::Complexes) {
                    return Result<bool>::success(true);
                }
                if (negative) return Result<bool>::success(false);
                if (domain == NumberDomain::Reals) {
                    return Result<bool>::success(true);
                }
                const bool rational_root =
                    exact->get_numerator().is_perfect_square() &&
                    exact->get_denominator().is_perfect_square();
                if (domain == NumberDomain::Rationals) {
                    return Result<bool>::success(rational_root);
                }
                if (domain == NumberDomain::Integers) {
                    return Result<bool>::success(
                        rational_root && exact->get_denominator() == BigInt(1));
                }
            }
        }
        if (function->arguments().size() == 1 &&
            (function->type() == FunctionNode::FuncType::Sin ||
             function->type() == FunctionNode::FuncType::Cos ||
             function->type() == FunctionNode::FuncType::Exp ||
             function->type() == FunctionNode::FuncType::Sinh ||
             function->type() == FunctionNode::FuncType::Cosh ||
             function->type() == FunctionNode::FuncType::Tanh ||
             function->type() == FunctionNode::FuncType::Abs)) {
            if (domain == NumberDomain::Integers ||
                domain == NumberDomain::Rationals) {
                return Result<bool>::success(false);
            }
            return domain_contains_node(
                domain == NumberDomain::Complexes
                    ? NumberDomain::Complexes
                    : NumberDomain::Reals,
                function->arguments()[0]);
        }
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

} // namespace

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

} // namespace lamina::lsr
