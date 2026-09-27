#include "polynomial_conversion.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/facts_query.hpp"
#include "internal/normalization_utils.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/visitors/normalization_visitor.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>

namespace LMCAS {

template <>
Result<SymbolicPolyCoeff> extract_coeff_value<SymbolicPolyCoeff>(
    const std::shared_ptr<SymbolicExpr>& coefficient) {
    if (!coefficient || !detail::node(coefficient)) {
        return Result<SymbolicPolyCoeff>::failure(CasErrc::InvalidArgument,
            "Coefficient must have an expression", "polynomial.coefficient");
    }
    return SymbolicPolyCoeff(coefficient);
}

template <>
Result<BigInt> extract_coeff_value<BigInt>(
    const std::shared_ptr<SymbolicExpr>& coefficient) {
    if (!coefficient || !detail::node(coefficient)) {
        return Result<BigInt>::failure(CasErrc::InvalidArgument,
            "Coefficient must have an expression", "polynomial.coefficient");
    }
    auto simplified = coefficient->simplify();
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(
            LMCAS::detail::node(simplified))) {
        if (std::holds_alternative<BigInt>(number->value())) {
            return std::get<BigInt>(number->value());
        }
        if (std::holds_alternative<Rational>(number->value())) {
            const Rational& value = std::get<Rational>(number->value());
            if (value.is_integer()) return value.to_bigint();
        }
        if (std::holds_alternative<lmmc_real_t>(number->value())) {
            const lmmc_real_t value = std::get<lmmc_real_t>(number->value());
            if (std::isfinite(value) && value == std::floor(value)) {
                return Rational::from_double(value).to_bigint();
            }
        }
    }
    return Result<BigInt>::failure(CasErrc::UnsupportedExpression,
        "Coefficient is not an integer", "polynomial.coefficient");
}

template <>
Result<Rational> extract_coeff_value<Rational>(
    const std::shared_ptr<SymbolicExpr>& coefficient) {
    if (!coefficient || !detail::node(coefficient)) {
        return Result<Rational>::failure(CasErrc::InvalidArgument,
            "Coefficient must have an expression", "polynomial.coefficient");
    }
    auto simplified = coefficient->simplify();
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(
            LMCAS::detail::node(simplified))) {
        if (std::holds_alternative<Rational>(number->value())) {
            return std::get<Rational>(number->value());
        }
        if (std::holds_alternative<BigInt>(number->value())) {
            return Rational(std::get<BigInt>(number->value()));
        }
        if (std::holds_alternative<lmmc_real_t>(number->value())) {
            const auto value = std::get<lmmc_real_t>(number->value());
            if (std::isfinite(value)) return Rational::from_double(value);
        }
    }
    if (simplified->is_one()) return Rational(1);
    return Result<Rational>::failure(CasErrc::UnsupportedExpression,
        "Coefficient is not representable as a rational", "polynomial.coefficient");
}

bool contains(const SymbolicExpr& expression, const std::string& variable) {
    return expression_depends_on_variable(
        LMCAS::detail::node(expression), variable);
}

template <typename T>
bool polynomial_product_fits(const Polynomial<T>& left, const Polynomial<T>& right) {
    if (left.is_zero() || right.is_zero()) return true;
    const auto max_terms = static_cast<std::size_t>(std::numeric_limits<int>::max());
    return left.coeffs.size() <= max_terms &&
        right.coeffs.size() - 1 <= max_terms - left.coeffs.size();
}

namespace {

void checked_coefficient_count(std::size_t count, detail::RewriteBudget& budget) {
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw CasError{CasErrc::ResourceLimit,
            "Polynomial coefficient storage limit exceeded", "polynomial.convert"};
    }
    auto terms = budget.context().require_expansion_terms(count, "polynomial.convert");
    if (!terms) { throw terms.error(); }
    budget.require_nodes(count);
}

std::size_t checked_product_count(std::size_t left, std::size_t right,
                                 detail::RewriteBudget& budget) {
    if (left == 0 || right == 0) return 0;
    const auto maximum = static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (left > maximum || right - 1 > maximum - left) {
        throw CasError{CasErrc::ResourceLimit,
            "Polynomial coefficient storage limit exceeded", "polynomial.convert"};
    }
    const auto count = left + right - 1;
    checked_coefficient_count(count, budget);
    return count;
}

detail::SymbolicNodePtr normalize_coefficient(
    const detail::SymbolicNodePtr& node, detail::RewriteBudget& budget) {
    auto step = budget.context().consume_steps(1, "polynomial.coefficient");
    if (!step) { throw step.error(); }
    NormalizationVisitor visitor(budget.context(), detail::no_facts(), Domain::Real, &budget);
    node->accept(visitor);
    return visitor.get_result();
}

template <typename Node>
detail::SymbolicNodePtr coefficient_arithmetic(
    const detail::SymbolicNodePtr& left, const detail::SymbolicNodePtr& right,
    detail::RewriteBudget& budget) {
    normalization_check_arithmetic<Node>(&budget, 0, left, right);
    return normalize_coefficient(detail::make_node<Node>(
        std::vector<detail::SymbolicNodePtr>{left, right}), budget);
}

void trim_checked_coefficients(Polynomial<SymbolicPolyCoeff>& polynomial,
                               detail::RewriteBudget& budget) {
    while (!polynomial.coeffs.empty()) {
        auto step = budget.context().consume_steps(1, "polynomial.coefficient");
        if (!step) { throw step.error(); }
        if (!detail::node(polynomial.coeffs.back().val)->is_zero()) break;
        polynomial.coeffs.pop_back();
    }
}

Polynomial<SymbolicPolyCoeff> add_checked_coefficients(
    const Polynomial<SymbolicPolyCoeff>& left, const Polynomial<SymbolicPolyCoeff>& right,
    detail::RewriteBudget& budget) {
    const auto count = std::max(left.coeffs.size(), right.coeffs.size());
    checked_coefficient_count(count, budget);
    Polynomial<SymbolicPolyCoeff> result(left.variable_name);
    result.coeffs.reserve(count);
    std::size_t nodes = 0;
    for (std::size_t index = 0; index < count; ++index) {
        detail::SymbolicNodePtr coefficient;
        if (index < left.coeffs.size() && index < right.coeffs.size()) {
            coefficient = coefficient_arithmetic<AddNode>(
                detail::node(left.coeffs[index].val), detail::node(right.coeffs[index].val), budget);
        } else {
            const auto& source = index < left.coeffs.size() ? left : right;
            coefficient = detail::node(source.coeffs[index].val);
        }
        nodes = budget.append_size(nodes, budget.measure(coefficient));
        result.coeffs.emplace_back(detail::make_expression_ptr(std::move(coefficient)));
    }
    trim_checked_coefficients(result, budget);
    return result;
}

Polynomial<SymbolicPolyCoeff> multiply_checked_coefficients(
    const Polynomial<SymbolicPolyCoeff>& left, const Polynomial<SymbolicPolyCoeff>& right,
    detail::RewriteBudget& budget) {
    const auto count = checked_product_count(left.coeffs.size(), right.coeffs.size(), budget);
    Polynomial<SymbolicPolyCoeff> result(left.variable_name);
    result.coeffs.reserve(count);
    std::size_t nodes = 0;
    for (std::size_t degree = 0; degree < count; ++degree) {
        detail::SymbolicNodePtr coefficient;
        std::size_t coefficient_nodes = 0;
        const auto first = degree < right.coeffs.size() ? 0 : degree - right.coeffs.size() + 1;
        const auto last = std::min(degree, left.coeffs.size() - 1);
        for (std::size_t index = first; index <= last; ++index) {
            auto product = coefficient_arithmetic<MultiplyNode>(
                detail::node(left.coeffs[index].val),
                detail::node(right.coeffs[degree - index].val), budget);
            coefficient = coefficient
                ? coefficient_arithmetic<AddNode>(coefficient, product, budget)
                : std::move(product);
            coefficient_nodes = budget.measure(coefficient);
            budget.append_size(nodes, coefficient_nodes);
        }
        nodes = budget.append_size(nodes, coefficient_nodes);
        result.coeffs.emplace_back(detail::make_expression_ptr(std::move(coefficient)));
    }
    trim_checked_coefficients(result, budget);
    return result;
}

Result<void> check_zero_power_base(
    const detail::SymbolicNodePtr& base, ComputationContext& context) {
    const auto& facts = detail::no_facts();
    auto real_defined = detail::query_definedness(base, facts, Domain::Real, context);
    if (!real_defined) { return Result<void>::failure(real_defined.error()); }
    auto complex_defined = detail::query_definedness(base, facts, Domain::Complex, context);
    if (!complex_defined) { return Result<void>::failure(complex_defined.error()); }
    auto nonzero = detail::query_nonzero_value(base, facts, Domain::Real, context);
    if (!nonzero) { return Result<void>::failure(nonzero.error()); }
    if (real_defined.value() != Tribool::True ||
        complex_defined.value() != Tribool::True || nonzero.value() != Tribool::True) {
        return Result<void>::failure(CasErrc::UnsupportedExpression,
            "Zero power requires a defined nonzero base", "polynomial.convert");
    }
    return Result<void>::success();
}

Result<void> check_power_coefficient_count(
    std::size_t base_count, int exponent, detail::RewriteBudget* budget) {
    if (!budget || base_count == 0) { return Result<void>::success(); }
    const auto degree = base_count - 1;
    const auto maximum = static_cast<std::size_t>(std::numeric_limits<int>::max());
    const auto count = static_cast<std::size_t>(exponent);
    if (degree > (maximum - 1) / count) {
        return Result<void>::failure(CasErrc::ResourceLimit,
            "Polynomial coefficient storage limit exceeded", "polynomial.convert");
    }
    checked_coefficient_count(degree * count + 1, *budget);
    return Result<void>::success();
}

template <typename T>
Result<Polynomial<T>> constant_to_poly(
    const detail::SymbolicNodePtr& node, const std::string& variable,
    detail::RewriteBudget* budget) {
    if constexpr (std::is_same_v<T, SymbolicPolyCoeff>) {
        if (budget) {
            checked_coefficient_count(1, *budget);
            auto coefficient = normalize_coefficient(node, *budget);
            budget->append_size(0, budget->measure(coefficient));
            Polynomial<T> result(variable);
            result.coeffs.emplace_back(detail::make_expression_ptr(std::move(coefficient)));
            trim_checked_coefficients(result, *budget);
            return result;
        }
    }
    auto coefficient = extract_coeff_value<T>(detail::make_expression_ptr(node));
    if (!coefficient) { return Result<Polynomial<T>>::failure(coefficient.error()); }
    return Polynomial<T>({std::move(coefficient.value())}, variable);
}

} // namespace

template <typename T>
Result<Polynomial<T>> symbolic_to_poly_recursive(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& variable, ComputationContext* context, detail::RewriteBudget* budget);

template <typename T>
Result<Polynomial<T>> symbolic_power_to_poly(
    const PowerNode& power, const std::string& variable,
    ComputationContext* context, detail::RewriteBudget* budget) {
    detail::RewriteScope scope(budget);
    BigInt exponent;
    if (!try_get_integer_value(
            std::dynamic_pointer_cast<const NumberNode>(power.exponent()), exponent) ||
        exponent.is_negative()) {
        return Result<Polynomial<T>>::failure(CasErrc::UnsupportedExpression,
            "Polynomial powers require nonnegative integer exponents", "polynomial.convert");
    }
    const auto bounded_exponent = exponent.try_to_int64();
    if (!bounded_exponent || *bounded_exponent >= 1000) {
        return Result<Polynomial<T>>::failure(CasErrc::ResourceLimit,
            "Polynomial power expansion limit exceeded", "polynomial.convert");
    }
    const int exponent_value = static_cast<int>(*bounded_exponent);

    Polynomial<T> result(variable);
    if (budget) {
        checked_coefficient_count(1, *budget);
        normalization_check_children(budget, 1);
    }
    result.coeffs.emplace_back(1);
    if (exponent_value == 0) {
        auto defined = check_zero_power_base(power.base(), *context);
        if (!defined) { return Result<Polynomial<T>>::failure(defined.error()); }
        return result;
    }
    const auto base = symbolic_to_poly_recursive<T>(power.base(), variable, context, budget);
    if (!base) { return Result<Polynomial<T>>::failure(base.error()); }
    auto count = check_power_coefficient_count(base.value().coeffs.size(), exponent_value, budget);
    if (!count) { return Result<Polynomial<T>>::failure(count.error()); }
    for (int index = 0; index < exponent_value; ++index) {
        if constexpr (std::is_same_v<T, SymbolicPolyCoeff>) {
            if (budget) {
                result = multiply_checked_coefficients(result, base.value(), *budget);
                continue;
            }
        }
        if (!polynomial_product_fits(result, base.value())) {
            return Result<Polynomial<T>>::failure(CasErrc::ResourceLimit,
                "Polynomial coefficient storage limit exceeded", "polynomial.convert");
        }
        result = result * base.value();
    }
    return result;
}

template <typename T>
Result<Polynomial<T>> symbolic_sum_to_poly(
    const AddNode& add, const std::string& variable,
    ComputationContext* context, detail::RewriteBudget* budget) {
    detail::RewriteScope scope(budget);
    Polynomial<T> result(variable);
    for (const auto& operand : add.operands()) {
        auto child = symbolic_to_poly_recursive<T>(operand, variable, context, budget);
        if (!child) { return Result<Polynomial<T>>::failure(child.error()); }
        if constexpr (std::is_same_v<T, SymbolicPolyCoeff>) {
            if (budget) {
                result = add_checked_coefficients(result, child.value(), *budget);
                continue;
            }
        }
        result = result + child.value();
    }
    return result;
}

template <typename T>
Result<Polynomial<T>> symbolic_product_to_poly(
    const MultiplyNode& multiply, const std::string& variable,
    ComputationContext* context, detail::RewriteBudget* budget) {
    detail::RewriteScope scope(budget);
    Polynomial<T> result(variable);
    if (budget) {
        checked_coefficient_count(1, *budget);
        normalization_check_children(budget, 1);
    }
    result.coeffs.emplace_back(1);
    for (const auto& operand : multiply.operands()) {
        auto child = symbolic_to_poly_recursive<T>(operand, variable, context, budget);
        if (!child) { return Result<Polynomial<T>>::failure(child.error()); }
        if constexpr (std::is_same_v<T, SymbolicPolyCoeff>) {
            if (budget) {
                result = multiply_checked_coefficients(result, child.value(), *budget);
                continue;
            }
        }
        if (!polynomial_product_fits(result, child.value())) {
            return Result<Polynomial<T>>::failure(CasErrc::ResourceLimit,
                "Polynomial coefficient storage limit exceeded", "polynomial.convert");
        }
        result = result * child.value();
    }
    return result;
}

template <typename T>
Result<Polynomial<T>> symbolic_to_poly_recursive(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& variable, ComputationContext* context, detail::RewriteBudget* budget) {
    detail::RewriteScope scope(budget);
    if (!node) {
        return Result<Polynomial<T>>::failure(CasErrc::InvalidArgument,
            "Expression node must not be null", "polynomial.convert");
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        BigInt exponent;
        if (try_get_integer_value(
                std::dynamic_pointer_cast<const NumberNode>(power->exponent()), exponent) &&
            exponent.is_zero()) {
            return symbolic_power_to_poly<T>(*power, variable, context, budget);
        }
    }
    if (!expression_depends_on_variable(node, variable, budget)) {
        return constant_to_poly<T>(node, variable, budget);
    }

    if (auto symbol = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (symbol->name() == variable) {
            if constexpr (std::is_same_v<T, SymbolicPolyCoeff>) {
                if (budget) {
                    checked_coefficient_count(2, *budget);
                    normalization_check_children(budget, 2);
                    Polynomial<T> result(variable);
                    result.coeffs.reserve(2);
                    result.coeffs.emplace_back(0);
                    result.coeffs.emplace_back(1);
                    return result;
                }
            }
            return Polynomial<T>({T(0), T(1)}, variable);
        }
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return symbolic_sum_to_poly<T>(*add, variable, context, budget);
    }

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return symbolic_product_to_poly<T>(*multiply, variable, context, budget);
    }

    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return symbolic_power_to_poly<T>(*power, variable, context, budget);
    }

    return Result<Polynomial<T>>::failure(CasErrc::UnsupportedExpression,
        "Expression is not a polynomial in the requested variable", "polynomial.convert");
}

template <typename T>
LMCAS_API Result<Polynomial<T>> symbolic_to_poly(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable) {
    if (!expression || !LMCAS::detail::node(expression) || variable.empty()) {
        return Result<Polynomial<T>>::failure(CasErrc::InvalidArgument,
            "Expression and variable must not be empty", "polynomial.convert");
    }
    try {
        ComputationContext context;
        return symbolic_to_poly_recursive<T>(LMCAS::detail::node(expression), variable, &context, nullptr);
    } catch (const CasError& error) {
        return Result<Polynomial<T>>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<Polynomial<T>>::failure(CasErrc::ResourceLimit,
            "Polynomial allocation failed", "polynomial.convert");
    } catch (const std::length_error& error) {
        return Result<Polynomial<T>>::failure(CasErrc::ResourceLimit,
            error.what(), "polynomial.convert");
    }
}

Result<Polynomial<SymbolicPolyCoeff>> detail::symbolic_to_poly_checked(
    const SymbolicExpr& expression, const std::string& variable, ComputationContext& context) {
    using ConversionResult = Result<Polynomial<SymbolicPolyCoeff>>;
    if (!detail::node(expression) || variable.empty()) {
        return ConversionResult::failure(CasErrc::InvalidArgument,
            "Expression and variable must not be empty", "polynomial.convert");
    }
    try {
        auto step = context.consume_steps(1, "polynomial.convert");
        if (!step) return ConversionResult::failure(step.error());
        detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
            context.limits().max_ast_nodes, "polynomial.convert");
        budget.measure(detail::node(expression));
        return symbolic_to_poly_recursive<SymbolicPolyCoeff>(
            detail::node(expression), variable, &context, &budget);
    } catch (const CasError& error) {
        return ConversionResult::failure(error);
    } catch (const std::bad_alloc&) {
        return ConversionResult::failure(CasErrc::ResourceLimit,
            "Polynomial allocation failed", "polynomial.convert");
    } catch (const std::length_error& error) {
        return ConversionResult::failure(CasErrc::ResourceLimit,
            error.what(), "polynomial.convert");
    }
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
                variable, SymbolicExpr::number(BigInt(static_cast<std::uint64_t>(degree))));
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

template LMCAS_API Result<Polynomial<BigInt>> symbolic_to_poly<BigInt>(
    const std::shared_ptr<SymbolicExpr>&, const std::string&);
template LMCAS_API Result<Polynomial<Rational>> symbolic_to_poly<Rational>(
    const std::shared_ptr<SymbolicExpr>&, const std::string&);
template LMCAS_API Result<Polynomial<SymbolicPolyCoeff>> symbolic_to_poly<SymbolicPolyCoeff>(
    const std::shared_ptr<SymbolicExpr>&, const std::string&);

template LMCAS_API std::shared_ptr<SymbolicExpr> poly_to_symbolic<BigInt>(
    const Polynomial<BigInt>&);
template LMCAS_API std::shared_ptr<SymbolicExpr> poly_to_symbolic<Rational>(
    const Polynomial<Rational>&);
template LMCAS_API std::shared_ptr<SymbolicExpr> poly_to_symbolic<SymbolicPolyCoeff>(
    const Polynomial<SymbolicPolyCoeff>&);

} // namespace LMCAS
