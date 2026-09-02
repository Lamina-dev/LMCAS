#include "poly_utils.hpp"

#include "multivariate_factor.hpp"
#include "symbolic_ast.hpp"

#include <algorithm>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <new>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace lamina {
namespace {

constexpr const char* kOperation = "symbolic_polynomial_gcd";

class RecursionFrame {
public:
    explicit RecursionFrame(ComputationContext& context) : context_(context) {}
    ~RecursionFrame() { context_.leave_recursion(); }

    RecursionFrame(const RecursionFrame&) = delete;
    RecursionFrame& operator=(const RecursionFrame&) = delete;

private:
    ComputationContext& context_;
};

Result<void> collect_polynomial_variables(
    const std::shared_ptr<const SymbolicNode>& node,
    std::set<std::string>& variables,
    ComputationContext& context) {
    auto entered = context.enter_recursion(kOperation);
    if (!entered) return entered;
    RecursionFrame frame(context);

    if (!node) {
        return Result<void>::failure(
            CasErrc::InternalInvariant, "polynomial contains a null node",
            kOperation);
    }
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        if (std::holds_alternative<lmmc_real_t>(number->value())) {
            return Result<void>::failure(
                CasErrc::UnsupportedExpression,
                "approximate coefficients are not supported by exact GCD",
                kOperation);
        }
        return Result<void>::success();
    }
    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        variables.insert(variable->name());
        return Result<void>::success();
    }
    if (auto addition = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& operand : addition->operands()) {
            auto result = collect_polynomial_variables(operand, variables, context);
            if (!result) return result;
        }
        return Result<void>::success();
    }
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& operand : product->operands()) {
            auto result = collect_polynomial_variables(operand, variables, context);
            if (!result) return result;
        }
        return Result<void>::success();
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
        if (!exponent || std::holds_alternative<lmmc_real_t>(exponent->value())) {
            return Result<void>::failure(
                CasErrc::UnsupportedExpression,
                "polynomial exponents must be exact nonnegative integers",
                kOperation);
        }
        return collect_polynomial_variables(power->base(), variables, context);
    }
    return Result<void>::failure(
        CasErrc::UnsupportedExpression,
        "expression is outside the exact rational polynomial domain",
        kOperation);
}

Result<std::size_t> polynomial_exponent(
    const std::shared_ptr<const SymbolicNode>& node) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number || std::holds_alternative<lmmc_real_t>(number->value())) {
        return Result<std::size_t>::failure(
            CasErrc::UnsupportedExpression,
            "polynomial exponents must be exact nonnegative integers",
            kOperation);
    }

    BigInt integer;
    if (std::holds_alternative<BigInt>(number->value())) {
        integer = std::get<BigInt>(number->value());
    } else {
        const Rational& rational = std::get<Rational>(number->value());
        if (!rational.is_integer()) {
            return Result<std::size_t>::failure(
                CasErrc::UnsupportedExpression,
                "polynomial exponents must be integers", kOperation);
        }
        integer = rational.to_BigInt();
    }

    const auto value = integer.try_to_int64();
    if (!value || *value < 0) {
        return Result<std::size_t>::failure(
            CasErrc::UnsupportedExpression,
            "polynomial exponents must be representable nonnegative integers",
            kOperation);
    }
    return Result<std::size_t>::success(static_cast<std::size_t>(*value));
}

Result<void> check_term_budget(std::size_t lhs_terms,
                               std::size_t rhs_terms,
                               ComputationContext& context) {
    const std::size_t limit = context.limits().max_expansion_terms;
    if (lhs_terms != 0 && rhs_terms > limit / lhs_terms) {
        return Result<void>::failure(
            CasErrc::ResourceLimit,
            "polynomial expansion term budget exhausted", kOperation);
    }
    return context.consume_steps(lhs_terms * rhs_terms + 1, kOperation);
}

Result<MultiPoly> multiply_checked(const MultiPoly& lhs,
                                   const MultiPoly& rhs,
                                   ComputationContext& context) {
    auto budget = check_term_budget(
        static_cast<std::size_t>(lhs.num_terms()),
        static_cast<std::size_t>(rhs.num_terms()), context);
    if (!budget) return Result<MultiPoly>::failure(budget.error());
    MultiPoly result = lhs * rhs;
    if (static_cast<std::size_t>(result.num_terms()) >
        context.limits().max_expansion_terms) {
        return Result<MultiPoly>::failure(
            CasErrc::ResourceLimit,
            "polynomial expansion term budget exhausted", kOperation);
    }
    return Result<MultiPoly>::success(std::move(result));
}

Result<MultiPoly> node_to_multivariate_polynomial(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::vector<std::string>& variables,
    ComputationContext& context) {
    auto entered = context.enter_recursion(kOperation);
    if (!entered) return Result<MultiPoly>::failure(entered.error());
    RecursionFrame frame(context);

    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        if (std::holds_alternative<BigInt>(number->value())) {
            return Result<MultiPoly>::success(MultiPoly(
                Rational(std::get<BigInt>(number->value())), variables));
        }
        if (std::holds_alternative<Rational>(number->value())) {
            return Result<MultiPoly>::success(MultiPoly(
                std::get<Rational>(number->value()), variables));
        }
        return Result<MultiPoly>::failure(
            CasErrc::UnsupportedExpression,
            "approximate coefficients are not supported by exact GCD",
            kOperation);
    }
    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        Monomial monomial(variables.size(), 0);
        auto position = std::find(variables.begin(), variables.end(),
                                  variable->name());
        if (position == variables.end()) {
            return Result<MultiPoly>::failure(
                CasErrc::InternalInvariant,
                "polynomial variable collection is inconsistent", kOperation);
        }
        monomial[static_cast<std::size_t>(position - variables.begin())] = 1;
        return Result<MultiPoly>::success(MultiPoly(
            {{std::move(monomial), Rational(1)}}, variables));
    }
    if (auto addition = std::dynamic_pointer_cast<const AddNode>(node)) {
        MultiPoly result(Rational(0), variables);
        for (const auto& operand : addition->operands()) {
            auto converted = node_to_multivariate_polynomial(
                operand, variables, context);
            if (!converted) return converted;
            auto step = context.consume_steps(
                static_cast<std::size_t>(result.num_terms() +
                                         converted.value().num_terms()) + 1,
                kOperation);
            if (!step) return Result<MultiPoly>::failure(step.error());
            result = result + converted.value();
            if (static_cast<std::size_t>(result.num_terms()) >
                context.limits().max_expansion_terms) {
                return Result<MultiPoly>::failure(
                    CasErrc::ResourceLimit,
                    "polynomial expansion term budget exhausted", kOperation);
            }
        }
        return Result<MultiPoly>::success(std::move(result));
    }
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        MultiPoly result(Rational(1), variables);
        for (const auto& operand : product->operands()) {
            auto converted = node_to_multivariate_polynomial(
                operand, variables, context);
            if (!converted) return converted;
            auto multiplied = multiply_checked(result, converted.value(), context);
            if (!multiplied) return multiplied;
            result = std::move(multiplied.value());
        }
        return Result<MultiPoly>::success(std::move(result));
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent = polynomial_exponent(power->exponent());
        if (!exponent) return Result<MultiPoly>::failure(exponent.error());
        auto base = node_to_multivariate_polynomial(
            power->base(), variables, context);
        if (!base) return base;

        MultiPoly result(Rational(1), variables);
        MultiPoly factor = std::move(base.value());
        std::size_t remaining = exponent.value();
        while (remaining != 0) {
            if ((remaining & 1U) != 0) {
                auto multiplied = multiply_checked(result, factor, context);
                if (!multiplied) return multiplied;
                result = std::move(multiplied.value());
            }
            remaining >>= 1U;
            if (remaining != 0) {
                auto squared = multiply_checked(factor, factor, context);
                if (!squared) return squared;
                factor = std::move(squared.value());
            }
        }
        return Result<MultiPoly>::success(std::move(result));
    }
    return Result<MultiPoly>::failure(
        CasErrc::UnsupportedExpression,
        "expression is outside the exact rational polynomial domain",
        kOperation);
}

MultiPoly make_monic(const MultiPoly& polynomial) {
    if (polynomial.is_zero()) return polynomial;
    const Rational leading = polynomial.terms().front().second;
    return polynomial * (Rational(1) / leading);
}

Result<std::shared_ptr<SymbolicExpr>> multivariate_polynomial_to_expression(
    const MultiPoly& polynomial,
    ComputationContext& context) {
    if (polynomial.is_zero()) {
        auto nodes = context.reserve_nodes(1, kOperation);
        if (!nodes) return Result<std::shared_ptr<SymbolicExpr>>::failure(nodes.error());
        return Result<std::shared_ptr<SymbolicExpr>>::success(
            SymbolicExpr::number(0));
    }

    std::vector<std::shared_ptr<SymbolicExpr>> terms;
    terms.reserve(polynomial.terms().size());
    for (const auto& [monomial, coefficient] : polynomial.terms()) {
        auto step = context.consume_steps(1, kOperation);
        if (!step) return Result<std::shared_ptr<SymbolicExpr>>::failure(step.error());

        std::vector<std::shared_ptr<SymbolicExpr>> factors;
        const bool constant_term = total_degree(monomial) == 0;
        if (coefficient != Rational(1) || constant_term) {
            factors.push_back(SymbolicExpr::number(coefficient));
        }
        for (std::size_t index = 0;
             index < polynomial.variables().size() && index < monomial.size();
             ++index) {
            if (monomial[index] == 0) continue;
            auto variable = SymbolicExpr::variable(polynomial.variables()[index]);
            factors.push_back(monomial[index] == 1
                ? variable
                : SymbolicExpr::power(
                      variable, SymbolicExpr::number(monomial[index])));
        }

        auto reserve = context.reserve_nodes(factors.size() * 2 + 1, kOperation);
        if (!reserve) {
            return Result<std::shared_ptr<SymbolicExpr>>::failure(reserve.error());
        }
        std::shared_ptr<SymbolicExpr> term = factors.empty()
            ? SymbolicExpr::number(1)
            : factors.front();
        for (std::size_t index = 1; index < factors.size(); ++index) {
            term = SymbolicExpr::multiply(term, factors[index]);
        }
        terms.push_back(std::move(term));
    }

    auto reserve = context.reserve_nodes(terms.size(), kOperation);
    if (!reserve) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(reserve.error());
    }
    auto expression = terms.front();
    for (std::size_t index = 1; index < terms.size(); ++index) {
        expression = SymbolicExpr::add(expression, terms[index]);
    }
    return Result<std::shared_ptr<SymbolicExpr>>::success(std::move(expression));
}

} // namespace

SymbolicGcdResult symbolic_polynomial_gcd(
    const SymbolicExpr& lhs,
    const SymbolicExpr& rhs,
    ComputationContext& context) {
    try {
        std::set<std::string> variable_set;
        auto lhs_variables = collect_polynomial_variables(
            detail::node(lhs), variable_set, context);
        if (!lhs_variables) return SymbolicGcdResult::failure(lhs_variables.error());
        auto rhs_variables = collect_polynomial_variables(
            detail::node(rhs), variable_set, context);
        if (!rhs_variables) return SymbolicGcdResult::failure(rhs_variables.error());
        std::vector<std::string> variables(variable_set.begin(), variable_set.end());

        auto lhs_polynomial = node_to_multivariate_polynomial(
            detail::node(lhs), variables, context);
        if (!lhs_polynomial) return SymbolicGcdResult::failure(lhs_polynomial.error());
        auto rhs_polynomial = node_to_multivariate_polynomial(
            detail::node(rhs), variables, context);
        if (!rhs_polynomial) return SymbolicGcdResult::failure(rhs_polynomial.error());

        auto step = context.consume_steps(
            static_cast<std::size_t>(lhs_polynomial.value().num_terms() +
                                     rhs_polynomial.value().num_terms()) + 1,
            kOperation);
        if (!step) return SymbolicGcdResult::failure(step.error());

        MultiPoly gcd = make_monic(multivariate_gcd(
            lhs_polynomial.value(), rhs_polynomial.value()));
        if (!gcd.is_zero()) {
            MultiPoly lhs_quotient;
            MultiPoly rhs_quotient;
            try {
                lhs_quotient = lhs_polynomial.value().exact_div(gcd);
                rhs_quotient = rhs_polynomial.value().exact_div(gcd);
            } catch (const std::bad_alloc&) {
                return SymbolicGcdResult::failure(
                    CasErrc::ResourceLimit,
                    "symbolic GCD verification allocation failed", kOperation);
            } catch (const std::runtime_error&) {
                return SymbolicGcdResult::failure(
                    CasErrc::Inconclusive,
                    "GCD candidate failed exact divisibility verification",
                    kOperation);
            } catch (const std::exception& error) {
                return SymbolicGcdResult::failure(
                    CasErrc::InternalInvariant, error.what(), kOperation);
            }
            MultiPoly residual = make_monic(multivariate_gcd(
                lhs_quotient, rhs_quotient));
            if (!residual.is_constant()) {
                return SymbolicGcdResult::failure(
                    CasErrc::Inconclusive,
                    "GCD candidate is not maximal", kOperation);
            }
        }
        return multivariate_polynomial_to_expression(gcd, context);
    } catch (const std::bad_alloc&) {
        return SymbolicGcdResult::failure(
            CasErrc::ResourceLimit, "symbolic GCD allocation failed", kOperation);
    } catch (const std::exception& error) {
        return SymbolicGcdResult::failure(
            CasErrc::InternalInvariant, error.what(), kOperation);
    }
}

Result<Rational> symbolic_polynomial_content(
    const SymbolicExpr& expression,
    ComputationContext& context) {
    try {
        std::set<std::string> variable_set;
        auto variables_result = collect_polynomial_variables(
            detail::node(expression), variable_set, context);
        if (!variables_result) {
            return Result<Rational>::failure(variables_result.error());
        }
        std::vector<std::string> variables(
            variable_set.begin(), variable_set.end());
        auto polynomial = node_to_multivariate_polynomial(
            detail::node(expression), variables, context);
        if (!polynomial) {
            return Result<Rational>::failure(polynomial.error());
        }
        return Result<Rational>::success(polynomial.value().numeric_content());
    } catch (const std::bad_alloc&) {
        return Result<Rational>::failure(
            CasErrc::ResourceLimit, "polynomial content allocation failed",
            kOperation);
    } catch (const std::exception& error) {
        return Result<Rational>::failure(
            CasErrc::InternalInvariant, error.what(), kOperation);
    }
}

} // namespace lamina
