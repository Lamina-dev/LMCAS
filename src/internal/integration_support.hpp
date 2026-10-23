#pragma once

#include "assumption_context.hpp"
#include "integration.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/expression_transform.hpp"
#include "internal/normalization_utils.hpp"
#include "residual_verification.hpp"
#include "polynomial.hpp"
#include "solve_polynomial.hpp"
#include "internal/symbolic_ast.hpp"

#include "lmmc/config.h"
#include "lmmc/numeric.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <optional>
#include <variant>

namespace LMCAS {

inline bool depends_on_integration_variable(
    const SymbolicExpr& expression,
    const std::string& variable) {
    return expression_depends_on_variable(detail::node(expression), variable);
}

inline std::shared_ptr<SymbolicExpr> sym_sub(
    const SymbolicExpr& lhs,
    const SymbolicExpr& rhs) {
    auto negated_rhs = SymbolicExpr::multiply(
        SymbolicExpr::number(-1), detail::make_expression_ptr(rhs));
    return SymbolicExpr::add(detail::make_expression_ptr(lhs), negated_rhs);
}

inline std::shared_ptr<SymbolicExpr> sym_rational(
    const BigInt& numerator, const BigInt& denominator) {
    return SymbolicExpr::number(Rational(numerator, denominator));
}
inline std::optional<int> integration_bounded_exponent(
    const NumberNode& number, int minimum_exponent, int maximum_exponent) {
    const auto narrow = [minimum_exponent, maximum_exponent](const BigInt& exponent)
        -> std::optional<int> {
        if (exponent < BigInt(minimum_exponent) || exponent > BigInt(maximum_exponent)) {
            return std::nullopt;
        }
        const auto bounded_exponent = exponent.try_to_int64();
        if (!bounded_exponent) { return std::nullopt; }
        return static_cast<int>(*bounded_exponent);
    };
    const auto& value = number.value();
    if (const auto* exponent = std::get_if<BigInt>(&value)) { return narrow(*exponent); }
    if (const auto* exponent = std::get_if<Rational>(&value)) {
        if (!exponent->is_integer()) { return std::nullopt; }
        return narrow(exponent->to_bigint());
    }
    if (const auto* exponent = std::get_if<lmmc_real_t>(&value)) {
        if (!std::isfinite(*exponent)) { return std::nullopt; }
        const Rational exact_exponent = Rational::from_double(*exponent);
        if (!exact_exponent.is_integer()) { return std::nullopt; }
        return narrow(exact_exponent.to_bigint());
    }
    return std::nullopt;
}


inline std::shared_ptr<SymbolicExpr> make_unary_function(
    FunctionNode::FuncType type,
    const std::shared_ptr<SymbolicExpr>& operand) {
    return detail::make_expression_ptr(detail::make_node<FunctionNode>(
        type,
        std::vector<std::shared_ptr<const SymbolicNode>>{
            detail::node(operand)}));
}

inline std::shared_ptr<SymbolicExpr> make_arctan(
    const std::shared_ptr<SymbolicExpr>& operand) {
    return make_unary_function(FunctionNode::FuncType::ArcTan, operand);
}

inline bool contains_unevaluated_integral(
    const std::shared_ptr<const SymbolicNode>& node) {
    return detail::contains_node_type<IntegralNode>(node);
}

inline bool certify_integral(
    const SymbolicExpr& candidate, const SymbolicExpr& integrand,
    const std::string& variable, detail::RewriteBudget& budget) {
    budget.measure(detail::node(integrand));
    budget.measure(detail::node(candidate));
    if (contains_unevaluated_integral(detail::node(candidate)) ||
        contains_inexact_number(detail::node(integrand)) ||
        contains_inexact_number(detail::node(candidate))) {
        return false;
    }
    try {
        auto derivative = detail::make_expression_ptr(detail::differentiate_node(
            detail::node(candidate), variable, budget));
        auto residual = sym_sub(*derivative, integrand);
        auto proof = check_zero_residual(residual, budget.context());
        if (!proof) { throw proof.error(); }
        return std::holds_alternative<ProvedZeroResidual>(proof.value());
    } catch (const CasError& error) {
        if (error.code == CasErrc::UnsupportedExpression) { return false; }
        throw;
    }
}

template <typename Transform>
std::shared_ptr<SymbolicExpr> lookup_integral(
    const IntegrationTable& table, const SymbolicExpr& matched,
    const std::string& matched_variable, const SymbolicExpr& integrand,
    const std::string& variable, detail::RewriteBudget& budget, Transform&& transform) {
    for (const auto* entry : table.get_all_sorted()) {
        budget.consume();
        MatchMap bindings;
        if (!Matcher::match(entry->pattern, matched, entry->wildcards, bindings)) { continue; }
        // Rewrite-rule partial matches must not discard integration operands.
        if (bindings.find("__Add_REST__") != bindings.end() ||
            bindings.find("__Mul_REST__") != bindings.end()) {
            continue;
        }
        const bool complete = std::all_of(entry->wildcards.begin(), entry->wildcards.end(),
            [&](const std::string& name) { return bindings.find(name) != bindings.end(); });
        if (!complete || (entry->condition && !entry->condition(bindings, matched_variable))) {
            continue;
        }
        auto instantiated = detail::make_expression_ptr(
            Matcher::replace(entry->result, bindings, false));
        budget.measure(detail::node(instantiated));
        if (contains_inexact_number(detail::node(instantiated)) ||
            contains_unevaluated_integral(detail::node(instantiated))) {
            continue;
        }
        auto candidate = transform(std::move(instantiated));
        if (!candidate || !certify_integral(*candidate, integrand, variable, budget)) { continue; }
        return detail::make_expression_ptr(detail::normalize_node(detail::node(candidate), budget));
    }
    return nullptr;
}

} // namespace LMCAS
