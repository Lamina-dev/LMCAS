#pragma once

#include "assumption_context.hpp"
#include "integration.hpp"
#include "numeric_probe.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
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

LMCAS_API Result<SymbolicExpr> apply_assumption_simplifications(
    const SymbolicExpr& expr, const std::string& var,
    const AssumptionContext* context);

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

inline std::optional<double> try_checked_numeric_constant(
    const SymbolicExpr& expr,
    ComputationContext& context) {
    return detail::try_finite_numeric(expr, &context);
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

} // namespace LMCAS
