#pragma once

#include "assumption_context.hpp"
#include "integration.hpp"
#include "numeric_evaluation.hpp"
#include "poly_utils.hpp"
#include "expression_analysis.hpp"
#include "polynomial.hpp"
#include "solve_polynomial.hpp"
#include "symbolic_ast.hpp"

#include "lmmc/config.h"
#include "lmmc/numeric.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <optional>
#include <variant>

namespace lamina {

inline std::shared_ptr<SymbolicExpr> make_expr_ptr(const SymbolicExpr& expr) {
    return detail::make_expression_ptr(expr);
}

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

inline std::shared_ptr<SymbolicExpr> sym_rational(long long numerator, long long denominator) {
    return SymbolicExpr::number(Rational(BigInt(numerator), BigInt(denominator)));
}

inline std::optional<double> try_checked_numeric_constant(const SymbolicExpr& expr) {
    ComputationContext context;
    auto evaluated = evaluate_numeric(expr, NumericBindings{}, context);
    if (!evaluated || !evaluated.value().is_finite() ||
        !std::isfinite(evaluated.value().value)) {
        return std::nullopt;
    }
    return evaluated.value().value;
}

inline std::shared_ptr<SymbolicExpr> make_arctan(
    const std::shared_ptr<SymbolicExpr>& operand) {
    return detail::make_expression_ptr(detail::make_node<FunctionNode>(
        FunctionNode::FuncType::ArcTan,
        std::vector<std::shared_ptr<const SymbolicNode>>{detail::node(operand)}));
}

inline bool contains_unevaluated_integral(
    const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) {
        return false;
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (function->type() == FunctionNode::FuncType::Calculus_Integral) {
            return true;
        }
        for (const auto& argument : function->arguments()) {
            if (contains_unevaluated_integral(argument)) {
                return true;
            }
        }
    } else if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& operand : add->operands()) {
            if (contains_unevaluated_integral(operand)) {
                return true;
            }
        }
    } else if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& operand : multiply->operands()) {
            if (contains_unevaluated_integral(operand)) {
                return true;
            }
        }
    } else if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return contains_unevaluated_integral(power->base()) ||
               contains_unevaluated_integral(power->exponent());
    }
    return false;
}

} // namespace lamina
