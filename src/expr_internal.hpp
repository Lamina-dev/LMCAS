/** @file expr_internal.hpp */
#pragma once

#include <optional>
#include <string>
#include <utility>

#include "expr.hpp"
#include "symbolic_ast.hpp"

namespace LMCAS {

// Operation tags shared by two or more of the expr_*.cpp translation units.
inline constexpr const char* kEquivalentOperation = "LMCAS.equivalent_core";
inline constexpr const char* kEquivalentProfileOperation =
    "LMCAS.equivalent_core.profile";
inline constexpr const char* kExprSetOperation = "LMCAS.expr_set";
inline constexpr const char* kSolveExprSetOperation = "LMCAS.solve_expr_set";
inline constexpr const char* kEvalComplexOperation = "LMCAS.eval_complex";

inline ExprResult expression_failure(CasErrc code, std::string message,
                                     const char* operation) {
    return ExprResult::failure(code, std::move(message), operation);
}

inline ExprSetResult expr_set_failure(CasErrc code, std::string message,
                                      const char* operation) {
    return ExprSetResult::failure(code, std::move(message), operation);
}

inline bool is_imaginary_unit_name(const std::string& name) {
    return name == "I";
}

inline std::optional<int> exact_small_integer_node(
    const std::shared_ptr<const SymbolicNode>& node,
    int min_value,
    int max_value) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return std::nullopt;

    BigInt value;
    if (std::holds_alternative<BigInt>(number->value())) {
        value = std::get<BigInt>(number->value());
    } else if (std::holds_alternative<Rational>(number->value())) {
        const Rational& rational = std::get<Rational>(number->value());
        if (!rational.is_integer()) return std::nullopt;
        value = rational.to_BigInt();
    } else {
        return std::nullopt;
    }

    if (value < BigInt(min_value) || value > BigInt(max_value)) {
        return std::nullopt;
    }
    return value.to_int();
}

} // namespace LMCAS
