/** @file expr_internal.hpp */
#pragma once

#include <optional>
#include <string>
#include <utility>

#include "expr.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/equivalence_engine.hpp"
#include "internal/expression_analysis.hpp"

namespace LMCAS {
inline constexpr const char* kExprSetOperation = "LMCAS.expr_set";
inline constexpr const char* kSolveExprSetOperation = "LMCAS.solve_expr_set";
inline constexpr const char* kEvalComplexOperation = "LMCAS.eval_complex";

inline ExprSetResult expr_set_failure(CasErrc code, std::string message,
                                      const char* operation) {
    return ExprSetResult::failure(code, std::move(message), operation);
}


}
