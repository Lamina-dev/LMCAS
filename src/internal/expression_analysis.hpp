#pragma once

#include "internal/symbolic_ast.hpp"
#include "computation_context.hpp"
#include "result.hpp"
#include "polynomial_conversion.hpp"

#include <optional>
#include <set>
#include <string>

namespace LMCAS {
LMCAS_API std::optional<int> exact_small_integer_node(
    const detail::SymbolicNodePtr& node, int min_value, int max_value);

namespace detail {
std::optional<Rational> exact_rational_value(const NumberNode& number);
std::optional<Rational> exact_rational_value(const SymbolicNodePtr& node);

constexpr RelationOp reversed_relation(RelationOp op) noexcept {
    switch (op) {
        case RelationOp::GT: return RelationOp::LT;
        case RelationOp::GEQ: return RelationOp::LEQ;
        case RelationOp::LT: return RelationOp::GT;
        case RelationOp::LEQ: return RelationOp::GEQ;
        default: return op;
    }
}

Result<Polynomial<SymbolicPolyCoeff>> symbolic_to_poly_checked(
    const SymbolicExpr& expression, const std::string& variable,
    ComputationContext& context);

struct AffineForm {
    std::shared_ptr<SymbolicExpr> slope;
    std::shared_ptr<SymbolicExpr> offset;
};
Result<std::optional<AffineForm>> recognize_affine(
    const SymbolicExpr& expression, const std::string& variable,
    ComputationContext& context);

bool is_imaginary_unit_name(const std::string& name);
std::set<std::string> all_variable_names(const SymbolicNodePtr& expression);
}
LMCAS_API std::set<std::string> free_variables(
    const detail::SymbolicNodePtr& expression,
    detail::RewriteBudget* budget = nullptr);
LMCAS_API bool expression_depends_on_variable(
    const detail::SymbolicNodePtr& expression,
    const std::string& variable,
    detail::RewriteBudget* budget = nullptr);
LMCAS_API detail::SymbolicNodePtr substitute_free(
    const detail::SymbolicNodePtr& expression,
    const std::string& variable,
    const detail::SymbolicNodePtr& replacement,
    detail::RewriteBudget* budget = nullptr);

}
