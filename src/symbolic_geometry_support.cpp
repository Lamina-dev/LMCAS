#include "internal/symbolic_geometry_support.hpp"
#include "internal/symbolic_ast.hpp"
#include "numeric_evaluation.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace LMCAS::geometry_detail {

Result<void> validate_symbolic_vector(
    const std::vector<std::shared_ptr<SymbolicExpr>>& vector,
    const std::string& name,
    const std::string& operation)
{
    for (size_t i = 0; i < vector.size(); ++i) {
        if (!vector[i] || !LMCAS::detail::node(vector[i])) {
            return Result<void>::failure(
                CasErrc::InvalidArgument,
                name + " contains a null component at index " + std::to_string(i),
                operation);
        }
    }
    return Result<void>::success();
}

Result<void> validate_same_dimension_vectors(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (a.empty() || b.empty()) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, "vectors cannot be empty", operation);
    }
    if (a.size() != b.size()) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, "vector dimensions do not match", operation);
    }
    auto a_check = validate_symbolic_vector(a, "left vector", operation);
    if (!a_check) return a_check;
    return validate_symbolic_vector(b, "right vector", operation);
}

Result<void> validate_surface_point(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!surf.F || !LMCAS::detail::node(surf.F)) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, "surface equation cannot be null", operation);
    }
    if (surf.vars.empty()) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, "surface variables cannot be empty", operation);
    }
    if (surf.vars.size() != point.size()) {
        return Result<void>::failure(
            CasErrc::InvalidArgument,
            "surface variable count must match point dimension",
            operation);
    }
    for (const auto& var : surf.vars) {
        if (var.empty()) {
            return Result<void>::failure(
                CasErrc::InvalidArgument, "surface variable names cannot be empty",
                operation);
        }
    }
    return validate_symbolic_vector(point, "surface point", operation);
}

ExpressionResult simplify_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& operation,
    const std::string& message)
{
    if (!expr || !LMCAS::detail::node(expr)) {
        return ExpressionResult::failure(CasErrc::Inconclusive, message, operation);
    }
    auto simplified = expr->simplify();
    if (!simplified || !LMCAS::detail::node(simplified)) {
        return ExpressionResult::failure(CasErrc::Inconclusive, message, operation);
    }
    return ExpressionResult::success(std::move(simplified));
}

Result<double> numeric_vector_scale_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& expressions,
    ComputationContext& context,
    const std::string& operation,
    std::vector<double>* values)
{
    double max_abs = 0.0;
    if (values) {
        values->clear();
        values->reserve(expressions.size());
    }
    for (const auto& expression : expressions) {
        if (!expression || !LMCAS::detail::node(expression)) {
            return Result<double>::failure(
                CasErrc::InvalidArgument,
                "numeric vector contains a null expression",
                operation);
        }
        auto evaluated = evaluate_numeric(
            *expression, NumericBindings{}, context);
        if (!evaluated) {
            if (evaluated.error().code == CasErrc::Cancelled ||
                evaluated.error().code == CasErrc::ResourceLimit) {
                return Result<double>::failure(evaluated.error());
            }
            return Result<double>::failure(
                CasErrc::NumericFailure,
                "vector is not finite numeric in the supported domain",
                operation);
        }
        if (!evaluated.value().is_finite() ||
            !std::isfinite(evaluated.value().value)) {
            return Result<double>::failure(
                CasErrc::NumericFailure,
                "vector is not finite numeric in the supported domain",
                operation);
        }
        max_abs = std::max(max_abs, std::abs(evaluated.value().value));
        if (values) {
            values->push_back(evaluated.value().value);
        }
    }
    return Result<double>::success(max_abs);
}

Result<void> checked_nonzero_numeric_vector(
    const std::vector<std::shared_ptr<SymbolicExpr>>& expressions,
    ComputationContext& context,
    const std::string& operation,
    const std::string& domain_message,
    const std::string& inconclusive_message)
{
    auto scaled = numeric_vector_scale_checked(
        expressions, context, operation);
    if (!scaled) {
        if (scaled.error().code == CasErrc::Cancelled ||
            scaled.error().code == CasErrc::ResourceLimit) {
            return Result<void>::failure(scaled.error());
        }
        return Result<void>::failure(
            CasErrc::Inconclusive, inconclusive_message, operation);
    }
    if (scaled.value() == 0.0) {
        return Result<void>::failure(
            CasErrc::DomainError, domain_message, operation);
    }
    return Result<void>::success();
}

Result<void> checked_nonzero_numeric_or_exact(
    const std::shared_ptr<SymbolicExpr>& expr,
    ComputationContext& context,
    const std::string& operation,
    const std::string& domain_message,
    const std::string& inconclusive_message)
{
    if (!expr || !LMCAS::detail::node(expr)) {
        return Result<void>::failure(CasErrc::Inconclusive,
                                     inconclusive_message, operation);
    }
    auto simplified = expr->simplify();
    if (!simplified || !LMCAS::detail::node(simplified)) {
        return Result<void>::failure(CasErrc::Inconclusive,
                                     inconclusive_message, operation);
    }
    if (simplified->is_zero()) {
        return Result<void>::failure(CasErrc::DomainError,
                                     domain_message, operation);
    }
    auto numeric = evaluate_numeric(*simplified, NumericBindings{}, context);
    if (!numeric) {
        if (numeric.error().code == CasErrc::Cancelled ||
            numeric.error().code == CasErrc::ResourceLimit) {
            return Result<void>::failure(numeric.error());
        }
        return Result<void>::failure(CasErrc::Inconclusive,
                                     inconclusive_message, operation);
    }
    if (!numeric.value().is_finite() ||
        !std::isfinite(numeric.value().value)) {
        return Result<void>::failure(CasErrc::NumericFailure,
                                     inconclusive_message, operation);
    }
    if (numeric.value().value == 0.0) {
        return Result<void>::failure(CasErrc::DomainError,
                                     domain_message, operation);
    }
    return Result<void>::success();
}

Result<void> validate_line_checked(const LineSymbolic& line,
                                   ComputationContext& context,
                                   const std::string& operation)
{
    auto valid = validate_same_dimension_vectors(line.point, line.direction,
                                                 context, operation);
    if (!valid) return valid;
    if (line.point.size() != 3) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "line geometry requires 3-dimensional data",
                                     operation);
    }
    return checked_nonzero_numeric_vector(
        line.direction, context, operation,
        "line direction cannot be zero",
        "line direction nonzero condition cannot be verified");
}

Result<void> validate_plane_checked(const PlaneSymbolic& plane,
                                    ComputationContext& context,
                                    const std::string& operation)
{
    if (!plane.d || !LMCAS::detail::node(plane.d)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "plane constant cannot be null",
                                     operation);
    }
    if (plane.normal.size() != 3) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "plane geometry requires a 3-dimensional normal",
                                     operation);
    }
    auto normal_valid = validate_symbolic_vector(plane.normal, "plane normal",
                                                 operation);
    if (!normal_valid) return normal_valid;
    return checked_nonzero_numeric_vector(
        plane.normal, context, operation,
        "plane normal cannot be zero",
        "plane normal nonzero condition cannot be verified");
}

}
