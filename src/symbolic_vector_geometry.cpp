#include "../include/symbolic_vector_geometry.hpp"
#include "symbolic_quadric_classification.hpp"
#include "symbolic_ast.hpp"
#include "../include/numeric_evaluation.hpp"
#include "../include/symbolic.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace LMCAS {

namespace {

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

ExpressionResult differentiate_surface_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const std::string& operation)
{
    try {
        auto derivative = expr->differentiate(var);
        return simplify_checked(
            derivative, operation,
            "surface derivative is outside the supported domain");
    } catch (const std::exception&) {
        return ExpressionResult::failure(
            CasErrc::Inconclusive,
            "surface derivative is outside the supported domain",
            operation);
    }
}

VectorExprListResult surface_gradient_at_point_checked(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    const std::string& operation)
{
    std::vector<std::shared_ptr<SymbolicExpr>> gradient;
    gradient.reserve(surf.vars.size());
    for (size_t i = 0; i < surf.vars.size(); ++i) {
        auto derivative = differentiate_surface_checked(surf.F, surf.vars[i], operation);
        if (!derivative) return VectorExprListResult::failure(derivative.error());
        auto substituted = derivative.value();
        for (size_t j = 0; j < surf.vars.size(); ++j) {
            substituted = substituted->substitute(surf.vars[j], point[j]);
            if (!substituted || !LMCAS::detail::node(substituted)) {
                return VectorExprListResult::failure(
                    CasErrc::Inconclusive,
                    "surface derivative substitution is outside the supported domain",
                    operation);
            }
        }
        auto simplified = simplify_checked(
            substituted, operation,
            "surface derivative substitution is outside the supported domain");
        if (!simplified) return VectorExprListResult::failure(simplified.error());
        gradient.push_back(std::move(simplified.value()));
    }
    return VectorExprListResult::success(std::move(gradient));
}


Result<double> numeric_vector_scale_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& expressions,
    ComputationContext& context,
    const std::string& operation,
    std::vector<double>* values = nullptr)
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

} // namespace

Result<double> symbolic_vector_finite_numeric(
    const std::shared_ptr<SymbolicExpr>& expr,
    ComputationContext& context,
    const std::string& operation)
{
    if (!expr || !LMCAS::detail::node(expr)) {
        return Result<double>::failure(
            CasErrc::InvalidArgument,
            "numeric expression cannot be null",
            operation);
    }
    auto evaluated = evaluate_numeric(*expr, NumericBindings{}, context);
    if (!evaluated) {
        if (evaluated.error().code == CasErrc::Cancelled ||
            evaluated.error().code == CasErrc::ResourceLimit) {
            return Result<double>::failure(evaluated.error());
        }
        return Result<double>::failure(
            CasErrc::NumericFailure,
            "expression is not finite numeric in the supported domain",
            operation);
    }
    if (!evaluated.value().is_finite() ||
        !std::isfinite(evaluated.value().value)) {
        return Result<double>::failure(
            CasErrc::NumericFailure,
            "expression is not finite numeric in the supported domain",
            operation);
    }
    return Result<double>::success(evaluated.value().value);
}


std::shared_ptr<SymbolicExpr> vector_dot(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
) {
    auto checked = vector_dot_checked(a, b);
    if (!checked) {
        throw std::invalid_argument("vector_dot: " + checked.error().message);
    }
    return checked.value();
}

ExpressionResult vector_dot_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b,
    ComputationContext& context
) {
    const std::string operation = "vector_dot";
    auto input = validate_same_dimension_vectors(a, b, context, operation);
    if (!input) return ExpressionResult::failure(input.error());
    std::vector<std::shared_ptr<const SymbolicNode>> sum_terms;
    for (size_t i = 0; i < a.size(); ++i) {
        auto step = context.consume_steps(1, operation);
        if (!step) return ExpressionResult::failure(step.error());
        sum_terms.push_back(LMCAS::detail::node(SymbolicExpr::multiply(a[i], b[i])));
    }
    return ExpressionResult::success(
        LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<AddNode>(sum_terms)));
}

ExpressionResult vector_dot_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
) {
    ComputationContext context;
    return vector_dot_checked(a, b, context);
}

std::vector<std::shared_ptr<SymbolicExpr>> vector_cross(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
) {
    auto checked = vector_cross_checked(a, b);
    if (!checked) {
        throw std::invalid_argument("vector_cross: " + checked.error().message);
    }
    return checked.value();
}

VectorExprListResult vector_cross_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b,
    ComputationContext& context
) {
    const std::string operation = "vector_cross";
    auto step = context.consume_steps(1, operation);
    if (!step) return VectorExprListResult::failure(step.error());
    if (a.size() != 3 || b.size() != 3) {
        return VectorExprListResult::failure(
            CasErrc::InvalidArgument,
            "cross product requires 3-dimensional vectors",
            operation);
    }
    auto a_check = validate_symbolic_vector(a, "left vector", operation);
    if (!a_check) return VectorExprListResult::failure(a_check.error());
    auto b_check = validate_symbolic_vector(b, "right vector", operation);
    if (!b_check) return VectorExprListResult::failure(b_check.error());
    auto x = SymbolicExpr::add(SymbolicExpr::multiply(a[1], b[2]), SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(a[2], b[1])));
    auto y = SymbolicExpr::add(SymbolicExpr::multiply(a[2], b[0]), SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(a[0], b[2])));
    auto z = SymbolicExpr::add(SymbolicExpr::multiply(a[0], b[1]), SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::multiply(a[1], b[0])));
    return VectorExprListResult::success({x, y, z});
}

VectorExprListResult vector_cross_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
) {
    ComputationContext context;
    return vector_cross_checked(a, b, context);
}


VectorAngleResult vector_angle_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b,
    ComputationContext& context
) {
    const std::string operation = "vector_angle";
    auto input = validate_same_dimension_vectors(a, b, context, operation);
    if (!input) return VectorAngleResult::failure(input.error());
    double scale_a = 0.0;
    double scale_b = 0.0;
    double norm_a_scaled = 0.0;
    double norm_b_scaled = 0.0;
    double dot_scaled = 0.0;
    double dot_correction = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        auto step = context.consume_steps(1, operation);
        if (!step) return VectorAngleResult::failure(step.error());
        auto na_result = symbolic_vector_finite_numeric(a[i], context, operation);
        if (!na_result) return VectorAngleResult::failure(na_result.error());
        auto nb_result = symbolic_vector_finite_numeric(b[i], context, operation);
        if (!nb_result) return VectorAngleResult::failure(nb_result.error());
        const double na = na_result.value();
        const double nb = nb_result.value();
        const double abs_a = std::abs(na);
        const double abs_b = std::abs(nb);

        if (abs_a > scale_a) {
            const double ratio = scale_a == 0.0 ? 0.0 : scale_a / abs_a;
            norm_a_scaled *= ratio * ratio;
            dot_scaled *= ratio;
            dot_correction *= ratio;
            scale_a = abs_a;
        }
        if (abs_b > scale_b) {
            const double ratio = scale_b == 0.0 ? 0.0 : scale_b / abs_b;
            norm_b_scaled *= ratio * ratio;
            dot_scaled *= ratio;
            dot_correction *= ratio;
            scale_b = abs_b;
        }
        if (scale_a != 0.0) {
            const double normalized = na / scale_a;
            norm_a_scaled += normalized * normalized;
        }
        if (scale_b != 0.0) {
            const double normalized = nb / scale_b;
            norm_b_scaled += normalized * normalized;
        }
        if (scale_a != 0.0 && scale_b != 0.0) {
            const double product = (na / scale_a) * (nb / scale_b);
            const double corrected = product - dot_correction;
            const double next = dot_scaled + corrected;
            dot_correction = (next - dot_scaled) - corrected;
            dot_scaled = next;
        }
    }
    if (scale_a == 0.0 || scale_b == 0.0) {
        return VectorAngleResult::failure(
            CasErrc::DomainError,
            "angle is undefined for zero-length vectors",
            operation);
    }
    const double denominator =
        std::sqrt(norm_a_scaled) * std::sqrt(norm_b_scaled);
    double cosv = dot_scaled / denominator;
    if (!std::isfinite(cosv)) {
        return VectorAngleResult::failure(
            CasErrc::NumericFailure,
            "vector angle normalization produced a non-finite result",
            operation);
    }
    // Clamp into [-1, 1] to absorb final-rounding excursions.
    if (cosv > 1.0) cosv = 1.0;
    else if (cosv < -1.0) cosv = -1.0;
    return VectorAngleResult::success(std::acos(cosv));
}

VectorAngleResult vector_angle_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b
) {
    ComputationContext context;
    return vector_angle_checked(a, b, context);
}


VectorExprListResult line_plane_intersection_checked(
    const LineSymbolic& line,
    const PlaneSymbolic& plane,
    ComputationContext& context
) {
    const std::string operation = "line_plane_intersection";
    try {
        auto line_valid = validate_line_checked(line, context, operation);
        if (!line_valid) return VectorExprListResult::failure(line_valid.error());
        auto plane_valid = validate_plane_checked(plane, context, operation);
        if (!plane_valid) return VectorExprListResult::failure(plane_valid.error());
        auto step = context.consume_steps(6, operation);
        if (!step) return VectorExprListResult::failure(step.error());

        std::vector<double> normal_values;
        std::vector<double> direction_values;
        auto normal_scale = numeric_vector_scale_checked(
            plane.normal, context, operation, &normal_values);
        if (!normal_scale) {
            return VectorExprListResult::failure(normal_scale.error());
        }
        auto direction_scale = numeric_vector_scale_checked(
            line.direction, context, operation, &direction_values);
        if (!direction_scale) {
            return VectorExprListResult::failure(direction_scale.error());
        }
        std::vector<double> point_values;
        auto point_scale = numeric_vector_scale_checked(
            line.point, context, operation, &point_values);
        const bool numeric_point = point_scale.has_value();
        if (!point_scale &&
            (point_scale.error().code == CasErrc::Cancelled ||
             point_scale.error().code == CasErrc::ResourceLimit)) {
            return VectorExprListResult::failure(point_scale.error());
        }
        auto plane_offset = evaluate_numeric(
            *plane.d, NumericBindings{}, context);
        const bool numeric_plane_offset =
            plane_offset.has_value() &&
            plane_offset.value().is_finite() &&
            std::isfinite(plane_offset.value().value);
        if (!plane_offset &&
            (plane_offset.error().code == CasErrc::Cancelled ||
             plane_offset.error().code == CasErrc::ResourceLimit)) {
            return VectorExprListResult::failure(plane_offset.error());
        }

        const double safe_high =
            std::numeric_limits<double>::max() / 3.0;
        const double safe_low = std::numeric_limits<double>::min();
        const bool scale_direction_products =
            normal_scale.value() > safe_high / direction_scale.value() ||
            normal_scale.value() < safe_low / direction_scale.value();
        const bool protect_point_dot =
            numeric_point && point_scale.value() > 0.0 &&
            normal_scale.value() >
                std::numeric_limits<double>::max() /
                    point_scale.value() / 3.0;

        const std::vector<std::shared_ptr<SymbolicExpr>>* effective_normal =
            &plane.normal;
        const std::vector<std::shared_ptr<SymbolicExpr>>* effective_direction =
            &line.direction;
        std::vector<std::shared_ptr<SymbolicExpr>> scaled_normal;
        std::vector<std::shared_ptr<SymbolicExpr>> scaled_direction;
        std::shared_ptr<SymbolicExpr> effective_d = plane.d;
        if (scale_direction_products || protect_point_dot) {
            scaled_normal.reserve(3);
            for (size_t i = 0; i < 3; ++i) {
                double component =
                    normal_values[i] / normal_scale.value();
                if (protect_point_dot) component /= 3.0;
                scaled_normal.push_back(SymbolicExpr::number(component));
            }
            effective_normal = &scaled_normal;
            effective_d = SymbolicExpr::divide(
                plane.d, SymbolicExpr::number(normal_scale.value()));
            if (protect_point_dot) {
                effective_d = SymbolicExpr::divide(
                    effective_d, SymbolicExpr::number(3.0));
            }
        }
        if (scale_direction_products) {
            scaled_direction.reserve(3);
            for (size_t i = 0; i < 3; ++i) {
                scaled_direction.push_back(SymbolicExpr::number(
                    direction_values[i] / direction_scale.value()));
            }
            effective_direction = &scaled_direction;
        }

        const bool protect_coordinate_update =
            numeric_point && numeric_plane_offset &&
            point_scale.value() >
                std::numeric_limits<double>::max() / 3.0;
        if (protect_coordinate_update) {
            const double plane_location =
                plane_offset.value().value / normal_scale.value();
            if (!std::isfinite(plane_location)) {
                return VectorExprListResult::failure(
                    CasErrc::NumericFailure,
                    "line-plane intersection coordinate is outside the supported domain",
                    operation);
            }
            double spatial_scale = std::max(
                point_scale.value(), std::abs(plane_location));
            if (spatial_scale == 0.0) spatial_scale = 1.0;

            double unit_normal[3];
            double unit_direction[3];
            double scaled_point[3];
            for (size_t i = 0; i < 3; ++i) {
                unit_normal[i] =
                    normal_values[i] / normal_scale.value();
                unit_direction[i] =
                    direction_values[i] / direction_scale.value();
                scaled_point[i] = point_values[i] / spatial_scale;
            }
            const double denominator = std::fma(
                unit_normal[0], unit_direction[0],
                std::fma(
                    unit_normal[1], unit_direction[1],
                    unit_normal[2] * unit_direction[2]));
            if (denominator == 0.0) {
                return VectorExprListResult::failure(
                    CasErrc::DomainError,
                    "line is parallel to plane, so no unique intersection exists",
                    operation);
            }
            const double point_dot = std::fma(
                unit_normal[0], scaled_point[0],
                std::fma(
                    unit_normal[1], scaled_point[1],
                    unit_normal[2] * scaled_point[2]));
            const double scaled_offset =
                plane_location / spatial_scale;
            const double parameter =
                (scaled_offset - point_dot) / denominator;

            std::vector<std::shared_ptr<SymbolicExpr>> intersection;
            intersection.reserve(3);
            for (size_t i = 0; i < 3; ++i) {
                const double scaled_coordinate = std::fma(
                    parameter, unit_direction[i], scaled_point[i]);
                const double coordinate =
                    scaled_coordinate * spatial_scale;
                if (!std::isfinite(coordinate)) {
                    return VectorExprListResult::failure(
                        CasErrc::NumericFailure,
                        "line-plane intersection coordinate is outside the supported domain",
                        operation);
                }
                intersection.push_back(
                    SymbolicExpr::number(coordinate));
            }
            return VectorExprListResult::success(
                std::move(intersection));
        }

        auto n_dot_a = vector_dot_checked(
            *effective_normal, line.point, context);
        if (!n_dot_a) return VectorExprListResult::failure(n_dot_a.error());
        auto n_dot_b = vector_dot_checked(
            *effective_normal, *effective_direction, context);
        if (!n_dot_b) return VectorExprListResult::failure(n_dot_b.error());
        auto denom = n_dot_b.value()->simplify();
        auto nonzero = checked_nonzero_numeric_or_exact(
            denom, context, operation,
            "line is parallel to plane, so no unique intersection exists",
            "line-plane denominator nonzero condition cannot be verified");
        if (!nonzero) return VectorExprListResult::failure(nonzero.error());

        auto t = SymbolicExpr::divide(
            SymbolicExpr::add(
                effective_d,
                SymbolicExpr::multiply(SymbolicExpr::number(-1), n_dot_a.value())),
            denom);
        std::vector<std::shared_ptr<SymbolicExpr>> intersection;
        intersection.reserve(line.point.size());
        for (size_t i = 0; i < line.point.size(); ++i) {
            auto coord = SymbolicExpr::add(
                line.point[i],
                SymbolicExpr::multiply(t, (*effective_direction)[i]));
            auto simplified = simplify_checked(
                coord, operation,
                "line-plane intersection coordinate is outside the supported domain");
            if (!simplified) {
                return VectorExprListResult::failure(simplified.error());
            }
            intersection.push_back(std::move(simplified.value()));
        }
        return VectorExprListResult::success(std::move(intersection));
    } catch (const std::bad_alloc&) {
        return VectorExprListResult::failure(CasErrc::ResourceLimit,
                                             "vector geometry allocation failed",
                                             operation);
    } catch (const std::exception& e) {
        return VectorExprListResult::failure(CasErrc::InternalInvariant,
                                             e.what(), operation);
    }
}

VectorExprListResult line_plane_intersection_checked(
    const LineSymbolic& line,
    const PlaneSymbolic& plane
) {
    ComputationContext context;
    return line_plane_intersection_checked(line, plane, context);
}


ExpressionResult point_plane_distance_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    const PlaneSymbolic& plane,
    ComputationContext& context
) {
    const std::string operation = "point_plane_distance";
    try {
        auto plane_valid = validate_plane_checked(plane, context, operation);
        if (!plane_valid) return ExpressionResult::failure(plane_valid.error());
        if (point.size() != 3) {
            return ExpressionResult::failure(
                CasErrc::InvalidArgument,
                "point-plane distance requires a 3-dimensional point",
                operation);
        }
        auto point_valid = validate_symbolic_vector(point, "point", operation);
        if (!point_valid) return ExpressionResult::failure(point_valid.error());
        auto step = context.consume_steps(4, operation);
        if (!step) return ExpressionResult::failure(step.error());

        std::vector<double> normal_values;
        auto normal_scale = numeric_vector_scale_checked(
            plane.normal, context, operation, &normal_values);
        if (!normal_scale) {
            return ExpressionResult::failure(normal_scale.error());
        }
        auto point_scale = numeric_vector_scale_checked(
            point, context, operation);
        const bool numeric_point = point_scale.has_value();
        if (!point_scale &&
            (point_scale.error().code == CasErrc::Cancelled ||
             point_scale.error().code == CasErrc::ResourceLimit)) {
            return ExpressionResult::failure(point_scale.error());
        }

        const double safe_low =
            std::sqrt(std::numeric_limits<double>::min());
        const double safe_high =
            std::sqrt(std::numeric_limits<double>::max() / 3.0);
        const bool protect_point_dot =
            numeric_point && point_scale.value() > 0.0 &&
            normal_scale.value() >
                std::numeric_limits<double>::max() /
                    point_scale.value() / 3.0;
        const std::vector<std::shared_ptr<SymbolicExpr>>* effective_normal =
            &plane.normal;
        std::vector<std::shared_ptr<SymbolicExpr>> scaled_normal;
        std::shared_ptr<SymbolicExpr> effective_d = plane.d;
        if (normal_scale.value() < safe_low ||
            normal_scale.value() > safe_high ||
            protect_point_dot) {
            scaled_normal.reserve(normal_values.size());
            for (double component : normal_values) {
                double scaled_component =
                    component / normal_scale.value();
                if (protect_point_dot) scaled_component /= 3.0;
                scaled_normal.push_back(
                    SymbolicExpr::number(scaled_component));
            }
            effective_normal = &scaled_normal;
            effective_d = SymbolicExpr::divide(
                plane.d, SymbolicExpr::number(normal_scale.value()));
            if (protect_point_dot) {
                effective_d = SymbolicExpr::divide(
                    effective_d, SymbolicExpr::number(3.0));
            }
        }

        auto n_dot_r = vector_dot_checked(
            *effective_normal, point, context);
        if (!n_dot_r) return ExpressionResult::failure(n_dot_r.error());
        auto norm_sq = vector_dot_checked(
            *effective_normal, *effective_normal, context);
        if (!norm_sq) return ExpressionResult::failure(norm_sq.error());

        auto diff = SymbolicExpr::add(
            n_dot_r.value(),
            SymbolicExpr::multiply(
                SymbolicExpr::number(-1), effective_d));
        auto norm_n = SymbolicExpr::sqrt(norm_sq.value());
        auto abs_diff = LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Abs,
                std::vector<std::shared_ptr<const SymbolicNode>>{
                    LMCAS::detail::node(diff)}));
        return simplify_checked(
            SymbolicExpr::divide(abs_diff, norm_n),
            operation,
            "point-plane distance is outside the supported domain");
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                         "vector geometry allocation failed",
                                         operation);
    } catch (const std::exception& e) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                         e.what(), operation);
    }
}

ExpressionResult point_plane_distance_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    const PlaneSymbolic& plane
) {
    ComputationContext context;
    return point_plane_distance_checked(point, plane, context);
}


ExpressionResult skew_lines_distance_checked(
    const LineSymbolic& l1,
    const LineSymbolic& l2,
    ComputationContext& context
) {
    const std::string operation = "skew_lines_distance";
    try {
        auto l1_valid = validate_line_checked(l1, context, operation);
        if (!l1_valid) return ExpressionResult::failure(l1_valid.error());
        auto l2_valid = validate_line_checked(l2, context, operation);
        if (!l2_valid) return ExpressionResult::failure(l2_valid.error());
        auto step = context.consume_steps(8, operation);
        if (!step) return ExpressionResult::failure(step.error());

        std::vector<double> direction1_values;
        std::vector<double> direction2_values;
        auto scale1 = numeric_vector_scale_checked(
            l1.direction, context, operation, &direction1_values);
        if (!scale1) return ExpressionResult::failure(scale1.error());
        auto scale2 = numeric_vector_scale_checked(
            l2.direction, context, operation, &direction2_values);
        if (!scale2) return ExpressionResult::failure(scale2.error());

        std::vector<std::shared_ptr<SymbolicExpr>> direction1;
        std::vector<std::shared_ptr<SymbolicExpr>> direction2;
        direction1.reserve(3);
        direction2.reserve(3);
        for (size_t i = 0; i < 3; ++i) {
            direction1.push_back(SymbolicExpr::number(
                direction1_values[i] / scale1.value()));
            direction2.push_back(SymbolicExpr::number(
                direction2_values[i] / scale2.value()));
        }

        std::vector<double> point1_values;
        std::vector<double> point2_values;
        auto point1_scale = numeric_vector_scale_checked(
            l1.point, context, operation, &point1_values);
        auto point2_scale = numeric_vector_scale_checked(
            l2.point, context, operation, &point2_values);
        const bool numeric_points =
            point1_scale.has_value() && point2_scale.has_value();
        if (!point1_scale &&
            (point1_scale.error().code == CasErrc::Cancelled ||
             point1_scale.error().code == CasErrc::ResourceLimit)) {
            return ExpressionResult::failure(point1_scale.error());
        }
        if (!point2_scale &&
            (point2_scale.error().code == CasErrc::Cancelled ||
             point2_scale.error().code == CasErrc::ResourceLimit)) {
            return ExpressionResult::failure(point2_scale.error());
        }

        bool scale_point_difference = false;
        double point_scale = 0.0;
        if (numeric_points) {
            point_scale = std::max(
                point1_scale.value(), point2_scale.value());
            for (size_t i = 0; i < 3; ++i) {
                if (!std::isfinite(point2_values[i] - point1_values[i])) {
                    scale_point_difference = true;
                    break;
                }
            }
        }

        std::vector<std::shared_ptr<SymbolicExpr>> a2_minus_a1;
        a2_minus_a1.reserve(3);
        for (size_t i = 0; i < l1.point.size(); ++i) {
            if (scale_point_difference) {
                a2_minus_a1.push_back(SymbolicExpr::number(
                    point2_values[i] / point_scale -
                    point1_values[i] / point_scale));
                continue;
            }
            a2_minus_a1.push_back(SymbolicExpr::add(
                l2.point[i],
                SymbolicExpr::multiply(
                    SymbolicExpr::number(-1), l1.point[i])));
        }
        auto cross = vector_cross_checked(direction1, direction2, context);
        if (!cross) return ExpressionResult::failure(cross.error());
        auto cross_nonzero = checked_nonzero_numeric_vector(
            cross.value(), context, operation,
            "skew line distance requires non-parallel directions",
            "cross-product nonzero condition cannot be verified");
        if (!cross_nonzero) {
            return ExpressionResult::failure(cross_nonzero.error());
        }
        auto cross_norm_sq = vector_dot_checked(
            cross.value(), cross.value(), context);
        if (!cross_norm_sq) {
            return ExpressionResult::failure(cross_norm_sq.error());
        }
        auto cross_norm = SymbolicExpr::sqrt(cross_norm_sq.value());
        auto numerator = vector_dot_checked(
            a2_minus_a1, cross.value(), context);
        if (!numerator) return ExpressionResult::failure(numerator.error());
        auto abs_num = LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Abs,
                std::vector<std::shared_ptr<const SymbolicNode>>{
                    LMCAS::detail::node(numerator.value())}));
        auto distance = SymbolicExpr::divide(abs_num, cross_norm);
        if (scale_point_difference) {
            distance = SymbolicExpr::multiply(
                distance, SymbolicExpr::number(point_scale));
        }
        return simplify_checked(
            distance, operation,
            "skew line distance is outside the supported domain");
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                         "vector geometry allocation failed",
                                         operation);
    } catch (const std::exception& e) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                         e.what(), operation);
    }
}

ExpressionResult skew_lines_distance_checked(
    const LineSymbolic& l1,
    const LineSymbolic& l2
) {
    ComputationContext context;
    return skew_lines_distance_checked(l1, l2, context);
}


LineSymbolicResult line_from_two_points_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2,
    ComputationContext& context) {
    const std::string operation = "line_from_two_points";
    try {
        auto valid = validate_same_dimension_vectors(p1, p2, context, operation);
        if (!valid) return LineSymbolicResult::failure(valid.error());
        if (p1.size() != 3) {
            return LineSymbolicResult::failure(
                CasErrc::InvalidArgument,
                "line construction requires 3-dimensional points",
                operation);
        }
        auto step = context.consume_steps(4, operation);
        if (!step) return LineSymbolicResult::failure(step.error());

        std::vector<double> p1_values;
        std::vector<double> p2_values;
        auto p1_scale = numeric_vector_scale_checked(
            p1, context, operation, &p1_values);
        auto p2_scale = numeric_vector_scale_checked(
            p2, context, operation, &p2_values);
        const bool numeric_points = p1_scale.has_value() && p2_scale.has_value();
        if (!p1_scale &&
            (p1_scale.error().code == CasErrc::Cancelled ||
             p1_scale.error().code == CasErrc::ResourceLimit)) {
            return LineSymbolicResult::failure(p1_scale.error());
        }
        if (!p2_scale &&
            (p2_scale.error().code == CasErrc::Cancelled ||
             p2_scale.error().code == CasErrc::ResourceLimit)) {
            return LineSymbolicResult::failure(p2_scale.error());
        }

        bool scale_difference = false;
        double point_scale = 0.0;
        if (numeric_points) {
            point_scale = std::max(p1_scale.value(), p2_scale.value());
            for (size_t i = 0; i < 3; ++i) {
                if (!std::isfinite(p2_values[i] - p1_values[i])) {
                    scale_difference = true;
                    break;
                }
            }
        }

        LineSymbolic line;
        line.point = p1;
        line.direction.reserve(3);
        for (size_t i = 0; i < p1.size(); ++i) {
            if (scale_difference) {
                line.direction.push_back(SymbolicExpr::number(
                    p2_values[i] / point_scale -
                    p1_values[i] / point_scale));
                continue;
            }
            auto component = SymbolicExpr::add(
                p2[i],
                SymbolicExpr::multiply(SymbolicExpr::number(-1), p1[i]));
            auto simplified = simplify_checked(
                component, operation,
                "line direction component is outside the supported domain");
            if (!simplified) return LineSymbolicResult::failure(simplified.error());
            line.direction.push_back(std::move(simplified.value()));
        }
        auto line_valid = validate_line_checked(line, context, operation);
        if (!line_valid) return LineSymbolicResult::failure(line_valid.error());
        return LineSymbolicResult::success(std::move(line));
    } catch (const std::bad_alloc&) {
        return LineSymbolicResult::failure(CasErrc::ResourceLimit,
                                           "vector geometry allocation failed",
                                           operation);
    } catch (const std::exception& e) {
        return LineSymbolicResult::failure(CasErrc::InternalInvariant,
                                           e.what(), operation);
    }
}

LineSymbolicResult line_from_two_points_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2) {
    ComputationContext context;
    return line_from_two_points_checked(p1, p2, context);
}


PlaneSymbolicResult plane_from_three_points_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p3,
    ComputationContext& context) {
    const std::string operation = "plane_from_three_points";
    try {
        auto v12 = validate_same_dimension_vectors(p1, p2, context, operation);
        if (!v12) return PlaneSymbolicResult::failure(v12.error());
        auto v13 = validate_same_dimension_vectors(p1, p3, context, operation);
        if (!v13) return PlaneSymbolicResult::failure(v13.error());
        if (p1.size() != 3) {
            return PlaneSymbolicResult::failure(
                CasErrc::InvalidArgument,
                "plane construction requires 3-dimensional points",
                operation);
        }
        auto step = context.consume_steps(8, operation);
        if (!step) return PlaneSymbolicResult::failure(step.error());

        std::vector<double> p1_values;
        std::vector<double> p2_values;
        std::vector<double> p3_values;
        auto p1_scale = numeric_vector_scale_checked(
            p1, context, operation, &p1_values);
        auto p2_scale = numeric_vector_scale_checked(
            p2, context, operation, &p2_values);
        auto p3_scale = numeric_vector_scale_checked(
            p3, context, operation, &p3_values);
        if (!p1_scale &&
            (p1_scale.error().code == CasErrc::Cancelled ||
             p1_scale.error().code == CasErrc::ResourceLimit)) {
            return PlaneSymbolicResult::failure(p1_scale.error());
        }
        if (!p2_scale &&
            (p2_scale.error().code == CasErrc::Cancelled ||
             p2_scale.error().code == CasErrc::ResourceLimit)) {
            return PlaneSymbolicResult::failure(p2_scale.error());
        }
        if (!p3_scale &&
            (p3_scale.error().code == CasErrc::Cancelled ||
             p3_scale.error().code == CasErrc::ResourceLimit)) {
            return PlaneSymbolicResult::failure(p3_scale.error());
        }

        const bool numeric_points =
            p1_scale.has_value() && p2_scale.has_value() && p3_scale.has_value();
        bool scale_differences = false;
        double point_scale = 0.0;
        if (numeric_points) {
            point_scale = std::max(
                p1_scale.value(), std::max(p2_scale.value(), p3_scale.value()));
            for (size_t i = 0; i < 3; ++i) {
                if (!std::isfinite(p2_values[i] - p1_values[i]) ||
                    !std::isfinite(p3_values[i] - p1_values[i])) {
                    scale_differences = true;
                    break;
                }
            }
        }

        std::vector<std::shared_ptr<SymbolicExpr>> v1, v2;
        v1.reserve(3);
        v2.reserve(3);
        for (size_t i = 0; i < 3; ++i) {
            if (scale_differences) {
                v1.push_back(SymbolicExpr::number(
                    p2_values[i] / point_scale -
                    p1_values[i] / point_scale));
                v2.push_back(SymbolicExpr::number(
                    p3_values[i] / point_scale -
                    p1_values[i] / point_scale));
                continue;
            }
            auto left = simplify_checked(
                SymbolicExpr::add(
                    p2[i],
                    SymbolicExpr::multiply(SymbolicExpr::number(-1), p1[i])),
                operation,
                "plane edge vector is outside the supported domain");
            if (!left) return PlaneSymbolicResult::failure(left.error());
            auto right = simplify_checked(
                SymbolicExpr::add(
                    p3[i],
                    SymbolicExpr::multiply(SymbolicExpr::number(-1), p1[i])),
                operation,
                "plane edge vector is outside the supported domain");
            if (!right) return PlaneSymbolicResult::failure(right.error());
            v1.push_back(std::move(left.value()));
            v2.push_back(std::move(right.value()));
        }
        std::vector<double> v1_values;
        std::vector<double> v2_values;
        auto v1_scale = numeric_vector_scale_checked(
            v1, context, operation, &v1_values);
        auto v2_scale = numeric_vector_scale_checked(
            v2, context, operation, &v2_values);
        const bool numeric_edges = v1_scale.has_value() && v2_scale.has_value();
        if (!numeric_edges) {
            if ((!v1_scale &&
                 (v1_scale.error().code == CasErrc::Cancelled ||
                  v1_scale.error().code == CasErrc::ResourceLimit)) ||
                (!v2_scale &&
                 (v2_scale.error().code == CasErrc::Cancelled ||
                  v2_scale.error().code == CasErrc::ResourceLimit))) {
                return PlaneSymbolicResult::failure(
                    !v1_scale ? v1_scale.error() : v2_scale.error());
            }
        } else if (v1_scale.value() != 0.0 && v2_scale.value() != 0.0) {
            const double safe_high =
                std::sqrt(std::numeric_limits<double>::max() / 12.0);
            const double safe_low =
                std::sqrt(std::numeric_limits<double>::min());
            const bool scale_cross =
                v1_scale.value() > safe_high / v2_scale.value() ||
                v1_scale.value() < safe_low / v2_scale.value();
            if (scale_cross) {
                for (size_t i = 0; i < 3; ++i) {
                    v1[i] = SymbolicExpr::number(
                        v1_values[i] / v1_scale.value());
                    v2[i] = SymbolicExpr::number(
                        v2_values[i] / v2_scale.value());
                }
            }
        }

        auto normal = vector_cross_checked(v1, v2, context);
        if (!normal) return PlaneSymbolicResult::failure(normal.error());
        Result<void> nonzero = numeric_edges
            ? checked_nonzero_numeric_vector(
                  normal.value(), context, operation,
                  "three points are collinear, so no unique plane exists",
                  "plane normal nonzero condition cannot be verified")
            : Result<void>::failure(
                  CasErrc::Inconclusive,
                  "plane normal nonzero condition cannot be verified",
                  operation);
        if (!numeric_edges) {
            auto norm_sq = vector_dot_checked(
                normal.value(), normal.value(), context);
            if (!norm_sq) return PlaneSymbolicResult::failure(norm_sq.error());
            nonzero = checked_nonzero_numeric_or_exact(
                norm_sq.value(), context, operation,
                "three points are collinear, so no unique plane exists",
                "plane normal nonzero condition cannot be verified");
        }
        if (!nonzero) return PlaneSymbolicResult::failure(nonzero.error());
        auto d = vector_dot_checked(normal.value(), p1, context);
        if (!d) return PlaneSymbolicResult::failure(d.error());
        auto d_simplified = simplify_checked(
            d.value(), operation,
            "plane constant is outside the supported domain");
        if (!d_simplified) return PlaneSymbolicResult::failure(d_simplified.error());

        PlaneSymbolic plane;
        plane.normal = std::move(normal.value());
        for (auto& component : plane.normal) {
            auto simplified = simplify_checked(
                component, operation,
                "plane normal component is outside the supported domain");
            if (!simplified) return PlaneSymbolicResult::failure(simplified.error());
            component = std::move(simplified.value());
        }
        plane.d = std::move(d_simplified.value());
        return PlaneSymbolicResult::success(std::move(plane));
    } catch (const std::bad_alloc&) {
        return PlaneSymbolicResult::failure(CasErrc::ResourceLimit,
                                            "vector geometry allocation failed",
                                            operation);
    } catch (const std::exception& e) {
        return PlaneSymbolicResult::failure(CasErrc::InternalInvariant,
                                            e.what(), operation);
    }
}

PlaneSymbolicResult plane_from_three_points_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p3) {
    ComputationContext context;
    return plane_from_three_points_checked(p1, p2, p3, context);
}


ExpressionResult dihedral_angle_checked(
    const PlaneSymbolic& p1,
    const PlaneSymbolic& p2,
    ComputationContext& context) {
    const std::string operation = "dihedral_angle";
    try {
        auto p1_valid = validate_plane_checked(p1, context, operation);
        if (!p1_valid) return ExpressionResult::failure(p1_valid.error());
        auto p2_valid = validate_plane_checked(p2, context, operation);
        if (!p2_valid) return ExpressionResult::failure(p2_valid.error());

        auto angle = vector_angle_checked(p1.normal, p2.normal, context);
        if (!angle) return ExpressionResult::failure(angle.error());
        const double pi = std::acos(-1.0);
        const double acute_angle =
            std::min(angle.value(), pi - angle.value());
        return ExpressionResult::success(
            SymbolicExpr::number(acute_angle));
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                         "vector geometry allocation failed",
                                         operation);
    } catch (const std::exception& e) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                         e.what(), operation);
    }
}

ExpressionResult dihedral_angle_checked(
    const PlaneSymbolic& p1,
    const PlaneSymbolic& p2) {
    ComputationContext context;
    return dihedral_angle_checked(p1, p2, context);
}


VectorStringResult classify_quadric_checked(
    const SurfaceSymbolic& surf,
    ComputationContext& context) {
    const std::string operation = "classify_quadric";
    try {
        auto valid = validate_surface_point(
            surf,
            std::vector<std::shared_ptr<SymbolicExpr>>(
                surf.vars.size(), SymbolicExpr::number(0)),
            context,
            operation);
        if (!valid) return VectorStringResult::failure(valid.error());
        if (surf.vars.size() != 3) {
            return VectorStringResult::failure(
                CasErrc::InvalidArgument,
                "quadric classification requires exactly three variables",
                operation);
        }
        auto step = context.consume_steps(12, operation);
        if (!step) return VectorStringResult::failure(step.error());
        return detail::classify_quadric_impl(surf, context, operation);
    } catch (const std::bad_alloc&) {
        return VectorStringResult::failure(CasErrc::ResourceLimit,
                                           "vector geometry allocation failed",
                                           operation);
    } catch (const std::exception& e) {
        return VectorStringResult::failure(CasErrc::Inconclusive,
                                           e.what(), operation);
    }
}

VectorStringResult classify_quadric_checked(const SurfaceSymbolic& surf) {
    ComputationContext context;
    return classify_quadric_checked(surf, context);
}


VectorExprListResult surface_normal_checked(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    ComputationContext& context) {
    const std::string operation = "surface_normal";
    auto valid = validate_surface_point(surf, point, context, operation);
    if (!valid) return VectorExprListResult::failure(valid.error());
    auto step = context.consume_steps(surf.vars.size() * 4 + 4, operation);
    if (!step) return VectorExprListResult::failure(step.error());

    auto gradient = surface_gradient_at_point_checked(surf, point, operation);
    if (!gradient) return gradient;
    std::vector<double> gradient_values;
    auto scaled = numeric_vector_scale_checked(
        gradient.value(), context, operation, &gradient_values);
    if (!scaled) {
        if (scaled.error().code == CasErrc::Cancelled ||
            scaled.error().code == CasErrc::ResourceLimit) {
            return VectorExprListResult::failure(scaled.error());
        }
        return VectorExprListResult::failure(
            CasErrc::Inconclusive,
            "surface gradient nonzero condition cannot be verified",
            operation);
    }
    if (scaled.value() == 0.0) {
        return VectorExprListResult::failure(
            CasErrc::DomainError,
            "surface normal is undefined at a singular point",
            operation);
    }

    double scaled_norm_sq = 0.0;
    for (double component : gradient_values) {
        const double ratio = component / scaled.value();
        scaled_norm_sq += ratio * ratio;
    }
    const double scaled_norm = std::sqrt(scaled_norm_sq);

    std::vector<std::shared_ptr<SymbolicExpr>> result;
    result.reserve(gradient_values.size());
    for (double component : gradient_values) {
        result.push_back(SymbolicExpr::number(
            (component / scaled.value()) / scaled_norm));
    }
    return VectorExprListResult::success(std::move(result));
}

VectorExprListResult surface_normal_checked(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point) {
    ComputationContext context;
    return surface_normal_checked(surf, point, context);
}


PlaneSymbolicResult tangent_plane_checked(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    ComputationContext& context) {
    const std::string operation = "tangent_plane";
    auto valid = validate_surface_point(surf, point, context, operation);
    if (!valid) return PlaneSymbolicResult::failure(valid.error());
    auto step = context.consume_steps(surf.vars.size() * 4 + 4, operation);
    if (!step) return PlaneSymbolicResult::failure(step.error());

    auto gradient = surface_gradient_at_point_checked(surf, point, operation);
    if (!gradient) return PlaneSymbolicResult::failure(gradient.error());
    std::vector<double> gradient_values;
    auto scaled = numeric_vector_scale_checked(
        gradient.value(), context, operation, &gradient_values);
    if (!scaled) {
        if (scaled.error().code == CasErrc::Cancelled ||
            scaled.error().code == CasErrc::ResourceLimit) {
            return PlaneSymbolicResult::failure(scaled.error());
        }
        return PlaneSymbolicResult::failure(
            CasErrc::Inconclusive,
            "surface gradient nonzero condition cannot be verified",
            operation);
    }
    if (scaled.value() == 0.0) {
        return PlaneSymbolicResult::failure(
            CasErrc::DomainError,
            "tangent plane is undefined at a singular point",
            operation);
    }

    std::vector<double> point_values;
    auto point_scale = numeric_vector_scale_checked(
        point, context, operation, &point_values);
    if (!point_scale &&
        (point_scale.error().code == CasErrc::Cancelled ||
         point_scale.error().code == CasErrc::ResourceLimit)) {
        return PlaneSymbolicResult::failure(point_scale.error());
    }
    if (point_scale && point_scale.value() != 0.0 &&
        scaled.value() >
            std::numeric_limits<double>::max() /
                point_scale.value() /
                static_cast<double>(point.size())) {
        double normalized_dot = 0.0;
        for (size_t i = 0; i < point.size(); ++i) {
            normalized_dot = std::fma(
                gradient_values[i] / scaled.value(),
                point_values[i] / point_scale.value(),
                normalized_dot);
        }
        const double plane_constant = normalized_dot * point_scale.value();
        if (!std::isfinite(plane_constant)) {
            return PlaneSymbolicResult::failure(
                CasErrc::NumericFailure,
                "tangent plane constant is outside the supported domain",
                operation);
        }

        PlaneSymbolic plane;
        plane.normal.reserve(gradient_values.size());
        for (double component : gradient_values) {
            plane.normal.push_back(
                SymbolicExpr::number(component / scaled.value()));
        }
        plane.d = SymbolicExpr::number(plane_constant);
        return PlaneSymbolicResult::success(std::move(plane));
    }

    auto dot = vector_dot_checked(gradient.value(), point, context);
    if (!dot) return PlaneSymbolicResult::failure(dot.error());
    auto d = simplify_checked(
        dot.value(), operation,
        "tangent plane constant is outside the supported domain");
    if (!d) return PlaneSymbolicResult::failure(d.error());

    PlaneSymbolic plane;
    plane.normal = std::move(gradient.value());
    plane.d = std::move(d.value());
    return PlaneSymbolicResult::success(std::move(plane));
}

PlaneSymbolicResult tangent_plane_checked(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point) {
    ComputationContext context;
    return tangent_plane_checked(surf, point, context);
}

}
