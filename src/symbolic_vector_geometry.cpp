#include "../include/symbolic_vector_geometry.hpp"
#include "symbolic_ast.hpp"
#include "../include/numeric_evaluation.hpp"
#include "../include/symbolic.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace lamina {

namespace {

Result<void> validate_symbolic_vector(
    const std::vector<std::shared_ptr<SymbolicExpr>>& vector,
    const std::string& name,
    const std::string& operation)
{
    for (size_t i = 0; i < vector.size(); ++i) {
        if (!vector[i] || !lamina::detail::node(vector[i])) {
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
    if (!surf.F || !lamina::detail::node(surf.F)) {
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
    if (!expr || !lamina::detail::node(expr)) {
        return ExpressionResult::failure(CasErrc::Inconclusive, message, operation);
    }
    auto simplified = expr->simplify();
    if (!simplified || !lamina::detail::node(simplified)) {
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
            if (!substituted || !lamina::detail::node(substituted)) {
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

ExpressionResult gradient_norm_sq_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& gradient,
    const std::string& operation)
{
    std::shared_ptr<SymbolicExpr> norm_sq = SymbolicExpr::number(0);
    for (const auto& component : gradient) {
        auto square = SymbolicExpr::multiply(component, component);
        norm_sq = SymbolicExpr::add(norm_sq, square);
        auto simplified = simplify_checked(
            norm_sq, operation,
            "surface gradient norm construction is outside the supported domain");
        if (!simplified) return simplified;
        norm_sq = std::move(simplified.value());
    }
    return ExpressionResult::success(std::move(norm_sq));
}

Result<void> checked_nonzero_numeric_or_exact(
    const std::shared_ptr<SymbolicExpr>& expr,
    ComputationContext& context,
    const std::string& operation,
    const std::string& domain_message,
    const std::string& inconclusive_message)
{
    if (!expr || !lamina::detail::node(expr)) {
        return Result<void>::failure(CasErrc::Inconclusive,
                                     inconclusive_message, operation);
    }
    auto simplified = expr->simplify();
    if (!simplified || !lamina::detail::node(simplified)) {
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
    auto norm = vector_dot_checked(line.direction, line.direction, context);
    if (!norm) return Result<void>::failure(norm.error());
    return checked_nonzero_numeric_or_exact(
        norm.value(), context, operation,
        "line direction cannot be zero",
        "line direction nonzero condition cannot be verified");
}

Result<void> validate_plane_checked(const PlaneSymbolic& plane,
                                    ComputationContext& context,
                                    const std::string& operation)
{
    if (!plane.d || !lamina::detail::node(plane.d)) {
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
    auto norm = vector_dot_checked(plane.normal, plane.normal, context);
    if (!norm) return Result<void>::failure(norm.error());
    return checked_nonzero_numeric_or_exact(
        norm.value(), context, operation,
        "plane normal cannot be zero",
        "plane normal nonzero condition cannot be verified");
}

} // namespace

Result<double> symbolic_vector_finite_numeric(
    const std::shared_ptr<SymbolicExpr>& expr,
    ComputationContext& context,
    const std::string& operation)
{
    if (!expr || !lamina::detail::node(expr)) {
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
        sum_terms.push_back(lamina::detail::node(SymbolicExpr::multiply(a[i], b[i])));
    }
    return ExpressionResult::success(
        lamina::detail::make_expression_ptr(lamina::detail::make_node<AddNode>(sum_terms)));
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
    double norm_a = 0, norm_b = 0, dot = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        auto step = context.consume_steps(1, operation);
        if (!step) return VectorAngleResult::failure(step.error());
        auto na_result = symbolic_vector_finite_numeric(a[i], context, operation);
        if (!na_result) return VectorAngleResult::failure(na_result.error());
        auto nb_result = symbolic_vector_finite_numeric(b[i], context, operation);
        if (!nb_result) return VectorAngleResult::failure(nb_result.error());
        double na = na_result.value();
        double nb = nb_result.value();
        norm_a += na * na;
        norm_b += nb * nb;
        dot += na * nb;
    }
    if (norm_a == 0.0 || norm_b == 0.0) {
        return VectorAngleResult::failure(
            CasErrc::DomainError,
            "angle is undefined for zero-length vectors",
            operation);
    }
    double cosv = dot / (std::sqrt(norm_a) * std::sqrt(norm_b));
    // Clamp into [-1, 1] to avoid std::acos domain errors caused by rounding.
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

        auto n_dot_a = vector_dot_checked(plane.normal, line.point, context);
        if (!n_dot_a) return VectorExprListResult::failure(n_dot_a.error());
        auto n_dot_b = vector_dot_checked(plane.normal, line.direction, context);
        if (!n_dot_b) return VectorExprListResult::failure(n_dot_b.error());
        auto denom = n_dot_b.value()->simplify();
        auto nonzero = checked_nonzero_numeric_or_exact(
            denom, context, operation,
            "line is parallel to plane, so no unique intersection exists",
            "line-plane denominator nonzero condition cannot be verified");
        if (!nonzero) return VectorExprListResult::failure(nonzero.error());

        auto t = SymbolicExpr::divide(
            SymbolicExpr::add(
                plane.d,
                SymbolicExpr::multiply(SymbolicExpr::number(-1), n_dot_a.value())),
            denom);
        std::vector<std::shared_ptr<SymbolicExpr>> intersection;
        intersection.reserve(line.point.size());
        for (size_t i = 0; i < line.point.size(); ++i) {
            auto coord = SymbolicExpr::add(
                line.point[i],
                SymbolicExpr::multiply(t, line.direction[i]));
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

        auto n_dot_r = vector_dot_checked(plane.normal, point, context);
        if (!n_dot_r) return ExpressionResult::failure(n_dot_r.error());
        auto norm_sq = vector_dot_checked(plane.normal, plane.normal, context);
        if (!norm_sq) return ExpressionResult::failure(norm_sq.error());
        auto norm_nonzero = checked_nonzero_numeric_or_exact(
            norm_sq.value(), context, operation,
            "plane normal cannot be zero",
            "plane normal nonzero condition cannot be verified");
        if (!norm_nonzero) return ExpressionResult::failure(norm_nonzero.error());

        auto diff = SymbolicExpr::add(
            n_dot_r.value(),
            SymbolicExpr::multiply(SymbolicExpr::number(-1), plane.d));
        auto norm_n = SymbolicExpr::sqrt(norm_sq.value());
        auto abs_diff = lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Abs,
                std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(diff)}));
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

        std::vector<std::shared_ptr<SymbolicExpr>> a2_minus_a1;
        a2_minus_a1.reserve(3);
        for (size_t i = 0; i < l1.point.size(); ++i) {
            a2_minus_a1.push_back(SymbolicExpr::add(
                l2.point[i],
                SymbolicExpr::multiply(SymbolicExpr::number(-1), l1.point[i])));
        }
        auto cross = vector_cross_checked(l1.direction, l2.direction, context);
        if (!cross) return ExpressionResult::failure(cross.error());
        auto cross_norm_sq = vector_dot_checked(cross.value(), cross.value(), context);
        if (!cross_norm_sq) return ExpressionResult::failure(cross_norm_sq.error());
        auto cross_nonzero = checked_nonzero_numeric_or_exact(
            cross_norm_sq.value(), context, operation,
            "skew line distance requires non-parallel directions",
            "cross-product nonzero condition cannot be verified");
        if (!cross_nonzero) return ExpressionResult::failure(cross_nonzero.error());
        auto cross_norm = SymbolicExpr::sqrt(cross_norm_sq.value());
        auto numerator = vector_dot_checked(a2_minus_a1, cross.value(), context);
        if (!numerator) return ExpressionResult::failure(numerator.error());
        auto abs_num = lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Abs,
                std::vector<std::shared_ptr<const SymbolicNode>>{
                    lamina::detail::node(numerator.value())}));
        return simplify_checked(
            SymbolicExpr::divide(abs_num, cross_norm),
            operation,
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

        LineSymbolic line;
        line.point = p1;
        line.direction.reserve(3);
        for (size_t i = 0; i < p1.size(); ++i) {
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

        std::vector<std::shared_ptr<SymbolicExpr>> v1, v2;
        v1.reserve(3);
        v2.reserve(3);
        for (size_t i = 0; i < 3; ++i) {
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
        auto normal = vector_cross_checked(v1, v2, context);
        if (!normal) return PlaneSymbolicResult::failure(normal.error());
        auto norm_sq = vector_dot_checked(normal.value(), normal.value(), context);
        if (!norm_sq) return PlaneSymbolicResult::failure(norm_sq.error());
        auto nonzero = checked_nonzero_numeric_or_exact(
            norm_sq.value(), context, operation,
            "three points are collinear, so no unique plane exists",
            "plane normal nonzero condition cannot be verified");
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
        auto step = context.consume_steps(6, operation);
        if (!step) return ExpressionResult::failure(step.error());

        auto dot = vector_dot_checked(p1.normal, p2.normal, context);
        if (!dot) return ExpressionResult::failure(dot.error());
        auto n1_sq = vector_dot_checked(p1.normal, p1.normal, context);
        if (!n1_sq) return ExpressionResult::failure(n1_sq.error());
        auto n2_sq = vector_dot_checked(p2.normal, p2.normal, context);
        if (!n2_sq) return ExpressionResult::failure(n2_sq.error());
        auto n1_nonzero = checked_nonzero_numeric_or_exact(
            n1_sq.value(), context, operation,
            "first plane normal cannot be zero",
            "first plane normal nonzero condition cannot be verified");
        if (!n1_nonzero) return ExpressionResult::failure(n1_nonzero.error());
        auto n2_nonzero = checked_nonzero_numeric_or_exact(
            n2_sq.value(), context, operation,
            "second plane normal cannot be zero",
            "second plane normal nonzero condition cannot be verified");
        if (!n2_nonzero) return ExpressionResult::failure(n2_nonzero.error());

        auto n1 = SymbolicExpr::sqrt(n1_sq.value());
        auto n2 = SymbolicExpr::sqrt(n2_sq.value());
        auto abs_dot = lamina::detail::make_expression_ptr(
            lamina::detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Abs,
                std::vector<std::shared_ptr<const SymbolicNode>>{
                    lamina::detail::node(dot.value())}));
        auto cos_theta = SymbolicExpr::divide(abs_dot, SymbolicExpr::multiply(n1, n2));
        auto arccos = lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::ArcCos,
            std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(cos_theta)});
        return simplify_checked(
            lamina::detail::make_expression_ptr(arccos),
            operation,
            "dihedral angle is outside the supported domain");
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

        const auto& vars = surf.vars;
        auto F = surf.F->expand();
        if (!F || !lamina::detail::node(F)) F = surf.F;

        auto numeric_coeff = [&](const std::shared_ptr<SymbolicExpr>& expr,
                                 const std::string& message)
            -> Result<double> {
            auto simplified = expr ? expr->simplify() : nullptr;
            if (!simplified || !lamina::detail::node(simplified)) {
                return Result<double>::failure(CasErrc::Inconclusive,
                                               message, operation);
            }
            auto value = symbolic_vector_finite_numeric(
                simplified, context, operation);
            if (!value) {
                if (value.error().code == CasErrc::Cancelled ||
                    value.error().code == CasErrc::ResourceLimit) {
                    return Result<double>::failure(value.error());
                }
                return Result<double>::failure(CasErrc::Inconclusive,
                                               message, operation);
            }
            return Result<double>::success(value.value());
        };

        auto coeff_of_sq = [&](const std::string& v) -> Result<double> {
            auto d2 = F->differentiate(v);
            if (!d2 || !lamina::detail::node(d2)) {
                return Result<double>::failure(
                    CasErrc::Inconclusive,
                    "quadric second derivative is outside the supported domain",
                    operation);
            }
            d2 = d2->differentiate(v);
            auto value = numeric_coeff(
                d2,
                "quadric second derivative coefficient cannot be evaluated");
            if (!value) return value;
            return Result<double>::success(value.value() / 2.0);
        };
        auto coeff_of_lin = [&](const std::string& v) -> Result<double> {
            auto d = F->differentiate(v);
            if (!d || !lamina::detail::node(d)) {
                return Result<double>::failure(
                    CasErrc::Inconclusive,
                    "quadric first derivative is outside the supported domain",
                    operation);
            }
            for (const auto& w : vars) {
                d = d->substitute(w, SymbolicExpr::number(0));
                if (!d || !lamina::detail::node(d)) {
                    return Result<double>::failure(
                        CasErrc::Inconclusive,
                        "quadric linear coefficient cannot be substituted",
                        operation);
                }
            }
            return numeric_coeff(
                d,
                "quadric linear coefficient cannot be evaluated");
        };

        auto a_result = coeff_of_sq(vars[0]);
        if (!a_result) return VectorStringResult::failure(a_result.error());
        auto b_result = coeff_of_sq(vars[1]);
        if (!b_result) return VectorStringResult::failure(b_result.error());
        auto c_result = coeff_of_sq(vars[2]);
        if (!c_result) return VectorStringResult::failure(c_result.error());
        auto lz_result = coeff_of_lin(vars[2]);
        if (!lz_result) return VectorStringResult::failure(lz_result.error());

        double a = a_result.value();
        double b = b_result.value();
        double c = c_result.value();
        double lz = lz_result.value();

        int n_sq = (a != 0) + (b != 0) + (c != 0);
        int n_pos = (a > 0) + (b > 0) + (c > 0);
        int n_neg = (a < 0) + (b < 0) + (c < 0);

        if (n_sq == 2 && lz != 0 && c == 0) {
            return VectorStringResult::success("paraboloid");
        }
        if (n_sq == 2) return VectorStringResult::success("cylinder");
        if (n_sq == 3) {
            if (n_neg == 0) {
                if (a == b && b == c) {
                    return VectorStringResult::success("sphere");
                }
                return VectorStringResult::success("ellipsoid");
            }
            if (n_pos == 0) return VectorStringResult::success("ellipsoid");

            auto constant = F;
            for (const auto& w : vars) {
                constant = constant->substitute(w, SymbolicExpr::number(0));
                if (!constant || !lamina::detail::node(constant)) {
                    return VectorStringResult::failure(
                        CasErrc::Inconclusive,
                        "quadric constant term cannot be substituted",
                        operation);
                }
            }
            auto cst_value = numeric_coeff(
                constant,
                "quadric constant term cannot be evaluated");
            if (!cst_value) return VectorStringResult::failure(cst_value.error());
            if (std::abs(cst_value.value()) < 1e-9) {
                return VectorStringResult::success("cone");
            }
            return VectorStringResult::success("hyperboloid");
        }
        return VectorStringResult::failure(
            CasErrc::Inconclusive,
            "quadric is outside the currently classified support domain",
            operation);
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

    auto norm_sq = gradient_norm_sq_checked(gradient.value(), operation);
    if (!norm_sq) return VectorExprListResult::failure(norm_sq.error());
    auto norm_sq_simplified = norm_sq.value()->simplify();
    if (!norm_sq_simplified || !lamina::detail::node(norm_sq_simplified)) {
        return VectorExprListResult::failure(
            CasErrc::Inconclusive,
            "surface gradient norm cannot be verified",
            operation);
    }
    if (norm_sq_simplified->is_zero()) {
        return VectorExprListResult::failure(
            CasErrc::DomainError,
            "surface normal is undefined at a singular point",
            operation);
    }
    auto numeric_norm_sq = symbolic_vector_finite_numeric(
        norm_sq_simplified, context, operation);
    if (!numeric_norm_sq) {
        if (numeric_norm_sq.error().code == CasErrc::Cancelled ||
            numeric_norm_sq.error().code == CasErrc::ResourceLimit) {
            return VectorExprListResult::failure(numeric_norm_sq.error());
        }
        return VectorExprListResult::failure(
            CasErrc::Inconclusive,
            "surface gradient nonzero condition cannot be verified",
            operation);
    }
    if (numeric_norm_sq.value() <= 0.0) {
        return VectorExprListResult::failure(
            CasErrc::DomainError,
            "surface normal is undefined at a singular point",
            operation);
    }

    auto norm = SymbolicExpr::sqrt(norm_sq_simplified);
    auto norm_checked = simplify_checked(
        norm, operation, "surface gradient norm construction is outside the supported domain");
    if (!norm_checked) return VectorExprListResult::failure(norm_checked.error());

    std::vector<std::shared_ptr<SymbolicExpr>> result;
    result.reserve(gradient.value().size());
    for (const auto& component : gradient.value()) {
        auto normalized = SymbolicExpr::divide(component, norm_checked.value());
        auto checked_component = simplify_checked(
            normalized, operation,
            "surface normal component is outside the supported domain");
        if (!checked_component) {
            return VectorExprListResult::failure(checked_component.error());
        }
        result.push_back(std::move(checked_component.value()));
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

    auto norm_sq = gradient_norm_sq_checked(gradient.value(), operation);
    if (!norm_sq) return PlaneSymbolicResult::failure(norm_sq.error());
    auto norm_sq_simplified = norm_sq.value()->simplify();
    if (!norm_sq_simplified || !lamina::detail::node(norm_sq_simplified)) {
        return PlaneSymbolicResult::failure(
            CasErrc::Inconclusive,
            "surface gradient norm cannot be verified",
            operation);
    }
    if (norm_sq_simplified->is_zero()) {
        return PlaneSymbolicResult::failure(
            CasErrc::DomainError,
            "tangent plane is undefined at a singular point",
            operation);
    }
    auto numeric_norm_sq = symbolic_vector_finite_numeric(
        norm_sq_simplified, context, operation);
    if (!numeric_norm_sq) {
        if (numeric_norm_sq.error().code == CasErrc::Cancelled ||
            numeric_norm_sq.error().code == CasErrc::ResourceLimit) {
            return PlaneSymbolicResult::failure(numeric_norm_sq.error());
        }
        return PlaneSymbolicResult::failure(
            CasErrc::Inconclusive,
            "surface gradient nonzero condition cannot be verified",
            operation);
    }
    if (numeric_norm_sq.value() <= 0.0) {
        return PlaneSymbolicResult::failure(
            CasErrc::DomainError,
            "tangent plane is undefined at a singular point",
            operation);
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
