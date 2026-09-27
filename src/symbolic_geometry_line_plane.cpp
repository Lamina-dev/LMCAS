#include "symbolic_geometry_line_plane.hpp"
#include "internal/symbolic_geometry_support.hpp"
#include "internal/symbolic_ast.hpp"
#include "numeric_evaluation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace LMCAS {
using namespace geometry_detail;

namespace {

Result<void> validate_intersection_inputs(
    const LineSymbolic& line, const PlaneSymbolic& plane,
    ComputationContext& context, const std::string& operation
) {
    auto line_valid = validate_line_checked(line, context, operation);
    if (!line_valid) return line_valid;
    auto plane_valid = validate_plane_checked(plane, context, operation);
    if (!plane_valid) return plane_valid;
    return context.consume_steps(6, operation);
}

bool interrupts_optional_numeric_evaluation(const CasError& error) {
    return error.code == CasErrc::Cancelled ||
           error.code == CasErrc::ResourceLimit;
}

bool finite_plane_offset(const Result<ApproxReal>& offset) {
    return offset.has_value() && offset.value().is_finite() &&
           std::isfinite(offset.value().value);
}

bool direction_products_need_scaling(double normal_scale, double direction_scale) {
    const double safe_high = std::numeric_limits<double>::max() / 3.0;
    const double safe_low = std::numeric_limits<double>::min();
    return normal_scale > safe_high / direction_scale ||
           normal_scale < safe_low / direction_scale;
}

bool point_dot_needs_protection(
    double normal_scale, const Result<double>& point_scale
) {
    return point_scale.has_value() && point_scale.value() > 0.0 &&
           normal_scale >
               std::numeric_limits<double>::max() / point_scale.value() / 3.0;
}

bool coordinate_update_needs_protection(
    const Result<double>& point_scale, bool numeric_plane_offset
) {
    return point_scale.has_value() && numeric_plane_offset &&
           point_scale.value() > std::numeric_limits<double>::max() / 3.0;
}

void scale_symbolic_components(
    const std::vector<double>& values, double scale, bool protect_point_dot,
    std::vector<std::shared_ptr<SymbolicExpr>>& scaled
) {
    scaled.reserve(values.size());
    for (double component : values) {
        double scaled_component = component / scale;
        if (protect_point_dot) scaled_component /= 3.0;
        scaled.push_back(SymbolicExpr::number(scaled_component));
    }
}

std::shared_ptr<SymbolicExpr> scale_plane_offset(
    const std::shared_ptr<SymbolicExpr>& offset,
    double normal_scale, bool protect_point_dot
) {
    auto scaled = SymbolicExpr::divide(
        offset, SymbolicExpr::number(normal_scale));
    if (protect_point_dot) {
        scaled = SymbolicExpr::divide(scaled, SymbolicExpr::number(3.0));
    }
    return scaled;
}

double fused_dot_three(const double (&a)[3], const double (&b)[3]) {
    return std::fma(a[0], b[0], std::fma(a[1], b[1], a[2] * b[2]));
}

VectorExprListResult protected_coordinate_intersection(
    const std::vector<double>& normal_values, double normal_scale,
    const std::vector<double>& direction_values, double direction_scale,
    const std::vector<double>& point_values, double point_scale,
    double plane_offset, const std::string& operation
) {
    const double plane_location = plane_offset / normal_scale;
    if (!std::isfinite(plane_location)) {
        return VectorExprListResult::failure(
            CasErrc::NumericFailure,
            "line-plane intersection coordinate is outside the supported domain",
            operation);
    }
    double spatial_scale = std::max(point_scale, std::abs(plane_location));
    if (spatial_scale == 0.0) spatial_scale = 1.0;

    double unit_normal[3];
    double unit_direction[3];
    double scaled_point[3];
    for (size_t i = 0; i < 3; ++i) {
        unit_normal[i] = normal_values[i] / normal_scale;
        unit_direction[i] = direction_values[i] / direction_scale;
        scaled_point[i] = point_values[i] / spatial_scale;
    }
    const double denominator = fused_dot_three(unit_normal, unit_direction);
    if (denominator == 0.0) {
        return VectorExprListResult::failure(
            CasErrc::DomainError,
            "line is parallel to plane, so no unique intersection exists",
            operation);
    }
    const double point_dot = fused_dot_three(unit_normal, scaled_point);
    const double scaled_offset = plane_location / spatial_scale;
    const double parameter = (scaled_offset - point_dot) / denominator;

    std::vector<std::shared_ptr<SymbolicExpr>> intersection;
    intersection.reserve(3);
    for (size_t i = 0; i < 3; ++i) {
        const double scaled_coordinate = std::fma(
            parameter, unit_direction[i], scaled_point[i]);
        const double coordinate = scaled_coordinate * spatial_scale;
        if (!std::isfinite(coordinate)) {
            return VectorExprListResult::failure(
                CasErrc::NumericFailure,
                "line-plane intersection coordinate is outside the supported domain",
                operation);
        }
        intersection.push_back(SymbolicExpr::number(coordinate));
    }
    return VectorExprListResult::success(std::move(intersection));
}

VectorExprListResult symbolic_line_plane_intersection(
    const std::vector<std::shared_ptr<SymbolicExpr>>& normal,
    const std::vector<std::shared_ptr<SymbolicExpr>>& direction,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    const std::shared_ptr<SymbolicExpr>& offset,
    ComputationContext& context, const std::string& operation
) {
    auto n_dot_a = vector_dot_checked(normal, point, context);
    if (!n_dot_a) return VectorExprListResult::failure(n_dot_a.error());
    auto n_dot_b = vector_dot_checked(normal, direction, context);
    if (!n_dot_b) return VectorExprListResult::failure(n_dot_b.error());
    auto denom = n_dot_b.value()->simplify();
    auto nonzero = checked_nonzero_numeric_or_exact(
        denom, context, operation,
        "line is parallel to plane, so no unique intersection exists",
        "line-plane denominator nonzero condition cannot be verified");
    if (!nonzero) return VectorExprListResult::failure(nonzero.error());

    auto t = SymbolicExpr::divide(
        SymbolicExpr::add(
            offset,
            SymbolicExpr::multiply(SymbolicExpr::number(-1), n_dot_a.value())),
        denom);
    std::vector<std::shared_ptr<SymbolicExpr>> intersection;
    intersection.reserve(point.size());
    for (size_t i = 0; i < point.size(); ++i) {
        auto coord = SymbolicExpr::add(
            point[i], SymbolicExpr::multiply(t, direction[i]));
        auto simplified = simplify_checked(
            coord, operation,
            "line-plane intersection coordinate is outside the supported domain");
        if (!simplified) {
            return VectorExprListResult::failure(simplified.error());
        }
        intersection.push_back(std::move(simplified.value()));
    }
    return VectorExprListResult::success(std::move(intersection));
}

ExpressionResult symbolic_point_plane_distance(
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    const std::vector<std::shared_ptr<SymbolicExpr>>& normal,
    const std::shared_ptr<SymbolicExpr>& offset,
    ComputationContext& context, const std::string& operation
) {
    auto n_dot_r = vector_dot_checked(normal, point, context);
    if (!n_dot_r) return ExpressionResult::failure(n_dot_r.error());
    auto norm_sq = vector_dot_checked(normal, normal, context);
    if (!norm_sq) return ExpressionResult::failure(norm_sq.error());

    auto diff = SymbolicExpr::add(
        n_dot_r.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), offset));
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
}

}

VectorExprListResult line_plane_intersection_checked(
    const LineSymbolic& line,
    const PlaneSymbolic& plane,
    ComputationContext& context
) {
    const std::string operation = "line_plane_intersection";
    try {
        auto valid = validate_intersection_inputs(line, plane, context, operation);
        if (!valid) {
            return VectorExprListResult::failure(valid.error());
        }

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
        if (!point_scale && interrupts_optional_numeric_evaluation(point_scale.error())) {
            return VectorExprListResult::failure(point_scale.error());
        }
        auto plane_offset = evaluate_numeric(*plane.d, NumericBindings{}, context);
        const bool numeric_plane_offset = finite_plane_offset(plane_offset);
        if (!plane_offset && interrupts_optional_numeric_evaluation(plane_offset.error())) {
            return VectorExprListResult::failure(plane_offset.error());
        }

        const bool scale_direction_products = direction_products_need_scaling(
            normal_scale.value(), direction_scale.value());
        const bool protect_point_dot =
            point_dot_needs_protection(normal_scale.value(), point_scale);
        const std::vector<std::shared_ptr<SymbolicExpr>>* effective_normal = &plane.normal;
        const std::vector<std::shared_ptr<SymbolicExpr>>* effective_direction = &line.direction;
        std::vector<std::shared_ptr<SymbolicExpr>> scaled_normal;
        std::vector<std::shared_ptr<SymbolicExpr>> scaled_direction;
        std::shared_ptr<SymbolicExpr> effective_d = plane.d;
        if (scale_direction_products || protect_point_dot) {
            scale_symbolic_components(
                normal_values, normal_scale.value(), protect_point_dot, scaled_normal);
            effective_normal = &scaled_normal;
            effective_d = scale_plane_offset(plane.d, normal_scale.value(), protect_point_dot);
        }
        if (scale_direction_products) {
            scale_symbolic_components(
                direction_values, direction_scale.value(), false, scaled_direction);
            effective_direction = &scaled_direction;
        }

        if (coordinate_update_needs_protection(point_scale, numeric_plane_offset)) {
            return protected_coordinate_intersection(
                normal_values, normal_scale.value(),
                direction_values, direction_scale.value(),
                point_values, point_scale.value(), plane_offset.value().value, operation);
        }
        return symbolic_line_plane_intersection(
            *effective_normal, *effective_direction, line.point, effective_d,
            context, operation);
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
        if (!point_scale && interrupts_optional_numeric_evaluation(point_scale.error())) {
            return ExpressionResult::failure(point_scale.error());
        }

        const double safe_low =
            std::sqrt(std::numeric_limits<double>::min());
        const double safe_high =
            std::sqrt(std::numeric_limits<double>::max() / 3.0);
        const bool protect_point_dot =
            point_dot_needs_protection(normal_scale.value(), point_scale);
        const std::vector<std::shared_ptr<SymbolicExpr>>* effective_normal =
            &plane.normal;
        std::vector<std::shared_ptr<SymbolicExpr>> scaled_normal;
        std::shared_ptr<SymbolicExpr> effective_d = plane.d;
        if (normal_scale.value() < safe_low ||
            normal_scale.value() > safe_high ||
            protect_point_dot) {
            scale_symbolic_components(
                normal_values, normal_scale.value(), protect_point_dot, scaled_normal);
            effective_normal = &scaled_normal;
            effective_d = scale_plane_offset(plane.d, normal_scale.value(), protect_point_dot);
        }

        return symbolic_point_plane_distance(
            point, *effective_normal, effective_d, context, operation);
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

}
