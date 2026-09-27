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

bool is_plane_budget_error(const Result<double>& result) {
    return !result &&
        (result.error().code == CasErrc::Cancelled ||
         result.error().code == CasErrc::ResourceLimit);
}

bool plane_point_differences_overflow(
    const std::vector<double>& p1,
    const std::vector<double>& p2,
    const std::vector<double>& p3) {
    for (size_t i = 0; i < 3; ++i) {
        if (!std::isfinite(p2[i] - p1[i]) ||
            !std::isfinite(p3[i] - p1[i])) {
            return true;
        }
    }
    return false;
}

void append_scaled_plane_edges(
    const std::vector<double>& p1,
    const std::vector<double>& p2,
    const std::vector<double>& p3,
    double point_scale,
    std::vector<std::shared_ptr<SymbolicExpr>>& v1,
    std::vector<std::shared_ptr<SymbolicExpr>>& v2) {
    for (size_t i = 0; i < 3; ++i) {
        v1.push_back(SymbolicExpr::number(
            p2[i] / point_scale - p1[i] / point_scale));
        v2.push_back(SymbolicExpr::number(
            p3[i] / point_scale - p1[i] / point_scale));
    }
}

Result<void> append_symbolic_plane_edges(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p3,
    std::vector<std::shared_ptr<SymbolicExpr>>& v1,
    std::vector<std::shared_ptr<SymbolicExpr>>& v2,
    const std::string& operation) {
    for (size_t i = 0; i < 3; ++i) {
        auto left = simplify_checked(
            SymbolicExpr::add(
                p2[i],
                SymbolicExpr::multiply(SymbolicExpr::number(-1), p1[i])),
            operation,
            "plane edge vector is outside the supported domain");
        if (!left) return Result<void>::failure(left.error());
        auto right = simplify_checked(
            SymbolicExpr::add(
                p3[i],
                SymbolicExpr::multiply(SymbolicExpr::number(-1), p1[i])),
            operation,
            "plane edge vector is outside the supported domain");
        if (!right) return Result<void>::failure(right.error());
        v1.push_back(std::move(left.value()));
        v2.push_back(std::move(right.value()));
    }
    return Result<void>::success();
}

Result<void> build_plane_edges(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p3,
    std::vector<std::shared_ptr<SymbolicExpr>>& v1,
    std::vector<std::shared_ptr<SymbolicExpr>>& v2,
    ComputationContext& context,
    const std::string& operation) {
    std::vector<double> p1_values;
    std::vector<double> p2_values;
    std::vector<double> p3_values;
    auto p1_scale = numeric_vector_scale_checked(
        p1, context, operation, &p1_values);
    auto p2_scale = numeric_vector_scale_checked(
        p2, context, operation, &p2_values);
    auto p3_scale = numeric_vector_scale_checked(
        p3, context, operation, &p3_values);
    if (is_plane_budget_error(p1_scale)) {
        return Result<void>::failure(p1_scale.error());
    }
    if (is_plane_budget_error(p2_scale)) {
        return Result<void>::failure(p2_scale.error());
    }
    if (is_plane_budget_error(p3_scale)) {
        return Result<void>::failure(p3_scale.error());
    }

    const bool numeric_points =
        p1_scale.has_value() && p2_scale.has_value() && p3_scale.has_value();
    bool scale_differences = false;
    double point_scale = 0.0;
    if (numeric_points) {
        point_scale = std::max(
            p1_scale.value(), std::max(p2_scale.value(), p3_scale.value()));
        scale_differences = plane_point_differences_overflow(
            p1_values, p2_values, p3_values);
    }
    v1.reserve(3);
    v2.reserve(3);
    if (scale_differences) {
        append_scaled_plane_edges(
            p1_values, p2_values, p3_values, point_scale, v1, v2);
        return Result<void>::success();
    }
    return append_symbolic_plane_edges(p1, p2, p3, v1, v2, operation);
}

Result<bool> scale_plane_edges_for_cross(
    std::vector<std::shared_ptr<SymbolicExpr>>& v1,
    std::vector<std::shared_ptr<SymbolicExpr>>& v2,
    ComputationContext& context,
    const std::string& operation) {
    std::vector<double> v1_values;
    std::vector<double> v2_values;
    auto v1_scale = numeric_vector_scale_checked(
        v1, context, operation, &v1_values);
    auto v2_scale = numeric_vector_scale_checked(
        v2, context, operation, &v2_values);
    const bool numeric_edges = v1_scale.has_value() && v2_scale.has_value();
    if (!numeric_edges) {
        if (is_plane_budget_error(v1_scale) || is_plane_budget_error(v2_scale)) {
            return Result<bool>::failure(
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
    return Result<bool>::success(numeric_edges);
}

Result<void> check_plane_normal_nonzero(
    const std::vector<std::shared_ptr<SymbolicExpr>>& normal,
    bool numeric_edges,
    ComputationContext& context,
    const std::string& operation) {
    Result<void> nonzero = numeric_edges
        ? checked_nonzero_numeric_vector(
              normal, context, operation,
              "three points are collinear, so no unique plane exists",
              "plane normal nonzero condition cannot be verified")
        : Result<void>::failure(
              CasErrc::Inconclusive,
              "plane normal nonzero condition cannot be verified",
              operation);
    if (!numeric_edges) {
        auto norm_sq = vector_dot_checked(normal, normal, context);
        if (!norm_sq) return Result<void>::failure(norm_sq.error());
        nonzero = checked_nonzero_numeric_or_exact(
            norm_sq.value(), context, operation,
            "three points are collinear, so no unique plane exists",
            "plane normal nonzero condition cannot be verified");
    }
    return nonzero;
}

PlaneSymbolicResult assemble_plane(
    std::vector<std::shared_ptr<SymbolicExpr>>&& normal,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    ComputationContext& context,
    const std::string& operation) {
    auto d = vector_dot_checked(normal, point, context);
    if (!d) return PlaneSymbolicResult::failure(d.error());
    auto d_simplified = simplify_checked(
        d.value(), operation,
        "plane constant is outside the supported domain");
    if (!d_simplified) return PlaneSymbolicResult::failure(d_simplified.error());

    PlaneSymbolic plane;
    plane.normal = std::move(normal);
    for (auto& component : plane.normal) {
        auto simplified = simplify_checked(
            component, operation,
            "plane normal component is outside the supported domain");
        if (!simplified) return PlaneSymbolicResult::failure(simplified.error());
        component = std::move(simplified.value());
    }
    plane.d = std::move(d_simplified.value());
    return PlaneSymbolicResult::success(std::move(plane));
}

}

PlaneSymbolicResult plane_from_three_points_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p3,
    ComputationContext& context) {
    const std::string operation = "plane_from_three_points";
    try {
        auto v12 = validate_same_dimension_vectors(p1, p2, context, operation);
        if (!v12) {
            return PlaneSymbolicResult::failure(v12.error());
        }
        auto v13 = validate_same_dimension_vectors(p1, p3, context, operation);
        if (!v13) {
            return PlaneSymbolicResult::failure(v13.error());
        }
        if (p1.size() != 3) {
            return PlaneSymbolicResult::failure(
                CasErrc::InvalidArgument,
                "plane construction requires 3-dimensional points",
                operation);
        }
        auto step = context.consume_steps(8, operation);
        if (!step) {
            return PlaneSymbolicResult::failure(step.error());
        }

        std::vector<std::shared_ptr<SymbolicExpr>> v1, v2;
        auto edges = build_plane_edges(p1, p2, p3, v1, v2, context, operation);
        if (!edges) {
            return PlaneSymbolicResult::failure(edges.error());
        }
        auto numeric_edges = scale_plane_edges_for_cross(v1, v2, context, operation);
        if (!numeric_edges) {
            return PlaneSymbolicResult::failure(numeric_edges.error());
        }

        auto normal = vector_cross_checked(v1, v2, context);
        if (!normal) {
            return PlaneSymbolicResult::failure(normal.error());
        }
        auto nonzero = check_plane_normal_nonzero(
            normal.value(), numeric_edges.value(), context, operation);
        if (!nonzero) {
            return PlaneSymbolicResult::failure(nonzero.error());
        }
        return assemble_plane(std::move(normal.value()), p1, context, operation);
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

}
