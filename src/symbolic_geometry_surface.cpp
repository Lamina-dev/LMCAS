#include "internal/expression_transform.hpp"
#include "symbolic_geometry_surface.hpp"
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

ExpressionResult differentiate_surface_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const std::string& operation, ComputationContext& context)
try {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    try {
        auto derivative = detail::make_expression_ptr(detail::differentiate_node(detail::node(expr), var, transform_budget));
        return simplify_checked(
            derivative, operation,
            "surface derivative is outside the supported domain", context);
    } catch (const std::length_error& error) { return ExpressionResult::failure(CasErrc::ResourceLimit, error.what(), operation); } catch (const std::exception&) {
        return ExpressionResult::failure(
            CasErrc::InternalInvariant,
            "surface derivative is outside the supported domain",
            operation);
    }
} catch (const CasError& error) {
    return ExpressionResult::failure(error);
} catch (const std::bad_alloc&) {
    return ExpressionResult::failure(CasErrc::ResourceLimit, "allocation failed", operation);
} catch (const std::length_error& error) {
    return ExpressionResult::failure(CasErrc::ResourceLimit, error.what(), operation);
} catch (const std::exception& error) {
    return ExpressionResult::failure(CasErrc::InternalInvariant, error.what(), operation);
}

VectorExprListResult surface_gradient_at_point_checked(
    const SurfaceSymbolic& surf,
    const std::vector<std::shared_ptr<SymbolicExpr>>& point,
    const std::string& operation, ComputationContext& context)
try {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    std::vector<std::shared_ptr<SymbolicExpr>> gradient;
    gradient.reserve(surf.vars.size());
    for (size_t i = 0; i < surf.vars.size(); ++i) {
        auto derivative = differentiate_surface_checked(surf.F, surf.vars[i], operation, context);
        if (!derivative) return VectorExprListResult::failure(derivative.error());
        auto substituted = derivative.value();
        for (size_t j = 0; j < surf.vars.size(); ++j) {
            substituted = detail::make_expression_ptr(detail::substitute_node(detail::node(substituted), surf.vars[j], detail::node(point[j]), transform_budget));
            if (!substituted || !LMCAS::detail::node(substituted)) {
                return VectorExprListResult::failure(
                    CasErrc::Inconclusive,
                    "surface derivative substitution is outside the supported domain",
                    operation);
            }
        }
        auto simplified = simplify_checked(
            substituted, operation,
            "surface derivative substitution is outside the supported domain", context);
        if (!simplified) return VectorExprListResult::failure(simplified.error());
        gradient.push_back(std::move(simplified.value()));
    }
    return VectorExprListResult::success(std::move(gradient));
} catch (const CasError& error) {
    return VectorExprListResult::failure(error);
} catch (const std::bad_alloc&) {
    return VectorExprListResult::failure(CasErrc::ResourceLimit, "allocation failed", operation);
} catch (const std::length_error& error) {
    return VectorExprListResult::failure(CasErrc::ResourceLimit, error.what(), operation);
} catch (const std::exception& error) {
    return VectorExprListResult::failure(CasErrc::InternalInvariant, error.what(), operation);
}

Result<double> surface_gradient_scale_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& gradient,
    ComputationContext& context,
    const std::string& operation,
    std::vector<double>& values,
    const std::string& singular_message)
{
    auto scaled = numeric_vector_scale_checked(
        gradient, context, operation, &values);
    if (!scaled) {
        if (scaled.error().code == CasErrc::Cancelled ||
            scaled.error().code == CasErrc::ResourceLimit) {
            return Result<double>::failure(scaled.error());
        }
        return Result<double>::failure(
            CasErrc::Inconclusive,
            "surface gradient nonzero condition cannot be verified",
            operation);
    }
    if (scaled.value() == 0.0) {
        return Result<double>::failure(
            CasErrc::DomainError, singular_message, operation);
    }
    return scaled;
}

PlaneSymbolicResult scaled_tangent_plane(
    const std::vector<double>& gradient_values,
    double gradient_scale,
    const std::vector<double>& point_values,
    double point_scale,
    const std::string& operation)
{
    double normalized_dot = 0.0;
    for (size_t i = 0; i < point_values.size(); ++i) {
        normalized_dot = std::fma(
            gradient_values[i] / gradient_scale,
            point_values[i] / point_scale,
            normalized_dot);
    }
    const double plane_constant = normalized_dot * point_scale;
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
            SymbolicExpr::number(component / gradient_scale));
    }
    plane.d = SymbolicExpr::number(plane_constant);
    return PlaneSymbolicResult::success(std::move(plane));
}

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

    auto gradient = surface_gradient_at_point_checked(surf, point, operation, context);
    if (!gradient) return gradient;
    std::vector<double> gradient_values;
    auto scaled = surface_gradient_scale_checked(
        gradient.value(), context, operation, gradient_values,
        "surface normal is undefined at a singular point");
    if (!scaled) return VectorExprListResult::failure(scaled.error());

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
    if (!valid) {
        return PlaneSymbolicResult::failure(valid.error());
    }
    auto step = context.consume_steps(surf.vars.size() * 4 + 4, operation);
    if (!step) {
        return PlaneSymbolicResult::failure(step.error());
    }

    auto gradient = surface_gradient_at_point_checked(surf, point, operation, context);
    if (!gradient) {
        return PlaneSymbolicResult::failure(gradient.error());
    }
    std::vector<double> gradient_values;
    auto scaled = surface_gradient_scale_checked(
        gradient.value(), context, operation, gradient_values,
        "tangent plane is undefined at a singular point");
    if (!scaled) {
        return PlaneSymbolicResult::failure(scaled.error());
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
        return scaled_tangent_plane(
            gradient_values, scaled.value(), point_values,
            point_scale.value(), operation);
    }

    auto dot = vector_dot_checked(gradient.value(), point, context);
    if (!dot) {
        return PlaneSymbolicResult::failure(dot.error());
    }
    auto d = simplify_checked(
        dot.value(), operation,
        "tangent plane constant is outside the supported domain", context);
    if (!d) {
        return PlaneSymbolicResult::failure(d.error());
    }

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
