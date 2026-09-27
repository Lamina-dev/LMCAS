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

Result<double> point_difference_scale_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& p1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& p2,
    ComputationContext& context,
    const std::string& operation,
    std::vector<double>& p1_values,
    std::vector<double>& p2_values)
{
    auto p1_scale = numeric_vector_scale_checked(
        p1, context, operation, &p1_values);
    auto p2_scale = numeric_vector_scale_checked(
        p2, context, operation, &p2_values);
    if (!p1_scale &&
        (p1_scale.error().code == CasErrc::Cancelled ||
         p1_scale.error().code == CasErrc::ResourceLimit)) {
        return Result<double>::failure(p1_scale.error());
    }
    if (!p2_scale &&
        (p2_scale.error().code == CasErrc::Cancelled ||
         p2_scale.error().code == CasErrc::ResourceLimit)) {
        return Result<double>::failure(p2_scale.error());
    }
    if (!p1_scale || !p2_scale) {
        return Result<double>::success(0.0);
    }
    const double scale = std::max(p1_scale.value(), p2_scale.value());
    for (size_t i = 0; i < 3; ++i) {
        if (!std::isfinite(p2_values[i] - p1_values[i])) {
            return Result<double>::success(scale);
        }
    }
    return Result<double>::success(0.0);
}

ExpressionResult skew_distance_from_directions(
    const std::vector<std::shared_ptr<SymbolicExpr>>& direction1,
    const std::vector<std::shared_ptr<SymbolicExpr>>& direction2,
    const std::vector<std::shared_ptr<SymbolicExpr>>& difference,
    double point_scale,
    ComputationContext& context,
    const std::string& operation)
{
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
    auto numerator = vector_dot_checked(difference, cross.value(), context);
    if (!numerator) return ExpressionResult::failure(numerator.error());
    auto abs_num = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<const SymbolicNode>>{
                LMCAS::detail::node(numerator.value())}));
    auto distance = SymbolicExpr::divide(abs_num, cross_norm);
    if (point_scale != 0.0) {
        distance = SymbolicExpr::multiply(
            distance, SymbolicExpr::number(point_scale));
    }
    return simplify_checked(
        distance, operation,
        "skew line distance is outside the supported domain");
}

}

ExpressionResult skew_lines_distance_checked(
    const LineSymbolic& l1,
    const LineSymbolic& l2,
    ComputationContext& context
) {
    const std::string operation = "skew_lines_distance";
    try {
        auto l1_valid = validate_line_checked(l1, context, operation);
        if (!l1_valid) {
            return ExpressionResult::failure(l1_valid.error());
        }
        auto l2_valid = validate_line_checked(l2, context, operation);
        if (!l2_valid) {
            return ExpressionResult::failure(l2_valid.error());
        }
        auto step = context.consume_steps(8, operation);
        if (!step) {
            return ExpressionResult::failure(step.error());
        }

        std::vector<double> direction1_values;
        std::vector<double> direction2_values;
        auto scale1 = numeric_vector_scale_checked(
            l1.direction, context, operation, &direction1_values);
        if (!scale1) {
            return ExpressionResult::failure(scale1.error());
        }
        auto scale2 = numeric_vector_scale_checked(
            l2.direction, context, operation, &direction2_values);
        if (!scale2) {
            return ExpressionResult::failure(scale2.error());
        }

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
        auto difference_scale = point_difference_scale_checked(
            l1.point, l2.point, context, operation, point1_values, point2_values);
        if (!difference_scale) {
            return ExpressionResult::failure(difference_scale.error());
        }
        const double point_scale = difference_scale.value();

        std::vector<std::shared_ptr<SymbolicExpr>> a2_minus_a1;
        a2_minus_a1.reserve(3);
        for (size_t i = 0; i < l1.point.size(); ++i) {
            if (point_scale != 0.0) {
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
        return skew_distance_from_directions(
            direction1, direction2, a2_minus_a1, point_scale, context, operation);
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
        auto difference_scale = point_difference_scale_checked(
            p1, p2, context, operation, p1_values, p2_values);
        if (!difference_scale) {
            return LineSymbolicResult::failure(difference_scale.error());
        }
        const double point_scale = difference_scale.value();

        LineSymbolic line;
        line.point = p1;
        line.direction.reserve(3);
        for (size_t i = 0; i < p1.size(); ++i) {
            if (point_scale != 0.0) {
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

}
