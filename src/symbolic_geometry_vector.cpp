#include "symbolic_geometry_vector.hpp"
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

void rescale_angle_component(
    double magnitude, double& scale, double& norm_squared,
    double& dot, double& correction)
{
    if (magnitude <= scale) return;
    const double ratio = scale == 0.0 ? 0.0 : scale / magnitude;
    norm_squared *= ratio * ratio;
    dot *= ratio;
    correction *= ratio;
    scale = magnitude;
}

void accumulate_angle_product(
    double left, double right, double& dot, double& correction)
{
    const double product = left * right;
    const double corrected = product - correction;
    const double next = dot + corrected;
    correction = (next - dot) - corrected;
    dot = next;
}

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

        rescale_angle_component(
            abs_a, scale_a, norm_a_scaled, dot_scaled, dot_correction);
        rescale_angle_component(
            abs_b, scale_b, norm_b_scaled, dot_scaled, dot_correction);
        if (scale_a != 0.0) {
            const double normalized = na / scale_a;
            norm_a_scaled += normalized * normalized;
        }
        if (scale_b != 0.0) {
            const double normalized = nb / scale_b;
            norm_b_scaled += normalized * normalized;
        }
        if (scale_a != 0.0 && scale_b != 0.0) {
            accumulate_angle_product(
                na / scale_a, nb / scale_b, dot_scaled, dot_correction);
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

}
