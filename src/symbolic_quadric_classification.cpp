#include "internal/symbolic_quadric_classification.hpp"

#include "numeric_evaluation.hpp"
#include "internal/symbolic_ast.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace LMCAS::detail {
namespace {

struct SymmetricEigen3 {
    double values[3];
    double vectors[3][3];
    double coefficient_scale;
    double zero_tolerance;
    bool unresolved_rank;
};

double largest_off_diagonal(const double matrix[3][3], size_t& p, size_t& q)
{
    p = 0;
    q = 1;
    double off_diagonal = std::abs(matrix[p][q]);
    for (size_t row = 0; row < 3; ++row) {
        for (size_t column = row + 1; column < 3; ++column) {
            const double candidate = std::abs(matrix[row][column]);
            if (candidate > off_diagonal) {
                off_diagonal = candidate;
                p = row;
                q = column;
            }
        }
    }
    return off_diagonal;
}

void rotate_principal_axes(
    double matrix[3][3], double vectors[3][3], size_t p, size_t q)
{
    const double angle = 0.5 * std::atan2(
        2.0 * matrix[p][q], matrix[q][q] - matrix[p][p]);
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    const double pp = matrix[p][p];
    const double pq = matrix[p][q];
    const double qq = matrix[q][q];

    for (size_t axis = 0; axis < 3; ++axis) {
        if (axis == p || axis == q) continue;
        const double axis_p = matrix[axis][p];
        const double axis_q = matrix[axis][q];
        matrix[axis][p] = cosine * axis_p - sine * axis_q;
        matrix[p][axis] = matrix[axis][p];
        matrix[axis][q] = sine * axis_p + cosine * axis_q;
        matrix[q][axis] = matrix[axis][q];
    }
    matrix[p][p] = cosine * cosine * pp -
                   2.0 * sine * cosine * pq +
                   sine * sine * qq;
    matrix[q][q] = sine * sine * pp +
                   2.0 * sine * cosine * pq +
                   cosine * cosine * qq;
    matrix[p][q] = 0.0;
    matrix[q][p] = 0.0;

    for (size_t axis = 0; axis < 3; ++axis) {
        const double vector_p = vectors[axis][p];
        const double vector_q = vectors[axis][q];
        vectors[axis][p] = cosine * vector_p - sine * vector_q;
        vectors[axis][q] = sine * vector_p + cosine * vector_q;
    }
}

SymmetricEigen3 diagonalize_symmetric_3x3(const double input[3][3])
{
    double matrix[3][3];
    SymmetricEigen3 result{{0.0, 0.0, 0.0},
                           {{1.0, 0.0, 0.0},
                            {0.0, 1.0, 0.0},
                            {0.0, 0.0, 1.0}},
                           0.0,
                           0.0,
                           false};
    for (size_t row = 0; row < 3; ++row) {
        for (size_t column = 0; column < 3; ++column) {
            result.coefficient_scale = std::max(
                result.coefficient_scale, std::abs(input[row][column]));
        }
    }
    if (result.coefficient_scale == 0.0) return result;
    for (size_t row = 0; row < 3; ++row) {
        for (size_t column = 0; column < 3; ++column) {
            matrix[row][column] =
                input[row][column] / result.coefficient_scale;
        }
    }
    result.zero_tolerance =
        256.0 * std::numeric_limits<double>::epsilon();

    for (size_t iteration = 0; iteration < 64; ++iteration) {
        size_t p;
        size_t q;
        if (largest_off_diagonal(matrix, p, q) <= result.zero_tolerance) break;
        rotate_principal_axes(matrix, result.vectors, p, q);
    }
    for (size_t axis = 0; axis < 3; ++axis) {
        result.values[axis] = matrix[axis][axis];
        if (std::abs(result.values[axis]) <= result.zero_tolerance) {
            result.unresolved_rank =
                result.unresolved_rank || result.values[axis] != 0.0;
            result.values[axis] = 0.0;
        }
    }
    return result;
}

Result<double> finite_numeric_coeff(
    const std::shared_ptr<SymbolicExpr>& expression,
    ComputationContext& context,
    const std::string& operation,
    const std::string& message)
{
    auto simplified = expression ? expression->simplify() : nullptr;
    if (!simplified || !node(simplified)) {
        return Result<double>::failure(
            CasErrc::Inconclusive, message, operation);
    }
    auto evaluated = evaluate_numeric(
        *simplified, NumericBindings{}, context);
    if (!evaluated) {
        if (evaluated.error().code == CasErrc::Cancelled ||
            evaluated.error().code == CasErrc::ResourceLimit) {
            return Result<double>::failure(evaluated.error());
        }
        return Result<double>::failure(
            CasErrc::Inconclusive, message, operation);
    }
    if (!evaluated.value().is_finite() ||
        !std::isfinite(evaluated.value().value)) {
        return Result<double>::failure(
            CasErrc::Inconclusive, message, operation);
    }
    return Result<double>::success(evaluated.value().value);
}

Result<double> coefficient_of_square(
    const std::shared_ptr<SymbolicExpr>& expanded,
    const std::string& variable,
    ComputationContext& context,
    const std::string& operation)
{
    auto derivative = expanded->differentiate(variable);
    if (!derivative || !node(derivative)) {
        return Result<double>::failure(
            CasErrc::Inconclusive,
            "quadric second derivative is outside the supported domain",
            operation);
    }
    derivative = derivative->differentiate(variable);
    auto value = finite_numeric_coeff(
        derivative, context, operation,
        "quadric second derivative coefficient cannot be evaluated");
    if (!value) return value;
    return Result<double>::success(value.value() / 2.0);
}

Result<double> coefficient_of_linear(
    const std::shared_ptr<SymbolicExpr>& expanded,
    const std::string& variable,
    const std::vector<std::string>& vars,
    ComputationContext& context,
    const std::string& operation)
{
    auto derivative = expanded->differentiate(variable);
    if (!derivative || !node(derivative)) {
        return Result<double>::failure(
            CasErrc::Inconclusive,
            "quadric first derivative is outside the supported domain",
            operation);
    }
    for (const auto& coordinate : vars) {
        derivative = derivative->substitute(
            coordinate, SymbolicExpr::number(0));
        if (!derivative || !node(derivative)) {
            return Result<double>::failure(
                CasErrc::Inconclusive,
                "quadric linear coefficient cannot be substituted",
                operation);
        }
    }
    return finite_numeric_coeff(
        derivative, context, operation,
        "quadric linear coefficient cannot be evaluated");
}

bool extract_axis_coefficients(
    const std::shared_ptr<SymbolicExpr>& expanded,
    const std::vector<std::string>& vars,
    ComputationContext& context,
    const std::string& operation,
    double quadratic_matrix[3][3],
    double linear[3],
    CasError& error)
{
    Result<double> diagonal[3] = {
        coefficient_of_square(expanded, vars[0], context, operation),
        coefficient_of_square(expanded, vars[1], context, operation),
        coefficient_of_square(expanded, vars[2], context, operation)};
    Result<double> linear_coefficients[3] = {
        coefficient_of_linear(expanded, vars[0], vars, context, operation),
        coefficient_of_linear(expanded, vars[1], vars, context, operation),
        coefficient_of_linear(expanded, vars[2], vars, context, operation)};
    for (size_t axis = 0; axis < 3; ++axis) {
        if (!diagonal[axis]) {
            error = std::move(diagonal[axis].error());
            return false;
        }
        if (!linear_coefficients[axis]) {
            error = std::move(linear_coefficients[axis].error());
            return false;
        }
    }
    for (size_t axis = 0; axis < 3; ++axis) {
        quadratic_matrix[axis][axis] = diagonal[axis].value();
        linear[axis] = linear_coefficients[axis].value();
    }
    return true;
}

bool extract_mixed_coefficients(
    const std::shared_ptr<SymbolicExpr>& expanded,
    const std::vector<std::string>& vars,
    ComputationContext& context,
    const std::string& operation,
    double quadratic_matrix[3][3],
    CasError& error)
{
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i + 1; j < 3; ++j) {
            auto mixed = expanded->differentiate(vars[i]);
            if (!mixed || !node(mixed)) {
                error = {
                    CasErrc::Inconclusive,
                    "quadric mixed derivative is outside the supported domain",
                    operation};
                return false;
            }
            mixed = mixed->differentiate(vars[j]);
            auto mixed_value = finite_numeric_coeff(
                mixed, context, operation,
                "quadric mixed coefficient cannot be evaluated");
            if (!mixed_value) {
                error = std::move(mixed_value.error());
                return false;
            }
            quadratic_matrix[i][j] = mixed_value.value() * 0.5;
            quadratic_matrix[j][i] = quadratic_matrix[i][j];
        }
    }
    return true;
}

bool project_linear_coefficients(
    const SymmetricEigen3& principal,
    const double coefficients[3],
    double linear[3],
    double linear_tolerance[3],
    const std::string& operation,
    CasError& error)
{
    double projection_scale[3] = {0.0, 0.0, 0.0};
    for (size_t row = 0; row < 3; ++row) {
        const double normalized =
            coefficients[row] / principal.coefficient_scale;
        if (!std::isfinite(normalized)) {
            error = {
                CasErrc::Inconclusive,
                "quadric linear part is outside the supported scale",
                operation};
            return false;
        }
        for (size_t axis = 0; axis < 3; ++axis) {
            const double term = principal.vectors[row][axis] * normalized;
            projection_scale[axis] =
                std::max(projection_scale[axis], std::abs(term));
            linear[axis] = std::fma(
                principal.vectors[row][axis], normalized, linear[axis]);
        }
    }
    for (size_t axis = 0; axis < 3; ++axis) {
        linear_tolerance[axis] =
            256.0 * std::numeric_limits<double>::epsilon() *
            projection_scale[axis];
    }
    return true;
}

bool has_null_axis_linear_term(
    const SymmetricEigen3& principal,
    const double linear[3],
    const double linear_tolerance[3])
{
    for (size_t axis = 0; axis < 3; ++axis) {
        if (principal.values[axis] == 0.0 &&
            std::abs(linear[axis]) > linear_tolerance[axis]) {
            return true;
        }
    }
    return false;
}

Result<double> coefficient_of_constant(
    const std::shared_ptr<SymbolicExpr>& expanded,
    const std::vector<std::string>& vars,
    ComputationContext& context,
    const std::string& operation)
{
    auto constant = expanded;
    for (const auto& coordinate : vars) {
        constant = constant->substitute(
            coordinate, SymbolicExpr::number(0));
        if (!constant || !node(constant)) {
            return Result<double>::failure(
                CasErrc::Inconclusive,
                "quadric constant term cannot be substituted",
                operation);
        }
    }
    return finite_numeric_coeff(
        constant, context, operation,
        "quadric constant term cannot be evaluated");
}

struct CenteredQuadric {
    double constant;
    double zero_tolerance;
};

Result<CenteredQuadric> center_quadric(
    const SymmetricEigen3& principal,
    const double linear[3],
    double constant,
    const std::string& operation)
{
    double centered_constant = constant / principal.coefficient_scale;
    if (!std::isfinite(centered_constant)) {
        return Result<CenteredQuadric>::failure(
            CasErrc::Inconclusive,
            "quadric constant is outside the supported scale",
            operation);
    }
    double centered_scale = std::max(1.0, std::abs(centered_constant));
    for (size_t axis = 0; axis < 3; ++axis) {
        if (principal.values[axis] == 0.0) continue;
        const double factor =
            (linear[axis] / principal.values[axis]) * 0.25;
        const double shift = factor * linear[axis];
        if (!std::isfinite(factor) || !std::isfinite(shift)) {
            return Result<CenteredQuadric>::failure(
                CasErrc::Inconclusive,
                "quadric translation is outside the supported domain",
                operation);
        }
        centered_constant = std::fma(
            -factor, linear[axis], centered_constant);
        centered_scale = std::max(centered_scale, std::abs(shift));
    }
    const double zero_tolerance =
        256.0 * std::numeric_limits<double>::epsilon() * centered_scale;
    return Result<CenteredQuadric>::success({centered_constant, zero_tolerance});
}

VectorStringResult classify_rank_two_locus(
    int positive,
    int negative,
    double centered_constant,
    double zero_tolerance,
    const std::string& operation)
{
    if (std::abs(centered_constant) <= zero_tolerance) {
        return VectorStringResult::failure(
            CasErrc::Inconclusive,
            "rank-two quadric has a degenerate real locus",
            operation);
    }
    const bool indefinite = positive == 1 && negative == 1;
    const bool real_definite =
        (positive == 2 && centered_constant < -zero_tolerance) ||
        (negative == 2 && centered_constant > zero_tolerance);
    if (indefinite || real_definite) {
        return VectorStringResult::success("cylinder");
    }
    return VectorStringResult::failure(
        CasErrc::Inconclusive,
        "rank-two quadric has no real cylinder locus",
        operation);
}

VectorStringResult classify_definite_locus(
    const SymmetricEigen3& principal,
    int positive,
    int negative,
    double centered_constant,
    double zero_tolerance,
    const std::string& operation)
{
    const bool real_ellipsoid =
        (positive == 3 && centered_constant < -zero_tolerance) ||
        (negative == 3 && centered_constant > zero_tolerance);
    if (!real_ellipsoid) {
        return VectorStringResult::failure(
            CasErrc::Inconclusive,
            std::abs(centered_constant) <= zero_tolerance
                ? "definite quadric degenerates to a point"
                : "definite quadric has no real locus",
            operation);
    }
    const double min_value = std::min(
        principal.values[0],
        std::min(principal.values[1], principal.values[2]));
    const double max_value = std::max(
        principal.values[0],
        std::max(principal.values[1], principal.values[2]));
    if (max_value - min_value <= principal.zero_tolerance) {
        return VectorStringResult::success("sphere");
    }
    return VectorStringResult::success("ellipsoid");
}

VectorStringResult classify_centered_locus(
    const SymmetricEigen3& principal,
    int rank,
    double centered_constant,
    double zero_tolerance,
    const std::string& operation)
{
    const int positive =
        (principal.values[0] > 0.0) +
        (principal.values[1] > 0.0) +
        (principal.values[2] > 0.0);
    const int negative =
        (principal.values[0] < 0.0) +
        (principal.values[1] < 0.0) +
        (principal.values[2] < 0.0);
    if (rank == 2) {
        return classify_rank_two_locus(
            positive, negative, centered_constant, zero_tolerance, operation);
    }
    if (rank == 3) {
        if (negative == 0 || positive == 0) {
            return classify_definite_locus(
                principal, positive, negative, centered_constant,
                zero_tolerance, operation);
        }
        if (std::abs(centered_constant) <= zero_tolerance) {
            return VectorStringResult::success("cone");
        }
        return VectorStringResult::success("hyperboloid");
    }
    return VectorStringResult::failure(
        CasErrc::Inconclusive,
        "quadric is outside the currently classified support domain",
        operation);
}

VectorStringResult classify_translated_quadric(
    const SymmetricEigen3& principal,
    int rank,
    const double linear[3],
    const std::shared_ptr<SymbolicExpr>& expanded,
    const std::vector<std::string>& vars,
    ComputationContext& context,
    const std::string& operation)
{
    auto constant = coefficient_of_constant(expanded, vars, context, operation);
    if (!constant) return VectorStringResult::failure(std::move(constant.error()));
    auto centered = center_quadric(
        principal, linear, constant.value(), operation);
    if (!centered) return VectorStringResult::failure(std::move(centered.error()));
    const auto& value = centered.value();
    return classify_centered_locus(
        principal, rank, value.constant, value.zero_tolerance, operation);
}

}

VectorStringResult classify_quadric_impl(
    const SurfaceSymbolic& surf,
    ComputationContext& context,
    const std::string& operation)
{
    const auto& vars = surf.vars;
    auto expanded = surf.F->expand();
    if (!expanded || !node(expanded)) expanded = surf.F;

    double quadratic_matrix[3][3] = {};
    double linear_coefficients[3];
    CasError error;
    if (!extract_axis_coefficients(
            expanded, vars, context, operation, quadratic_matrix,
            linear_coefficients, error)) {
        return VectorStringResult::failure(std::move(error));
    }
    if (!extract_mixed_coefficients(
            expanded, vars, context, operation, quadratic_matrix, error)) {
        return VectorStringResult::failure(std::move(error));
    }

    const SymmetricEigen3 principal =
        diagonalize_symmetric_3x3(quadratic_matrix);
    if (principal.coefficient_scale == 0.0) {
        return VectorStringResult::failure(
            CasErrc::Inconclusive,
            "quadric has no nonzero quadratic part",
            operation);
    }
    if (principal.unresolved_rank) {
        return VectorStringResult::failure(
            CasErrc::Inconclusive,
            "quadric rank is unresolved within the eigensolver error bound",
            operation);
    }

    double linear[3] = {0.0, 0.0, 0.0};
    double linear_tolerance[3];
    if (!project_linear_coefficients(
            principal, linear_coefficients, linear, linear_tolerance,
            operation, error)) {
        return VectorStringResult::failure(std::move(error));
    }
    const int rank =
        (principal.values[0] != 0.0) +
        (principal.values[1] != 0.0) +
        (principal.values[2] != 0.0);

    if (rank == 1) {
        if (has_null_axis_linear_term(principal, linear, linear_tolerance)) {
            return VectorStringResult::success("cylinder");
        }
        return VectorStringResult::failure(
            CasErrc::Inconclusive,
            "rank-one quadric is degenerate or has no real locus",
            operation);
    }
    if (rank == 2 &&
        has_null_axis_linear_term(principal, linear, linear_tolerance)) {
        return VectorStringResult::success("paraboloid");
    }
    return classify_translated_quadric(
        principal, rank, linear, expanded, vars, context, operation);
}

}
