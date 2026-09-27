#include "internal/symbolic_matrix_support.hpp"

#include <cmath>
#include <new>

namespace LMCAS {
using namespace detail::matrix_api;
namespace {
ExpressionResult unsupported_matrix_dimension(const std::string& operation)
{
    return ExpressionResult::failure(
        CasErrc::UnsupportedExpression,
        "only two-dimensional transformation matrices are supported", operation);
}

ExpressionResult invalid_matrix_dimension(const std::string& operation)
{
    return ExpressionResult::failure(
        CasErrc::DimensionMismatch,
        "only two-dimensional transformation matrices are supported", operation);
}
}

ExpressionResult matrix_rotation_checked(double theta, int dim,
                                         ComputationContext& context) {
    const std::string operation = "matrix_rotation";
    auto step = context.consume_steps(1, operation);
    if (!step) return ExpressionResult::failure(step.error());
    if (dim != 2) return invalid_matrix_dimension(operation);
    if (!std::isfinite(theta)) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument,
            "rotation angle must be finite", operation);
    }
    auto c = SymbolicExpr::number(std::cos(theta));
    auto s = SymbolicExpr::number(std::sin(theta));
    return require_matrix_result(
        SymbolicExpr::matrix(
            {{c, SymbolicExpr::multiply(SymbolicExpr::number(-1), s)},
             {s, c}}),
        operation);
}

ExpressionResult matrix_rotation_checked(double theta, int dim) {
    ComputationContext context;
    return matrix_rotation_checked(theta, dim, context);
}


ExpressionResult matrix_reflection_checked(double angle, int dim,
                                           ComputationContext& context) {
    const std::string operation = "matrix_reflection";
    auto step = context.consume_steps(1, operation);
    if (!step) return ExpressionResult::failure(step.error());
    if (dim == 2) {
        if (!std::isfinite(angle)) {
            return ExpressionResult::failure(CasErrc::InvalidArgument,
                "reflection angle must be finite", operation);
        }
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);
        try {
            auto c = SymbolicExpr::number(cosine * cosine - sine * sine);
            auto s = SymbolicExpr::number(2.0 * sine * cosine);
            return require_matrix_result(
                SymbolicExpr::matrix({{c, s}, {s, SymbolicExpr::multiply(SymbolicExpr::number(-1), c)}}),
                operation);
        } catch (const CasError& error) {
            return ExpressionResult::failure(error);
        } catch (const std::bad_alloc&) {
            return ExpressionResult::failure(CasErrc::ResourceLimit,
                "allocation failed while constructing reflection matrix", operation);
        } catch (const std::exception& error) {
            return ExpressionResult::failure(CasErrc::InternalInvariant,
                error.what(), operation);
        }
    }
    return unsupported_matrix_dimension(operation);
}

ExpressionResult matrix_reflection_checked(double angle, int dim) {
    ComputationContext context;
    return matrix_reflection_checked(angle, dim, context);
}


ExpressionResult matrix_scaling_checked(double sx, double sy, int dim,
                                        ComputationContext& context) {
    const std::string operation = "matrix_scaling";
    auto step = context.consume_steps(1, operation);
    if (!step) return ExpressionResult::failure(step.error());
    if (dim != 2) return invalid_matrix_dimension(operation);
    if (!std::isfinite(sx) || !std::isfinite(sy)) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument,
            "scaling factors must be finite", operation);
    }
    return require_matrix_result(
        SymbolicExpr::matrix(
            {{SymbolicExpr::number(sx), SymbolicExpr::number(0)},
             {SymbolicExpr::number(0), SymbolicExpr::number(sy)}}),
        operation);
}

ExpressionResult matrix_scaling_checked(double sx, double sy, int dim) {
    ComputationContext context;
    return matrix_scaling_checked(sx, sy, dim, context);
}

}
