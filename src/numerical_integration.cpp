/**
 * @file numerical_integration.cpp
 * @brief 数值积分实现.
 */
#include "numerical_integration.hpp"

#include "internal/lmmc_lifecycle.hpp"
#include "lmmc/quadrature.h"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <string>

namespace LMCAS {
namespace {

Result<double> checked_finite_value(
    Result<ApproxReal> value,
    const std::string& operation) {
    if (!value) return Result<double>::failure(value.error());
    if (!value.value().is_finite() || !std::isfinite(value.value().value)) {
        return Result<double>::failure(
            CasErrc::NumericFailure,
            "numeric evaluation did not produce a finite double",
            operation);
    }
    return Result<double>::success(value.value().value);
}

Result<double> eval_bound_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    ComputationContext& context,
    const std::string& operation) {
    if (!expr) {
        return Result<double>::failure(
            CasErrc::InvalidArgument, "integration bound is null", operation);
    }
    return checked_finite_value(
        evaluate_numeric(*expr, NumericBindings{}, context), operation);
}

Result<void> validate_quadrature_args(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    const std::string& operation) {
    if (!function || !lower || !upper) {
        return Result<void>::failure(
            CasErrc::InvalidArgument,
            "integrand and bounds must be non-null",
            operation);
    }
    if (variable.empty()) {
        return Result<void>::failure(
            CasErrc::InvalidArgument,
            "integration variable must not be empty",
            operation);
    }
    return Result<void>::success();
}

struct NormalizedInterval {
    double lower = 0.0;
    double upper = 0.0;
    double sign = 1.0;
};

Result<NormalizedInterval> evaluate_interval(
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    ComputationContext& context,
    const std::string& operation) {
    auto lower_value = eval_bound_checked(lower, context, operation);
    if (!lower_value) {
        return Result<NormalizedInterval>::failure(lower_value.error());
    }
    auto upper_value = eval_bound_checked(upper, context, operation);
    if (!upper_value) {
        return Result<NormalizedInterval>::failure(upper_value.error());
    }
    NormalizedInterval interval{lower_value.value(), upper_value.value(), 1.0};
    if (interval.lower > interval.upper) {
        std::swap(interval.lower, interval.upper);
        interval.sign = -1.0;
    }
    return Result<NormalizedInterval>::success(interval);
}

struct QuadratureCallback {
    const std::shared_ptr<SymbolicExpr>* function = nullptr;
    const std::string* variable = nullptr;
    ComputationContext* context = nullptr;
    const char* operation = nullptr;
    std::optional<CasError> error;
};

double evaluate_quadrature_callback(double x, void* user_data) {
    auto& callback = *static_cast<QuadratureCallback*>(user_data);
    if (callback.error) return std::numeric_limits<double>::quiet_NaN();
    NumericBindings bindings;
    bindings.emplace(*callback.variable, x);
    auto value = checked_finite_value(
        evaluate_numeric(**callback.function, bindings, *callback.context),
        callback.operation);
    if (!value) {
        callback.error = value.error();
        return std::numeric_limits<double>::quiet_NaN();
    }
    return value.value();
}

CasErrc quadrature_error_code(lmmc_status_t status) {
    if (status == LMMC_STATUS_INVALID_ARGUMENT) return CasErrc::InvalidArgument;
    if (status == LMMC_STATUS_ALLOCATION_FAILED ||
        status == LMMC_STATUS_REFERENCE_LIMIT ||
        status == LMMC_STATUS_WARNING_MAX_DEPTH) {
        return CasErrc::ResourceLimit;
    }
    return CasErrc::NumericFailure;
}

template <typename T>
Result<T> quadrature_failure(lmmc_status_t status, const std::string& operation) {
    return Result<T>::failure(
        quadrature_error_code(status), lmmc_status_string(status), operation);
}

QuadratureEstimate zero_quadrature_result() {
    return {};
}

Result<QuadratureEstimate> estimated_quadrature_result(
    double value,
    double estimated_error,
    double sign,
    const std::string& operation) {
    value *= sign;
    estimated_error = std::abs(estimated_error);
    if (!std::isfinite(value) || !std::isfinite(estimated_error)) {
        return Result<QuadratureEstimate>::failure(
            CasErrc::NumericFailure,
            "quadrature accumulation or error estimate is non-finite",
            operation);
    }
    QuadratureEstimate result;
    result.value = value;
    result.estimated_error = estimated_error;
    return Result<QuadratureEstimate>::success(result);
}

}

Result<QuadratureEstimate> quadrature_simpson_numeric(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    ComputationContext& context,
    int subdivisions) {
    constexpr const char* operation = "quadrature_simpson";
    auto valid = validate_quadrature_args(
        function, variable, lower, upper, operation);
    if (!valid) return Result<QuadratureEstimate>::failure(valid.error());
    if (subdivisions <= 0 || subdivisions % 2 != 0) {
        return Result<QuadratureEstimate>::failure(
            CasErrc::InvalidArgument,
            "Simpson subinterval count must be positive and even",
            operation);
    }
    if (subdivisions > std::numeric_limits<int>::max() / 2) {
        return Result<QuadratureEstimate>::failure(
            CasErrc::ResourceLimit,
            "Simpson refinement count overflows",
            operation);
    }
    auto interval = evaluate_interval(lower, upper, context, operation);
    if (!interval) return Result<QuadratureEstimate>::failure(interval.error());
    if (interval.value().lower == interval.value().upper) {
        return Result<QuadratureEstimate>::success(zero_quadrature_result());
    }
    const auto refined_subdivisions =
        static_cast<std::size_t>(subdivisions) * 2;
    const auto samples =
        static_cast<std::size_t>(subdivisions) + 1 +
        refined_subdivisions + 1;
    auto budget = context.consume_steps(samples, operation);
    if (!budget) return Result<QuadratureEstimate>::failure(budget.error());

    detail::ensure_lmmc_lifecycle();
    QuadratureCallback callback{
        &function, &variable, &context, operation, std::nullopt};
    double value = 0.0;
    auto status = lmmc_quad_simpson(
        evaluate_quadrature_callback, &callback,
        interval.value().lower, interval.value().upper,
        static_cast<std::size_t>(subdivisions), &value);
    if (callback.error) return Result<QuadratureEstimate>::failure(*callback.error);
    if (status != LMMC_STATUS_OK) {
        return quadrature_failure<QuadratureEstimate>(status, operation);
    }

    double refined_value = 0.0;
    status = lmmc_quad_simpson(
        evaluate_quadrature_callback, &callback,
        interval.value().lower, interval.value().upper,
        refined_subdivisions, &refined_value);
    if (callback.error) return Result<QuadratureEstimate>::failure(*callback.error);
    if (status != LMMC_STATUS_OK) {
        return quadrature_failure<QuadratureEstimate>(status, operation);
    }
    const double correction = (refined_value - value) / 15.0;
    const double extrapolated_value = refined_value + correction;
    return estimated_quadrature_result(
        extrapolated_value, std::abs(correction),
        interval.value().sign, operation);
}

Result<QuadratureEstimate> quadrature_simpson_numeric(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    int subdivisions) {
    ComputationContext context;
    return quadrature_simpson_numeric(
        function, variable, lower, upper, context, subdivisions);
}

Result<QuadratureEstimate> quadrature_gaussian_numeric(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    ComputationContext& context,
    int order) {
    constexpr const char* operation = "quadrature_gaussian";
    auto valid = validate_quadrature_args(
        function, variable, lower, upper, operation);
    if (!valid) return Result<QuadratureEstimate>::failure(valid.error());
    if (order < 1 || order > 20) {
        return Result<QuadratureEstimate>::failure(
            CasErrc::InvalidArgument,
            "Gaussian quadrature order must be in [1, 20]",
            operation);
    }
    auto interval = evaluate_interval(lower, upper, context, operation);
    if (!interval) return Result<QuadratureEstimate>::failure(interval.error());
    if (interval.value().lower == interval.value().upper) {
        return Result<QuadratureEstimate>::success(zero_quadrature_result());
    }
    const int comparison_order = order == 20 ? 19 : order + 1;
    auto budget = context.consume_steps(
        static_cast<std::size_t>(order + comparison_order + 2), operation);
    if (!budget) return Result<QuadratureEstimate>::failure(budget.error());

    detail::ensure_lmmc_lifecycle();
    QuadratureCallback callback{
        &function, &variable, &context, operation, std::nullopt};
    double value = 0.0;
    auto status = lmmc_quad_gauss_legendre(
        evaluate_quadrature_callback, &callback,
        interval.value().lower, interval.value().upper,
        static_cast<std::size_t>(order), &value);
    if (callback.error) return Result<QuadratureEstimate>::failure(*callback.error);
    if (status != LMMC_STATUS_OK) {
        return quadrature_failure<QuadratureEstimate>(status, operation);
    }

    double comparison_value = 0.0;
    status = lmmc_quad_gauss_legendre(
        evaluate_quadrature_callback, &callback,
        interval.value().lower, interval.value().upper,
        static_cast<std::size_t>(comparison_order), &comparison_value);
    if (callback.error) return Result<QuadratureEstimate>::failure(*callback.error);
    if (status != LMMC_STATUS_OK) {
        return quadrature_failure<QuadratureEstimate>(status, operation);
    }
    return estimated_quadrature_result(
        value, std::abs(comparison_value - value),
        interval.value().sign, operation);
}

Result<QuadratureEstimate> quadrature_gaussian_numeric(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    int order) {
    ComputationContext context;
    return quadrature_gaussian_numeric(
        function, variable, lower, upper, context, order);
}

Result<QuadratureEstimate> adaptive_simpson_numeric(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    ComputationContext& context,
    double tolerance,
    int max_depth) {
    constexpr const char* operation = "adaptive_simpson";
    auto valid = validate_quadrature_args(
        function, variable, lower, upper, operation);
    if (!valid) return Result<QuadratureEstimate>::failure(valid.error());
    if (!(tolerance > 0.0) || !std::isfinite(tolerance)) {
        return Result<QuadratureEstimate>::failure(
            CasErrc::InvalidArgument,
            "tolerance must be finite and positive",
            operation);
    }
    if (max_depth < 0) {
        return Result<QuadratureEstimate>::failure(
            CasErrc::InvalidArgument,
            "maximum recursion depth must be non-negative",
            operation);
    }
    auto interval = evaluate_interval(lower, upper, context, operation);
    if (!interval) return Result<QuadratureEstimate>::failure(interval.error());
    if (interval.value().lower == interval.value().upper) {
        return Result<QuadratureEstimate>::success(zero_quadrature_result());
    }
    auto recursion = context.enter_recursion(operation);
    if (!recursion) return Result<QuadratureEstimate>::failure(recursion.error());
    struct RecursionGuard {
        ComputationContext& context;
        ~RecursionGuard() { context.leave_recursion(); }
    } recursion_guard{context};
    const std::size_t remaining_depth =
        context.limits().max_recursion_depth > context.recursion_depth()
            ? context.limits().max_recursion_depth - context.recursion_depth()
            : 0;
    const std::size_t effective_depth = std::min(
        static_cast<std::size_t>(max_depth), remaining_depth);

    detail::ensure_lmmc_lifecycle();
    QuadratureCallback callback{
        &function, &variable, &context, operation, std::nullopt};
    lmmc_quad_result_t backend_result{};
    const auto status = lmmc_quad_adaptive(
        evaluate_quadrature_callback, &callback,
        interval.value().lower, interval.value().upper,
        tolerance, 0.0, effective_depth, &backend_result);
    if (callback.error) return Result<QuadratureEstimate>::failure(*callback.error);
    if (status == LMMC_STATUS_WARNING_MAX_DEPTH) {
        return Result<QuadratureEstimate>::failure(
            CasErrc::ResourceLimit,
            "adaptive Simpson recursion limit reached before meeting tolerance",
            operation);
    }
    if (status != LMMC_STATUS_OK) {
        return quadrature_failure<QuadratureEstimate>(status, operation);
    }
    QuadratureEstimate result;
    result.value = interval.value().sign * backend_result.value;
    result.estimated_error = std::abs(backend_result.error) +
        std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::abs(result.value)) * 16.0;
    if (!std::isfinite(result.value) || !std::isfinite(result.estimated_error)) {
        return Result<QuadratureEstimate>::failure(
            CasErrc::NumericFailure,
            "adaptive Simpson result is not finite",
            operation);
    }
    return Result<QuadratureEstimate>::success(result);
}

Result<QuadratureEstimate> adaptive_simpson_numeric(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& lower,
    const std::shared_ptr<SymbolicExpr>& upper,
    double tolerance,
    int max_depth) {
    ComputationContext context;
    return adaptive_simpson_numeric(
        function, variable, lower, upper, context, tolerance, max_depth);
}

}
