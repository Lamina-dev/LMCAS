#include "newton_raphson.hpp"

#include "internal/exact_sturm.hpp"
#include "internal/lmmc_lifecycle.hpp"
#include "lmmc/nonlinear.h"
#include "numeric_evaluation.hpp"
#include "poly_utils.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <optional>

namespace LMCAS {
namespace {

constexpr const char* kBisectionOperation = "bisection";
constexpr const char* kNewtonOperation = "newton_raphson";

double finite_midpoint(double lower, double upper) {
    if (lower < 0.0 && upper > 0.0) {
        return lower * 0.5 + upper * 0.5;
    }
    return lower + (upper - lower) * 0.5;
}

Result<double> evaluate_root_function(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable,
    lmmc_real_t value,
    ComputationContext& context,
    const char* operation) {
    if (!expression) {
        return Result<double>::failure(
            CasErrc::InvalidArgument,
            "root function expression cannot be null",
            operation);
    }
    if (variable.empty()) {
        return Result<double>::failure(
            CasErrc::InvalidArgument,
            "root variable cannot be empty",
            operation);
    }
    auto evaluated = evaluate_numeric(
        *expression, NumericBindings{{variable, static_cast<double>(value)}},
        context);
    if (!evaluated) return Result<double>::failure(evaluated.error());
    if (!evaluated.value().is_finite() ||
        !std::isfinite(evaluated.value().value)) {
        return Result<double>::failure(
            CasErrc::NumericFailure,
            "root function evaluation is not finite",
            operation);
    }
    return Result<double>::success(evaluated.value().value);
}

NumericRootResult invalid_root_options(
    const SolveOptions& options,
    const char* operation) {
    if (!std::isfinite(options.tolerance) || options.tolerance <= 0.0) {
        return NumericRootResult::failure(
            CasErrc::InvalidArgument,
            "root tolerance must be finite and positive",
            operation);
    }
    if (options.max_newton_iterations <= 0) {
        return NumericRootResult::failure(
            CasErrc::InvalidArgument,
            "maximum iteration count must be positive",
            operation);
    }
    return NumericRootResult::success(std::nullopt);
}

struct RootCallbacks {
    const std::shared_ptr<SymbolicExpr>* function = nullptr;
    const std::shared_ptr<SymbolicExpr>* derivative = nullptr;
    const std::string* variable = nullptr;
    ComputationContext* context = nullptr;
    const char* operation = nullptr;
    bool enforce_bracket = false;
    double bracket_lower = 0.0;
    double bracket_upper = 0.0;
    bool left_bracket = false;
    std::optional<CasError> error;
};

double evaluate_root_callback(
    const std::shared_ptr<SymbolicExpr>& expression,
    double x,
    RootCallbacks& callbacks) {
    if (callbacks.error) return std::numeric_limits<double>::quiet_NaN();
    if (callbacks.enforce_bracket &&
        (x < callbacks.bracket_lower || x > callbacks.bracket_upper)) {
        callbacks.left_bracket = true;
        return std::numeric_limits<double>::quiet_NaN();
    }
    auto value = evaluate_root_function(
        expression, *callbacks.variable, x, *callbacks.context,
        callbacks.operation);
    if (!value) {
        callbacks.error = value.error();
        return std::numeric_limits<double>::quiet_NaN();
    }
    return value.value();
}

double root_function_callback(double x, void* user_data) {
    auto& callbacks = *static_cast<RootCallbacks*>(user_data);
    return evaluate_root_callback(*callbacks.function, x, callbacks);
}

double root_derivative_callback(double x, void* user_data) {
    auto& callbacks = *static_cast<RootCallbacks*>(user_data);
    return evaluate_root_callback(*callbacks.derivative, x, callbacks);
}

lmmc_nonlinear_config_t root_config(
    const SolveOptions& options,
    std::size_t max_iterations) {
    lmmc_nonlinear_config_t config{};
    if (lmmc_nonlinear_default_config(&config) != LMMC_STATUS_OK) {
        std::terminate();
    }
    config.abs_tol = options.tolerance;
    config.rel_tol = 0.0;
    config.max_iter = max_iterations;
    config.min_derivative = 1.0e-15;
    config.min_step = 0.0;
    return config;
}

NumericRootResult backend_root_result(
    const lmmc_nonlinear_result_t& backend,
    double residual_limit) {
    if (!backend.converged || !std::isfinite(backend.root) ||
        !std::isfinite(backend.residual_norm) ||
        backend.residual_norm > residual_limit) {
        return NumericRootResult::success(std::nullopt);
    }
    const int iterations =
        backend.num_iter > static_cast<std::size_t>(std::numeric_limits<int>::max())
            ? std::numeric_limits<int>::max()
            : static_cast<int>(backend.num_iter);
    return NumericRootResult::success(
        NumericRoot{backend.root, backend.residual_norm, iterations});
}

NumericRootResult run_bisection_backend(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::string& variable,
    double lower,
    double upper,
    ComputationContext& context,
    const SolveOptions& options,
    const char* operation) {
    if (lower == upper) {
        auto value = evaluate_root_function(
            function, variable, lower, context, operation);
        if (!value) return NumericRootResult::failure(value.error());
        if (std::abs(value.value()) <= options.tolerance) {
            return NumericRootResult::success(
                NumericRoot{lower, std::abs(value.value()), 0});
        }
        return NumericRootResult::success(std::nullopt);
    }

    detail::ensure_lmmc_lifecycle();
    RootCallbacks callbacks{
        &function, nullptr, &variable, &context, operation,
        false, lower, upper, false, std::nullopt};
    const std::size_t max_iterations =
        static_cast<std::size_t>(options.max_newton_iterations) * 3;
    auto config = root_config(options, max_iterations);
    lmmc_nonlinear_result_t backend{};
    const auto status = lmmc_bisection_solve(
        root_function_callback, &callbacks, lower, upper, &config, &backend);
    if (callbacks.error) return NumericRootResult::failure(*callbacks.error);
    if (status == LMMC_STATUS_INVALID_ARGUMENT &&
        backend.failure_reason == LMMC_NONLINEAR_FAILURE_INVALID_BRACKET) {
        return NumericRootResult::success(std::nullopt);
    }
    return backend_root_result(backend, options.tolerance * 100.0);
}

NumericRootResult run_newton_backend(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& variable,
    double initial,
    ComputationContext& context,
    const SolveOptions& options,
    bool enforce_bracket,
    double bracket_lower,
    double bracket_upper) {
    detail::ensure_lmmc_lifecycle();
    RootCallbacks callbacks{
        &function, &derivative, &variable, &context, kNewtonOperation,
        enforce_bracket, bracket_lower, bracket_upper, false, std::nullopt};
    auto config = root_config(
        options, static_cast<std::size_t>(options.max_newton_iterations));
    lmmc_nonlinear_result_t backend{};
    (void)lmmc_newton_solve(
        root_function_callback, root_derivative_callback, &callbacks, initial,
        &config, &backend);
    if (callbacks.error) return NumericRootResult::failure(*callbacks.error);
    return backend_root_result(backend, options.tolerance);
}

} // namespace

Result<std::vector<std::pair<Rational, Rational>>>
isolate_real_roots_checked(
    const Polynomial<Rational>& polynomial,
    ComputationContext& context) {
    return detail::isolate_real_roots_exact(
        polynomial, context, "isolate_real_roots");
}

Result<std::vector<std::pair<Rational, Rational>>>
isolate_real_roots_checked(const Polynomial<Rational>& polynomial) {
    ComputationContext context;
    return isolate_real_roots_checked(polynomial, context);
}

NumericRootResult bisection_checked(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::string& variable,
    lmmc_real_t lower,
    lmmc_real_t upper,
    ComputationContext& context,
    const SolveOptions& options) {
    auto valid = invalid_root_options(options, kBisectionOperation);
    if (!valid) return valid;
    if (!function || variable.empty() || !std::isfinite(lower) ||
        !std::isfinite(upper) || lower > upper) {
        return NumericRootResult::failure(
            CasErrc::InvalidArgument,
            "bisection function, variable, and ordered finite interval are required",
            kBisectionOperation);
    }
    return run_bisection_backend(
        function, variable, lower, upper, context, options,
        kBisectionOperation);
}

NumericRootResult bisection_checked(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::string& variable,
    lmmc_real_t lower,
    lmmc_real_t upper,
    const SolveOptions& options) {
    ComputationContext context;
    return bisection_checked(
        function, variable, lower, upper, context, options);
}

NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& variable,
    lmmc_real_t initial,
    lmmc_real_t bracket_lower,
    lmmc_real_t bracket_upper,
    ComputationContext& context,
    const SolveOptions& options) {
    auto valid = invalid_root_options(options, kNewtonOperation);
    if (!valid) return valid;
    if (!function || !derivative || variable.empty() ||
        !std::isfinite(initial) || !std::isfinite(bracket_lower) ||
        !std::isfinite(bracket_upper) || bracket_lower > bracket_upper ||
        initial < bracket_lower || initial > bracket_upper) {
        return NumericRootResult::failure(
            CasErrc::InvalidArgument,
            "Newton function, derivative, variable, bracket, and initial value are invalid",
            kNewtonOperation);
    }
    auto newton = run_newton_backend(
        function, derivative, variable, initial, context, options, true,
        bracket_lower, bracket_upper);
    if (!newton || newton.value()) return newton;
    return run_bisection_backend(
        function, variable, bracket_lower, bracket_upper, context, options,
        kBisectionOperation);
}

NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& variable,
    lmmc_real_t initial,
    lmmc_real_t bracket_lower,
    lmmc_real_t bracket_upper,
    const SolveOptions& options) {
    ComputationContext context;
    return newton_raphson_checked(
        function, derivative, variable, initial, bracket_lower, bracket_upper,
        context, options);
}

NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& variable,
    lmmc_real_t initial,
    ComputationContext& context,
    const SolveOptions& options) {
    auto valid = invalid_root_options(options, kNewtonOperation);
    if (!valid) return valid;
    if (!function || !derivative || variable.empty() ||
        !std::isfinite(initial)) {
        return NumericRootResult::failure(
            CasErrc::InvalidArgument,
            "Newton function, derivative, variable, and finite initial value are required",
            kNewtonOperation);
    }
    return run_newton_backend(
        function, derivative, variable, initial, context, options, false, 0.0,
        0.0);
}

NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& function,
    const std::shared_ptr<SymbolicExpr>& derivative,
    const std::string& variable,
    lmmc_real_t initial,
    const SolveOptions& options) {
    ComputationContext context;
    return newton_raphson_checked(
        function, derivative, variable, initial, context, options);
}


NumericRootsResult solve_numeric_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts)
{
    ComputationContext context;
    return solve_numeric_checked(expr, var, context, opts);
}

NumericRootsResult solve_numeric_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    ComputationContext& context,
    const SolveOptions& opts)
{
    constexpr const char* operation = "solve_numeric";
    std::vector<NumericRoot> results;
    auto options = invalid_root_options(opts, operation);
    if (!options) return NumericRootsResult::failure(options.error());
    if (!expr) {
        return NumericRootsResult::failure(CasErrc::InvalidArgument,
                                           "expression cannot be null",
                                           operation);
    }
    if (var.empty()) {
        return NumericRootsResult::failure(CasErrc::InvalidArgument,
                                           "solve variable cannot be empty",
                                           operation);
    }
    if (opts.max_roots < -1) {
        return NumericRootsResult::failure(CasErrc::InvalidArgument,
                                           "maximum root count must be -1 or non-negative",
                                           operation);
    }
    if (opts.max_roots == 0) {
        return NumericRootsResult::success(std::move(results));
    }

    auto initial_step = context.consume_steps(1, operation);
    if (!initial_step) return NumericRootsResult::failure(initial_step.error());

    auto recognized_poly = recognize_rational_polynomial(*expr, var, context);
    if (!recognized_poly) return NumericRootsResult::failure(recognized_poly.error());

    if (recognized_poly.value() && !recognized_poly.value()->is_zero() &&
        recognized_poly.value()->degree() >= 1) {
        const Polynomial<Rational>& poly = *recognized_poly.value();
        auto df_expr = expr->differentiate(var);
        auto isolated = isolate_real_roots_checked(poly, context);
        if (!isolated) return NumericRootsResult::failure(isolated.error());
        for (const auto& [lo_rat, hi_rat] : isolated.value()) {
            if (opts.max_roots > 0 &&
                static_cast<int>(results.size()) >= opts.max_roots) {
                break;
            }

            auto interval_step = context.consume_steps(1, operation);
            if (!interval_step) return NumericRootsResult::failure(interval_step.error());

            lmmc_real_t lo = lo_rat.to_double();
            lmmc_real_t hi = hi_rat.to_double();
            lmmc_real_t x0 = finite_midpoint(lo, hi);
            auto root_result = newton_raphson_checked(
                expr, df_expr, var, x0, lo, hi, context, opts);
            if (!root_result) return NumericRootsResult::failure(root_result.error());
            if (!root_result.value()) continue;

            NumericRoot root = *root_result.value();
            auto verified = evaluate_root_function(
                expr, var, root.value, context, operation);
            if (!verified) return NumericRootsResult::failure(verified.error());
            root.residual = std::abs(verified.value());
            if (root.residual <= opts.tolerance * 100.0) {
                results.push_back(root);
            }
        }
    } else {
        lmmc_real_t x0 = opts.has_initial_guess ? opts.initial_guess : 0.0;
        auto df_expr = expr->differentiate(var);
        auto root_result = newton_raphson_checked(
            expr, df_expr, var, x0, context, opts);
        if (!root_result) return NumericRootsResult::failure(root_result.error());
        if (root_result.value()) {
            NumericRoot root = *root_result.value();
            auto verified = evaluate_root_function(
                expr, var, root.value, context, operation);
            if (!verified) return NumericRootsResult::failure(verified.error());
            root.residual = std::abs(verified.value());
            if (root.residual <= opts.tolerance) {
                results.push_back(root);
            }
        }
    }

    return NumericRootsResult::success(std::move(results));
}

}
