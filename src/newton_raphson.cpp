#include "newton_raphson.hpp"
#include "internal/exact_sturm.hpp"
#include "numeric_evaluation.hpp"
#include "poly_utils.hpp"
#include <algorithm>
#include <cmath>

namespace lamina {

namespace {

constexpr const char* kBisectionOperation = "bisection";
constexpr const char* kNewtonOperation = "newton_raphson";

Result<double> evaluate_root_function(const std::shared_ptr<SymbolicExpr>& expression,
                                      const std::string& variable,
                                      lmmc_real_t value,
                                      ComputationContext& context,
                                      const char* operation) {
    if (!expression) {
        return Result<double>::failure(CasErrc::InvalidArgument,
                                       "root function expression cannot be null",
                                       operation);
    }
    if (variable.empty()) {
        return Result<double>::failure(CasErrc::InvalidArgument,
                                       "root variable cannot be empty",
                                       operation);
    }

    auto evaluated = evaluate_numeric(
        *expression, NumericBindings{{variable, static_cast<double>(value)}}, context);
    if (!evaluated) return Result<double>::failure(evaluated.error());
    if (!evaluated.value().is_finite()) {
        return Result<double>::failure(CasErrc::NumericFailure,
                                       "root function evaluation is not finite",
                                       operation);
    }
    return Result<double>::success(evaluated.value().value);
}

NumericRootResult invalid_root_options(const SolveOptions& opts,
                                       const char* operation) {
    if (!std::isfinite(opts.tolerance) || opts.tolerance <= 0.0) {
        return NumericRootResult::failure(CasErrc::InvalidArgument,
                                          "root tolerance must be finite and positive",
                                          operation);
    }
    if (opts.max_newton_iterations <= 0) {
        return NumericRootResult::failure(CasErrc::InvalidArgument,
                                          "maximum iteration count must be positive",
                                          operation);
    }
    return NumericRootResult::success(std::nullopt);
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
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    lmmc_real_t lo,
    lmmc_real_t hi,
    ComputationContext& context,
    const SolveOptions& opts)
{
    auto options = invalid_root_options(opts, kBisectionOperation);
    if (!options) return options;
    if (!std::isfinite(lo) || !std::isfinite(hi) || lo > hi) {
        return NumericRootResult::failure(CasErrc::InvalidArgument,
                                          "bisection interval must be finite and ordered",
                                          kBisectionOperation);
    }

    auto f_lo_result = evaluate_root_function(f, var, lo, context, kBisectionOperation);
    if (!f_lo_result) return NumericRootResult::failure(f_lo_result.error());
    auto f_hi_result = evaluate_root_function(f, var, hi, context, kBisectionOperation);
    if (!f_hi_result) return NumericRootResult::failure(f_hi_result.error());
    lmmc_real_t f_lo = f_lo_result.value();
    lmmc_real_t f_hi = f_hi_result.value();

    if (std::abs(f_lo) < opts.tolerance) {
        return NumericRootResult::success(NumericRoot{lo, std::abs(f_lo), 0});
    }
    if (std::abs(f_hi) < opts.tolerance) {
        return NumericRootResult::success(NumericRoot{hi, std::abs(f_hi), 0});
    }

    if (f_lo * f_hi > 0) {
        return NumericRootResult::success(std::nullopt);
    }

    lmmc_real_t best_x = std::abs(f_lo) <= std::abs(f_hi) ? lo : hi;
    lmmc_real_t best_residual = std::min(std::abs(f_lo), std::abs(f_hi));
    int best_iteration = 0;
    int max_iter = opts.max_newton_iterations * 3;
    for (int i = 1; i <= max_iter; ++i) {
        auto step = context.consume_steps(1, kBisectionOperation);
        if (!step) return NumericRootResult::failure(step.error());
        lmmc_real_t mid = (lo + hi) * 0.5;
        auto f_mid_result = evaluate_root_function(f, var, mid, context, kBisectionOperation);
        if (!f_mid_result) return NumericRootResult::failure(f_mid_result.error());
        lmmc_real_t f_mid = f_mid_result.value();
        lmmc_real_t residual = std::abs(f_mid);
        if (residual < best_residual) {
            best_x = mid;
            best_residual = residual;
            best_iteration = i;
        }

        if (residual < opts.tolerance) {
            return NumericRootResult::success(NumericRoot{mid, residual, i});
        }

        if (mid == lo || mid == hi) {
            break;
        }

        if (f_lo * f_mid < 0) {
            hi = mid;
            f_hi = f_mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }

    if (best_residual < opts.tolerance * 100.0) {
        return NumericRootResult::success(
            NumericRoot{best_x, best_residual, best_iteration});
    }
    return NumericRootResult::success(std::nullopt);
}


NumericRootResult bisection_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    lmmc_real_t lo,
    lmmc_real_t hi,
    const SolveOptions& opts)
{
    ComputationContext context;
    return bisection_checked(f, var, lo, hi, context, opts);
}

NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    lmmc_real_t bracket_lo,
    lmmc_real_t bracket_hi,
    ComputationContext& context,
    const SolveOptions& opts)
{
    auto options = invalid_root_options(opts, kNewtonOperation);
    if (!options) return options;
    if (!f || !df) {
        return NumericRootResult::failure(CasErrc::InvalidArgument,
                                          "function and derivative cannot be null",
                                          kNewtonOperation);
    }
    if (!std::isfinite(x0) || !std::isfinite(bracket_lo) ||
        !std::isfinite(bracket_hi) || bracket_lo > bracket_hi ||
        x0 < bracket_lo || x0 > bracket_hi) {
        return NumericRootResult::failure(CasErrc::InvalidArgument,
                                          "Newton bracket and initial value must be finite and ordered",
                                          kNewtonOperation);
    }
    lmmc_real_t x = x0;
    for (int i = 1; i <= opts.max_newton_iterations; ++i) {
        auto step = context.consume_steps(1, kNewtonOperation);
        if (!step) return NumericRootResult::failure(step.error());

        auto fx_result = evaluate_root_function(f, var, x, context, kNewtonOperation);
        if (!fx_result) return NumericRootResult::failure(fx_result.error());
        auto dfx_result = evaluate_root_function(df, var, x, context, kNewtonOperation);
        if (!dfx_result) return NumericRootResult::failure(dfx_result.error());
        lmmc_real_t fx = fx_result.value();
        lmmc_real_t dfx = dfx_result.value();

        if (std::abs(fx) < opts.tolerance) {
            return NumericRootResult::success(NumericRoot{x, std::abs(fx), i});
        }

        if (std::abs(dfx) < 1e-15) {
            return bisection_checked(f, var, bracket_lo, bracket_hi, context, opts);
        }

        lmmc_real_t x_new = x - fx / dfx;
        if (!std::isfinite(x_new)) {
            return NumericRootResult::failure(CasErrc::NumericFailure,
                                              "Newton update produced a non-finite value",
                                              kNewtonOperation);
        }

        if (i > 1 && std::abs(x_new - x) > 2.0 * std::abs(x - x0)) {
            x_new = x - 0.5 * fx / dfx;
        }

        if (x_new < bracket_lo || x_new > bracket_hi) {
            return bisection_checked(f, var, bracket_lo, bracket_hi, context, opts);
        }

        x = x_new;
    }

    return NumericRootResult::success(std::nullopt);
}


NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    lmmc_real_t bracket_lo,
    lmmc_real_t bracket_hi,
    const SolveOptions& opts)
{
    ComputationContext context;
    return newton_raphson_checked(
        f, df, var, x0, bracket_lo, bracket_hi, context, opts);
}

NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    ComputationContext& context,
    const SolveOptions& opts)
{
    auto options = invalid_root_options(opts, kNewtonOperation);
    if (!options) return options;
    if (!f || !df) {
        return NumericRootResult::failure(CasErrc::InvalidArgument,
                                          "function and derivative cannot be null",
                                          kNewtonOperation);
    }
    if (!std::isfinite(x0)) {
        return NumericRootResult::failure(CasErrc::InvalidArgument,
                                          "Newton initial value must be finite",
                                          kNewtonOperation);
    }
    lmmc_real_t x = x0;
    for (int i = 1; i <= opts.max_newton_iterations; ++i) {
        auto step = context.consume_steps(1, kNewtonOperation);
        if (!step) return NumericRootResult::failure(step.error());

        auto fx_result = evaluate_root_function(f, var, x, context, kNewtonOperation);
        if (!fx_result) return NumericRootResult::failure(fx_result.error());
        auto dfx_result = evaluate_root_function(df, var, x, context, kNewtonOperation);
        if (!dfx_result) return NumericRootResult::failure(dfx_result.error());
        lmmc_real_t fx = fx_result.value();
        lmmc_real_t dfx = dfx_result.value();

        if (std::abs(fx) < opts.tolerance) {
            return NumericRootResult::success(NumericRoot{x, std::abs(fx), i});
        }

        if (std::abs(dfx) < 1e-15) {
            return NumericRootResult::success(std::nullopt);
        }

        lmmc_real_t x_new = x - fx / dfx;
        if (!std::isfinite(x_new)) {
            return NumericRootResult::failure(CasErrc::NumericFailure,
                                              "Newton update produced a non-finite value",
                                              kNewtonOperation);
        }

        if (i > 1 && std::abs(x_new - x) > 2.0 * std::abs(x - x0)) {
            x_new = x - 0.5 * fx / dfx;
        }

        x = x_new;
    }

    return NumericRootResult::success(std::nullopt);
}


NumericRootResult newton_raphson_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    const SolveOptions& opts)
{
    ComputationContext context;
    return newton_raphson_checked(f, df, var, x0, context, opts);
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
            lmmc_real_t x0 = (lo + hi) * 0.5;
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
