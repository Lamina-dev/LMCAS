/**
 * @file numerical_integration.cpp
 * @brief 数值积分实现。
 */
#include "numerical_integration.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>

namespace lamina {


namespace {

Result<double> checked_finite_value(Result<ApproxReal> value, const std::string& operation) {
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
        return Result<double>::failure(CasErrc::InvalidArgument,
                                       "integration bound is null", operation);
    }

    return checked_finite_value(evaluate_numeric(*expr, NumericBindings{}, context),
                                operation);
}

/// 在数值点处求 f 的 double 值；失败返回 CasError。
Result<double> eval_f_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    double xval,
    ComputationContext& context) {
    if (!f) {
        return Result<double>::failure(CasErrc::InvalidArgument,
                                       "integrand is null",
                                       "adaptive_simpson");
    }
    NumericBindings bindings;
    bindings.emplace(var, xval);
    return checked_finite_value(evaluate_numeric(*f, bindings, context),
                                "adaptive_simpson");
}

Result<double> eval_f_checked_for(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    double xval,
    ComputationContext& context,
    const std::string& operation) {
    if (!f) {
        return Result<double>::failure(CasErrc::InvalidArgument,
                                       "integrand is null",
                                       operation);
    }
    NumericBindings bindings;
    bindings.emplace(var, xval);
    return checked_finite_value(evaluate_numeric(*f, bindings, context),
                                operation);
}

Result<void> validate_fixed_quadrature_args(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n,
    const std::string& operation) {
    if (!f || !a || !b) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "integrand and bounds must be non-null",
                                     operation);
    }
    if (var.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "integration variable must not be empty",
                                     operation);
    }
    if (n <= 0) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "quadrature sample count must be positive",
                                     operation);
    }
    return Result<void>::success();
}

Result<ApproxReal> make_finite_quadrature_result(
    double value,
    int samples,
    const std::string& operation) {
    if (!std::isfinite(value)) {
        return Result<ApproxReal>::failure(CasErrc::NumericFailure,
                                           "quadrature accumulation produced a non-finite value",
                                           operation);
    }
    ApproxReal result;
    result.value = value;
    result.absolute_error = std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(value)) *
        static_cast<double>(std::max(samples, 1) * 16);
    result.status = NumericStatus::Finite;
    return Result<ApproxReal>::success(result);
}

/// 单段辛普森
double simpson_seg(double a, double b, double fa, double fb, double fm) {
    double h = b - a;
    return (h / 6.0) * (fa + 4.0 * fm + fb);
}

struct IntegrationEstimate {
    double value = 0.0;
    double absolute_error = 0.0;
};

Result<IntegrationEstimate> adaptive_rec_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    double a,
    double b,
    double fa,
    double fb,
    double fm,
    double whole,
    double tol,
    int depth,
    ComputationContext& context) {
    auto recursion = context.enter_recursion("adaptive_simpson");
    if (!recursion) return Result<IntegrationEstimate>::failure(recursion.error());
    struct RecursionGuard {
        ComputationContext& context;
        ~RecursionGuard() { context.leave_recursion(); }
    } guard{context};

    double m = 0.5 * a + 0.5 * b;
    double lm = 0.5 * a + 0.5 * m;
    double rm = 0.5 * m + 0.5 * b;
    auto flm = eval_f_checked(f, var, lm, context);
    if (!flm) return Result<IntegrationEstimate>::failure(flm.error());
    auto frm = eval_f_checked(f, var, rm, context);
    if (!frm) return Result<IntegrationEstimate>::failure(frm.error());
    double left = simpson_seg(a, m, fa, fm, flm.value());
    double right = simpson_seg(m, b, fm, fb, frm.value());
    if (!std::isfinite(left) || !std::isfinite(right)) {
        return Result<IntegrationEstimate>::failure(
            CasErrc::NumericFailure,
            "adaptive Simpson segment overflowed",
            "adaptive_simpson");
    }
    const double delta = left + right - whole;
    const double truncation_error = std::abs(delta) / 15.0;
    if (truncation_error <= tol) {
        const double corrected = left + right + delta / 15.0;
        if (!std::isfinite(corrected)) {
            return Result<IntegrationEstimate>::failure(
                CasErrc::NumericFailure,
                "adaptive Simpson correction produced a non-finite value",
                "adaptive_simpson");
        }
        const double rounding_error = std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::abs(corrected)) * 16.0;
        return Result<IntegrationEstimate>::success(
            IntegrationEstimate{corrected, truncation_error + rounding_error});
    }
    if (depth <= 0) {
        return Result<IntegrationEstimate>::failure(
            CasErrc::ResourceLimit,
            "adaptive Simpson recursion limit reached before meeting tolerance",
            "adaptive_simpson");
    }
    auto left_result = adaptive_rec_checked(
        f, var, a, m, fa, fm, flm.value(), left, tol / 2.0, depth - 1, context);
    if (!left_result) return Result<IntegrationEstimate>::failure(left_result.error());
    auto right_result = adaptive_rec_checked(
        f, var, m, b, fm, fb, frm.value(), right, tol / 2.0, depth - 1, context);
    if (!right_result) return Result<IntegrationEstimate>::failure(right_result.error());
    const double combined_value =
        left_result.value().value + right_result.value().value;
    const double combined_error =
        left_result.value().absolute_error + right_result.value().absolute_error;
    if (!std::isfinite(combined_value) || !std::isfinite(combined_error)) {
        return Result<IntegrationEstimate>::failure(
            CasErrc::NumericFailure,
            "adaptive Simpson accumulation produced a non-finite value",
            "adaptive_simpson");
    }
    return Result<IntegrationEstimate>::success(
        IntegrationEstimate{combined_value, combined_error});
}

} // anonymous namespace

Result<ApproxReal> quadrature_simpson_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    int n) {
    constexpr const char* operation = "quadrature_simpson";
    auto valid = validate_fixed_quadrature_args(f, var, a, b, n, operation);
    if (!valid) return Result<ApproxReal>::failure(valid.error());

    if (n % 2 != 0) {
        return Result<ApproxReal>::failure(
            CasErrc::InvalidArgument,
            "Simpson subinterval count must be even",
            operation);
    }
    const auto sample_count = static_cast<std::size_t>(n) + 1;
    auto budget = context.consume_steps(sample_count, operation);
    if (!budget) return Result<ApproxReal>::failure(budget.error());

    auto av = eval_bound_checked(a, context, operation);
    if (!av) return Result<ApproxReal>::failure(av.error());
    auto bv = eval_bound_checked(b, context, operation);
    if (!bv) return Result<ApproxReal>::failure(bv.error());
    const double width = bv.value() - av.value();
    if (!std::isfinite(width)) {
        return Result<ApproxReal>::failure(
            CasErrc::NumericFailure,
            "integration interval width is not representable as a finite double",
            operation);
    }
    const double h = width / static_cast<double>(n);
    if (!std::isfinite(h)) {
        return Result<ApproxReal>::failure(
            CasErrc::NumericFailure,
            "Simpson step size is not finite",
            operation);
    }

    auto fa = eval_f_checked_for(f, var, av.value(), context, operation);
    if (!fa) return Result<ApproxReal>::failure(fa.error());
    auto fb = eval_f_checked_for(f, var, bv.value(), context, operation);
    if (!fb) return Result<ApproxReal>::failure(fb.error());
    double sum = fa.value() + fb.value();
    if (!std::isfinite(sum)) {
        return Result<ApproxReal>::failure(CasErrc::NumericFailure,
                                           "Simpson endpoint accumulation is not finite",
                                           operation);
    }

    for (int i = 1; i < n; ++i) {
        const double xi = av.value() + h * static_cast<double>(i);
        if (!std::isfinite(xi)) {
            return Result<ApproxReal>::failure(CasErrc::NumericFailure,
                                               "Simpson sample point is not finite",
                                               operation);
        }
        auto fxi = eval_f_checked_for(f, var, xi, context, operation);
        if (!fxi) return Result<ApproxReal>::failure(fxi.error());
        const double coeff = (i % 2 == 1) ? 4.0 : 2.0;
        sum += coeff * fxi.value();
        if (!std::isfinite(sum)) {
            return Result<ApproxReal>::failure(CasErrc::NumericFailure,
                                               "Simpson accumulation is not finite",
                                               operation);
        }
    }

    const double value = (h / 3.0) * sum;
    return make_finite_quadrature_result(value, n + 1, operation);
}

Result<ApproxReal> quadrature_simpson_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n) {
    ComputationContext context;
    return quadrature_simpson_numeric(f, var, a, b, context, n);
}

Result<ApproxReal> quadrature_gaussian_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    int n) {
    constexpr const char* operation = "quadrature_gaussian";
    auto valid = validate_fixed_quadrature_args(f, var, a, b, n, operation);
    if (!valid) return Result<ApproxReal>::failure(valid.error());
    if (n > 3) {
        if (n > std::numeric_limits<int>::max() / 2) {
            return Result<ApproxReal>::failure(CasErrc::ResourceLimit,
                                               "Gaussian fallback sample count overflows",
                                               operation);
        }
        return quadrature_simpson_numeric(f, var, a, b, context, n * 2);
    }

    auto budget = context.consume_steps(static_cast<std::size_t>(n) + 2, operation);
    if (!budget) return Result<ApproxReal>::failure(budget.error());
    auto av = eval_bound_checked(a, context, operation);
    if (!av) return Result<ApproxReal>::failure(av.error());
    auto bv = eval_bound_checked(b, context, operation);
    if (!bv) return Result<ApproxReal>::failure(bv.error());
    const double width = bv.value() - av.value();
    const double half_diff = 0.5 * width;
    const double half_sum = 0.5 * (av.value() + bv.value());
    if (!std::isfinite(width) || !std::isfinite(half_diff) ||
        !std::isfinite(half_sum)) {
        return Result<ApproxReal>::failure(
            CasErrc::NumericFailure,
            "Gaussian interval transform is not finite",
            operation);
    }

    std::array<double, 3> roots{};
    std::array<double, 3> weights{};
    int count = n;
    if (n <= 1) {
        count = 1;
        roots = {0.0, 0.0, 0.0};
        weights = {2.0, 0.0, 0.0};
    } else if (n == 2) {
        const double r = 1.0 / std::sqrt(3.0);
        roots = {-r, r, 0.0};
        weights = {1.0, 1.0, 0.0};
    } else {
        const double r = std::sqrt(3.0 / 5.0);
        roots = {-r, 0.0, r};
        weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    }

    double sum = 0.0;
    for (int i = 0; i < count; ++i) {
        const double x = half_diff * roots[static_cast<std::size_t>(i)] + half_sum;
        if (!std::isfinite(x)) {
            return Result<ApproxReal>::failure(CasErrc::NumericFailure,
                                               "Gaussian sample point is not finite",
                                               operation);
        }
        auto fx = eval_f_checked_for(f, var, x, context, operation);
        if (!fx) return Result<ApproxReal>::failure(fx.error());
        sum += weights[static_cast<std::size_t>(i)] * fx.value();
        if (!std::isfinite(sum)) {
            return Result<ApproxReal>::failure(CasErrc::NumericFailure,
                                               "Gaussian accumulation is not finite",
                                               operation);
        }
    }

    return make_finite_quadrature_result(half_diff * sum, count, operation);
}

Result<ApproxReal> quadrature_gaussian_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n) {
    ComputationContext context;
    return quadrature_gaussian_numeric(f, var, a, b, context, n);
}

Result<ApproxReal> adaptive_simpson_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    double tol,
    int max_depth) {
    if (!f || !a || !b) {
        return Result<ApproxReal>::failure(CasErrc::InvalidArgument,
                                           "integrand and bounds must be non-null",
                                           "adaptive_simpson");
    }
    if (var.empty()) {
        return Result<ApproxReal>::failure(CasErrc::InvalidArgument,
                                           "integration variable must not be empty",
                                           "adaptive_simpson");
    }
    if (!(tol > 0.0) || !std::isfinite(tol)) {
        return Result<ApproxReal>::failure(CasErrc::InvalidArgument,
                                           "tolerance must be finite and positive",
                                           "adaptive_simpson");
    }
    if (max_depth < 0) {
        return Result<ApproxReal>::failure(CasErrc::InvalidArgument,
                                           "maximum recursion depth must be non-negative",
                                           "adaptive_simpson");
    }

    auto av = eval_bound_checked(a, context, "adaptive_simpson");
    if (!av) return Result<ApproxReal>::failure(av.error());
    auto bv = eval_bound_checked(b, context, "adaptive_simpson");
    if (!bv) return Result<ApproxReal>::failure(bv.error());
    if (!std::isfinite(bv.value() - av.value())) {
        return Result<ApproxReal>::failure(
            CasErrc::NumericFailure,
            "integration interval width is not representable as a finite double",
            "adaptive_simpson");
    }

    auto fa = eval_f_checked(f, var, av.value(), context);
    if (!fa) return Result<ApproxReal>::failure(fa.error());
    auto fb = eval_f_checked(f, var, bv.value(), context);
    if (!fb) return Result<ApproxReal>::failure(fb.error());
    double m = 0.5 * av.value() + 0.5 * bv.value();
    auto fm = eval_f_checked(f, var, m, context);
    if (!fm) return Result<ApproxReal>::failure(fm.error());

    double whole = simpson_seg(av.value(), bv.value(),
                               fa.value(), fb.value(), fm.value());
    if (!std::isfinite(whole)) {
        return Result<ApproxReal>::failure(
            CasErrc::NumericFailure,
            "initial Simpson estimate is not finite",
            "adaptive_simpson");
    }
    auto numeric_result = adaptive_rec_checked(
        f, var, av.value(), bv.value(), fa.value(), fb.value(), fm.value(),
        whole, tol, max_depth, context);
    if (!numeric_result) return Result<ApproxReal>::failure(numeric_result.error());

    ApproxReal approx;
    approx.value = numeric_result.value().value;
    approx.absolute_error = numeric_result.value().absolute_error;
    approx.status = NumericStatus::Finite;
    return Result<ApproxReal>::success(approx);
}

Result<ApproxReal> adaptive_simpson_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    double tol,
    int max_depth) {
    ComputationContext context;
    return adaptive_simpson_numeric(f, var, a, b, context, tol, max_depth);
}


Result<ApproxReal> numerical_integrate_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    int n) {
    return quadrature_simpson_numeric(f, var, a, b, context, n);
}

Result<ApproxReal> numerical_integrate_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n) {
    ComputationContext context;
    return numerical_integrate_numeric(f, var, a, b, context, n);
}

} // namespace lamina
