#include "internal/exact_root.hpp"
#include "internal/exact_sturm.hpp"


#include <cmath>
#include <algorithm>
#include <limits>

namespace LMCAS::detail {
namespace {

bool valid_options(const NumericEvaluationOptions& options) {
    return options.absolute_tolerance > 0.0 &&
           options.relative_tolerance > 0.0 &&
           std::isfinite(options.absolute_tolerance) &&
           std::isfinite(options.relative_tolerance);
}

double conversion_ulp(double value) {
    if (!std::isfinite(value)) {
        return std::numeric_limits<double>::infinity();
    }
    const double upward = std::nextafter(
        value, std::numeric_limits<double>::infinity());
    const double downward = std::nextafter(
        value, -std::numeric_limits<double>::infinity());
    double spacing = 0.0;
    if (std::isfinite(upward)) {
        spacing = std::max(spacing, std::abs(upward - value));
    }
    if (std::isfinite(downward)) {
        spacing = std::max(spacing, std::abs(value - downward));
    }
    return spacing;
}

Result<ApproxReal> approximate_real_isolation(
    ExactRealAlgebraic isolation,
    const NumericEvaluationOptions& options,
    ComputationContext& context,
    const std::string& operation) {
    while (isolation.lower != isolation.upper) {
        const bool positive = isolation.lower > Rational(0) &&
            isolation.upper <= Rational(2) * isolation.lower;
        const bool negative = isolation.upper < Rational(0) &&
            isolation.lower >= Rational(2) * isolation.upper;
        if (positive || negative) break;
        auto narrowed = refine_exact_real_algebraic(isolation, context, operation);
        if (!narrowed) return Result<ApproxReal>::failure(narrowed.error());
    }
    const double scale = std::min(
        std::abs(isolation.lower.to_double()), std::abs(isolation.upper.to_double()));
    const double absolute_tolerance = std::min(
        options.absolute_tolerance,
        std::max(std::numeric_limits<double>::denorm_min(),
                 options.relative_tolerance * scale));
    auto refined = refine_exact_real_algebraic_to_tolerance(
        isolation, absolute_tolerance, options.relative_tolerance,
        context, operation);
    if (!refined) return Result<ApproxReal>::failure(refined.error());
    const double lower = isolation.lower.to_double();
    const double upper = isolation.upper.to_double();
    const double value = lower + (upper - lower) * 0.5;
    const double error = std::abs(upper - lower) * 0.5 +
        conversion_ulp(lower) + conversion_ulp(upper) +
        conversion_ulp(value);
    if (!std::isfinite(value) || !std::isfinite(error) ||
        (value == 0.0 &&
         (isolation.lower > Rational(0) || isolation.upper < Rational(0)))) {
        return Result<ApproxReal>::failure(
            CasErrc::NumericFailure,
            "certified root enclosure could not be represented as finite double",
            operation);
    }
    return Result<ApproxReal>::success(
        ApproxReal{value, error, NumericStatus::Finite});
}

ApproxReal approximate_component(
    const Rational& lower,
    const Rational& upper) {
    const double lower_value = lower.to_double();
    const double upper_value = upper.to_double();
    const double value =
        lower_value + (upper_value - lower_value) * 0.5;
    const double error = std::abs(upper_value - lower_value) * 0.5 +
        conversion_ulp(lower_value) + conversion_ulp(upper_value) +
        conversion_ulp(value);
    return ApproxReal{value, error, NumericStatus::Finite};
}

} // namespace


Result<ExactRootId> make_exact_root_id(
    Polynomial<Rational> polynomial,
    std::size_t index,
    ComputationContext& context,
    const std::string& operation) {
    auto access = context.consume_steps(1, operation);
    if (!access) return Result<ExactRootId>::failure(access.error());
    try {
        polynomial.trim();
        if (polynomial.degree() <= 0) {
            return Result<ExactRootId>::failure(
                CasErrc::InvalidArgument,
                "RootOf requires a nonconstant exact rational polynomial",
                operation);
        }
        polynomial = polynomial.square_free_part().make_monic();
        polynomial.variable_name = "_root";
        if (polynomial.degree() <= 0 ||
            index >= static_cast<std::size_t>(polynomial.degree())) {
            return Result<ExactRootId>::failure(
                CasErrc::InvalidArgument,
                "RootOf index exceeds the number of distinct roots",
                operation);
        }
        return Result<ExactRootId>::success(
            ExactRootId{std::move(polynomial), index});
    } catch (const std::bad_alloc&) {
        return Result<ExactRootId>::failure(
            CasErrc::ResourceLimit,
            "RootOf canonicalization allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<ExactRootId>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

Result<ApproxReal> evaluate_root_real(
    const ExactRootId& root,
    const NumericEvaluationOptions& options,
    ComputationContext& context) {
    constexpr const char* operation = "rootof.evaluate_real";
    if (!valid_options(options)) {
        return Result<ApproxReal>::failure(
            CasErrc::InvalidArgument,
            "root evaluation tolerances must be finite and positive",
            operation);
    }
    if (root.index >= static_cast<std::size_t>(root.polynomial.degree())) {
        return Result<ApproxReal>::failure(
            CasErrc::InvalidArgument,
            "RootOf index exceeds the number of distinct roots",
            operation);
    }
    auto real_roots = isolate_real_roots_exact(root.polynomial, context, operation);
    if (!real_roots) return Result<ApproxReal>::failure(real_roots.error());
    if (root.index >= real_roots.value().size()) {
        return Result<ApproxReal>::failure(
            CasErrc::DomainError,
            "selected RootOf is not real",
            operation);
    }
    const auto& interval = real_roots.value()[root.index];
    return approximate_real_isolation(
        ExactRealAlgebraic{root.polynomial, interval.first, interval.second, root.index, 1},
        options, context, operation);
}

Result<ApproxComplex> evaluate_root_complex(
    const ExactRootId& root,
    const NumericEvaluationOptions& options,
    ComputationContext& context) {
    constexpr const char* operation = "rootof.evaluate_complex";
    if (!valid_options(options)) {
        return Result<ApproxComplex>::failure(
            CasErrc::InvalidArgument,
            "root evaluation tolerances must be finite and positive",
            operation);
    }
    auto isolated = isolate_exact_root(root, context, operation);
    if (!isolated) return Result<ApproxComplex>::failure(isolated.error());
    if (auto* real = std::get_if<RealIsolation>(&isolated.value())) {
        auto approximation = approximate_real_isolation(
            std::move(real->value), options, context, operation);
        if (!approximation) {
            return Result<ApproxComplex>::failure(approximation.error());
        }
        return Result<ApproxComplex>::success(ApproxComplex{
            approximation.value(),
            ApproxReal{0.0, 0.0, NumericStatus::Finite}});
    }
    auto complex = std::get<ComplexIsolation>(std::move(isolated.value()));
    auto refined = refine_complex_isolation(
        root.polynomial, complex, options, context, operation);
    if (!refined) return Result<ApproxComplex>::failure(refined.error());
    auto real_part = approximate_component(
        complex.real_lower, complex.real_upper);
    auto imaginary_part = approximate_component(
        complex.imaginary_lower, complex.imaginary_upper);
    if (!std::isfinite(real_part.value) ||
        !std::isfinite(real_part.absolute_error) ||
        !std::isfinite(imaginary_part.value) ||
        !std::isfinite(imaginary_part.absolute_error)) {
        return Result<ApproxComplex>::failure(
            CasErrc::NumericFailure,
            "certified complex root enclosure could not be represented as finite double",
            operation);
    }
    return Result<ApproxComplex>::success(
        ApproxComplex{real_part, imaginary_part});
}

} // namespace LMCAS::detail
