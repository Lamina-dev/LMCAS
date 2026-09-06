#include "internal/exact_algebraic.hpp"

#include "internal/exact_sturm.hpp"
#include "solve_polynomial.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <utility>
#include <vector>

namespace LMCAS::detail {
namespace {

constexpr const char* kOperation = "exact_algebraic";

Polynomial<Rational> monic(Polynomial<Rational> polynomial) {
    polynomial.trim();
    if (polynomial.is_zero()) return polynomial;
    const Rational leading = polynomial.lead_coeff();
    for (auto& coefficient : polynomial.coeffs) coefficient = coefficient / leading;
    polynomial.trim();
    return polynomial;
}

bool same_polynomial(const Polynomial<Rational>& lhs,
                     const Polynomial<Rational>& rhs) {
    return lhs.variable_name == rhs.variable_name && lhs.coeffs == rhs.coeffs;
}

bool overlaps(const Rational& lhs_lower, const Rational& lhs_upper,
              const Rational& rhs_lower, const Rational& rhs_upper) {
    return !(lhs_upper < rhs_lower) && !(rhs_upper < lhs_lower);
}

int sign(const Rational& value) {
    if (value < Rational(0)) return -1;
    if (value > Rational(0)) return 1;
    return 0;
}

} // namespace

ExactRealAlgebraicResult make_exact_real_algebraic(
    Polynomial<Rational> polynomial,
    std::size_t root_index,
    std::size_t multiplicity,
    ComputationContext& context) {
    if (polynomial.degree() < 1) {
        return ExactRealAlgebraicResult::failure(
            CasErrc::InvalidArgument,
            "an algebraic value requires a non-constant polynomial",
            kOperation);
    }
    if (multiplicity == 0) {
        return ExactRealAlgebraicResult::failure(
            CasErrc::InvalidArgument,
            "algebraic root multiplicity must be positive",
            kOperation);
    }
    auto step = context.consume_steps(1, kOperation);
    if (!step) return ExactRealAlgebraicResult::failure(step.error());

    try {
        polynomial = monic(std::move(polynomial));
        auto rational_roots = find_rational_roots(polynomial);
        std::sort(rational_roots.begin(), rational_roots.end());
        rational_roots.erase(
            std::unique(rational_roots.begin(), rational_roots.end()),
            rational_roots.end());
        if (rational_roots.size() ==
            static_cast<std::size_t>(polynomial.degree())) {
            if (root_index >= rational_roots.size()) {
                return ExactRealAlgebraicResult::failure(
                    CasErrc::InvalidArgument,
                    "real algebraic root index is out of range",
                    kOperation);
            }
            return ExactRealAlgebraicResult::success(ExactRealAlgebraic{
                std::move(polynomial), rational_roots[root_index],
                rational_roots[root_index], root_index, multiplicity});
        }
        auto isolated = isolate_real_roots_exact(
            polynomial, context, kOperation);
        if (!isolated) {
            return ExactRealAlgebraicResult::failure(isolated.error());
        }
        const auto& intervals = isolated.value();
        if (root_index >= intervals.size()) {
            return ExactRealAlgebraicResult::failure(
                CasErrc::InvalidArgument,
                "real algebraic root index is out of range",
                kOperation);
        }
        return ExactRealAlgebraicResult::success(ExactRealAlgebraic{
            std::move(polynomial), intervals[root_index].first,
            intervals[root_index].second, root_index, multiplicity});
    } catch (const std::bad_alloc&) {
        return ExactRealAlgebraicResult::failure(
            CasErrc::ResourceLimit,
            "algebraic root isolation allocation failed",
            kOperation);
    } catch (const std::exception& error) {
        return ExactRealAlgebraicResult::failure(
            CasErrc::InternalInvariant, error.what(), kOperation);
    }
}

Result<void> refine_exact_real_algebraic(
    ExactRealAlgebraic& value,
    ComputationContext& context,
    const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (value.lower == value.upper) return Result<void>::success();

    const Rational lower_value = value.polynomial.eval(value.lower);
    if (lower_value == Rational(0)) {
        value.upper = value.lower;
        return Result<void>::success();
    }
    const Rational upper_value = value.polynomial.eval(value.upper);
    if (upper_value == Rational(0)) {
        value.lower = value.upper;
        return Result<void>::success();
    }

    const Rational midpoint = (value.lower + value.upper) / Rational(2);
    const Rational midpoint_value = value.polynomial.eval(midpoint);
    if (midpoint_value == Rational(0)) {
        value.lower = midpoint;
        value.upper = midpoint;
        return Result<void>::success();
    }

    if (sign(lower_value) != sign(midpoint_value)) {
        value.upper = midpoint;
    } else {
        value.lower = midpoint;
    }
    return Result<void>::success();
}

Result<void> refine_exact_real_algebraic_to_tolerance(
    ExactRealAlgebraic& value,
    double absolute_tolerance,
    double relative_tolerance,
    ComputationContext& context,
    const std::string& operation) {
    if (!(absolute_tolerance > 0.0) ||
        !(relative_tolerance > 0.0) ||
        !std::isfinite(absolute_tolerance) ||
        !std::isfinite(relative_tolerance)) {
        return Result<void>::failure(
            CasErrc::InvalidArgument,
            "algebraic refinement tolerances must be finite and positive",
            operation);
    }
    if (value.lower == value.upper) return Result<void>::success();

    try {

        const Rational base = value.lower;
        const Rational span = value.upper - value.lower;
        const double base_double = base.to_double();
        const double upper_double = value.upper.to_double();
        auto evaluate_double = [&](double point) {
            double result = 0.0;
            for (std::size_t position = value.polynomial.coeffs.size();
                 position-- > 0;) {
                result = result * point +
                    value.polynomial.coeffs[position].to_double();
            }
            return result;
        };
        double approximate_lower = base_double;
        double approximate_upper = upper_double;
        double lower_value = evaluate_double(approximate_lower);
        for (int iteration = 0; iteration < 128; ++iteration) {
            const double midpoint =
                (approximate_lower + approximate_upper) * 0.5;
            const double midpoint_value = evaluate_double(midpoint);
            if (!std::isfinite(midpoint_value) || midpoint_value == 0.0) {
                approximate_lower = midpoint;
                approximate_upper = midpoint;
                break;
            }
            if ((lower_value < 0.0) != (midpoint_value < 0.0)) {
                approximate_upper = midpoint;
            } else {
                approximate_lower = midpoint;
                lower_value = midpoint_value;
            }
        }
        const double approximation =
            (approximate_lower + approximate_upper) * 0.5;
        const double target_tolerance = absolute_tolerance +
            relative_tolerance * std::abs(approximation);
        const double span_double = span.to_double();
        int initial_depth = 1;
        if (span_double > 2.0 * target_tolerance) {
            initial_depth = static_cast<int>(std::ceil(std::log2(
                span_double / (2.0 * target_tolerance)))) + 2;
        }
        if (initial_depth > 60) {
            return Result<void>::failure(
                CasErrc::NumericFailure,
                "requested algebraic tolerance is finer than the double evaluation boundary",
                operation);
        }

        const double fraction = std::clamp(
            (approximation - base_double) / span_double, 0.0, 1.0);
        for (int depth = std::max(1, initial_depth); depth <= 60; ++depth) {
            auto step = context.consume_steps(1, operation);
            if (!step) return step;
            const unsigned long long scale_value = 1ULL << depth;
            unsigned long long center = static_cast<unsigned long long>(
                std::floor(fraction * static_cast<double>(scale_value)));
            if (center >= scale_value) center = scale_value - 1;
            for (int radius = 0; radius <= 16; ++radius) {
                for (int direction : {radius == 0 ? 0 : -1,
                                      radius == 0 ? 0 : 1}) {
                    if (radius == 0 && direction != 0) continue;
                    const long long candidate_signed =
                        static_cast<long long>(center) +
                        static_cast<long long>(direction * radius);
                    if (candidate_signed < 0 ||
                        static_cast<unsigned long long>(candidate_signed) >=
                            scale_value) {
                        continue;
                    }
                    const auto candidate =
                        static_cast<unsigned long long>(candidate_signed);
                    const BigInt scale_integer(scale_value);
                    const Rational lower = base + span * Rational(
                        BigInt(candidate), scale_integer);
                    const Rational upper = base + span * Rational(
                        BigInt(candidate + 1), scale_integer);
                    auto exact_count = count_real_roots_exact(
                        value.polynomial, lower, upper,
                        context, operation);
                    if (!exact_count) {
                        return Result<void>::failure(exact_count.error());
                    }
                    const bool contains_root = exact_count.value() == 1;
                    if (contains_root) {
                        value.lower = lower;
                        value.upper = upper;
                        return Result<void>::success();
                    }
                    if (radius == 0) break;
                }
            }
        }
        return Result<void>::failure(
            CasErrc::InternalInvariant,
            "exact projection refinement could not certify the suggested bracket",
            operation);
    } catch (const std::bad_alloc&) {
        return Result<void>::failure(
            CasErrc::ResourceLimit,
            "algebraic refinement allocation failed",
            operation);
    } catch (const std::exception& error) {
        return Result<void>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

Result<bool> equal_exact_real_algebraic(
    ExactRealAlgebraic lhs,
    ExactRealAlgebraic rhs,
    ComputationContext& context) {
    auto step = context.consume_steps(1, "exact_algebraic.equal");
    if (!step) return Result<bool>::failure(step.error());
    if (!overlaps(lhs.lower, lhs.upper, rhs.lower, rhs.upper)) {
        return Result<bool>::success(false);
    }
    if (same_polynomial(lhs.polynomial, rhs.polynomial)) {
        return Result<bool>::success(lhs.root_index == rhs.root_index);
    }

    const Rational overlap_lower =
        lhs.lower < rhs.lower ? rhs.lower : lhs.lower;
    const Rational overlap_upper =
        lhs.upper < rhs.upper ? lhs.upper : rhs.upper;
    if (overlap_lower == overlap_upper &&
        (lhs.polynomial.eval(overlap_lower) != Rational(0) ||
         rhs.polynomial.eval(overlap_lower) != Rational(0))) {
        return Result<bool>::success(false);
    }

    try {
        auto common = monic(Polynomial<Rational>::gcd(
            lhs.polynomial, rhs.polynomial));
        if (common.degree() < 1) return Result<bool>::success(false);
        auto isolated = isolate_real_roots_exact(
            common, context, "exact_algebraic.equal");
        if (!isolated) return Result<bool>::failure(isolated.error());
        for (const auto& interval : isolated.value()) {
            if (overlaps(lhs.lower, lhs.upper, interval.first, interval.second) &&
                overlaps(rhs.lower, rhs.upper, interval.first, interval.second)) {
                return Result<bool>::success(true);
            }
        }
        return Result<bool>::success(false);
    } catch (const std::bad_alloc&) {
        return Result<bool>::failure(
            CasErrc::ResourceLimit,
            "algebraic equality allocation failed",
            "exact_algebraic.equal");
    } catch (const std::exception& error) {
        return Result<bool>::failure(
            CasErrc::InternalInvariant,
            error.what(),
            "exact_algebraic.equal");
    }
}

Result<int> compare_exact_real_algebraic(
    ExactRealAlgebraic lhs,
    ExactRealAlgebraic rhs,
    ComputationContext& context) {
    auto equal = equal_exact_real_algebraic(lhs, rhs, context);
    if (!equal) return Result<int>::failure(equal.error());
    if (equal.value()) return Result<int>::success(0);

    while (true) {
        if (lhs.upper < rhs.lower) return Result<int>::success(-1);
        if (rhs.upper < lhs.lower) return Result<int>::success(1);
        const double lhs_width = (lhs.upper - lhs.lower).to_double();
        const double rhs_width = (rhs.upper - rhs.lower).to_double();
        const double width = std::max(lhs_width, rhs_width);
        if (!(width > 0.0) || !std::isfinite(width)) {
            return Result<int>::failure(
                CasErrc::InternalInvariant,
                "distinct algebraic values retained an inseparable enclosure",
                "exact_algebraic.compare");
        }
        if (lhs_width > 0.0) {
            auto refined = refine_exact_real_algebraic_to_tolerance(
                lhs, width / 16.0, 1e-15, context,
                "exact_algebraic.compare");
            if (!refined) return Result<int>::failure(refined.error());
        }
        if (rhs_width > 0.0) {
            auto refined = refine_exact_real_algebraic_to_tolerance(
                rhs, width / 16.0, 1e-15, context,
                "exact_algebraic.compare");
            if (!refined) return Result<int>::failure(refined.error());
        }
    }
}

} // namespace LMCAS::detail
