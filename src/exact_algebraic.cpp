#include "internal/exact_algebraic.hpp"

#include "internal/exact_sturm.hpp"
#include "rational_polynomial.hpp"

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

Result<bool> common_root_in_enclosures(
    const ExactRealAlgebraic& lhs, const ExactRealAlgebraic& rhs,
    ComputationContext& context) {
    auto common = monic(Polynomial<Rational>::gcd(lhs.polynomial, rhs.polynomial));
    if (common.degree() < 1) {
        return Result<bool>::success(false);
    }
    auto isolated = isolate_real_roots_exact(
        common, context, "exact_algebraic.equal");
    if (!isolated) {
        return Result<bool>::failure(isolated.error());
    }
    for (const auto& interval : isolated.value()) {
        if (overlaps(lhs.lower, lhs.upper, interval.first, interval.second) &&
            overlaps(rhs.lower, rhs.upper, interval.first, interval.second)) {
            return Result<bool>::success(true);
        }
    }
    return Result<bool>::success(false);
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


Result<bool> equal_exact_real_algebraic(
    ExactRealAlgebraic lhs,
    ExactRealAlgebraic rhs,
    ComputationContext& context) {
    auto step = context.consume_steps(1, "exact_algebraic.equal");
    if (!step) {
        return Result<bool>::failure(step.error());
    }
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
        return common_root_in_enclosures(lhs, rhs, context);
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
    if (!equal) {
        return Result<int>::failure(equal.error());
    }
    if (equal.value()) {
        return Result<int>::success(0);
    }

    while (true) {
        if (lhs.upper < rhs.lower) {
            return Result<int>::success(-1);
        }
        if (rhs.upper < lhs.lower) {
            return Result<int>::success(1);
        }
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
            if (!refined) {
                return Result<int>::failure(refined.error());
            }
        }
        if (rhs_width > 0.0) {
            auto refined = refine_exact_real_algebraic_to_tolerance(
                rhs, width / 16.0, 1e-15, context,
                "exact_algebraic.compare");
            if (!refined) {
                return Result<int>::failure(refined.error());
            }
        }
    }
}

} // namespace LMCAS::detail
