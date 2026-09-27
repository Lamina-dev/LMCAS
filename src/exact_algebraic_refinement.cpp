#include "internal/exact_algebraic.hpp"
#include "internal/exact_sturm.hpp"

#include <algorithm>
#include <cmath>
#include <exception>

namespace LMCAS::detail {
namespace {

struct DyadicProposal {
    Rational base;
    Rational span;
    Rational fraction;
    int initial_depth;
};

double evaluate_polynomial_double(const Polynomial<Rational>& polynomial,
                                  double point) {
    double result = 0.0;
    for (std::size_t position = polynomial.coeffs.size(); position-- > 0;) {
        result = result * point + polynomial.coeffs[position].to_double();
    }
    return result;
}

double approximate_root(const ExactRealAlgebraic& value) {
    double lower = value.lower.to_double();
    double upper = value.upper.to_double();
    double lower_value = evaluate_polynomial_double(value.polynomial, lower);
    for (int iteration = 0; iteration < 128; ++iteration) {
        const double midpoint = (lower + upper) * 0.5;
        const double midpoint_value = evaluate_polynomial_double(value.polynomial, midpoint);
        if (!std::isfinite(midpoint_value) || midpoint_value == 0.0) {
            lower = midpoint;
            upper = midpoint;
            break;
        }
        if ((lower_value < 0.0) != (midpoint_value < 0.0)) {
            upper = midpoint;
        } else {
            lower = midpoint;
            lower_value = midpoint_value;
        }
    }
    return (lower + upper) * 0.5;
}

Result<DyadicProposal> propose_dyadic_bracket(
    const ExactRealAlgebraic& value, double absolute_tolerance,
    double relative_tolerance, const std::string& operation) {
    const Rational base = value.lower;
    const Rational span = value.upper - value.lower;
    const double approximation = approximate_root(value);
    const double target_tolerance = absolute_tolerance +
        relative_tolerance * std::abs(approximation);
    const double span_double = span.to_double();
    double requested_depth = 1.0;
    if (span_double > 2.0 * target_tolerance) {
        requested_depth = std::ceil(std::log2(
            span_double / (2.0 * target_tolerance))) + 2.0;
    }
    if (!std::isfinite(approximation) || !std::isfinite(target_tolerance) ||
        !std::isfinite(span_double) || !std::isfinite(requested_depth) ||
        requested_depth > 60.0) {
        return Result<DyadicProposal>::failure(
            CasErrc::NumericFailure,
            "requested algebraic tolerance is finer than the double evaluation boundary",
            operation);
    }
    const int initial_depth = static_cast<int>(std::max(1.0, requested_depth));
    const Rational fraction = std::clamp(
        (Rational::from_double(approximation) - base) / span,
        Rational(0), Rational(1));
    return Result<DyadicProposal>::success(
        DyadicProposal{base, span, fraction, initial_depth});
}

Result<bool> certify_dyadic_candidate(
    ExactRealAlgebraic& value, const DyadicProposal& proposal,
    const BigInt& candidate, const BigInt& scale,
    ComputationContext& context, const std::string& operation) {
    const Rational lower = proposal.base + proposal.span * Rational(candidate, scale);
    const Rational upper = proposal.base + proposal.span * Rational(candidate + BigInt(1), scale);
    auto count = count_real_roots_exact(value.polynomial, lower, upper, context, operation);
    if (!count) return Result<bool>::failure(count.error());
    if (count.value() != 1) return Result<bool>::success(false);
    value.lower = lower;
    value.upper = upper;
    return Result<bool>::success(true);
}

Result<bool> search_dyadic_scale(
    ExactRealAlgebraic& value, const DyadicProposal& proposal, int depth,
    ComputationContext& context, const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) {
        return Result<bool>::failure(step.error());
    }
    auto bits = context.require_integer_bits(static_cast<std::size_t>(depth) + 1, operation);
    if (!bits) {
        return Result<bool>::failure(bits.error());
    }
    const BigInt scale = BigInt(1) << depth;
    BigInt center = (proposal.fraction.get_numerator() * scale) /
        proposal.fraction.get_denominator();
    if (center >= scale) {
        center = scale - BigInt(1);
    }
    for (int radius = 0; radius <= 16; ++radius) {
        for (int direction : {radius == 0 ? 0 : -1, radius == 0 ? 0 : 1}) {
            const BigInt candidate = center + BigInt(direction * radius);
            if (candidate.is_negative() || candidate >= scale) {
                continue;
            }
            auto certified = certify_dyadic_candidate(
                value, proposal, candidate, scale, context, operation);
            if (!certified || certified.value()) {
                return certified;
            }
            if (radius == 0) {
                break;
            }
        }
    }
    return Result<bool>::success(false);
}

}

Result<void> refine_exact_real_algebraic_to_tolerance(
    ExactRealAlgebraic& value, double absolute_tolerance, double relative_tolerance,
    ComputationContext& context, const std::string& operation) {
    if (!(absolute_tolerance > 0.0) || !(relative_tolerance > 0.0) ||
        !std::isfinite(absolute_tolerance) || !std::isfinite(relative_tolerance)) {
        return Result<void>::failure(
            CasErrc::InvalidArgument,
            "algebraic refinement tolerances must be finite and positive", operation);
    }
    if (value.lower == value.upper) {
        return Result<void>::success();
    }
    try {
        auto proposal = propose_dyadic_bracket(
            value, absolute_tolerance, relative_tolerance, operation);
        if (!proposal) {
            return Result<void>::failure(proposal.error());
        }
        for (int depth = std::max(1, proposal.value().initial_depth); depth <= 60; ++depth) {
            auto certified = search_dyadic_scale(value, proposal.value(), depth, context, operation);
            if (!certified) {
                return Result<void>::failure(certified.error());
            }
            if (certified.value()) {
                return Result<void>::success();
            }
        }
        return Result<void>::failure(
            CasErrc::InternalInvariant,
            "exact projection refinement could not certify the suggested bracket", operation);
    } catch (const std::bad_alloc&) {
        return Result<void>::failure(
            CasErrc::ResourceLimit, "algebraic refinement allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<void>::failure(CasErrc::InternalInvariant, error.what(), operation);
    }
}

}
