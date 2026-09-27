#include "internal/exact_sturm.hpp"

#include <algorithm>
#include <exception>

namespace LMCAS::detail {
namespace {

using SturmSequence = std::vector<Polynomial<Rational>>;

int sign_changes(const std::vector<int>& signs) {
    int changes = 0;
    int previous = 0;
    for (const int current : signs) {
        if (current == 0) continue;
        if (previous != 0 && current != previous) ++changes;
        previous = current;
    }
    return changes;
}

int polynomial_sign_at(
    const Polynomial<Rational>& polynomial,
    const Rational& point) {
    if (polynomial.is_zero()) return 0;
    const int degree = polynomial.degree();
    const BigInt numerator = point.get_numerator();
    const BigInt denominator = point.get_denominator();
    BigInt coefficient_denominator(1);
    for (const auto& coefficient : polynomial.coeffs) {
        coefficient_denominator *= coefficient.get_denominator();
    }
    std::vector<BigInt> numerator_powers(
        static_cast<std::size_t>(degree) + 1, BigInt(1));
    std::vector<BigInt> denominator_powers(
        static_cast<std::size_t>(degree) + 1, BigInt(1));
    for (int exponent = 1; exponent <= degree; ++exponent) {
        numerator_powers[static_cast<std::size_t>(exponent)] =
            numerator_powers[static_cast<std::size_t>(exponent - 1)] *
            numerator;
        denominator_powers[static_cast<std::size_t>(exponent)] =
            denominator_powers[static_cast<std::size_t>(exponent - 1)] *
            denominator;
    }
    BigInt value(0);
    for (int exponent = 0; exponent <= degree; ++exponent) {
        const auto& coefficient =
            polynomial.coeffs[static_cast<std::size_t>(exponent)];
        if (coefficient == Rational(0)) continue;
        value += coefficient.get_numerator() *
            (coefficient_denominator / coefficient.get_denominator()) *
            numerator_powers[static_cast<std::size_t>(exponent)] *
            denominator_powers[
                static_cast<std::size_t>(degree - exponent)];
    }
    if (value == BigInt(0)) return 0;
    return value < BigInt(0) ? -1 : 1;
}

Result<SturmSequence> make_sturm_sequence(
    Polynomial<Rational> polynomial,
    ComputationContext& context,
    const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) return Result<SturmSequence>::failure(step.error());

    polynomial.trim();
    if (polynomial.degree() <= 0) {
        return Result<SturmSequence>::success({});
    }
    polynomial = polynomial.square_free_part().make_monic();

    SturmSequence sequence;
    sequence.push_back(polynomial);
    sequence.push_back(polynomial.differentiate());
    while (!sequence.back().is_zero()) {
        step = context.consume_steps(1, operation);
        if (!step) return Result<SturmSequence>::failure(step.error());
        const auto& dividend = sequence[sequence.size() - 2];
        const auto& divisor = sequence.back();
        auto division = dividend.div_mod(divisor);
        auto remainder = std::move(division.second);
        if (remainder.is_zero()) break;
        for (auto& coefficient : remainder.coeffs) coefficient = -coefficient;
        remainder.trim();
        sequence.push_back(std::move(remainder));
    }
    if (!sequence.empty() && sequence.back().is_zero()) sequence.pop_back();
    return Result<SturmSequence>::success(std::move(sequence));
}

Result<int> variations_at(
    const SturmSequence& sequence,
    const Rational& point,
    ComputationContext& context,
    const std::string& operation) {
    auto step = context.consume_steps(sequence.size() + 1, operation);
    if (!step) return Result<int>::failure(step.error());
    std::vector<int> signs;
    signs.reserve(sequence.size());
    for (const auto& member : sequence) {
        signs.push_back(polynomial_sign_at(member, point));
    }
    return Result<int>::success(sign_changes(signs));
}

Result<std::size_t> count_open_roots(
    const Polynomial<Rational>& square_free,
    const SturmSequence& sequence,
    const Rational& lower,
    const Rational& upper,
    ComputationContext& context,
    const std::string& operation) {
    if (!(lower < upper) || sequence.empty()) {
        return Result<std::size_t>::success(0);
    }
    auto lower_variations = variations_at(sequence, lower, context, operation);
    if (!lower_variations) {
        return Result<std::size_t>::failure(lower_variations.error());
    }
    auto upper_variations = variations_at(sequence, upper, context, operation);
    if (!upper_variations) {
        return Result<std::size_t>::failure(upper_variations.error());
    }
    int count = lower_variations.value() - upper_variations.value();
    if (polynomial_sign_at(square_free, upper) == 0) --count;
    if (count < 0) {
        return Result<std::size_t>::failure(
            CasErrc::InternalInvariant,
            "Sturm variation produced a negative open-interval count",
            operation);
    }
    return Result<std::size_t>::success(static_cast<std::size_t>(count));
}


struct PendingInterval {
    Rational lower;
    Rational upper;
    std::size_t count;
};

Result<void> clear_isolation_endpoints(
    const Polynomial<Rational>& square_free, const SturmSequence& sequence,
    PendingInterval& interval, ComputationContext& context,
    const std::string& operation) {
    while (polynomial_sign_at(square_free, interval.lower) == 0 ||
           polynomial_sign_at(square_free, interval.upper) == 0) {
        const Rational midpoint = (interval.lower + interval.upper) / Rational(2);
        if (polynomial_sign_at(square_free, midpoint) == 0) {
            interval.lower = midpoint;
            interval.upper = midpoint;
            break;
        }
        const bool lower_is_root = polynomial_sign_at(square_free, interval.lower) == 0;
        auto half = lower_is_root
            ? count_open_roots(square_free, sequence, midpoint, interval.upper,
                               context, operation)
            : count_open_roots(square_free, sequence, interval.lower, midpoint,
                               context, operation);
        if (!half) return Result<void>::failure(half.error());
        if (lower_is_root == (half.value() == 1)) {
            interval.lower = midpoint;
        } else {
            interval.upper = midpoint;
        }
    }
    return Result<void>::success();
}

Result<void> subdivide_isolation(
    const Polynomial<Rational>& square_free, const SturmSequence& sequence,
    const PendingInterval& current, std::vector<PendingInterval>& pending,
    std::vector<RationalInterval>& isolated, ComputationContext& context,
    const std::string& operation) {
    const Rational midpoint = (current.lower + current.upper) / Rational(2);
    const bool midpoint_is_root = polynomial_sign_at(square_free, midpoint) == 0;
    auto left = count_open_roots(
        square_free, sequence, current.lower, midpoint, context, operation);
    if (!left) return Result<void>::failure(left.error());
    auto right = count_open_roots(
        square_free, sequence, midpoint, current.upper, context, operation);
    if (!right) return Result<void>::failure(right.error());
    const std::size_t accounted = left.value() + right.value() +
        (midpoint_is_root ? 1U : 0U);
    if (accounted != current.count) {
        return Result<void>::failure(
            CasErrc::InternalInvariant,
            "Sturm subdivision did not preserve the root count", operation);
    }
    if (right.value() != 0) {
        pending.push_back({midpoint, current.upper, right.value()});
    }
    if (midpoint_is_root) isolated.emplace_back(midpoint, midpoint);
    if (left.value() != 0) {
        pending.push_back({current.lower, midpoint, left.value()});
    }
    return Result<void>::success();
}

Result<std::vector<RationalInterval>> isolate_sturm_sequence(
    const Polynomial<Rational>& square_free, const SturmSequence& sequence,
    ComputationContext& context, const std::string& operation) {
    const Rational bound = strict_cauchy_root_bound(square_free);
    auto total = count_open_roots(square_free, sequence, -bound, bound, context, operation);
    if (!total) return Result<std::vector<RationalInterval>>::failure(total.error());
    std::vector<PendingInterval> pending;
    std::vector<RationalInterval> isolated;
    if (total.value() != 0) pending.push_back({-bound, bound, total.value()});
    isolated.reserve(total.value());
    while (!pending.empty()) {
        auto step = context.consume_steps(1, operation);
        if (!step) return Result<std::vector<RationalInterval>>::failure(step.error());
        PendingInterval current = std::move(pending.back());
        pending.pop_back();
        if (current.count == 1) {
            auto cleared = clear_isolation_endpoints(
                square_free, sequence, current, context, operation);
            if (!cleared) return Result<std::vector<RationalInterval>>::failure(cleared.error());
            isolated.emplace_back(current.lower, current.upper);
            continue;
        }
        auto subdivided = subdivide_isolation(
            square_free, sequence, current, pending, isolated, context, operation);
        if (!subdivided) return Result<std::vector<RationalInterval>>::failure(subdivided.error());
    }
    std::sort(isolated.begin(), isolated.end(),
        [](const RationalInterval& left, const RationalInterval& right) {
            if (left.first != right.first) return left.first < right.first;
            return left.second < right.second;
        });
    if (isolated.size() != total.value()) {
        return Result<std::vector<RationalInterval>>::failure(
            CasErrc::InternalInvariant,
            "Sturm isolation lost a distinct real root", operation);
    }
    return Result<std::vector<RationalInterval>>::success(std::move(isolated));
}

} // namespace

Result<std::size_t> count_real_roots_exact(
    const Polynomial<Rational>& polynomial,
    const Rational& lower,
    const Rational& upper,
    ComputationContext& context,
    const std::string& operation) {
    if (upper < lower) {
        return Result<std::size_t>::failure(
            CasErrc::InvalidArgument,
            "real-root interval bounds must be ordered",
            operation);
    }
    try {
        auto step = context.consume_steps(1, operation);
        if (!step) return Result<std::size_t>::failure(step.error());
        if (polynomial.degree() <= 0) {
            return Result<std::size_t>::success(0);
        }
        auto square_free = polynomial.square_free_part().make_monic();
        if (lower == upper) {
            return Result<std::size_t>::success(
                polynomial_sign_at(square_free, lower) == 0 ? 1U : 0U);
        }
        auto sequence = make_sturm_sequence(square_free, context, operation);
        if (!sequence) return Result<std::size_t>::failure(sequence.error());
        auto interior = count_open_roots(
            square_free, sequence.value(), lower, upper, context, operation);
        if (!interior) return interior;
        std::size_t count = interior.value();
        if (polynomial_sign_at(square_free, lower) == 0) ++count;
        if (polynomial_sign_at(square_free, upper) == 0) ++count;
        return Result<std::size_t>::success(count);
    } catch (const std::bad_alloc&) {
        return Result<std::size_t>::failure(
            CasErrc::ResourceLimit,
            "Sturm root counting allocation failed",
            operation);
    } catch (const std::exception& error) {
        return Result<std::size_t>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

Result<std::vector<RationalInterval>> isolate_real_roots_exact(
    const Polynomial<Rational>& polynomial, ComputationContext& context,
    const std::string& operation) {
    try {
        auto step = context.consume_steps(1, operation);
        if (!step) return Result<std::vector<RationalInterval>>::failure(step.error());
        if (polynomial.degree() <= 0) {
            return Result<std::vector<RationalInterval>>::success({});
        }
        if (polynomial.degree() == 1) {
            const Rational root = -polynomial.coeffs[0] / polynomial.coeffs[1];
            return Result<std::vector<RationalInterval>>::success({{root, root}});
        }
        auto square_free = polynomial.square_free_part().make_monic();
        auto sequence = make_sturm_sequence(square_free, context, operation);
        if (!sequence) {
            return Result<std::vector<RationalInterval>>::failure(sequence.error());
        }
        return isolate_sturm_sequence(square_free, sequence.value(), context, operation);
    } catch (const std::bad_alloc&) {
        return Result<std::vector<RationalInterval>>::failure(
            CasErrc::ResourceLimit,
            "Sturm root isolation allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<std::vector<RationalInterval>>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

} // namespace LMCAS::detail
