#include "internal/interval_endpoint.hpp"

namespace LMCAS::detail {
namespace {

int rational_sign(const Rational& value) {
    if (value == Rational(0)) {
        return 0;
    }
    return value.get_numerator().is_negative() ? -1 : 1;
}

int compare_rationals(const Rational& left, const Rational& right) {
    if (left < right) {
        return -1;
    }
    if (left > right) {
        return 1;
    }
    return 0;
}

int compare_algebraic_to_rational(
    const ExactRealAlgebraic& algebraic, const Rational& rational) {
            if (algebraic.upper < rational) {
                return -1;
            }
            if (rational < algebraic.lower) {
                return 1;
            }
            const Rational value = algebraic.polynomial.eval(rational);
            if (value == Rational(0)) {
                return 0;
            }
            const Rational lower_value =
                algebraic.polynomial.eval(algebraic.lower);
            const bool root_before_rational =
                rational_sign(lower_value) != rational_sign(value);
            return root_before_rational ? -1 : 1;
}

ExactRealAlgebraicResult to_algebraic(
    const ComparableEndpoint& endpoint, ComputationContext& context) {
        if (endpoint.algebraic) {
            return LMCAS::detail::ExactRealAlgebraicResult::success(
                *endpoint.algebraic);
        }
        if (endpoint.radical_coefficient == Rational(0)) {
            return LMCAS::detail::make_exact_real_algebraic(
                Polynomial<Rational>({
                    Rational(0) - endpoint.rational,
                    Rational(1)
                }),
                0, 1, context);
        }
        const Rational a = endpoint.rational;
        const Rational b = endpoint.radical_coefficient;
        const Rational d = endpoint.radicand;
        auto value = LMCAS::detail::make_exact_real_algebraic(
            Polynomial<Rational>({
                a * a - b * b * d,
                Rational(-2) * a,
                Rational(1)
            }),
            b > Rational(0) ? 1 : 0, 1, context);
        return value;
}

}

Result<int> compare_comparable(const ComparableEndpoint& left,
                               const ComparableEndpoint& right,
                               ComputationContext& context) {
    if (left.infinity < right.infinity) {
        return Result<int>::success(-1);
    }
    if (left.infinity > right.infinity) {
        return Result<int>::success(1);
    }
    if (left.infinity != 0) {
        return Result<int>::success(0);
    }

    const bool left_radical =
        left.radical_coefficient != Rational(0) || left.algebraic.has_value();
    const bool right_radical =
        right.radical_coefficient != Rational(0) || right.algebraic.has_value();
    if (!left_radical && !right_radical) {
        return Result<int>::success(compare_rationals(left.rational, right.rational));
    }
    if (left.algebraic && !right_radical) {
        return Result<int>::success(
            compare_algebraic_to_rational(*left.algebraic, right.rational));
    }
    if (right.algebraic && !left_radical) {
        return Result<int>::success(
            -compare_algebraic_to_rational(*right.algebraic, left.rational));
    }
    auto left_algebraic = to_algebraic(left, context);
    if (!left_algebraic) {
        return Result<int>::failure(left_algebraic.error());
    }
    auto right_algebraic = to_algebraic(right, context);
    if (!right_algebraic) {
        return Result<int>::failure(right_algebraic.error());
    }
    return LMCAS::detail::compare_exact_real_algebraic(
        left_algebraic.value(), right_algebraic.value(), context);
}

}
