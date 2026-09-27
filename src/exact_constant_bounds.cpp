#include "internal/exact_constant_bounds.hpp"
#include <new>
#include <stdexcept>

namespace LMCAS::detail {
namespace {
using Bounds = ExactConstantInterval;

/**
 * @brief 用交错部分和夹逼 atan(x)。
 *
 * |x|≤1/2 时相邻项绝对值递减；相邻交错部分和夹住 atan(x)，差值恰为首个省略项。
 */
Bounds small_atan(const Rational& x, std::size_t precision,
                  const ExactBoundArithmetic& arithmetic) {
    const Rational target(BigInt(1), arithmetic.dyadic_scale(precision));
    if (x.is_zero()) return {Rational(0), Rational(0)};
    if (x.get_numerator().is_negative()) {
        auto bounds = small_atan(-x, precision, arithmetic);
        return {-bounds.second, -bounds.first};
    }
    const Rational square = arithmetic.mul(x, x);
    Rational power = x;
    Rational sum = x;
    BigInt odd(1);
    bool subtract = true;
    for (;;) {
        arithmetic.step();
        arithmetic.growth(odd.bit_length(), 2);
        odd += BigInt(2);
        power = arithmetic.mul(power, square);
        const Rational omitted = arithmetic.div(power, Rational(odd));
        const Rational adjacent = subtract ? arithmetic.sub(sum, omitted)
                                           : arithmetic.add(sum, omitted);
        if (!arithmetic.less(target, omitted)) {
            return subtract ? Bounds{adjacent, sum} : Bounds{sum, adjacent};
        }
        sum = adjacent;
        subtract = !subtract;
    }
}

Bounds pi_bounds(std::size_t precision, const ExactBoundArithmetic& arithmetic) {
    /**
     * @brief 以小参数 atan 区间构造 π 界。
     *
     * 每个 atan 区间宽度≤2^-(p+5)，组合宽度≤(16+4)*2^-(p+5)<2^-p。
     */
    const auto internal = arithmetic.growth(precision, 5);
    auto fifth = small_atan(Rational(1, 5), internal, arithmetic);
    auto two_thirty_ninth = small_atan(Rational(1, 239), internal, arithmetic);
    return {arithmetic.sub(arithmetic.mul(Rational(16), fifth.first),
                           arithmetic.mul(Rational(4), two_thirty_ninth.second)),
            arithmetic.sub(arithmetic.mul(Rational(16), fifth.second),
                           arithmetic.mul(Rational(4), two_thirty_ninth.first))};
}

Bounds atan_bounds(const Rational& x, std::size_t precision,
                   const ExactBoundArithmetic& arithmetic) {
    arithmetic.step();
    if (x.get_numerator().is_negative()) {
        auto bounds = atan_bounds(-x, precision, arithmetic);
        return {-bounds.second, -bounds.first};
    }
    if (!arithmetic.less(Rational(1, 2), x))
        return small_atan(x, precision, arithmetic);
    const auto internal = arithmetic.growth(precision, 2);
    auto pi = pi_bounds(internal, arithmetic);
    if (arithmetic.less(Rational(1), x)) {
        auto reciprocal = atan_bounds(arithmetic.div(Rational(1), x), internal, arithmetic);
        return {arithmetic.sub(arithmetic.div(pi.first, Rational(2)), reciprocal.second),
                arithmetic.sub(arithmetic.div(pi.second, Rational(2)), reciprocal.first)};
    }
    const Rational transformed = arithmetic.div(arithmetic.sub(x, Rational(1)),
                                                  arithmetic.add(x, Rational(1)));
    auto correction = small_atan(transformed, internal, arithmetic);
    return {arithmetic.add(arithmetic.div(pi.first, Rational(4)), correction.first),
            arithmetic.add(arithmetic.div(pi.second, Rational(4)), correction.second)};
}

template<class Function>
Result<Bounds> checked_bounds(ComputationContext& context, const char* operation,
                              Function&& function) {
    try {
        ExactBoundArithmetic arithmetic(context, operation);
        arithmetic.step();
        return function(arithmetic);
    } catch (const CasError& error) {
        return Result<Bounds>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<Bounds>::failure(CasErrc::ResourceLimit,
                                       "constant enclosure allocation failed", operation);
    } catch (const std::length_error&) {
        return Result<Bounds>::failure(CasErrc::ResourceLimit,
                                       "constant enclosure size exceeded", operation);
    }
}
}

Result<ExactConstantInterval> exact_atan_bounds(
    const Rational& value, std::size_t precision_bits, ComputationContext& context) {
    return checked_bounds(context, "exact.atan.bounds", [&](const ExactBoundArithmetic& arithmetic) {
        return atan_bounds(value, precision_bits, arithmetic);
    });
}

Result<ExactConstantInterval> exact_pi_bounds(
    std::size_t precision_bits, ComputationContext& context) {
    return checked_bounds(context, "exact.pi.bounds", [&](const ExactBoundArithmetic& arithmetic) {
        return pi_bounds(precision_bits, arithmetic);
    });
}

Result<ExactConstantInterval> exact_e_bounds(
    std::size_t precision_bits, ComputationContext& context) {
    return checked_bounds(context, "exact.e.bounds", [&](const ExactBoundArithmetic& arithmetic) {
        const Rational target(BigInt(1), arithmetic.dyadic_scale(precision_bits));
        BigInt n(1);
        Rational term(1); /**< 级数项 1/n!。 */
        Rational sum(2); /**< 含 k=0 项的部分和 S_n。 */
        for (;;) {
            arithmetic.step();
            const Rational tail = arithmetic.div(term, Rational(n));
            if (!arithmetic.less(target, tail))
                return Bounds{sum, arithmetic.add(sum, tail)};
            arithmetic.growth(n.bit_length(), 1);
            n += BigInt(1);
            term = arithmetic.div(term, Rational(n));
            sum = arithmetic.add(sum, term);
        }
    });
}
}
