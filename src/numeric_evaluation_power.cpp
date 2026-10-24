#include "internal/numeric_evaluation_support.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

namespace LMCAS::detail {
namespace {

struct ExponentParity {
    bool integral;
    bool odd;
};

ExponentParity exponent_parity(const PowerNode& power, double value) {
    const auto* number = dynamic_cast<const NumberNode*>(power.exponent().get());
    const auto* integer = number ? std::get_if<BigInt>(&number->value()) : nullptr;
    if (integer) {
        return {true, integer->is_odd()};
    }
    const auto* rational = number ? std::get_if<Rational>(&number->value()) : nullptr;
    if (rational) {
        const bool integral = rational->is_integer();
        return {integral, integral && rational->get_numerator().is_odd()};
    }
    const bool integral = std::trunc(value) == value;
    return {integral, integral && std::fmod(value, 2.0) != 0.0};
}

}

Result<ApproxReal> NumericEvaluator::square_root(const NumberNode& node) {
    auto step = context_.consume_steps(1, "evaluate_numeric");
    if (!step) return Result<ApproxReal>::failure(step.error());
    if (node.is_zero()) return numeric_approximation(0.0);
    BigInt numerator;
    BigInt denominator(1);
    if (const auto* integer = std::get_if<BigInt>(&node.value())) {
        numerator = *integer;
    } else {
        const auto& rational = std::get<Rational>(node.value());
        numerator = rational.get_numerator();
        denominator = rational.get_denominator();
    }
    if (numerator.is_negative()) {
        return numeric_failure(CasErrc::DomainError,
                               "real square root requires a non-negative argument");
    }
    const auto num_bits = numerator.bit_length();
    const auto den_bits = denominator.bit_length();
    auto bits = context_.require_integer_bits(std::max(num_bits, den_bits), "evaluate_numeric");
    if (!bits) return Result<ApproxReal>::failure(bits.error());
    const bool reciprocal_scale = num_bits < den_bits;
    const auto difference = reciprocal_scale ? den_bits - num_bits : num_bits - den_bits;
    constexpr int max_scale = 2 * std::numeric_limits<double>::max_exponent;
    constexpr int min_scale = 2 * (std::numeric_limits<double>::digits -
                                   std::numeric_limits<double>::min_exponent + 1);
    if (difference > static_cast<std::size_t>(reciprocal_scale ? min_scale : max_scale)) {
        return numeric_failure(CasErrc::NumericFailure, "numeric square root is outside the floating-point range");
    }
    const auto shift = difference - difference % 2;
    if (reciprocal_scale) numerator <<= shift;
    else denominator <<= shift;
    const int scale = static_cast<int>(shift / 2) * (reciprocal_scale ? -1 : 1);
    const double value = std::ldexp(std::sqrt(Rational(numerator, denominator).to_double()), scale);
    if (value == 0.0) {
        return numeric_failure(CasErrc::NumericFailure, "numeric square root underflow");
    }
    return numeric_approximation(value);
}

Result<ApproxReal> NumericEvaluator::power(const PowerNode& node) {
    const auto* base_number = dynamic_cast<const NumberNode*>(node.base().get());
    const auto* exponent_number = dynamic_cast<const NumberNode*>(node.exponent().get());
    if (base_number && exponent_number) {
        const auto* exponent = std::get_if<Rational>(&exponent_number->value());
        if (exponent && *exponent == Rational(1, 2)) {
            auto step = context_.consume_steps(1, "evaluate_numeric");
            if (!step) return Result<ApproxReal>::failure(step.error());
            return square_root(*base_number);
        }
        if (exponent_number->is_negative_one()) {
            auto step = context_.consume_steps(2, "evaluate_numeric");
            if (!step) return Result<ApproxReal>::failure(step.error());
            if (base_number->is_zero()) {
                return numeric_failure(CasErrc::DomainError, "zero cannot be raised to a non-positive power");
            }
            try {
                if (const auto* integer = std::get_if<BigInt>(&base_number->value())) {
                    return numeric_approximation(Rational(BigInt(1), *integer).to_double());
                }
                const auto& rational = std::get<Rational>(base_number->value());
                return numeric_approximation(
                    Rational(rational.get_denominator(), rational.get_numerator()).to_double());
            } catch (const std::range_error& error) {
                return numeric_failure(CasErrc::NumericFailure, error.what());
            } catch (const std::overflow_error& error) {
                return numeric_failure(CasErrc::NumericFailure, error.what());
            } catch (const std::underflow_error& error) {
                return numeric_failure(CasErrc::NumericFailure, error.what());
            }
        }
    }
    auto base = evaluate(node.base());
    if (!base) {
        return base;
    }
    auto exponent = evaluate(node.exponent());
    if (!exponent) {
        return exponent;
    }
    const double base_value = base.value().value;
    const double exponent_value = exponent.value().value;
    if (base_value == 0.0 && exponent_value <= 0.0) {
        return numeric_failure(CasErrc::DomainError, "zero cannot be raised to a non-positive power");
    }
    if (std::signbit(base_value) && std::isfinite(exponent_value)) {
        const auto parity = exponent_parity(node, exponent_value);
        if (base_value < 0.0 && std::isfinite(base_value) && !parity.integral) {
            return numeric_failure(CasErrc::DomainError,
                                   "a negative real base requires an integer exponent");
        }
        if (parity.integral) {
            const double magnitude = std::pow(std::abs(base_value), exponent_value);
            if (magnitude == 0.0 && base_value != 0.0) {
                return numeric_failure(CasErrc::NumericFailure, "numeric power underflow");
            }
            return numeric_approximation(parity.odd ? -magnitude : magnitude);
        }
    }
    const double value = std::pow(base_value, exponent_value);
    if (value == 0.0 && base_value != 0.0) {
        return numeric_failure(CasErrc::NumericFailure, "numeric power underflow");
    }
    return numeric_approximation(value);
}

}
