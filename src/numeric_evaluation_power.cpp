#include "internal/numeric_evaluation_support.hpp"
#include <cmath>

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

Result<ApproxReal> NumericEvaluator::power(const PowerNode& node) {
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
            return numeric_approximation(parity.odd ? -magnitude : magnitude);
        }
    }
    return numeric_approximation(std::pow(base_value, exponent_value));
}

}
