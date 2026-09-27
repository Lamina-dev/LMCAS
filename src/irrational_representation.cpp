#include "irrational.hpp"
#include "lmmc/numeric.h"
#include <cmath>
#include <ostream>

namespace LMCAS {
namespace {

std::string named_power(const std::string& name, const BigInt& exponent) {
    if (exponent.is_zero()) return {};
    return exponent == BigInt(1) ? name : name + "^(" + exponent.to_string() + ")";
}

void append_factor(std::string& text, const std::string& factor) {
    if (factor.empty()) return;
    if (!text.empty()) text += "*";
    text += factor;
}

}

std::string Irrational::Basis::to_string() const {
    std::string text;
    append_factor(text, named_power("e", e_power));
    append_factor(text, named_power("π", pi_power));
    if (radicand != BigInt(1)) append_factor(text, "√" + radicand.to_string());
    return text;
}

lmmc_real_t Irrational::Basis::to_double() const {
    lmmc_real_t value = 1.0;
    if (!e_power.is_zero()) value *= std::exp(e_power.to_double());
    if (!pi_power.is_zero()) value *= std::pow(LMMC_PI, pi_power.to_double());
    if (radicand != BigInt(1)) value *= std::sqrt(radicand.to_double());
    return value;
}

std::string Irrational::terms_to_string(const Terms& terms) {
    std::string result;
    for (const auto& term : terms) {
        const bool negative = term.second < Rational(0);
        if (result.empty()) {
            if (negative) result = "-";
        } else {
            result += negative ? " - " : " + ";
        }
        const Rational magnitude = term.second.abs();
        if (term.first.is_constant()) {
            result += magnitude.to_string();
        } else {
            if (magnitude != Rational(1)) {
                result += magnitude.is_integer() ? magnitude.to_string() :
                          "(" + magnitude.to_string() + ")";
            }
            result += term.first.to_string();
        }
    }
    return result.empty() ? "0" : result;
}

lmmc_real_t Irrational::terms_to_double(const Terms& terms) {
    lmmc_real_t result = 0.0;
    for (const auto& term : terms) result += term.second.to_double() * term.first.to_double();
    return result;
}

lmmc_real_t Irrational::to_double() const {
    const auto rational = as_rational_checked();
    if (rational) return rational.value().to_double();
    return terms_to_double(numerator_) / terms_to_double(denominator_);
}

std::string Irrational::to_string() const {
    const auto rational = as_rational_checked();
    if (rational) return rational.value().to_string();
    const std::string numerator = terms_to_string(numerator_);
    if (denominator_.size() == 1 && denominator_.begin()->first.is_constant() &&
        denominator_.begin()->second == Rational(1)) return numerator;
    return "(" + numerator + ")/(" + terms_to_string(denominator_) + ")";
}

std::ostream& operator<<(std::ostream& os, const Irrational& value) {
    return os << value.to_string();
}

}
