#include "irrational.hpp"
#include "internal/exact_constant_bounds.hpp"
#include <limits>
#include <new>
#include <optional>
#include <stdexcept>

namespace LMCAS {

void Irrational::add_term(Terms& terms, const Basis& basis, const Rational& coefficient) {
    if (coefficient.is_zero()) return;
    auto inserted = terms.emplace(basis, coefficient);
    if (inserted.second) return;
    auto& value = inserted.first->second;
    value = value + coefficient;
    if (value.is_zero()) terms.erase(inserted.first);
}

Irrational::Terms Irrational::multiply_terms(const Terms& left, const Terms& right) {
    Terms result;
    for (const auto& a : left) {
        for (const auto& b : right) {
            const BigInt common = BigInt::gcd(a.first.radicand, b.first.radicand);
            Basis basis;
            basis.radicand = (a.first.radicand / common) * (b.first.radicand / common);
            basis.pi_power = a.first.pi_power + b.first.pi_power;
            basis.e_power = a.first.e_power + b.first.e_power;
            add_term(result, basis, a.second * b.second * Rational(common));
        }
    }
    return result;
}

Irrational Irrational::operator+(const Irrational& other) const {
    Irrational result;
    if (denominator_ == other.denominator_) {
        result.numerator_ = numerator_;
        result.denominator_ = denominator_;
        for (const auto& term : other.numerator_) {
            add_term(result.numerator_, term.first, term.second);
        }
    } else {
        result.numerator_ = multiply_terms(numerator_, other.denominator_);
        const Terms right = multiply_terms(other.numerator_, denominator_);
        for (const auto& term : right) add_term(result.numerator_, term.first, term.second);
        result.denominator_ = multiply_terms(denominator_, other.denominator_);
    }
    result.simplify();
    return result;
}

Irrational Irrational::operator-(const Irrational& other) const {
    return *this + (-other);
}

Irrational Irrational::operator*(const Irrational& other) const {
    Irrational result;
    result.numerator_ = multiply_terms(numerator_, other.numerator_);
    result.denominator_ = multiply_terms(denominator_, other.denominator_);
    result.simplify();
    return result;
}

Irrational Irrational::operator*(const Rational& scalar) const {
    Irrational result = *this;
    if (scalar.is_zero()) {
        result.numerator_.clear();
    } else {
        for (auto& term : result.numerator_) term.second = term.second * scalar;
    }
    return result;
}

Irrational Irrational::operator*(lmmc_real_t scalar) const {
    return *this * Rational::from_double(scalar);
}

Irrational Irrational::operator/(const Irrational& other) const {
    if (other.is_zero()) throw std::domain_error("Irrational: division by zero");
    Irrational result;
    result.numerator_ = multiply_terms(numerator_, other.denominator_);
    result.denominator_ = multiply_terms(denominator_, other.numerator_);
    result.simplify();
    return result;
}

Irrational Irrational::operator-() const {
    return *this * Rational(-1);
}

Irrational Irrational::pow(const BigInt& exponent) const {
    BigInt remaining = exponent.abs();
    Irrational result = constant(Rational(1));
    Irrational base = *this;
    while (!remaining.is_zero()) {
        if (remaining.is_odd()) result = result * base;
        remaining >>= 1;
        if (!remaining.is_zero()) base = base * base;
    }
    return exponent.is_negative() ? constant(Rational(1)) / result : result;
}

bool Irrational::operator==(const Irrational& other) const {
    if (denominator_ == other.denominator_) return numerator_ == other.numerator_;
    return multiply_terms(numerator_, other.denominator_) ==
           multiply_terms(other.numerator_, denominator_);
}

bool Irrational::operator<(const Irrational& other) const {
    auto left = as_rational_checked();
    auto right = other.as_rational_checked();
    if (left && right) return left.value() < right.value();
    if (*this == other) return false;
    return to_double() < other.to_double();
}

bool Irrational::operator<=(const Irrational& other) const {
    return *this < other || *this == other;
}

bool Irrational::operator>(const Irrational& other) const {
    return other < *this;
}

bool Irrational::operator>=(const Irrational& other) const {
    return other <= *this;
}

bool Irrational::is_positive() const {
    return constant(Rational(0)) < *this;
}

bool Irrational::is_negative() const {
    return *this < constant(Rational(0));
}

struct Irrational::EnclosureRound {
    using Bounds = detail::ExactConstantInterval;
    ComputationContext& context;
    std::size_t precision;
    detail::ExactBoundArithmetic arithmetic;
    std::map<BigInt, Bounds> radicals;
    std::optional<Bounds> pi;
    std::optional<Bounds> e;

    EnclosureRound(ComputationContext& context, std::size_t precision)
        : context(context), precision(precision),
          arithmetic(context, "Irrational::sign_checked") {}

    const Bounds& radical(const BigInt& radicand) {
        auto found = radicals.find(radicand);
        if (found != radicals.end()) return found->second;
        arithmetic.step();
        if (radicand.is_negative() || radicand.is_zero())
            throw CasError{CasErrc::InternalInvariant,
                           "nonpositive normalized radicand", "Irrational::sign_checked"};
        if (radicand == BigInt(1))
            return radicals.emplace(radicand, Bounds{Rational(1), Rational(1)}).first->second;
        const BigInt scale = arithmetic.dyadic_scale(precision);
        const auto shift = arithmetic.growth(precision, precision);
        arithmetic.growth(radicand.bit_length(), shift);
        if (shift > std::numeric_limits<mp_size_t>::max())
            throw CasError{CasErrc::ResourceLimit,
                           "radical shift size exceeded", "Irrational::sign_checked"};
        const BigInt scaled = radicand << static_cast<mp_size_t>(shift);
        arithmetic.step();
        const BigInt lower = scaled.sqrt();
        arithmetic.growth(lower.bit_length(), lower.bit_length());
        const bool exact = lower * lower == scaled;
        arithmetic.growth(lower.bit_length(), 1);
        const BigInt upper = exact ? lower : lower + BigInt(1);
        return radicals.emplace(radicand, Bounds{Rational(lower, scale),
                                                  Rational(upper, scale)}).first->second;
    }

    const Bounds& constant(bool use_pi) {
        auto& cached = use_pi ? pi : e;
        if (!cached) {
            auto bounds = use_pi ? detail::exact_pi_bounds(precision, context)
                                 : detail::exact_e_bounds(precision, context);
            if (!bounds) throw bounds.error();
            cached = std::move(bounds.value());
        }
        return *cached;
    }

    Bounds multiply_positive(const Bounds& a, const Bounds& b) {
        return {arithmetic.mul(a.first, b.first), arithmetic.mul(a.second, b.second)};
    }

    Bounds power(Bounds base, const BigInt& exponent) {
        arithmetic.bits(exponent.bit_length());
        BigInt remaining = exponent.abs();
        Bounds result{Rational(1), Rational(1)};
        while (!remaining.is_zero()) {
            arithmetic.step();
            if (remaining.is_odd()) result = multiply_positive(result, base);
            remaining >>= 1;
            if (!remaining.is_zero()) base = multiply_positive(base, base);
        }
        if (exponent.is_negative())
            return {arithmetic.div(Rational(1), result.second),
                    arithmetic.div(Rational(1), result.first)};
        return result;
    }

    static bool validate_terms(const Terms& terms,
        detail::ExactBoundArithmetic& arithmetic, bool& transcendental) {
        bool zero = true;
        for (const auto& term : terms) {
            arithmetic.step();
            arithmetic.bits(term.second.get_numerator().bit_length());
            arithmetic.bits(term.second.get_denominator().bit_length());
            arithmetic.bits(term.first.radicand.bit_length());
            arithmetic.bits(term.first.pi_power.bit_length());
            arithmetic.bits(term.first.e_power.bit_length());
            if (term.second.is_zero()) { continue; }
            zero = false;
            transcendental |= !term.first.pi_power.is_zero() ||
                              !term.first.e_power.is_zero();
        }
        return zero;
    }

    static int separated_sign(const Bounds& bounds) {
        const auto lower = bounds.first.get_numerator();
        if (!lower.is_zero() && !lower.is_negative()) { return 1; }
        if (bounds.second.get_numerator().is_negative()) { return -1; }
        return 0;
    }

    static Result<Sign> refine_sign(const Terms& numerator, const Terms& denominator,
        bool numerator_zero, bool transcendental, ComputationContext& context,
        detail::ExactBoundArithmetic& arithmetic, const char* operation) {
        /**
         * @brief 不同规范无平方因子根式在 Q 上线性无关，故非空纯根式和非零，细化过程在数学上终止；
         * 实际执行仍受调用方资源预算限制。
         */
        for (std::size_t precision = 32;;) {
            arithmetic.step();
            EnclosureRound round(context, precision);
            const auto denominator_bounds = terms_bounds(denominator, round);
            const int denominator_sign = separated_sign(denominator_bounds);
            if (denominator_sign != 0) {
                if (numerator_zero) { return Sign::Zero; }
                const auto numerator_bounds = terms_bounds(numerator, round);
                const int numerator_sign = separated_sign(numerator_bounds);
                if (numerator_sign != 0) {
                    return numerator_sign == denominator_sign ? Sign::Positive : Sign::Negative;
                }
            }
            if (transcendental && precision == 4096) {
                return Result<Sign>::failure(CasErrc::Inconclusive,
                    "pi/e enclosure did not separate zero at 4096 bits", operation);
            }
            if (precision > std::numeric_limits<std::size_t>::max() / 2) {
                return Result<Sign>::failure(CasErrc::ResourceLimit,
                                             "refinement size exceeded", operation);
            }
            precision *= 2;
        }
    }
};

std::pair<Rational, Rational> Irrational::terms_bounds(
    const Terms& terms, EnclosureRound& round) {
    auto& arithmetic = round.arithmetic;
    EnclosureRound::Bounds sum{Rational(0), Rational(0)};
    for (const auto& term : terms) {
        arithmetic.step();
        if (term.second.is_zero()) continue;
        auto basis = round.radical(term.first.radicand);
        if (!term.first.pi_power.is_zero())
            basis = round.multiply_positive(
                basis, round.power(round.constant(true), term.first.pi_power));
        if (!term.first.e_power.is_zero())
            basis = round.multiply_positive(
                basis, round.power(round.constant(false), term.first.e_power));
        if (term.second.get_numerator().is_negative()) std::swap(basis.first, basis.second);
        sum.first = arithmetic.add(sum.first, arithmetic.mul(term.second, basis.first));
        sum.second = arithmetic.add(sum.second, arithmetic.mul(term.second, basis.second));
    }
    return sum;
}

Result<Sign> Irrational::sign_checked(ComputationContext& context) const {
    const char* operation = "Irrational::sign_checked";
    try {
        detail::ExactBoundArithmetic arithmetic(context, operation);
        arithmetic.step();
        arithmetic.require(context.require_expansion_terms(numerator_.size(), operation));
        arithmetic.require(context.require_expansion_terms(denominator_.size(), operation));
        bool transcendental = false;
        const bool numerator_zero = EnclosureRound::validate_terms(
            numerator_, arithmetic, transcendental);
        const bool denominator_zero = EnclosureRound::validate_terms(
            denominator_, arithmetic, transcendental);
        if (denominator_zero) {
            return Result<Sign>::failure(CasErrc::DomainError, "zero denominator", operation);
        }
        return EnclosureRound::refine_sign(numerator_, denominator_, numerator_zero,
            transcendental, context, arithmetic, operation);
    } catch (const CasError& error) {
        return Result<Sign>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<Sign>::failure(CasErrc::ResourceLimit, "sign enclosure allocation failed", operation);
    } catch (const std::length_error&) {
        return Result<Sign>::failure(CasErrc::ResourceLimit, "sign enclosure size exceeded", operation);
    }
}

Result<Sign> Irrational::sign_checked() const {
    ComputationContext context;
    return sign_checked(context);
}

Result<Irrational> Irrational::abs_checked(ComputationContext& context) const {
    auto sign = sign_checked(context);
    if (!sign) return Result<Irrational>::failure(sign.error());
    try {
        if (sign.value() == Sign::Zero) return constant(Rational(0));
        Irrational result = *this;
        if (sign.value() == Sign::Negative) {
            for (auto& term : result.numerator_) {
                auto step = context.consume_steps(1, "Irrational::abs_checked");
                if (!step) return Result<Irrational>::failure(step.error());
                term.second = -term.second; /**< 取负保持分子、分母整数的位数不变。 */
            }
        }
        return result;
    } catch (const std::bad_alloc&) {
        return Result<Irrational>::failure(CasErrc::ResourceLimit,
            "absolute value allocation failed", "Irrational::abs_checked");
    } catch (const std::length_error&) {
        return Result<Irrational>::failure(CasErrc::ResourceLimit,
            "absolute value size exceeded", "Irrational::abs_checked");
    }
}

Result<Irrational> Irrational::abs_checked() const {
    ComputationContext context;
    return abs_checked(context);
}

Irrational Irrational::abs() const {
    auto result = abs_checked();
    if (result) return std::move(result.value());
    const auto& error = result.error();
    const std::string message = error.operation + ": " + error.message;
    switch (error.code) {
        case CasErrc::DomainError: throw std::domain_error(message);
        case CasErrc::ResourceLimit: throw std::length_error(message);
        case CasErrc::InternalInvariant: throw std::logic_error(message);
        default: throw std::runtime_error(message);
    }
}

}
