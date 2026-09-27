#include "irrational.hpp"
#include <cmath>
#include <new>
#include <utility>

namespace LMCAS {

bool Irrational::Basis::operator<(const Basis& other) const {
    if (radicand != other.radicand) return radicand < other.radicand;
    if (pi_power != other.pi_power) return pi_power < other.pi_power;
    return e_power < other.e_power;
}

bool Irrational::Basis::operator==(const Basis& other) const {
    return radicand == other.radicand && pi_power == other.pi_power &&
           e_power == other.e_power;
}

bool Irrational::Basis::is_constant() const {
    return radicand == BigInt(1) && pi_power.is_zero() && e_power.is_zero();
}

Result<std::pair<BigInt, BigInt>> Irrational::extract_square_factors(
    const BigInt& n, ComputationContext& context) {
    using Factors = std::pair<BigInt, BigInt>;
    if (n.is_zero()) return Factors{BigInt(0), BigInt(1)};
    const BigInt root = n.sqrt();
    if (root * root == n) return Factors{root, BigInt(1)};
    BigInt perfect(1);
    BigInt remainder = n;
    const BigInt one(1);
    for (BigInt divisor(2); divisor <= remainder / divisor; divisor += one) {
        auto step = context.consume_steps(1, "Irrational::sqrt_checked");
        if (!step) return Result<Factors>::failure(step.error());
        const BigInt square = divisor * divisor;
        while ((remainder % square).is_zero()) {
            step = context.consume_steps(1, "Irrational::sqrt_checked");
            if (!step) return Result<Factors>::failure(step.error());
            perfect *= divisor;
            remainder /= square;
        }
    }
    return Factors{std::move(perfect), std::move(remainder)};
}

Result<Irrational> Irrational::sqrt_checked(const BigInt& n, lmmc_real_t coefficient,
                                           ComputationContext& context) {
    const std::string operation = "Irrational::sqrt_checked";
    auto step = context.consume_steps(1, operation);
    if (!step) {
        return Result<Irrational>::failure(step.error());
    }
    if (n.is_negative()) {
        return Result<Irrational>::failure(CasErrc::DomainError,
                                           "negative square root radicand", operation);
    }
    if (!std::isfinite(coefficient)) {
        return Result<Irrational>::failure(CasErrc::InvalidArgument,
                                           "square root coefficient must be finite", operation);
    }
    auto bits = context.require_integer_bits(n.bit_length(), operation);
    if (!bits) {
        return Result<Irrational>::failure(bits.error());
    }
    try {
        Rational exact = Rational::from_double(coefficient);
        bits = context.require_integer_bits(exact.get_numerator().bit_length(), operation);
        if (!bits) {
            return Result<Irrational>::failure(bits.error());
        }
        bits = context.require_integer_bits(exact.get_denominator().bit_length(), operation);
        if (!bits) {
            return Result<Irrational>::failure(bits.error());
        }
        Irrational result;
        result.type_ = Type::SQRT;
        if (exact.is_zero() || n.is_zero()) {
            return result;
        }
        auto factors = extract_square_factors(n, context);
        if (!factors) {
            return Result<Irrational>::failure(factors.error());
        }
        exact = exact * Rational(factors.value().first);
        bits = context.require_integer_bits(exact.get_numerator().bit_length(), operation);
        if (!bits) {
            return Result<Irrational>::failure(bits.error());
        }
        Basis basis;
        basis.radicand = std::move(factors.value().second);
        result.numerator_.emplace(std::move(basis), std::move(exact));
        return result;
    } catch (const std::bad_alloc&) {
        return Result<Irrational>::failure(CasErrc::ResourceLimit,
                                           "square root allocation failed", operation);
    }
}

Result<Irrational> Irrational::sqrt_checked(const BigInt& n, lmmc_real_t coefficient) {
    ComputationContext context;
    return sqrt_checked(n, coefficient, context);
}

Irrational Irrational::pi(lmmc_real_t coefficient) {
    Irrational result;
    result.type_ = Type::PI;
    Basis basis;
    basis.pi_power = BigInt(1);
    add_term(result.numerator_, basis, Rational::from_double(coefficient));
    return result;
}

Irrational Irrational::e(lmmc_real_t coefficient) {
    Irrational result;
    result.type_ = Type::E;
    Basis basis;
    basis.e_power = BigInt(1);
    add_term(result.numerator_, basis, Rational::from_double(coefficient));
    return result;
}

Irrational Irrational::constant(lmmc_real_t value) {
    return constant(Rational::from_double(value));
}

Irrational Irrational::constant(const Rational& value) {
    Irrational result;
    add_term(result.numerator_, Basis{}, value);
    return result;
}

bool Irrational::is_zero() const {
    return numerator_.empty();
}

Result<Rational> Irrational::as_rational_checked() const {
    if (is_zero()) return Rational(0);
    if (numerator_.size() == denominator_.size()) {
        const Rational ratio = numerator_.begin()->second / denominator_.begin()->second;
        auto numerator = numerator_.begin();
        auto denominator = denominator_.begin();
        for (; numerator != numerator_.end(); ++numerator, ++denominator) {
            if (!(numerator->first == denominator->first) ||
                numerator->second != denominator->second * ratio) break;
        }
        if (numerator == numerator_.end()) return ratio;
    }
    return Result<Rational>::failure(CasErrc::InvalidArgument,
                                     "irrational value is not rational",
                                     "Irrational::as_rational_checked");
}

bool Irrational::is_rational() const {
    return as_rational_checked().has_value();
}

void Irrational::simplify() {
    /**
     * @brief 纯根式非空和由 Q 线性无关性保证非零；规范单项式因各基值为正而非零。
     * 含 pi/e 的多项式不采用该线性无关性结论；即使分子为零或与分母成比例，
     * 也保留分母，由 sign_checked 证明其非零后再返回值。
     */
    if (denominator_.size() > 1) {
        for (const auto& term : denominator_) {
            if (!term.first.pi_power.is_zero() || !term.first.e_power.is_zero())
                return;
        }
    }
    auto rational = as_rational_checked();
    if (rational) {
        numerator_.clear();
        denominator_.clear();
        denominator_.emplace(Basis{}, Rational(1));
        add_term(numerator_, Basis{}, rational.value());
    } else if (denominator_.size() == 1 && denominator_.begin()->first.is_constant()) {
        const Rational divisor = denominator_.begin()->second;
        for (auto& term : numerator_) term.second = term.second / divisor;
        denominator_.begin()->second = Rational(1);
    }
}

}
