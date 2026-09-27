#include "rational.hpp"
#include <limits>
#include <stdexcept>
#include <utility>

namespace LMCAS {

bool Rational::power_leq(BigInt base, unsigned long exponent,
                      const BigInt& limit) {
    BigInt result(1);
    while (exponent != 0) {
        if ((exponent & 1UL) != 0) {
            result *= base;
            if (result > limit) return false;
        }
        exponent >>= 1UL;
        if (exponent != 0) {
            base *= base;
            if (base > limit) base = limit + BigInt(1);
        }
    }
    return result <= limit;
}

BigInt Rational::integer_nth_root(const BigInt& value, unsigned long degree) {
    if (degree == 0) {
        throw std::domain_error("Zeroth root is undefined");
    }
    if (value < BigInt(0)) {
        throw std::domain_error("Root input must be non-negative");
    }
    if (value <= BigInt(1) || degree == 1) {
        return value;
    }

    BigInt low(0);
    BigInt high(1);
    while (power_leq(high, degree, value)) {
        high *= BigInt(2);
    }

    while (high - low > BigInt(1)) {
        BigInt mid = (low + high) / BigInt(2);
        if (power_leq(mid, degree, value)) {
            low = std::move(mid);
        } else {
            high = std::move(mid);
        }
    }
    return low;
}

void Rational::sqrt_self(std::int64_t n) {
    if (n < 0) throw std::invalid_argument("Square-root precision must be non-negative");
    if (numerator < BigInt(0)) throw std::domain_error("Square root of a negative Rational");
    if (numerator == BigInt(0)) {
        denominator = BigInt(1);
        return;
    }
    if (n > static_cast<std::int64_t>(max_decimal_scale / 2)) {
        throw std::length_error("Square-root precision exceeds safety limit");
    }

    BigInt numerator_root = numerator.sqrt();
    BigInt denominator_root = denominator.sqrt();
    if (numerator_root * numerator_root == numerator &&
        denominator_root * denominator_root == denominator) {
        numerator = std::move(numerator_root);
        denominator = std::move(denominator_root);
        simplify();
        return;
    }

    const BigInt scale = pow10(static_cast<std::size_t>(n));
    const BigInt scaled = numerator * scale * scale / denominator;
    numerator = scaled.sqrt();
    denominator = scale;
    simplify();
}

Rational Rational::sqrt(std::int64_t n) const{
    Rational re = *this;
    re.sqrt_self(n);
    return re;
}

void Rational::radicand_self(const BigInt& radical, std::int64_t n) {
    if (radical <= BigInt(0)) throw std::domain_error("Root degree must be positive");
    if (radical > BigInt(1024)) throw std::length_error("Root degree exceeds safety limit");
    if (n < 0) throw std::invalid_argument("Root precision must be non-negative");

    const unsigned long degree = static_cast<unsigned long>(*radical.try_to_uint64());
    const bool negative = numerator < BigInt(0);
    if (negative && (degree % 2UL) == 0) {
        throw std::domain_error("Even root of a negative Rational");
    }
    if (n > static_cast<std::int64_t>(max_decimal_scale / degree)) {
        throw std::length_error("Root precision exceeds safety limit");
    }

    const BigInt scale = pow10(static_cast<std::size_t>(n));
    const BigInt scaled = numerator.abs() * scale.power(degree) / denominator;
    numerator = integer_nth_root(scaled, degree);
    if (negative) numerator = -numerator;
    denominator = scale;
    simplify();
}

Rational Rational::radicand(const BigInt& radical, std::int64_t n) {
    Rational re = *this;
    re.radicand_self(radical, n);
    return re;
}

void Rational::power_self(const Rational& exponent, std::size_t n) {
    Rational result = power(exponent.numerator);
    if (n > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
        throw std::length_error("Root precision exceeds safety limit");
    }
    result.radicand_self(exponent.denominator, static_cast<std::int64_t>(n));
    *this = std::move(result);
}

Rational Rational::power(const Rational& exponent, std::size_t n) {
    Rational re = *this;
    re.power_self(exponent, n);
    return re;
}

}
