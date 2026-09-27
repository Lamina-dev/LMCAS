#include "rational.hpp"
#include <stdexcept>

namespace LMCAS {

void Rational::simplify() {
    if (!denominator) {
        throw std::runtime_error("Denominator cannot be zero");
    }

    if (denominator.is_negative()) {
        numerator = -numerator;
        denominator = denominator.abs();
    }

    BigInt g = BigInt::gcd(numerator, denominator);
    if (g != 0 && g != 1) {
        numerator = numerator / g;
        denominator = denominator / g;
    }
}

Rational::Rational(const BigInt& num, const BigInt& den) : numerator(num), denominator(den) {
    simplify();
}

Rational::Rational(int num, int den) : numerator(num), denominator(den) {
    simplify();
}

Rational Rational::operator+(const Rational& other) const {
    BigInt new_num = numerator * other.denominator + other.numerator * denominator;
    BigInt new_den = denominator * other.denominator;
    return Rational(new_num, new_den);
}

Rational Rational::operator-(const Rational& other) const {
    BigInt new_num = numerator * other.denominator - other.numerator * denominator;
    BigInt new_den = denominator * other.denominator;
    return Rational(new_num, new_den);
}

Rational Rational::operator*(const Rational& other) const {
    BigInt new_num = numerator * other.numerator;
    BigInt new_den = denominator * other.denominator;
    return Rational(new_num, new_den);
}

Rational Rational::operator/(const Rational& other) const {
    if (other.is_zero()) {
        throw std::runtime_error("Division by zero");
    }
    BigInt new_num = numerator * other.denominator;
    BigInt new_den = denominator * other.numerator;
    return Rational(new_num, new_den);
}

Rational Rational::power(const BigInt& exponent) const {
    if (exponent < BigInt(0)) {

        if (is_zero()) {
            throw std::runtime_error("Cannot raise zero to negative power");
        }
        BigInt pos_exp = -exponent;
        return Rational(denominator.power(pos_exp), numerator.power(pos_exp));
    }

    if (!exponent) {
        if (is_zero()) {
            throw std::runtime_error("0^0 is undefined");
        }
        return Rational(1);
    }

    return Rational(numerator.power(exponent), denominator.power(exponent));
}

bool Rational::operator==(const Rational& other) const {
    return numerator * other.denominator == other.numerator * denominator;
}

bool Rational::operator!=(const Rational& other) const {
    return !(*this == other);
}

bool Rational::operator<(const Rational& other) const {
    BigInt left = numerator * other.denominator;
    BigInt right = other.numerator * denominator;
    return left < right;
}

bool Rational::operator<=(const Rational& other) const {
    return *this < other || *this == other;
}

bool Rational::operator>(const Rational& other) const {
    return !(*this <= other);
}

bool Rational::operator>=(const Rational& other) const {
    return !(*this < other);
}

Rational Rational::reciprocal() const {
    if (is_zero()) {
        throw std::runtime_error("Cannot take reciprocal of zero");
    }
    return Rational(denominator, numerator);
}

Rational Rational::abs() const {
    Rational result = *this;
    result.numerator = result.numerator.abs();
    return result;
}

Rational Rational::operator-() const {
    Rational result = *this;
    result.numerator = -result.numerator;
    return result;
}

}
