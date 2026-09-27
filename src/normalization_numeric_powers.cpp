#include "internal/normalization_utils.hpp"
#include <cstdint>

namespace LMCAS {
namespace {
std::shared_ptr<const NumberNode> reciprocal_exact_number(const NumberNode& number) {
    if (const auto* integer = std::get_if<BigInt>(&number.value())) {
        return detail::make_node<NumberNode>(Rational(BigInt(1), *integer));
    }
    if (const auto* rational = std::get_if<Rational>(&number.value())) {
        if (!rational->get_numerator().is_zero()) {
            return detail::make_node<NumberNode>(
                Rational(rational->get_denominator(), rational->get_numerator()));
        }
    }
    return nullptr;
}

std::shared_ptr<const NumberNode> small_exact_power(const NumberNode& base,
                                                   const BigInt& exponent) {
    const auto count = exponent.abs().try_to_uint64();
    if (const auto* integer = std::get_if<BigInt>(&base.value())) {
        BigInt value(1);
        for (std::uint64_t k = 0; k < *count; ++k) value = value * *integer;
        return detail::make_node<NumberNode>(value);
    }
    if (const auto* rational = std::get_if<Rational>(&base.value())) {
        Rational value(1);
        for (std::uint64_t k = 0; k < *count; ++k) value = value * *rational;
        return detail::make_node<NumberNode>(value);
    }
    return nullptr;
}

std::shared_ptr<const SymbolicNode> integer_numeric_power(
    const std::shared_ptr<const NumberNode>& base,
    const std::shared_ptr<const NumberNode>& exponent_node, const BigInt& exponent,
    detail::RewriteBudget* budget) {
    const bool small = exponent > BigInt(-64) && exponent < BigInt(64);
    const auto* real = std::get_if<lmmc_real_t>(&base->value());
    if (real && small) {
        const double value = std::pow(*real, exponent.to_double());
        if (std::isfinite(value)) {
            return detail::make_node<NumberNode>(value);
        }
        normalization_check_children(budget, 1, base, exponent_node);
        return detail::make_node<PowerNode>(base, exponent_node);
    }
    if (exponent == BigInt(-1)) {
        return reciprocal_exact_number(*base);
    }
    if (exponent.is_zero()) {
        return detail::make_node<NumberNode>(BigInt(1));
    }
    if (!small) {
        return nullptr;
    }
    auto powered = small_exact_power(*base, exponent);
    if (!powered) {
        return nullptr;
    }
    if (exponent > BigInt(0)) {
        return powered;
    }
    return reciprocal_exact_number(*powered);
}

bool is_half_exponent(const NumberNode& exponent) {
    if (const auto* real = std::get_if<lmmc_real_t>(&exponent.value())) {
        return *real == 0.5;
    }
    if (const auto* rational = std::get_if<Rational>(&exponent.value())) {
        return *rational == Rational(1, 2);
    }
    return false;
}

std::shared_ptr<const SymbolicNode> extract_small_square_factor(
    const BigInt& value, const BigInt& root, detail::RewriteBudget* budget) {
    if (value <= BigInt(0) || value >= BigInt(1000000)) {
        return nullptr;
    }
    for (BigInt factor = root; factor >= BigInt(2); factor -= BigInt(1)) {
        const BigInt square = factor * factor;
        if (!(value % square).is_zero()) {
            continue;
        }
        if (budget) { budget->require_nodes(5); }
        auto coefficient = detail::make_node<NumberNode>(factor);
        auto radicand = detail::make_node<NumberNode>(value / square);
        auto half = detail::make_node<NumberNode>(Rational(1, 2));
        auto radical = detail::make_node<PowerNode>(radicand, half);
        return SymbolicFactory::create_multiply({coefficient, radical}, budget);
    }
    return nullptr;
}

std::shared_ptr<const SymbolicNode> half_numeric_power(
    const NumberNode& base, detail::RewriteBudget* budget) {
    if (const auto* integer = std::get_if<BigInt>(&base.value())) {
        if (*integer < BigInt(0)) {
            return nullptr;
        }
        const BigInt root = integer->sqrt();
        if (root * root == *integer) {
            return detail::make_node<NumberNode>(root);
        }
        return extract_small_square_factor(*integer, root, budget);
    }
    if (auto root = normalization_exact_square_root(base)) {
        return root;
    }
    if (const auto* real = std::get_if<lmmc_real_t>(&base.value())) {
        if (*real >= 0) {
            return detail::make_node<NumberNode>(std::sqrt(*real));
        }
    }
    return nullptr;
}
}

std::shared_ptr<const NumberNode> normalization_exact_square_root(const NumberNode& number) {
    if (const auto* integer = std::get_if<BigInt>(&number.value())) {
        if (*integer < BigInt(0)) {
            return nullptr;
        }
        const BigInt root = integer->sqrt();
        if (root * root == *integer) {
            return detail::make_node<NumberNode>(root);
        }
        return nullptr;
    }
    const auto* rational = std::get_if<Rational>(&number.value());
    if (!rational || *rational < Rational(0)) {
        return nullptr;
    }
    const BigInt numerator = rational->get_numerator().sqrt();
    const BigInt denominator = rational->get_denominator().sqrt();
    if (numerator * numerator != rational->get_numerator()) {
        return nullptr;
    }
    if (denominator * denominator != rational->get_denominator()) {
        return nullptr;
    }
    return detail::make_node<NumberNode>(Rational(numerator, denominator));
}

std::shared_ptr<const SymbolicNode> normalization_imaginary_power(
    const BigInt& exponent, detail::RewriteBudget* budget) {
    BigInt remainder = exponent % BigInt(4);
    if (remainder.is_negative()) remainder += BigInt(4);
    if (budget) { budget->require_nodes(3); }
    if (remainder.is_zero()) {
        return detail::make_node<NumberNode>(BigInt(1));
    }
    if (remainder == BigInt(2)) {
        return detail::make_node<NumberNode>(BigInt(-1));
    }
    auto imaginary = detail::make_node<PowerNode>(
        detail::make_node<NumberNode>(BigInt(-1)),
        detail::make_node<NumberNode>(Rational(1, 2)));
    if (remainder == BigInt(1)) {
        return imaginary;
    }
    return make_normalized_multiply_node(
        {detail::make_node<NumberNode>(BigInt(-1)), imaginary}, budget);
}

std::shared_ptr<const SymbolicNode> normalization_numeric_power(
    const std::shared_ptr<const NumberNode>& base,
    const std::shared_ptr<const NumberNode>& exponent, detail::RewriteBudget* budget) {
    if (budget) { budget->require_nodes(1); }
    BigInt integer;
    if (try_get_integer_value(exponent, integer)) {
        if (auto value = integer_numeric_power(base, exponent, integer, budget)) {
            return value;
        }
    } else if (is_half_exponent(*exponent)) {
        if (auto value = half_numeric_power(*base, budget)) {
            return value;
        }
    }
    const auto* rational = std::get_if<Rational>(&exponent->value());
    if (!rational || !base->is_negative_one()) {
        return nullptr;
    }
    if (rational->get_denominator() != BigInt(2)) {
        return nullptr;
    }
    return normalization_imaginary_power(rational->get_numerator(), budget);
}
}
