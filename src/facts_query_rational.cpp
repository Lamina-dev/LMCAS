#include "internal/facts_query.hpp"
#include "computation_context.hpp"
#include "internal/exact_constant_bounds.hpp"
#include "internal/normalization_utils.hpp"
#include "internal/symbolic_ast.hpp"

#include <algorithm>
#include <cmath>

namespace LMCAS::detail {
using Node = std::shared_ptr<const SymbolicNode>;

std::optional<Rational> closed_rational_value(
    const Node& node, Domain domain, ComputationContext& context);

namespace {

std::optional<Rational> closed_rational_power(
    const PowerNode& power, Domain domain, ComputationContext& context,
    ExactBoundArithmetic& arithmetic) {
    auto base = closed_rational_value(power.base(), domain, context);
    if (!base) { return std::nullopt; }
    auto exponent = closed_rational_value(power.exponent(), domain, context);
    if (!exponent || !exponent->is_integer()) { return std::nullopt; }
    BigInt count = exponent->to_bigint();
    if (base->is_zero() && (count.is_negative() ||
        (count.is_zero() && domain == Domain::Real))) { return std::nullopt; }
    if (count.is_negative()) {
        *base = arithmetic.div(Rational(1), *base);
        count = -count;
    }
    Rational value(1);
    while (!count.is_zero()) {
        arithmetic.step();
        if (count.is_odd()) { value = arithmetic.mul(value, *base); }
        count = count >> 1;
        if (!count.is_zero()) { *base = arithmetic.mul(*base, *base); }
    }
    return value;
}

}

// Fold closed rational arithmetic without normalization or fact queries.
std::optional<Rational> closed_rational_value(
    const std::shared_ptr<const SymbolicNode>& node, Domain domain, ComputationContext& context) {
    ExactBoundArithmetic arithmetic(context, "facts.constant_exponent");
    arithmetic.require(context.enter_recursion("facts.constant_exponent"));
    struct Leave {
        ComputationContext& context;
        ~Leave() { context.leave_recursion(); }
    } leave{context};
    if (const auto* number = dynamic_cast<const NumberNode*>(node.get())) {
        if (const auto* value = std::get_if<double>(&number->value());
            value && !std::isfinite(*value)) { return std::nullopt; }
        auto value = exact_number_as_rational(*number);
        arithmetic.bits(std::max(value.get_numerator().bit_length(),
                                 value.get_denominator().bit_length()));
        return value;
    }
    const auto* sum = dynamic_cast<const AddNode*>(node.get());
    const auto* product = dynamic_cast<const MultiplyNode*>(node.get());
    if (sum || product) {
        Rational value(product ? 1 : 0);
        for (const auto& operand : sum ? sum->operands() : product->operands()) {
            auto child = closed_rational_value(operand, domain, context);
            if (!child) { return std::nullopt; }
            value = product ? arithmetic.mul(value, *child) : arithmetic.add(value, *child);
        }
        return value;
    }
    const auto* power = dynamic_cast<const PowerNode*>(node.get());
    if (!power) { return std::nullopt; }
    return closed_rational_power(*power, domain, context, arithmetic);
}

}
