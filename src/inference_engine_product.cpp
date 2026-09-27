#define _USE_MATH_DEFINES
#include "internal/inference_engine_impl.hpp"
#include "internal/assumption_facts.hpp"
#include <cstdint>

namespace LMCAS {
namespace {
constexpr std::uint8_t negative_sign = 1;
constexpr std::uint8_t zero_sign = 2;
constexpr std::uint8_t positive_sign = 4;
constexpr std::uint8_t all_signs = 7;

std::uint8_t sign_set(Sign target) {
    switch (target) {
    case Sign::Positive: return positive_sign;
    case Sign::Negative: return negative_sign;
    case Sign::Zero: return zero_sign;
    case Sign::NonNegative: return zero_sign | positive_sign;
    case Sign::NonPositive: return negative_sign | zero_sign;
    case Sign::NonZero: return negative_sign | positive_sign;
    default: return all_signs;
    }
}

std::uint8_t multiply_sign_sets(std::uint8_t left, std::uint8_t right) {
    std::uint8_t result = 0;
    for (std::uint8_t a : {negative_sign, zero_sign, positive_sign}) {
        for (std::uint8_t b : {negative_sign, zero_sign, positive_sign}) {
            if (!(left & a) || !(right & b)) continue;
            result |= a == zero_sign || b == zero_sign ? zero_sign :
                a == b ? positive_sign : negative_sign;
        }
    }
    return result;
}
}

InferenceTriboolResult InferenceEngine::Impl::accumulate_product_sign(
    const InferenceEngine& engine,
    const std::shared_ptr<const SymbolicNode>& operand,
    std::uint8_t& possible, ComputationContext& context) const {
    const auto expression = detail::expression_from_node(operand);
    std::uint8_t operand_signs = all_signs;
    for (Sign sign : {Sign::Positive, Sign::Negative, Sign::Zero,
                      Sign::NonNegative, Sign::NonPositive, Sign::NonZero}) {
        auto fact = engine.query_sign_of_checked(expression, sign, context);
        if (!fact) return fact;
        if (fact.value() == Tribool::True) operand_signs &= sign_set(sign);
    }
    possible = multiply_sign_sets(possible, operand_signs);
    return Tribool::True;
}

InferenceTriboolResult InferenceEngine::Impl::product_parity_sign(
    std::uint8_t possible, Sign target) {
    if (!possible) return Tribool::Unknown;
    const auto wanted = sign_set(target);
    if ((possible & wanted) == possible) return Tribool::True;
    if (!(possible & wanted)) return Tribool::False;
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::product_sign(
    const InferenceEngine& engine, const MultiplyNode& node, Sign target, ComputationContext& context) const {
    auto defined = real_definedness(engine, node, context);
    if (!defined || defined.value() != Tribool::True) { return defined; }
    const detail::AssumptionFacts facts(engine);
    bool all_real = true;
    for (const auto& operand : node.operands()) {
        auto real = detail::query_real_value(operand, facts, context);
        if (!real) { return real; }
        all_real = all_real && real.value() == Tribool::True;
    }
    if (!all_real) { return Tribool::Unknown; }
    std::uint8_t possible = positive_sign;
    for (const auto& operand : node.operands()) {
        auto result = accumulate_product_sign(engine, operand, possible, context);
        if (!result) { return result; }
    }
    return product_parity_sign(possible, target);
}

InferenceTriboolResult InferenceEngine::Impl::quotient_sign(
    const InferenceEngine& engine, const SymbolicExpr& numerator,
    const SymbolicExpr& denominator, Sign target, ComputationContext& context) const {
    const detail::AssumptionFacts facts(engine);
    bool all_defined = true;
    bool all_real = true;
    for (const auto* expression : {&numerator, &denominator}) {
        auto defined = real_definedness(engine, *detail::node(*expression), context);
        if (!defined) { return defined; }
        all_defined = all_defined && defined.value() == Tribool::True;
        auto real = detail::query_real_value(detail::node(*expression), facts, context);
        if (!real) { return real; }
        all_real = all_real && real.value() == Tribool::True;
    }
    auto nonzero = detail::query_nonzero_value(detail::node(denominator), facts, Domain::Real, context);
    if (!nonzero) { return nonzero; }
    if (nonzero.value() == Tribool::False) {
        return InferenceTriboolResult::failure(
            CasErrc::DomainError, "division denominator is zero", "inference.division");
    }
    if (!all_defined || !all_real || nonzero.value() != Tribool::True)
        { return Tribool::Unknown; }
    std::uint8_t possible = positive_sign;
    auto first = accumulate_product_sign(engine, detail::node(numerator), possible, context);
    if (!first) { return first; }
    auto second = accumulate_product_sign(engine, detail::node(denominator), possible, context);
    if (!second) { return second; }
    return product_parity_sign(possible, target);
}

InferenceTriboolResult InferenceEngine::infer_division_sign_checked(
    const MultiplyNode& node, Sign target, ComputationContext& context) const {
    try {
        if (node.operands().size() != 2) {
            return Tribool::Unknown;
        }
        std::shared_ptr<const SymbolicNode> numerator;
        std::shared_ptr<const SymbolicNode> denominator;
        for (const auto& operand : node.operands()) {
            auto step = context.consume_steps(1, "inference.product");
            if (!step) {
                throw step.error();
            }
            const auto power = std::dynamic_pointer_cast<const PowerNode>(operand);
            const auto exponent = power
                ? std::dynamic_pointer_cast<const NumberNode>(
                      power->exponent())
                : nullptr;
            if (exponent && exponent->is_negative_one() && !denominator) {
                denominator = power->base();
            } else {
                numerator = operand;
            }
        }
        if (!numerator || !denominator) {
            return Tribool::Unknown;
        }
        return impl_->quotient_sign(*this, detail::expression_from_node(numerator),
            detail::expression_from_node(denominator), target, context);
    } catch (const CasError& error) {
        return InferenceTriboolResult::failure(error);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(CasErrc::ResourceLimit,
            "division sign inference allocation failed", "infer_division_sign_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(CasErrc::InternalInvariant,
            ex.what(), "infer_division_sign_checked");
    }
}

InferenceTriboolResult InferenceEngine::infer_multiply_sign_checked(
    const MultiplyNode& node, Sign target, ComputationContext& context) const {
    try {
        if (node.operands().size() == 2) {
            for (const auto& operand : node.operands()) {
                auto step = context.consume_steps(1, "inference.product");
                if (!step) throw step.error();
                const auto power =
                    std::dynamic_pointer_cast<const PowerNode>(operand);
                const auto exponent = power
                    ? std::dynamic_pointer_cast<const NumberNode>(
                          power->exponent())
                    : nullptr;
                if (exponent && exponent->is_negative_one())
                    return infer_division_sign_checked(node, target, context);
            }
        }
        return impl_->product_sign(*this, node, target, context);
    } catch (const CasError& error) {
        return InferenceTriboolResult::failure(error);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(CasErrc::ResourceLimit,
            "multiplication sign inference allocation failed", "infer_multiply_sign_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(CasErrc::InternalInvariant,
            ex.what(), "infer_multiply_sign_checked");
    }
}

}
