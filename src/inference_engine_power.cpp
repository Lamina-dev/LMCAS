#define _USE_MATH_DEFINES
#include "internal/inference_engine_impl.hpp"
#include "internal/assumption_facts.hpp"

namespace LMCAS {

InferenceTriboolResult InferenceEngine::Impl::zero_power_sign(
    const InferenceEngine& engine,
    const SymbolicExpr& base_expr,
    const NumberNode* exp_num,
    Sign target, ComputationContext& context) const {
    if (!exp_num || !is_zero_number(*exp_num)) {
        return Tribool::Unknown;
    }
    const detail::AssumptionFacts facts(engine);
    auto base_nonzero = detail::query_nonzero_value(detail::node(base_expr), facts, Domain::Real, context);
    if (!base_nonzero) return base_nonzero;
    if (base_nonzero.value() == Tribool::True) {
        return known_sign(1, target);
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::positive_power_sign(
    const InferenceEngine& engine,
    const SymbolicExpr& base_expr,
    const SymbolicExpr& exp_expr,
    Sign target, ComputationContext& context) const {
    const detail::AssumptionFacts facts(engine);
    auto base_pos = detail::query_positive_value(detail::node(base_expr), facts, context);
    if (!base_pos) {
        return base_pos;
    }
    if (base_pos.value() != Tribool::True) {
        return Tribool::Unknown;
    }
    auto exp_real = detail::query_real_value(detail::node(exp_expr), facts, context);
    if (!exp_real) {
        return exp_real;
    }
    if (exp_real.value() == Tribool::True) {
        return known_sign(1, target);
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::even_power_sign(
    const InferenceEngine& engine,
    const SymbolicExpr& base_expr,
    const NumberNode& exponent,
    Sign target, ComputationContext& context) const {
    if (!is_even_integer_number(exponent)) {
        return Tribool::Unknown;
    }
    const detail::AssumptionFacts facts(engine);
    auto base = detail::query_real_value(detail::node(base_expr), facts, context);
    if (!base) {
        return base;
    }
    if (base.value() == Tribool::True) {
        auto nonzero = detail::query_nonzero_value(detail::node(base_expr), facts, Domain::Real, context);
        if (!nonzero) return nonzero;
        if (nonzero.value() == Tribool::True) return known_sign(1, target);
        if (nonzero.value() == Tribool::False) return known_sign(0, target);
        return nonnegative_sign(target);
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::positive_integer_power_sign(
    const InferenceEngine& engine,
    const SymbolicExpr& base_expr,
    const NumberNode& exponent,
    Sign target, ComputationContext& context) const {
    if (!is_positive_integer_number(exponent)) {
        return Tribool::Unknown;
    }
    auto base = engine.query_sign_of_checked(base_expr, Sign::NonNegative, context);
    if (!base) {
        return base;
    }
    if (base.value() == Tribool::True) {
        return nonnegative_sign(target);
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::nonzero_power_sign(
    const InferenceEngine& engine,
    const SymbolicExpr& base_expr,
    const NumberNode& exponent,
    Sign target, ComputationContext& context) const {
    if (!is_integer_number(exponent)) {
        return Tribool::Unknown;
    }
    auto base = engine.query_sign_of_checked(base_expr, Sign::NonZero, context);
    if (!base) {
        return base;
    }
    if (base.value() == Tribool::True) {
        return target == Sign::NonZero ? Tribool::True : target == Sign::Zero ? Tribool::False : Tribool::Unknown;
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::integer_power_sign(
    const InferenceEngine& engine,
    const SymbolicExpr& base_expr,
    const NumberNode& exponent,
    Sign target, ComputationContext& context) const {
    if (is_integer_number(exponent) && !is_even_integer_number(exponent)) {
        const detail::AssumptionFacts facts(engine);
        auto real = detail::query_real_value(detail::node(base_expr), facts, context);
        if (!real) return real;
        if (real.value() != Tribool::True)
            return nonzero_power_sign(engine, base_expr, exponent, target, context);
        return engine.query_sign_of_checked(base_expr, target, context);
    }
    auto even = even_power_sign(engine, base_expr, exponent, target, context);
    if (!even) {
        return even;
    }
    if (even.value() != Tribool::Unknown) {
        return even;
    }
    auto positive = positive_integer_power_sign(engine, base_expr, exponent, target, context);
    if (!positive) {
        return positive;
    }
    if (positive.value() != Tribool::Unknown) {
        return positive;
    }
    return nonzero_power_sign(engine, base_expr, exponent, target, context);
}

InferenceTriboolResult InferenceEngine::infer_power_property_checked(
    const PowerNode& node, Sign target, ComputationContext& context) const {
    try {
        auto defined = Impl::real_definedness(*this, node, context);
        if (!defined || defined.value() != Tribool::True) return defined;
        auto base_expr = LMCAS::detail::expression_from_node(node.base());
        auto exp_expr = LMCAS::detail::expression_from_node(node.exponent());
        auto exp_num = std::dynamic_pointer_cast<const NumberNode>(node.exponent());

        auto zero = impl_->zero_power_sign(*this, base_expr, exp_num.get(), target, context);
        if (!zero) {
            return zero;
        }
        if (zero.value() != Tribool::Unknown) {
            return zero;
        }

        auto positive = impl_->positive_power_sign(*this, base_expr, exp_expr, target, context);
        if (!positive) {
            return positive;
        }
        if (positive.value() != Tribool::Unknown) {
            return positive;
        }

        if (exp_num) {
            return impl_->integer_power_sign(*this, base_expr, *exp_num, target, context);
        }
        return InferenceTriboolResult::success(Tribool::Unknown);
    } catch (const CasError& error) {
        return InferenceTriboolResult::failure(error);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "power sign inference allocation failed",
            "infer_power_property_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_power_property_checked");
    }
}

}
