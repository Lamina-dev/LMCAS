#define _USE_MATH_DEFINES
#include "internal/inference_engine_impl.hpp"
#include "internal/assumption_facts.hpp"

namespace LMCAS {



InferenceTriboolResult InferenceEngine::infer_add_domain_checked(
    const AddNode& node, Domain target, ComputationContext& context) const {
    try {
        if (node.operands().empty()) {
            return InferenceTriboolResult::success(Tribool::Unknown);
        }
        for (const auto& operand : node.operands()) {
            auto op_expr = LMCAS::detail::expression_from_node(operand);
            Tribool op_result = Tribool::Unknown;
            switch (target) {
                case Domain::Integer: {
                    auto queried = query_integer_checked(op_expr, context);
                    if (!queried) {
                        return queried;
                    }
                    op_result = queried.value();
                    break;
                }
                case Domain::Real: {
                    auto queried = query_real_checked(op_expr, context);
                    if (!queried) {
                        return queried;
                    }
                    op_result = queried.value();
                    break;
                }
                default:
                    return InferenceTriboolResult::success(Tribool::Unknown);
            }

            if (op_result != Tribool::True) {
                if (target == Domain::Real) {
                    auto integer = query_integer_checked(op_expr, context);
                    if (!integer) {
                        return integer;
                    }
                    if (integer.value() == Tribool::True) {
                        continue;
                    }
                }
                return InferenceTriboolResult::success(Tribool::Unknown);
            }
        }

        return InferenceTriboolResult::success(Tribool::True);
    } catch (const CasError& error) {
        return InferenceTriboolResult::failure(error);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "addition domain inference allocation failed",
            "infer_add_domain_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_add_domain_checked");
    }
}


InferenceTriboolResult InferenceEngine::infer_multiply_domain_checked(
    const MultiplyNode& node, Domain target, ComputationContext& context) const {
    try {
        auto defined = Impl::real_definedness(*this, node, context);
        if (!defined || defined.value() != Tribool::True) { return defined; }
        if (node.operands().empty()) {
            return InferenceTriboolResult::success(Tribool::Unknown);
        }

        return impl_->multiply_domain(*this, node, target, context);
    } catch (const CasError& error) {
        return InferenceTriboolResult::failure(error);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "multiplication domain inference allocation failed",
            "infer_multiply_domain_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_multiply_domain_checked");
    }
}

InferenceTriboolResult InferenceEngine::Impl::multiply_domain(
    const InferenceEngine& engine, const MultiplyNode& node, Domain target, ComputationContext& context) const {
    if (target != Domain::Integer && target != Domain::Real) { return Tribool::Unknown; }
    for (const auto& operand : node.operands()) {
        const auto expression = detail::expression_from_node(operand);
        auto result = target == Domain::Integer ? engine.query_integer_checked(expression, context)
                                                : real_or_integer(engine, expression, context);
        if (!result) { return result; }
        if (result.value() != Tribool::True) { return Tribool::Unknown; }
    }
    return Tribool::True;
}
InferenceTriboolResult InferenceEngine::infer_power_domain_checked(
    const PowerNode& node, Domain target, ComputationContext& context) const {
    try {
        auto defined = Impl::real_definedness(*this, node, context);
        if (!defined || defined.value() != Tribool::True) { return defined; }
        if (target == Domain::Real) {
            const detail::AssumptionFacts facts(*this);
            const std::shared_ptr<const SymbolicNode> borrowed(
                std::shared_ptr<const SymbolicNode>{}, &node);
            auto real = detail::query_real_value(borrowed, facts, context);
            if (!real) { return real; }
            return real.value() == Tribool::True ? Tribool::True : Tribool::Unknown;
        }
        auto base_expr = LMCAS::detail::expression_from_node(node.base());
        auto exp_num = std::dynamic_pointer_cast<const NumberNode>(node.exponent());

        if (exp_num && is_integer_number(*exp_num)) {
            return Impl::integral_power_domain(*this, base_expr, *exp_num, target, context);
        }

        return InferenceTriboolResult::success(Tribool::Unknown);
    } catch (const CasError& error) {
        return InferenceTriboolResult::failure(error);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "power domain inference allocation failed",
            "infer_power_domain_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_power_domain_checked");
    }
}

InferenceTriboolResult InferenceEngine::Impl::integral_power_domain(
    const InferenceEngine& engine, const SymbolicExpr& base,
    const NumberNode& exponent, Domain target, ComputationContext& context) {
    if (target == Domain::Integer) {
        auto base_int_result = engine.query_integer_checked(base, context);
        if (!base_int_result) {
            return base_int_result;
        }
        const Tribool base_int = base_int_result.value();
        if (base_int == Tribool::True && is_positive_integer_number(exponent)) {
            return InferenceTriboolResult::success(Tribool::True);
        }
        if (base_int == Tribool::True && is_zero_number(exponent)) {
            return InferenceTriboolResult::success(Tribool::True);
        }
    }
    return Tribool::Unknown;
}

}
