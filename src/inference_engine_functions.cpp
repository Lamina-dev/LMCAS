#define _USE_MATH_DEFINES
#include "internal/inference_engine_impl.hpp"
#include "internal/assumption_facts.hpp"
#include <type_traits>

namespace LMCAS {

InferenceTriboolResult InferenceEngine::Impl::real_definedness(
    const InferenceEngine& engine, const SymbolicNode& node, ComputationContext& context) {
    const detail::AssumptionFacts facts(engine);
    const std::shared_ptr<const SymbolicNode> borrowed(
        std::shared_ptr<const SymbolicNode>{}, &node);
    bool violated_real_restriction = false;
    auto defined = detail::query_definedness(
        borrowed, facts, Domain::Real, context, &violated_real_restriction);
    if (defined && defined.value() == Tribool::False) {
        if (!violated_real_restriction) {
            auto complex_defined =
                detail::query_definedness(borrowed, facts, Domain::Complex, context);
            if (!complex_defined) return InferenceTriboolResult::failure(complex_defined.error());
            if (complex_defined.value() != Tribool::False) {
                auto real_value = detail::query_real_value(borrowed, facts, context);
                if (!real_value) return InferenceTriboolResult::failure(real_value.error());
                if (real_value.value() == Tribool::Unknown) return Tribool::Unknown;
            }
        }
        return InferenceTriboolResult::failure(
            CasErrc::DomainError, "expression is undefined in the real domain",
            "inference.definedness");
    }
    return defined;
}

InferenceTriboolResult InferenceEngine::Impl::exponential_sign(
    const InferenceEngine& engine,
    const SymbolicExpr& arg_expr,
    Sign target, ComputationContext& context) const {
    const detail::AssumptionFacts facts(engine);
    auto real = detail::query_real_value(detail::node(arg_expr), facts, context);
    if (!real) {
        return real;
    }
    if (real.value() == Tribool::True) {
        return known_sign(1, target);
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::trigonometric_sign(
    const InferenceEngine& engine,
    const SymbolicExpr& arg_expr, ComputationContext& context) const {
    auto real = engine.query_domain_of_checked(arg_expr, Domain::Real, context);
    if (!real) {
        return real;
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::absolute_sign(
    const InferenceEngine& engine,
    const SymbolicExpr& arg_expr,
    Sign target, ComputationContext& context) const {
    auto arg_real_result = engine.query_domain_of_checked(arg_expr, Domain::Real, context);
    if (!arg_real_result) {
        return arg_real_result;
    }
    const Tribool arg_real = arg_real_result.value();
    if (arg_real != Tribool::True) {
        return InferenceTriboolResult::success(Tribool::Unknown);
    }

    auto arg_pos_result = engine.query_sign_of_checked(arg_expr, Sign::Positive, context);
    if (!arg_pos_result) {
        return arg_pos_result;
    }
    const Tribool arg_pos = arg_pos_result.value();
    auto arg_neg_result = engine.query_sign_of_checked(arg_expr, Sign::Negative, context);
    if (!arg_neg_result) {
        return arg_neg_result;
    }
    const Tribool arg_neg = arg_neg_result.value();
    auto arg_nz_result = engine.query_sign_of_checked(arg_expr, Sign::NonZero, context);
    if (!arg_nz_result) {
        return arg_nz_result;
    }
    const Tribool arg_nz = arg_nz_result.value();
    bool abs_is_positive = (arg_pos == Tribool::True) ||
                           (arg_neg == Tribool::True) ||
                           (arg_nz == Tribool::True);

    if (abs_is_positive) {
        return known_sign(1, target);
    }
    return nonnegative_sign(target);
}

InferenceTriboolResult InferenceEngine::Impl::logarithm_sign(
    const InferenceEngine&,
    const SymbolicExpr& arg_expr,
    Sign target, ComputationContext& context) const {
    if (const auto* number = dynamic_cast<const NumberNode*>(detail::node(arg_expr).get())) {
        const int sign = std::visit([](const auto& value) {
            using T = std::decay_t<decltype(value)>;
            return value > T(1) ? 1 : value < T(1) ? -1 : 0;
        }, number->value());
        return known_sign(sign, target);
    }
    auto step = context.consume_steps(1, "inference.logarithm");
    if (!step) {
        return InferenceTriboolResult::failure(step.error());
    }
    const auto one = detail::expression_from_node(detail::make_node<NumberNode>(BigInt(1)));
    auto greater = ctx.has_relation_checked(arg_expr, one, RelationOp::GT, context);
    if (!greater) {
        return InferenceTriboolResult::failure(greater.error());
    }
    if (greater.value()) {
        return known_sign(1, target);
    }
    auto equal = ctx.has_relation_checked(arg_expr, one, RelationOp::EQ, context);
    if (!equal) {
        return InferenceTriboolResult::failure(equal.error());
    }
    if (equal.value()) {
        return known_sign(0, target);
    }
    auto less = ctx.has_relation_checked(arg_expr, one, RelationOp::LT, context);
    if (!less) {
        return InferenceTriboolResult::failure(less.error());
    }
    if (less.value()) {
        return known_sign(-1, target);
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::square_root_sign(
    const InferenceEngine& engine,
    const SymbolicExpr& arg_expr,
    Sign target, ComputationContext& context) const {
    auto nonnegative = engine.query_sign_of_checked(arg_expr, Sign::NonNegative, context);
    if (!nonnegative) {
        return nonnegative;
    }
    if (nonnegative.value() == Tribool::True) {
        return nonnegative_sign(target);
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::real_or_integer(
    const InferenceEngine& engine,
    const SymbolicExpr& arg_expr, ComputationContext& context) const {
    auto arg_real_result = engine.query_real_checked(arg_expr, context);
    if (!arg_real_result) {
        return arg_real_result;
    }
    const Tribool arg_real = arg_real_result.value();
    if (arg_real == Tribool::True) {
        return InferenceTriboolResult::success(Tribool::True);
    }

    auto arg_int_result = engine.query_integer_checked(arg_expr, context);
    if (!arg_int_result) {
        return arg_int_result;
    }
    const Tribool arg_int = arg_int_result.value();
    if (arg_int == Tribool::True) {
        return InferenceTriboolResult::success(Tribool::True);
    }
    return InferenceTriboolResult::success(Tribool::Unknown);
}

InferenceTriboolResult InferenceEngine::Impl::exponential_domain(
    const InferenceEngine& engine,
    const SymbolicExpr& arg_expr,
    Domain target, ComputationContext& context) const {
    if (target != Domain::Real) {
        return Tribool::Unknown;
    }
    const detail::AssumptionFacts facts(engine);
    auto real = detail::query_real_value(detail::node(arg_expr), facts, context);
    if (!real) return real;
    return real.value() == Tribool::True ? Tribool::True : Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::absolute_domain(
    const InferenceEngine& engine,
    const SymbolicExpr& arg_expr,
    Domain target, ComputationContext& context) const {
    if (target == Domain::Real) {
        return real_or_integer(engine, arg_expr, context);
    }
    if (target == Domain::Integer) {
        auto integer = engine.query_integer_checked(arg_expr, context);
        if (!integer) {
            return integer;
        }
        if (integer.value() == Tribool::True) {
            return Tribool::True;
        }
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::logarithm_domain(
    const InferenceEngine& engine,
    const SymbolicExpr& arg_expr,
    Domain target, ComputationContext& context) const {
    if (target != Domain::Real) return Tribool::Unknown;
    const detail::AssumptionFacts facts(engine);
    auto positive = detail::query_positive_value(detail::node(arg_expr), facts, context);
    if (!positive) return positive;
    return positive.value() == Tribool::True ? Tribool::True : Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::square_root_domain(
    const InferenceEngine& engine,
    const SymbolicExpr& arg_expr,
    Domain target, ComputationContext& context) const {
    if (target == Domain::Real) {
        auto arg_nonneg_result = engine.query_sign_of_checked(arg_expr, Sign::NonNegative, context);
        if (!arg_nonneg_result) {
            return arg_nonneg_result;
        }
        const Tribool arg_nonneg = arg_nonneg_result.value();
        if (arg_nonneg == Tribool::True) {
            return InferenceTriboolResult::success(Tribool::True);
        }
    }
    return InferenceTriboolResult::success(Tribool::Unknown);
}

InferenceTriboolResult InferenceEngine::Impl::arctangent_domain(
    const InferenceEngine& engine,
    const SymbolicExpr& arg_expr,
    Domain target, ComputationContext& context) const {
    if (target == Domain::Real) {
        auto arg_real_result = engine.query_real_checked(arg_expr, context);
        if (!arg_real_result) {
            return arg_real_result;
        }
        const Tribool arg_real = arg_real_result.value();
        if (arg_real == Tribool::True) {
            return InferenceTriboolResult::success(Tribool::True);
        }
    }
    return InferenceTriboolResult::success(Tribool::Unknown);
}

InferenceTriboolResult InferenceEngine::infer_function_property_checked(
    const FunctionNode& node, Sign target, ComputationContext& context) const {
    try {
        if (node.arguments().empty()) {
            return InferenceTriboolResult::success(Tribool::Unknown);
        }
        if (node.type() == FunctionNode::FuncType::Ln ||
            node.type() == FunctionNode::FuncType::Exp) {
            auto defined = Impl::real_definedness(*this, node, context);
            if (!defined || defined.value() != Tribool::True) { return defined; }
        }

        auto arg_expr = LMCAS::detail::expression_from_node(node.arguments()[0]);
        return impl_->function_sign(*this, node.type(), arg_expr, target, context);
    } catch (const CasError& error) {
        return InferenceTriboolResult::failure(error);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "function sign inference allocation failed",
            "infer_function_property_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_function_property_checked");
    }
}

InferenceTriboolResult InferenceEngine::Impl::function_sign(
    const InferenceEngine& engine, FunctionNode::FuncType type,
    const SymbolicExpr& arg_expr, Sign target, ComputationContext& context) const {
    switch (type) {
        case FunctionNode::FuncType::Exp: {
            return exponential_sign(engine, arg_expr, target, context);
        }
        case FunctionNode::FuncType::Sin:
        case FunctionNode::FuncType::Cos:
        case FunctionNode::FuncType::Tan: {
            return trigonometric_sign(engine, arg_expr, context);
        }
        case FunctionNode::FuncType::Abs: {
            return absolute_sign(engine, arg_expr, target, context);
        }
        case FunctionNode::FuncType::Ln: {
            return logarithm_sign(engine, arg_expr, target, context);
        }
        case FunctionNode::FuncType::Sqrt: {
            return square_root_sign(engine, arg_expr, target, context);
        }
        case FunctionNode::FuncType::ArcTan: {
            return engine.query_sign_of_checked(arg_expr, target, context);
        }
        default:
            return InferenceTriboolResult::success(Tribool::Unknown);
    }
}

InferenceTriboolResult InferenceEngine::infer_function_domain_checked(
    const FunctionNode& node, Domain target, ComputationContext& context) const {
    try {
        if (node.arguments().empty()) {
            return InferenceTriboolResult::success(Tribool::Unknown);
        }
        if (node.type() == FunctionNode::FuncType::Ln ||
            node.type() == FunctionNode::FuncType::Exp) {
            auto defined = Impl::real_definedness(*this, node, context);
            if (!defined || defined.value() != Tribool::True) { return defined; }
        }

        auto arg_expr = LMCAS::detail::expression_from_node(node.arguments()[0]);
        return impl_->function_domain(*this, node.type(), arg_expr, target, context);
    } catch (const CasError& error) {
        return InferenceTriboolResult::failure(error);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "function domain inference allocation failed",
            "infer_function_domain_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_function_domain_checked");
    }
}

InferenceTriboolResult InferenceEngine::Impl::function_domain(
    const InferenceEngine& engine, FunctionNode::FuncType type,
    const SymbolicExpr& arg_expr, Domain target, ComputationContext& context) const {
    switch (type) {
        case FunctionNode::FuncType::Exp: {
            return exponential_domain(engine, arg_expr, target, context);
        }
        case FunctionNode::FuncType::Sin:
        case FunctionNode::FuncType::Cos:
        case FunctionNode::FuncType::Tan: {
            if (target == Domain::Real) {
                return real_or_integer(engine, arg_expr, context);
            }
            return InferenceTriboolResult::success(Tribool::Unknown);
        }
        case FunctionNode::FuncType::Abs: {
            return absolute_domain(engine, arg_expr, target, context);
        }
        case FunctionNode::FuncType::Ln: {
            return logarithm_domain(engine, arg_expr, target, context);
        }
        case FunctionNode::FuncType::Sqrt: {
            return square_root_domain(engine, arg_expr, target, context);
        }
        case FunctionNode::FuncType::ArcTan: {
            return arctangent_domain(engine, arg_expr, target, context);
        }
        default:
            return InferenceTriboolResult::success(Tribool::Unknown);
    }
}

}
