#define _USE_MATH_DEFINES
#include "internal/inference_engine_impl.hpp"

namespace LMCAS {

InferenceTriboolResult InferenceEngine::query_positive_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return query_positive_checked(expr, context);
}

InferenceTriboolResult InferenceEngine::query_positive_checked(const SymbolicExpr& expr, ComputationContext& context) const {
    return checked_inference_result<Tribool>(expr, "query_positive_checked", context,
        [&]() -> InferenceTriboolResult {
            const int infinity = infinity_sign(expr, context);
            if (infinity != 0) {
                return InferenceTriboolResult::success(
                    infinity > 0 ? Tribool::True : Tribool::False);
            }
            DepthGuard guard(*this, *LMCAS::detail::node(expr));
            if (guard.should_abort()) {
                return Tribool::Unknown;
            }

            if (auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
                return Impl::number_positive(*num);
            }
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
                const auto& props = impl_->ctx;
                if (props.has_sign(var->name(), Sign::Positive)) {
                    return Tribool::True;
                }
                if (props.has_sign(var->name(), Sign::Negative) ||
                    props.has_sign(var->name(), Sign::Zero) ||
                    props.has_sign(var->name(), Sign::NonPositive)) {
                    return Tribool::False;
                }
                if (props.has_domain(var->name(), Domain::PositiveInt)) {
                    return Tribool::True;
                }
                return Tribool::Unknown;
            }

            return impl_->composite_positive(*this, expr, context);
        });
}

InferenceTriboolResult InferenceEngine::query_negative_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return query_negative_checked(expr, context);
}

InferenceTriboolResult InferenceEngine::query_negative_checked(const SymbolicExpr& expr, ComputationContext& context) const {
    return checked_inference_result<Tribool>(expr, "query_negative_checked", context,
        [&]() -> InferenceTriboolResult {
            const int infinity = infinity_sign(expr, context);
            if (infinity != 0) {
                return InferenceTriboolResult::success(
                    infinity < 0 ? Tribool::True : Tribool::False);
            }
            DepthGuard guard(*this, *LMCAS::detail::node(expr));
            if (guard.should_abort()) {
                return Tribool::Unknown;
            }

            if (auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
                return Impl::number_negative(*num);
            }

            if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
                const auto& props = impl_->ctx;
                if (props.has_sign(var->name(), Sign::Negative)) {
                    return Tribool::True;
                }
                if (props.has_sign(var->name(), Sign::Positive) ||
                    props.has_sign(var->name(), Sign::Zero) ||
                    props.has_sign(var->name(), Sign::NonNegative)) {
                    return Tribool::False;
                }
                return Tribool::Unknown;
            }

            return impl_->composite_negative(*this, expr, context);
        });
}

InferenceTriboolResult InferenceEngine::query_nonnegative_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return query_nonnegative_checked(expr, context);
}

InferenceTriboolResult InferenceEngine::query_nonnegative_checked(const SymbolicExpr& expr, ComputationContext& context) const {
    return checked_inference_result<Tribool>(expr, "query_nonnegative_checked", context,
        [&]() -> InferenceTriboolResult {
            const int infinity = infinity_sign(expr, context);
            if (infinity != 0) {
                return InferenceTriboolResult::success(
                    infinity > 0 ? Tribool::True : Tribool::False);
            }
            DepthGuard guard(*this, *LMCAS::detail::node(expr));
            if (guard.should_abort()) {
                return Tribool::Unknown;
            }

            if (auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
                return Impl::number_nonnegative(*num);
            }

            if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
                const auto& props = impl_->ctx;
                if (props.has_sign(var->name(), Sign::NonNegative)) {
                    return Tribool::True;
                }
                if (props.has_sign(var->name(), Sign::Negative)) {
                    return Tribool::False;
                }
                return Tribool::Unknown;
            }

            return impl_->composite_nonnegative(*this, expr, context);
        });
}

InferenceTriboolResult InferenceEngine::query_nonpositive_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return query_nonpositive_checked(expr, context);
}

InferenceTriboolResult InferenceEngine::query_nonpositive_checked(const SymbolicExpr& expr, ComputationContext& context) const {
    return checked_inference_result<Tribool>(expr, "query_nonpositive_checked", context,
        [&]() -> InferenceTriboolResult {
            DepthGuard guard(*this, *LMCAS::detail::node(expr));
            if (guard.should_abort()) {
                return Tribool::Unknown;
            }

            if (auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
                return Impl::number_nonpositive(*num);
            }

            if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
                const auto& props = impl_->ctx;
                if (props.has_sign(var->name(), Sign::NonPositive)) {
                    return Tribool::True;
                }
                if (props.has_sign(var->name(), Sign::Positive)) {
                    return Tribool::False;
                }
                return Tribool::Unknown;
            }

            return impl_->composite_nonpositive(*this, expr, context);
        });
}

InferenceTriboolResult InferenceEngine::query_nonzero_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return query_nonzero_checked(expr, context);
}

InferenceTriboolResult InferenceEngine::query_nonzero_checked(const SymbolicExpr& expr, ComputationContext& context) const {
    return checked_inference_result<Tribool>(expr, "query_nonzero_checked", context,
        [&]() -> InferenceTriboolResult {
            if (infinity_sign(expr, context) != 0) {
                return InferenceTriboolResult::success(Tribool::True);
            }
            DepthGuard guard(*this, *LMCAS::detail::node(expr));
            if (guard.should_abort()) {
                return Tribool::Unknown;
            }

            if (auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
                return Impl::number_nonzero(*num);
            }

            if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
                const auto& props = impl_->ctx;
                if (props.has_sign(var->name(), Sign::NonZero)) {
                    return Tribool::True;
                }
                if (props.has_sign(var->name(), Sign::Zero)) {
                    return Tribool::False;
                }
                return Tribool::Unknown;
            }

            return impl_->composite_nonzero(*this, expr, context);
        });
}

}
