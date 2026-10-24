#include "internal/inference_engine_impl.hpp"

namespace LMCAS {

Tribool InferenceEngine::Impl::number_positive(const NumberNode& num) {
    return num.is_positive() ? Tribool::True : Tribool::False;
}

Tribool InferenceEngine::Impl::number_negative(const NumberNode& num) {
    return !num.is_positive() && !num.is_zero() ? Tribool::True : Tribool::False;
}

Tribool InferenceEngine::Impl::number_nonnegative(const NumberNode& num) {
    return num.is_positive() || num.is_zero() ? Tribool::True : Tribool::False;
}

Tribool InferenceEngine::Impl::number_nonpositive(const NumberNode& num) {
    return num.is_positive() ? Tribool::False : Tribool::True;
}

Tribool InferenceEngine::Impl::number_integer(const NumberNode& num) {
    return is_integer_number(num) ? Tribool::True : Tribool::False;
}

Tribool InferenceEngine::Impl::number_nonzero(const NumberNode& num) {
    return num.is_zero() ? Tribool::False : Tribool::True;
}

InferenceTriboolResult InferenceEngine::Impl::query_rational(
    const InferenceEngine& engine,
    const SymbolicExpr& expr, ComputationContext& context) const {
    if (!LMCAS::detail::node(expr)) {
        return InferenceTriboolResult::failure(
            CasErrc::InvalidArgument,
            "inference expression must not be null",
            "query_domain_of_checked");
    }
    try {
        auto integer = engine.query_integer_checked(expr, context);
        if (!integer) {
            return InferenceTriboolResult::failure(integer.error());
        }
        if (integer.value() == Tribool::True) {
            return InferenceTriboolResult::success(Tribool::True);
        }
        if (auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
            if (std::holds_alternative<Rational>(num->value())) {
                return InferenceTriboolResult::success(Tribool::True);
            }
        }
        if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
            if (var->is_constant()) return InferenceTriboolResult::success(Tribool::False);
            const auto& props = ctx;
            if (props.has_domain(var->name(), Domain::Rational)) {
                return InferenceTriboolResult::success(Tribool::True);
            }
        }
        return InferenceTriboolResult::success(Tribool::Unknown);
    } catch (const CasError& error) {
        return InferenceTriboolResult::failure(error);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "inference query allocation failed",
            "query_domain_of_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "query_domain_of_checked");
    }
}

InferenceTriboolResult InferenceEngine::Impl::query_natural(const SymbolicExpr& expr, ComputationContext& context) const {
    return checked_inference_result<Tribool>(expr, "query_domain_of_checked", context,
        [&]() -> InferenceTriboolResult {
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
                return natural_number(*num);
            }
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
                if (var->is_constant()) return Tribool::False;
                const auto& props = ctx;
                if (props.has_domain(var->name(), Domain::Natural)) {
                    return Tribool::True;
                }
            }
            return Tribool::Unknown;
        });
}

Tribool InferenceEngine::Impl::natural_number(const NumberNode& num) {
    return is_integer_number(num) && (num.is_positive() || num.is_zero())
        ? Tribool::True : Tribool::False;
}


InferenceTriboolResult InferenceEngine::query_real_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return query_real_checked(expr, context);
}

InferenceTriboolResult InferenceEngine::query_real_checked(const SymbolicExpr& expr, ComputationContext& context) const {
    return checked_inference_result<Tribool>(expr, "query_real_checked", context,
        [&]() -> InferenceTriboolResult {
            DepthGuard guard(*this, *LMCAS::detail::node(expr));
            if (guard.should_abort()) {
                return Tribool::Unknown;
            }
            if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
                return Tribool::True;
            }

            if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
                if (var->is_constant()) return Tribool::True;
                const auto& props = impl_->ctx;
                if (props.has_domain(var->name(), Domain::Real)) {
                    return Tribool::True;
                }
                for (Sign sign : {Sign::Positive, Sign::Negative, Sign::Zero,
                                  Sign::NonNegative, Sign::NonPositive}) {
                    auto step = context.consume_steps(1, "query_real_checked");
                    if (!step) return InferenceTriboolResult::failure(step.error());
                    if (props.has_sign(var->name(), sign)) { return Tribool::True; }
                }
                return Tribool::Unknown;
            }

            return impl_->composite_real(*this, expr, context);
        });
}

InferenceTriboolResult InferenceEngine::query_integer_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return query_integer_checked(expr, context);
}

InferenceTriboolResult InferenceEngine::query_integer_checked(const SymbolicExpr& expr, ComputationContext& context) const {
    return checked_inference_result<Tribool>(expr, "query_integer_checked", context,
        [&]() -> InferenceTriboolResult {
            if (infinity_sign(expr, context) != 0) {
                return InferenceTriboolResult::success(Tribool::False);
            }
            DepthGuard guard(*this, *LMCAS::detail::node(expr));
            if (guard.should_abort()) {
                return Tribool::Unknown;
            }
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
                return Impl::number_integer(*num);
            }

            if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
                if (var->is_constant()) return Tribool::False;
                const auto& props = impl_->ctx;
                if (props.has_domain(var->name(), Domain::Integer)) {
                    return Tribool::True;
                }
                return Tribool::Unknown;
            }

            return impl_->composite_integer(*this, expr, context);
        });
}



InferenceTriboolResult InferenceEngine::query_algebraic_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return query_algebraic_checked(expr, context);
}

InferenceTriboolResult InferenceEngine::query_algebraic_checked(const SymbolicExpr& expr, ComputationContext& context) const {
    return checked_inference_result<Tribool>(expr, "query_algebraic_checked", context,
        [&]() -> InferenceTriboolResult {
            DepthGuard guard(*this, *LMCAS::detail::node(expr));
            if (guard.should_abort()) {
                return Tribool::Unknown;
            }
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
                if (std::holds_alternative<BigInt>(num->value())) {
                    return Tribool::True;
                }
                if (std::holds_alternative<Rational>(num->value())) {
                    return Tribool::True;
                }
                return Tribool::Unknown;
            }

            if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
                if (var->is_constant()) return var->name() == "phi" ? Tribool::True : Tribool::False;
                const auto& props = impl_->ctx;
                if (props.has_domain(var->name(), Domain::Algebraic)) {
                    return Tribool::True;
                }
                if (props.is_transcendental(var->name())) {
                    return Tribool::False;
                }
                return Tribool::Unknown;
            }
            return Tribool::Unknown;
        });
}


InferenceTriboolResult InferenceEngine::query_transcendental_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return query_transcendental_checked(expr, context);
}

InferenceTriboolResult InferenceEngine::query_transcendental_checked(const SymbolicExpr& expr, ComputationContext& context) const {
    return checked_inference_result<Tribool>(expr, "query_transcendental_checked", context,
        [&]() -> InferenceTriboolResult {
            DepthGuard guard(*this, *LMCAS::detail::node(expr));
            if (guard.should_abort()) {
                return Tribool::Unknown;
            }
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
                if (std::holds_alternative<BigInt>(num->value())) {
                    return Tribool::False;
                }
                if (std::holds_alternative<Rational>(num->value())) {
                    return Tribool::False;
                }
                return Tribool::Unknown;
            }

            if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
                if (var->is_constant()) return var->name() == "phi" ? Tribool::False : Tribool::True;
                const auto& props = impl_->ctx;
                if (props.is_transcendental(var->name())) {
                    return Tribool::True;
                }
                if (props.has_domain(var->name(), Domain::Algebraic)) {
                    return Tribool::False;
                }
                return Tribool::Unknown;
            }
            return Tribool::Unknown;
        });
}


InferenceTriboolResult InferenceEngine::query_finite_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return query_finite_checked(expr, context);
}

InferenceTriboolResult InferenceEngine::query_finite_checked(const SymbolicExpr& expr, ComputationContext& context) const {
    return checked_inference_result<Tribool>(expr, "query_finite_checked", context,
        [&]() -> InferenceTriboolResult {
            DepthGuard guard(*this, *LMCAS::detail::node(expr));
            if (guard.should_abort()) {
                return Tribool::Unknown;
            }
            if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
                return Tribool::True;
            }

            if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
                if (var->is_constant()) return Tribool::True;
                const auto& props = impl_->ctx;
                Finiteness f = props.get_finiteness(var->name());
                if (f == Finiteness::Finite) {
                    return Tribool::True;
                }
                if (f == Finiteness::Divergent) {
                    return Tribool::False;
                }
                return Tribool::Unknown;
            }
            return Tribool::Unknown;
        });
}


InferenceTriboolResult InferenceEngine::query_divergent_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return query_divergent_checked(expr, context);
}

InferenceTriboolResult InferenceEngine::query_divergent_checked(const SymbolicExpr& expr, ComputationContext& context) const {
    return checked_inference_result<Tribool>(expr, "query_divergent_checked", context,
        [&]() -> InferenceTriboolResult {
            DepthGuard guard(*this, *LMCAS::detail::node(expr));
            if (guard.should_abort()) {
                return Tribool::Unknown;
            }
            if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
                return Tribool::False;
            }

            if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
                if (var->is_constant()) return Tribool::False;
                const auto& props = impl_->ctx;
                Finiteness f = props.get_finiteness(var->name());
                if (f == Finiteness::Divergent) {
                    return Tribool::True;
                }
                if (f == Finiteness::Finite) {
                    return Tribool::False;
                }
                return Tribool::Unknown;
            }
            return Tribool::Unknown;
        });
}

bool is_integer_number(const NumberNode& num) {
    if (std::holds_alternative<BigInt>(num.value())) {
        return true;
    }
    return std::get<Rational>(num.value()).is_integer();
}
bool is_even_integer_number(const NumberNode& num) {
    if (!is_integer_number(num)) {
        return false;
    }

    if (std::holds_alternative<BigInt>(num.value())) {
        return std::get<BigInt>(num.value()).is_even();
    }
    return std::get<Rational>(num.value()).get_numerator().is_even();
}
bool is_positive_integer_number(const NumberNode& num) {
    return is_integer_number(num) && num.is_positive();
}
bool is_zero_number(const NumberNode& num) {
    return num.is_zero();
}

}
