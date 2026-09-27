#define _USE_MATH_DEFINES
#include "internal/inference_engine_impl.hpp"

namespace LMCAS {

InferenceEngine::InferenceEngine(const AssumptionContext& ctx)
    : impl_(std::make_unique<Impl>(ctx)) {}

InferenceEngine::~InferenceEngine() = default;
InferenceEngine::InferenceEngine(InferenceEngine&&) noexcept = default;
InferenceEngine& InferenceEngine::operator=(InferenceEngine&&) noexcept = default;

void InferenceEngine::set_max_depth(int depth) {
    if (depth > 0) {
        impl_->max_depth = depth;
    }
}

int InferenceEngine::get_max_depth() const {
    return impl_->max_depth;
}

InferenceEngine::DepthGuard::DepthGuard(
    const InferenceEngine& engine, const SymbolicNode& node)
    : engine_(engine), node_(&node) {
    engine_.impl_->current_depth++;

    if (engine_.impl_->current_depth > engine_.impl_->max_depth) {
        abort_ = true;
        return;
    }

    // Insert and detect cycles atomically. If allocation throws, restore the
    // depth counter because a throwing constructor has no matching destructor.
    try {
        inserted_ = engine_.impl_->visited.insert(node_).second;
    } catch (...) {
        engine_.impl_->current_depth--;
        if (engine_.impl_->current_depth == 0) {
            engine_.impl_->visited.clear();
        }
        throw;
    }
    if (!inserted_) {
        abort_ = true;
    }
}

InferenceEngine::DepthGuard::~DepthGuard() {
    if (inserted_) {
        engine_.impl_->visited.erase(node_);
    }

    engine_.impl_->current_depth--;

    if (engine_.impl_->current_depth == 0) {
        engine_.impl_->visited.clear();
    }
}


InferenceTriboolResult InferenceEngine::query_sign_of_checked(const SymbolicExpr& expr, Sign sign, ComputationContext& context) const {
    switch (sign) {
        case Sign::Positive:    return query_positive_checked(expr, context);
        case Sign::Negative:    return query_negative_checked(expr, context);
        case Sign::NonNegative: return query_nonnegative_checked(expr, context);
        case Sign::NonPositive: return query_nonpositive_checked(expr, context);
        case Sign::Zero:
            // Zero means both NonNegative and NonPositive
            return impl_->query_zero(*this, expr, context);
        case Sign::NonZero:     return query_nonzero_checked(expr, context);
    }
    return InferenceTriboolResult::success(Tribool::Unknown);
}
InferenceTriboolResult InferenceEngine::Impl::composite_positive(
    const InferenceEngine& engine,
    const SymbolicExpr& expr, ComputationContext& context) const {
    if (auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr))) {
        auto result =
            engine.infer_add_sign_checked(*add, Sign::Positive, context);
        if (!result) {
            return result;
        }
        if (result.value() != Tribool::Unknown) {
            return result;
        }
        return engine.infer_sign_from_relations_checked(expr, Sign::Positive, context);
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        auto result =
            engine.infer_multiply_sign_checked(*mul, Sign::Positive, context);
        if (!result) {
            return result;
        }
        if (result.value() != Tribool::Unknown) {
            return result;
        }
        return engine.infer_sign_from_relations_checked(expr, Sign::Positive, context);
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(expr))) {
        return engine.infer_power_property_checked(*pow, Sign::Positive, context);
    }
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(expr))) {
        return engine.infer_function_property_checked(*func, Sign::Positive, context);
    }
    return engine.infer_sign_from_relations_checked(expr, Sign::Positive, context);
}

InferenceTriboolResult InferenceEngine::Impl::composite_negative(
    const InferenceEngine& engine,
    const SymbolicExpr& expr, ComputationContext& context) const {
    if (auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr))) {
        return engine.infer_add_sign_checked(*add, Sign::Negative, context);
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        return engine.infer_multiply_sign_checked(*mul, Sign::Negative, context);
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(expr))) {
        return engine.infer_power_property_checked(*pow, Sign::Negative, context);
    }
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(expr))) {
        return engine.infer_function_property_checked(*func, Sign::Negative, context);
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::composite_nonnegative(
    const InferenceEngine& engine,
    const SymbolicExpr& expr, ComputationContext& context) const {
    if (auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr))) {
        auto result =
            engine.infer_add_sign_checked(*add, Sign::NonNegative, context);
        if (!result) {
            return result;
        }
        if (result.value() != Tribool::Unknown) {
            return result;
        }
        return engine.infer_sign_from_relations_checked(expr, Sign::NonNegative, context);
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        auto result =
            engine.infer_multiply_sign_checked(*mul, Sign::NonNegative, context);
        if (!result) {
            return result;
        }
        if (result.value() != Tribool::Unknown) {
            return result;
        }
        return engine.infer_sign_from_relations_checked(expr, Sign::NonNegative, context);
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(expr))) {
        return engine.infer_power_property_checked(*pow, Sign::NonNegative, context);
    }
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(expr))) {
        return engine.infer_function_property_checked(*func, Sign::NonNegative, context);
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::composite_nonpositive(
    const InferenceEngine& engine,
    const SymbolicExpr& expr, ComputationContext& context) const {
    if (auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr))) {
        return engine.infer_add_sign_checked(*add, Sign::NonPositive, context);
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        return engine.infer_multiply_sign_checked(*mul, Sign::NonPositive, context);
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(expr))) {
        return engine.infer_power_property_checked(*pow, Sign::NonPositive, context);
    }
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(expr))) {
        return engine.infer_function_property_checked(*func, Sign::NonPositive, context);
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::composite_real(
    const InferenceEngine& engine,
    const SymbolicExpr& expr, ComputationContext& context) const {
    if (auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr))) {
        return engine.infer_add_domain_checked(*add, Domain::Real, context);
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        return engine.infer_multiply_domain_checked(*mul, Domain::Real, context);
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(expr))) {
        return engine.infer_power_domain_checked(*pow, Domain::Real, context);
    }
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(expr))) {
        return engine.infer_function_domain_checked(*func, Domain::Real, context);
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::composite_integer(
    const InferenceEngine& engine,
    const SymbolicExpr& expr, ComputationContext& context) const {
    if (auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr))) {
        return engine.infer_add_domain_checked(*add, Domain::Integer, context);
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        return engine.infer_multiply_domain_checked(*mul, Domain::Integer, context);
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(expr))) {
        return engine.infer_power_domain_checked(*pow, Domain::Integer, context);
    }
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(expr))) {
        return engine.infer_function_domain_checked(*func, Domain::Integer, context);
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::sum_nonzero(
    const InferenceEngine& engine,
    const AddNode& add, ComputationContext& context) const {
    auto positive =
        engine.infer_add_sign_checked(add, Sign::Positive, context);
    if (!positive) {
        return positive;
    }
    if (positive.value() == Tribool::True) {
        return Tribool::True;
    }
    auto negative =
        engine.infer_add_sign_checked(add, Sign::Negative, context);
    if (!negative) {
        return negative;
    }
    if (negative.value() == Tribool::True) {
        return Tribool::True;
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::composite_nonzero(
    const InferenceEngine& engine,
    const SymbolicExpr& expr, ComputationContext& context) const {
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        return engine.infer_multiply_sign_checked(*mul, Sign::NonZero, context);
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(expr))) {
        return engine.infer_power_property_checked(*pow, Sign::NonZero, context);
    }
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(expr))) {
        return engine.infer_function_property_checked(*func, Sign::NonZero, context);
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr))) {
        return sum_nonzero(engine, *add, context);
    }
    return Tribool::Unknown;
}

InferenceTriboolResult InferenceEngine::Impl::query_zero(
    const InferenceEngine& engine,
    const SymbolicExpr& expr, ComputationContext& context) const {
    auto nn = engine.query_nonnegative_checked(expr, context);
    if (!nn) {
        return nn;
    }
    auto np = engine.query_nonpositive_checked(expr, context);
    if (!np) {
        return np;
    }
    if (nn.value() == Tribool::True && np.value() == Tribool::True) {
        return InferenceTriboolResult::success(Tribool::True);
    }
    if (nn.value() == Tribool::False || np.value() == Tribool::False) {
        return InferenceTriboolResult::success(Tribool::False);
    }
    auto nonzero = engine.query_nonzero_checked(expr, context);
    if (!nonzero) return nonzero;
    if (nonzero.value() == Tribool::True) return Tribool::False;
    return InferenceTriboolResult::success(Tribool::Unknown);
}

InferenceTriboolResult InferenceEngine::query_domain_of_checked(const SymbolicExpr& expr, Domain domain, ComputationContext& context) const {
    switch (domain) {
        case Domain::Integer:  return query_integer_checked(expr, context);
        case Domain::Real:     return query_real_checked(expr, context);
        case Domain::Rational: {
            return impl_->query_rational(*this, expr, context);
        }
        case Domain::Natural: {
            return impl_->query_natural(expr, context);
        }
        default:
            return InferenceTriboolResult::success(Tribool::Unknown);
    }
}

} // namespace LMCAS
