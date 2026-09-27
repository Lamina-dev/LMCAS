#include "query_interface.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/query_support.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"

namespace LMCAS {

// Construction

QueryInterface::QueryInterface(const AssumptionContext& ctx)
    : ctx_(ctx), observed_generation_(ctx.cache_generation()) {}

// Cache management

void QueryInterface::invalidate_cache() const {
    cache_.clear();
    observed_generation_ = ctx_.cache_generation();
}


QueryTriboolResult QueryInterface::cached_query_checked(
    const SymbolicExpr& expr,
    PropType prop,
    const std::string& operation,
    const std::function<QueryTriboolResult()>& compute,
    const std::string& variable) const {
    if (!LMCAS::detail::node(expr)) {
        return QueryTriboolResult::failure(
            CasErrc::InvalidArgument, "query expression must not be null", operation);
    }

    // Check if the context has been mutated since we last validated the cache.
    // If so, invalidate and update our observed generation.
    uint64_t current_gen = ctx_.cache_generation();
    if (current_gen != observed_generation_) {
        cache_.clear();
        observed_generation_ = current_gen;
    }

    CacheKey key{LMCAS::detail::node(expr)->hash(), prop, variable};
    auto it = cache_.find(key);
    if (it != cache_.end()) {
        for (const auto& entry : it->second) {
            if (LMCAS::detail::node(entry.expression)->equals(*LMCAS::detail::node(expr))) {
                return QueryTriboolResult::success(entry.result);
            }
        }
    }

    auto result = compute();
    if (!result) return result;
    cache_[key].push_back(CacheEntry{expr, result.value()});
    return result;
}
// Public query methods

QueryTriboolResult QueryInterface::query_positive(const SymbolicExpr& expr) const {
    return query_positive_checked(expr);
}

QueryTriboolResult QueryInterface::query_positive_checked(const SymbolicExpr& expr) const {
    return query_detail::checked_expression_result(expr, "query_positive_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Positive, "query_positive_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        int sign = get_infinity_sign(expr);
                        if (sign > 0) return Tribool::True;
                        if (sign < 0) return Tribool::False;
                        return Tribool::Unknown;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_positive_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_negative(const SymbolicExpr& expr) const {
    return query_negative_checked(expr);
}

QueryTriboolResult QueryInterface::query_negative_checked(const SymbolicExpr& expr) const {
    return query_detail::checked_expression_result(expr, "query_negative_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Negative, "query_negative_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        int sign = get_infinity_sign(expr);
                        if (sign < 0) return Tribool::True;
                        if (sign > 0) return Tribool::False;
                        return Tribool::Unknown;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_negative_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_nonnegative(const SymbolicExpr& expr) const {
    return query_nonnegative_checked(expr);
}

QueryTriboolResult QueryInterface::query_nonnegative_checked(const SymbolicExpr& expr) const {
    return query_detail::checked_expression_result(expr, "query_nonnegative_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::NonNegative, "query_nonnegative_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        int sign = get_infinity_sign(expr);
                        if (sign > 0) return Tribool::True;
                        if (sign < 0) return Tribool::False;
                        return Tribool::Unknown;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_nonnegative_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_real(const SymbolicExpr& expr) const {
    return query_real_checked(expr);
}

QueryTriboolResult QueryInterface::query_real_checked(const SymbolicExpr& expr) const {
    return query_detail::checked_expression_result(expr, "query_real_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Real, "query_real_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::Unknown;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_real_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_integer(const SymbolicExpr& expr) const {
    return query_integer_checked(expr);
}

QueryTriboolResult QueryInterface::query_integer_checked(const SymbolicExpr& expr) const {
    return query_detail::checked_expression_result(expr, "query_integer_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Integer, "query_integer_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::False;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::False;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_integer_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_nonzero(const SymbolicExpr& expr) const {
    return query_nonzero_checked(expr);
}

QueryTriboolResult QueryInterface::query_nonzero_checked(const SymbolicExpr& expr) const {
    return query_detail::checked_expression_result(expr, "query_nonzero_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::NonZero, "query_nonzero_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::True;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_nonzero_checked(expr);
                });
            return result;
        });
}

// Extended property queries

QueryTriboolResult QueryInterface::query_algebraic(const SymbolicExpr& expr) const {
    return query_algebraic_checked(expr);
}

QueryTriboolResult QueryInterface::query_algebraic_checked(const SymbolicExpr& expr) const {
    return query_detail::checked_expression_result(expr, "query_algebraic_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Algebraic, "query_algebraic_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::False;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_algebraic_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_transcendental(const SymbolicExpr& expr) const {
    return query_transcendental_checked(expr);
}

QueryTriboolResult QueryInterface::query_transcendental_checked(const SymbolicExpr& expr) const {
    return query_detail::checked_expression_result(expr, "query_transcendental_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Transcendental, "query_transcendental_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::False;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_transcendental_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_finite(const SymbolicExpr& expr) const {
    return query_finite_checked(expr);
}

QueryTriboolResult QueryInterface::query_finite_checked(const SymbolicExpr& expr) const {
    return query_detail::checked_expression_result(expr, "query_finite_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Finite, "query_finite_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::False;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_finite_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_divergent(const SymbolicExpr& expr) const {
    return query_divergent_checked(expr);
}

QueryTriboolResult QueryInterface::query_divergent_checked(const SymbolicExpr& expr) const {
    return query_detail::checked_expression_result(expr, "query_divergent_checked",
        [&]() {
            auto result = cached_query_checked(
                expr, PropType::Divergent, "query_divergent_checked", [&]() -> QueryTriboolResult {
                    if (is_unhandled_type(expr)) {
                        return Tribool::Unknown;
                    }
                    if (std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr)) &&
                        is_nan_number(expr)) {
                        return Tribool::Unknown;
                    }
                    if (is_infinity_node(expr)) {
                        return Tribool::True;
                    }
                    InferenceEngine engine(ctx_);
                    return engine.query_divergent_checked(expr);
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_periodic(
    const SymbolicExpr& expr, const std::string& variable) const {
    return query_periodic_checked(expr, variable);
}

QueryTriboolResult QueryInterface::query_periodic_checked(
    const SymbolicExpr& expr, const std::string& variable) const {
    if (variable.empty()) {
        return QueryTriboolResult::failure(CasErrc::InvalidArgument,
            "independent variable must not be empty", "query_periodic_checked");
    }
    return query_detail::checked_expression_result(expr, "query_periodic_checked", [&]() {
        return cached_query_checked(expr, PropType::Periodic, "query_periodic_checked",
            [&]() -> QueryTriboolResult {
                InferenceEngine engine(ctx_);
                return engine.query_periodic_checked(expr, variable);
            }, variable);
    });
}

QueryPeriodResult QueryInterface::get_period(
    const SymbolicExpr& expr, const std::string& variable) const {
    return get_period_checked(expr, variable);
}

QueryPeriodResult QueryInterface::get_period_checked(
    const SymbolicExpr& expr, const std::string& variable) const {
    InferenceEngine engine(ctx_);
    return engine.infer_period_checked(expr, variable);
}

QueryTriboolResult QueryInterface::query_positive_definite(
    const SymbolicExpr& expr) const {
    return query_positive_definite_checked(expr);
}

QueryTriboolResult QueryInterface::query_positive_definite_checked(const SymbolicExpr& expr) const {
    return query_detail::checked_expression_result(expr, "query_positive_definite_checked",
        [&]() {
            auto result = cached_query_checked(
                expr,
                PropType::PositiveDefinite,
                "query_positive_definite_checked",
                [&]() -> QueryTriboolResult {
                    if (auto var = std::dynamic_pointer_cast<const VariableNode>(
                            LMCAS::detail::node(expr))) {
                        const auto& props = ctx_;
                        Definiteness d = props.get_definiteness(var->name());
                        if (d == Definiteness::PositiveDefinite) return Tribool::True;
                        if (d == Definiteness::NegativeDefinite ||
                            d == Definiteness::NegativeSemiDefinite ||
                            d == Definiteness::Indefinite) return Tribool::False;
                        return Tribool::Unknown;
                    }
                    return Tribool::Unknown;
                });
            return result;
        });
}

QueryTriboolResult QueryInterface::query_positive_semidefinite(
    const SymbolicExpr& expr) const {
    return query_positive_semidefinite_checked(expr);
}

QueryTriboolResult QueryInterface::query_positive_semidefinite_checked(const SymbolicExpr& expr) const {
    return query_detail::checked_expression_result(expr, "query_positive_semidefinite_checked",
        [&]() {
            auto result = cached_query_checked(
                expr,
                PropType::PositiveSemiDefinite,
                "query_positive_semidefinite_checked",
                [&]() -> QueryTriboolResult {
                    if (auto var = std::dynamic_pointer_cast<const VariableNode>(
                            LMCAS::detail::node(expr))) {
                        const auto& props = ctx_;
                        Definiteness d = props.get_definiteness(var->name());
                        if (d == Definiteness::PositiveDefinite ||
                            d == Definiteness::PositiveSemiDefinite) return Tribool::True;
                        if (d == Definiteness::NegativeDefinite ||
                            d == Definiteness::NegativeSemiDefinite ||
                            d == Definiteness::Indefinite) return Tribool::False;
                        return Tribool::Unknown;
                    }
                    return Tribool::Unknown;
                });
            return result;
        });
}

}
