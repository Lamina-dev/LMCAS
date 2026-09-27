#include "assumption_context.hpp"
#include "internal/symbolic_ast.hpp"
#include "inference_engine.hpp"
#include "computation_context.hpp"

namespace LMCAS {


AssumptionTriboolResult AssumptionContext::is_positive(const SymbolicExpr& expr) const {
    return is_positive_checked(expr);
}

AssumptionTriboolResult AssumptionContext::is_positive_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return is_positive_checked(expr, context);
}

AssumptionTriboolResult AssumptionContext::is_positive_checked(
    const SymbolicExpr& expr, ComputationContext& context) const {
    InferenceEngine inference(*this);
    return inference.query_positive_checked(expr, context);
}

AssumptionTriboolResult AssumptionContext::is_negative(const SymbolicExpr& expr) const {
    return is_negative_checked(expr);
}

AssumptionTriboolResult AssumptionContext::is_negative_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return is_negative_checked(expr, context);
}

AssumptionTriboolResult AssumptionContext::is_negative_checked(
    const SymbolicExpr& expr, ComputationContext& context) const {
    InferenceEngine inference(*this);
    return inference.query_negative_checked(expr, context);
}

AssumptionTriboolResult AssumptionContext::is_nonnegative(const SymbolicExpr& expr) const {
    return is_nonnegative_checked(expr);
}

AssumptionTriboolResult AssumptionContext::is_nonnegative_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return is_nonnegative_checked(expr, context);
}

AssumptionTriboolResult AssumptionContext::is_nonnegative_checked(
    const SymbolicExpr& expr, ComputationContext& context) const {
    InferenceEngine inference(*this);
    return inference.query_nonnegative_checked(expr, context);
}

AssumptionTriboolResult AssumptionContext::is_real(const SymbolicExpr& expr) const {
    return is_real_checked(expr);
}

AssumptionTriboolResult AssumptionContext::is_real_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return is_real_checked(expr, context);
}

AssumptionTriboolResult AssumptionContext::is_real_checked(
    const SymbolicExpr& expr, ComputationContext& context) const {
    InferenceEngine inference(*this);
    return inference.query_real_checked(expr, context);
}

AssumptionTriboolResult AssumptionContext::is_integer(const SymbolicExpr& expr) const {
    return is_integer_checked(expr);
}

AssumptionTriboolResult AssumptionContext::is_integer_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return is_integer_checked(expr, context);
}

AssumptionTriboolResult AssumptionContext::is_integer_checked(
    const SymbolicExpr& expr, ComputationContext& context) const {
    InferenceEngine inference(*this);
    return inference.query_integer_checked(expr, context);
}

AssumptionTriboolResult AssumptionContext::is_nonzero(const SymbolicExpr& expr) const {
    return is_nonzero_checked(expr);
}

AssumptionTriboolResult AssumptionContext::is_nonzero_checked(const SymbolicExpr& expr) const {
    ComputationContext context;
    return is_nonzero_checked(expr, context);
}

AssumptionTriboolResult AssumptionContext::is_nonzero_checked(
    const SymbolicExpr& expr, ComputationContext& context) const {
    InferenceEngine inference(*this);
    return inference.query_nonzero_checked(expr, context);
}
AssumptionTriboolResult AssumptionContext::is_continuous(
    const std::string& symbol,
    const Interval& interval) const {
    return is_continuous_checked(symbol, interval);
}

AssumptionTriboolResult AssumptionContext::is_continuous_checked(
    const std::string& symbol,
    const Interval& interval) const {
    constexpr const char* operation = "is_continuous";
    if (symbol.empty()) {
        return AssumptionTriboolResult::failure(
            CasErrc::InvalidArgument, "symbol name must not be empty", operation);
    }
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        auto result = it->properties.is_continuous_checked(symbol, interval);
        if (!result) {
            return AssumptionTriboolResult::failure(result.error());
        }
        if (result.value()) {
            return AssumptionTriboolResult::success(Tribool::True);
        }
    }
    return AssumptionTriboolResult::success(Tribool::Unknown);
}

AssumptionTriboolResult AssumptionContext::is_differentiable(
    const std::string& symbol,
    const Interval& interval) const {
    return is_differentiable_checked(symbol, interval);
}

AssumptionTriboolResult AssumptionContext::is_differentiable_checked(
    const std::string& symbol,
    const Interval& interval) const {
    constexpr const char* operation = "is_differentiable";
    if (symbol.empty()) {
        return AssumptionTriboolResult::failure(
            CasErrc::InvalidArgument, "symbol name must not be empty", operation);
    }
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        auto result = it->properties.is_differentiable_checked(symbol, interval);
        if (!result) {
            return AssumptionTriboolResult::failure(result.error());
        }
        if (result.value()) {
            return AssumptionTriboolResult::success(Tribool::True);
        }
    }
    return AssumptionTriboolResult::success(Tribool::Unknown);
}

AssumptionTriboolResult AssumptionContext::is_positive_definite(
    const std::string& symbol) const {
    return is_positive_definite_checked(symbol);
}

AssumptionTriboolResult AssumptionContext::is_positive_definite_checked(
    const std::string& symbol) const {
    constexpr const char* operation = "is_positive_definite";
    if (symbol.empty()) {
        return AssumptionTriboolResult::failure(
            CasErrc::InvalidArgument, "symbol name must not be empty", operation);
    }
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        Definiteness d = it->properties.get_definiteness(symbol);
        if (d != Definiteness::Unknown) {
            return AssumptionTriboolResult::success(
                (d == Definiteness::PositiveDefinite) ? Tribool::True : Tribool::False);
        }
    }
    return AssumptionTriboolResult::success(Tribool::Unknown);
}

AssumptionTriboolResult AssumptionContext::is_positive_semidefinite(
    const std::string& symbol) const {
    return is_positive_semidefinite_checked(symbol);
}

AssumptionTriboolResult AssumptionContext::is_positive_semidefinite_checked(
    const std::string& symbol) const {
    constexpr const char* operation = "is_positive_semidefinite";
    if (symbol.empty()) {
        return AssumptionTriboolResult::failure(
            CasErrc::InvalidArgument, "symbol name must not be empty", operation);
    }
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        Definiteness d = it->properties.get_definiteness(symbol);
        if (d != Definiteness::Unknown) {
            return AssumptionTriboolResult::success(
                (d == Definiteness::PositiveDefinite ||
                 d == Definiteness::PositiveSemiDefinite) ? Tribool::True : Tribool::False);
        }
    }
    return AssumptionTriboolResult::success(Tribool::Unknown);
}

}
