/**
 * @file assumption_context.hpp
 * @brief AssumptionContext class with scoped push/pop management.
 *
 * Provides a scope stack where each scope has its own PropertyStore and
 * RelationStore. Queries read through from the top scope down to the root,
 * with child scope declarations shadowing parent declarations.
 */
#pragma once

#include "assumption.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "result.hpp"
#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_set>
#include <optional>
#include <sstream>
#include <type_traits>
#include <initializer_list>
#include <variant>

namespace LMCAS {

// Forward declarations
class InferenceEngine;
class ComputationContext;
struct Interval;

using AssumptionVoidResult = Result<void>;
using AssumptionTriboolResult = Result<Tribool>;

/**
 * @brief Scoped assumption context with push/pop management.
 *
 * Maintains a stack of scopes, each containing a PropertyStore and RelationStore.
 * Queries use read-through semantics: searching from the top scope down to root,
 * 查询按子作用域到父作用域的顺序定位首个符号声明.
 * 子作用域声明覆盖查询结果,同时保留父作用域状态;弹出作用域后恢复先前结果.
 */
class LMCAS_API AssumptionContext {
public:
    AssumptionContext();

    // --- Scope management ---

    /// Push a new scope onto the stack.
    void push();

    /// Pop the current scope.
    AssumptionVoidResult pop();

    /// Return the current nesting depth (number of scopes on the stack).
    int depth() const;

    // --- Direct access to current (top) scope stores ---

    /// Access the current scope's PropertyStore (mutable, for declarations).
    PropertyStore& current_properties();

    /// Access the current scope's PropertyStore (const).
    const PropertyStore& current_properties() const;

    /// Access the current scope's RelationStore (mutable).
    RelationStore& current_relations();

    /// Access the current scope's RelationStore (const).
    const RelationStore& current_relations() const;

    // --- Read-through query methods ---
    // These search from the top scope down to root, returning the first
    // scope's value that has a declaration for the symbol. This provides
    // shadowing semantics: child declarations override parent declarations.

    /// Check if a symbol has a specific sign (read-through all scopes).
    bool has_sign(const std::string& symbol, Sign sign) const;

    /// Check if a symbol has at least the given domain specificity (read-through).
    bool has_domain(const std::string& symbol, Domain domain) const;

    /// Get the most specific domain for a symbol (read-through).
    Domain get_domain(const std::string& symbol) const;

    /// Get all signs for a symbol (read-through).
    std::unordered_set<Sign, SignHash> get_signs(const std::string& symbol) const;

    /// Get parity for a symbol (read-through).
    Parity get_parity(const std::string& symbol) const;

    /// Get boundedness for a symbol (read-through).
    Boundedness get_boundedness(const std::string& symbol) const;

    /// Get interval bounds for a symbol (read-through).
    std::optional<Interval> get_bounds(const std::string& symbol) const;

    /// Declare domain for a variable in the current scope.
    /// Returns an explicit error if the variable name is empty or contradictory.
    AssumptionVoidResult assume_domain(const std::string& variable, Domain domain);

    /** @brief Declares a domain and reports invalid or contradictory input. */
    AssumptionVoidResult assume_domain_checked(const std::string& variable, Domain domain);

    /// Declare sign for a variable in the current scope.
    /// Returns an explicit error if the variable name is empty or contradictory.
    AssumptionVoidResult assume_sign(const std::string& variable, Sign sign);

    /** @brief Declares a sign and reports invalid or contradictory input. */
    AssumptionVoidResult assume_sign_checked(const std::string& variable, Sign sign);

    /// Store a relational constraint. The expression's root must be a RelationalNode.
    /// Invalid or contradictory relations are returned as an explicit error.
    AssumptionVoidResult assume(const SymbolicExpr& relation);

    /** @brief Stores a relation and reports invalid or contradictory input. */
    AssumptionVoidResult assume_checked(const SymbolicExpr& relation);

    /**
     * @brief Store a conditional assumption: if condition holds, conclude conclusion.
     *
     * The conditional is stored in the current scope and discarded on scope pop.
     * When a query depends on the conclusion, the system checks if the condition
     * is satisfied by the current assumption state. If satisfied, the conclusion
     * is used; if not verifiable, the conditional is inactive.
     *
     * @param condition  A relational expression (e.g., x > 1) serving as the guard
     * @param conclusion A relational expression (e.g., ln(x) > 0) to use when condition holds
     * Invalid or contradictory relations are returned as an explicit error.
     */
    AssumptionVoidResult assume_conditional(const SymbolicExpr& condition, const SymbolicExpr& conclusion);

    /** @brief Stores a guarded conclusion and reports invalid input. */
    AssumptionVoidResult assume_conditional_checked(
        const SymbolicExpr& condition,
        const SymbolicExpr& conclusion);

    /**
     * @brief A conditional assumption: "if condition then conclusion".
     */
    struct ConditionalAssumption {
        SymbolicExpr condition;
        SymbolicExpr conclusion;
    };

    /**
     * @brief Retrieve all active conditional assumptions across all scopes (top to bottom).
     *
     * Used by the QueryInterface/InferenceEngine to check if any conditional
     * conclusions are applicable given the current assumption state.
     *
     * @return Vector of all conditional assumptions from all scopes (most recent first)
     */
    std::vector<ConditionalAssumption> get_active_conditionals() const;

    /**
     * @brief Evaluate whether a condition expression is satisfied by the current state.
     *
     * Checks if the condition (a RelationalNode) can be verified against the
     * current RelationStore and PropertyStore. Returns True if the condition
     * is provably satisfied, False if provably unsatisfied, Unknown otherwise.
     *
     * @param condition The condition expression to evaluate
     * @return Tribool indicating whether the condition is satisfied
     */
    Tribool evaluate_condition(const SymbolicExpr& condition) const;

    // --- Expression property queries ---

    /// Query whether the expression is positive (> 0).
    AssumptionTriboolResult is_positive(const SymbolicExpr& expr) const;

    /** @brief Queries positivity and propagates computation failures. */
    AssumptionTriboolResult is_positive_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is negative (< 0).
    AssumptionTriboolResult is_negative(const SymbolicExpr& expr) const;

    /** @brief Queries negativity and propagates computation failures. */
    AssumptionTriboolResult is_negative_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is non-negative (>= 0).
    AssumptionTriboolResult is_nonnegative(const SymbolicExpr& expr) const;

    /** @brief Queries non-negativity and propagates computation failures. */
    AssumptionTriboolResult is_nonnegative_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is real.
    AssumptionTriboolResult is_real(const SymbolicExpr& expr) const;

    /** @brief Queries real-domain membership and propagates computation failures. */
    AssumptionTriboolResult is_real_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is an integer.
    AssumptionTriboolResult is_integer(const SymbolicExpr& expr) const;

    /** @brief Queries integer-domain membership and propagates computation failures. */
    AssumptionTriboolResult is_integer_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is non-zero (!= 0).
    AssumptionTriboolResult is_nonzero(const SymbolicExpr& expr) const;

    /** @brief Queries nonzero status and propagates computation failures. */
    AssumptionTriboolResult is_nonzero_checked(const SymbolicExpr& expr) const;

    /**
     * @brief Serialize the entire AssumptionContext to a human-readable string.
     *
     * Format is line-oriented with SCOPE delimiters. Each scope's properties,
     * relations, and conditionals are serialized one per line. The output ends
     * with "END".
     *
     * @return Serialized string representation of all scopes.
     */
    std::string serialize() const;

    /**
     * @brief Deserialize a string into a fully reconstructed AssumptionContext.
     *
     * Parses the line-oriented format produced by serialize().
     *
     * @param data Serialized string (as produced by serialize()).
     * @return Reconstructed context or a parse/contradiction error.
     */
    static Result<AssumptionContext> deserialize(const std::string& data);

    /// 面向不可信输入的受检反序列化接口,通过 CasError 传递格式与约束诊断.
    static Result<AssumptionContext> deserialize_checked(const std::string& data);
    /// Checked deserialization using the caller's resource and cancellation policy.
    static Result<AssumptionContext> deserialize_checked(
        const std::string& data, ComputationContext& context);

    /**
     * @brief Query whether a symbol is continuous on a given interval (read-through all scopes).
     * @param symbol Symbol name to query.
     * @param interval Interval on which to check continuity.
     * @return True if continuous, False if not, Unknown if undetermined.
     */
    AssumptionTriboolResult is_continuous(const std::string& symbol, const Interval& interval) const;

    /** @brief Queries interval continuity and propagates computation failures. */
    AssumptionTriboolResult is_continuous_checked(
        const std::string& symbol,
        const Interval& interval) const;

    /**
     * @brief Query whether a symbol is differentiable on a given interval (read-through all scopes).
     * @param symbol Symbol name to query.
     * @param interval Interval on which to check differentiability.
     * @return True if differentiable, False if not, Unknown if undetermined.
     */
    AssumptionTriboolResult is_differentiable(const std::string& symbol, const Interval& interval) const;

    /** @brief Queries interval differentiability and propagates computation failures. */
    AssumptionTriboolResult is_differentiable_checked(
        const std::string& symbol,
        const Interval& interval) const;

    /**
     * @brief Query whether a symbol (matrix) is positive definite (read-through all scopes).
     * @param symbol Symbol name to query.
     * @return True if positive definite, False if known not, Unknown if undetermined.
     */
    AssumptionTriboolResult is_positive_definite(const std::string& symbol) const;

    /** @brief Queries positive definiteness and propagates computation failures. */
    AssumptionTriboolResult is_positive_definite_checked(const std::string& symbol) const;

    /**
     * @brief Query whether a symbol (matrix) is positive semidefinite (read-through all scopes).
     * @param symbol Symbol name to query.
     * @return True if positive semidefinite, False if known not, Unknown if undetermined.
     */
    AssumptionTriboolResult is_positive_semidefinite(const std::string& symbol) const;

    /** @brief Queries positive semidefiniteness and propagates computation failures. */
    AssumptionTriboolResult is_positive_semidefinite_checked(const std::string& symbol) const;

    /**
     * @brief Set the maximum recursion depth for inference queries.
     *
     * This depth limit is wired to the InferenceEngine used by QueryInterface.
     * @param depth Maximum depth (must be > 0). Default is 32.
     */
    void set_max_query_depth(int depth);

    /**
     * @brief Get the current maximum recursion depth for inference queries.
     * @return The configured maximum depth.
     */
    int get_max_query_depth() const;

private:
    static AssumptionContext deserialize_impl(const std::string& data);

    struct Scope {
        PropertyStore properties;
        RelationStore relations;
        std::vector<ConditionalAssumption> conditionals;
    };

    std::vector<Scope> scope_stack_;

    /// Maximum recursion depth for inference queries (default 32).
    int max_query_depth_ = 32;

    /// Generation counter incremented on every mutation (push/pop/assume).
    /// Used by QueryInterface to detect stale caches.
    mutable uint64_t cache_generation_ = 0;

public:
    /// Return the current cache generation counter.
    /// Incremented on every state mutation (push, pop, assume_domain, assume_sign, assume).
    uint64_t cache_generation() const { return cache_generation_; }
};

/**
 * @brief Declaration of an assumption to apply within a with_assumptions block.
 *
 * Represents a single domain, sign, or relational assumption that can be
 * applied to an AssumptionContext scope.
 */
class AssumptionDecl {
public:
    /// The type of assumption being declared.
    enum class Type { Domain, Sign, Relation };

    /// Construct a domain assumption declaration.
    static AssumptionDecl make_domain(const std::string& sym, Domain d) {
        return AssumptionDecl(DomainDeclaration{sym, d});
    }

    /// Construct a sign assumption declaration.
    static AssumptionDecl make_sign(const std::string& sym, Sign s) {
        return AssumptionDecl(SignDeclaration{sym, s});
    }

    /// Construct a relational assumption declaration.
    static AssumptionDecl make_relation(const SymbolicExpr& rel) {
        return AssumptionDecl(rel);
    }

    Type type() const noexcept {
        if (std::holds_alternative<DomainDeclaration>(payload_)) return Type::Domain;
        if (std::holds_alternative<SignDeclaration>(payload_)) return Type::Sign;
        return Type::Relation;
    }

    const std::string& symbol() const {
        if (const auto* domain = std::get_if<DomainDeclaration>(&payload_)) {
            return domain->symbol;
        }
        if (const auto* sign = std::get_if<SignDeclaration>(&payload_)) {
            return sign->symbol;
        }
        throw std::logic_error("relational assumption has no symbol field");
    }

    Domain domain() const {
        return std::get<DomainDeclaration>(payload_).domain;
    }

    Sign sign() const {
        return std::get<SignDeclaration>(payload_).sign;
    }

    const SymbolicExpr& relation() const {
        return std::get<SymbolicExpr>(payload_);
    }

private:
    struct DomainDeclaration {
        std::string symbol;
        Domain domain;
    };

    struct SignDeclaration {
        std::string symbol;
        Sign sign;
    };

    using Payload = std::variant<DomainDeclaration, SignDeclaration, SymbolicExpr>;

    explicit AssumptionDecl(DomainDeclaration declaration)
        : payload_(std::move(declaration)) {}
    explicit AssumptionDecl(SignDeclaration declaration)
        : payload_(std::move(declaration)) {}
    explicit AssumptionDecl(SymbolicExpr relation)
        : payload_(std::move(relation)) {}

    Payload payload_;
};

inline AssumptionVoidResult apply_assumption_decl_checked(
    AssumptionContext& ctx,
    const AssumptionDecl& decl) {
    switch (decl.type()) {
        case AssumptionDecl::Type::Domain:
            return ctx.assume_domain_checked(decl.symbol(), decl.domain());
        case AssumptionDecl::Type::Sign:
            return ctx.assume_sign_checked(decl.symbol(), decl.sign());
        case AssumptionDecl::Type::Relation:
            return ctx.assume_checked(decl.relation());
    }
    return AssumptionVoidResult::failure(
        CasErrc::InternalInvariant, "unknown assumption declaration type", "with_assumptions");
}

/**
 * @brief Execute a callable within a temporary assumption scope.
 *
 * Declaration and callable failures are returned as `CasError`. The temporary
 * scope is removed before every return.
 */
template<typename F>
auto with_assumptions(AssumptionContext& ctx,
                      const std::vector<AssumptionDecl>& decls,
                      F&& callable)
    -> std::enable_if_t<!std::is_void_v<decltype(callable())>, Result<decltype(callable())>>
{
    using ReturnT = decltype(callable());
    ctx.push();
    try {
        for (const auto& decl : decls) {
            auto applied = apply_assumption_decl_checked(ctx, decl);
            if (!applied) {
                auto popped = ctx.pop();
                if (!popped) return Result<ReturnT>::failure(popped.error());
                return Result<ReturnT>::failure(applied.error());
            }
        }
        auto result = callable();
        auto popped = ctx.pop();
        if (!popped) return Result<ReturnT>::failure(popped.error());
        return Result<ReturnT>::success(std::move(result));
    } catch (const std::bad_alloc&) {
        auto popped = ctx.pop();
        if (!popped) return Result<ReturnT>::failure(popped.error());
        return Result<ReturnT>::failure(
            CasErrc::ResourceLimit, "with_assumptions allocation failed", "with_assumptions");
    } catch (const std::exception& ex) {
        auto popped = ctx.pop();
        if (!popped) return Result<ReturnT>::failure(popped.error());
        return Result<ReturnT>::failure(
            CasErrc::InternalInvariant, ex.what(), "with_assumptions");
    }
}

template<typename F>
auto with_assumptions(AssumptionContext& ctx,
                      const std::vector<AssumptionDecl>& decls,
                      F&& callable)
    -> std::enable_if_t<std::is_void_v<decltype(callable())>, AssumptionVoidResult>
{
    ctx.push();
    try {
        for (const auto& decl : decls) {
            auto applied = apply_assumption_decl_checked(ctx, decl);
            if (!applied) {
                auto popped = ctx.pop();
                if (!popped) return AssumptionVoidResult::failure(popped.error());
                return AssumptionVoidResult::failure(applied.error());
            }
        }
        callable();
        auto popped = ctx.pop();
        if (!popped) return AssumptionVoidResult::failure(popped.error());
        return AssumptionVoidResult::success();
    } catch (const std::bad_alloc&) {
        auto popped = ctx.pop();
        if (!popped) return AssumptionVoidResult::failure(popped.error());
        return AssumptionVoidResult::failure(
            CasErrc::ResourceLimit, "with_assumptions allocation failed", "with_assumptions");
    } catch (const std::exception& ex) {
        auto popped = ctx.pop();
        if (!popped) return AssumptionVoidResult::failure(popped.error());
        return AssumptionVoidResult::failure(
            CasErrc::InternalInvariant, ex.what(), "with_assumptions");
    }
}

template<typename F>
auto with_assumptions(AssumptionContext& ctx,
                      std::initializer_list<AssumptionDecl> decls,
                      F&& callable)
    -> decltype(with_assumptions(
        ctx, std::vector<AssumptionDecl>(decls), std::forward<F>(callable)))
{
    return with_assumptions(
        ctx, std::vector<AssumptionDecl>(decls), std::forward<F>(callable));
}

} // namespace LMCAS
