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
#include "query_interface.hpp"
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

namespace lamina {

// Forward declarations
class InferenceEngine;
struct Interval;

using AssumptionVoidResult = Result<void>;
using AssumptionTriboolResult = Result<Tribool>;

/**
 * @brief Scoped assumption context with push/pop management.
 *
 * Maintains a stack of scopes, each containing a PropertyStore and RelationStore.
 * Queries use read-through semantics: searching from the top scope down to root,
 * returning the first scope that has a declaration for the queried symbol.
 * Child scope declarations shadow parent declarations without modifying them.
 * Popping a scope restores all query results to their pre-push values.
 */
class LAMINA_API AssumptionContext {
public:
    AssumptionContext();

    // --- Scope management ---

    /// Push a new scope onto the stack.
    void push();

    /// Pop the current scope. Throws std::runtime_error if at root scope.
    void pop();

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

    // --- Convenience declaration API ---

    /// Declare domain for a variable in the current scope.
    /// Throws std::invalid_argument if variable name is empty.
    void assume_domain(const std::string& variable, Domain domain);

    /// Checked domain declaration for migration away from exception-only APIs.
    AssumptionVoidResult assume_domain_checked(const std::string& variable, Domain domain);

    /// Declare sign for a variable in the current scope.
    /// Throws std::invalid_argument if variable name is empty.
    void assume_sign(const std::string& variable, Sign sign);

    /// Checked sign declaration for migration away from exception-only APIs.
    AssumptionVoidResult assume_sign_checked(const std::string& variable, Sign sign);

    /// Store a relational constraint. The expression's root must be a RelationalNode.
    /// Throws std::invalid_argument if expression is null/empty or root is not RelationalNode.
    void assume(const SymbolicExpr& relation);

    /// Checked relation declaration for migration away from exception-only APIs.
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
     * @throws std::invalid_argument if condition or conclusion is not relational
     */
    void assume_conditional(const SymbolicExpr& condition, const SymbolicExpr& conclusion);

    /// Checked conditional declaration for migration away from exception-only APIs.
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

    // --- Convenience query API (delegates to QueryInterface) ---

    /// Query whether the expression is positive (> 0).
    Tribool is_positive(const SymbolicExpr& expr) const;

    /// Checked positive query for migration away from exception-only APIs.
    AssumptionTriboolResult is_positive_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is negative (< 0).
    Tribool is_negative(const SymbolicExpr& expr) const;

    /// Checked negative query for migration away from exception-only APIs.
    AssumptionTriboolResult is_negative_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is non-negative (>= 0).
    Tribool is_nonnegative(const SymbolicExpr& expr) const;

    /// Checked non-negative query for migration away from exception-only APIs.
    AssumptionTriboolResult is_nonnegative_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is real.
    Tribool is_real(const SymbolicExpr& expr) const;

    /// Checked real-domain query for migration away from exception-only APIs.
    AssumptionTriboolResult is_real_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is an integer.
    Tribool is_integer(const SymbolicExpr& expr) const;

    /// Checked integer-domain query for migration away from exception-only APIs.
    AssumptionTriboolResult is_integer_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is non-zero (!= 0).
    Tribool is_nonzero(const SymbolicExpr& expr) const;

    /// Checked non-zero query for migration away from exception-only APIs.
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
     * Parses the line-oriented format produced by serialize(). Throws
     * std::invalid_argument with line number on malformed input.
     *
     * @param data Serialized string (as produced by serialize()).
     * @return Reconstructed AssumptionContext.
     * @throws std::invalid_argument on malformed input with line number info.
     */
    static AssumptionContext deserialize(const std::string& data);

    /// Checked deserialization for untrusted input. Malformed or contradictory
    /// data returns CasError instead of throwing.
    static Result<AssumptionContext> deserialize_checked(const std::string& data);

    /**
     * @brief Query whether a symbol is continuous on a given interval (read-through all scopes).
     * @param symbol Symbol name to query.
     * @param interval Interval on which to check continuity.
     * @return True if continuous, False if not, Unknown if undetermined.
     */
    Tribool is_continuous(const std::string& symbol, const Interval& interval) const;

    /// Checked continuity query for migration away from bare Tribool APIs.
    AssumptionTriboolResult is_continuous_checked(
        const std::string& symbol,
        const Interval& interval) const;

    /**
     * @brief Query whether a symbol is differentiable on a given interval (read-through all scopes).
     * @param symbol Symbol name to query.
     * @param interval Interval on which to check differentiability.
     * @return True if differentiable, False if not, Unknown if undetermined.
     */
    Tribool is_differentiable(const std::string& symbol, const Interval& interval) const;

    /// Checked differentiability query for migration away from bare Tribool APIs.
    AssumptionTriboolResult is_differentiable_checked(
        const std::string& symbol,
        const Interval& interval) const;

    /**
     * @brief Query whether a symbol (matrix) is positive definite (read-through all scopes).
     * @param symbol Symbol name to query.
     * @return True if positive definite, False if known not, Unknown if undetermined.
     */
    Tribool is_positive_definite(const std::string& symbol) const;

    /// Checked positive-definite query for migration away from bare Tribool APIs.
    AssumptionTriboolResult is_positive_definite_checked(const std::string& symbol) const;

    /**
     * @brief Query whether a symbol (matrix) is positive semidefinite (read-through all scopes).
     * @param symbol Symbol name to query.
     * @return True if positive semidefinite, False if known not, Unknown if undetermined.
     */
    Tribool is_positive_semidefinite(const std::string& symbol) const;

    /// Checked positive-semidefinite query for migration away from bare Tribool APIs.
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

// with_assumptions: scoped assumption application with RAII safety

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
 * Pushes a new scope on the context, applies all declarations from the
 * provided set, invokes the callable, pops the scope, and returns the result.
 * Exception-safe: the scope is popped even if the callable throws.
 *
 * @tparam F Callable type (must be invocable with no arguments)
 * @param ctx AssumptionContext to operate on
 * @param decls Set of assumption declarations to apply in the new scope
 * @param callable Function to invoke under the temporary assumptions
 * @return The result of invoking callable()
 *
 * @note Supports domain, sign, and relational assumptions (Req 17.4).
 * @note Preserves context depth before and after, even on exception (Req 17.2, 17.3).
 */
template<typename F>
auto with_assumptions(AssumptionContext& ctx,
                      std::initializer_list<AssumptionDecl> decls,
                      F&& callable)
    -> std::enable_if_t<!std::is_void_v<decltype(callable())>, decltype(callable())>
{
    ctx.push();
    try {
        for (const auto& decl : decls) {
            switch (decl.type()) {
                case AssumptionDecl::Type::Domain:
                    ctx.assume_domain(decl.symbol(), decl.domain());
                    break;
                case AssumptionDecl::Type::Sign:
                    ctx.assume_sign(decl.symbol(), decl.sign());
                    break;
                case AssumptionDecl::Type::Relation:
                    ctx.assume(decl.relation());
                    break;
            }
        }
        auto result = callable();
        ctx.pop();
        return result;
    } catch (...) {
        ctx.pop();
        throw;
    }
}

/**
 * @brief Checked scoped assumption application for non-void callables.
 *
 * Applies declarations with checked APIs and returns CasError on declaration or
 * callable failure. The temporary scope is always popped before returning.
 */
template<typename F>
auto with_assumptions_checked(AssumptionContext& ctx,
                              const std::vector<AssumptionDecl>& decls,
                              F&& callable)
    -> std::enable_if_t<!std::is_void_v<decltype(callable())>, Result<decltype(callable())>>
{
    using ReturnT = decltype(callable());
    ctx.push();
    try {
        for (const auto& decl : decls) {
            auto applied = apply_assumption_decl_checked(ctx, decl);
            if (!applied.has_value()) {
                ctx.pop();
                return Result<ReturnT>::failure(applied.error());
            }
        }
        auto result = callable();
        ctx.pop();
        return Result<ReturnT>::success(std::move(result));
    } catch (const std::bad_alloc&) {
        ctx.pop();
        return Result<ReturnT>::failure(
            CasErrc::ResourceLimit, "with_assumptions allocation failed", "with_assumptions");
    } catch (const std::exception& ex) {
        ctx.pop();
        return Result<ReturnT>::failure(
            CasErrc::InternalInvariant, ex.what(), "with_assumptions");
    }
}

/**
 * @brief Checked scoped assumption application for void callables.
 */
template<typename F>
auto with_assumptions_checked(AssumptionContext& ctx,
                              const std::vector<AssumptionDecl>& decls,
                              F&& callable)
    -> std::enable_if_t<std::is_void_v<decltype(callable())>, AssumptionVoidResult>
{
    ctx.push();
    try {
        for (const auto& decl : decls) {
            auto applied = apply_assumption_decl_checked(ctx, decl);
            if (!applied.has_value()) {
                ctx.pop();
                return AssumptionVoidResult::failure(applied.error());
            }
        }
        callable();
        ctx.pop();
        return AssumptionVoidResult::success();
    } catch (const std::bad_alloc&) {
        ctx.pop();
        return AssumptionVoidResult::failure(
            CasErrc::ResourceLimit, "with_assumptions allocation failed", "with_assumptions");
    } catch (const std::exception& ex) {
        ctx.pop();
        return AssumptionVoidResult::failure(
            CasErrc::InternalInvariant, ex.what(), "with_assumptions");
    }
}

template<typename F>
auto with_assumptions_checked(AssumptionContext& ctx,
                              std::initializer_list<AssumptionDecl> decls,
                              F&& callable)
    -> decltype(with_assumptions_checked(
        ctx, std::vector<AssumptionDecl>(decls), std::forward<F>(callable)))
{
    return with_assumptions_checked(
        ctx, std::vector<AssumptionDecl>(decls), std::forward<F>(callable));
}

/**
 * @brief Execute a void-returning callable within a temporary assumption scope.
 *
 * Overload for callables that return void. Pushes a new scope, applies
 * declarations, invokes the callable, and pops the scope.
 * Exception-safe: the scope is popped even if the callable throws.
 *
 * @tparam F Callable type (must be invocable with no arguments, returning void)
 * @param ctx AssumptionContext to operate on
 * @param decls Set of assumption declarations to apply in the new scope
 * @param callable Function to invoke under the temporary assumptions
 */
template<typename F>
auto with_assumptions(AssumptionContext& ctx,
                      std::initializer_list<AssumptionDecl> decls,
                      F&& callable)
    -> std::enable_if_t<std::is_void_v<decltype(callable())>, void>
{
    ctx.push();
    try {
        for (const auto& decl : decls) {
            switch (decl.type()) {
                case AssumptionDecl::Type::Domain:
                    ctx.assume_domain(decl.symbol(), decl.domain());
                    break;
                case AssumptionDecl::Type::Sign:
                    ctx.assume_sign(decl.symbol(), decl.sign());
                    break;
                case AssumptionDecl::Type::Relation:
                    ctx.assume(decl.relation());
                    break;
            }
        }
        callable();
        ctx.pop();
    } catch (...) {
        ctx.pop();
        throw;
    }
}

/**
 * @brief Execute a callable within a temporary assumption scope (vector overload).
 *
 * Same as the initializer_list overload but accepts a std::vector of declarations.
 *
 * @tparam F Callable type (must be invocable with no arguments)
 * @param ctx AssumptionContext to operate on
 * @param decls Vector of assumption declarations to apply in the new scope
 * @param callable Function to invoke under the temporary assumptions
 * @return The result of invoking callable()
 */
template<typename F>
auto with_assumptions(AssumptionContext& ctx,
                      const std::vector<AssumptionDecl>& decls,
                      F&& callable)
    -> std::enable_if_t<!std::is_void_v<decltype(callable())>, decltype(callable())>
{
    ctx.push();
    try {
        for (const auto& decl : decls) {
            switch (decl.type()) {
                case AssumptionDecl::Type::Domain:
                    ctx.assume_domain(decl.symbol(), decl.domain());
                    break;
                case AssumptionDecl::Type::Sign:
                    ctx.assume_sign(decl.symbol(), decl.sign());
                    break;
                case AssumptionDecl::Type::Relation:
                    ctx.assume(decl.relation());
                    break;
            }
        }
        auto result = callable();
        ctx.pop();
        return result;
    } catch (...) {
        ctx.pop();
        throw;
    }
}

/**
 * @brief Execute a void-returning callable within a temporary assumption scope (vector overload).
 *
 * Same as the initializer_list overload but accepts a std::vector of declarations.
 *
 * @tparam F Callable type (must be invocable with no arguments, returning void)
 * @param ctx AssumptionContext to operate on
 * @param decls Vector of assumption declarations to apply in the new scope
 * @param callable Function to invoke under the temporary assumptions
 */
template<typename F>
auto with_assumptions(AssumptionContext& ctx,
                      const std::vector<AssumptionDecl>& decls,
                      F&& callable)
    -> std::enable_if_t<std::is_void_v<decltype(callable())>, void>
{
    ctx.push();
    try {
        for (const auto& decl : decls) {
            switch (decl.type()) {
                case AssumptionDecl::Type::Domain:
                    ctx.assume_domain(decl.symbol(), decl.domain());
                    break;
                case AssumptionDecl::Type::Sign:
                    ctx.assume_sign(decl.symbol(), decl.sign());
                    break;
                case AssumptionDecl::Type::Relation:
                    ctx.assume(decl.relation());
                    break;
            }
        }
        callable();
        ctx.pop();
    } catch (...) {
        ctx.pop();
        throw;
    }
}

} // namespace lamina
