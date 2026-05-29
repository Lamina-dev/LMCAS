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
#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_set>
#include <optional>

namespace lamina {

// Forward declarations
class InferenceEngine;
struct Interval;

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

    /// Declare sign for a variable in the current scope.
    /// Throws std::invalid_argument if variable name is empty.
    void assume_sign(const std::string& variable, Sign sign);

    /// Store a relational constraint. The expression's root must be a RelationalNode.
    /// Throws std::invalid_argument if expression is null/empty or root is not RelationalNode.
    void assume(const SymbolicExpr& relation);

    // --- Convenience query API (delegates to QueryInterface) ---

    /// Query whether the expression is positive (> 0).
    /// Throws std::invalid_argument if expression has null root.
    Tribool is_positive(const SymbolicExpr& expr) const;

    /// Query whether the expression is negative (< 0).
    /// Throws std::invalid_argument if expression has null root.
    Tribool is_negative(const SymbolicExpr& expr) const;

    /// Query whether the expression is non-negative (>= 0).
    /// Throws std::invalid_argument if expression has null root.
    Tribool is_nonnegative(const SymbolicExpr& expr) const;

    /// Query whether the expression is real.
    /// Throws std::invalid_argument if expression has null root.
    Tribool is_real(const SymbolicExpr& expr) const;

    /// Query whether the expression is an integer.
    /// Throws std::invalid_argument if expression has null root.
    Tribool is_integer(const SymbolicExpr& expr) const;

    /// Query whether the expression is non-zero (!= 0).
    /// Throws std::invalid_argument if expression has null root.
    Tribool is_nonzero(const SymbolicExpr& expr) const;

private:
    struct Scope {
        PropertyStore properties;
        RelationStore relations;
    };

    std::vector<Scope> scope_stack_;
};

} // namespace lamina
