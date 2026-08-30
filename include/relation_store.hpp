/**
 * @file relation_store.hpp
 * @brief RelationStore class for storing relational constraints between symbolic expressions.
 *
 * Stores inequality relations (GT, LT, GEQ, LEQ, NEQ) between SymbolicExpr instances.
 * When a simple "variable op 0" pattern is detected, the corresponding sign property
 * is automatically propagated to the PropertyStore.
 */
#pragma once

#include <vector>
#include <string>
#include <memory>
#include "symbolic.hpp"
#include "assumption.hpp"
#include "result.hpp"

namespace lamina {

// Forward declaration — PropertyStore may not yet be compiled
class PropertyStore;

using RelationStoreResult = Result<void>;

/**
 * @brief A stored relational constraint between two symbolic expressions.
 */
struct Relation {
    SymbolicExpr lhs;          ///< Left-hand side expression
    SymbolicExpr rhs;          ///< Right-hand side expression
    RelationOp op;             ///< Relational operator (GT, LT, GEQ, LEQ, NEQ, EQ)
};

/**
 * @brief Stores relational constraints and derives sign properties for simple patterns.
 *
 * When a relation of the form `variable op 0` is added, the RelationStore notifies
 * the PropertyStore to mark the variable with the corresponding sign:
 *   - GT  → Positive
 *   - GEQ → NonNegative
 *   - LT  → Negative
 *   - LEQ → NonPositive
 *   - NEQ → NonZero
 *
 * Composite relations (multi-variable LHS or non-zero RHS) are stored without
 * decomposition for later use by the InferenceEngine.
 */
class LAMINA_API RelationStore {
public:
    /**
     * @brief Store a relation and optionally derive sign properties.
     *
     * If the relation compares one variable with exact zero, the corresponding
     * sign property is declared on the PropertyStore.
     *
     * @param lhs Left-hand side expression
     * @param rhs Right-hand side expression
     * @param op  Relational operator
     * @param prop_store PropertyStore to notify for simple variable > 0 patterns
     */
    RelationStoreResult add_relation(
        const SymbolicExpr& lhs, const SymbolicExpr& rhs,
        RelationOp op, PropertyStore& prop_store);

    /**
     * @brief Checked relation insertion with explicit failure reporting.
     *
     * Applies the relation plus derived property declarations transactionally.
     * On failure, neither this store nor the PropertyStore is modified.
     */
    RelationStoreResult add_relation_checked(const SymbolicExpr& lhs,
                                             const SymbolicExpr& rhs,
                                             RelationOp op,
                                             PropertyStore& prop_store);

    /**
     * @brief Retrieve all stored relations.
     * @return Const reference to the vector of stored relations
     */
    const std::vector<Relation>& get_relations() const;

    /**
     * @brief Check if a specific relation is stored.
     *
     * Compares LHS, RHS, and operator using structural equality of the AST nodes.
     *
     * @param lhs Left-hand side expression
     * @param rhs Right-hand side expression
     * @param op  Relational operator
     * @return true if the relation is found in the store
     */
    bool has_relation(const SymbolicExpr& lhs, const SymbolicExpr& rhs,
                      RelationOp op) const;

    /**
     * @brief Clear all stored relations (used during scope pop).
     */
    void clear();

private:
    std::vector<Relation> relations_;

    RelationStoreResult add_relation_unchecked(
        const SymbolicExpr& lhs,
        const SymbolicExpr& rhs,
        RelationOp op,
        PropertyStore& prop_store);

    /// Maximum number of new relations deduced per add_relation call via transitive closure.
    static constexpr int MAX_TRANSITIVE_DEDUCTIONS = 64;

    /**
     * @brief Detect reversed "0 op variable" patterns and derive sign properties.
     *
     * When exact zero appears on the left and a variable on the right, the
     * operator semantics are reversed to derive the variable's sign:
     *   - 0 LT  var → var is Positive   (0 < var means var > 0)
     *   - 0 GT  var → var is Negative   (0 > var means var < 0)
     *   - 0 GEQ var → var is NonPositive (0 >= var means var <= 0)
     *   - 0 LEQ var → var is NonNegative (0 <= var means var >= 0)
     *   - 0 NEQ var → var is NonZero
     *
     * @param lhs Left-hand side expression (expected to be zero)
     * @param rhs Right-hand side expression (expected to be a variable)
     * @param op  Relational operator
     * @param prop_store PropertyStore to update with derived sign
     */
    RelationStoreResult detect_reversed_pattern(
        const SymbolicExpr& lhs, const SymbolicExpr& rhs,
        RelationOp op, PropertyStore& prop_store);

    /**
     * @brief Compute transitive closure after adding a new relation.
     *
     * BFS from the new relation, combining GT/GEQ operators transitively:
     *   GT+GT→GT, GT+GEQ→GT, GEQ+GT→GT, GEQ+GEQ→GEQ.
     * Only GT and GEQ participate in transitive closure.
     * Stops after MAX_TRANSITIVE_DEDUCTIONS new relations are deduced.
     *
     * @param new_rel The newly added relation that triggers closure computation
     * @param prop_store PropertyStore for sign derivation of deduced relations
     */
    RelationStoreResult compute_transitive_closure(
        const Relation& new_rel, PropertyStore& prop_store);
};

} // namespace lamina
