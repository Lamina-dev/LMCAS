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
#include "symbolic_ast.hpp"
#include "assumption.hpp"

namespace lamina {

// Forward declaration — PropertyStore may not yet be compiled
class PropertyStore;

/**
 * @brief A stored relational constraint between two symbolic expressions.
 */
struct Relation {
    SymbolicExpr lhs;          ///< Left-hand side expression
    SymbolicExpr rhs;          ///< Right-hand side expression
    RelationalNode::Op op;     ///< Relational operator (GT, LT, GEQ, LEQ, NEQ, EQ)
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
     * If the LHS is a single VariableNode and the RHS is a NumberNode with value 0,
     * the corresponding sign property is declared on the PropertyStore.
     *
     * @param lhs Left-hand side expression
     * @param rhs Right-hand side expression
     * @param op  Relational operator
     * @param prop_store PropertyStore to notify for simple variable > 0 patterns
     */
    void add_relation(const SymbolicExpr& lhs, const SymbolicExpr& rhs,
                      RelationalNode::Op op, PropertyStore& prop_store);

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
                      RelationalNode::Op op) const;

    /**
     * @brief Clear all stored relations (used during scope pop).
     */
    void clear();

private:
    std::vector<Relation> relations_;
};

} // namespace lamina
