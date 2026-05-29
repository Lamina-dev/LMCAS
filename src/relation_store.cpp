/**
 * @file relation_store.cpp
 * @brief Implementation of RelationStore for storing relational constraints.
 */

#include "relation_store.hpp"
#include "property_store.hpp"

namespace lamina {

void RelationStore::add_relation(const SymbolicExpr& lhs, const SymbolicExpr& rhs,
                                 RelationalNode::Op op, PropertyStore& prop_store) {
    // Store the relation regardless of pattern
    relations_.push_back(Relation{lhs, rhs, op});

    // Detect simple "variable op 0" pattern for sign property derivation.
    // LHS must be a single VariableNode and RHS must be a NumberNode with value 0.
    if (!lhs.root || !rhs.root) {
        return;
    }

    auto var_node = std::dynamic_pointer_cast<VariableNode>(lhs.root);
    if (!var_node) {
        return;  // LHS is not a single variable — composite relation, store without decomposition
    }

    auto num_node = std::dynamic_pointer_cast<NumberNode>(rhs.root);
    if (!num_node || !num_node->is_zero()) {
        return;  // RHS is not zero — store without decomposition
    }

    // Map operator to sign property
    // GT  → Positive
    // GEQ → NonNegative
    // LT  → Negative
    // LEQ → NonPositive
    // NEQ → NonZero
    switch (op) {
        case RelationalNode::Op::GT:
            prop_store.declare_sign(var_node->name, Sign::Positive);
            break;
        case RelationalNode::Op::GEQ:
            prop_store.declare_sign(var_node->name, Sign::NonNegative);
            break;
        case RelationalNode::Op::LT:
            prop_store.declare_sign(var_node->name, Sign::Negative);
            break;
        case RelationalNode::Op::LEQ:
            prop_store.declare_sign(var_node->name, Sign::NonPositive);
            break;
        case RelationalNode::Op::NEQ:
            prop_store.declare_sign(var_node->name, Sign::NonZero);
            break;
        case RelationalNode::Op::EQ:
            // EQ with zero could imply Zero sign, but the design only specifies
            // GT, GEQ, LT, LEQ, NEQ mappings. Store without sign derivation.
            break;
    }
}

const std::vector<Relation>& RelationStore::get_relations() const {
    return relations_;
}

bool RelationStore::has_relation(const SymbolicExpr& lhs, const SymbolicExpr& rhs,
                                 RelationalNode::Op op) const {
    for (const auto& rel : relations_) {
        if (rel.op != op) {
            continue;
        }
        // Compare LHS structurally
        if (!rel.lhs.root && !lhs.root) {
            // Both null — match on LHS
        } else if (!rel.lhs.root || !lhs.root) {
            continue;  // One null, one not — no match
        } else if (!rel.lhs.root->equals(*lhs.root)) {
            continue;
        }
        // Compare RHS structurally
        if (!rel.rhs.root && !rhs.root) {
            return true;  // Both null — full match
        } else if (!rel.rhs.root || !rhs.root) {
            continue;  // One null, one not — no match
        } else if (rel.rhs.root->equals(*rhs.root)) {
            return true;
        }
    }
    return false;
}

void RelationStore::clear() {
    relations_.clear();
}

} // namespace lamina
