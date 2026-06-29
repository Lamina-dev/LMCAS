/**
 * @file relation_store.cpp
 * @brief Implementation of RelationStore for storing relational constraints.
 */

#include "relation_store.hpp"
#include "property_store.hpp"
#include <queue>

namespace lamina {

namespace {

/**
 * @brief Check if an operator participates in transitive closure.
 * Only GT and GEQ form transitive chains.
 */
bool is_transitive_op(RelationalNode::Op op) {
    return op == RelationalNode::Op::GT || op == RelationalNode::Op::GEQ;
}

/**
 * @brief Combine two transitive operators according to the combination rules:
 *   GT+GT→GT, GT+GEQ→GT, GEQ+GT→GT, GEQ+GEQ→GEQ.
 */
RelationalNode::Op combine_ops(RelationalNode::Op op1, RelationalNode::Op op2) {
    if (op1 == RelationalNode::Op::GEQ && op2 == RelationalNode::Op::GEQ) {
        return RelationalNode::Op::GEQ;
    }
    return RelationalNode::Op::GT;
}

/**
 * @brief Check structural equality of two expression roots.
 */
bool expr_equals(const SymbolicExpr& a, const SymbolicExpr& b) {
    if (!a.root && !b.root) return true;
    if (!a.root || !b.root) return false;
    return a.root->equals(*b.root);
}

} // anonymous namespace

void RelationStore::add_relation(const SymbolicExpr& lhs, const SymbolicExpr& rhs,
                                 RelationalNode::Op op, PropertyStore& prop_store) {
    // Store the relation regardless of pattern
    relations_.push_back(Relation{lhs, rhs, op});

    if (!lhs.root || !rhs.root) {
        return;
    }

    // Detect simple "variable op 0" pattern for sign property derivation.
    // LHS must be a single VariableNode and RHS must be a NumberNode with value 0.
    auto var_node = std::dynamic_pointer_cast<VariableNode>(lhs.root);
    auto num_node = std::dynamic_pointer_cast<NumberNode>(rhs.root);

    if (var_node && num_node && num_node->is_zero()) {
        // Map operator to sign property
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
                break;
        }
    } else {
        // Detect reversed "0 op variable" pattern
        detect_reversed_pattern(lhs, rhs, op, prop_store);
    }

    // Compute transitive closure for GT/GEQ relations
    if (is_transitive_op(op)) {
        compute_transitive_closure(relations_.back(), prop_store);
    }
}

void RelationStore::compute_transitive_closure(const Relation& new_rel, PropertyStore& prop_store) {
    // BFS queue: each entry is a deduced relation to explore further
    struct QueueEntry {
        SymbolicExpr lhs;
        SymbolicExpr rhs;
        RelationalNode::Op op;
    };

    std::queue<QueueEntry> bfs_queue;
    bfs_queue.push({new_rel.lhs, new_rel.rhs, new_rel.op});

    int deductions = 0;

    while (!bfs_queue.empty() && deductions < MAX_TRANSITIVE_DEDUCTIONS) {
        auto current = bfs_queue.front();
        bfs_queue.pop();

        // Snapshot the current relation count to avoid iterating over newly added relations.
        // We copy relevant relations to avoid invalidation from vector reallocation.
        const size_t relation_count = relations_.size();

        // Forward chaining: current is (A op B), find existing (B op2 C) → deduce (A combined_op C)
        for (size_t i = 0; i < relation_count && deductions < MAX_TRANSITIVE_DEDUCTIONS; ++i) {
            // Access by index each iteration since vector may have grown
            if (!is_transitive_op(relations_[i].op)) continue;

            // Check if current.rhs matches existing.lhs (forward chain)
            if (expr_equals(current.rhs, relations_[i].lhs)) {
                RelationalNode::Op combined = combine_ops(current.op, relations_[i].op);
                SymbolicExpr deduced_lhs = current.lhs;
                SymbolicExpr deduced_rhs = relations_[i].rhs;

                // Only add if not already stored
                if (!has_relation(deduced_lhs, deduced_rhs, combined)) {
                    relations_.push_back(Relation{deduced_lhs, deduced_rhs, combined});
                    ++deductions;

                    // Derive sign properties for the new deduced relation
                    if (deduced_lhs.root && deduced_rhs.root) {
                        detect_reversed_pattern(deduced_lhs, deduced_rhs, combined, prop_store);
                        // Also check "variable op 0" pattern
                        auto var_node = std::dynamic_pointer_cast<VariableNode>(deduced_lhs.root);
                        auto num_node = std::dynamic_pointer_cast<NumberNode>(deduced_rhs.root);
                        if (var_node && num_node && num_node->is_zero()) {
                            switch (combined) {
                                case RelationalNode::Op::GT:
                                    prop_store.declare_sign(var_node->name, Sign::Positive);
                                    break;
                                case RelationalNode::Op::GEQ:
                                    prop_store.declare_sign(var_node->name, Sign::NonNegative);
                                    break;
                                default:
                                    break;
                            }
                        }
                    }

                    // Enqueue for further BFS exploration
                    bfs_queue.push({deduced_lhs, deduced_rhs, combined});
                }
            }
        }

        // Backward chaining: current is (A op B), find existing (C op2 A) → deduce (C combined_op B)
        for (size_t i = 0; i < relation_count && deductions < MAX_TRANSITIVE_DEDUCTIONS; ++i) {
            if (!is_transitive_op(relations_[i].op)) continue;

            // Check if existing.rhs matches current.lhs (backward chain)
            if (expr_equals(relations_[i].rhs, current.lhs)) {
                RelationalNode::Op combined = combine_ops(relations_[i].op, current.op);
                SymbolicExpr deduced_lhs = relations_[i].lhs;
                SymbolicExpr deduced_rhs = current.rhs;

                // Only add if not already stored
                if (!has_relation(deduced_lhs, deduced_rhs, combined)) {
                    relations_.push_back(Relation{deduced_lhs, deduced_rhs, combined});
                    ++deductions;

                    // Derive sign properties for the new deduced relation
                    if (deduced_lhs.root && deduced_rhs.root) {
                        detect_reversed_pattern(deduced_lhs, deduced_rhs, combined, prop_store);
                        // Also check "variable op 0" pattern
                        auto var_node = std::dynamic_pointer_cast<VariableNode>(deduced_lhs.root);
                        auto num_node = std::dynamic_pointer_cast<NumberNode>(deduced_rhs.root);
                        if (var_node && num_node && num_node->is_zero()) {
                            switch (combined) {
                                case RelationalNode::Op::GT:
                                    prop_store.declare_sign(var_node->name, Sign::Positive);
                                    break;
                                case RelationalNode::Op::GEQ:
                                    prop_store.declare_sign(var_node->name, Sign::NonNegative);
                                    break;
                                default:
                                    break;
                            }
                        }
                    }

                    // Enqueue for further BFS exploration
                    bfs_queue.push({deduced_lhs, deduced_rhs, combined});
                }
            }
        }
    }
}

void RelationStore::detect_reversed_pattern(const SymbolicExpr& lhs, const SymbolicExpr& rhs,
                                            RelationalNode::Op op, PropertyStore& prop_store) {
    // LHS must be a NumberNode with value 0
    auto num_node = std::dynamic_pointer_cast<NumberNode>(lhs.root);
    if (!num_node || !num_node->is_zero()) {
        return;
    }

    // RHS must be a single VariableNode
    auto var_node = std::dynamic_pointer_cast<VariableNode>(rhs.root);
    if (!var_node) {
        return;
    }

    // Reversed semantics:
    //   0 LT  var → 0 < var → var > 0 → Positive
    //   0 GT  var → 0 > var → var < 0 → Negative
    //   0 GEQ var → 0 >= var → var <= 0 → NonPositive
    //   0 LEQ var → 0 <= var → var >= 0 → NonNegative
    //   0 NEQ var → 0 != var → NonZero
    switch (op) {
        case RelationalNode::Op::LT:
            prop_store.declare_sign(var_node->name, Sign::Positive);
            break;
        case RelationalNode::Op::GT:
            prop_store.declare_sign(var_node->name, Sign::Negative);
            break;
        case RelationalNode::Op::GEQ:
            prop_store.declare_sign(var_node->name, Sign::NonPositive);
            break;
        case RelationalNode::Op::LEQ:
            prop_store.declare_sign(var_node->name, Sign::NonNegative);
            break;
        case RelationalNode::Op::NEQ:
            prop_store.declare_sign(var_node->name, Sign::NonZero);
            break;
        case RelationalNode::Op::EQ:
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
