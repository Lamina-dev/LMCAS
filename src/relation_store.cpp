/**
 * @file relation_store.cpp
 * @brief Implementation of RelationStore for storing relational constraints.
 */

#include "relation_store.hpp"
#include "symbolic_ast.hpp"
#include "property_store.hpp"
#include <queue>

namespace lamina {

namespace {

/**
 * @brief Check if an operator participates in transitive closure.
 * Only GT and GEQ form transitive chains.
 */
bool is_transitive_op(RelationOp op) {
    return op == RelationOp::GT || op == RelationOp::GEQ;
}

/**
 * @brief Combine two transitive operators according to the combination rules:
 *   GT+GT→GT, GT+GEQ→GT, GEQ+GT→GT, GEQ+GEQ→GEQ.
 */
RelationOp combine_ops(RelationOp op1, RelationOp op2) {
    if (op1 == RelationOp::GEQ && op2 == RelationOp::GEQ) {
        return RelationOp::GEQ;
    }
    return RelationOp::GT;
}

/**
 * @brief Check structural equality of two expression roots.
 */
bool expr_equals(const SymbolicExpr& a, const SymbolicExpr& b) {
    if (!lamina::detail::node(a) && !lamina::detail::node(b)) return true;
    if (!lamina::detail::node(a) || !lamina::detail::node(b)) return false;
    return lamina::detail::node(a)->equals(*lamina::detail::node(b));
}

void require_sign_declaration(PropertyStore& store,
                              const std::string& symbol,
                              Sign sign) {
    auto result = store.declare_sign_checked(symbol, sign);
    if (result) {
        return;
    }
    if (result.error().code == CasErrc::InvalidArgument) {
        throw std::invalid_argument(result.error().message);
    }
    if (result.error().code == CasErrc::ResourceLimit) {
        throw std::bad_alloc();
    }
    throw std::runtime_error(result.error().message);
}

} // anonymous namespace

void RelationStore::add_relation_unchecked(const SymbolicExpr& lhs,
                                           const SymbolicExpr& rhs,
                                           RelationOp op,
                                           PropertyStore& prop_store) {
    // Store the relation regardless of pattern
    relations_.push_back(Relation{lhs, rhs, op});

    if (!lamina::detail::node(lhs) || !lamina::detail::node(rhs)) {
        return;
    }

    // Detect simple "variable op 0" pattern for sign property derivation.
    // LHS must be a single VariableNode and RHS must be a NumberNode with value 0.
    auto var_node = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(lhs));
    auto num_node = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(rhs));

    if (var_node && num_node && num_node->is_zero()) {
        // Map operator to sign property
        switch (op) {
            case RelationOp::GT:
                require_sign_declaration(prop_store, var_node->name(), Sign::Positive);
                break;
            case RelationOp::GEQ:
                require_sign_declaration(prop_store, var_node->name(), Sign::NonNegative);
                break;
            case RelationOp::LT:
                require_sign_declaration(prop_store, var_node->name(), Sign::Negative);
                break;
            case RelationOp::LEQ:
                require_sign_declaration(prop_store, var_node->name(), Sign::NonPositive);
                break;
            case RelationOp::NEQ:
                require_sign_declaration(prop_store, var_node->name(), Sign::NonZero);
                break;
            case RelationOp::EQ:
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

void RelationStore::add_relation(const SymbolicExpr& lhs,
                                 const SymbolicExpr& rhs,
                                 RelationOp op,
                                 PropertyStore& prop_store) {
    auto result = add_relation_checked(lhs, rhs, op, prop_store);
    if (result) {
        return;
    }

    const auto& error = result.error();
    if (error.code == CasErrc::InvalidArgument) {
        throw std::invalid_argument(error.message);
    }
    if (error.code == CasErrc::ResourceLimit) {
        throw std::bad_alloc();
    }
    throw std::runtime_error(error.message);
}

RelationStoreResult RelationStore::add_relation_checked(
    const SymbolicExpr& lhs,
    const SymbolicExpr& rhs,
    RelationOp op,
    PropertyStore& prop_store) {
    if (!lamina::detail::node(lhs)) {
        return RelationStoreResult::failure(
            CasErrc::InvalidArgument, "relation lhs must not be null", "add_relation");
    }
    if (!lamina::detail::node(rhs)) {
        return RelationStoreResult::failure(
            CasErrc::InvalidArgument, "relation rhs must not be null", "add_relation");
    }

    try {
        RelationStore relation_candidate = *this;
        PropertyStore property_candidate = prop_store;
        relation_candidate.add_relation_unchecked(lhs, rhs, op, property_candidate);
        *this = std::move(relation_candidate);
        prop_store = std::move(property_candidate);
    } catch (const std::bad_alloc&) {
        return RelationStoreResult::failure(
            CasErrc::ResourceLimit, "relation-store allocation failed", "add_relation");
    } catch (const std::invalid_argument& ex) {
        return RelationStoreResult::failure(CasErrc::InvalidArgument, ex.what(), "add_relation");
    } catch (const std::exception& ex) {
        return RelationStoreResult::failure(
            CasErrc::InternalInvariant, ex.what(), "add_relation");
    }

    return RelationStoreResult::success();
}

void RelationStore::compute_transitive_closure(const Relation& new_rel, PropertyStore& prop_store) {
    // BFS queue: each entry is a deduced relation to explore further
    struct QueueEntry {
        SymbolicExpr lhs;
        SymbolicExpr rhs;
        RelationOp op;
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
                RelationOp combined = combine_ops(current.op, relations_[i].op);
                SymbolicExpr deduced_lhs = current.lhs;
                SymbolicExpr deduced_rhs = relations_[i].rhs;

                // Only add if not already stored
                if (!has_relation(deduced_lhs, deduced_rhs, combined)) {
                    relations_.push_back(Relation{deduced_lhs, deduced_rhs, combined});
                    ++deductions;

                    // Derive sign properties for the new deduced relation
                    if (lamina::detail::node(deduced_lhs) && lamina::detail::node(deduced_rhs)) {
                        detect_reversed_pattern(deduced_lhs, deduced_rhs, combined, prop_store);
                        // Also check "variable op 0" pattern
                        auto var_node = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(deduced_lhs));
                        auto num_node = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(deduced_rhs));
                        if (var_node && num_node && num_node->is_zero()) {
                            switch (combined) {
                                case RelationOp::GT:
                                    require_sign_declaration(prop_store, var_node->name(), Sign::Positive);
                                    break;
                                case RelationOp::GEQ:
                                    require_sign_declaration(prop_store, var_node->name(), Sign::NonNegative);
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
                RelationOp combined = combine_ops(relations_[i].op, current.op);
                SymbolicExpr deduced_lhs = relations_[i].lhs;
                SymbolicExpr deduced_rhs = current.rhs;

                // Only add if not already stored
                if (!has_relation(deduced_lhs, deduced_rhs, combined)) {
                    relations_.push_back(Relation{deduced_lhs, deduced_rhs, combined});
                    ++deductions;

                    // Derive sign properties for the new deduced relation
                    if (lamina::detail::node(deduced_lhs) && lamina::detail::node(deduced_rhs)) {
                        detect_reversed_pattern(deduced_lhs, deduced_rhs, combined, prop_store);
                        // Also check "variable op 0" pattern
                        auto var_node = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(deduced_lhs));
                        auto num_node = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(deduced_rhs));
                        if (var_node && num_node && num_node->is_zero()) {
                            switch (combined) {
                                case RelationOp::GT:
                                    require_sign_declaration(prop_store, var_node->name(), Sign::Positive);
                                    break;
                                case RelationOp::GEQ:
                                    require_sign_declaration(prop_store, var_node->name(), Sign::NonNegative);
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
                                            RelationOp op, PropertyStore& prop_store) {
    // LHS must be a NumberNode with value 0
    auto num_node = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(lhs));
    if (!num_node || !num_node->is_zero()) {
        return;
    }

    // RHS must be a single VariableNode
    auto var_node = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(rhs));
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
        case RelationOp::LT:
            require_sign_declaration(prop_store, var_node->name(), Sign::Positive);
            break;
        case RelationOp::GT:
            require_sign_declaration(prop_store, var_node->name(), Sign::Negative);
            break;
        case RelationOp::GEQ:
            require_sign_declaration(prop_store, var_node->name(), Sign::NonPositive);
            break;
        case RelationOp::LEQ:
            require_sign_declaration(prop_store, var_node->name(), Sign::NonNegative);
            break;
        case RelationOp::NEQ:
            require_sign_declaration(prop_store, var_node->name(), Sign::NonZero);
            break;
        case RelationOp::EQ:
            break;
    }
}

const std::vector<Relation>& RelationStore::get_relations() const {
    return relations_;
}

bool RelationStore::has_relation(const SymbolicExpr& lhs, const SymbolicExpr& rhs,
                                 RelationOp op) const {
    for (const auto& rel : relations_) {
        if (rel.op != op) {
            continue;
        }
        // Compare LHS structurally
        if (!lamina::detail::node(rel.lhs) && !lamina::detail::node(lhs)) {
            // Both null — match on LHS
        } else if (!lamina::detail::node(rel.lhs) || !lamina::detail::node(lhs)) {
            continue;  // One null, one not — no match
        } else if (!lamina::detail::node(rel.lhs)->equals(*lamina::detail::node(lhs))) {
            continue;
        }
        // Compare RHS structurally
        if (!lamina::detail::node(rel.rhs) && !lamina::detail::node(rhs)) {
            return true;  // Both null — full match
        } else if (!lamina::detail::node(rel.rhs) || !lamina::detail::node(rhs)) {
            continue;  // One null, one not — no match
        } else if (lamina::detail::node(rel.rhs)->equals(*lamina::detail::node(rhs))) {
            return true;
        }
    }
    return false;
}

void RelationStore::clear() {
    relations_.clear();
}

} // namespace lamina
