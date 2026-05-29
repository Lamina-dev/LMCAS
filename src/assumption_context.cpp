/**
 * @file assumption_context.cpp
 * @brief Implementation of the AssumptionContext class.
 *
 * Provides scoped push/pop management with read-through query semantics.
 * Each scope has its own PropertyStore and RelationStore. Queries search
 * from the top scope down to root, with child declarations shadowing parent.
 */

#include "assumption_context.hpp"
#include "query_interface.hpp"
#include "interval.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"

namespace lamina {

// ============================================================
// Construction
// ============================================================

AssumptionContext::AssumptionContext() {
    // Start with root scope
    scope_stack_.emplace_back();
}

// ============================================================
// Scope management
// ============================================================

void AssumptionContext::push() {
    scope_stack_.emplace_back();
}

void AssumptionContext::pop() {
    if (scope_stack_.size() <= 1) {
        throw std::runtime_error("Cannot pop root scope");
    }
    scope_stack_.pop_back();
}

int AssumptionContext::depth() const {
    return static_cast<int>(scope_stack_.size());
}

// ============================================================
// Direct access to current (top) scope stores
// ============================================================

PropertyStore& AssumptionContext::current_properties() {
    return scope_stack_.back().properties;
}

const PropertyStore& AssumptionContext::current_properties() const {
    return scope_stack_.back().properties;
}

RelationStore& AssumptionContext::current_relations() {
    return scope_stack_.back().relations;
}

const RelationStore& AssumptionContext::current_relations() const {
    return scope_stack_.back().relations;
}

// ============================================================
// Read-through query methods
// ============================================================

bool AssumptionContext::has_sign(const std::string& symbol, Sign sign) const {
    // Search from top scope down to root.
    // Return the result from the first scope that has sign info for this symbol.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (!it->properties.get_signs(symbol).empty()) {
            return it->properties.has_sign(symbol, sign);
        }
    }
    // No scope has sign info for this symbol
    return false;
}

bool AssumptionContext::has_domain(const std::string& symbol, Domain domain) const {
    // Search from top scope down to root.
    // Return the result from the first scope that has domain info for this symbol.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->properties.get_domain(symbol) != Domain::Complex) {
            return it->properties.has_domain(symbol, domain);
        }
    }
    // No scope has domain info — default is Complex.
    // has_domain checks if the symbol has at least the given specificity.
    return domain == Domain::Complex;
}

Domain AssumptionContext::get_domain(const std::string& symbol) const {
    // Search from top scope down to root.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->properties.get_domain(symbol) != Domain::Complex) {
            return it->properties.get_domain(symbol);
        }
    }
    return Domain::Complex;
}

std::unordered_set<Sign, SignHash> AssumptionContext::get_signs(const std::string& symbol) const {
    // Search from top scope down to root.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (!it->properties.get_signs(symbol).empty()) {
            return it->properties.get_signs(symbol);
        }
    }
    return {};
}

Parity AssumptionContext::get_parity(const std::string& symbol) const {
    // Search from top scope down to root.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->properties.get_parity(symbol) != Parity::Unknown) {
            return it->properties.get_parity(symbol);
        }
    }
    return Parity::Unknown;
}

Boundedness AssumptionContext::get_boundedness(const std::string& symbol) const {
    // Search from top scope down to root.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->properties.get_boundedness(symbol) != Boundedness::Unknown) {
            return it->properties.get_boundedness(symbol);
        }
    }
    return Boundedness::Unknown;
}

std::optional<Interval> AssumptionContext::get_bounds(const std::string& symbol) const {
    // Search from top scope down to root.
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->properties.get_bounds(symbol).has_value()) {
            return it->properties.get_bounds(symbol);
        }
    }
    return std::nullopt;
}

// ============================================================
// Convenience declaration API
// ============================================================

void AssumptionContext::assume_domain(const std::string& variable, Domain domain) {
    if (variable.empty()) {
        throw std::invalid_argument("assume_domain: variable name must not be empty");
    }
    scope_stack_.back().properties.declare_domain(variable, domain);
}

void AssumptionContext::assume_sign(const std::string& variable, Sign sign) {
    if (variable.empty()) {
        throw std::invalid_argument("assume_sign: variable name must not be empty");
    }
    scope_stack_.back().properties.declare_sign(variable, sign);
}

void AssumptionContext::assume(const SymbolicExpr& relation) {
    if (!relation.root) {
        throw std::invalid_argument("assume: expression must not be null/empty");
    }
    auto rel_node = std::dynamic_pointer_cast<RelationalNode>(relation.root);
    if (!rel_node) {
        throw std::invalid_argument("assume: expression root must be a RelationalNode");
    }
    // Extract lhs, rhs, and op from the RelationalNode and store in RelationStore
    SymbolicExpr lhs(rel_node->left);
    SymbolicExpr rhs(rel_node->right);
    scope_stack_.back().relations.add_relation(lhs, rhs, rel_node->op,
                                               scope_stack_.back().properties);
}

// ============================================================
// Convenience query API (delegates to QueryInterface)
// ============================================================

Tribool AssumptionContext::is_positive(const SymbolicExpr& expr) const {
    if (!expr.root) {
        throw std::invalid_argument("is_positive: expression must not be null/empty");
    }
    QueryInterface qi(*this);
    return qi.query_positive(expr);
}

Tribool AssumptionContext::is_negative(const SymbolicExpr& expr) const {
    if (!expr.root) {
        throw std::invalid_argument("is_negative: expression must not be null/empty");
    }
    QueryInterface qi(*this);
    return qi.query_negative(expr);
}

Tribool AssumptionContext::is_nonnegative(const SymbolicExpr& expr) const {
    if (!expr.root) {
        throw std::invalid_argument("is_nonnegative: expression must not be null/empty");
    }
    QueryInterface qi(*this);
    return qi.query_nonnegative(expr);
}

Tribool AssumptionContext::is_real(const SymbolicExpr& expr) const {
    if (!expr.root) {
        throw std::invalid_argument("is_real: expression must not be null/empty");
    }
    QueryInterface qi(*this);
    return qi.query_real(expr);
}

Tribool AssumptionContext::is_integer(const SymbolicExpr& expr) const {
    if (!expr.root) {
        throw std::invalid_argument("is_integer: expression must not be null/empty");
    }
    QueryInterface qi(*this);
    return qi.query_integer(expr);
}

Tribool AssumptionContext::is_nonzero(const SymbolicExpr& expr) const {
    if (!expr.root) {
        throw std::invalid_argument("is_nonzero: expression must not be null/empty");
    }
    QueryInterface qi(*this);
    return qi.query_nonzero(expr);
}

} // namespace lamina
