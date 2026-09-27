#include "assumption_context.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include "computation_context.hpp"
#include "internal/rewrite_budget.hpp"
#include <limits>
#include <set>

namespace LMCAS {

namespace {

bool has_value_declaration(const PropertyStore& properties, const std::string& symbol) {
    return !properties.get_signs(symbol).empty() ||
        properties.get_domain(symbol) != Domain::Complex ||
        properties.get_parity(symbol) != Parity::Unknown ||
        properties.get_boundedness(symbol) != Boundedness::Unknown ||
        properties.get_finiteness(symbol) != Finiteness::Unknown;
}

bool has_visible_proof(
    const std::vector<std::array<std::size_t, 2>>& alternatives,
    const std::vector<bool>& visible, detail::RewriteBudget& budget) {
    for (const auto& premises : alternatives) {
        auto step = budget.context().consume_steps(1, "assumption.relations");
        if (!step) throw step.error();
        if (visible[premises[0]] && visible[premises[1]]) return true;
    }
    return false;
}

}


AssumptionContext::AssumptionContext() {
    scope_stack_.emplace_back();
}

AssumptionContext::AssumptionContext(AssumptionContext&& other)
    : AssumptionContext() {
    scope_stack_.swap(other.scope_stack_);
    max_query_depth_ = other.max_query_depth_;
    ++other.cache_generation_;
}

AssumptionContext& AssumptionContext::operator=(const AssumptionContext& other) {
    if (this != &other) {
        auto candidate = other.scope_stack_;
        scope_stack_.swap(candidate);
        max_query_depth_ = other.max_query_depth_;
        ++cache_generation_;
    }
    return *this;
}

AssumptionContext& AssumptionContext::operator=(AssumptionContext&& other) noexcept {
    if (this != &other) {
        scope_stack_.swap(other.scope_stack_);
        std::swap(max_query_depth_, other.max_query_depth_);
        ++cache_generation_;
        ++other.cache_generation_;
    }
    return *this;
}

uint64_t AssumptionContext::cache_generation() const {
    bool changed = false;
    for (const auto& scope : scope_stack_) {
        const auto properties = scope.properties.revision();
        const auto relations = scope.relations.revision();
        if (scope.observed_property_revision != properties ||
            scope.observed_relation_revision != relations) {
            scope.observed_property_revision = properties;
            scope.observed_relation_revision = relations;
            changed = true;
        }
    }
    if (changed) ++cache_generation_;
    return cache_generation_;
}

// Scope management

void AssumptionContext::push() {
    scope_stack_.emplace_back();
    ++cache_generation_;
}

AssumptionVoidResult AssumptionContext::pop() {
    if (scope_stack_.size() <= 1) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "cannot pop the root assumption scope", "assumption.pop");
    }
    scope_stack_.pop_back();
    ++cache_generation_;
    return AssumptionVoidResult::success();
}

int AssumptionContext::depth() const {
    return static_cast<int>(scope_stack_.size());
}

// Direct access to current (top) scope stores

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

// Read-through query methods

bool AssumptionContext::has_sign(const std::string& symbol, Sign sign) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        auto signs = it->properties.get_signs(symbol);
        if (!signs.empty()) return signs.count(sign) > 0;
    }
    return false;
}

bool AssumptionContext::has_domain(const std::string& symbol, Domain domain) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        const auto declared_domain = it->properties.get_domain(symbol);
        if (declared_domain != Domain::Complex) {
            return it->properties.has_domain(symbol, domain);
        }
    }
    return domain == Domain::Complex;
}

Domain AssumptionContext::get_domain(const std::string& symbol) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        const auto domain = it->properties.get_domain(symbol);
        if (domain != Domain::Complex) return domain;
    }
    return Domain::Complex;
}

std::unordered_set<Sign, SignHash> AssumptionContext::get_signs(
    const std::string& symbol) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        auto signs = it->properties.get_signs(symbol);
        if (!signs.empty()) return signs;
    }
    return {};
}

Parity AssumptionContext::get_parity(const std::string& symbol) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        const auto parity = it->properties.get_parity(symbol);
        if (parity != Parity::Unknown) return parity;
    }
    return Parity::Unknown;
}

Boundedness AssumptionContext::get_boundedness(const std::string& symbol) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        const auto boundedness = it->properties.get_boundedness(symbol);
        if (boundedness != Boundedness::Unknown) return boundedness;
    }
    return Boundedness::Unknown;
}

std::optional<Interval> AssumptionContext::get_bounds(const std::string& symbol) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        auto bounds = it->properties.get_bounds(symbol);
        if (bounds) return bounds;
    }
    return std::nullopt;
}

bool AssumptionContext::is_transcendental(const std::string& symbol) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (it->properties.is_transcendental(symbol)) return true;
        if (it->properties.has_domain(symbol, Domain::Algebraic)) return false;
    }
    return false;
}

Finiteness AssumptionContext::get_finiteness(const std::string& symbol) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        auto value = it->properties.get_finiteness(symbol);
        if (value != Finiteness::Unknown) return value;
    }
    return Finiteness::Unknown;
}

Definiteness AssumptionContext::get_definiteness(const std::string& symbol) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        auto value = it->properties.get_definiteness(symbol);
        if (value != Definiteness::Unknown) return value;
    }
    return Definiteness::Unknown;
}

std::optional<SymbolicExpr> AssumptionContext::get_period(
    const std::string& symbol, const std::string& variable) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        auto value = it->properties.get_period(symbol, variable);
        if (value) return value;
    }
    return std::nullopt;
}

bool AssumptionContext::is_periodic(const std::string& symbol, const std::string& variable) const {
    return get_period(symbol, variable).has_value();
}

bool AssumptionContext::has_period_declarations(const std::string& symbol) const {
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        if (!it->properties.get_period_decls(symbol).empty()) return true;
    }
    return false;
}

Result<Monotonicity> AssumptionContext::get_monotonicity_checked(
    const std::string& symbol, const std::string& variable,
    const Interval& interval) const {
    ComputationContext context;
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        auto value = it->properties.get_monotonicity_checked(
            symbol, variable, interval, context);
        if (!value || value.value() != Monotonicity::Unknown) return value;
    }
    return Result<Monotonicity>::success(Monotonicity::Unknown);
}

struct AssumptionContext::RelationShadowing {
    std::set<std::string> symbols;
    std::vector<SymbolicExpr> expressions;
};

AssumptionContext::RelationShadowing AssumptionContext::collect_relation_shadowing(
    std::size_t scope, detail::RewriteBudget& budget) const {
    RelationShadowing shadowing;
    for (std::size_t i = scope + 1; i < scope_stack_.size(); ++i) {
        auto step = budget.context().consume_steps(1, "assumption.relations");
        if (!step) throw step.error();
        const auto& nearer = scope_stack_[i];
        const auto& store = nearer.relations;
        for (std::size_t j = 0; j < store.relations_.size(); ++j) {
            step = budget.context().consume_steps(1, "assumption.relations");
            if (!step) throw step.error();
            if (!store.proofs_[j].declared) continue;
            const auto& declared = store.relations_[j];
            const SymbolicExpr* subject = nullptr;
            if (!free_variables(detail::node(declared.lhs), &budget).empty()) {
                subject = &declared.lhs;
            } else if (!free_variables(detail::node(declared.rhs), &budget).empty()) {
                subject = &declared.rhs;
            }
            if (!subject) continue;
            if (const auto* variable =
                    dynamic_cast<const VariableNode*>(detail::node(*subject).get())) {
                shadowing.symbols.insert(variable->name());
            } else {
                shadowing.expressions.push_back(*subject);
            }
        }
        for (const auto& symbol : nearer.properties.get_all_symbols()) {
            step = budget.context().consume_steps(1, "assumption.relations");
            if (!step) throw step.error();
            if (has_value_declaration(nearer.properties, symbol)) {
                shadowing.symbols.insert(symbol);
            }
        }
    }
    return shadowing;
}

bool AssumptionContext::relation_shadowed(
    const Relation& relation, const RelationShadowing& shadowing,
    detail::RewriteBudget& budget) {
    const auto& lhs = detail::node(relation.lhs);
    const auto& rhs = detail::node(relation.rhs);
    for (const auto& symbol : shadowing.symbols) {
        auto step = budget.context().consume_steps(1, "assumption.relations");
        if (!step) { throw step.error(); }
        if (expression_depends_on_variable(lhs, symbol, &budget) ||
            expression_depends_on_variable(rhs, symbol, &budget)) { return true; }
    }
    if (shadowing.expressions.empty()) { return false; }
    budget.measure(lhs);
    budget.measure(rhs);
    for (const auto& expression : shadowing.expressions) {
        auto step = budget.context().consume_steps(1, "assumption.relations");
        if (!step) { throw step.error(); }
        const auto& subject = detail::node(expression);
        budget.measure(subject);
        if ((lhs && lhs->equals(*subject)) || (rhs && rhs->equals(*subject))) {
            return true;
        }
    }
    return false;
}

std::vector<bool> AssumptionContext::relation_visibility(
    std::size_t scope, detail::RewriteBudget& budget) const {
    // An empty mask denotes that every relation in the top scope is visible.
    if (scope + 1 == scope_stack_.size()) return {};
    const auto& store = scope_stack_[scope].relations;
    const auto shadowing = collect_relation_shadowing(scope, budget);
    std::vector<bool> unshadowed(store.relations_.size(), false);
    std::vector<bool> visible(store.relations_.size(), false);
    for (std::size_t i = 0; i < store.relations_.size(); ++i) {
        auto step = budget.context().consume_steps(1, "assumption.relations");
        if (!step) throw step.error();
        unshadowed[i] = !relation_shadowed(store.relations_[i], shadowing, budget);
        visible[i] = unshadowed[i] && store.proofs_[i].declared;
    }
    bool changed;
    do {
        auto step = budget.context().consume_steps(1, "assumption.relations");
        if (!step) throw step.error();
        changed = false;
        for (std::size_t i = 0; i < store.relations_.size(); ++i) {
            step = budget.context().consume_steps(1, "assumption.relations");
            if (!step) throw step.error();
            if (!unshadowed[i] || visible[i]) continue;
            if (!has_visible_proof(store.proofs_[i].alternatives, visible, budget)) continue;
            visible[i] = true;
            changed = true;
        }
    } while (changed);
    return visible;
}

Result<bool> AssumptionContext::has_relation_in_scope(
    std::size_t scope, const SymbolicExpr& lhs, const SymbolicExpr& rhs,
    RelationOp op, detail::RewriteBudget& budget) const {
    const auto& relations = scope_stack_[scope].relations.get_relations();
    const auto visible = relation_visibility(scope, budget);
    const auto& left = detail::node(lhs);
    const auto& right = detail::node(rhs);
    for (std::size_t j = 0; j < relations.size(); ++j) {
        auto step = budget.context().consume_steps(1, "assumption.has_relation");
        if (!step) { return Result<bool>::failure(step.error()); }
        const auto& relation = relations[j];
        if (relation.op != op || (!visible.empty() && !visible[j])) { continue; }
        const auto& relation_lhs = detail::node(relation.lhs);
        const auto& relation_rhs = detail::node(relation.rhs);
        if (!relation_lhs || !relation_rhs) { continue; }
        budget.measure(relation_lhs);
        budget.measure(relation_rhs);
        if (left->equals(*relation_lhs) && right->equals(*relation_rhs)) { return true; }
    }
    return false;
}

bool AssumptionContext::has_relation(
    const SymbolicExpr& lhs, const SymbolicExpr& rhs, RelationOp op) const {
    ComputationContext context;
    auto result = has_relation_checked(lhs, rhs, op, context);
    if (!result) throw result.error();
    return result.value();
}

Result<bool> AssumptionContext::has_relation_checked(
    const SymbolicExpr& lhs, const SymbolicExpr& rhs, RelationOp op,
    ComputationContext& context) const {
    try {
        detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
            std::numeric_limits<std::size_t>::max(), "assumption.has_relation");
        auto step = context.consume_steps(1, "assumption.has_relation");
        if (!step) { return Result<bool>::failure(step.error()); }
        const auto& left = detail::node(lhs);
        const auto& right = detail::node(rhs);
        if (!left || !right) { return false; }
        budget.measure(left);
        budget.measure(right);
        for (std::size_t i = scope_stack_.size(); i-- > 0;) {
            step = context.consume_steps(1, "assumption.has_relation");
            if (!step) { return Result<bool>::failure(step.error()); }
            auto match = has_relation_in_scope(i, lhs, rhs, op, budget);
            if (!match || match.value()) { return match; }
        }
        return false;
    } catch (const CasError& error) {
        return Result<bool>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<bool>::failure(CasErrc::ResourceLimit,
            "relation query allocation failed", "assumption.has_relation");
    }
}

std::vector<Relation> AssumptionContext::get_visible_relations() const {
    ComputationContext context;
    auto result = get_visible_relations_checked(context);
    if (!result) throw result.error();
    return std::move(result.value());
}

Result<std::vector<Relation>> AssumptionContext::get_visible_relations_checked(
    ComputationContext& context) const {
    try {
        detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
            std::numeric_limits<std::size_t>::max(), "assumption.visible_relations");
        std::vector<Relation> visible;
        for (std::size_t i = scope_stack_.size(); i-- > 0;) {
            auto step = context.consume_steps(1, "assumption.visible_relations");
            if (!step) return Result<std::vector<Relation>>::failure(step.error());
            const auto& relations = scope_stack_[i].relations.get_relations();
            const auto visibility = relation_visibility(i, budget);
            for (std::size_t j = 0; j < relations.size(); ++j) {
                step = context.consume_steps(1, "assumption.visible_relations");
                if (!step) return Result<std::vector<Relation>>::failure(step.error());
                if (visibility.empty() || visibility[j]) visible.push_back(relations[j]);
            }
        }
        return visible;
    } catch (const CasError& error) {
        return Result<std::vector<Relation>>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<std::vector<Relation>>::failure(CasErrc::ResourceLimit,
            "relation query allocation failed", "assumption.visible_relations");
    }
}

AssumptionVoidResult AssumptionContext::assume_domain(
    const std::string& variable,
    Domain domain) {
    return assume_domain_checked(variable, domain);
}

AssumptionVoidResult AssumptionContext::assume_domain_checked(
    const std::string& variable,
    Domain domain) {
    constexpr const char* operation = "assume_domain";
    if (variable.empty()) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "variable name must not be empty", operation);
    }
    auto result = scope_stack_.back().properties.declare_domain_checked(variable, domain);
    if (!result) {
        return AssumptionVoidResult::failure(
            result.error().code, result.error().message, operation);
    }
    return AssumptionVoidResult::success();
}

AssumptionVoidResult AssumptionContext::assume_sign(
    const std::string& variable,
    Sign sign) {
    return assume_sign_checked(variable, sign);
}

AssumptionVoidResult AssumptionContext::assume_sign_checked(
    const std::string& variable,
    Sign sign) {
    constexpr const char* operation = "assume_sign";
    if (variable.empty()) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "variable name must not be empty", operation);
    }
    auto result = scope_stack_.back().properties.declare_sign_checked(variable, sign);
    if (!result) {
        return AssumptionVoidResult::failure(
            result.error().code, result.error().message, operation);
    }
    return AssumptionVoidResult::success();
}

AssumptionVoidResult AssumptionContext::assume(const SymbolicExpr& relation) {
    return assume_checked(relation);
}

AssumptionVoidResult AssumptionContext::assume_checked(const SymbolicExpr& relation) {
    constexpr const char* operation = "assume";
    if (!LMCAS::detail::node(relation)) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "relation expression must not be null", operation);
    }
    if (!std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(relation))) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "relation expression root must be relational", operation);
    }
    try {
        auto rel_node = std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(relation));
        auto lhs = LMCAS::detail::expression_from_node(rel_node->left());
        auto rhs = LMCAS::detail::expression_from_node(rel_node->right());
        auto result = scope_stack_.back().relations.add_relation_checked(
            lhs, rhs, rel_node->op(), scope_stack_.back().properties);
        if (!result.has_value()) {
            return AssumptionVoidResult::failure(result.error());
        }
    } catch (const std::bad_alloc&) {
        return AssumptionVoidResult::failure(
            CasErrc::ResourceLimit, "assumption allocation failed", operation);
    } catch (const std::invalid_argument& ex) {
        return AssumptionVoidResult::failure(CasErrc::InvalidArgument, ex.what(), operation);
    } catch (const std::exception& ex) {
        return AssumptionVoidResult::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
    return AssumptionVoidResult::success();
}
// Depth limit configuration

void AssumptionContext::set_max_query_depth(int depth) {
    if (depth <= 0 || depth == max_query_depth_) return;
    max_query_depth_ = depth;
    ++cache_generation_;
}

int AssumptionContext::get_max_query_depth() const {
    return max_query_depth_;
}

}
