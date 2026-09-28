#include "relation_store.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include "property_store.hpp"
#include <algorithm>
#include <queue>

namespace LMCAS {

RelationStore::RelationStore(RelationStore&& other) noexcept
    : relations_(std::move(other.relations_)), proofs_(std::move(other.proofs_)),
      declaration_order_(std::move(other.declaration_order_)) {
    other.relations_.clear();
    other.proofs_.clear();
    other.declaration_order_.clear();
    ++other.revision_;
}

RelationStore& RelationStore::operator=(const RelationStore& other) {
    if (this != &other) {
        RelationStore candidate(other);
        relations_.swap(candidate.relations_);
        proofs_.swap(candidate.proofs_);
        declaration_order_.swap(candidate.declaration_order_);
        ++revision_;
    }
    return *this;
}

RelationStore& RelationStore::operator=(RelationStore&& other) noexcept {
    if (this != &other) {
        relations_ = std::move(other.relations_);
        proofs_ = std::move(other.proofs_);
        declaration_order_ = std::move(other.declaration_order_);
        other.relations_.clear();
        other.proofs_.clear();
        other.declaration_order_.clear();
        ++revision_;
        ++other.revision_;
    }
    return *this;
}

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
    if (!LMCAS::detail::node(a) && !LMCAS::detail::node(b)) return true;
    if (!LMCAS::detail::node(a) || !LMCAS::detail::node(b)) return false;
    return LMCAS::detail::node(a)->equals(*LMCAS::detail::node(b));
}

RelationStoreResult declare_comparison_sign(const std::string& symbol,
                                          RelationOp op, PropertyStore& store) {
    struct SignRule {
        RelationOp op;
        Sign sign;
    };
    static constexpr SignRule rules[] = {
        {RelationOp::GT, Sign::Positive},
        {RelationOp::GEQ, Sign::NonNegative},
        {RelationOp::LT, Sign::Negative},
        {RelationOp::LEQ, Sign::NonPositive},
        {RelationOp::NEQ, Sign::NonZero}
    };
    for (const auto& rule : rules) {
        if (rule.op == op) return store.declare_sign(symbol, rule.sign);
    }
    return RelationStoreResult::success();
}


RelationStoreResult derive_comparison_sign(const SymbolicExpr& lhs,
                                         const SymbolicExpr& rhs,
                                         RelationOp op, PropertyStore& store) {
    const auto variable = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(lhs));
    const auto number = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(rhs));
    if (variable && !variable->is_constant() && number && number->is_zero()) {
        return declare_comparison_sign(variable->name(), op, store);
    }
    const auto zero = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(lhs));
    const auto reversed_variable = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(rhs));
    if (zero && zero->is_zero() && reversed_variable &&
        !reversed_variable->is_constant()) {
        return declare_comparison_sign(
            reversed_variable->name(), detail::reversed_relation(op), store);
    }
    return RelationStoreResult::success();
}

}

RelationStoreResult RelationStore::add_relation_unchecked(
    const SymbolicExpr& lhs,
    const SymbolicExpr& rhs,
    RelationOp op,
    PropertyStore& prop_store) {
    auto index = find_relation(lhs, rhs, op);
    if (!index) {
        index = relations_.size();
        relations_.push_back(Relation{lhs, rhs, op});
        proofs_.push_back(RelationProofs{});
    }
    const bool already_declared = proofs_[*index].declared;
    if (!already_declared) {
        proofs_[*index].declared = true;
        declaration_order_.push_back(*index);
    }

    if (!LMCAS::detail::node(lhs) || !LMCAS::detail::node(rhs)) {
        return RelationStoreResult::success();
    }

    auto derived = derive_comparison_sign(lhs, rhs, op, prop_store);
    if (!derived) return derived;

    // Compute transitive closure for GT/GEQ relations
    if (is_transitive_op(op)) {
        auto closure = compute_transitive_closure(*index, prop_store);
        if (!closure) return RelationStoreResult::failure(closure.error());
        if (already_declared && closure.value()) {
            // 重复声明可以推进受限闭包，序列化需重放这次操作。
            declaration_order_.push_back(*index);
        }
    }
    return RelationStoreResult::success();
}

RelationStoreResult RelationStore::add_relation(
    const SymbolicExpr& lhs,
    const SymbolicExpr& rhs,
    RelationOp op,
    PropertyStore& prop_store) {
    return add_relation_checked(lhs, rhs, op, prop_store);
}

RelationStoreResult RelationStore::add_relation_checked(
    const SymbolicExpr& lhs,
    const SymbolicExpr& rhs,
    RelationOp op,
    PropertyStore& prop_store) {
    if (!LMCAS::detail::node(lhs)) {
        return RelationStoreResult::failure(
            CasErrc::InvalidArgument, "relation lhs must not be null", "add_relation");
    }
    if (!LMCAS::detail::node(rhs)) {
        return RelationStoreResult::failure(
            CasErrc::InvalidArgument, "relation rhs must not be null", "add_relation");
    }

    try {
        RelationStore relation_candidate = *this;
        PropertyStore property_candidate = prop_store;
        auto inserted = relation_candidate.add_relation_unchecked(
            lhs, rhs, op, property_candidate);
        if (!inserted) return inserted;
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
struct RelationStore::TransitiveWork {
    std::queue<std::size_t> queue;
    int deductions = 0;
    bool changed = false;
};

RelationStoreResult RelationStore::enqueue_transitive_deduction(
    const Relation& deduced, Premises premises,
    PropertyStore& prop_store, TransitiveWork& work) {
    if (premises[1] < premises[0]) std::swap(premises[0], premises[1]);
    if (const auto index = find_relation(deduced.lhs, deduced.rhs, deduced.op)) {
        if (premises[0] == *index || premises[1] == *index) {
            return RelationStoreResult::success();
        }
        auto& alternatives = proofs_[*index].alternatives;
        if (std::find(alternatives.begin(), alternatives.end(), premises) ==
            alternatives.end()) {
            alternatives.push_back(premises);
            work.changed = true;
        }
        return RelationStoreResult::success();
    }
    if (work.deductions == MAX_TRANSITIVE_DEDUCTIONS) {
        return RelationStoreResult::success();
    }
    relations_.push_back(deduced);
    proofs_.push_back(RelationProofs{false, {premises}});
    ++work.deductions;
    work.changed = true;
    auto derived = derive_comparison_sign(
        deduced.lhs, deduced.rhs, deduced.op, prop_store);
    if (!derived) return derived;
    work.queue.push(relations_.size() - 1);
    return RelationStoreResult::success();
}

Result<bool> RelationStore::compute_transitive_closure(
    std::size_t new_relation, PropertyStore& prop_store) {
    TransitiveWork work;
    work.queue.push(new_relation);
    while (!work.queue.empty()) {
        const auto current_index = work.queue.front();
        const auto current = relations_[current_index];
        work.queue.pop();
        /**
         * @brief 双向推导共用同一快照，先正向后反向。
         * @note 插入可能使向量引用失效，跨插入仅保留索引或副本。
         */
        const std::size_t relation_count = relations_.size();
        for (std::size_t i = 0; i < relation_count; ++i) {
            if (!is_transitive_op(relations_[i].op)) continue;
            if (!expr_equals(current.rhs, relations_[i].lhs)) continue;
            auto result = enqueue_transitive_deduction(
                {current.lhs, relations_[i].rhs, combine_ops(current.op, relations_[i].op)},
                {current_index, i}, prop_store, work);
            if (!result) return Result<bool>::failure(result.error());
        }
        for (std::size_t i = 0; i < relation_count; ++i) {
            if (!is_transitive_op(relations_[i].op)) continue;
            if (!expr_equals(relations_[i].rhs, current.lhs)) continue;
            auto result = enqueue_transitive_deduction(
                {relations_[i].lhs, current.rhs, combine_ops(relations_[i].op, current.op)},
                {i, current_index}, prop_store, work);
            if (!result) return Result<bool>::failure(result.error());
        }
    }
    return Result<bool>::success(work.changed);
}

const std::vector<Relation>& RelationStore::get_relations() const {
    return relations_;
}

std::optional<std::size_t> RelationStore::find_relation(
    const SymbolicExpr& lhs, const SymbolicExpr& rhs, RelationOp op) const {
    for (std::size_t i = 0; i < relations_.size(); ++i) {
        const auto& relation = relations_[i];
        if (relation.op == op && expr_equals(relation.lhs, lhs) &&
            expr_equals(relation.rhs, rhs)) {
            return i;
        }
    }
    return std::nullopt;
}

bool RelationStore::has_relation(const SymbolicExpr& lhs, const SymbolicExpr& rhs,
                                 RelationOp op) const {
    return find_relation(lhs, rhs, op).has_value();
}

void RelationStore::clear() {
    if (!relations_.empty()) {
        relations_.clear();
        proofs_.clear();
        declaration_order_.clear();
        ++revision_;
    }
}

}
