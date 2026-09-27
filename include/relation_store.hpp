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
#include <cstdint>
#include <array>
#include <cstddef>
#include <optional>

namespace LMCAS {

// Forward declaration - PropertyStore may not yet be compiled
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
 *   - GT  -> Positive
 *   - GEQ -> NonNegative
 *   - LT  -> Negative
 *   - LEQ -> NonPositive
 *   - NEQ -> NonZero
 *
 * 复合关系(多变量左端或任意右端)按原结构存储,供 InferenceEngine 后续推导.
 */
class LMCAS_API RelationStore {
public:
    RelationStore() = default;
    RelationStore(const RelationStore&) = default;
    RelationStore(RelationStore&& other) noexcept;
    RelationStore& operator=(const RelationStore& other);
    RelationStore& operator=(RelationStore&& other) noexcept;

    /** @brief 成功提交（含赋值）均更新此存储的修订号。 */
    std::uint64_t revision() const noexcept { return revision_; }

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
    friend class AssumptionContext;
    using Premises = std::array<std::size_t, 2>;
    struct RelationProofs {
        bool declared = false;
        std::vector<Premises> alternatives;
    };
    std::vector<Relation> relations_;
    std::vector<RelationProofs> proofs_;
    std::vector<std::size_t> declaration_order_;
    std::uint64_t revision_ = 0;

    std::optional<std::size_t> find_relation(
        const SymbolicExpr& lhs, const SymbolicExpr& rhs, RelationOp op) const;

    RelationStoreResult add_relation_unchecked(
        const SymbolicExpr& lhs,
        const SymbolicExpr& rhs,
        RelationOp op,
        PropertyStore& prop_store);

    /// Maximum number of new relations deduced per add_relation call via transitive closure.
    static constexpr int MAX_TRANSITIVE_DEDUCTIONS = 64;

    struct TransitiveWork;
    RelationStoreResult enqueue_transitive_deduction(
        const Relation& deduced, Premises premises,
        PropertyStore& prop_store, TransitiveWork& work);

    /**
     * @brief 新增关系后计算传递闭包。
     *
     * 从新关系开始 BFS，仅对 GT、GEQ 进行传递组合：
     * GT+GT→GT、GT+GEQ→GT、GEQ+GT→GT、GEQ+GEQ→GEQ。
     * 最多新增 MAX_TRANSITIVE_DEDUCTIONS 个结论，之后仍登记已有结论的替代证明。
     * @param new_relation 触发本轮闭包计算的关系索引
     * @param prop_store 用于推导关系符号的属性存储
     * @return 本轮是否新增结论或替代证明。
     */
    Result<bool> compute_transitive_closure(
        std::size_t new_relation, PropertyStore& prop_store);
};

} // namespace LMCAS
