/**
 * @file query_interface.hpp
 * @brief QueryInterface class - unified tri-state query API for property questions on arbitrary expressions.
 *
 * The QueryInterface is the single public entry point for property queries on SymbolicExpr trees.
 * It handles special cases (NaN, Infinity, Matrix, Relational, Logical nodes) and
 * delegates to the InferenceEngine for NumberNode, VariableNode, and composite node queries.
 */
#pragma once

#include "assumption.hpp"
#include "symbolic.hpp"
#include "relation_store.hpp"
#include "result.hpp"
#include <unordered_map>
#include <functional>
#include <vector>
#include <utility>
#include <optional>

namespace lamina {

// Forward declaration
class AssumptionContext;
class InferenceEngine;

using QueryTriboolResult = Result<Tribool>;
using QueryPeriodResult = Result<std::optional<SymbolicExpr>>;

/**
 * @brief Identifies which property is being queried, used as part of the cache key.
 */
enum class PropType {
    Positive,
    Negative,
    NonNegative,
    NonPositive,
    Real,
    Integer,
    NonZero,
    Algebraic,
    Transcendental,
    Finite,
    Divergent,
    Periodic,
    PositiveDefinite,
    PositiveSemiDefinite
};

/**
 * @brief Cache key combining an expression's structural hash with the property type.
 *
 * Uses the expression structural hash - two structurally identical
 * expressions will produce the same hash regardless of pointer identity.
 */
struct CacheKey {
    std::size_t expression_hash;
    PropType property;

    bool operator==(const CacheKey& other) const {
        return expression_hash == other.expression_hash && property == other.property;
    }
};

/**
 * @brief Hash function for CacheKey, combining expression hash and property type.
 */
struct CacheKeyHash {
    std::size_t operator()(const CacheKey& key) const {
        std::size_t seed = key.expression_hash;
        seed ^= static_cast<std::size_t>(key.property) + 0x9e3779b9U
            + (seed << 6U) + (seed >> 2U);
        return seed;
    }
};

/**
 * @brief Unified query API for property questions on arbitrary SymbolicExpr trees.
 *
 * Provides query_positive, query_negative, query_nonnegative, query_real,
 * query_integer, and query_nonzero methods that return Tribool results.
 *
 * 结果先按(表达式哈希,属性类型)分桶,再以结构相等性确认表达式身份.
 * AssumptionContext 发生 push,pop 或 assume 变更时,缓存世代触发整体失效.
 *
 * Dispatch logic:
 *   - Null root node -> Unknown
 *   - MatrixNode, RelationalNode, LogicalNode -> Unknown
 *   - NumberNode with NaN -> False for integer, Unknown for sign queries
 *   - FunctionNode::Infinity -> derive sign from infinity's sign, False for integer
 *   - All other nodes -> delegate to InferenceEngine
 */
class LAMINA_API QueryInterface {
public:
    /**
     * @brief Construct a QueryInterface bound to an AssumptionContext.
     * @param ctx The AssumptionContext used for property lookups and inference
     */
    explicit QueryInterface(const AssumptionContext& ctx);

    /// @name Public query methods
    /// @{

    /// Query whether the expression is positive (> 0).
    QueryTriboolResult query_positive(const SymbolicExpr& expr) const;

    QueryTriboolResult query_positive_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is negative (< 0).
    QueryTriboolResult query_negative(const SymbolicExpr& expr) const;

    QueryTriboolResult query_negative_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is non-negative (>= 0).
    QueryTriboolResult query_nonnegative(const SymbolicExpr& expr) const;

    QueryTriboolResult query_nonnegative_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is real.
    QueryTriboolResult query_real(const SymbolicExpr& expr) const;

    QueryTriboolResult query_real_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is an integer.
    QueryTriboolResult query_integer(const SymbolicExpr& expr) const;

    QueryTriboolResult query_integer_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is non-zero (!= 0).
    QueryTriboolResult query_nonzero(const SymbolicExpr& expr) const;

    QueryTriboolResult query_nonzero_checked(const SymbolicExpr& expr) const;

    /// @}

    /// @name Extended property queries
    /// @{

    /// Query whether the expression is algebraic (root of a polynomial with rational coefficients).
    QueryTriboolResult query_algebraic(const SymbolicExpr& expr) const;

    QueryTriboolResult query_algebraic_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is transcendental (real but not algebraic).
    QueryTriboolResult query_transcendental(const SymbolicExpr& expr) const;

    QueryTriboolResult query_transcendental_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression has a finite value/limit.
    QueryTriboolResult query_finite(const SymbolicExpr& expr) const;

    QueryTriboolResult query_finite_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression diverges.
    QueryTriboolResult query_divergent(const SymbolicExpr& expr) const;

    QueryTriboolResult query_divergent_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression is periodic.
    QueryTriboolResult query_periodic(const SymbolicExpr& expr) const;

    QueryTriboolResult query_periodic_checked(const SymbolicExpr& expr) const;

    /**
     * @brief Get the period of an expression, if known.
     * @param expr The expression to query
     * @return The period as a SymbolicExpr, or std::nullopt if not periodic or unknown
     */
    QueryPeriodResult get_period(const SymbolicExpr& expr) const;

    QueryPeriodResult get_period_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression (matrix symbol) is positive definite.
    QueryTriboolResult query_positive_definite(const SymbolicExpr& expr) const;

    QueryTriboolResult query_positive_definite_checked(const SymbolicExpr& expr) const;

    /// Query whether the expression (matrix symbol) is positive semidefinite.
    QueryTriboolResult query_positive_semidefinite(const SymbolicExpr& expr) const;

    QueryTriboolResult query_positive_semidefinite_checked(const SymbolicExpr& expr) const;

    /// @}

    /**
     * @brief Invalidate the entire query result cache.
     *
     * AssumptionContext 执行 push_scope,pop_scope,assume_domain,
     * assume_sign 或关系插入后调用本函数.
     */
    void invalidate_cache() const;

    /**
     * @brief A set of sufficient conditions for a property to hold on an expression.
     *
     * Each ConditionSet represents one alternative set of conditions. If all conditions
     * within a single ConditionSet are satisfied, the target property holds.
     */
    struct ConditionSet {
        /// Sign conditions: pairs of (variable_name, required_sign)
        std::vector<std::pair<std::string, Sign>> sign_conditions;
        /// Domain conditions: pairs of (variable_name, required_domain)
        std::vector<std::pair<std::string, Domain>> domain_conditions;
        /// Relational conditions between expressions
        std::vector<Relation> relational_conditions;
    };

    using QueryConditionSetsResult = Result<std::vector<ConditionSet>>;

    /**
     * @brief Query sufficient conditions for a target sign property to hold on an expression.
     *
     * Analyzes the expression structure and returns a list of alternative condition sets.
     * Each condition set, if fully satisfied, guarantees the target sign property holds.
     *
     * - Single variable: returns the direct sign assumption needed (e.g., {x: Positive})
     * - Composite expression (e.g., x - y): derives condition sets from sub-expression
     *   analysis (e.g., {x: Positive, y: Negative} or {x GT y, y: NonNegative})
     * - Returns empty list when no sufficient conditions can be determined
     *
     * @param expr The expression to analyze
     * @param target The target sign property (e.g., Sign::Positive)
     * @return A list of alternative sufficient condition sets (empty if undetermined)
     */
    QueryConditionSetsResult query_conditions(
        const SymbolicExpr& expr, Sign target) const;

    QueryConditionSetsResult query_conditions_checked(const SymbolicExpr& expr, Sign target) const;

private:
    struct CacheEntry {
        SymbolicExpr expression;
        Tribool result;
    };

    const AssumptionContext& ctx_;

    /// The cache generation at the time of last cache validation.
    /// If ctx_.cache_generation() differs, the cache is stale and must be cleared.
    mutable uint64_t observed_generation_;

    /// Query result cache: (expression_hash, property_type) -> collision bucket.
    /// Mutable because queries are logically const but populate the cache.
    mutable std::unordered_map<CacheKey, std::vector<CacheEntry>, CacheKeyHash> cache_;

    /// Check if the root node is a special case that should return Unknown immediately.
    /// Returns true if the node is null, MatrixNode, RelationalNode, or LogicalNode.
    bool is_unhandled_type(const SymbolicExpr& expression) const;

    /// Check if a NumberNode holds NaN.
    bool is_nan_number(const SymbolicExpr& expression) const;

    /// Check if a FunctionNode represents infinity (FuncType::Infinity).
    bool is_infinity_node(const SymbolicExpr& expression) const;

    /// Determine the sign of an infinity expression.
    /// Positive infinity: FunctionNode::Infinity directly.
    /// Negative infinity: MultiplyNode(-1, FunctionNode::Infinity).
    /// Returns +1 for positive infinity, -1 for negative infinity, 0 if indeterminate.
    int get_infinity_sign(const SymbolicExpr& expression) const;

    /**
     * @brief Look up a cached result or compute, cache, and return it.
     * @param expr The expression to query
     * @param prop The property type being queried
     * @param compute Function that performs the actual computation
     * @return Cached or freshly computed Tribool result
     */
    Tribool cached_query(const SymbolicExpr& expr, PropType prop,
                         const std::function<Tribool()>& compute) const;

    QueryTriboolResult cached_query_checked(
        const SymbolicExpr& expr,
        PropType prop,
        const std::string& operation,
        const std::function<Tribool()>& compute) const;

    std::vector<ConditionSet> query_conditions_impl(
        const SymbolicExpr& expr,
        Sign target) const;
};

} // namespace lamina
