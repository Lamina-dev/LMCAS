/**
 * @file inference_engine.hpp
 * @brief InferenceEngine class for deriving properties of composite expressions.
 *
 * Analyzes arithmetic and function expressions to infer sign, domain, and
 * boundedness properties from their sub-expressions.
 */
#pragma once

#include "assumption.hpp"
#include "symbolic.hpp"
#include "relation_store.hpp"
#include "result.hpp"
#include <memory>
#include <optional>

namespace LMCAS {

class AddNode;
class MultiplyNode;
class PowerNode;
class FunctionNode;

// Forward declaration for Interval (avoid circular include with interval.hpp)


struct Interval;
}

namespace LMCAS {

// Forward declaration - AssumptionContext is implemented in a later task
class AssumptionContext;
class ComputationContext;

using InferenceTriboolResult = Result<Tribool>;
using InferencePeriodResult = Result<std::optional<SymbolicExpr>>;

/**
 * @brief 递归分析子表达式，推导复合表达式的符号、定义域和有界性。
 *
 * 符号规则：
 * - 加法：各项同号时和同号；混合符号或未知符号返回 Unknown。
 * - 乘法：所有因子均有定义且为实数时，合并可能符号集合。
 * - 幂：先满足原定义域；零次幂和负整数次幂要求底数非零。
 * - 函数：exp 的参数有定义且为实数时，结果为正。
 *
 * 定义域规则：
 * - 加法：全为 Integer 时得 Integer，全为 Real 时得 Real。
 * - 乘法：全为 Integer 时得 Integer，混合 Real/Integer 时得 Real。
 * - 整数幂：要求底数有定义，零次幂和负次幂还要求底数非零。
 * - 对数：参数有定义且为正时可证明 Real，仅有 Integer 属性不足。
 * - 无定义的实表达式返回 DomainError，前提未证实时保留 Unknown。
 *
 * @par 线程安全
 * 递归保护和活动查询状态可变，同一实例限单线程调用。
 * 绑定不可变上下文的独立实例可并发使用；空闲实例可在下次调用前迁移线程。
 */
class LMCAS_API InferenceEngine {
public:
    /**
     * @brief Construct an InferenceEngine bound to an AssumptionContext.
     * @param ctx The AssumptionContext used to query sub-expression properties
     */
    explicit InferenceEngine(const AssumptionContext& ctx);
    ~InferenceEngine();

    InferenceEngine(const InferenceEngine&) = delete;
    InferenceEngine& operator=(const InferenceEngine&) = delete;
    InferenceEngine(InferenceEngine&&) noexcept;
    InferenceEngine& operator=(InferenceEngine&&) noexcept;

    /// @name Public query methods
    /// These query sign/domain properties for arbitrary expressions by dispatching
    /// to the appropriate inference method based on the expression's root node type.
    /// Context overloads share cancellation and resource limits with recursive queries;
    /// proof-depth exhaustion still returns Unknown, independently of resource errors.
    /// @{


    InferenceTriboolResult query_positive_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_positive_checked(const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult query_negative_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_negative_checked(const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult query_nonnegative_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_nonnegative_checked(const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult query_nonpositive_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_nonpositive_checked(const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult query_real_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_real_checked(const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult query_integer_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_integer_checked(const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult query_nonzero_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_nonzero_checked(const SymbolicExpr& expr, ComputationContext& context) const;

    /**
     * @brief Query whether an expression is algebraic.
     *
     * For variables, checks if the domain is Algebraic or more specific
     * (Rational, Integer, Natural, PositiveInt). For composite expressions,
     * returns Unknown (inference rules for composite algebraic expressions
     * can be added later).
     *
     * @param expr The expression to query
     * @return True if algebraic, False if transcendental, Unknown otherwise
     */

    InferenceTriboolResult query_algebraic_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_algebraic_checked(const SymbolicExpr& expr, ComputationContext& context) const;

    /**
     * @brief Query whether an expression is transcendental.
     *
     * For variables, checks if the transcendental flag is set in the PropertyStore.
     * For composite expressions, returns Unknown.
     *
     * @param expr The expression to query
     * @return True if transcendental, False if algebraic, Unknown otherwise
     */

    InferenceTriboolResult query_transcendental_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_transcendental_checked(const SymbolicExpr& expr, ComputationContext& context) const;

    /**
     * @brief Query whether an expression has a finite value/limit.
     *
     * For variables, checks if Finiteness::Finite is declared in the PropertyStore.
     * For composite expressions, returns Unknown.
     *
     * @param expr The expression to query
     * @return True if finite, False if divergent, Unknown otherwise
     */

    InferenceTriboolResult query_finite_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_finite_checked(const SymbolicExpr& expr, ComputationContext& context) const;

    /**
     * @brief Query whether an expression diverges.
     *
     * For variables, checks if Finiteness::Divergent is declared in the PropertyStore.
     * For composite expressions, returns Unknown.
     *
     * @param expr The expression to query
     * @return True if divergent, False if finite, Unknown otherwise
     */

    InferenceTriboolResult query_divergent_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_divergent_checked(const SymbolicExpr& expr, ComputationContext& context) const;

    /// @}

    /// @name Depth limit configuration
    /// @{

    /**
     * @brief Set the maximum recursion depth for property queries.
     * @param depth Maximum depth (must be > 0). Default is 32.
     */
    void set_max_depth(int depth);

    /**
     * @brief Get the current maximum recursion depth.
     * @return The configured maximum depth
     */
    int get_max_depth() const;

    /// @}

    /**
     * @brief 查询表达式相对于指定自变量的周期性。
     *
     * sin/cos/tan 的实仿射参数要求系数有定义且斜率非零；
     * 非线性参数或前提未证实时返回 Unknown。
     *
     * @param expr 待查询表达式
     * @param variable 非空自变量名
     * @return 周期为 Tribool::True，已知非周期为 False，其余为 Unknown
     */

    InferenceTriboolResult query_periodic_checked(
        const SymbolicExpr& expr, const std::string& variable) const;
    InferenceTriboolResult query_periodic_checked(
        const SymbolicExpr& expr, const std::string& variable, ComputationContext& context) const;

    /**
     * @brief 推导表达式相对于指定自变量的已证实正周期。
     *
     * 仿射斜率按精确符号 pi 缩放周期，周期声明按变量索引。
     * 已证实的常量具有周期性，但无最小正周期。
     *
     * @param expr 待查询表达式
     * @param variable 非空自变量名
     * @return 已证实的周期；常量或周期性未知时返回 nullopt
     */

    InferencePeriodResult infer_period_checked(
        const SymbolicExpr& expr, const std::string& variable) const;
    InferencePeriodResult infer_period_checked(
        const SymbolicExpr& expr, const std::string& variable, ComputationContext& context) const;

    /**
     * @brief Infer the monotonicity of an expression with respect to a variable on an interval.
     *
     * Auto-infers monotonicity for known functions:
     * - exp: Increasing on all of R
     * - ln: Increasing on R+ (positive reals)
     * - negation (multiply by -1): reverses monotonicity
     *
     * Also checks PropertyStore for user-declared monotonicity.
     *
     * @param expr The expression to analyze
     * @param var The variable with respect to which monotonicity is queried
     * @param interval The interval on which to check monotonicity
     * @return The inferred Monotonicity classification
     */
    Monotonicity infer_monotonicity(const SymbolicExpr& expr, const std::string& var,
                                    const Interval& interval) const;

    /**
     * @brief 沿表达式树传播区间界。
     *
     * 保留精确 Number 点和声明的端点开闭性。
     * 算术运算返回有理数向外包络，实超越运算也使用已证实包络；结果未必最紧。
     *
     * @param expr 用于传播边界的表达式
     * @return 推导区间；事实不足或资源失败时返回 std::nullopt，不返回部分区间
     */
    std::optional<Interval> propagate_bounds(const SymbolicExpr& expr) const;

    /**
     * @brief Apply monotonicity deduction rules when a new relation is added.
     *
     * For relations like x > y where both operands have appropriate domain assumptions,
     * deduces new relations (e.g., ln(x) > ln(y) when both are Positive).
     *
     * @param rel The newly added relation
     * @param store The RelationStore to add deduced relations to
     * @param prop_store The PropertyStore for checking domain assumptions
     * @param depth Current recursion depth (stops at MAX_MONOTONICITY_DEPTH)
     */
    void apply_monotonicity_rules(const Relation& rel, RelationStore& store,
                                  PropertyStore& prop_store, int depth = 0);

    /**
     * @brief Applies monotonicity deductions transactionally.
     *
     * Deduced relations are inserted through RelationStore::add_relation_checked.
     * 派生关系与既有属性冲突时,插入操作以 CasError 原子地返回诊断,
     * 关系存储保持提交前状态.
     */
    Result<void> apply_monotonicity_rules_checked(const Relation& rel,
                                                  RelationStore& store,
                                                  PropertyStore& prop_store,
                                                  int depth = 0);

private:
    static constexpr int MAX_MONOTONICITY_DEPTH = 8;

    struct Impl;
    class DepthGuard;
    std::unique_ptr<Impl> impl_;

    /// @name Sign inference for specific node types
    /// @{
    InferenceTriboolResult infer_add_sign_checked(const AddNode& node, Sign target, ComputationContext& context) const;
    InferenceTriboolResult infer_multiply_sign_checked(const MultiplyNode& node, Sign target, ComputationContext& context) const;
    InferenceTriboolResult infer_power_property_checked(const PowerNode& node, Sign target, ComputationContext& context) const;
    InferenceTriboolResult infer_function_property_checked(const FunctionNode& node, Sign target, ComputationContext& context) const;
    /// @}

    /**
     * @brief Infer sign of an arithmetic expression by checking relational constraints.
     *
     * Examines the RelationStore for GT/GEQ zero patterns on operands of sums
     * and products:
     *   - sum: all operands GT 0 -> Positive; all GEQ 0 -> NonNegative
     *   - product: all operands GT 0 -> Positive
     *   - Also checks: x GT y with y non-negative -> x is Positive
     *
     * @param expr The expression to check
     * @param target The sign property being queried
     * @return True if the target sign can be inferred from relations, Unknown otherwise
     */
    InferenceTriboolResult infer_sign_from_relations_checked(const SymbolicExpr& expr, Sign target, ComputationContext& context) const;

    /// @name Division and subtraction sign inference
    /// @{

    /**
     * @brief Infer the sign of an internally represented division expression.
     *
     * Applies the sign multiplication table for division:
     *   positive / positive -> positive, negative / negative -> positive,
     *   positive / negative -> negative, negative / positive -> negative.
     * Returns Unknown when denominator sign is unknown or zero.
    */
    InferenceTriboolResult infer_division_sign_checked(const MultiplyNode& node, Sign target, ComputationContext& context) const;

    /**
     * @brief Infer the sign of an internally represented subtraction expression.
     *
     * Detects a negated subtrahend and applies subtraction sign rules, such as
     * positive minus negative producing a positive result.
     */
    InferenceTriboolResult infer_subtraction_sign_checked(const AddNode& node, Sign target, ComputationContext& context) const;

    /// @}

    /// @name Domain inference for specific node types
    /// @{
    InferenceTriboolResult infer_add_domain_checked(const AddNode& node, Domain target, ComputationContext& context) const;
    InferenceTriboolResult infer_multiply_domain_checked(const MultiplyNode& node, Domain target, ComputationContext& context) const;
    InferenceTriboolResult infer_power_domain_checked(const PowerNode& node, Domain target, ComputationContext& context) const;
    InferenceTriboolResult infer_function_domain_checked(const FunctionNode& node, Domain target, ComputationContext& context) const;
    /// @}

    /// @name Helper methods for querying sub-expression properties
    /// @{
    InferenceTriboolResult query_sign_of_checked(const SymbolicExpr& expr, Sign sign, ComputationContext& context) const;
    InferenceTriboolResult query_domain_of_checked(const SymbolicExpr& expr, Domain domain, ComputationContext& context) const;
    /// @}
};

} // namespace LMCAS
