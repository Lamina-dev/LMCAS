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

// Forward declaration for Interval (avoid circular include with interval.hpp)
namespace lamina {
struct Interval;
}

namespace lamina {

// Forward declaration — AssumptionContext is implemented in a later task
class AssumptionContext;

using InferenceTriboolResult = Result<Tribool>;
using InferencePeriodResult = Result<std::optional<SymbolicExpr>>;

/**
 * @brief Derives properties of composite expressions by recursively analyzing sub-expressions.
 *
 * The InferenceEngine implements arithmetic inference rules for addition, multiplication,
 * power, and function expressions. It determines sign, domain, and boundedness properties
 * based on the properties of operands/arguments.
 *
 * Sign inference rules:
 *   - Addition: uniform sign → same sign on sum; mixed/unknown → Unknown
 *   - Multiplication: zero → Zero; parity of negatives determines sign; unknown → Unknown
 *   - Power: positive base + real exponent → Positive; even exponent → NonNegative; etc.
 *   - Functions: exp → Positive; abs → NonNegative; sin/cos → Bounded[-1,1]; etc.
 *
 * Domain inference rules:
 *   - Addition: all Integer → Integer; all Real → Real
 *   - Multiplication: all Integer → Integer; all Real/Integer → Real
 *   - Power: Real base + integer exponent → Real
 *   - Functions: exp/sin/cos/abs/sqrt/ln/tan with Real argument → Real
 */
class LAMINA_API InferenceEngine {
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
    /// @{

    Tribool query_positive(const SymbolicExpr& expr) const;
    Tribool query_negative(const SymbolicExpr& expr) const;
    Tribool query_nonnegative(const SymbolicExpr& expr) const;
    Tribool query_nonpositive(const SymbolicExpr& expr) const;
    Tribool query_real(const SymbolicExpr& expr) const;
    Tribool query_integer(const SymbolicExpr& expr) const;
    Tribool query_nonzero(const SymbolicExpr& expr) const;

    InferenceTriboolResult query_positive_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_negative_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_nonnegative_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_nonpositive_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_real_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_integer_checked(const SymbolicExpr& expr) const;
    InferenceTriboolResult query_nonzero_checked(const SymbolicExpr& expr) const;

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
    Tribool query_algebraic(const SymbolicExpr& expr) const;

    /// Checked algebraic query for direct InferenceEngine callers.
    InferenceTriboolResult query_algebraic_checked(const SymbolicExpr& expr) const;

    /**
     * @brief Query whether an expression is transcendental.
     *
     * For variables, checks if the transcendental flag is set in the PropertyStore.
     * For composite expressions, returns Unknown.
     *
     * @param expr The expression to query
     * @return True if transcendental, False if algebraic, Unknown otherwise
     */
    Tribool query_transcendental(const SymbolicExpr& expr) const;

    /// Checked transcendental query for direct InferenceEngine callers.
    InferenceTriboolResult query_transcendental_checked(const SymbolicExpr& expr) const;

    /**
     * @brief Query whether an expression has a finite value/limit.
     *
     * For variables, checks if Finiteness::Finite is declared in the PropertyStore.
     * For composite expressions, returns Unknown.
     *
     * @param expr The expression to query
     * @return True if finite, False if divergent, Unknown otherwise
     */
    Tribool query_finite(const SymbolicExpr& expr) const;

    /// Checked finite-value query for direct InferenceEngine callers.
    InferenceTriboolResult query_finite_checked(const SymbolicExpr& expr) const;

    /**
     * @brief Query whether an expression diverges.
     *
     * For variables, checks if Finiteness::Divergent is declared in the PropertyStore.
     * For composite expressions, returns Unknown.
     *
     * @param expr The expression to query
     * @return True if divergent, False if finite, Unknown otherwise
     */
    Tribool query_divergent(const SymbolicExpr& expr) const;

    /// Checked divergent-value query for direct InferenceEngine callers.
    InferenceTriboolResult query_divergent_checked(const SymbolicExpr& expr) const;

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
     * @brief Query whether an expression is periodic.
     *
     * Returns True for known periodic functions (sin, cos, tan) and symbols
     * declared periodic in the PropertyStore.
     *
     * @param expr The expression to query
     * @return Tribool::True if periodic, False if known non-periodic, Unknown otherwise
     */
    Tribool query_periodic(const SymbolicExpr& expr) const;

    /// Checked periodic query for direct InferenceEngine callers.
    InferenceTriboolResult query_periodic_checked(const SymbolicExpr& expr) const;

    /**
     * @brief Infer the period of an expression.
     *
     * For known periodic functions: sin/cos → 2π, tan → π.
     * For symbols declared periodic, returns the stored period.
     *
     * @param expr The expression to query
     * @return The period as a SymbolicExpr, or std::nullopt if not periodic or unknown
     */
    std::optional<SymbolicExpr> infer_period(const SymbolicExpr& expr) const;

    /// Checked period inference for direct InferenceEngine callers.
    InferencePeriodResult infer_period_checked(const SymbolicExpr& expr) const;

    /**
     * @brief Infer the monotonicity of an expression with respect to a variable on an interval.
     *
     * Auto-infers monotonicity for known functions:
     * - exp: Increasing on all of ℝ
     * - ln: Increasing on ℝ⁺ (positive reals)
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
     * @brief Propagate interval bounds through an expression tree.
     * @param expr The expression to propagate bounds for
     * @return The propagated interval, or std::nullopt if bounds cannot be determined
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
     * @brief Checked monotonicity deduction with explicit relation-store failures.
     *
     * Deduced relations are inserted through RelationStore::add_relation_checked.
     * If any derived relation contradicts stored properties, the failing insert
     * returns CasError instead of partially applying that relation.
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
    Tribool infer_add_sign(const void* node, Sign target) const;
    InferenceTriboolResult infer_add_sign_checked(const void* node, Sign target) const;
    Tribool infer_multiply_sign(const void* node, Sign target) const;
    InferenceTriboolResult infer_multiply_sign_checked(const void* node, Sign target) const;
    Tribool infer_power_property(const void* node, Sign target) const;
    InferenceTriboolResult infer_power_property_checked(const void* node, Sign target) const;
    Tribool infer_function_property(const void* node, Sign target) const;
    InferenceTriboolResult infer_function_property_checked(const void* node, Sign target) const;
    /// @}

    /**
     * @brief Infer sign of an arithmetic expression by checking relational constraints.
     *
     * Examines the RelationStore for GT/GEQ zero patterns on operands of sums
     * and products:
     *   - sum: all operands GT 0 -> Positive; all GEQ 0 -> NonNegative
     *   - product: all operands GT 0 -> Positive
     *   - Also checks: x GT y with y non-negative → x is Positive
     *
     * @param expr The expression to check
     * @param target The sign property being queried
     * @return True if the target sign can be inferred from relations, Unknown otherwise
     */
    Tribool infer_sign_from_relations(const SymbolicExpr& expr, Sign target) const;
    InferenceTriboolResult infer_sign_from_relations_checked(const SymbolicExpr& expr, Sign target) const;

    /// @name Division and subtraction sign inference
    /// @{

    /**
     * @brief Infer the sign of an internally represented division expression.
     *
     * Applies the sign multiplication table for division:
     *   positive ÷ positive → positive, negative ÷ negative → positive,
     *   positive ÷ negative → negative, negative ÷ positive → negative.
     * Returns Unknown when denominator sign is unknown or zero.
    */
    Tribool infer_division_sign(const void* node, Sign target) const;
    InferenceTriboolResult infer_division_sign_checked(const void* node, Sign target) const;

    /**
     * @brief Infer the sign of an internally represented subtraction expression.
     *
     * Detects a negated subtrahend and applies subtraction sign rules, such as
     * positive minus negative producing a positive result.
     */
    Tribool infer_subtraction_sign(const void* node, Sign target) const;
    InferenceTriboolResult infer_subtraction_sign_checked(const void* node, Sign target) const;

    /// @}

    /// @name Domain inference for specific node types
    /// @{
    Tribool infer_add_domain(const void* node, Domain target) const;
    InferenceTriboolResult infer_add_domain_checked(const void* node, Domain target) const;
    Tribool infer_multiply_domain(const void* node, Domain target) const;
    InferenceTriboolResult infer_multiply_domain_checked(const void* node, Domain target) const;
    Tribool infer_power_domain(const void* node, Domain target) const;
    InferenceTriboolResult infer_power_domain_checked(const void* node, Domain target) const;
    Tribool infer_function_domain(const void* node, Domain target) const;
    InferenceTriboolResult infer_function_domain_checked(const void* node, Domain target) const;
    /// @}

    /// @name Helper methods for querying sub-expression properties
    /// @{
    Tribool query_sign_of(const SymbolicExpr& expr, Sign sign) const;
    Tribool query_domain_of(const SymbolicExpr& expr, Domain domain) const;
    InferenceTriboolResult query_sign_of_checked(const SymbolicExpr& expr, Sign sign) const;
    InferenceTriboolResult query_domain_of_checked(const SymbolicExpr& expr, Domain domain) const;
    /// @}
};

} // namespace lamina
