/**
 * @file inference_engine.hpp
 * @brief InferenceEngine class for deriving properties of composite expressions.
 *
 * Analyzes AddNode, MultiplyNode, PowerNode, and FunctionNode expressions to infer
 * sign, domain, and boundedness properties from the properties of their sub-expressions.
 */
#pragma once

#include "assumption.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "relation_store.hpp"
#include <optional>
#include <unordered_set>

// Forward declaration for Interval (avoid circular include with interval.hpp)
namespace lamina {
struct Interval;
}

namespace lamina {

// Forward declaration — AssumptionContext is implemented in a later task
class AssumptionContext;

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

private:
    static constexpr int MAX_MONOTONICITY_DEPTH = 8;

    const AssumptionContext& ctx_;

    /// Maximum recursion depth for property queries (default 32)
    int max_depth_ = 32;

    /// Set of node pointers currently being visited (cycle detection).
    /// Uses raw pointers for identity-based comparison (not structural equality).
    mutable std::unordered_set<const SymbolicNode*> visited_;

    /// Current recursion depth counter for depth limiting.
    mutable int current_depth_ = 0;

    /**
     * @brief RAII guard for recursive query depth tracking and cycle detection.
     *
     * On construction: increments current_depth_ and inserts the node pointer into visited_.
     * On destruction: decrements current_depth_ and removes the node pointer from visited_.
     * When depth returns to 0, the visited set is cleared (top-level query completion).
     */
    class DepthGuard {
    public:
        DepthGuard(const InferenceEngine& engine, const SymbolicNode* node);
        ~DepthGuard();

        /// Returns true if the query should be aborted (cycle or depth exceeded)
        bool should_abort() const { return abort_; }

        // Non-copyable, non-movable
        DepthGuard(const DepthGuard&) = delete;
        DepthGuard& operator=(const DepthGuard&) = delete;

    private:
        const InferenceEngine& engine_;
        const SymbolicNode* node_;
        bool abort_ = false;
        bool inserted_ = false; ///< Whether the node was newly inserted into visited_
    };

    /// @name Sign inference for specific node types
    /// @{
    Tribool infer_add_sign(const AddNode& node, Sign target) const;
    Tribool infer_multiply_sign(const MultiplyNode& node, Sign target) const;
    Tribool infer_power_property(const PowerNode& node, Sign target) const;
    Tribool infer_function_property(const FunctionNode& node, Sign target) const;
    /// @}

    /**
     * @brief Infer sign of an arithmetic expression by checking relational constraints.
     *
     * Examines the RelationStore for GT/GEQ zero patterns on operands of AddNode
     * and MultiplyNode expressions:
     *   - AddNode: all operands GT 0 → Positive; all GEQ 0 → NonNegative
     *   - MultiplyNode: all operands GT 0 → Positive
     *   - Also checks: x GT y with y non-negative → x is Positive
     *
     * @param expr The expression to check
     * @param target The sign property being queried
     * @return True if the target sign can be inferred from relations, Unknown otherwise
     */
    Tribool infer_sign_from_relations(const SymbolicExpr& expr, Sign target) const;

    /// @name Division and subtraction sign inference
    /// @{

    /**
     * @brief Infer sign of a division expression represented as MultiplyNode([num, PowerNode(den, -1)]).
     *
     * Applies the sign multiplication table for division:
     *   positive ÷ positive → positive, negative ÷ negative → positive,
     *   positive ÷ negative → negative, negative ÷ positive → negative.
     * Returns Unknown when denominator sign is unknown or zero.
     */
    Tribool infer_division_sign(const MultiplyNode& node, Sign target) const;

    /**
     * @brief Infer sign of a subtraction expression represented as AddNode with negated operands.
     *
     * Detects operands of the form MultiplyNode([NumberNode(-1), subtrahend]) and applies
     * subtraction sign rules (e.g., positive minus negative → positive).
     */
    Tribool infer_subtraction_sign(const AddNode& node, Sign target) const;

    /// @}

    /// @name Domain inference for specific node types
    /// @{
    Tribool infer_add_domain(const AddNode& node, Domain target) const;
    Tribool infer_multiply_domain(const MultiplyNode& node, Domain target) const;
    Tribool infer_power_domain(const PowerNode& node, Domain target) const;
    Tribool infer_function_domain(const FunctionNode& node, Domain target) const;
    /// @}

    /// @name Helper methods for querying sub-expression properties
    /// @{
    Tribool query_sign_of(const SymbolicExpr& expr, Sign sign) const;
    Tribool query_domain_of(const SymbolicExpr& expr, Domain domain) const;
    /// @}
};

} // namespace lamina
