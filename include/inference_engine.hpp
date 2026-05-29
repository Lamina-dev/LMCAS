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

    /// @}

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

    /// @name Sign inference for specific node types
    /// @{
    Tribool infer_add_sign(const AddNode& node, Sign target) const;
    Tribool infer_multiply_sign(const MultiplyNode& node, Sign target) const;
    Tribool infer_power_property(const PowerNode& node, Sign target) const;
    Tribool infer_function_property(const FunctionNode& node, Sign target) const;
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
