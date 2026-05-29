/**
 * @file query_interface.hpp
 * @brief QueryInterface class — unified tri-state query API for property questions on arbitrary expressions.
 *
 * The QueryInterface is the single public entry point for property queries on SymbolicExpr trees.
 * It handles special cases (null root, NaN, Infinity, Matrix, Relational, Logical nodes) and
 * delegates to the InferenceEngine for NumberNode, VariableNode, and composite node queries.
 */
#pragma once

#include "assumption.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"

namespace lamina {

// Forward declaration
class AssumptionContext;
class InferenceEngine;

/**
 * @brief Unified query API for property questions on arbitrary SymbolicExpr trees.
 *
 * Provides query_positive, query_negative, query_nonnegative, query_real,
 * query_integer, and query_nonzero methods that return Tribool results.
 *
 * Dispatch logic:
 *   - Null root node → Unknown
 *   - MatrixNode, RelationalNode, LogicalNode → Unknown
 *   - NumberNode with NaN → False for integer, Unknown for sign queries
 *   - FunctionNode::Infinity → derive sign from infinity's sign, False for integer
 *   - All other nodes → delegate to InferenceEngine
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
    Tribool query_positive(const SymbolicExpr& expr) const;

    /// Query whether the expression is negative (< 0).
    Tribool query_negative(const SymbolicExpr& expr) const;

    /// Query whether the expression is non-negative (>= 0).
    Tribool query_nonnegative(const SymbolicExpr& expr) const;

    /// Query whether the expression is real.
    Tribool query_real(const SymbolicExpr& expr) const;

    /// Query whether the expression is an integer.
    Tribool query_integer(const SymbolicExpr& expr) const;

    /// Query whether the expression is non-zero (!= 0).
    Tribool query_nonzero(const SymbolicExpr& expr) const;

    /// @}

private:
    const AssumptionContext& ctx_;

    /// Check if the root node is a special case that should return Unknown immediately.
    /// Returns true if the node is null, MatrixNode, RelationalNode, or LogicalNode.
    bool is_unhandled_type(const std::shared_ptr<SymbolicNode>& node) const;

    /// Check if a NumberNode holds NaN.
    bool is_nan_number(const NumberNode& node) const;

    /// Check if a FunctionNode represents infinity (FuncType::Infinity).
    bool is_infinity_node(const std::shared_ptr<SymbolicNode>& node) const;

    /// Determine the sign of an infinity expression.
    /// Positive infinity: FunctionNode::Infinity directly.
    /// Negative infinity: MultiplyNode(-1, FunctionNode::Infinity).
    /// Returns +1 for positive infinity, -1 for negative infinity, 0 if indeterminate.
    int get_infinity_sign(const std::shared_ptr<SymbolicNode>& node) const;
};

} // namespace lamina
