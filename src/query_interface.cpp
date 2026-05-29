/**
 * @file query_interface.cpp
 * @brief Implementation of the QueryInterface class.
 *
 * The QueryInterface is the single public entry point for property queries.
 * It handles special cases (null, NaN, Infinity, Matrix, Relational, Logical)
 * and delegates to the InferenceEngine for all other node types.
 */

#include "query_interface.hpp"
#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include <cmath>

namespace lamina {

// ============================================================
// Construction
// ============================================================

QueryInterface::QueryInterface(const AssumptionContext& ctx)
    : ctx_(ctx) {}

// ============================================================
// Private helpers
// ============================================================

bool QueryInterface::is_unhandled_type(const std::shared_ptr<SymbolicNode>& node) const {
    if (!node) return true;
    if (std::dynamic_pointer_cast<MatrixNode>(node)) return true;
    if (std::dynamic_pointer_cast<RelationalNode>(node)) return true;
    if (std::dynamic_pointer_cast<LogicalNode>(node)) return true;
    return false;
}

bool QueryInterface::is_nan_number(const NumberNode& node) const {
    if (std::holds_alternative<lmmc_real_t>(node.value)) {
        lmmc_real_t v = std::get<lmmc_real_t>(node.value);
        return std::isnan(v);
    }
    return false;
}

bool QueryInterface::is_infinity_node(const std::shared_ptr<SymbolicNode>& node) const {
    // Direct infinity: FunctionNode with FuncType::Infinity
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(node)) {
        return func->type == FunctionNode::FuncType::Infinity;
    }
    // Negative infinity: MultiplyNode(-1, Infinity)
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        if (mul->operands.size() == 2) {
            for (const auto& op : mul->operands) {
                if (auto func = std::dynamic_pointer_cast<FunctionNode>(op)) {
                    if (func->type == FunctionNode::FuncType::Infinity) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

int QueryInterface::get_infinity_sign(const std::shared_ptr<SymbolicNode>& node) const {
    // Direct infinity node → positive infinity
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(node)) {
        if (func->type == FunctionNode::FuncType::Infinity) {
            return +1;
        }
    }
    // Negative infinity: MultiplyNode containing -1 and Infinity
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        bool has_infinity = false;
        bool has_neg_one = false;
        for (const auto& op : mul->operands) {
            if (auto func = std::dynamic_pointer_cast<FunctionNode>(op)) {
                if (func->type == FunctionNode::FuncType::Infinity) {
                    has_infinity = true;
                }
            }
            if (auto num = std::dynamic_pointer_cast<NumberNode>(op)) {
                // Check if it's -1
                if (std::holds_alternative<BigInt>(num->value)) {
                    if (std::get<BigInt>(num->value) == BigInt(-1)) {
                        has_neg_one = true;
                    }
                } else if (std::holds_alternative<Rational>(num->value)) {
                    if (std::get<Rational>(num->value) == Rational(-1)) {
                        has_neg_one = true;
                    }
                } else if (std::holds_alternative<lmmc_real_t>(num->value)) {
                    lmmc_real_t v = std::get<lmmc_real_t>(num->value);
                    if (v == -1.0) {
                        has_neg_one = true;
                    }
                }
            }
        }
        if (has_infinity && has_neg_one) return -1;
        if (has_infinity) return +1; // Infinity multiplied by something else positive
    }
    return 0; // indeterminate
}

// ============================================================
// Public query methods
// ============================================================

Tribool QueryInterface::query_positive(const SymbolicExpr& expr) const {
    // Null root or unhandled types → Unknown (Req 10.10)
    if (is_unhandled_type(expr.root)) {
        return Tribool::Unknown;
    }

    // Check for NaN (Req 10.11): Unknown for sign queries
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (is_nan_number(*num)) {
            return Tribool::Unknown;
        }
    }

    // Check for Infinity (Req 10.11): derive sign from infinity's sign
    if (is_infinity_node(expr.root)) {
        int sign = get_infinity_sign(expr.root);
        if (sign > 0) return Tribool::True;
        if (sign < 0) return Tribool::False;
        return Tribool::Unknown;
    }

    // Delegate to InferenceEngine for all other node types
    // (NumberNode, VariableNode, AddNode, MultiplyNode, PowerNode, FunctionNode)
    InferenceEngine engine(ctx_);
    return engine.query_positive(expr);
}

Tribool QueryInterface::query_negative(const SymbolicExpr& expr) const {
    // Null root or unhandled types → Unknown (Req 10.10)
    if (is_unhandled_type(expr.root)) {
        return Tribool::Unknown;
    }

    // Check for NaN (Req 10.11): Unknown for sign queries
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (is_nan_number(*num)) {
            return Tribool::Unknown;
        }
    }

    // Check for Infinity (Req 10.11): derive sign from infinity's sign
    if (is_infinity_node(expr.root)) {
        int sign = get_infinity_sign(expr.root);
        if (sign < 0) return Tribool::True;
        if (sign > 0) return Tribool::False;
        return Tribool::Unknown;
    }

    // Delegate to InferenceEngine
    InferenceEngine engine(ctx_);
    return engine.query_negative(expr);
}

Tribool QueryInterface::query_nonnegative(const SymbolicExpr& expr) const {
    // Null root or unhandled types → Unknown (Req 10.10)
    if (is_unhandled_type(expr.root)) {
        return Tribool::Unknown;
    }

    // Check for NaN (Req 10.11): Unknown for sign queries
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (is_nan_number(*num)) {
            return Tribool::Unknown;
        }
    }

    // Check for Infinity (Req 10.11): derive from infinity's sign
    if (is_infinity_node(expr.root)) {
        int sign = get_infinity_sign(expr.root);
        if (sign > 0) return Tribool::True;   // +∞ is nonnegative
        if (sign < 0) return Tribool::False;   // -∞ is not nonnegative
        return Tribool::Unknown;
    }

    // Delegate to InferenceEngine
    InferenceEngine engine(ctx_);
    return engine.query_nonnegative(expr);
}

Tribool QueryInterface::query_real(const SymbolicExpr& expr) const {
    // Null root or unhandled types → Unknown (Req 10.10)
    if (is_unhandled_type(expr.root)) {
        return Tribool::Unknown;
    }

    // Check for NaN (Req 10.11): NaN is not a real number in the usual sense,
    // but the requirement says Unknown for sign queries on NaN.
    // For query_real on NaN, return Unknown (it's not a well-defined real value).
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (is_nan_number(*num)) {
            return Tribool::Unknown;
        }
    }

    // Infinity is not a real number (it's extended real)
    if (is_infinity_node(expr.root)) {
        return Tribool::Unknown;
    }

    // Delegate to InferenceEngine
    InferenceEngine engine(ctx_);
    return engine.query_real(expr);
}

Tribool QueryInterface::query_integer(const SymbolicExpr& expr) const {
    // Null root or unhandled types → Unknown (Req 10.10)
    if (is_unhandled_type(expr.root)) {
        return Tribool::Unknown;
    }

    // Check for NaN (Req 10.11): False for integer
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (is_nan_number(*num)) {
            return Tribool::False;
        }
    }

    // Infinity (Req 10.11): False for integer
    if (is_infinity_node(expr.root)) {
        return Tribool::False;
    }

    // Delegate to InferenceEngine
    InferenceEngine engine(ctx_);
    return engine.query_integer(expr);
}

Tribool QueryInterface::query_nonzero(const SymbolicExpr& expr) const {
    // Null root or unhandled types → Unknown (Req 10.10)
    if (is_unhandled_type(expr.root)) {
        return Tribool::Unknown;
    }

    // Check for NaN (Req 10.11): Unknown for sign queries
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (is_nan_number(*num)) {
            return Tribool::Unknown;
        }
    }

    // Infinity is nonzero
    if (is_infinity_node(expr.root)) {
        return Tribool::True;
    }

    // Delegate to InferenceEngine
    InferenceEngine engine(ctx_);
    return engine.query_nonzero(expr);
}

} // namespace lamina
