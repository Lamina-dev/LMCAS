/**
 * @file inference_engine.cpp
 * @brief Implementation of the InferenceEngine class.
 *
 * Implements arithmetic inference rules for addition, multiplication,
 * power expressions, and built-in functions.
 */

#include "inference_engine.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "interval.hpp"
#include <algorithm>
#include <cmath>
#include <climits>

namespace lamina {

// ============================================================
// Construction
// ============================================================

InferenceEngine::InferenceEngine(const AssumptionContext& ctx)
    : ctx_(ctx) {}

// ============================================================
// Helper: query sign/domain of a sub-expression
// ============================================================

// These helpers will eventually delegate through the AssumptionContext's
// QueryInterface. For now, they provide the interface that the inference
// logic uses. Once AssumptionContext is implemented, these will call
// ctx_.is_positive(expr), ctx_.is_integer(expr), etc.

Tribool InferenceEngine::query_sign_of(const SymbolicExpr& expr, Sign sign) const {
    // Delegate to the public query methods which will be wired through
    // the AssumptionContext/QueryInterface in a later task.
    // For now, this provides the internal dispatch.
    switch (sign) {
        case Sign::Positive:    return query_positive(expr);
        case Sign::Negative:    return query_negative(expr);
        case Sign::NonNegative: return query_nonnegative(expr);
        case Sign::NonPositive: return query_nonpositive(expr);
        case Sign::Zero:
            // Zero means both NonNegative and NonPositive
            {
                auto nn = query_nonnegative(expr);
                auto np = query_nonpositive(expr);
                if (nn == Tribool::True && np == Tribool::True)
                    return Tribool::True;
                if (nn == Tribool::False || np == Tribool::False)
                    return Tribool::False;
                return Tribool::Unknown;
            }
        case Sign::NonZero:     return query_nonzero(expr);
    }
    return Tribool::Unknown;
}

Tribool InferenceEngine::query_domain_of(const SymbolicExpr& expr, Domain domain) const {
    switch (domain) {
        case Domain::Integer:  return query_integer(expr);
        case Domain::Real:     return query_real(expr);
        default:               return Tribool::Unknown;
    }
}

// ============================================================
// Public query methods — dispatch based on node type
// ============================================================

// These public methods are called by the QueryInterface for composite nodes.
// They inspect the root node type and dispatch to the appropriate inference method.

Tribool InferenceEngine::query_positive(const SymbolicExpr& expr) const {
    if (!expr.root) return Tribool::Unknown;

    // Handle NumberNode: determine sign directly from numeric value
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (num->is_positive()) return Tribool::True;
        if (num->is_zero()) return Tribool::False;
        // Check if negative — if not negative and not zero, it must be positive
        if (std::holds_alternative<BigInt>(num->value)) {
            return std::get<BigInt>(num->value).IsNegative() ? Tribool::False : Tribool::True;
        }
        if (std::holds_alternative<Rational>(num->value)) {
            return std::get<Rational>(num->value) < Rational(0) ? Tribool::False : Tribool::True;
        }
        if (std::holds_alternative<lmmc_real_t>(num->value)) {
            lmmc_real_t v = std::get<lmmc_real_t>(num->value);
            if (!std::isfinite(v)) return Tribool::Unknown;
            if (v > 0.0) return Tribool::True;
            if (v < 0.0) return Tribool::False;
            return Tribool::False; // zero
        }
        return Tribool::Unknown;
    }

    // Handle VariableNode: check PropertyStore in the AssumptionContext
    if (auto var = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
        const auto& props = ctx_.current_properties();
        if (props.has_sign(var->name, Sign::Positive)) return Tribool::True;
        // If NonPositive or Negative or Zero, then not positive
        if (props.has_sign(var->name, Sign::Negative) ||
            props.has_sign(var->name, Sign::Zero) ||
            props.has_sign(var->name, Sign::NonPositive)) return Tribool::False;
        return Tribool::Unknown;
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(expr.root)) {
        return infer_add_sign(*add, Sign::Positive);
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        return infer_multiply_sign(*mul, Sign::Positive);
    }
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        return infer_power_property(*pow, Sign::Positive);
    }
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(expr.root)) {
        return infer_function_property(*func, Sign::Positive);
    }
    return Tribool::Unknown;
}

Tribool InferenceEngine::query_negative(const SymbolicExpr& expr) const {
    if (!expr.root) return Tribool::Unknown;

    // Handle NumberNode
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (num->is_zero()) return Tribool::False;
        if (std::holds_alternative<BigInt>(num->value)) {
            return std::get<BigInt>(num->value).IsNegative() ? Tribool::True : Tribool::False;
        }
        if (std::holds_alternative<Rational>(num->value)) {
            const auto& r = std::get<Rational>(num->value);
            if (r < Rational(0)) return Tribool::True;
            return Tribool::False;
        }
        if (std::holds_alternative<lmmc_real_t>(num->value)) {
            lmmc_real_t v = std::get<lmmc_real_t>(num->value);
            if (!std::isfinite(v)) return Tribool::Unknown;
            if (v < 0.0) return Tribool::True;
            return Tribool::False;
        }
        return Tribool::Unknown;
    }

    // Handle VariableNode
    if (auto var = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
        const auto& props = ctx_.current_properties();
        if (props.has_sign(var->name, Sign::Negative)) return Tribool::True;
        if (props.has_sign(var->name, Sign::Positive) ||
            props.has_sign(var->name, Sign::Zero) ||
            props.has_sign(var->name, Sign::NonNegative)) return Tribool::False;
        return Tribool::Unknown;
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(expr.root)) {
        return infer_add_sign(*add, Sign::Negative);
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        return infer_multiply_sign(*mul, Sign::Negative);
    }
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        return infer_power_property(*pow, Sign::Negative);
    }
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(expr.root)) {
        return infer_function_property(*func, Sign::Negative);
    }
    return Tribool::Unknown;
}

Tribool InferenceEngine::query_nonnegative(const SymbolicExpr& expr) const {
    if (!expr.root) return Tribool::Unknown;

    // Handle NumberNode
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (num->is_zero()) return Tribool::True;
        if (num->is_positive()) return Tribool::True;
        if (std::holds_alternative<BigInt>(num->value)) {
            return std::get<BigInt>(num->value).IsNegative() ? Tribool::False : Tribool::True;
        }
        if (std::holds_alternative<Rational>(num->value)) {
            return std::get<Rational>(num->value) < Rational(0) ? Tribool::False : Tribool::True;
        }
        if (std::holds_alternative<lmmc_real_t>(num->value)) {
            lmmc_real_t v = std::get<lmmc_real_t>(num->value);
            if (!std::isfinite(v)) return Tribool::Unknown;
            return v >= 0.0 ? Tribool::True : Tribool::False;
        }
        return Tribool::Unknown;
    }

    // Handle VariableNode
    if (auto var = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
        const auto& props = ctx_.current_properties();
        if (props.has_sign(var->name, Sign::NonNegative)) return Tribool::True;
        if (props.has_sign(var->name, Sign::Negative)) return Tribool::False;
        return Tribool::Unknown;
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(expr.root)) {
        return infer_add_sign(*add, Sign::NonNegative);
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        return infer_multiply_sign(*mul, Sign::NonNegative);
    }
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        return infer_power_property(*pow, Sign::NonNegative);
    }
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(expr.root)) {
        return infer_function_property(*func, Sign::NonNegative);
    }
    return Tribool::Unknown;
}

Tribool InferenceEngine::query_nonpositive(const SymbolicExpr& expr) const {
    if (!expr.root) return Tribool::Unknown;

    // Handle NumberNode
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (num->is_zero()) return Tribool::True;
        if (num->is_positive()) return Tribool::False;
        if (std::holds_alternative<BigInt>(num->value)) {
            return std::get<BigInt>(num->value).IsNegative() ? Tribool::True : Tribool::False;
        }
        if (std::holds_alternative<Rational>(num->value)) {
            return std::get<Rational>(num->value) > Rational(0) ? Tribool::False : Tribool::True;
        }
        if (std::holds_alternative<lmmc_real_t>(num->value)) {
            lmmc_real_t v = std::get<lmmc_real_t>(num->value);
            if (!std::isfinite(v)) return Tribool::Unknown;
            return v <= 0.0 ? Tribool::True : Tribool::False;
        }
        return Tribool::Unknown;
    }

    // Handle VariableNode
    if (auto var = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
        const auto& props = ctx_.current_properties();
        if (props.has_sign(var->name, Sign::NonPositive)) return Tribool::True;
        if (props.has_sign(var->name, Sign::Positive)) return Tribool::False;
        return Tribool::Unknown;
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(expr.root)) {
        return infer_add_sign(*add, Sign::NonPositive);
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        return infer_multiply_sign(*mul, Sign::NonPositive);
    }
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        return infer_power_property(*pow, Sign::NonPositive);
    }
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(expr.root)) {
        return infer_function_property(*func, Sign::NonPositive);
    }
    return Tribool::Unknown;
}

Tribool InferenceEngine::query_real(const SymbolicExpr& expr) const {
    if (!expr.root) return Tribool::Unknown;

    // Handle NumberNode: finite numbers are Real
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (std::holds_alternative<BigInt>(num->value)) return Tribool::True;
        if (std::holds_alternative<Rational>(num->value)) return Tribool::True;
        if (std::holds_alternative<lmmc_real_t>(num->value)) {
            lmmc_real_t v = std::get<lmmc_real_t>(num->value);
            return std::isfinite(v) ? Tribool::True : Tribool::Unknown;
        }
        return Tribool::Unknown;
    }

    // Handle VariableNode
    if (auto var = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
        const auto& props = ctx_.current_properties();
        if (props.has_domain(var->name, Domain::Real)) return Tribool::True;
        return Tribool::Unknown;
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(expr.root)) {
        return infer_add_domain(*add, Domain::Real);
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        return infer_multiply_domain(*mul, Domain::Real);
    }
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        return infer_power_domain(*pow, Domain::Real);
    }
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(expr.root)) {
        return infer_function_domain(*func, Domain::Real);
    }
    return Tribool::Unknown;
}

Tribool InferenceEngine::query_integer(const SymbolicExpr& expr) const {
    if (!expr.root) return Tribool::Unknown;

    // Handle NumberNode: BigInt values are Integer
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (std::holds_alternative<BigInt>(num->value)) return Tribool::True;
        if (std::holds_alternative<Rational>(num->value)) {
            const auto& r = std::get<Rational>(num->value);
            // Integer if denominator is 1
            return r.is_integer() ? Tribool::True : Tribool::False;
        }
        if (std::holds_alternative<lmmc_real_t>(num->value)) {
            lmmc_real_t v = std::get<lmmc_real_t>(num->value);
            if (!std::isfinite(v)) return Tribool::False;
            return (v == std::floor(v)) ? Tribool::True : Tribool::False;
        }
        return Tribool::Unknown;
    }

    // Handle VariableNode
    if (auto var = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
        const auto& props = ctx_.current_properties();
        if (props.has_domain(var->name, Domain::Integer)) return Tribool::True;
        return Tribool::Unknown;
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(expr.root)) {
        return infer_add_domain(*add, Domain::Integer);
    }
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        return infer_multiply_domain(*mul, Domain::Integer);
    }
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        return infer_power_domain(*pow, Domain::Integer);
    }
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(expr.root)) {
        return infer_function_domain(*func, Domain::Integer);
    }
    return Tribool::Unknown;
}

Tribool InferenceEngine::query_nonzero(const SymbolicExpr& expr) const {
    if (!expr.root) return Tribool::Unknown;

    // Handle NumberNode
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        if (num->is_zero()) return Tribool::False;
        if (num->is_positive()) return Tribool::True;
        // Check if it's a non-zero number
        if (std::holds_alternative<BigInt>(num->value)) {
            return (std::get<BigInt>(num->value) == BigInt(0)) ? Tribool::False : Tribool::True;
        }
        if (std::holds_alternative<Rational>(num->value)) {
            return (std::get<Rational>(num->value) == Rational(0)) ? Tribool::False : Tribool::True;
        }
        if (std::holds_alternative<lmmc_real_t>(num->value)) {
            lmmc_real_t v = std::get<lmmc_real_t>(num->value);
            if (!std::isfinite(v)) return Tribool::Unknown;
            return (v != 0.0) ? Tribool::True : Tribool::False;
        }
        return Tribool::Unknown;
    }

    // Handle VariableNode
    if (auto var = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
        const auto& props = ctx_.current_properties();
        if (props.has_sign(var->name, Sign::NonZero)) return Tribool::True;
        if (props.has_sign(var->name, Sign::Zero)) return Tribool::False;
        return Tribool::Unknown;
    }

    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        return infer_multiply_sign(*mul, Sign::NonZero);
    }
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        return infer_power_property(*pow, Sign::NonZero);
    }
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(expr.root)) {
        return infer_function_property(*func, Sign::NonZero);
    }
    // For addition, NonZero is hard to determine in general
    // (positive + positive = nonzero, but that's covered by positive inference)
    if (auto add = std::dynamic_pointer_cast<AddNode>(expr.root)) {
        // If the sum is positive or negative, it's nonzero
        auto pos = infer_add_sign(*add, Sign::Positive);
        if (pos == Tribool::True) return Tribool::True;
        auto neg = infer_add_sign(*add, Sign::Negative);
        if (neg == Tribool::True) return Tribool::True;
        return Tribool::Unknown;
    }
    return Tribool::Unknown;
}

// ============================================================
// Addition sign inference
// ============================================================

Tribool InferenceEngine::infer_add_sign(const AddNode& node, Sign target) const {
    if (node.operands.empty()) return Tribool::Unknown;

    // For sign inference on addition:
    // - If all operands have the target sign property → sum has that property
    // - If any operand has Unknown → result is Unknown
    // - If operands have mixed definite signs → result is Unknown

    bool all_have_property = true;

    for (const auto& operand : node.operands) {
        SymbolicExpr op_expr;
        op_expr.root = operand;

        Tribool op_result = Tribool::Unknown;
        switch (target) {
            case Sign::Positive:    op_result = query_positive(op_expr); break;
            case Sign::Negative:    op_result = query_negative(op_expr); break;
            case Sign::NonNegative: op_result = query_nonnegative(op_expr); break;
            case Sign::NonPositive: op_result = query_nonpositive(op_expr); break;
            default:                return Tribool::Unknown;
        }

        if (op_result == Tribool::Unknown) {
            return Tribool::Unknown;
        }
        if (op_result == Tribool::False) {
            all_have_property = false;
        }
    }

    if (all_have_property) {
        return Tribool::True;
    }

    return Tribool::Unknown;
}

// ============================================================
// Addition domain inference
// ============================================================

Tribool InferenceEngine::infer_add_domain(const AddNode& node, Domain target) const {
    if (node.operands.empty()) return Tribool::Unknown;

    // For domain inference on addition:
    // - All Integer → Integer (closure under addition)
    // - All Real → Real (closure under addition)
    // - Integer is a subset of Real, so all Integer also implies all Real

    for (const auto& operand : node.operands) {
        SymbolicExpr op_expr;
        op_expr.root = operand;

        Tribool op_result = Tribool::Unknown;
        switch (target) {
            case Domain::Integer: op_result = query_integer(op_expr); break;
            case Domain::Real:    op_result = query_real(op_expr); break;
            default:              return Tribool::Unknown;
        }

        if (op_result != Tribool::True) {
            // For Real target: if operand is Integer, it's also Real
            if (target == Domain::Real) {
                Tribool int_result = query_integer(op_expr);
                if (int_result == Tribool::True) {
                    continue; // Integer implies Real
                }
            }
            return Tribool::Unknown;
        }
    }

    return Tribool::True;
}

// ============================================================
// Multiplication sign inference
// ============================================================

Tribool InferenceEngine::infer_multiply_sign(const MultiplyNode& node, Sign target) const {
    if (node.operands.empty()) return Tribool::Unknown;

    // Step 1: Check if any operand is Zero → product is Zero
    bool has_zero = false;
    for (const auto& operand : node.operands) {
        SymbolicExpr op_expr;
        op_expr.root = operand;

        // Check if operand is zero (both nonnegative and nonpositive)
        Tribool nn = query_nonnegative(op_expr);
        Tribool np = query_nonpositive(op_expr);
        if (nn == Tribool::True && np == Tribool::True) {
            has_zero = true;
            break;
        }
        // Also check the node directly
        if (operand->is_zero()) {
            has_zero = true;
            break;
        }
    }

    if (has_zero) {
        // Product is Zero
        switch (target) {
            case Sign::Zero:        return Tribool::True;
            case Sign::NonNegative: return Tribool::True;
            case Sign::NonPositive: return Tribool::True;
            case Sign::Positive:    return Tribool::False;
            case Sign::Negative:    return Tribool::False;
            case Sign::NonZero:     return Tribool::False;
        }
    }

    // Step 2: Handle NonZero query — all NonZero → product is NonZero
    if (target == Sign::NonZero) {
        for (const auto& operand : node.operands) {
            SymbolicExpr op_expr;
            op_expr.root = operand;
            Tribool nz = query_nonzero(op_expr);
            if (nz != Tribool::True) {
                // Also check if operand is positive or negative (both imply nonzero)
                Tribool pos = query_positive(op_expr);
                if (pos == Tribool::True) continue;
                Tribool neg = query_negative(op_expr);
                if (neg == Tribool::True) continue;
                return Tribool::Unknown;
            }
        }
        return Tribool::True;
    }

    // Step 3: Determine sign by counting negatives
    // We need each operand to have a definite sign (Positive or Negative)
    // or at least NonNegative/NonPositive for weaker results.
    int negative_count = 0;
    bool all_definite = true;       // all are Positive or Negative
    bool all_nonneg_or_nonpos = true; // all are at least NonNeg or NonPos
    bool has_unknown = false;

    for (const auto& operand : node.operands) {
        SymbolicExpr op_expr;
        op_expr.root = operand;

        Tribool pos = query_positive(op_expr);
        Tribool neg = query_negative(op_expr);

        if (pos == Tribool::True) {
            // Positive operand — doesn't change sign
            continue;
        } else if (neg == Tribool::True) {
            // Negative operand — flips sign
            negative_count++;
            continue;
        }

        // Not definitively Positive or Negative
        all_definite = false;

        // Check weaker properties
        Tribool nn = query_nonnegative(op_expr);
        Tribool np = query_nonpositive(op_expr);

        if (nn == Tribool::True) {
            // NonNegative but not Positive (could be zero)
            continue;
        } else if (np == Tribool::True) {
            // NonPositive but not Negative (could be zero)
            negative_count++;
            continue;
        }

        // Unknown sign
        all_nonneg_or_nonpos = false;
        has_unknown = true;
    }

    if (has_unknown) {
        return Tribool::Unknown;
    }

    // Determine result based on parity of negative count
    bool even_negatives = (negative_count % 2 == 0);

    if (all_definite) {
        // All operands are definitively Positive or Negative
        switch (target) {
            case Sign::Positive:
                return even_negatives ? Tribool::True : Tribool::False;
            case Sign::Negative:
                return even_negatives ? Tribool::False : Tribool::True;
            case Sign::NonNegative:
                return even_negatives ? Tribool::True : Tribool::False;
            case Sign::NonPositive:
                return even_negatives ? Tribool::False : Tribool::True;
            case Sign::NonZero:
                return Tribool::True; // product of nonzero values is nonzero
            default:
                return Tribool::Unknown;
        }
    }

    if (all_nonneg_or_nonpos) {
        // All operands are at least NonNegative or NonPositive (some may be zero)
        // Req 6.8: even negatives + all remaining NonNeg (with at least one Pos/Neg) → NonNeg
        // Req 6.9: odd negatives + all remaining NonNeg (with at least one Pos/Neg) → NonPos
        switch (target) {
            case Sign::NonNegative:
                return even_negatives ? Tribool::True : Tribool::False;
            case Sign::NonPositive:
                return even_negatives ? Tribool::False : Tribool::True;
            case Sign::Positive:
                return Tribool::Unknown; // could be zero
            case Sign::Negative:
                return Tribool::Unknown; // could be zero
            default:
                return Tribool::Unknown;
        }
    }

    return Tribool::Unknown;
}

// ============================================================
// Multiplication domain inference
// ============================================================

Tribool InferenceEngine::infer_multiply_domain(const MultiplyNode& node, Domain target) const {
    if (node.operands.empty()) return Tribool::Unknown;

    if (target == Domain::Integer) {
        // All Integer → product is Integer
        for (const auto& operand : node.operands) {
            SymbolicExpr op_expr;
            op_expr.root = operand;
            Tribool result = query_integer(op_expr);
            if (result != Tribool::True) {
                return Tribool::Unknown;
            }
        }
        return Tribool::True;
    }

    if (target == Domain::Real) {
        // All Real (or Integer, since Integer ⊂ Real) → product is Real
        for (const auto& operand : node.operands) {
            SymbolicExpr op_expr;
            op_expr.root = operand;
            Tribool real_result = query_real(op_expr);
            if (real_result == Tribool::True) continue;
            // Integer implies Real
            Tribool int_result = query_integer(op_expr);
            if (int_result == Tribool::True) continue;
            return Tribool::Unknown;
        }
        return Tribool::True;
    }

    return Tribool::Unknown;
}

// ============================================================
// Power expression inference
// ============================================================

/**
 * @brief Helper: check if a NumberNode holds an integer value.
 * Works for BigInt (always integer), Rational (if denominator == 1),
 * and double (if finite and equal to its floor).
 */
static bool is_integer_number(const NumberNode& num) {
    if (std::holds_alternative<BigInt>(num.value)) {
        return true;
    }
    if (std::holds_alternative<Rational>(num.value)) {
        return std::get<Rational>(num.value).is_integer();
    }
    // double case
    double v = std::get<lmmc_real_t>(num.value);
    return std::isfinite(v) && v == std::floor(v);
}

/**
 * @brief Helper: check if a NumberNode holds an even integer value.
 */
static bool is_even_integer_number(const NumberNode& num) {
    if (!is_integer_number(num)) return false;

    if (std::holds_alternative<BigInt>(num.value)) {
        return std::get<BigInt>(num.value).is_even();
    }
    if (std::holds_alternative<Rational>(num.value)) {
        // Rational with denominator 1 — check numerator
        BigInt n = std::get<Rational>(num.value).get_numerator();
        return n.is_even();
    }
    // double case — guard against overflow for very large values
    double v = std::get<lmmc_real_t>(num.value);
    if (std::fabs(v) > static_cast<double>(LLONG_MAX)) {
        return std::fmod(v, 2.0) == 0.0;
    }
    long long iv = static_cast<long long>(v);
    return (iv % 2) == 0;
}

/**
 * @brief Helper: check if a NumberNode holds a positive integer value (> 0).
 */
static bool is_positive_integer_number(const NumberNode& num) {
    if (!is_integer_number(num)) return false;

    if (std::holds_alternative<BigInt>(num.value)) {
        const auto& b = std::get<BigInt>(num.value);
        return !b.IsNegative() && !(b == BigInt(0));
    }
    if (std::holds_alternative<Rational>(num.value)) {
        BigInt n = std::get<Rational>(num.value).get_numerator();
        return !n.IsNegative() && !(n == BigInt(0));
    }
    // double case
    double v = std::get<lmmc_real_t>(num.value);
    return v > 0.0;
}

/**
 * @brief Helper: check if a NumberNode holds the value 0.
 */
static bool is_zero_number(const NumberNode& num) {
    return num.is_zero();
}

Tribool InferenceEngine::infer_power_property(const PowerNode& node, Sign target) const {
    // Wrap base and exponent as SymbolicExpr for querying
    SymbolicExpr base_expr;
    base_expr.root = node.base;
    SymbolicExpr exp_expr;
    exp_expr.root = node.exponent;

    // Check if exponent is a NumberNode (needed for several rules)
    auto exp_num = std::dynamic_pointer_cast<NumberNode>(node.exponent);

    // Rule (Req 7.6): NonNegative base + exponent == 0 → Positive
    // x^0 = 1 for any non-negative x (including 0^0 = 1 by convention here)
    if (exp_num && is_zero_number(*exp_num)) {
        Tribool base_nn = query_nonnegative(base_expr);
        if (base_nn == Tribool::True) {
            // Result is Positive (value is 1)
            switch (target) {
                case Sign::Positive:    return Tribool::True;
                case Sign::NonNegative: return Tribool::True;
                case Sign::NonZero:     return Tribool::True;
                case Sign::Negative:    return Tribool::False;
                case Sign::NonPositive: return Tribool::False;
                case Sign::Zero:        return Tribool::False;
            }
        }
    }

    // Rule (Req 7.1): Positive base + Real exponent → Positive
    {
        Tribool base_pos = query_positive(base_expr);
        if (base_pos == Tribool::True) {
            Tribool exp_real = query_real(exp_expr);
            if (exp_real == Tribool::True) {
                switch (target) {
                    case Sign::Positive:    return Tribool::True;
                    case Sign::NonNegative: return Tribool::True;
                    case Sign::NonZero:     return Tribool::True;
                    case Sign::Negative:    return Tribool::False;
                    case Sign::NonPositive: return Tribool::False;
                    case Sign::Zero:        return Tribool::False;
                }
            }
        }
    }

    // Rule (Req 7.3): Real base + even integer exponent → NonNegative
    if (exp_num && is_even_integer_number(*exp_num)) {
        Tribool base_real = query_real(base_expr);
        if (base_real == Tribool::True) {
            switch (target) {
                case Sign::NonNegative: return Tribool::True;
                case Sign::Negative:    return Tribool::False;
                // Cannot determine Positive (base could be 0)
                // Cannot determine NonPositive (result is >= 0, not necessarily <= 0)
                default: break;
            }
        }
    }

    // Rule (Req 7.2): NonNegative base + positive integer exponent → NonNegative
    if (exp_num && is_positive_integer_number(*exp_num)) {
        Tribool base_nn = query_nonnegative(base_expr);
        if (base_nn == Tribool::True) {
            switch (target) {
                case Sign::NonNegative: return Tribool::True;
                case Sign::Negative:    return Tribool::False;
                default: break;
            }
        }
    }

    // Rule (Req 7.5): NonZero base + integer exponent → NonZero
    if (exp_num && is_integer_number(*exp_num)) {
        Tribool base_nz = query_nonzero(base_expr);
        if (base_nz == Tribool::True) {
            switch (target) {
                case Sign::NonZero: return Tribool::True;
                case Sign::Zero:    return Tribool::False;
                default: break;
            }
        }
    }

    // Default (Req 7.7): cannot determine property
    return Tribool::Unknown;
}

Tribool InferenceEngine::infer_power_domain(const PowerNode& node, Domain target) const {
    // Wrap base as SymbolicExpr for querying
    SymbolicExpr base_expr;
    base_expr.root = node.base;

    // Check if exponent is a NumberNode with integer value
    auto exp_num = std::dynamic_pointer_cast<NumberNode>(node.exponent);

    // Rule (Req 7.4): Real base + integer exponent → Real
    if (exp_num && is_integer_number(*exp_num)) {
        Tribool base_real = query_real(base_expr);
        if (base_real == Tribool::True) {
            switch (target) {
                case Domain::Real: return Tribool::True;
                // Integer base + integer exponent doesn't necessarily give integer
                // (e.g., 2^(-1) = 0.5), so we only infer Real here
                default: break;
            }
        }

        // If base is Integer and exponent is a non-negative integer, result is Integer
        if (target == Domain::Integer) {
            Tribool base_int = query_integer(base_expr);
            if (base_int == Tribool::True && is_positive_integer_number(*exp_num)) {
                return Tribool::True;
            }
            // Also handle exponent == 0: x^0 = 1 which is integer
            if (base_int == Tribool::True && is_zero_number(*exp_num)) {
                return Tribool::True;
            }
        }
    }

    // Default (Req 7.7): cannot determine domain
    return Tribool::Unknown;
}

// ============================================================
// Function inference
// ============================================================

Tribool InferenceEngine::infer_function_property(const FunctionNode& node, Sign target) const {
    // Need at least one argument for all built-in function rules
    if (node.arguments.empty()) return Tribool::Unknown;

    SymbolicExpr arg_expr;
    arg_expr.root = node.arguments[0];

    switch (node.type) {
        case FunctionNode::FuncType::Exp: {
            // exp(Real) → Positive
            // exp is always positive for real arguments
            Tribool arg_real = query_real(arg_expr);
            if (arg_real != Tribool::True) return Tribool::Unknown;

            switch (target) {
                case Sign::Positive:    return Tribool::True;
                case Sign::NonNegative: return Tribool::True;  // Positive implies NonNegative
                case Sign::NonZero:     return Tribool::True;  // Positive implies NonZero
                case Sign::Negative:    return Tribool::False;
                case Sign::NonPositive: return Tribool::False;
                case Sign::Zero:        return Tribool::False;
            }
            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::Sin: {
            // sin(Real) → Real, Bounded[-1,1]
            // sin can be positive, negative, or zero — cannot determine sign in general
            Tribool arg_real = query_real(arg_expr);
            if (arg_real != Tribool::True) return Tribool::Unknown;

            // sin(Real) is bounded in [-1,1], so it's both NonNegative? No.
            // sin can be negative, so we can't say NonNegative or Positive.
            // We can't determine sign without knowing the specific argument value.
            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::Cos: {
            // cos(Real) → Real, Bounded[-1,1]
            // cos can be positive, negative, or zero — cannot determine sign in general
            Tribool arg_real = query_real(arg_expr);
            if (arg_real != Tribool::True) return Tribool::Unknown;

            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::Abs: {
            // abs(Real) → NonNegative, Real
            Tribool arg_real = query_real(arg_expr);
            if (arg_real != Tribool::True) return Tribool::Unknown;

            switch (target) {
                case Sign::NonNegative: return Tribool::True;
                case Sign::Negative:    return Tribool::False;
                case Sign::NonPositive: return Tribool::Unknown; // abs could be zero
                case Sign::Positive:    return Tribool::Unknown; // abs could be zero
                case Sign::Zero:        return Tribool::Unknown; // depends on argument
                case Sign::NonZero:     return Tribool::Unknown; // depends on argument
            }
            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::Ln: {
            // ln(Positive) → Real
            // ln can be positive (x>1), negative (0<x<1), or zero (x=1)
            // Cannot determine sign in general
            Tribool arg_positive = query_positive(arg_expr);
            if (arg_positive != Tribool::True) return Tribool::Unknown;

            // ln(Positive) has no definite sign property
            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::Sqrt: {
            // sqrt(NonNegative) → NonNegative, Real
            Tribool arg_nonneg = query_nonnegative(arg_expr);
            if (arg_nonneg != Tribool::True) return Tribool::Unknown;

            switch (target) {
                case Sign::NonNegative: return Tribool::True;
                case Sign::Negative:    return Tribool::False;
                case Sign::NonPositive: return Tribool::Unknown; // sqrt could be zero
                case Sign::Positive:    return Tribool::Unknown; // sqrt(0) = 0
                case Sign::Zero:        return Tribool::Unknown; // depends on argument
                case Sign::NonZero:     return Tribool::Unknown; // depends on argument
            }
            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::Tan: {
            // tan(Real) → Real
            // tan can be any real value — no sign determination possible
            Tribool arg_real = query_real(arg_expr);
            if (arg_real != Tribool::True) return Tribool::Unknown;

            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::ArcSin:
        case FunctionNode::FuncType::ArcCos:
        case FunctionNode::FuncType::ArcTan: {
            // Recognized but no specific sign rules defined in requirements
            return Tribool::Unknown;
        }

        default:
            // Unrecognized function → Unknown for all properties
            return Tribool::Unknown;
    }
}

Tribool InferenceEngine::infer_function_domain(const FunctionNode& node, Domain target) const {
    // Need at least one argument for all built-in function rules
    if (node.arguments.empty()) return Tribool::Unknown;

    SymbolicExpr arg_expr;
    arg_expr.root = node.arguments[0];

    switch (node.type) {
        case FunctionNode::FuncType::Exp: {
            // exp(Real) → Real
            if (target == Domain::Real) {
                Tribool arg_real = query_real(arg_expr);
                if (arg_real == Tribool::True) return Tribool::True;
                // Integer implies Real
                Tribool arg_int = query_integer(arg_expr);
                if (arg_int == Tribool::True) return Tribool::True;
            }
            // exp does not produce integers in general
            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::Sin: {
            // sin(Real) → Real
            if (target == Domain::Real) {
                Tribool arg_real = query_real(arg_expr);
                if (arg_real == Tribool::True) return Tribool::True;
                Tribool arg_int = query_integer(arg_expr);
                if (arg_int == Tribool::True) return Tribool::True;
            }
            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::Cos: {
            // cos(Real) → Real
            if (target == Domain::Real) {
                Tribool arg_real = query_real(arg_expr);
                if (arg_real == Tribool::True) return Tribool::True;
                Tribool arg_int = query_integer(arg_expr);
                if (arg_int == Tribool::True) return Tribool::True;
            }
            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::Abs: {
            // abs(Real) → Real
            if (target == Domain::Real) {
                Tribool arg_real = query_real(arg_expr);
                if (arg_real == Tribool::True) return Tribool::True;
                Tribool arg_int = query_integer(arg_expr);
                if (arg_int == Tribool::True) return Tribool::True;
            }
            // abs(Integer) → Integer (non-negative integer)
            if (target == Domain::Integer) {
                Tribool arg_int = query_integer(arg_expr);
                if (arg_int == Tribool::True) return Tribool::True;
            }
            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::Ln: {
            // ln(Positive) → Real
            if (target == Domain::Real) {
                Tribool arg_positive = query_positive(arg_expr);
                if (arg_positive == Tribool::True) return Tribool::True;
            }
            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::Sqrt: {
            // sqrt(NonNegative) → Real
            if (target == Domain::Real) {
                Tribool arg_nonneg = query_nonnegative(arg_expr);
                if (arg_nonneg == Tribool::True) return Tribool::True;
            }
            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::Tan: {
            // tan(Real) → Real
            if (target == Domain::Real) {
                Tribool arg_real = query_real(arg_expr);
                if (arg_real == Tribool::True) return Tribool::True;
                Tribool arg_int = query_integer(arg_expr);
                if (arg_int == Tribool::True) return Tribool::True;
            }
            return Tribool::Unknown;
        }

        case FunctionNode::FuncType::ArcSin:
        case FunctionNode::FuncType::ArcCos:
        case FunctionNode::FuncType::ArcTan: {
            // Recognized inverse trig functions — produce Real from appropriate input
            // but no specific requirements defined, return Unknown
            return Tribool::Unknown;
        }

        default:
            // Unrecognized function → Unknown for all properties
            return Tribool::Unknown;
    }
}

// ============================================================
// Interval propagation
// ============================================================

/**
 * @brief Helper: extract a numeric double value from an Endpoint.
 * Returns std::nullopt if the endpoint is not a finite numeric value.
 */
static std::optional<double> endpoint_to_double(const Endpoint& ep) {
    if (ep.is_neg_infinity || ep.is_pos_infinity) return std::nullopt;
    if (!ep.value || !ep.value->root) return std::nullopt;
    if (!ep.value->is_number()) return std::nullopt;
    double v = ep.value->to_numeric();
    if (!std::isfinite(v)) return std::nullopt;
    return v;
}

/**
 * @brief Helper: create a closed Endpoint from a double value.
 */
static Endpoint make_closed_endpoint(double val) {
    auto expr = std::make_shared<SymbolicExpr>(
        std::make_shared<NumberNode>(static_cast<lmmc_real_t>(val)));
    return Endpoint::closed(expr);
}

/**
 * @brief Helper: create an Interval from two double values [lo, hi].
 */
static Interval make_interval(double lo, double hi) {
    Interval result;
    result.lower = make_closed_endpoint(lo);
    result.upper = make_closed_endpoint(hi);
    return result;
}

/**
 * @brief Helper: check if an exponent NumberNode is exactly 2 (for squaring detection).
 */
static bool is_exponent_two(const std::shared_ptr<SymbolicNode>& node) {
    auto num = std::dynamic_pointer_cast<NumberNode>(node);
    if (!num) return false;
    if (std::holds_alternative<BigInt>(num->value)) {
        return std::get<BigInt>(num->value) == BigInt(2);
    }
    if (std::holds_alternative<Rational>(num->value)) {
        return std::get<Rational>(num->value) == Rational(2);
    }
    if (std::holds_alternative<lmmc_real_t>(num->value)) {
        lmmc_real_t v = std::get<lmmc_real_t>(num->value);
        return std::isfinite(v) && v == 2.0;
    }
    return false;
}

/**
 * @brief Helper: check if an exponent NumberNode is exactly -1 (for division detection).
 */
static bool is_exponent_neg_one(const std::shared_ptr<SymbolicNode>& node) {
    auto num = std::dynamic_pointer_cast<NumberNode>(node);
    if (!num) return false;
    if (std::holds_alternative<BigInt>(num->value)) {
        return std::get<BigInt>(num->value) == BigInt(-1);
    }
    if (std::holds_alternative<Rational>(num->value)) {
        return std::get<Rational>(num->value) == Rational(-1);
    }
    if (std::holds_alternative<lmmc_real_t>(num->value)) {
        lmmc_real_t v = std::get<lmmc_real_t>(num->value);
        return std::isfinite(v) && v == -1.0;
    }
    return false;
}

std::optional<Interval> InferenceEngine::propagate_bounds(const SymbolicExpr& expr) const {
    if (!expr.root) return std::nullopt;

    // --- NumberNode: point interval [n, n] ---
    if (auto num = std::dynamic_pointer_cast<NumberNode>(expr.root)) {
        double val = 0.0;
        if (std::holds_alternative<BigInt>(num->value)) {
            val = std::get<BigInt>(num->value).to_double();
        } else if (std::holds_alternative<Rational>(num->value)) {
            val = std::get<Rational>(num->value).to_double();
        } else {
            val = std::get<lmmc_real_t>(num->value);
        }
        if (!std::isfinite(val)) return std::nullopt;
        return make_interval(val, val);
    }

    // --- VariableNode: look up bounds from AssumptionContext ---
    if (auto var = std::dynamic_pointer_cast<VariableNode>(expr.root)) {
        auto bounds = ctx_.get_bounds(var->name);
        if (!bounds.has_value()) return std::nullopt;

        // Verify the bounds have numeric endpoints
        auto lo = endpoint_to_double(bounds->lower);
        auto hi = endpoint_to_double(bounds->upper);
        if (!lo.has_value() || !hi.has_value()) return std::nullopt;

        return bounds.value();
    }

    // --- AddNode: [a,b] + [c,d] → [a+c, b+d] for each operand ---
    if (auto add = std::dynamic_pointer_cast<AddNode>(expr.root)) {
        if (add->operands.empty()) return std::nullopt;

        // Start with the bounds of the first operand
        SymbolicExpr first_expr;
        first_expr.root = add->operands[0];
        auto result = propagate_bounds(first_expr);
        if (!result.has_value()) return std::nullopt;

        // Accumulate bounds for remaining operands
        for (size_t i = 1; i < add->operands.size(); ++i) {
            SymbolicExpr op_expr;
            op_expr.root = add->operands[i];
            auto op_bounds = propagate_bounds(op_expr);
            if (!op_bounds.has_value()) return std::nullopt;

            auto r_lo = endpoint_to_double(result->lower);
            auto r_hi = endpoint_to_double(result->upper);
            auto o_lo = endpoint_to_double(op_bounds->lower);
            auto o_hi = endpoint_to_double(op_bounds->upper);

            if (!r_lo || !r_hi || !o_lo || !o_hi) return std::nullopt;

            // Addition: [a,b] + [c,d] = [a+c, b+d]
            result = make_interval(*r_lo + *o_lo, *r_hi + *o_hi);
        }

        return result;
    }

    // --- MultiplyNode: handle multiplication and division ---
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr.root)) {
        if (mul->operands.empty()) return std::nullopt;

        // Separate operands into numerator terms and denominator terms (PowerNode with exp -1)
        std::vector<std::shared_ptr<SymbolicNode>> numerator_ops;
        std::vector<std::shared_ptr<SymbolicNode>> denominator_ops;

        for (const auto& operand : mul->operands) {
            auto pow_node = std::dynamic_pointer_cast<PowerNode>(operand);
            if (pow_node && is_exponent_neg_one(pow_node->exponent)) {
                // This is a division: base^(-1) means we divide by base
                denominator_ops.push_back(pow_node->base);
            } else {
                numerator_ops.push_back(operand);
            }
        }

        // Compute numerator product bounds
        std::optional<Interval> num_result;
        if (numerator_ops.empty()) {
            // All operands are denominators; numerator is implicitly 1
            num_result = make_interval(1.0, 1.0);
        } else {
            SymbolicExpr first_expr;
            first_expr.root = numerator_ops[0];
            num_result = propagate_bounds(first_expr);
            if (!num_result.has_value()) return std::nullopt;

            for (size_t i = 1; i < numerator_ops.size(); ++i) {
                SymbolicExpr op_expr;
                op_expr.root = numerator_ops[i];
                auto op_bounds = propagate_bounds(op_expr);
                if (!op_bounds.has_value()) return std::nullopt;

                auto a = endpoint_to_double(num_result->lower);
                auto b = endpoint_to_double(num_result->upper);
                auto c = endpoint_to_double(op_bounds->lower);
                auto d = endpoint_to_double(op_bounds->upper);

                if (!a || !b || !c || !d) return std::nullopt;

                // Multiplication: [a,b] * [c,d] = [min(ac,ad,bc,bd), max(ac,ad,bc,bd)]
                double ac = (*a) * (*c);
                double ad = (*a) * (*d);
                double bc = (*b) * (*c);
                double bd = (*b) * (*d);

                double lo = std::min({ac, ad, bc, bd});
                double hi = std::max({ac, ad, bc, bd});
                num_result = make_interval(lo, hi);
            }
        }

        // If no denominator terms, return the numerator result
        if (denominator_ops.empty()) {
            return num_result;
        }

        // Compute denominator product bounds, then divide
        // First compute the combined denominator interval
        SymbolicExpr den_first_expr;
        den_first_expr.root = denominator_ops[0];
        auto den_result = propagate_bounds(den_first_expr);
        if (!den_result.has_value()) return std::nullopt;

        for (size_t i = 1; i < denominator_ops.size(); ++i) {
            SymbolicExpr op_expr;
            op_expr.root = denominator_ops[i];
            auto op_bounds = propagate_bounds(op_expr);
            if (!op_bounds.has_value()) return std::nullopt;

            auto a = endpoint_to_double(den_result->lower);
            auto b = endpoint_to_double(den_result->upper);
            auto c = endpoint_to_double(op_bounds->lower);
            auto d = endpoint_to_double(op_bounds->upper);

            if (!a || !b || !c || !d) return std::nullopt;

            double ac = (*a) * (*c);
            double ad = (*a) * (*d);
            double bc = (*b) * (*c);
            double bd = (*b) * (*d);

            double lo = std::min({ac, ad, bc, bd});
            double hi = std::max({ac, ad, bc, bd});
            den_result = make_interval(lo, hi);
        }

        // Now divide numerator by denominator
        auto num_lo = endpoint_to_double(num_result->lower);
        auto num_hi = endpoint_to_double(num_result->upper);
        auto den_lo = endpoint_to_double(den_result->lower);
        auto den_hi = endpoint_to_double(den_result->upper);

        if (!num_lo || !num_hi || !den_lo || !den_hi) return std::nullopt;

        // Division where divisor interval contains zero → unbounded (no finite interval)
        if (*den_lo <= 0.0 && *den_hi >= 0.0) {
            return std::nullopt;
        }

        // Division by positive interval [c,d] (c>0, a>=0): [a,b]/[c,d] → [a/d, b/c]
        // General division: [a,b]/[c,d] where interval doesn't contain zero
        // Use the general formula: compute all combinations and take min/max
        double a = *num_lo, b = *num_hi, c = *den_lo, d = *den_hi;
        double ac_d = a / c;
        double ad_d = a / d;
        double bc_d = b / c;
        double bd_d = b / d;

        double lo = std::min({ac_d, ad_d, bc_d, bd_d});
        double hi = std::max({ac_d, ad_d, bc_d, bd_d});
        return make_interval(lo, hi);
    }

    // --- PowerNode: handle squaring ---
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(expr.root)) {
        // Squaring: x^2 where x >= 0 → [a², b²]
        if (is_exponent_two(pow->exponent)) {
            SymbolicExpr base_expr;
            base_expr.root = pow->base;
            auto base_bounds = propagate_bounds(base_expr);
            if (!base_bounds.has_value()) return std::nullopt;

            auto a = endpoint_to_double(base_bounds->lower);
            auto b = endpoint_to_double(base_bounds->upper);
            if (!a || !b) return std::nullopt;

            // If a >= 0 (non-negative base): [a², b²]
            if (*a >= 0.0) {
                return make_interval((*a) * (*a), (*b) * (*b));
            }

            // If b <= 0 (non-positive base): [b², a²]
            if (*b <= 0.0) {
                return make_interval((*b) * (*b), (*a) * (*a));
            }

            // If interval spans zero: [0, max(a², b²)]
            double a_sq = (*a) * (*a);
            double b_sq = (*b) * (*b);
            return make_interval(0.0, std::max(a_sq, b_sq));
        }

        // For other exponents, we don't propagate bounds in this implementation
        return std::nullopt;
    }

    // --- FunctionNode: sin/cos → [-1, 1] ---
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(expr.root)) {
        if (func->arguments.empty()) return std::nullopt;

        SymbolicExpr arg_expr;
        arg_expr.root = func->arguments[0];

        switch (func->type) {
            case FunctionNode::FuncType::Sin:
            case FunctionNode::FuncType::Cos: {
                // sin/cos of bounded input → [-1, 1]
                auto arg_bounds = propagate_bounds(arg_expr);
                if (!arg_bounds.has_value()) return std::nullopt;
                // If the argument has finite bounds, sin/cos is bounded by [-1, 1]
                return make_interval(-1.0, 1.0);
            }
            default:
                return std::nullopt;
        }
    }

    return std::nullopt;
}

// ============================================================
// Monotonicity deduction
// ============================================================

/**
 * @brief Helper: collect positive integer exponents that appear in PowerNode
 * expressions within the RelationStore involving the given variable names.
 */
static std::vector<int> collect_power_exponents(const RelationStore& store,
                                                const std::string& lhs_name,
                                                const std::string& rhs_name) {
    std::unordered_set<int> exponents;

    auto check_node = [&](const std::shared_ptr<SymbolicNode>& node) {
        if (!node) return;
        auto pow_node = std::dynamic_pointer_cast<PowerNode>(node);
        if (!pow_node) return;

        // Check if the base is one of our variables
        auto base_var = std::dynamic_pointer_cast<VariableNode>(pow_node->base);
        if (!base_var) return;
        if (base_var->name != lhs_name && base_var->name != rhs_name) return;

        // Check if the exponent is a positive integer
        auto exp_num = std::dynamic_pointer_cast<NumberNode>(pow_node->exponent);
        if (!exp_num) return;
        if (!is_positive_integer_number(*exp_num)) return;

        // Extract the integer value
        int n = 0;
        if (std::holds_alternative<BigInt>(exp_num->value)) {
            const auto& b = std::get<BigInt>(exp_num->value);
            // Only consider small exponents to avoid explosion
            if (b > BigInt(100)) return;
            n = static_cast<int>(b.IsNegative() ? 0 : b.is_zero() ? 0 : std::stoi(b.to_string()));
        } else if (std::holds_alternative<Rational>(exp_num->value)) {
            BigInt num_val = std::get<Rational>(exp_num->value).get_numerator();
            if (num_val > BigInt(100)) return;
            n = static_cast<int>(std::stoi(num_val.to_string()));
        } else {
            double v = std::get<lmmc_real_t>(exp_num->value);
            if (v <= 0.0 || v > 100.0 || v != std::floor(v)) return;
            n = static_cast<int>(v);
        }

        if (n > 0) {
            exponents.insert(n);
        }
    };

    // Scan all relations in the store for PowerNode expressions
    for (const auto& rel : store.get_relations()) {
        check_node(rel.lhs.root);
        check_node(rel.rhs.root);
    }

    return std::vector<int>(exponents.begin(), exponents.end());
}

void InferenceEngine::apply_monotonicity_rules(const Relation& rel, RelationStore& store,
                                               PropertyStore& prop_store, int depth) {
    // Stop recursion at maximum depth
    if (depth >= MAX_MONOTONICITY_DEPTH) return;

    // Only apply to GT (greater-than) relations
    if (rel.op != RelationalNode::Op::GT) return;

    // Extract LHS and RHS — both must be single VariableNodes
    auto lhs_var = std::dynamic_pointer_cast<VariableNode>(rel.lhs.root);
    auto rhs_var = std::dynamic_pointer_cast<VariableNode>(rel.rhs.root);
    if (!lhs_var || !rhs_var) return;

    const std::string& x_name = lhs_var->name;
    const std::string& y_name = rhs_var->name;

    // Helper to add a deduced relation and recursively apply monotonicity rules
    auto add_deduced = [&](const SymbolicExpr& new_lhs, const SymbolicExpr& new_rhs) {
        // Don't add if already present
        if (store.has_relation(new_lhs, new_rhs, RelationalNode::Op::GT)) return;

        store.add_relation(new_lhs, new_rhs, RelationalNode::Op::GT, prop_store);

        // Recursively apply monotonicity rules to the newly added relation
        Relation new_rel{new_lhs, new_rhs, RelationalNode::Op::GT};
        apply_monotonicity_rules(new_rel, store, prop_store, depth + 1);
    };

    // Check domain assumptions using the AssumptionContext (read-through all scopes)
    bool both_positive = ctx_.has_sign(x_name, Sign::Positive) &&
                         ctx_.has_sign(y_name, Sign::Positive);
    bool both_real = ctx_.has_domain(x_name, Domain::Real) &&
                     ctx_.has_domain(y_name, Domain::Real);
    bool both_nonnegative = ctx_.has_sign(x_name, Sign::NonNegative) &&
                            ctx_.has_sign(y_name, Sign::NonNegative);

    // Rule 1 (Req 15.1): ln rule — x > y, both Positive → ln(x) > ln(y)
    if (both_positive) {
        auto ln_x_node = std::make_shared<FunctionNode>(
            FunctionNode::FuncType::Ln,
            std::vector<std::shared_ptr<SymbolicNode>>{lhs_var->clone()});
        auto ln_y_node = std::make_shared<FunctionNode>(
            FunctionNode::FuncType::Ln,
            std::vector<std::shared_ptr<SymbolicNode>>{rhs_var->clone()});

        SymbolicExpr ln_x;
        ln_x.root = ln_x_node;
        SymbolicExpr ln_y;
        ln_y.root = ln_y_node;

        add_deduced(ln_x, ln_y);
    }

    // Rule 2 (Req 15.2): sqrt rule — x > y, both Positive → sqrt(x) > sqrt(y)
    if (both_positive) {
        auto sqrt_x_node = std::make_shared<FunctionNode>(
            FunctionNode::FuncType::Sqrt,
            std::vector<std::shared_ptr<SymbolicNode>>{lhs_var->clone()});
        auto sqrt_y_node = std::make_shared<FunctionNode>(
            FunctionNode::FuncType::Sqrt,
            std::vector<std::shared_ptr<SymbolicNode>>{rhs_var->clone()});

        SymbolicExpr sqrt_x;
        sqrt_x.root = sqrt_x_node;
        SymbolicExpr sqrt_y;
        sqrt_y.root = sqrt_y_node;

        add_deduced(sqrt_x, sqrt_y);
    }

    // Rule 3 (Req 15.3): exp rule — x > y, both Real → exp(x) > exp(y)
    if (both_real) {
        auto exp_x_node = std::make_shared<FunctionNode>(
            FunctionNode::FuncType::Exp,
            std::vector<std::shared_ptr<SymbolicNode>>{lhs_var->clone()});
        auto exp_y_node = std::make_shared<FunctionNode>(
            FunctionNode::FuncType::Exp,
            std::vector<std::shared_ptr<SymbolicNode>>{rhs_var->clone()});

        SymbolicExpr exp_x;
        exp_x.root = exp_x_node;
        SymbolicExpr exp_y;
        exp_y.root = exp_y_node;

        add_deduced(exp_x, exp_y);
    }

    // Rule 4 (Req 15.4): power rule — x > y, both NonNegative → x^n > y^n
    // for each positive integer n appearing in expressions involving x or y
    if (both_nonnegative) {
        std::vector<int> exponents = collect_power_exponents(store, x_name, y_name);

        for (int n : exponents) {
            auto exp_node = std::make_shared<NumberNode>(BigInt(n));
            auto pow_x_node = std::make_shared<PowerNode>(lhs_var->clone(), exp_node->clone());
            auto pow_y_node = std::make_shared<PowerNode>(rhs_var->clone(), exp_node->clone());

            SymbolicExpr pow_x;
            pow_x.root = pow_x_node;
            SymbolicExpr pow_y;
            pow_y.root = pow_y_node;

            add_deduced(pow_x, pow_y);
        }
    }
}

} // namespace lamina
