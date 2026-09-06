#define _USE_MATH_DEFINES
#include "internal/inference_engine_impl.hpp"

namespace LMCAS {

// Interval propagation

/**
 * @brief Helper: extract a numeric double value from an Endpoint.
 * Returns std::nullopt if the endpoint is not a finite numeric value.
 */
std::optional<double> endpoint_to_double(const Endpoint& ep) {
    if (ep.is_neg_infinity || ep.is_pos_infinity) return std::nullopt;
    if (!ep.value || !LMCAS::detail::node(ep.value)) return std::nullopt;
    ComputationContext context;
    auto evaluated = evaluate_numeric(*ep.value, NumericBindings{}, context);
    if (!evaluated || !evaluated.value().is_finite()) return std::nullopt;
    double v = evaluated.value().value;
    if (!std::isfinite(v)) return std::nullopt;
    return v;
}

/**
 * @brief Helper: create a closed Endpoint from a double value.
 */
static Endpoint make_closed_endpoint(double val) {
    auto expr = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(val)));
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
static bool is_exponent_two(const std::shared_ptr<const SymbolicNode>& node) {
    auto num = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!num) return false;
    if (std::holds_alternative<BigInt>(num->value())) {
        return std::get<BigInt>(num->value()) == BigInt(2);
    }
    if (std::holds_alternative<Rational>(num->value())) {
        return std::get<Rational>(num->value()) == Rational(2);
    }
    if (std::holds_alternative<lmmc_real_t>(num->value())) {
        lmmc_real_t v = std::get<lmmc_real_t>(num->value());
        return std::isfinite(v) && v == 2.0;
    }
    return false;
}

/**
 * @brief Helper: check if an exponent NumberNode is exactly -1 (for division detection).
 */
bool is_exponent_neg_one(const std::shared_ptr<const SymbolicNode>& node) {
    auto num = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!num) return false;
    if (std::holds_alternative<BigInt>(num->value())) {
        return std::get<BigInt>(num->value()) == BigInt(-1);
    }
    if (std::holds_alternative<Rational>(num->value())) {
        return std::get<Rational>(num->value()) == Rational(-1);
    }
    if (std::holds_alternative<lmmc_real_t>(num->value())) {
        lmmc_real_t v = std::get<lmmc_real_t>(num->value());
        return std::isfinite(v) && v == -1.0;
    }
    return false;
}

std::optional<Interval> InferenceEngine::propagate_bounds(const SymbolicExpr& expr) const {
    if (!LMCAS::detail::node(expr)) return std::nullopt;

    // --- NumberNode: point interval [n, n] ---
    if (auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr))) {
        double val = 0.0;
        if (std::holds_alternative<BigInt>(num->value())) {
            val = std::get<BigInt>(num->value()).to_double();
        } else if (std::holds_alternative<Rational>(num->value())) {
            val = std::get<Rational>(num->value()).to_double();
        } else {
            val = std::get<lmmc_real_t>(num->value());
        }
        if (!std::isfinite(val)) return std::nullopt;
        return make_interval(val, val);
    }

    // --- VariableNode: look up bounds from AssumptionContext ---
    if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
        auto bounds = impl_->ctx.get_bounds(var->name());
        if (!bounds.has_value()) return std::nullopt;

        // Verify the bounds have numeric endpoints
        auto lo = endpoint_to_double(bounds->lower);
        auto hi = endpoint_to_double(bounds->upper);
        if (!lo.has_value() || !hi.has_value()) return std::nullopt;

        return bounds.value();
    }

    // --- AddNode: [a,b] + [c,d] → [a+c, b+d] for each operand ---
    if (auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr))) {
        if (add->operands().empty()) return std::nullopt;

        // Start with the bounds of the first operand
        auto first_expr = LMCAS::detail::expression_from_node(add->operands()[0]);
        auto result = propagate_bounds(first_expr);
        if (!result.has_value()) return std::nullopt;

        // Accumulate bounds for remaining operands
        for (size_t i = 1; i < add->operands().size(); ++i) {
            auto op_expr = LMCAS::detail::expression_from_node(add->operands()[i]);
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
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        if (mul->operands().empty()) return std::nullopt;

        // Separate operands into numerator terms and denominator terms (PowerNode with exp -1)
        std::vector<std::shared_ptr<const SymbolicNode>> numerator_ops;
        std::vector<std::shared_ptr<const SymbolicNode>> denominator_ops;

        for (const auto& operand : mul->operands()) {
            auto pow_node = std::dynamic_pointer_cast<const PowerNode>(operand);
            if (pow_node && is_exponent_neg_one(pow_node->exponent())) {
                // This is a division: base^(-1) means we divide by base
                denominator_ops.push_back(pow_node->base());
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
            auto first_expr = LMCAS::detail::expression_from_node(numerator_ops[0]);
            num_result = propagate_bounds(first_expr);
            if (!num_result.has_value()) return std::nullopt;

            for (size_t i = 1; i < numerator_ops.size(); ++i) {
                auto op_expr = LMCAS::detail::expression_from_node(numerator_ops[i]);
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
        auto den_first_expr = LMCAS::detail::expression_from_node(denominator_ops[0]);
        auto den_result = propagate_bounds(den_first_expr);
        if (!den_result.has_value()) return std::nullopt;

        for (size_t i = 1; i < denominator_ops.size(); ++i) {
            auto op_expr = LMCAS::detail::expression_from_node(denominator_ops[i]);
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
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(expr))) {
        // Squaring: x^2 where x >= 0 → [a², b²]
        if (is_exponent_two(pow->exponent())) {
            auto base_expr = LMCAS::detail::expression_from_node(pow->base());
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
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(expr))) {
        if (func->arguments().empty()) return std::nullopt;

        auto arg_expr = LMCAS::detail::expression_from_node(func->arguments()[0]);
        switch (func->type()) {
            case FunctionNode::FuncType::Sin:
            case FunctionNode::FuncType::Cos: {
                // sin/cos of bounded input → [-1, 1]
                auto arg_bounds = propagate_bounds(arg_expr);
                if (!arg_bounds.has_value()) return std::nullopt;
                // If the argument has finite bounds, sin/cos is bounded by [-1, 1]
                return make_interval(-1.0, 1.0);
            }
            case FunctionNode::FuncType::ArcTan: {
                auto argument_bounds = propagate_bounds(arg_expr);
                if (argument_bounds) {
                    Interval result;
                    if (argument_bounds->lower.is_neg_infinity) {
                        result.lower = Endpoint::open(
                            SymbolicExpr::number(-LMMC_CONST_PI / 2.0));
                    } else if (auto lower =
                                   endpoint_to_double(
                                       argument_bounds->lower)) {
                        auto value = SymbolicExpr::number(std::atan(*lower));
                        result.lower = argument_bounds->lower.is_open
                            ? Endpoint::open(value)
                            : Endpoint::closed(value);
                    } else {
                        return std::nullopt;
                    }

                    if (argument_bounds->upper.is_pos_infinity) {
                        result.upper = Endpoint::open(
                            SymbolicExpr::number(LMMC_CONST_PI / 2.0));
                    } else if (auto upper =
                                   endpoint_to_double(
                                       argument_bounds->upper)) {
                        auto value = SymbolicExpr::number(std::atan(*upper));
                        result.upper = argument_bounds->upper.is_open
                            ? Endpoint::open(value)
                            : Endpoint::closed(value);
                    } else {
                        return std::nullopt;
                    }
                    return result;
                }

                auto real = query_real_checked(arg_expr);
                if (real && real.value() == Tribool::True) {
                    return Interval{
                        Endpoint::open(SymbolicExpr::number(
                            -LMMC_CONST_PI / 2.0)),
                        Endpoint::open(SymbolicExpr::number(
                            LMMC_CONST_PI / 2.0))};
                }
                return std::nullopt;
            }
            default:
                return std::nullopt;
        }
    }

    return std::nullopt;
}


} // namespace LMCAS
