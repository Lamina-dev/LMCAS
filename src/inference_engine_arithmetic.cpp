#define _USE_MATH_DEFINES
#include "internal/inference_engine_impl.hpp"

namespace LMCAS {

#define LMCAS_TRY_TRIBOOL(NAME, EXPRESSION) \
    auto NAME##_result = (EXPRESSION);       \
    if (!NAME##_result) return NAME##_result; \
    const Tribool NAME = NAME##_result.value()


InferenceTriboolResult InferenceEngine::infer_add_sign_checked(
    const AddNode& node, Sign target) const {
    try {
        if (node.operands().empty()) return InferenceTriboolResult::success(Tribool::Unknown);

        // First, try subtraction pattern detection for 2-operand AddNodes
        if (node.operands().size() == 2) {
            // Check if one operand is a negated term (MultiplyNode with -1 coefficient)
            bool has_negated = false;
            for (const auto& operand : node.operands()) {
                auto mul = std::dynamic_pointer_cast<const MultiplyNode>(operand);
                if (mul && mul->operands().size() == 2) {
                    for (const auto& mul_op : mul->operands()) {
                        auto num = std::dynamic_pointer_cast<const NumberNode>(mul_op);
                        if (num) {
                            bool is_neg_one = false;
                            if (std::holds_alternative<BigInt>(num->value())) {
                                is_neg_one = (std::get<BigInt>(num->value()) == BigInt(-1));
                            } else if (std::holds_alternative<Rational>(num->value())) {
                                is_neg_one = (std::get<Rational>(num->value()) == Rational(-1));
                            } else if (std::holds_alternative<lmmc_real_t>(num->value())) {
                                lmmc_real_t v = std::get<lmmc_real_t>(num->value());
                                is_neg_one = (std::isfinite(v) && v == -1.0);
                            }
                            if (is_neg_one) {
                                has_negated = true;
                                break;
                            }
                        }
                    }
                    if (has_negated) break;
                }
            }
            if (has_negated) {
                auto sub_result =
                    infer_subtraction_sign_checked(node, target);
                if (!sub_result) return sub_result;
                if (sub_result.value() != Tribool::Unknown) {
                    return sub_result;
                }
            }
        }

        // For sign inference on addition:
        // - If all operands have the target sign property -> sum has that property
        // - If any operand has Unknown -> result is Unknown
        // - If operands have mixed definite signs -> result is Unknown

        bool all_have_property = true;

        for (const auto& operand : node.operands()) {
            auto op_expr = LMCAS::detail::expression_from_node(operand);
            Tribool op_result = Tribool::Unknown;
            switch (target) {
                case Sign::Positive:
                case Sign::Negative:
                case Sign::NonNegative:
                case Sign::NonPositive:
                    {
                        auto queried = query_sign_of_checked(op_expr, target);
                        if (!queried) return queried;
                        op_result = queried.value();
                    }
                    break;
                default:
                    return InferenceTriboolResult::success(Tribool::Unknown);
            }

            if (op_result == Tribool::Unknown) {
                return InferenceTriboolResult::success(Tribool::Unknown);
            }
            if (op_result == Tribool::False) {
                all_have_property = false;
            }
        }

        if (all_have_property) {
            return InferenceTriboolResult::success(Tribool::True);
        }

        return InferenceTriboolResult::success(Tribool::Unknown);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "addition sign inference allocation failed",
            "infer_add_sign_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_add_sign_checked");
    }
}

// Sign inference from relational constraints


InferenceTriboolResult InferenceEngine::infer_sign_from_relations_checked(
    const SymbolicExpr& expr, Sign target) const {
    if (!LMCAS::detail::node(expr)) {
        return InferenceTriboolResult::failure(
            CasErrc::InvalidArgument,
            "relation sign inference expression must not be null",
            "infer_sign_from_relations_checked");
    }

    try {
        const auto& rel_store = impl_->ctx.current_relations();
        auto zero_expr = LMCAS::detail::expression_from_node(
            LMCAS::detail::make_node<NumberNode>(BigInt(0)));

        if (auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr))) {
            if (add->operands().empty()) return InferenceTriboolResult::success(Tribool::Unknown);

            bool all_gt_zero = true;
            bool all_geq_zero = true;

            for (const auto& operand : add->operands()) {
                auto op_expr = LMCAS::detail::expression_from_node(operand);
                bool op_gt = rel_store.has_relation(op_expr, zero_expr, RelationalNode::Op::GT);
                bool op_geq = rel_store.has_relation(op_expr, zero_expr, RelationalNode::Op::GEQ);

                if (!op_gt) {
                    if (auto var = std::dynamic_pointer_cast<const VariableNode>(operand)) {
                        op_gt = impl_->ctx.current_properties().has_sign(var->name(), Sign::Positive);
                    } else if (auto num = std::dynamic_pointer_cast<const NumberNode>(operand)) {
                        op_gt = num->is_positive();
                    }
                }
                if (!op_geq && !op_gt) {
                    if (auto var = std::dynamic_pointer_cast<const VariableNode>(operand)) {
                        op_geq = impl_->ctx.current_properties().has_sign(var->name(), Sign::NonNegative);
                    } else if (auto num = std::dynamic_pointer_cast<const NumberNode>(operand)) {
                        op_geq = num->is_zero() || num->is_positive();
                    }
                }

                if (!op_gt) all_gt_zero = false;
                if (!op_gt && !op_geq) all_geq_zero = false;
            }

            if (target == Sign::Positive && all_gt_zero) {
                return InferenceTriboolResult::success(Tribool::True);
            }
            if (target == Sign::NonNegative && all_geq_zero) {
                return InferenceTriboolResult::success(Tribool::True);
            }
        }

        if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
            if (mul->operands().empty()) return InferenceTriboolResult::success(Tribool::Unknown);

            if (target == Sign::Positive) {
                bool all_gt_zero = true;

                for (const auto& operand : mul->operands()) {
                    auto op_expr = LMCAS::detail::expression_from_node(operand);
                    bool op_gt = rel_store.has_relation(op_expr, zero_expr, RelationalNode::Op::GT);

                    if (!op_gt) {
                        if (auto var = std::dynamic_pointer_cast<const VariableNode>(operand)) {
                            op_gt = impl_->ctx.current_properties().has_sign(var->name(), Sign::Positive);
                        } else if (auto num = std::dynamic_pointer_cast<const NumberNode>(operand)) {
                            op_gt = num->is_positive();
                        }
                    }

                    if (!op_gt) {
                        all_gt_zero = false;
                        break;
                    }
                }

                if (all_gt_zero) return InferenceTriboolResult::success(Tribool::True);
            }
        }

        if (target == Sign::Positive) {
            const auto& relations = rel_store.get_relations();
            for (const auto& rel : relations) {
                if (rel.op != RelationalNode::Op::GT) continue;

                if (!LMCAS::detail::node(rel.lhs) || !LMCAS::detail::node(expr)) continue;
                if (!LMCAS::detail::node(rel.lhs)->equals(*LMCAS::detail::node(expr))) continue;
                if (!LMCAS::detail::node(rel.rhs)) continue;

                if (auto rhs_num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(rel.rhs))) {
                    if (rhs_num->is_zero() || rhs_num->is_positive()) {
                        return InferenceTriboolResult::success(Tribool::True);
                    }
                }

                if (auto rhs_var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(rel.rhs))) {
                    if (impl_->ctx.current_properties().has_sign(rhs_var->name(), Sign::NonNegative)) {
                        return InferenceTriboolResult::success(Tribool::True);
                    }
                }

                if (rel_store.has_relation(rel.rhs, zero_expr, RelationalNode::Op::GEQ) ||
                    rel_store.has_relation(rel.rhs, zero_expr, RelationalNode::Op::GT)) {
                    return InferenceTriboolResult::success(Tribool::True);
                }
            }
        }

        return InferenceTriboolResult::success(Tribool::Unknown);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "relation sign inference allocation failed",
            "infer_sign_from_relations_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_sign_from_relations_checked");
    }
}

// Addition domain inference


InferenceTriboolResult InferenceEngine::infer_add_domain_checked(
    const AddNode& node, Domain target) const {
    try {
        if (node.operands().empty()) return InferenceTriboolResult::success(Tribool::Unknown);

        // Addition is closed over Integer and Real. Integer operands also prove Real.
        for (const auto& operand : node.operands()) {
            auto op_expr = LMCAS::detail::expression_from_node(operand);
            Tribool op_result = Tribool::Unknown;
            switch (target) {
                case Domain::Integer: {
                    auto queried = query_integer_checked(op_expr);
                    if (!queried) return queried;
                    op_result = queried.value();
                    break;
                }
                case Domain::Real: {
                    auto queried = query_real_checked(op_expr);
                    if (!queried) return queried;
                    op_result = queried.value();
                    break;
                }
                default:
                    return InferenceTriboolResult::success(Tribool::Unknown);
            }

            if (op_result != Tribool::True) {
                if (target == Domain::Real) {
                    auto integer = query_integer_checked(op_expr);
                    if (!integer) return integer;
                    if (integer.value() == Tribool::True) continue;
                }
                return InferenceTriboolResult::success(Tribool::Unknown);
            }
        }

        return InferenceTriboolResult::success(Tribool::True);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "addition domain inference allocation failed",
            "infer_add_domain_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_add_domain_checked");
    }
}

// Division sign inference


InferenceTriboolResult InferenceEngine::infer_division_sign_checked(
    const MultiplyNode& node, Sign target) const {
    try {
        // Division pattern: MultiplyNode with exactly 2 children where one is PowerNode(den, -1)
        if (node.operands().size() != 2) return InferenceTriboolResult::success(Tribool::Unknown);

        // Find the PowerNode with exponent -1 (denominator) and the other operand (numerator)
        std::shared_ptr<const SymbolicNode> numerator_node;
        std::shared_ptr<const SymbolicNode> denominator_node;

        for (const auto& operand : node.operands()) {
            auto pow_node = std::dynamic_pointer_cast<const PowerNode>(operand);
            if (pow_node && is_exponent_neg_one(pow_node->exponent())) {
                denominator_node = pow_node->base();
            } else {
                numerator_node = operand;
            }
        }

        // Must have exactly one denominator and one numerator
        if (!denominator_node || !numerator_node) return InferenceTriboolResult::success(Tribool::Unknown);

        auto num_expr = LMCAS::detail::expression_from_node(numerator_node);
        auto den_expr = LMCAS::detail::expression_from_node(denominator_node);
        // Check if denominator is zero -> return Unknown for all sign queries
        auto den_nn_result =
            query_sign_of_checked(den_expr, Sign::NonNegative);
        if (!den_nn_result) return den_nn_result;
        const Tribool den_nn = den_nn_result.value();
        auto den_np_result =
            query_sign_of_checked(den_expr, Sign::NonPositive);
        if (!den_np_result) return den_np_result;
        const Tribool den_np = den_np_result.value();
        if (den_nn == Tribool::True && den_np == Tribool::True) {
            // Denominator is zero
            return InferenceTriboolResult::success(Tribool::Unknown);
        }

        auto num_pos_result =
            query_sign_of_checked(num_expr, Sign::Positive);
        if (!num_pos_result) return num_pos_result;
        const Tribool num_pos = num_pos_result.value();
        auto num_neg_result =
            query_sign_of_checked(num_expr, Sign::Negative);
        if (!num_neg_result) return num_neg_result;
        const Tribool num_neg = num_neg_result.value();

        auto den_pos_result =
            query_sign_of_checked(den_expr, Sign::Positive);
        if (!den_pos_result) return den_pos_result;
        const Tribool den_pos = den_pos_result.value();
        auto den_neg_result =
            query_sign_of_checked(den_expr, Sign::Negative);
        if (!den_neg_result) return den_neg_result;
        const Tribool den_neg = den_neg_result.value();

        // If denominator sign is unknown, return Unknown
        if (den_pos != Tribool::True && den_neg != Tribool::True) {
            return InferenceTriboolResult::success(Tribool::Unknown);
        }

        // If numerator sign is unknown, return Unknown
        if (num_pos != Tribool::True && num_neg != Tribool::True) {
            // Check if numerator is zero
            auto num_nn_result =
                query_sign_of_checked(num_expr, Sign::NonNegative);
            if (!num_nn_result) return num_nn_result;
            const Tribool num_nn = num_nn_result.value();
            auto num_np_result =
                query_sign_of_checked(num_expr, Sign::NonPositive);
            if (!num_np_result) return num_np_result;
            const Tribool num_np_check = num_np_result.value();
            if (num_nn == Tribool::True && num_np_check == Tribool::True) {
                // Numerator is zero -> result is zero
                switch (target) {
                    case Sign::Zero:        return InferenceTriboolResult::success(Tribool::True);
                    case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::True);
                    case Sign::NonPositive: return InferenceTriboolResult::success(Tribool::True);
                    case Sign::Positive:    return InferenceTriboolResult::success(Tribool::False);
                    case Sign::Negative:    return InferenceTriboolResult::success(Tribool::False);
                    case Sign::NonZero:     return InferenceTriboolResult::success(Tribool::False);
                }
            }
            return InferenceTriboolResult::success(Tribool::Unknown);
        }

        // Apply sign table for division
        // positive / positive -> positive
        // negative / negative -> positive
        // positive / negative -> negative
        // negative / positive -> negative
        bool result_positive = (num_pos == Tribool::True && den_pos == Tribool::True) ||
                               (num_neg == Tribool::True && den_neg == Tribool::True);
        bool result_negative = (num_pos == Tribool::True && den_neg == Tribool::True) ||
                               (num_neg == Tribool::True && den_pos == Tribool::True);

        if (result_positive) {
            switch (target) {
                case Sign::Positive:    return InferenceTriboolResult::success(Tribool::True);
                case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::True);
                case Sign::NonZero:     return InferenceTriboolResult::success(Tribool::True);
                case Sign::Negative:    return InferenceTriboolResult::success(Tribool::False);
                case Sign::NonPositive: return InferenceTriboolResult::success(Tribool::False);
                case Sign::Zero:        return InferenceTriboolResult::success(Tribool::False);
            }
        }

        if (result_negative) {
            switch (target) {
                case Sign::Negative:    return InferenceTriboolResult::success(Tribool::True);
                case Sign::NonPositive: return InferenceTriboolResult::success(Tribool::True);
                case Sign::NonZero:     return InferenceTriboolResult::success(Tribool::True);
                case Sign::Positive:    return InferenceTriboolResult::success(Tribool::False);
                case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::False);
                case Sign::Zero:        return InferenceTriboolResult::success(Tribool::False);
            }
        }

        return InferenceTriboolResult::success(Tribool::Unknown);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "division sign inference allocation failed",
            "infer_division_sign_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_division_sign_checked");
    }
}

// Multiplication sign inference


InferenceTriboolResult InferenceEngine::infer_multiply_sign_checked(
    const MultiplyNode& node, Sign target) const {
    try {
        if (node.operands().empty()) return InferenceTriboolResult::success(Tribool::Unknown);

        // Check for division pattern: exactly 2 operands with one being PowerNode(den, -1)
        if (node.operands().size() == 2) {
            bool has_power_neg_one = false;
            for (const auto& operand : node.operands()) {
                auto pow_node = std::dynamic_pointer_cast<const PowerNode>(operand);
                if (pow_node && is_exponent_neg_one(pow_node->exponent())) {
                    has_power_neg_one = true;
                    break;
                }
            }
            if (has_power_neg_one) {
                auto div_result =
                    infer_division_sign_checked(node, target);
                if (!div_result) return div_result;
                if (div_result.value() != Tribool::Unknown) {
                    return div_result;
                }
            }
        }

        // Step 1: Check if any operand is Zero -> product is Zero
        bool has_zero = false;
        for (const auto& operand : node.operands()) {
            auto op_expr = LMCAS::detail::expression_from_node(operand);
            // Check if operand is zero (both nonnegative and nonpositive)
            auto nonnegative =
                query_sign_of_checked(op_expr, Sign::NonNegative);
            if (!nonnegative) return nonnegative;
            auto nonpositive =
                query_sign_of_checked(op_expr, Sign::NonPositive);
            if (!nonpositive) return nonpositive;
            const Tribool nn = nonnegative.value();
            const Tribool np = nonpositive.value();
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
                case Sign::Zero:        return InferenceTriboolResult::success(Tribool::True);
                case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::True);
                case Sign::NonPositive: return InferenceTriboolResult::success(Tribool::True);
                case Sign::Positive:    return InferenceTriboolResult::success(Tribool::False);
                case Sign::Negative:    return InferenceTriboolResult::success(Tribool::False);
                case Sign::NonZero:     return InferenceTriboolResult::success(Tribool::False);
            }
        }

        // Step 2: Handle NonZero query -- all NonZero -> product is NonZero
        if (target == Sign::NonZero) {
            for (const auto& operand : node.operands()) {
                auto op_expr = LMCAS::detail::expression_from_node(operand);
                auto nonzero =
                    query_sign_of_checked(op_expr, Sign::NonZero);
                if (!nonzero) return nonzero;
                if (nonzero.value() != Tribool::True) {
                    auto positive =
                        query_sign_of_checked(op_expr, Sign::Positive);
                    if (!positive) return positive;
                    if (positive.value() == Tribool::True) continue;
                    auto negative =
                        query_sign_of_checked(op_expr, Sign::Negative);
                    if (!negative) return negative;
                    if (negative.value() == Tribool::True) continue;
                    return InferenceTriboolResult::success(Tribool::Unknown);
                }
            }
            return InferenceTriboolResult::success(Tribool::True);
        }

        // Step 3: Determine sign by counting negatives
        // We need each operand to have a definite sign (Positive or Negative)
        // or at least NonNegative/NonPositive for weaker results.
        int negative_count = 0;
        bool all_definite = true;       // all are Positive or Negative
        bool all_nonneg_or_nonpos = true; // all are at least NonNeg or NonPos
        bool has_unknown = false;

        for (const auto& operand : node.operands()) {
            auto op_expr = LMCAS::detail::expression_from_node(operand);
            auto positive =
                query_sign_of_checked(op_expr, Sign::Positive);
            if (!positive) return positive;
            const Tribool pos = positive.value();
            auto negative =
                query_sign_of_checked(op_expr, Sign::Negative);
            if (!negative) return negative;
            const Tribool neg = negative.value();

            if (pos == Tribool::True) {
                // Positive operand -- doesn't change sign
                continue;
            } else if (neg == Tribool::True) {
                // Negative operand -- flips sign
                negative_count++;
                continue;
            }

            // Not definitively Positive or Negative
            all_definite = false;

            // Check weaker properties
            auto nonnegative =
                query_sign_of_checked(op_expr, Sign::NonNegative);
            if (!nonnegative) return nonnegative;
            const Tribool nn = nonnegative.value();
            auto nonpositive =
                query_sign_of_checked(op_expr, Sign::NonPositive);
            if (!nonpositive) return nonpositive;
            const Tribool np = nonpositive.value();

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
            return InferenceTriboolResult::success(Tribool::Unknown);
        }

        // Determine result based on parity of negative count
        bool even_negatives = (negative_count % 2 == 0);

        if (all_definite) {
            // All operands are definitively Positive or Negative
            switch (target) {
                case Sign::Positive:
                    return InferenceTriboolResult::success(even_negatives ? Tribool::True : Tribool::False);
                case Sign::Negative:
                    return InferenceTriboolResult::success(even_negatives ? Tribool::False : Tribool::True);
                case Sign::NonNegative:
                    return InferenceTriboolResult::success(even_negatives ? Tribool::True : Tribool::False);
                case Sign::NonPositive:
                    return InferenceTriboolResult::success(even_negatives ? Tribool::False : Tribool::True);
                case Sign::NonZero:
                    return InferenceTriboolResult::success(Tribool::True); // product of nonzero values is nonzero
                default:
                    return InferenceTriboolResult::success(Tribool::Unknown);
            }
        }

        if (all_nonneg_or_nonpos) {
            switch (target) {
                case Sign::NonNegative:
                    return InferenceTriboolResult::success(even_negatives ? Tribool::True : Tribool::False);
                case Sign::NonPositive:
                    return InferenceTriboolResult::success(even_negatives ? Tribool::False : Tribool::True);
                case Sign::Positive:
                    return InferenceTriboolResult::success(Tribool::Unknown); // could be zero
                case Sign::Negative:
                    return InferenceTriboolResult::success(Tribool::Unknown); // could be zero
                default:
                    return InferenceTriboolResult::success(Tribool::Unknown);
            }
        }

        return InferenceTriboolResult::success(Tribool::Unknown);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "multiplication sign inference allocation failed",
            "infer_multiply_sign_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_multiply_sign_checked");
    }
}

// Multiplication domain inference


InferenceTriboolResult InferenceEngine::infer_multiply_domain_checked(
    const MultiplyNode& node, Domain target) const {
    try {
        if (node.operands().empty()) return InferenceTriboolResult::success(Tribool::Unknown);

        if (target == Domain::Integer) {
            for (const auto& operand : node.operands()) {
                auto op_expr = LMCAS::detail::expression_from_node(operand);
                auto queried = query_integer_checked(op_expr);
                if (!queried) return queried;
                if (queried.value() != Tribool::True) {
                    return InferenceTriboolResult::success(Tribool::Unknown);
                }
            }
            return InferenceTriboolResult::success(Tribool::True);
        }

        if (target == Domain::Real) {
            for (const auto& operand : node.operands()) {
                auto op_expr = LMCAS::detail::expression_from_node(operand);
                auto real = query_real_checked(op_expr);
                if (!real) return real;
                if (real.value() == Tribool::True) continue;

                auto integer = query_integer_checked(op_expr);
                if (!integer) return integer;
                if (integer.value() == Tribool::True) continue;

                return InferenceTriboolResult::success(Tribool::Unknown);
            }
            return InferenceTriboolResult::success(Tribool::True);
        }

        return InferenceTriboolResult::success(Tribool::Unknown);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "multiplication domain inference allocation failed",
            "infer_multiply_domain_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_multiply_domain_checked");
    }
}

// Power expression inference

/**
 * @brief Helper: check if a NumberNode holds an integer value.
 * Works for BigInt (always integer), Rational (if denominator == 1),
 * and double (if finite and equal to its floor).
 */
bool is_integer_number(const NumberNode& num) {
    if (std::holds_alternative<BigInt>(num.value())) {
        return true;
    }
    if (std::holds_alternative<Rational>(num.value())) {
        return std::get<Rational>(num.value()).is_integer();
    }
    // double case
    double v = std::get<lmmc_real_t>(num.value());
    return std::isfinite(v) && v == std::floor(v);
}

/**
 * @brief Helper: check if a NumberNode holds an even integer value.
 */
bool is_even_integer_number(const NumberNode& num) {
    if (!is_integer_number(num)) return false;

    if (std::holds_alternative<BigInt>(num.value())) {
        return std::get<BigInt>(num.value()).is_even();
    }
    if (std::holds_alternative<Rational>(num.value())) {
        // Rational with denominator 1 - check numerator
        BigInt n = std::get<Rational>(num.value()).get_numerator();
        return n.is_even();
    }
    // double case - guard against overflow for very large values
    double v = std::get<lmmc_real_t>(num.value());
    if (std::fabs(v) >= static_cast<double>(LLONG_MAX)) {
        return std::fmod(v, 2.0) == 0.0;
    }
    long long iv = static_cast<long long>(v);
    return (iv % 2) == 0;
}

/**
 * @brief Helper: check if a NumberNode holds a positive integer value (> 0).
 */
bool is_positive_integer_number(const NumberNode& num) {
    if (!is_integer_number(num)) return false;

    if (std::holds_alternative<BigInt>(num.value())) {
        const auto& b = std::get<BigInt>(num.value());
        return !b.IsNegative() && !(b == BigInt(0));
    }
    if (std::holds_alternative<Rational>(num.value())) {
        BigInt n = std::get<Rational>(num.value()).get_numerator();
        return !n.IsNegative() && !(n == BigInt(0));
    }
    // double case
    double v = std::get<lmmc_real_t>(num.value());
    return v > 0.0;
}

/**
 * @brief Helper: check if a NumberNode holds the value 0.
 */
bool is_zero_number(const NumberNode& num) {
    return num.is_zero();
}


InferenceTriboolResult InferenceEngine::infer_power_property_checked(
    const PowerNode& node, Sign target) const {
    try {
        // Wrap base and exponent as SymbolicExpr for querying
        auto base_expr = LMCAS::detail::expression_from_node(node.base());
        auto exp_expr = LMCAS::detail::expression_from_node(node.exponent());
        // Check if exponent is a NumberNode (needed for several rules)
        auto exp_num = std::dynamic_pointer_cast<const NumberNode>(node.exponent());

        if (exp_num && is_zero_number(*exp_num)) {
            LMCAS_TRY_TRIBOOL(base_nn, query_sign_of_checked(base_expr, Sign::NonNegative));
            if (base_nn == Tribool::True) {
                // Result is Positive (value is 1)
                switch (target) {
                    case Sign::Positive:    return InferenceTriboolResult::success(Tribool::True);
                    case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::True);
                    case Sign::NonZero:     return InferenceTriboolResult::success(Tribool::True);
                    case Sign::Negative:    return InferenceTriboolResult::success(Tribool::False);
                    case Sign::NonPositive: return InferenceTriboolResult::success(Tribool::False);
                    case Sign::Zero:        return InferenceTriboolResult::success(Tribool::False);
                }
            }
        }

        {
            LMCAS_TRY_TRIBOOL(base_pos, query_sign_of_checked(base_expr, Sign::Positive));
            if (base_pos == Tribool::True) {
                LMCAS_TRY_TRIBOOL(exp_real, query_domain_of_checked(exp_expr, Domain::Real));
                if (exp_real == Tribool::True) {
                    switch (target) {
                        case Sign::Positive:    return InferenceTriboolResult::success(Tribool::True);
                        case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::True);
                        case Sign::NonZero:     return InferenceTriboolResult::success(Tribool::True);
                        case Sign::Negative:    return InferenceTriboolResult::success(Tribool::False);
                        case Sign::NonPositive: return InferenceTriboolResult::success(Tribool::False);
                        case Sign::Zero:        return InferenceTriboolResult::success(Tribool::False);
                    }
                }
            }
        }

        if (exp_num && is_even_integer_number(*exp_num)) {
            LMCAS_TRY_TRIBOOL(base_real, query_domain_of_checked(base_expr, Domain::Real));
            if (base_real == Tribool::True) {
                switch (target) {
                    case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::True);
                    case Sign::Negative:    return InferenceTriboolResult::success(Tribool::False);
                    /// base 可能为零,因此 Positive 与 NonPositive 均保持 Unknown.
                    default: break;
                }
            }
        }

        if (exp_num && is_positive_integer_number(*exp_num)) {
            LMCAS_TRY_TRIBOOL(base_nn, query_sign_of_checked(base_expr, Sign::NonNegative));
            if (base_nn == Tribool::True) {
                switch (target) {
                    case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::True);
                    case Sign::Negative:    return InferenceTriboolResult::success(Tribool::False);
                    default: break;
                }
            }
        }

        if (exp_num && is_integer_number(*exp_num)) {
            LMCAS_TRY_TRIBOOL(base_nz, query_sign_of_checked(base_expr, Sign::NonZero));
            if (base_nz == Tribool::True) {
                switch (target) {
                    case Sign::NonZero: return InferenceTriboolResult::success(Tribool::True);
                    case Sign::Zero:    return InferenceTriboolResult::success(Tribool::False);
                    default: break;
                }
            }
        }

        return InferenceTriboolResult::success(Tribool::Unknown);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "power sign inference allocation failed",
            "infer_power_property_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_power_property_checked");
    }
}


InferenceTriboolResult InferenceEngine::infer_power_domain_checked(
    const PowerNode& node, Domain target) const {
    try {
        auto base_expr = LMCAS::detail::expression_from_node(node.base());
        auto exp_num = std::dynamic_pointer_cast<const NumberNode>(node.exponent());

        if (exp_num && is_integer_number(*exp_num)) {
            LMCAS_TRY_TRIBOOL(base_real, query_real_checked(base_expr));
            if (base_real == Tribool::True) {
                switch (target) {
                    case Domain::Real:
                        return InferenceTriboolResult::success(Tribool::True);
                    default:
                        break;
                }
            }

            if (target == Domain::Integer) {
                LMCAS_TRY_TRIBOOL(base_int, query_integer_checked(base_expr));
                if (base_int == Tribool::True && is_positive_integer_number(*exp_num)) {
                    return InferenceTriboolResult::success(Tribool::True);
                }
                if (base_int == Tribool::True && is_zero_number(*exp_num)) {
                    return InferenceTriboolResult::success(Tribool::True);
                }
            }
        }

        return InferenceTriboolResult::success(Tribool::Unknown);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "power domain inference allocation failed",
            "infer_power_domain_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_power_domain_checked");
    }
}

// Function inference


InferenceTriboolResult InferenceEngine::infer_function_property_checked(
    const FunctionNode& node, Sign target) const {
    try {
        // Need at least one argument for all built-in function rules
        if (node.arguments().empty()) return InferenceTriboolResult::success(Tribool::Unknown);

        auto arg_expr = LMCAS::detail::expression_from_node(node.arguments()[0]);
        switch (node.type()) {
            case FunctionNode::FuncType::Exp: {
                // exp(Real) -> Positive
                // exp is always positive for real arguments
                LMCAS_TRY_TRIBOOL(arg_real, query_domain_of_checked(arg_expr, Domain::Real));
                if (arg_real != Tribool::True) return InferenceTriboolResult::success(Tribool::Unknown);

                switch (target) {
                    case Sign::Positive:    return InferenceTriboolResult::success(Tribool::True);
                    case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::True);
                    case Sign::NonZero:     return InferenceTriboolResult::success(Tribool::True);
                    case Sign::Negative:    return InferenceTriboolResult::success(Tribool::False);
                    case Sign::NonPositive: return InferenceTriboolResult::success(Tribool::False);
                    case Sign::Zero:        return InferenceTriboolResult::success(Tribool::False);
                }
                return InferenceTriboolResult::success(Tribool::Unknown);
            }

            case FunctionNode::FuncType::Sin: {
                /// sin(Real) 属于 Real 且值域为 [-1,1];其符号由具体参数值决定.
                LMCAS_TRY_TRIBOOL(arg_real, query_domain_of_checked(arg_expr, Domain::Real));
                if (arg_real != Tribool::True) return InferenceTriboolResult::success(Tribool::Unknown);

                return InferenceTriboolResult::success(Tribool::Unknown);
            }

            case FunctionNode::FuncType::Cos: {
                /// cos(Real) 属于 Real 且值域为 [-1,1];其符号由具体参数值决定.
                LMCAS_TRY_TRIBOOL(arg_real, query_domain_of_checked(arg_expr, Domain::Real));
                if (arg_real != Tribool::True) return InferenceTriboolResult::success(Tribool::Unknown);

                return InferenceTriboolResult::success(Tribool::Unknown);
            }

            case FunctionNode::FuncType::Abs: {
                LMCAS_TRY_TRIBOOL(arg_real, query_domain_of_checked(arg_expr, Domain::Real));
                if (arg_real != Tribool::True) return InferenceTriboolResult::success(Tribool::Unknown);

                // Check if argument is known Positive, Negative, or NonZero
                LMCAS_TRY_TRIBOOL(arg_pos, query_sign_of_checked(arg_expr, Sign::Positive));
                LMCAS_TRY_TRIBOOL(arg_neg, query_sign_of_checked(arg_expr, Sign::Negative));
                LMCAS_TRY_TRIBOOL(arg_nz, query_sign_of_checked(arg_expr, Sign::NonZero));

                // If argument is Positive or Negative (or NonZero), abs is Positive
                bool abs_is_positive = (arg_pos == Tribool::True) ||
                                       (arg_neg == Tribool::True) ||
                                       (arg_nz == Tribool::True);

                switch (target) {
                    case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::True);
                    case Sign::Negative:    return InferenceTriboolResult::success(Tribool::False);
                    case Sign::Positive:
                        return InferenceTriboolResult::success(abs_is_positive ? Tribool::True : Tribool::Unknown);
                    case Sign::NonZero:
                        return InferenceTriboolResult::success(abs_is_positive ? Tribool::True : Tribool::Unknown);
                    case Sign::NonPositive:
                        return InferenceTriboolResult::success(abs_is_positive ? Tribool::False : Tribool::Unknown);
                    case Sign::Zero:
                        return InferenceTriboolResult::success(abs_is_positive ? Tribool::False : Tribool::Unknown);
                }
                return InferenceTriboolResult::success(Tribool::Unknown);
            }

            case FunctionNode::FuncType::Ln: {
                // ln(Positive) -> Real
                /// ln 在 x>1,0<x<1 与 x=1 时分别为正,负与零,通用符号保持 Unknown.
                LMCAS_TRY_TRIBOOL(arg_positive, query_sign_of_checked(arg_expr, Sign::Positive));
                if (arg_positive != Tribool::True) return InferenceTriboolResult::success(Tribool::Unknown);

                // ln(Positive) has no definite sign property
                return InferenceTriboolResult::success(Tribool::Unknown);
            }

            case FunctionNode::FuncType::Sqrt: {
                // sqrt(NonNegative) -> NonNegative, Real
                LMCAS_TRY_TRIBOOL(arg_nonneg, query_sign_of_checked(arg_expr, Sign::NonNegative));
                if (arg_nonneg != Tribool::True) return InferenceTriboolResult::success(Tribool::Unknown);

                switch (target) {
                    case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::True);
                    case Sign::Negative:    return InferenceTriboolResult::success(Tribool::False);
                    case Sign::NonPositive: return InferenceTriboolResult::success(Tribool::Unknown); // sqrt could be zero
                    case Sign::Positive:    return InferenceTriboolResult::success(Tribool::Unknown); // sqrt(0) = 0
                    case Sign::Zero:        return InferenceTriboolResult::success(Tribool::Unknown); // depends on argument
                    case Sign::NonZero:     return InferenceTriboolResult::success(Tribool::Unknown); // depends on argument
                }
                return InferenceTriboolResult::success(Tribool::Unknown);
            }

            case FunctionNode::FuncType::Tan: {
                // tan(Real) -> Real
                // tan can be any real value; no sign determination possible.
                LMCAS_TRY_TRIBOOL(arg_real, query_domain_of_checked(arg_expr, Domain::Real));
                if (arg_real != Tribool::True) return InferenceTriboolResult::success(Tribool::Unknown);

                return InferenceTriboolResult::success(Tribool::Unknown);
            }

            case FunctionNode::FuncType::ArcTan:
                // atan 在实数域严格递增且 atan(0)=0,六类符号性质均与参数一致.
                return query_sign_of_checked(arg_expr, target);

            case FunctionNode::FuncType::ArcSin:
            case FunctionNode::FuncType::ArcCos:
                return InferenceTriboolResult::success(Tribool::Unknown);

            default:
                // Unrecognized function -> Unknown for all properties
                return InferenceTriboolResult::success(Tribool::Unknown);
        }
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "function sign inference allocation failed",
            "infer_function_property_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_function_property_checked");
    }
}


InferenceTriboolResult InferenceEngine::infer_function_domain_checked(
    const FunctionNode& node, Domain target) const {
    try {
        if (node.arguments().empty()) return InferenceTriboolResult::success(Tribool::Unknown);

        auto arg_expr = LMCAS::detail::expression_from_node(node.arguments()[0]);
        switch (node.type()) {
            case FunctionNode::FuncType::Exp: {
                if (target == Domain::Real) {
                    LMCAS_TRY_TRIBOOL(arg_real, query_real_checked(arg_expr));
                    if (arg_real == Tribool::True) return InferenceTriboolResult::success(Tribool::True);

                    LMCAS_TRY_TRIBOOL(arg_int, query_integer_checked(arg_expr));
                    if (arg_int == Tribool::True) return InferenceTriboolResult::success(Tribool::True);

                    LMCAS_TRY_TRIBOOL(arg_rational, query_domain_of_checked(arg_expr, Domain::Rational));
                    if (arg_rational == Tribool::True) return InferenceTriboolResult::success(Tribool::True);
                }
                return InferenceTriboolResult::success(Tribool::Unknown);
            }

            case FunctionNode::FuncType::Sin:
            case FunctionNode::FuncType::Cos:
            case FunctionNode::FuncType::Tan: {
                if (target == Domain::Real) {
                    LMCAS_TRY_TRIBOOL(arg_real, query_real_checked(arg_expr));
                    if (arg_real == Tribool::True) return InferenceTriboolResult::success(Tribool::True);

                    LMCAS_TRY_TRIBOOL(arg_int, query_integer_checked(arg_expr));
                    if (arg_int == Tribool::True) return InferenceTriboolResult::success(Tribool::True);
                }
                return InferenceTriboolResult::success(Tribool::Unknown);
            }

            case FunctionNode::FuncType::Abs: {
                if (target == Domain::Real) {
                    LMCAS_TRY_TRIBOOL(arg_real, query_real_checked(arg_expr));
                    if (arg_real == Tribool::True) return InferenceTriboolResult::success(Tribool::True);

                    LMCAS_TRY_TRIBOOL(arg_int, query_integer_checked(arg_expr));
                    if (arg_int == Tribool::True) return InferenceTriboolResult::success(Tribool::True);
                }
                if (target == Domain::Integer) {
                    LMCAS_TRY_TRIBOOL(arg_int, query_integer_checked(arg_expr));
                    if (arg_int == Tribool::True) return InferenceTriboolResult::success(Tribool::True);
                }
                return InferenceTriboolResult::success(Tribool::Unknown);
            }

            case FunctionNode::FuncType::Ln: {
                if (target == Domain::Real) {
                    LMCAS_TRY_TRIBOOL(arg_positive, query_sign_of_checked(arg_expr, Sign::Positive));
                    if (arg_positive == Tribool::True) return InferenceTriboolResult::success(Tribool::True);

                    LMCAS_TRY_TRIBOOL(arg_int, query_integer_checked(arg_expr));
                    if (arg_int == Tribool::True) return InferenceTriboolResult::success(Tribool::True);
                }
                return InferenceTriboolResult::success(Tribool::Unknown);
            }

            case FunctionNode::FuncType::Sqrt: {
                if (target == Domain::Real) {
                    LMCAS_TRY_TRIBOOL(arg_nonneg, query_sign_of_checked(arg_expr, Sign::NonNegative));
                    if (arg_nonneg == Tribool::True) return InferenceTriboolResult::success(Tribool::True);
                }
                return InferenceTriboolResult::success(Tribool::Unknown);
            }

            case FunctionNode::FuncType::ArcTan:
                if (target == Domain::Real) {
                    LMCAS_TRY_TRIBOOL(arg_real, query_real_checked(arg_expr));
                    if (arg_real == Tribool::True) {
                        return InferenceTriboolResult::success(Tribool::True);
                    }
                }
                return InferenceTriboolResult::success(Tribool::Unknown);


            case FunctionNode::FuncType::ArcSin:
            case FunctionNode::FuncType::ArcCos:
                return InferenceTriboolResult::success(Tribool::Unknown);

            default:
                return InferenceTriboolResult::success(Tribool::Unknown);
        }
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "function domain inference allocation failed",
            "infer_function_domain_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_function_domain_checked");
    }
}

#undef LMCAS_TRY_TRIBOOL

} // namespace LMCAS
