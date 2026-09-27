#define _USE_MATH_DEFINES
#include "internal/inference_engine_impl.hpp"

namespace LMCAS {

Tribool InferenceEngine::Impl::known_sign(int value, Sign target) {
    switch (target) {
        case Sign::Positive: {
            return value > 0 ? Tribool::True : Tribool::False;
        }
        case Sign::Negative: {
            return value < 0 ? Tribool::True : Tribool::False;
        }
        case Sign::NonNegative: {
            return value >= 0 ? Tribool::True : Tribool::False;
        }
        case Sign::NonPositive: {
            return value <= 0 ? Tribool::True : Tribool::False;
        }
        case Sign::Zero: {
            return value == 0 ? Tribool::True : Tribool::False;
        }
        case Sign::NonZero: {
            return value != 0 ? Tribool::True : Tribool::False;
        }
    }
    return Tribool::Unknown;
}

Tribool InferenceEngine::Impl::nonnegative_sign(Sign target) {
    if (target == Sign::NonNegative) {
        return Tribool::True;
    }
    if (target == Sign::Negative) {
        return Tribool::False;
    }
    return Tribool::Unknown;
}

std::shared_ptr<const SymbolicNode> InferenceEngine::Impl::subtrahend(
    const std::shared_ptr<const SymbolicNode>& operand, ComputationContext& context) {
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(operand);
    if (!mul || mul->operands().size() != 2) {
        return nullptr;
    }
    bool found_neg_one = false;
    std::shared_ptr<const SymbolicNode> other_operand;
    for (const auto& mul_op : mul->operands()) {
        auto step = context.consume_steps(1, "inference.arithmetic");
        if (!step) throw step.error();
        const auto number =
            std::dynamic_pointer_cast<const NumberNode>(mul_op);
        if (number && number->is_negative_one()) {
            found_neg_one = true;
            continue;
        }
        other_operand = mul_op;
    }
    return found_neg_one ? other_operand : nullptr;
}

bool InferenceEngine::Impl::has_negated_operand(const AddNode& node, ComputationContext& context) {
    if (node.operands().size() != 2) {
        return false;
    }
    for (const auto& operand : node.operands()) {
        auto step = context.consume_steps(1, "inference.arithmetic");
        if (!step) throw step.error();
        auto mul = std::dynamic_pointer_cast<const MultiplyNode>(operand);
        if (!mul || mul->operands().size() != 2) {
            continue;
        }
        for (const auto& factor : mul->operands()) {
            auto step = context.consume_steps(1, "inference.arithmetic");
            if (!step) throw step.error();
            const auto number =
                std::dynamic_pointer_cast<const NumberNode>(factor);
            if (number && number->is_negative_one()) {
                return true;
            }
        }
    }
    return false;
}

InferenceTriboolResult InferenceEngine::Impl::sum_sign(
    const InferenceEngine& engine,
    const AddNode& node,
    Sign target, ComputationContext& context) const {
    bool all_have_property = true;

    for (const auto& operand : node.operands()) {
        auto step = context.consume_steps(1, "inference.arithmetic");
        if (!step) throw step.error();
        auto op_expr = LMCAS::detail::expression_from_node(operand);
        Tribool op_result = Tribool::Unknown;
        if (target != Sign::Positive && target != Sign::Negative &&
            target != Sign::NonNegative && target != Sign::NonPositive) {
            return InferenceTriboolResult::success(Tribool::Unknown);
        }
        auto queried = engine.query_sign_of_checked(op_expr, target, context);
        if (!queried) {
            return queried;
        }
        op_result = queried.value();

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
}

InferenceTriboolResult InferenceEngine::Impl::subtraction_sign(
    const InferenceEngine& engine,
    const SymbolicExpr& minuend_expr,
    const SymbolicExpr& subtrahend_expr,
    Sign target, ComputationContext& context) const {
    auto min_pos_result =
        engine.query_sign_of_checked(minuend_expr, Sign::Positive, context);
    if (!min_pos_result) {
        return min_pos_result;
    }
    const Tribool min_pos = min_pos_result.value();
    auto min_neg_result =
        engine.query_sign_of_checked(minuend_expr, Sign::Negative, context);
    if (!min_neg_result) {
        return min_neg_result;
    }
    const Tribool min_neg = min_neg_result.value();
    auto sub_pos_result =
        engine.query_sign_of_checked(subtrahend_expr, Sign::Positive, context);
    if (!sub_pos_result) {
        return sub_pos_result;
    }
    const Tribool sub_pos = sub_pos_result.value();
    auto sub_neg_result =
        engine.query_sign_of_checked(subtrahend_expr, Sign::Negative, context);
    if (!sub_neg_result) {
        return sub_neg_result;
    }
    const Tribool sub_neg = sub_neg_result.value();

    if (min_pos == Tribool::True && sub_neg == Tribool::True) {
        return known_sign(1, target);
    }

    if (min_neg == Tribool::True && sub_pos == Tribool::True) {
        return known_sign(-1, target);
    }

    return weak_subtraction_sign(engine, minuend_expr, subtrahend_expr, target, context);
}

InferenceTriboolResult InferenceEngine::Impl::weak_subtraction_sign(
    const InferenceEngine& engine,
    const SymbolicExpr& minuend_expr,
    const SymbolicExpr& subtrahend_expr,
    Sign target, ComputationContext& context) const {
    auto min_nn_result =
        engine.query_sign_of_checked(minuend_expr, Sign::NonNegative, context);
    if (!min_nn_result) {
        return min_nn_result;
    }
    const Tribool min_nn = min_nn_result.value();
    auto sub_np_result =
        engine.query_sign_of_checked(subtrahend_expr, Sign::NonPositive, context);
    if (!sub_np_result) {
        return sub_np_result;
    }
    const Tribool sub_np = sub_np_result.value();
    if (min_nn == Tribool::True && sub_np == Tribool::True) {
        switch (target) {
            case Sign::NonNegative: return InferenceTriboolResult::success(Tribool::True);
            case Sign::Negative:    return InferenceTriboolResult::success(Tribool::False);
            default: break;
        }
    }

    auto min_np_result =
        engine.query_sign_of_checked(minuend_expr, Sign::NonPositive, context);
    if (!min_np_result) {
        return min_np_result;
    }
    const Tribool min_np = min_np_result.value();
    auto sub_nn_result =
        engine.query_sign_of_checked(subtrahend_expr, Sign::NonNegative, context);
    if (!sub_nn_result) {
        return sub_nn_result;
    }
    const Tribool sub_nn = sub_nn_result.value();
    if (min_np == Tribool::True && sub_nn == Tribool::True) {
        switch (target) {
            case Sign::NonPositive: return InferenceTriboolResult::success(Tribool::True);
            case Sign::Positive:    return InferenceTriboolResult::success(Tribool::False);
            default: break;
        }
    }

    return InferenceTriboolResult::success(Tribool::Unknown);
}

InferenceTriboolResult InferenceEngine::infer_add_sign_checked(
    const AddNode& node, Sign target, ComputationContext& context) const {
    try {
        if (node.operands().empty()) {
            return InferenceTriboolResult::success(Tribool::Unknown);
        }
        if (Impl::has_negated_operand(node, context)) {
            auto sub_result = infer_subtraction_sign_checked(node, target, context);
            if (!sub_result) {
                return sub_result;
            }
            if (sub_result.value() != Tribool::Unknown) {
                return sub_result;
            }
        }

        /**
         * @brief 各加数均具备目标符号属性时，和也具备该属性。
         * @note 含 Unknown 或已知异号加数时返回 Unknown。
         */

        return impl_->sum_sign(*this, node, target, context);
    } catch (const CasError& error) {
        return InferenceTriboolResult::failure(error);
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

InferenceTriboolResult InferenceEngine::infer_subtraction_sign_checked(
    const AddNode& node, Sign target, ComputationContext& context) const {
    try {
        if (node.operands().size() != 2) {
            return InferenceTriboolResult::success(Tribool::Unknown);
        }

        std::shared_ptr<const SymbolicNode> minuend_node;
        std::shared_ptr<const SymbolicNode> subtrahend_node;

        for (const auto& operand : node.operands()) {
            auto step = context.consume_steps(1, "inference.arithmetic");
            if (!step) throw step.error();
            auto inner = Impl::subtrahend(operand, context);
            if (inner) {
                subtrahend_node = inner;
            } else {
                minuend_node = operand;
            }
        }

        if (!minuend_node || !subtrahend_node) {
            return InferenceTriboolResult::success(Tribool::Unknown);
        }

        auto minuend_expr = LMCAS::detail::expression_from_node(minuend_node);
        auto subtrahend_expr = LMCAS::detail::expression_from_node(subtrahend_node);

        return impl_->subtraction_sign(*this, minuend_expr, subtrahend_expr, target, context);
    } catch (const CasError& error) {
        return InferenceTriboolResult::failure(error);
    } catch (const std::bad_alloc&) {
        return InferenceTriboolResult::failure(
            CasErrc::ResourceLimit,
            "subtraction sign inference allocation failed",
            "infer_subtraction_sign_checked");
    } catch (const std::exception& ex) {
        return InferenceTriboolResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            "infer_subtraction_sign_checked");
    }
}

} // namespace LMCAS
