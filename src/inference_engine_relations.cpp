#define _USE_MATH_DEFINES
#include "internal/inference_engine_impl.hpp"
#include "internal/rewrite_budget.hpp"
#include <limits>

namespace LMCAS {

bool InferenceEngine::Impl::operand_sign(const std::shared_ptr<const SymbolicNode>& operand, Sign target, ComputationContext& context) const {
    auto step = context.consume_steps(0, "inference.relations");
    if (!step) throw step.error();
    if (auto var = std::dynamic_pointer_cast<const VariableNode>(operand)) {
        return var->is_constant()
            ? (target == Sign::Positive || target == Sign::NonNegative)
            : ctx.has_sign(var->name(), target);
    }
    if (auto num = std::dynamic_pointer_cast<const NumberNode>(operand)) {
        if (target == Sign::Positive) {
            return num->is_positive();
        }
        return num->is_zero() || num->is_positive();
    }
    return false;
}

Result<void> InferenceEngine::Impl::accumulate_sum_relation_sign(
    const std::shared_ptr<const SymbolicNode>& operand, const SymbolicExpr& zero_expr,
    bool& all_gt_zero, bool& all_geq_zero, ComputationContext& context) const {
    auto step = context.consume_steps(1, "inference.relations");
    if (!step) {
        throw step.error();
    }
    auto op_expr = detail::expression_from_node(operand);
    auto gt = ctx.has_relation_checked(
        op_expr, zero_expr, RelationalNode::Op::GT, context);
    if (!gt) {
        return Result<void>::failure(gt.error());
    }
    auto geq = ctx.has_relation_checked(
        op_expr, zero_expr, RelationalNode::Op::GEQ, context);
    if (!geq) {
        return Result<void>::failure(geq.error());
    }
    bool op_gt = gt.value();
    bool op_geq = geq.value();
    if (!op_gt) {
        op_gt = operand_sign(operand, Sign::Positive, context);
    }
    if (!op_geq && !op_gt) {
        op_geq = operand_sign(operand, Sign::NonNegative, context);
    }
    if (!op_gt) {
        all_gt_zero = false;
    }
    if (!op_gt && !op_geq) {
        all_geq_zero = false;
    }
    return Result<void>::success();
}

InferenceTriboolResult InferenceEngine::Impl::sum_relation_sign(
    const AddNode& add,
    const SymbolicExpr& zero_expr,
    Sign target, ComputationContext& context) const {
    if (add.operands().empty()) {
        return InferenceTriboolResult::success(Tribool::Unknown);
    }

    bool all_gt_zero = true;
    bool all_geq_zero = true;

    for (const auto& operand : add.operands()) {
        auto result = accumulate_sum_relation_sign(
            operand, zero_expr, all_gt_zero, all_geq_zero, context);
        if (!result) {
            return InferenceTriboolResult::failure(result.error());
        }
    }

    if (target == Sign::Positive && all_gt_zero) {
        return InferenceTriboolResult::success(Tribool::True);
    }
    if (target == Sign::NonNegative && all_geq_zero) {
        return InferenceTriboolResult::success(Tribool::True);
    }
    return InferenceTriboolResult::success(Tribool::Unknown);
}

InferenceTriboolResult InferenceEngine::Impl::product_relation_sign(
    const MultiplyNode& mul,
    const SymbolicExpr& zero_expr,
    Sign target, ComputationContext& context) const {
    const auto& rel_store = ctx;
    if (mul.operands().empty()) {
        return InferenceTriboolResult::success(Tribool::Unknown);
    }

    if (target == Sign::Positive) {
        bool all_gt_zero = true;

        for (const auto& operand : mul.operands()) {
            auto step = context.consume_steps(1, "inference.relations");
            if (!step) throw step.error();
            auto op_expr = LMCAS::detail::expression_from_node(operand);
            auto gt = rel_store.has_relation_checked(
                op_expr, zero_expr, RelationalNode::Op::GT, context);
            if (!gt) return InferenceTriboolResult::failure(gt.error());
            bool op_gt = gt.value();

            if (!op_gt) {
                op_gt = operand_sign(operand, Sign::Positive, context);
            }

            if (!op_gt) {
                all_gt_zero = false;
                break;
            }
        }

        if (all_gt_zero) {
            return InferenceTriboolResult::success(Tribool::True);
        }
    }
    return InferenceTriboolResult::success(Tribool::Unknown);
}

Result<bool> InferenceEngine::Impl::relation_rhs_nonnegative(const Relation& rel, const SymbolicExpr& zero_expr, ComputationContext& context) const {
    auto step = context.consume_steps(1, "inference.relations");
    if (!step) return Result<bool>::failure(step.error());
    const auto& rel_store = ctx;
    if (auto rhs_num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(rel.rhs))) {
        if (rhs_num->is_zero() || rhs_num->is_positive()) {
            return true;
        }
    }

    if (auto rhs_var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(rel.rhs))) {
        if (rhs_var->is_constant() ||
            ctx.has_sign(rhs_var->name(), Sign::NonNegative)) {
            return true;
        }
    }

    auto geq = rel_store.has_relation_checked(
        rel.rhs, zero_expr, RelationalNode::Op::GEQ, context);
    if (!geq || geq.value()) return geq;
    return rel_store.has_relation_checked(
        rel.rhs, zero_expr, RelationalNode::Op::GT, context);
}

InferenceTriboolResult InferenceEngine::Impl::positive_relation_chain(
    const SymbolicExpr& expr,
    const SymbolicExpr& zero_expr, ComputationContext& context) const {
    auto relations = ctx.get_visible_relations_checked(context);
    if (!relations) return InferenceTriboolResult::failure(relations.error());
    detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
        std::numeric_limits<std::size_t>::max(), "inference.relations");
    budget.measure(detail::node(expr));
    for (const auto& rel : relations.value()) {
        auto step = context.consume_steps(1, "inference.relations");
        if (!step) throw step.error();
        if (rel.op != RelationalNode::Op::GT) {
            continue;
        }

        if (!LMCAS::detail::node(rel.lhs) || !LMCAS::detail::node(expr)) {
            continue;
        }
        budget.measure(detail::node(rel.lhs));
        if (!LMCAS::detail::node(rel.lhs)->equals(*LMCAS::detail::node(expr))) {
            continue;
        }
        if (!LMCAS::detail::node(rel.rhs)) {
            continue;
        }

        auto nonnegative = relation_rhs_nonnegative(rel, zero_expr, context);
        if (!nonnegative) return InferenceTriboolResult::failure(nonnegative.error());
        if (nonnegative.value()) {
            return InferenceTriboolResult::success(Tribool::True);
        }
    }
    return InferenceTriboolResult::success(Tribool::Unknown);
}

InferenceTriboolResult InferenceEngine::infer_sign_from_relations_checked(
    const SymbolicExpr& expr, Sign target, ComputationContext& context) const {
    if (!LMCAS::detail::node(expr)) {
        return InferenceTriboolResult::failure(
            CasErrc::InvalidArgument,
            "relation sign inference expression must not be null",
            "infer_sign_from_relations_checked");
    }

    try {
        auto zero_expr = LMCAS::detail::expression_from_node(
            LMCAS::detail::make_node<NumberNode>(BigInt(0)));

        if (auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr))) {
            if (add->operands().empty()) {
                return InferenceTriboolResult::success(Tribool::Unknown);
            }
            auto result = impl_->sum_relation_sign(*add, zero_expr, target, context);
            if (!result || result.value() != Tribool::Unknown) {
                return result;
            }
        }

        if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
            if (mul->operands().empty()) {
                return InferenceTriboolResult::success(Tribool::Unknown);
            }
            auto result = impl_->product_relation_sign(*mul, zero_expr, target, context);
            if (!result || result.value() != Tribool::Unknown) {
                return result;
            }
        }

        if (target == Sign::Positive) {
            return impl_->positive_relation_chain(expr, zero_expr, context);
        }

        return InferenceTriboolResult::success(Tribool::Unknown);
    } catch (const CasError& error) {
        return InferenceTriboolResult::failure(error);
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

}
