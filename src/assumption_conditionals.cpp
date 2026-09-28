#include "assumption_context.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/normalization_utils.hpp"
#include "internal/expression_analysis.hpp"

namespace LMCAS {

namespace {


Tribool evaluate_symbol_sign(const AssumptionContext& ctx,
                            const std::string& name, RelationOp op) {
    struct SignRule {
        RelationOp op;
        std::initializer_list<Sign> true_signs;
        std::initializer_list<Sign> false_signs;
    };
    const SignRule rules[] = {
        {RelationOp::GT, {Sign::Positive},
         {Sign::Negative, Sign::NonPositive, Sign::Zero}},
        {RelationOp::GEQ, {Sign::NonNegative, Sign::Positive, Sign::Zero},
         {Sign::Negative}},
        {RelationOp::LT, {Sign::Negative},
         {Sign::Positive, Sign::NonNegative, Sign::Zero}},
        {RelationOp::LEQ, {Sign::NonPositive, Sign::Negative, Sign::Zero},
         {Sign::Positive}},
        {RelationOp::NEQ, {Sign::NonZero, Sign::Positive, Sign::Negative},
         {Sign::Zero}},
        {RelationOp::EQ, {Sign::Zero},
         {Sign::Positive, Sign::Negative, Sign::NonZero}}
    };
    for (const auto& rule : rules) {
        if (rule.op != op) continue;
        for (Sign sign : rule.true_signs) {
            if (ctx.has_sign(name, sign)) return Tribool::True;
        }
        for (Sign sign : rule.false_signs) {
            if (ctx.has_sign(name, sign)) return Tribool::False;
        }
        break;
    }
    return Tribool::Unknown;
}

Tribool evaluate_zero_comparison(const AssumptionContext& ctx,
                                const RelationalNode& relation) {
    const auto lhs_var = std::dynamic_pointer_cast<const VariableNode>(relation.left());
    const auto rhs_num = std::dynamic_pointer_cast<const NumberNode>(relation.right());
    if (lhs_var && !lhs_var->is_constant() && rhs_num && rhs_num->is_zero()) {
        return evaluate_symbol_sign(ctx, lhs_var->name(), relation.op());
    }
    const auto lhs_num = std::dynamic_pointer_cast<const NumberNode>(relation.left());
    const auto rhs_var = std::dynamic_pointer_cast<const VariableNode>(relation.right());
    if (lhs_num && lhs_num->is_zero() && rhs_var && !rhs_var->is_constant()) {
        return evaluate_symbol_sign(
            ctx, rhs_var->name(), detail::reversed_relation(relation.op()));
    }
    return Tribool::Unknown;
}

}


AssumptionVoidResult AssumptionContext::assume_conditional(
    const SymbolicExpr& condition,
    const SymbolicExpr& conclusion) {
    return assume_conditional_checked(condition, conclusion);
}

AssumptionVoidResult AssumptionContext::assume_conditional_checked(
    const SymbolicExpr& condition,
    const SymbolicExpr& conclusion) {
    constexpr const char* operation = "assume_conditional";
    if (!std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(condition))) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "condition expression must be relational", operation);
    }
    if (!std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(conclusion))) {
        return AssumptionVoidResult::failure(
            CasErrc::InvalidArgument, "conclusion expression must be relational", operation);
    }
    try {
        Tribool cond_result = evaluate_condition(condition);
        if (cond_result == Tribool::True &&
            evaluate_condition(conclusion) == Tribool::False) {
            return AssumptionVoidResult::failure(
                CasErrc::InvalidArgument,
                "condition is satisfied but conclusion contradicts the current assumption state",
                operation);
        }
        scope_stack_.back().conditionals.push_back({condition, conclusion});
        ++cache_generation_;
    } catch (const std::bad_alloc&) {
        return AssumptionVoidResult::failure(
            CasErrc::ResourceLimit, "assumption allocation failed", operation);
    } catch (const std::invalid_argument& ex) {
        return AssumptionVoidResult::failure(CasErrc::InvalidArgument, ex.what(), operation);
    } catch (const std::exception& ex) {
        return AssumptionVoidResult::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
    return AssumptionVoidResult::success();
}

std::vector<AssumptionContext::ConditionalAssumption> AssumptionContext::get_active_conditionals() const {
    std::vector<ConditionalAssumption> result;
    for (auto it = scope_stack_.rbegin(); it != scope_stack_.rend(); ++it) {
        for (const auto& cond : it->conditionals) {
            result.push_back(cond);
        }
    }
    return result;
}

Tribool AssumptionContext::evaluate_condition(const SymbolicExpr& condition) const {
    if (!LMCAS::detail::node(condition)) {
        return Tribool::Unknown;
    }
    auto rel_node = std::dynamic_pointer_cast<const RelationalNode>(LMCAS::detail::node(condition));
    if (!rel_node) {
        return Tribool::Unknown;
    }
    const auto* lhs_number = dynamic_cast<const NumberNode*>(rel_node->left().get());
    const auto* rhs_number = dynamic_cast<const NumberNode*>(rel_node->right().get());
    if (lhs_number && rhs_number) {
        const auto lhs = exact_number_as_rational(*lhs_number);
        const auto rhs = exact_number_as_rational(*rhs_number);
        bool value = false;
        switch (rel_node->op()) {
            case RelationOp::GT: { value = lhs > rhs; break; }
            case RelationOp::GEQ: { value = lhs >= rhs; break; }
            case RelationOp::LT: { value = lhs < rhs; break; }
            case RelationOp::LEQ: { value = lhs <= rhs; break; }
            case RelationOp::EQ: { value = lhs == rhs; break; }
            case RelationOp::NEQ: { value = lhs != rhs; break; }
        }
        return value ? Tribool::True : Tribool::False;
    }

    auto lhs = LMCAS::detail::expression_from_node(rel_node->left());
    auto rhs = LMCAS::detail::expression_from_node(rel_node->right());
    RelationalNode::Op op = rel_node->op();
    if (has_relation(lhs, rhs, op)) { return Tribool::True; }
    return evaluate_zero_comparison(*this, *rel_node);
}

}
