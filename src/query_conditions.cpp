#include "internal/query_support.hpp"

namespace LMCAS {

namespace {

using ConditionSet = QueryInterface::ConditionSet;
using ConditionSets = std::vector<ConditionSet>;

std::shared_ptr<const SymbolicNode> subtracted_operand(
    const std::shared_ptr<const SymbolicNode>& node) {
    const auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!mul || mul->operands().size() != 2) return nullptr;
    for (std::size_t j = 0; j < 2; ++j) {
        const auto coefficient =
            std::dynamic_pointer_cast<const NumberNode>(mul->operands()[j]);
        if (coefficient && coefficient->is_negative_one()) {
            return mul->operands()[1 - j];
        }
    }
    return nullptr;
}

ConditionSets difference_conditions(
    const std::shared_ptr<const SymbolicNode>& pos_operand,
    const std::shared_ptr<const SymbolicNode>& neg_operand, Sign target) {
    const auto lhs_var = std::dynamic_pointer_cast<const VariableNode>(pos_operand);
    const auto rhs_var = std::dynamic_pointer_cast<const VariableNode>(neg_operand);
    if (!lhs_var || !rhs_var) return {};
    if (target == Sign::Positive) {
        ConditionSet cs1;
        cs1.sign_conditions.emplace_back(lhs_var->name(), Sign::Positive);
        cs1.sign_conditions.emplace_back(rhs_var->name(), Sign::Negative);
        ConditionSet cs2;
        auto lhs_expr = LMCAS::detail::expression_from_node(pos_operand);
        auto rhs_expr = LMCAS::detail::expression_from_node(neg_operand);
        Relation rel{lhs_expr, rhs_expr, RelationalNode::Op::GT};
        cs2.relational_conditions.push_back(rel);
        cs2.sign_conditions.emplace_back(rhs_var->name(), Sign::NonNegative);

        return {cs1, cs2};
    } else if (target == Sign::Negative) {
        ConditionSet cs1;
        cs1.sign_conditions.emplace_back(lhs_var->name(), Sign::Negative);
        cs1.sign_conditions.emplace_back(rhs_var->name(), Sign::Positive);

        ConditionSet cs2;
        auto lhs_expr = LMCAS::detail::expression_from_node(pos_operand);
        auto rhs_expr = LMCAS::detail::expression_from_node(neg_operand);
        Relation rel{rhs_expr, lhs_expr, RelationalNode::Op::GT};
        cs2.relational_conditions.push_back(rel);
        cs2.sign_conditions.emplace_back(lhs_var->name(), Sign::NonNegative);

        return {cs1, cs2};
    } else if (target == Sign::NonNegative) {
        ConditionSet cs1;
        cs1.sign_conditions.emplace_back(lhs_var->name(), Sign::NonNegative);
        cs1.sign_conditions.emplace_back(rhs_var->name(), Sign::NonPositive);

        ConditionSet cs2;
        auto lhs_expr = LMCAS::detail::expression_from_node(pos_operand);
        auto rhs_expr = LMCAS::detail::expression_from_node(neg_operand);
        Relation rel{lhs_expr, rhs_expr, RelationalNode::Op::GEQ};
        cs2.relational_conditions.push_back(rel);

        return {cs1, cs2};
    }
    return {};
}

ConditionSets variable_operand_conditions(
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands, Sign target) {
    ConditionSet cs;
    for (const auto& operand : operands) {
        const auto variable = std::dynamic_pointer_cast<const VariableNode>(operand);
        if (!variable) return {};
        cs.sign_conditions.emplace_back(variable->name(), target);
    }
    if (cs.sign_conditions.empty()) return {};
    return {cs};
}

ConditionSets sum_conditions(const AddNode& add, Sign target) {
    if (add.operands().size() == 2) {
        for (std::size_t i = 0; i < 2; ++i) {
            const auto negative = subtracted_operand(add.operands()[i]);
            if (!negative) continue;
            auto result = difference_conditions(add.operands()[1 - i], negative, target);
            if (!result.empty()) return result;
            break;
        }
    }
    if (target == Sign::Positive || target == Sign::NonNegative) {
        return variable_operand_conditions(add.operands(), target);
    }
    return {};
}

ConditionSets quotient_conditions(
    const std::shared_ptr<const SymbolicNode>& numerator,
    const std::shared_ptr<const SymbolicNode>& denominator, Sign target) {
    const auto num_var = std::dynamic_pointer_cast<const VariableNode>(numerator);
    const auto den_var = std::dynamic_pointer_cast<const VariableNode>(denominator);
    if (!num_var || !den_var) return {};
    if (target == Sign::Positive) {
        ConditionSet cs1;
        cs1.sign_conditions.emplace_back(num_var->name(), Sign::Positive);
        cs1.sign_conditions.emplace_back(den_var->name(), Sign::Positive);

        ConditionSet cs2;
        cs2.sign_conditions.emplace_back(num_var->name(), Sign::Negative);
        cs2.sign_conditions.emplace_back(den_var->name(), Sign::Negative);

        return {cs1, cs2};
    } else if (target == Sign::Negative) {
        ConditionSet cs1;
        cs1.sign_conditions.emplace_back(num_var->name(), Sign::Positive);
        cs1.sign_conditions.emplace_back(den_var->name(), Sign::Negative);

        ConditionSet cs2;
        cs2.sign_conditions.emplace_back(num_var->name(), Sign::Negative);
        cs2.sign_conditions.emplace_back(den_var->name(), Sign::Positive);

        return {cs1, cs2};
    }
    return {};
}

ConditionSets product_conditions(const MultiplyNode& mul, Sign target) {
    if (mul.operands().size() == 2) {
        for (std::size_t i = 0; i < 2; ++i) {
            const auto power = std::dynamic_pointer_cast<const PowerNode>(mul.operands()[i]);
            const auto exponent =
                power ? std::dynamic_pointer_cast<const NumberNode>(
                            power->exponent())
                      : nullptr;
            if (exponent && exponent->is_negative_one()) {
                return quotient_conditions(mul.operands()[1 - i], power->base(), target);
            }
        }
    }
    if (target == Sign::Positive) {
        return variable_operand_conditions(mul.operands(), Sign::Positive);
    }
    return {};
}

}

QueryInterface::QueryConditionSetsResult QueryInterface::query_conditions(
    const SymbolicExpr& expr, Sign target) const {
    return query_conditions_checked(expr, target);
}

std::vector<QueryInterface::ConditionSet> QueryInterface::query_conditions_impl(
    const SymbolicExpr& expr, Sign target) const {
    if (!LMCAS::detail::node(expr)) {
        return {};
    }
    if (auto var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
        ConditionSet cs;
        cs.sign_conditions.emplace_back(var->name(), target);
        return {cs};
    }
    if (const auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr))) {
        return sum_conditions(*add, target);
    }
    if (const auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        return product_conditions(*mul, target);
    }
    return {};
}

QueryInterface::QueryConditionSetsResult QueryInterface::query_conditions_checked(
    const SymbolicExpr& expr, Sign target) const {
    return query_detail::checked_expression_result<std::vector<ConditionSet>>(
        expr, "query_conditions", [&]() {
            return query_conditions_impl(expr, target);
        });
}

}
