#include "internal/integration_support.hpp"

namespace LMCAS {

/**
 * @brief 在积分变量已知 Positive 时，递归将 |var| 替换为 var。
 * 仅替换 FunctionNode(Abs, [arg]) 中 arg 恰为积分变量的节点；
 * 调用方须由 AssumptionContext 确认该变量为正。
 */
static std::shared_ptr<const SymbolicNode> simplify_abs_positive(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& var);

static std::shared_ptr<const SymbolicNode> simplify_abs_function(
    const std::shared_ptr<const SymbolicNode>& node, const FunctionNode& fn,
    const std::string& var) {
        if (fn.type() == FunctionNode::FuncType::Abs && fn.arguments().size() == 1) {
            if (auto vn = std::dynamic_pointer_cast<const VariableNode>(fn.arguments()[0])) {
                if (!vn->is_constant() && vn->name() == var) {
                    return fn.arguments()[0];
                }
            }
        }
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        bool changed = false;
        for (auto& op : fn.arguments()) {
            auto new_op = simplify_abs_positive(op, var);
            if (new_op != op) changed = true;
            new_ops.push_back(new_op);
        }
        if (changed) {
            return LMCAS::detail::make_node<FunctionNode>(fn.type(), new_ops);
        }
        return node;
}

static std::shared_ptr<const SymbolicNode> simplify_abs_positive(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    if (!node) { return node; }

    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return simplify_abs_function(node, *fn, var);
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        bool changed = false;
        for (auto& op : add->operands()) {
            auto new_op = simplify_abs_positive(op, var);
            if (new_op != op) changed = true;
            new_ops.push_back(new_op);
        }
        if (changed) { return LMCAS::detail::make_node<AddNode>(new_ops); }
        return node;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        bool changed = false;
        for (auto& op : mul->operands()) {
            auto new_op = simplify_abs_positive(op, var);
            if (new_op != op) changed = true;
            new_ops.push_back(new_op);
        }
        if (changed) { return LMCAS::detail::make_node<MultiplyNode>(new_ops); }
        return node;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto new_base = simplify_abs_positive(pow->base(), var);
        auto new_exp = simplify_abs_positive(pow->exponent(), var);
        if (new_base != pow->base() || new_exp != pow->exponent()) {
            return LMCAS::detail::make_node<PowerNode>(new_base, new_exp);
        }
        return node;
    }

    return node;
}

/**
 * @brief 根据假设化简被积函数。
 * AssumptionContext 确认积分变量为 Positive 时，将各处 |var| 替换为 var；
 * 无适用化简时返回原表达式。
 */
Result<SymbolicExpr> apply_assumption_simplifications(
    const SymbolicExpr& expr, const std::string& var,
    const AssumptionContext* ctx) {
    if (!ctx) { return expr; }

    SymbolicExpr var_expr = *SymbolicExpr::variable(var);
    auto positive = ctx->is_positive(var_expr);
    if (!positive) {
        return Result<SymbolicExpr>::failure(positive.error());
    }

    if (positive.value() == Tribool::True) {
        auto new_root = simplify_abs_positive(
            LMCAS::detail::node(expr), var);
        if (new_root != LMCAS::detail::node(expr)) {
            return LMCAS::detail::expression_from_node(new_root);
        }
    }

    return expr;
}

}
