#include "transcendental_factor.hpp"
#include "internal/symbolic_ast.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/transcendental_support.hpp"

#include <string>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <limits>
#include <cstdint>

namespace LMCAS {

/**
 * @brief 判断函数类型是否支持超越函数换元。
 * @param[in] ft 函数类型枚举
 * @return 属于 Sin/Cos/Exp/Ln/Tan 之一返回 true
 * @internal
 */
bool tf_is_transcendental_type(FunctionNode::FuncType ft) {
    return ft == FunctionNode::FuncType::Sin ||
           ft == FunctionNode::FuncType::Cos ||
           ft == FunctionNode::FuncType::Exp ||
           ft == FunctionNode::FuncType::Ln  ||
           ft == FunctionNode::FuncType::Tan;
}

/**
 * @brief 按先序收集依赖 var 的超越函数子表达式。
 *
 * 嵌套函数按外层优先收集，如 sin(exp(x)) 先于 exp(x)。
 *
 * @param[in]  node       当前 AST 节点
 * @param[in]  var        目标变量名
 * @param[out] candidates 收集到的超越子表达式列表
 * @internal
 */
static void tf_collect_transcendental(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var,
    std::vector<std::shared_ptr<const SymbolicNode>>& candidates) {

    if (!node) return;

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (func->arguments().size() == 1 &&
            tf_is_transcendental_type(func->type()) &&
            expression_depends_on_variable(func->arguments()[0], var)) {
            candidates.push_back(node);
        }
        for (auto& arg : func->arguments()) {
            tf_collect_transcendental(arg, var, candidates);
        }
        return;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (auto& op : add->operands()) {
            tf_collect_transcendental(op, var, candidates);
        }
        return;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (auto& op : mul->operands()) {
            tf_collect_transcendental(op, var, candidates);
        }
        return;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        tf_collect_transcendental(pow->base(), var, candidates);
        tf_collect_transcendental(pow->exponent(), var, candidates);
        return;
    }

}

/**
 * @brief 按结构哈希和 equals 比较去重，保留首次出现的顺序。
 *
 * @param[in,out] candidates 候选列表,去重后仅保留唯一项
 * @internal
 */
static void tf_deduplicate(std::vector<std::shared_ptr<const SymbolicNode>>& candidates) {
    NodeSet seen;
    std::vector<std::shared_ptr<const SymbolicNode>> unique;

    for (auto& node : candidates) {
        if (seen.find(node) == seen.end()) {
            seen.insert(node);
            unique.push_back(node);
        }
    }

    candidates = std::move(unique);
}

/**
 * @brief 判断两个 AST 节点是否结构相等。
 * @param[in] a 第一个节点
 * @param[in] b 第二个节点
 * @return 结构相等返回 true
 * @internal
 */
static bool tf_nodes_equal(const std::shared_ptr<const SymbolicNode>& a,
                           const std::shared_ptr<const SymbolicNode>& b) {
    if (!a || !b) return a == b;
    return a->equals(*b);
}

/**
 * @brief 匹配 (-1)*arg 或 arg*(-1) 的二因子乘积结构。
 *
 * @param[in] node 待检测节点
 * @param[in] arg  参考参数节点
 * @return 若 node 表示 -arg 则返回 true
 * @internal
 */
static bool tf_is_negation_of(const std::shared_ptr<const SymbolicNode>& node,
                              const std::shared_ptr<const SymbolicNode>& arg) {
    if (!node || !arg) {
        return false;
    }

    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!mul || mul->operands().size() != 2) {
        return false;
    }

    for (int i = 0; i < 2; ++i) {
        auto num_node = std::dynamic_pointer_cast<const NumberNode>(mul->operands()[i]);
        if (!num_node) {
            continue;
        }

        bool is_neg_one = false;
        if (std::holds_alternative<BigInt>(num_node->value())) {
            is_neg_one = (std::get<BigInt>(num_node->value()) == BigInt(-1));
        } else if (std::holds_alternative<Rational>(num_node->value())) {
            is_neg_one = (std::get<Rational>(num_node->value()) == Rational(-1));
        } else if (std::holds_alternative<lmmc_real_t>(num_node->value())) {
            lmmc_real_t v = std::get<lmmc_real_t>(num_node->value());
            int eq;
            lmmc_double_nearly_equal(v, -1.0, &eq);
            is_neg_one = (eq != 0);
        }

        if (is_neg_one) {
            int other_idx = 1 - i;
            return tf_nodes_equal(mul->operands()[other_idx], arg);
        }
    }

    return false;
}

/**
 * @brief 按映射顺序将超越子表达式替换为不定元，递归重建未匹配节点。
 *
 * 外层映射优先：sin(exp(x)) 整体替换为 u0，exp(x) 保留在 u0 的定义中。
 *
 * @param[in] node     当前 AST 节点
 * @param[in] mappings 换元映射列表(按分配顺序)
 * @return 替换后的新节点
 * @internal
 */
static std::shared_ptr<const SymbolicNode> tf_substitute_expr(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::vector<TransSubstitution>& mappings) {

    if (!node) return nullptr;

    for (const auto& m : mappings) {
        if (m.trans_expr && LMCAS::detail::node(m.trans_expr) &&
            node->equals(*LMCAS::detail::node(m.trans_expr))) {
            return LMCAS::detail::make_node<VariableNode>(m.indeterminate);
        }
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        new_args.reserve(func->arguments().size());
        for (const auto& arg : func->arguments()) {
            new_args.push_back(tf_substitute_expr(arg, mappings));
        }
        return LMCAS::detail::make_node<FunctionNode>(func->type(), std::move(new_args));
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        new_ops.reserve(add->operands().size());
        for (const auto& op : add->operands()) {
            new_ops.push_back(tf_substitute_expr(op, mappings));
        }
        return LMCAS::detail::make_node<AddNode>(std::move(new_ops));
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        new_ops.reserve(mul->operands().size());
        for (const auto& op : mul->operands()) {
            new_ops.push_back(tf_substitute_expr(op, mappings));
        }
        return LMCAS::detail::make_node<MultiplyNode>(std::move(new_ops));
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto new_base = tf_substitute_expr(pow->base(), mappings);
        auto new_exp = tf_substitute_expr(pow->exponent(), mappings);
        return LMCAS::detail::make_node<PowerNode>(std::move(new_base), std::move(new_exp));
    }

    return node;
}

static bool tf_is_sine_cosine_pair(FunctionNode::FuncType first, FunctionNode::FuncType second) {
    if (first == FunctionNode::FuncType::Sin) {
        return second == FunctionNode::FuncType::Cos;
    }
    return first == FunctionNode::FuncType::Cos && second == FunctionNode::FuncType::Sin;
}

static std::shared_ptr<SymbolicExpr> tf_pair_constraint(
    const TransSubstitution& left, const TransSubstitution& right,
    const FunctionNode& first, const FunctionNode& second) {
    const auto first_type = first.type();
    const auto second_type = second.type();
    const bool sin_cos = tf_is_sine_cosine_pair(first_type, second_type);
    if (sin_cos && tf_nodes_equal(first.arguments()[0], second.arguments()[0])) {
        auto ui = SymbolicExpr::variable(left.indeterminate);
        auto uj = SymbolicExpr::variable(right.indeterminate);
        auto two = SymbolicExpr::number(2);
        auto sum = SymbolicExpr::add(SymbolicExpr::power(ui, two), SymbolicExpr::power(uj, two));
        return SymbolicExpr::add(sum, SymbolicExpr::number(-1));
    }
    if (first_type != FunctionNode::FuncType::Exp ||
        second_type != FunctionNode::FuncType::Exp) {
        return nullptr;
    }
    if (!tf_is_negation_of(second.arguments()[0], first.arguments()[0]) &&
        !tf_is_negation_of(first.arguments()[0], second.arguments()[0])) {
        return nullptr;
    }
    auto ui = SymbolicExpr::variable(left.indeterminate);
    auto uj = SymbolicExpr::variable(right.indeterminate);
    return SymbolicExpr::add(SymbolicExpr::multiply(ui, uj), SymbolicExpr::number(-1));
}

static void tf_detect_constraints(TransSubstitutionResult& result) {
    const auto& mappings = result.mappings;
    for (size_t i = 0; i < mappings.size(); ++i) {
        auto first = std::dynamic_pointer_cast<const FunctionNode>(
            detail::node(mappings[i].trans_expr));
        if (!first || first->arguments().size() != 1) continue;
        for (size_t j = i + 1; j < mappings.size(); ++j) {
            auto second = std::dynamic_pointer_cast<const FunctionNode>(
                detail::node(mappings[j].trans_expr));
            if (!second || second->arguments().size() != 1) continue;
            auto constraint = tf_pair_constraint(mappings[i], mappings[j], *first, *second);
            if (constraint) result.constraints.push_back(std::move(constraint));
        }
    }
}

/**
 * @brief 为依赖 var 的 Sin、Cos、Exp、Ln、Tan 子表达式构造换元结果。
 *
 * 按结构去重后，依次分配不定元 u0、u1、u2 等。
 *
 * @param[in] expr 待检测的符号表达式
 * @param[in] var  目标变量名
 * @return 换元映射、代数约束及替换后的表达式；无匹配时保留原表达式。
 */
TransSubstitutionResult detect_trans_substitutions(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    TransSubstitutionResult result;
    result.poly_expr = nullptr;

    if (!expr || !LMCAS::detail::node(expr)) {
        return result;
    }

    std::vector<std::shared_ptr<const SymbolicNode>> candidates;
    tf_collect_transcendental(LMCAS::detail::node(expr), var, candidates);

    tf_deduplicate(candidates);

    for (size_t i = 0; i < candidates.size(); ++i) {
        TransSubstitution mapping;
        mapping.trans_expr = LMCAS::detail::make_expression_ptr(candidates[i]);
        mapping.indeterminate = "u" + std::to_string(i);
        result.mappings.push_back(std::move(mapping));
    }

    tf_detect_constraints(result);

    if (!result.mappings.empty()) {
        auto substituted = tf_substitute_expr(LMCAS::detail::node(expr), result.mappings);
        result.poly_expr = LMCAS::detail::make_expression_ptr(substituted);
    } else {
        result.poly_expr = expr;
    }

    return result;
}

}
