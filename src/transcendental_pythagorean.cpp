#include "internal/transcendental_support.hpp"
#include "lmmc/numeric_scalar.h"

#include <string>
#include <vector>

namespace LMCAS {
/**
 * @brief 匹配 sin^2(f) 或 cos^2(f)，提取函数类型和参数。
 * @note 对应 PowerNode(FunctionNode(Sin/Cos, [f]), NumberNode(2))。
 *
 * @param[in]  node      待检测的 AST 节点
 * @param[out] func_type 输出函数类型(Sin 或 Cos)
 * @param[out] argument  输出函数参数节点
 * @return 匹配成功返回 true
 * @internal
 */
static bool tf_is_trig_squared(
    const std::shared_ptr<const SymbolicNode>& node,
    FunctionNode::FuncType& func_type,
    std::shared_ptr<const SymbolicNode>& argument) {

    auto pow = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!pow) {
        return false;
    }

    auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent());
    if (!exp_num) {
        return false;
    }

    bool is_two = false;
    if (std::holds_alternative<BigInt>(exp_num->value())) {
        is_two = (std::get<BigInt>(exp_num->value()) == BigInt(2));
    } else if (std::holds_alternative<Rational>(exp_num->value())) {
        is_two = (std::get<Rational>(exp_num->value()) == Rational(2));
    } else if (std::holds_alternative<lmmc_real_t>(exp_num->value())) {
        lmmc_real_t v = std::get<lmmc_real_t>(exp_num->value());
        int eq;
        lmmc_double_nearly_equal(v, 2.0, &eq);
        is_two = (eq != 0);
    }
    if (!is_two) {
        return false;
    }

    auto func = std::dynamic_pointer_cast<const FunctionNode>(pow->base());
    if (!func || func->arguments().size() != 1) {
        return false;
    }

    if (func->type() != FunctionNode::FuncType::Sin &&
        func->type() != FunctionNode::FuncType::Cos) {
        return false;
    }

    func_type = func->type();
    argument = func->arguments()[0];
    return true;
}

/**
 * @brief 从 a*sin^2(f) 或 a*cos^2(f) 提取系数和三角平方项。
 * @note 单独的三角平方项系数为 1。
 *
 * @param[in]  node       待分析的加法操作数节点
 * @param[out] coeff      输出系数节点(nullptr 表示系数为 1)
 * @param[out] func_type  输出函数类型(Sin 或 Cos)
 * @param[out] argument   输出函数参数节点
 * @return 匹配成功返回 true
 * @internal
 */
static bool tf_extract_coeff_trig_squared(
    const std::shared_ptr<const SymbolicNode>& node,
    std::shared_ptr<const SymbolicNode>& coeff,
    FunctionNode::FuncType& func_type,
    std::shared_ptr<const SymbolicNode>& argument) {

    if (tf_is_trig_squared(node, func_type, argument)) {
        coeff = nullptr;  /**< 表示系数为 1。 */
        return true;
    }

    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!mul || mul->operands().size() < 2) {
        return false;
    }

    for (size_t i = 0; i < mul->operands().size(); ++i) {
        if (tf_is_trig_squared(mul->operands()[i], func_type, argument)) {
            std::vector<std::shared_ptr<const SymbolicNode>> coeff_ops;
            for (size_t j = 0; j < mul->operands().size(); ++j) {
                if (j != i) {
                    coeff_ops.push_back(mul->operands()[j]);
                }
            }

            if (coeff_ops.size() == 1) {
                coeff = coeff_ops[0];
            } else {
                coeff = LMCAS::detail::make_node<MultiplyNode>(std::move(coeff_ops));
            }
            return true;
        }
    }

    return false;
}

static std::shared_ptr<const SymbolicNode> tf_simplify_pythagorean_node(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& var);

struct TrigSquare {
    std::shared_ptr<const SymbolicNode> coefficient;
    FunctionNode::FuncType type;
    std::shared_ptr<const SymbolicNode> argument;
};

static bool complementary_square(
    const std::shared_ptr<const SymbolicNode>& node, const TrigSquare& square) {
    TrigSquare other;
    if (!tf_extract_coeff_trig_squared(
            node, other.coefficient, other.type, other.argument)) {
        return false;
    }
    const auto target = square.type == FunctionNode::FuncType::Sin
        ? FunctionNode::FuncType::Cos : FunctionNode::FuncType::Sin;
    if (other.type != target || !square.argument || !other.argument ||
        !square.argument->equals(*other.argument)) {
        return false;
    }
    if (!square.coefficient || !other.coefficient) {
        return !square.coefficient && !other.coefficient;
    }
    return square.coefficient->equals(*other.coefficient);
}

static std::size_t find_complementary_square(
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands,
    const std::vector<bool>& consumed, std::size_t start, const TrigSquare& square) {
    for (std::size_t i = start; i < operands.size(); ++i) {
        if (!consumed[i] && complementary_square(operands[i], square)) {
            return i;
        }
    }
    return operands.size();
}

static std::shared_ptr<const SymbolicNode> simplify_pythagorean_sum(
    const AddNode& sum, const std::string& var) {
    std::vector<std::shared_ptr<const SymbolicNode>> operands;
    operands.reserve(sum.operands().size());
    for (const auto& operand : sum.operands()) {
        operands.push_back(tf_simplify_pythagorean_node(operand, var));
    }
    std::vector<bool> consumed(operands.size(), false);
    std::vector<std::shared_ptr<const SymbolicNode>> result;
    for (std::size_t i = 0; i < operands.size(); ++i) {
        if (consumed[i]) {
            continue;
        }
        TrigSquare square;
        if (!tf_extract_coeff_trig_squared(
                operands[i], square.coefficient, square.type, square.argument)) {
            result.push_back(operands[i]);
            continue;
        }
        const auto partner = find_complementary_square(operands, consumed, i + 1, square);
        if (partner == operands.size()) {
            result.push_back(operands[i]);
            continue;
        }
        consumed[i] = true;
        consumed[partner] = true;
        if (square.coefficient) {
            result.push_back(square.coefficient);
        } else {
            result.push_back(LMCAS::detail::make_node<NumberNode>(BigInt(1)));
        }
    }
    if (result.empty()) {
        return LMCAS::detail::make_node<NumberNode>(BigInt(0));
    }
    if (result.size() == 1) {
        return result[0];
    }
    return LMCAS::detail::make_node<AddNode>(std::move(result));
}

/**
 * @brief 递归化简 a*sin^2(f) + a*cos^2(f) 为 a，系数为 1 时结果为 1。
 * @note 在加法节点中配对，要求参数 f 与系数 a 分别结构相等。
 *
 * @param[in] node 当前 AST 节点
 * @param[in] var  目标变量名(用于限定化简范围)
 * @return 化简后的节点;若无可化简的模式则返回原节点
 * @internal
 */
static std::shared_ptr<const SymbolicNode> tf_simplify_pythagorean_node(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var) {

    if (!node) {
        return node;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return simplify_pythagorean_sum(*add, var);
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        new_ops.reserve(mul->operands().size());
        for (const auto& op : mul->operands()) {
            new_ops.push_back(tf_simplify_pythagorean_node(op, var));
        }
        return LMCAS::detail::make_node<MultiplyNode>(std::move(new_ops));
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto new_base = tf_simplify_pythagorean_node(pow->base(), var);
        auto new_exp = tf_simplify_pythagorean_node(pow->exponent(), var);
        return LMCAS::detail::make_node<PowerNode>(std::move(new_base), std::move(new_exp));
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        new_args.reserve(func->arguments().size());
        for (const auto& arg : func->arguments()) {
            new_args.push_back(tf_simplify_pythagorean_node(arg, var));
        }
        return LMCAS::detail::make_node<FunctionNode>(func->type(), std::move(new_args));
    }

    return node;
}

/**
 * @brief 在因式分解前用 sin^2(f) + cos^2(f) = 1 化简表达式。
 *
 * @param[in] expr 待化简的符号表达式
 * @param[in] var  目标变量名
 * @return 化简后的表达式;若无可化简模式则返回原表达式
 * @internal
 */
std::shared_ptr<SymbolicExpr> tf_simplify_pythagorean(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    if (!expr || !LMCAS::detail::node(expr)) {
        return expr;
    }

    auto simplified_root = tf_simplify_pythagorean_node(LMCAS::detail::node(expr), var);

    if (!simplified_root) {
        return expr;
    }
    if (simplified_root->equals(*LMCAS::detail::node(expr))) {
        return expr;
    }

    return LMCAS::detail::make_expression_ptr(simplified_root);
}

}
