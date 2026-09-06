/**
 * @file transcendental_factor.cpp
 * @brief 混合超越方程不可约因式分解:换元检测与主入口实现.
 *
 * 本文件实现 Phase 1(换元检测)的核心逻辑:遍历表达式 AST,
 * 收集依赖目标变量的超越子表达式,去重后分配代数不定元.
 */

#include "transcendental_factor.hpp"
#include "symbolic_ast.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"

#include <string>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <limits>

#include "internal/transcendental_support.hpp"

namespace LMCAS {

/**
 * @brief 将因子表达式中的不定元变量替换回原始超越子表达式.
 *
 * 对换元后的因子执行逆操作:遍历映射列表,将每个不定元(u0, u1, ...)
 * 替换为其对应的原始超越表达式(如 sin(x),cos(x) 等).
 * 利用 SymbolicExpr::substitute() 逐一执行变量替换.
 *
 * @param[in] factor_expr 以不定元表示的因子表达式
 * @param[in] mappings    换元映射列表(indeterminate -> trans_expr)
 * @return 替换后的符号表达式,以原始变量和超越函数表示
 * @internal
 */
std::shared_ptr<SymbolicExpr> tf_back_substitute(
    const std::shared_ptr<SymbolicExpr>& factor_expr,
    const std::vector<TransSubstitution>& mappings) {

    if (!factor_expr || !LMCAS::detail::node(factor_expr)) return factor_expr;
    if (mappings.empty()) return factor_expr;

    auto result = factor_expr;

    for (const auto& m : mappings) {
        if (!m.trans_expr || m.indeterminate.empty()) continue;

        /// 仅当表达式依赖该不定元时才执行替换
        if (expression_depends_on_variable(LMCAS::detail::node(result), m.indeterminate)) {
            result = result->substitute(m.indeterminate, m.trans_expr);
        }
    }

    return result;
}


/**
 * @brief 从 NumberNode 中提取有理数值.
 *
 * 将 BigInt,Rational,lmmc_real_t 统一转换为 Rational 表示.
 * 对于浮点数,仅当其为精确整数时才转换;否则返回失败.
 *
 * @param[in]  num_node 数值节点
 * @param[out] out      输出的有理数值
 * @return 提取成功返回 true
 * @internal
 */
static bool tf_extract_rational(const std::shared_ptr<const NumberNode>& num_node, Rational& out) {
    if (!num_node) return false;

    if (std::holds_alternative<BigInt>(num_node->value())) {
        out = Rational(std::get<BigInt>(num_node->value()));
        return true;
    }
    if (std::holds_alternative<Rational>(num_node->value())) {
        out = std::get<Rational>(num_node->value());
        return true;
    }
    if (std::holds_alternative<lmmc_real_t>(num_node->value())) {
        lmmc_real_t v = std::get<lmmc_real_t>(num_node->value());
        /// 仅处理精确整数浮点值
        if (std::isfinite(v) && v == std::floor(v) && std::abs(v) < 1e15) {
            out = Rational(BigInt(static_cast<long long>(v)));
            return true;
        }
        return false;
    }
    return false;
}

/**
 * @brief 对逆换元后的因子列表执行化简与常数乘子提取.
 *
 * 算法:
 * 1. 对每个因子调用 simplify() 进行规范化
 * 2. 若化简后的因子为纯数值(NumberNode),将其累乘到常数积中
 * 3. 若化简后的因子为乘积形式(MultiplyNode)且含数值前导系数,
 *    提取该系数并保留非常数部分
 * 4. 若累积常数 != 1,将其作为独立数值因子插入结果列表首位
 *
 * @param[in,out] factors 因子列表,就地修改为化简后的结果
 * @return 化简并提取常数后的因子列表
 */
std::vector<std::shared_ptr<SymbolicExpr>> tf_simplify_factors(
    std::vector<std::shared_ptr<SymbolicExpr>>& factors) {

    std::vector<std::shared_ptr<SymbolicExpr>> result;
    Rational constant_product(1);

    for (auto& factor : factors) {
        if (!factor || !LMCAS::detail::node(factor)) continue;

        /// 调用 simplify() 规范化因子
        auto simplified = factor->simplify();
        if (!simplified || !LMCAS::detail::node(simplified)) {
            simplified = factor;
        }

        /// 情形 1:因子为纯数值
        if (auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(simplified))) {
            Rational val;
            if (tf_extract_rational(num, val)) {
                if (val != Rational(0)) {
                    constant_product = constant_product * val;
                }
                /// 零因子不累乘,但保留(整个乘积为零)
                else {
                    result.clear();
                    result.push_back(SymbolicExpr::number(0));
                    return result;
                }
            } else {
                /// 当前有理数提取规则之外的数值节点保持原结构.
                result.push_back(simplified);
            }
            continue;
        }

        /// 情形 2:因子为乘积形式,检查是否含数值前导系数
        if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(simplified))) {
            std::vector<std::shared_ptr<const SymbolicNode>> numeric_ops;
            std::vector<std::shared_ptr<const SymbolicNode>> non_numeric_ops;

            for (const auto& op : mul->operands()) {
                if (auto num = std::dynamic_pointer_cast<const NumberNode>(op)) {
                    numeric_ops.push_back(op);
                } else {
                    non_numeric_ops.push_back(op);
                }
            }

            /// 提取所有数值操作数为常数乘子
            if (!numeric_ops.empty() && !non_numeric_ops.empty()) {
                for (const auto& nop : numeric_ops) {
                    auto num = std::dynamic_pointer_cast<const NumberNode>(nop);
                    Rational val;
                    if (tf_extract_rational(num, val)) {
                        constant_product = constant_product * val;
                    } else {
                        /// 当前提取规则之外的数值操作数保留在符号部分.
                        non_numeric_ops.push_back(nop);
                    }
                }

                /// 重建非常数部分
                if (non_numeric_ops.size() == 1) {
                    result.push_back(LMCAS::detail::make_expression_ptr(non_numeric_ops[0]));
                } else {
                    result.push_back(LMCAS::detail::make_expression_ptr(
                        LMCAS::detail::make_node<MultiplyNode>(std::move(non_numeric_ops))));
                }
            } else if (numeric_ops.empty()) {
                /// 无数值操作数,保留原因子
                result.push_back(simplified);
            } else {
                /// 全为数值操作数:整个因子为常数
                Rational val(1);
                for (const auto& nop : numeric_ops) {
                    auto num = std::dynamic_pointer_cast<const NumberNode>(nop);
                    Rational v;
                    if (tf_extract_rational(num, v)) {
                        val = val * v;
                    }
                }
                constant_product = constant_product * val;
            }
            continue;
        }

        /// 情形 3:非数值,非乘积形式,直接保留
        result.push_back(simplified);
    }

    /// 若累积常数 != 1,插入为首个因子
    if (constant_product != Rational(1)) {
        auto const_expr = SymbolicExpr::number(constant_product);
        result.insert(result.begin(), const_expr);
    }

    return result;
}


/**
 * @brief 判断换元后的表达式是否对所有不定元和原始变量均为线性.
 *
 * 换元表达式对每个不定元(u0, u1, ...)及原始变量的次数均小于等于 1 时,
 * 该表达式在超越多项式环中为整体元素.
 * 例如 a*sin(x) + b*x + c 映射为 a*u0 + b*x + c,对 u0 与 x 均为线性.
 *
 * @param[in] sub_result 换元结果(含 poly_expr 和 mappings)
 * @param[in] var        原始目标变量名
 * @return 表达式对所有变量均为线性返回 true
 * @internal
 */
bool tf_is_linear_irreducible(
    const TransSubstitutionResult& sub_result,
    const std::string& var) {

    if (!sub_result.poly_expr || !LMCAS::detail::node(sub_result.poly_expr)) return false;
    if (sub_result.mappings.empty()) return false;

    const auto& root = LMCAS::detail::node(sub_result.poly_expr);

    /// 检查每个不定元的次数是否 <= 1
    for (const auto& m : sub_result.mappings) {
        int deg = tf_degree_in(root, m.indeterminate);
        if (deg < 0 || deg > 1) return false;
    }

    /// 检查原始变量的次数是否 <= 1
    if (expression_depends_on_variable(root, var)) {
        int deg = tf_degree_in(root, var);
        if (deg < 0 || deg > 1) return false;
    }

    return true;
}

/**
 * @brief 检测表达式是否已为独立子表达式的乘积形式.
 *
 * MultiplyNode 的各操作数直接形成独立因子,沿乘法结构完成分解;
 * 数值常数单独累积,并在值异于 1 时形成常数因子.
 *
 * @param[in] expr 待检测的符号表达式
 * @return 因子列表;若表达式非乘积形式则返回空向量(表示无快速路径)
 * @internal
 */
std::vector<std::shared_ptr<SymbolicExpr>> tf_detect_multiplicative_structure(
    const std::shared_ptr<SymbolicExpr>& expr) {

    if (!expr || !LMCAS::detail::node(expr)) return {};

    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr));
    if (!mul || mul->operands().size() < 2) return {};

    std::vector<std::shared_ptr<SymbolicExpr>> factors;
    Rational constant_acc(1);

    for (const auto& op : mul->operands()) {
        if (!op) continue;

        /// 数值常数单独累积
        if (op->is_number()) {
            auto num = std::dynamic_pointer_cast<const NumberNode>(op);
            if (num) {
                if (std::holds_alternative<BigInt>(num->value())) {
                    constant_acc = constant_acc * Rational(std::get<BigInt>(num->value()));
                } else if (std::holds_alternative<Rational>(num->value())) {
                    constant_acc = constant_acc * std::get<Rational>(num->value());
                } else {
                    /// 浮点数值:作为独立因子保留
                    factors.push_back(LMCAS::detail::make_expression_ptr(op));
                }
            }
            continue;
        }

        /// 非数值操作数作为独立因子
        factors.push_back(LMCAS::detail::make_expression_ptr(op));
    }

    /// 仅当存在至少两个非常数因子(或一个非常数因子加一个非 1 常数)时才视为有效乘积分解
    if (factors.size() < 2 && (factors.empty() || constant_acc == Rational(1))) {
        return {};
    }

    /// 插入累积常数因子(若非 1)
    if (constant_acc != Rational(1)) {
        auto const_expr = SymbolicExpr::number(constant_acc);
        factors.insert(factors.begin(), const_expr);
    }

    return factors;
}

/**
 * @brief 从乘积项中提取指数函数因子.
 *
 * 若节点本身为 exp(f(x)) 形式,直接返回该节点.
 * 若节点为 MultiplyNode,遍历其操作数寻找 exp(f(x)) 因子.
 * 仅提取第一个匹配的指数函数因子.
 *
 * @param[in] node 待检测的 AST 节点
 * @param[in] var  目标变量名
 * @return 找到的 exp 因子节点;未找到返回 nullptr
 * @internal
 */
static std::shared_ptr<const SymbolicNode> tf_extract_exp_factor(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var) {

    if (!node) return nullptr;

    /// 直接为 exp(f(x)) 形式
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (func->type() == FunctionNode::FuncType::Exp &&
            func->arguments().size() == 1 &&
            expression_depends_on_variable(func->arguments()[0], var)) {
            return node;
        }
    }

    /// 乘积形式:遍历操作数寻找 exp 因子
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            auto func = std::dynamic_pointer_cast<const FunctionNode>(op);
            if (func &&
                func->type() == FunctionNode::FuncType::Exp &&
                func->arguments().size() == 1 &&
                expression_depends_on_variable(func->arguments()[0], var)) {
                return op;
            }
        }
    }

    return nullptr;
}

/**
 * @brief 从乘积项中移除指定的指数函数因子,返回剩余部分.
 *
 * 若节点本身即为该 exp 因子,返回数值 1.
 * 若节点为 MultiplyNode,移除匹配的 exp 操作数后重建乘积.
 *
 * @param[in] node       原始乘积项节点
 * @param[in] exp_factor 待移除的 exp 因子节点
 * @return 移除 exp 因子后的剩余节点
 * @internal
 */
static std::shared_ptr<const SymbolicNode> tf_remove_exp_factor(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::shared_ptr<const SymbolicNode>& exp_factor) {

    if (!node || !exp_factor) return node;

    /// 节点本身即为 exp 因子
    if (node->equals(*exp_factor)) {
        return LMCAS::detail::make_node<NumberNode>(BigInt(1));
    }

    /// 乘积形式:移除匹配的操作数
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> remaining_ops;
        bool removed = false;

        for (const auto& op : mul->operands()) {
            if (!removed && op->equals(*exp_factor)) {
                removed = true;
                continue;
            }
            remaining_ops.push_back(op);
        }

        if (!removed) return node;

        if (remaining_ops.empty()) {
            return LMCAS::detail::make_node<NumberNode>(BigInt(1));
        }
        if (remaining_ops.size() == 1) {
            return remaining_ops[0];
        }
        return LMCAS::detail::make_node<MultiplyNode>(std::move(remaining_ops));
    }

    return node;
}

/**
 * @brief 检测加法表达式中的公共指数因子并执行分离.
 *
 * 对于 AddNode 形式的表达式,检查所有加法项是否共享相同的 exp(f(x)) 因子.
 * 若是,则提取公因子:expr = exp(f(x)) * (t1' + t2' + ... + tn'),
 * 其中 ti' = ti / exp(f(x)).
 *
 * 典型用例:
 * - exp(x)*x + exp(x) -> [exp(x), x+1]
 * - exp(x)*x^2 + 2*exp(x)*x + exp(x) -> [exp(x), x^2+2x+1]
 *
 * @param[in] expr 待检测的符号表达式
 * @param[in] var  目标变量名
 * @return 因子列表 [exp(f(x)), remaining_sum];若无公共 exp 因子则返回空向量
 * @internal
 */
std::vector<std::shared_ptr<SymbolicExpr>> tf_detect_exponential_separation(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    if (!expr || !LMCAS::detail::node(expr)) return {};

    auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr));
    if (!add || add->operands().size() < 2) return {};

    /// 从第一个加法项中提取 exp 因子作为候选公因子
    std::shared_ptr<const SymbolicNode> common_exp = tf_extract_exp_factor(add->operands()[0], var);
    if (!common_exp) return {};

    /// 验证所有加法项均含有相同的 exp 因子
    for (size_t i = 1; i < add->operands().size(); ++i) {
        std::shared_ptr<const SymbolicNode> term_exp = tf_extract_exp_factor(add->operands()[i], var);
        if (!term_exp || !term_exp->equals(*common_exp)) {
            return {};
        }
    }

    /// 所有项共享相同的 exp(f(x)),执行分离
    /// 构造剩余和:对每个项移除 exp 因子
    std::vector<std::shared_ptr<const SymbolicNode>> remainder_terms;
    remainder_terms.reserve(add->operands().size());

    for (const auto& op : add->operands()) {
        auto remainder = tf_remove_exp_factor(op, common_exp);
        remainder_terms.push_back(remainder);
    }

    /// 构造结果
    auto exp_factor_expr = LMCAS::detail::make_expression_ptr(common_exp);

    std::shared_ptr<SymbolicExpr> sum_expr;
    if (remainder_terms.size() == 1) {
        sum_expr = LMCAS::detail::make_expression_ptr(remainder_terms[0]);
    } else {
        sum_expr = LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<AddNode>(std::move(remainder_terms)));
    }

    /// 化简剩余和
    auto simplified_sum = sum_expr->simplify();
    if (simplified_sum && LMCAS::detail::node(simplified_sum)) {
        sum_expr = simplified_sum;
    }

    return {exp_factor_expr, sum_expr};
}


/**
 * @brief 判断表达式 AST 中是否包含依赖指定变量的超越函数.
 *
 * 递归遍历 AST,若发现任何 FunctionNode 类型为 Sin/Cos/Exp/Ln/Tan
 * 且其参数依赖 var,则返回 true.
 *
 * @param[in] node 当前 AST 节点
 * @param[in] var  目标变量名
 * @return 包含超越函数返回 true
 * @internal
 */
bool tf_contains_transcendental(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var) {

    if (!node) return false;

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (func->arguments().size() == 1 &&
            tf_is_transcendental_type(func->type()) &&
            expression_depends_on_variable(func->arguments()[0], var)) {
            return true;
        }
        for (const auto& arg : func->arguments()) {
            if (tf_contains_transcendental(arg, var)) return true;
        }
        return false;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& op : add->operands()) {
            if (tf_contains_transcendental(op, var)) return true;
        }
        return false;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            if (tf_contains_transcendental(op, var)) return true;
        }
        return false;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return tf_contains_transcendental(pow->base(), var) ||
               tf_contains_transcendental(pow->exponent(), var);
    }

    return false;
}


/**
 * @brief 判断节点是否为 sin^2(f) 或 cos^2(f) 形式,并提取函数类型和参数.
 *
 * 匹配模式:PowerNode(FunctionNode(Sin/Cos, [f]), NumberNode(2))
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
    if (!pow) return false;

    /// 检查指数是否为 2
    auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent());
    if (!exp_num) return false;

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
    if (!is_two) return false;

    /// 检查底数是否为 sin 或 cos
    auto func = std::dynamic_pointer_cast<const FunctionNode>(pow->base());
    if (!func || func->arguments().size() != 1) return false;

    if (func->type() != FunctionNode::FuncType::Sin &&
        func->type() != FunctionNode::FuncType::Cos) {
        return false;
    }

    func_type = func->type();
    argument = func->arguments()[0];
    return true;
}

/**
 * @brief 从乘积项中提取系数和 sin^2/cos^2 核心部分.
 *
 * 对于形如 a*sin^2(f) 的项,提取系数 a 和 sin^2(f) 部分.
 * 若项本身即为 sin^2(f),系数为 1.
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

    /// 直接为 sin^2(f) 或 cos^2(f)
    if (tf_is_trig_squared(node, func_type, argument)) {
        coeff = nullptr;  // 系数为 1
        return true;
    }

    /// 乘积形式:a * sin^2(f) 或 sin^2(f) * a
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!mul || mul->operands().size() < 2) return false;

    /// 在操作数中寻找 sin^2(f) 或 cos^2(f) 部分
    for (size_t i = 0; i < mul->operands().size(); ++i) {
        if (tf_is_trig_squared(mul->operands()[i], func_type, argument)) {
            /// 提取剩余操作数作为系数
            std::vector<std::shared_ptr<const SymbolicNode>> coeff_ops;
            for (size_t j = 0; j < mul->operands().size(); ++j) {
                if (j != i) coeff_ops.push_back(mul->operands()[j]);
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

/**
 * @brief 对表达式执行毕达哥拉斯恒等式化简:sin^2(f) + cos^2(f) -> 1.
 *
 * 递归遍历 AST,在每个 AddNode 中扫描操作数对,寻找满足以下模式的配对:
 * - sin^2(f) + cos^2(f) -> 替换为 1
 * - a*sin^2(f) + a*cos^2(f) -> 替换为 a(公共系数)
 *
 * 要求 sin^2 和 cos^2 的参数 f 结构相等,且公共系数结构相等.
 *
 * @param[in] node 当前 AST 节点
 * @param[in] var  目标变量名(用于限定化简范围)
 * @return 化简后的节点;若无可化简的模式则返回原节点
 * @internal
 */
static std::shared_ptr<const SymbolicNode> tf_simplify_pythagorean_node(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var) {

    if (!node) return node;

    /// 递归处理子节点
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        /// 先递归化简每个操作数
        std::vector<std::shared_ptr<const SymbolicNode>> simplified_ops;
        simplified_ops.reserve(add->operands().size());
        for (const auto& op : add->operands()) {
            simplified_ops.push_back(tf_simplify_pythagorean_node(op, var));
        }

        /// 在化简后的操作数中寻找 sin^2(f) + cos^2(f) 配对
        std::vector<bool> consumed(simplified_ops.size(), false);
        std::vector<std::shared_ptr<const SymbolicNode>> result_ops;

        for (size_t i = 0; i < simplified_ops.size(); ++i) {
            if (consumed[i]) continue;

            std::shared_ptr<const SymbolicNode> coeff_i;
            FunctionNode::FuncType type_i;
            std::shared_ptr<const SymbolicNode> arg_i;

            if (!tf_extract_coeff_trig_squared(simplified_ops[i], coeff_i, type_i, arg_i)) {
                result_ops.push_back(simplified_ops[i]);
                continue;
            }

            /// 确定配对目标类型
            FunctionNode::FuncType target_type =
                (type_i == FunctionNode::FuncType::Sin)
                    ? FunctionNode::FuncType::Cos
                    : FunctionNode::FuncType::Sin;

            bool found_pair = false;
            for (size_t j = i + 1; j < simplified_ops.size(); ++j) {
                if (consumed[j]) continue;

                std::shared_ptr<const SymbolicNode> coeff_j;
                FunctionNode::FuncType type_j;
                std::shared_ptr<const SymbolicNode> arg_j;

                if (!tf_extract_coeff_trig_squared(simplified_ops[j], coeff_j, type_j, arg_j)) {
                    continue;
                }

                /// 检查类型互补且参数相同
                if (type_j != target_type) continue;
                if (!arg_i || !arg_j || !arg_i->equals(*arg_j)) continue;

                /// 检查系数相等
                bool coeffs_equal = false;
                if (!coeff_i && !coeff_j) {
                    /// 两者系数均为 1
                    coeffs_equal = true;
                } else if (coeff_i && coeff_j) {
                    coeffs_equal = coeff_i->equals(*coeff_j);
                }

                if (!coeffs_equal) continue;

                /// 找到配对:sin^2(f) + cos^2(f) -> 1,或 a*sin^2(f) + a*cos^2(f) -> a
                consumed[i] = true;
                consumed[j] = true;
                found_pair = true;

                if (!coeff_i) {
                    /// 系数为 1:替换为 NumberNode(1)
                    result_ops.push_back(LMCAS::detail::make_node<NumberNode>(BigInt(1)));
                } else {
                    /// 有公共系数:替换为系数本身
                    result_ops.push_back(coeff_i);
                }
                break;
            }

            if (!found_pair) {
                result_ops.push_back(simplified_ops[i]);
            }
        }

        /// 重建 AddNode
        if (result_ops.empty()) {
            return LMCAS::detail::make_node<NumberNode>(BigInt(0));
        }
        if (result_ops.size() == 1) {
            return result_ops[0];
        }
        return LMCAS::detail::make_node<AddNode>(std::move(result_ops));
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

    /// 叶节点(NumberNode,VariableNode)直接返回
    return node;
}

/**
 * @brief 对符号表达式执行毕达哥拉斯恒等式化简预处理.
 *
 * 在因式分解主流程之前调用,将表达式中的 sin^2(f) + cos^2(f) 子模式
 * 替换为 1,减少后续换元和分解的复杂度.
 *
 * @param[in] expr 待化简的符号表达式
 * @param[in] var  目标变量名
 * @return 化简后的表达式;若无可化简模式则返回原表达式
 * @internal
 */
std::shared_ptr<SymbolicExpr> tf_simplify_pythagorean(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    if (!expr || !LMCAS::detail::node(expr)) return expr;

    auto simplified_root = tf_simplify_pythagorean_node(LMCAS::detail::node(expr), var);

    if (!simplified_root) return expr;

    /// 化简结果与原表达式结构相同时复用原对象.
    if (simplified_root->equals(*LMCAS::detail::node(expr))) {
        return expr;
    }

    return LMCAS::detail::make_expression_ptr(simplified_root);
}

} // namespace LMCAS
