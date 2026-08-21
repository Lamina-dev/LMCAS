/**
 * @file transcendental_factor.cpp
 * @brief 混合超越方程不可约因式分解：换元检测与主入口实现。
 *
 * 本文件实现 Phase 1（换元检测）的核心逻辑：遍历表达式 AST，
 * 收集依赖目标变量的超越子表达式，去重后分配代数不定元。
 */

#include "transcendental_factor.hpp"
#include "symbolic_ast.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/transcendental_support.hpp"

#include <string>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <limits>

namespace lamina {


/**
 * @brief 判断 FunctionNode 的函数类型是否为需要换元的超越函数。
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
 * @brief 递归遍历 AST，收集所有依赖目标变量的超越函数子表达式。
 *
 * 对于每个 FunctionNode，若其类型为超越函数且参数依赖 var，
 * 则将该节点作为候选加入结果集。遍历采用先序策略：
 * 先检查当前节点是否为超越函数，若是则收集后继续递归其参数
 * （以支持嵌套超越函数如 sin(exp(x)) 的外层优先收集）。
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
        /// 继续递归参数，以收集嵌套的超越子表达式
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

    /// NumberNode、VariableNode 等叶节点无需递归
}

/**
 * @brief 对候选超越子表达式去重，保留结构唯一的节点。
 *
 * 使用 NodeSet（基于结构哈希和 equals 比较）进行去重，
 * 保持首次出现的顺序。
 *
 * @param[in,out] candidates 候选列表，去重后仅保留唯一项
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
 * @brief 判断节点是否为给定参数的取负形式。
 *
 * 检测 node 是否等价于 -arg，即结构为 MultiplyNode({NumberNode(-1), arg})
 * 或 MultiplyNode({arg, NumberNode(-1)})。
 *
 * @param[in] node 待检测节点
 * @param[in] arg  参考参数节点
 * @return 若 node 表示 -arg 则返回 true
 * @internal
 */
static bool tf_is_negation_of(const std::shared_ptr<const SymbolicNode>& node,
                              const std::shared_ptr<const SymbolicNode>& arg) {
    if (!node || !arg) return false;

    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!mul || mul->operands().size() != 2) return false;

    /// 检查两种排列：(-1)*arg 或 arg*(-1)
    for (int i = 0; i < 2; ++i) {
        auto num_node = std::dynamic_pointer_cast<const NumberNode>(mul->operands()[i]);
        if (!num_node) continue;

        /// 检查数值是否为 -1
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
 * @brief 递归替换 AST 中的超越子表达式为对应不定元变量节点。
 *
 * 对每个节点，按映射顺序（u0 优先于 u1）检查是否与某个超越子表达式结构相等。
 * 若匹配则替换为 VariableNode；否则递归处理子节点并重建当前节点。
 * 由于映射按外层优先顺序排列，先匹配外层可确保 sin(exp(x)) 整体被替换为 u0，
 * 而非先将内部 exp(x) 替换为 u1。
 *
 * @param[in] node     当前 AST 节点
 * @param[in] mappings 换元映射列表（按分配顺序）
 * @return 替换后的新节点
 * @internal
 */
static std::shared_ptr<const SymbolicNode> tf_substitute_expr(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::vector<TransSubstitution>& mappings) {

    if (!node) return nullptr;

    /// 按映射顺序检查当前节点是否匹配某个超越子表达式
    for (const auto& m : mappings) {
        if (m.trans_expr && lamina::detail::node(m.trans_expr) &&
            node->equals(*lamina::detail::node(m.trans_expr))) {
            return lamina::detail::make_node<VariableNode>(m.indeterminate);
        }
    }

    /// 对各节点类型递归处理子节点
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        new_args.reserve(func->arguments().size());
        for (const auto& arg : func->arguments()) {
            new_args.push_back(tf_substitute_expr(arg, mappings));
        }
        return lamina::detail::make_node<FunctionNode>(func->type(), std::move(new_args));
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        new_ops.reserve(add->operands().size());
        for (const auto& op : add->operands()) {
            new_ops.push_back(tf_substitute_expr(op, mappings));
        }
        return lamina::detail::make_node<AddNode>(std::move(new_ops));
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        new_ops.reserve(mul->operands().size());
        for (const auto& op : mul->operands()) {
            new_ops.push_back(tf_substitute_expr(op, mappings));
        }
        return lamina::detail::make_node<MultiplyNode>(std::move(new_ops));
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto new_base = tf_substitute_expr(pow->base(), mappings);
        auto new_exp = tf_substitute_expr(pow->exponent(), mappings);
        return lamina::detail::make_node<PowerNode>(std::move(new_base), std::move(new_exp));
    }

    /// 叶节点（NumberNode、VariableNode）直接返回
    return node;
}

/**
 * @brief 检测换元映射之间的代数约束关系。
 *
 * 遍历所有映射对，识别以下约束模式：
 * - sin(f) 与 cos(f) 同时出现（相同参数 f）→ u_sin² + u_cos² - 1 = 0
 * - exp(f) 与 exp(-f) 同时出现（互为取负参数）→ u_pos * u_neg - 1 = 0
 *
 * 约束以 SymbolicExpr 形式存储，表示等于零的多项式。
 *
 * @param[in,out] result 换元结果，constraints 字段将被填充
 * @internal
 */
static void tf_detect_constraints(TransSubstitutionResult& result) {
    const auto& mappings = result.mappings;
    size_t n = mappings.size();

    for (size_t i = 0; i < n; ++i) {
        auto node_i = lamina::detail::node(mappings[i].trans_expr);
        auto func_i = std::dynamic_pointer_cast<const FunctionNode>(node_i);
        if (!func_i || func_i->arguments().size() != 1) continue;

        for (size_t j = i + 1; j < n; ++j) {
            auto node_j = lamina::detail::node(mappings[j].trans_expr);
            auto func_j = std::dynamic_pointer_cast<const FunctionNode>(node_j);
            if (!func_j || func_j->arguments().size() != 1) continue;

            /// --- sin/cos Pythagorean 约束 ---
            bool is_sin_cos = (func_i->type() == FunctionNode::FuncType::Sin &&
                               func_j->type() == FunctionNode::FuncType::Cos) ||
                              (func_i->type() == FunctionNode::FuncType::Cos &&
                               func_j->type() == FunctionNode::FuncType::Sin);

            if (is_sin_cos &&
                tf_nodes_equal(func_i->arguments()[0], func_j->arguments()[0])) {
                /// 构造约束：u_i² + u_j² - 1 = 0
                auto ui = SymbolicExpr::variable(mappings[i].indeterminate);
                auto uj = SymbolicExpr::variable(mappings[j].indeterminate);
                auto two = SymbolicExpr::number(2);

                auto ui_sq = SymbolicExpr::power(ui, two);
                auto uj_sq = SymbolicExpr::power(uj, two);
                auto sum = SymbolicExpr::add(ui_sq, uj_sq);
                auto constraint = SymbolicExpr::add(
                    sum, SymbolicExpr::number(-1));

                result.constraints.push_back(constraint);
                continue;
            }

            /// --- exp/exp(-) 逆元约束 ---
            if (func_i->type() == FunctionNode::FuncType::Exp &&
                func_j->type() == FunctionNode::FuncType::Exp) {

                bool is_inverse_pair = false;

                /// 检查 arg_j == -arg_i 或 arg_i == -arg_j
                if (tf_is_negation_of(func_j->arguments()[0], func_i->arguments()[0]) ||
                    tf_is_negation_of(func_i->arguments()[0], func_j->arguments()[0])) {
                    is_inverse_pair = true;
                }

                if (is_inverse_pair) {
                    /// 构造约束：u_i * u_j - 1 = 0
                    auto ui = SymbolicExpr::variable(mappings[i].indeterminate);
                    auto uj = SymbolicExpr::variable(mappings[j].indeterminate);

                    auto product = SymbolicExpr::multiply(ui, uj);
                    auto constraint = SymbolicExpr::add(
                        product, SymbolicExpr::number(-1));

                    result.constraints.push_back(constraint);
                }
            }
        }
    }
}


/**
 * @brief 检测表达式中的超越函数换元模式。
 *
 * 遍历表达式 AST，识别所有依赖目标变量的超越子表达式（Sin、Cos、Exp、Ln、Tan），
 * 对结构相同的子表达式去重后，为每个分配唯一的代数不定元名（u0, u1, u2, ...）。
 *
 * @param[in] expr 待检测的符号表达式
 * @param[in] var  目标变量名
 * @return 换元结果，包含映射表；poly_expr 和 constraints 留空待后续任务填充
 */
TransSubstitutionResult detect_trans_substitutions(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    TransSubstitutionResult result;
    result.poly_expr = nullptr;

    if (!expr || !lamina::detail::node(expr)) {
        return result;
    }

    /// 收集所有依赖 var 的超越子表达式
    std::vector<std::shared_ptr<const SymbolicNode>> candidates;
    tf_collect_transcendental(lamina::detail::node(expr), var, candidates);

    /// 结构去重
    tf_deduplicate(candidates);

    /// 为每个唯一的超越子表达式分配不定元
    for (size_t i = 0; i < candidates.size(); ++i) {
        TransSubstitution mapping;
        mapping.trans_expr = lamina::detail::make_expression_ptr(candidates[i]);
        mapping.indeterminate = "u" + std::to_string(i);
        result.mappings.push_back(std::move(mapping));
    }

    /// 检测映射之间的代数约束
    tf_detect_constraints(result);

    /// 执行替换：将超越子表达式替换为对应不定元变量
    if (!result.mappings.empty()) {
        auto substituted = tf_substitute_expr(lamina::detail::node(expr), result.mappings);
        result.poly_expr = lamina::detail::make_expression_ptr(substituted);
    } else {
        /// 无超越子表达式时，poly_expr 即为原表达式
        result.poly_expr = expr;
    }

    return result;
}


/**
 * @brief 计算符号表达式关于指定变量的次数。
 *
 * 递归遍历 AST，计算表达式展开后关于 var 的最高次幂。
 * 对于不依赖 var 的子表达式返回 0；对于无法确定次数的情形返回 -1。
 *
 * @param[in] node 符号节点
 * @param[in] var  变量名
 * @return 关于 var 的次数；-1 表示无法确定（非多项式结构）
 * @internal
 */
int tf_degree_in(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    if (!node) return 0;

    if (!expression_depends_on_variable(node, var)) return 0;

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        return (v->name() == var) ? 1 : 0;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        int max_deg = 0;
        for (const auto& op : add->operands()) {
            int d = tf_degree_in(op, var);
            if (d < 0) return -1;
            max_deg = std::max(max_deg, d);
        }
        return max_deg;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        int total_deg = 0;
        for (const auto& op : mul->operands()) {
            int d = tf_degree_in(op, var);
            if (d < 0) return -1;
            total_deg += d;
        }
        return total_deg;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        /// 指数必须为非负整数常量
        auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent());
        if (!exp_num) return -1;

        int e_val = -1;
        if (std::holds_alternative<BigInt>(exp_num->value())) {
            const auto& bi = std::get<BigInt>(exp_num->value());
            if (!bi.IsNegative()) e_val = bi.to_int();
        } else if (std::holds_alternative<Rational>(exp_num->value())) {
            const auto& r = std::get<Rational>(exp_num->value());
            if (r.is_integer()) {
                BigInt bi = r.to_BigInt();
                if (!bi.IsNegative()) e_val = bi.to_int();
            }
        } else if (std::holds_alternative<lmmc_real_t>(exp_num->value())) {
            lmmc_real_t d = std::get<lmmc_real_t>(exp_num->value());
            if (std::isfinite(d) && d >= 0 && d == std::floor(d) && d < 1000.0) {
                e_val = static_cast<int>(d);
            }
        }

        if (e_val < 0) return -1;

        int base_deg = tf_degree_in(pow->base(), var);
        if (base_deg < 0) return -1;
        return base_deg * e_val;
    }

    /// FunctionNode 依赖 var 意味着非多项式结构
    if (std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return -1;
    }

    return 0;
}

/**
 * @brief 验证换元后的表达式确实是不定元的多项式。
 *
 * 检查表达式中是否仍残留超越函数依赖于任何不定元或原始变量。
 * 若表达式对所有候选变量的次数均可确定（非负整数），则视为有效多项式。
 *
 * @param[in] poly_expr      换元后的表达式
 * @param[in] all_variables  所有候选变量名（不定元 + 原始变量）
 * @return 表达式是有效多项式返回 true
 * @internal
 */
static bool tf_validate_polynomial(
    const std::shared_ptr<SymbolicExpr>& poly_expr,
    const std::vector<std::string>& all_variables) {

    if (!poly_expr || !lamina::detail::node(poly_expr)) return false;

    /// 对每个变量检查次数是否可确定
    for (const auto& var : all_variables) {
        if (expression_depends_on_variable(lamina::detail::node(poly_expr), var)) {
            int deg = tf_degree_in(lamina::detail::node(poly_expr), var);
            if (deg < 0) return false;
        }
    }
    return true;
}

/**
 * @brief 将换元后的表达式构造为有理系数多项式。
 *
 * 对换元后的多项式表达式，先展开再选择次数最高的不定元作为主变量，
 * 通过 symbolic_to_poly 转换为 Polynomial<Rational>。
 * 单变量情形（仅一个不定元且不含原始变量）直接转换；
 * 多变量情形选择次数最高的变量作为主变量。
 *
 * @param[in] poly_expr       换元后的符号表达式
 * @param[in] indeterminates  不定元名称列表（如 {"u0", "u1"}）
 * @param[in] original_var    原始目标变量名（如 "x"）
 * @return 多项式构造结果
 */
TfPolyBuildResult tf_build_polynomial(
    const std::shared_ptr<SymbolicExpr>& poly_expr,
    const std::vector<std::string>& indeterminates,
    const std::string& original_var) {

    TfPolyBuildResult result;
    result.success = false;

    if (!poly_expr || !lamina::detail::node(poly_expr)) {
        return result;
    }

    /// 先对原始表达式进行多项式验证（在展开之前）
    /// 这可以捕获分数指数等情形，避免 expand() 错误简化后遗漏
    {
        std::vector<std::string> pre_vars;
        for (const auto& ind : indeterminates) {
            if (expression_depends_on_variable(lamina::detail::node(poly_expr), ind)) {
                pre_vars.push_back(ind);
            }
        }
        if (expression_depends_on_variable(lamina::detail::node(poly_expr), original_var)) {
            pre_vars.push_back(original_var);
        }
        if (!pre_vars.empty() && !tf_validate_polynomial(poly_expr, pre_vars)) {
            return result;
        }
    }

    /// 展开表达式以确保多项式结构可见
    auto expanded = poly_expr->expand();
    if (!expanded || !lamina::detail::node(expanded)) {
        expanded = poly_expr;
    }

    /// 收集所有候选变量：不定元 + 原始变量（若表达式仍依赖它）
    std::vector<std::string> all_variables;
    for (const auto& ind : indeterminates) {
        if (expression_depends_on_variable(lamina::detail::node(expanded), ind)) {
            all_variables.push_back(ind);
        }
    }
    if (expression_depends_on_variable(lamina::detail::node(expanded), original_var)) {
        all_variables.push_back(original_var);
    }

    /// 无变量依赖：常数表达式
    if (all_variables.empty()) {
        Rational c = extract_coeff_value<Rational>(expanded);
        result.poly = Polynomial<Rational>({c}, "x");
        result.main_variable = "x";
        result.success = true;
        return result;
    }

    /// 验证表达式是否为有效多项式
    if (!tf_validate_polynomial(expanded, all_variables)) {
        return result;
    }

    /// 计算每个变量的次数，选择最高次的作为主变量
    std::string main_var;
    int max_degree = -1;

    for (const auto& var : all_variables) {
        int deg = tf_degree_in(lamina::detail::node(expanded), var);
        if (deg > max_degree) {
            max_degree = deg;
            main_var = var;
        }
    }

    if (main_var.empty() || max_degree < 0) {
        return result;
    }

    result.main_variable = main_var;

    /// 收集参数变量（非主变量）
    for (const auto& var : all_variables) {
        if (var != main_var) {
            result.param_variables.push_back(var);
        }
    }

    /// 单变量情形：仅主变量，无参数变量
    /// 此时 symbolic_to_poly<Rational> 可直接提取有理系数
    if (result.param_variables.empty()) {
        result.poly = symbolic_to_poly<Rational>(expanded, main_var);
        /// 验证转换结果非零（除非原表达式确实为零）
        if (result.poly.is_zero() && !expanded->is_zero()) {
            /// 转换失败：可能系数提取出错
            result.success = false;
            return result;
        }
        result.success = true;
        return result;
    }

    /// 多变量情形：主变量策略
    /// 使用 symbolic_to_poly<Rational> 以主变量转换。
    /// 注意：非主变量若出现在系数位置，extract_coeff_value<Rational> 会返回 0，
    /// 导致信息丢失。因此仅当参数变量不实际出现在系数中时才能成功。
    //
    /// 对于如 u0² - x² 这种情形，以 u0 为主变量时，
    /// 常数项为 -x²，无法表示为 Rational。
    /// 此时尝试：若所有参数变量的系数均为有理数（即参数变量仅以有理数倍出现），
    /// 则可以成功；否则标记为需要符号系数处理。
    result.poly = symbolic_to_poly<Rational>(expanded, main_var);

    /// 验证：重建多项式并与原表达式比较
    /// 简单验证：检查多项式次数是否与预期一致
    if (result.poly.degree() == max_degree) {
        result.success = true;
        return result;
    }

    /// 次数不匹配，说明系数提取有损失
    /// 尝试反转主变量选择：如果有其他变量能产生正确结果
    for (const auto& var : all_variables) {
        if (var == main_var) continue;

        int deg = tf_degree_in(lamina::detail::node(expanded), var);
        auto trial_poly = symbolic_to_poly<Rational>(expanded, var);
        if (trial_poly.degree() == deg && deg > 0) {
            result.main_variable = var;
            result.poly = trial_poly;
            result.param_variables.clear();
            for (const auto& v : all_variables) {
                if (v != var) result.param_variables.push_back(v);
            }
            result.success = true;
            return result;
        }
    }

    /// 所有尝试均失败：多变量情形无法用 Polynomial<Rational> 精确表示
    result.success = false;
    return result;
}


/**
 * @brief 计算多项式的无平方因子部分。
 *
 * 算法：
 *   1. 计算 f' = derivative(f)
 *   2. g = gcd(f, f')
 *   3. 若 deg(g) == 0，f 已为 square-free，直接返回
 *   4. 否则 square_free_part = f / g
 *
 * 用于 Berlekamp 模分解前的预处理，确保输入多项式无重因子。
 *
 * @param[in] poly 输入的有理系数多项式
 * @return 无平方因子预处理结果
 */
TfSquareFreeResult tf_square_free(const Polynomial<Rational>& poly) {
    TfSquareFreeResult result;

    /// 常数或零多项式：本身即为 square-free
    if (poly.degree() <= 0) {
        result.square_free = poly;
        result.repeated_factor = Polynomial<Rational>(Rational(1), poly.variable_name);
        result.had_repeated_factors = false;
        return result;
    }

    /// 线性多项式：始终 square-free
    if (poly.degree() == 1) {
        result.square_free = poly;
        result.repeated_factor = Polynomial<Rational>(Rational(1), poly.variable_name);
        result.had_repeated_factors = false;
        return result;
    }

    /// 计算形式导数
    Polynomial<Rational> deriv = poly.differentiate();

    /// 计算 gcd(f, f')
    Polynomial<Rational> g = Polynomial<Rational>::gcd(poly, deriv);

    /// 若 gcd 为常数（degree 0），f 已为 square-free
    if (g.degree() <= 0) {
        result.square_free = poly.make_monic();
        result.repeated_factor = Polynomial<Rational>(Rational(1), poly.variable_name);
        result.had_repeated_factors = false;
        return result;
    }

    /// 存在重复因子：square_free_part = f / gcd(f, f')
    auto [quotient, remainder] = poly.div_mod(g);
    result.square_free = quotient.make_monic();
    result.repeated_factor = g.make_monic();
    result.had_repeated_factors = true;
    return result;
}


/**
 * @brief 计算整数的 popcount（二进制中 1 的个数）。
 *
 * @param[in] mask 待计算的无符号整数
 * @return 二进制表示中 1 的个数
 * @internal
 */
/**
 * @brief 混合超越表达式不可约因式分解主入口。
 *
 * 对含超越函数（sin、cos、exp、ln、tan）的表达式执行因式分解。
 * 完整流程：毕达哥拉斯化简 → 换元检测 → 多项式构造 → 无平方因子预处理
 * → Berlekamp 模分解 → Hensel 提升 → Zassenhaus 因子组合 → 逆换元 → 化简。
 *
 * 若任一阶段失败（如表达式不可多项式化、无合适素数等），
 * 返回原表达式作为单一不可约因子。
 *
 * @param[in] expr 待分解的符号表达式
 * @param[in] var  目标变量名
 * @return 不可约因子列表（乘积等于原表达式，可能含常数因子）
 */
std::vector<std::shared_ptr<SymbolicExpr>> factor_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    if (!expr || !lamina::detail::node(expr)) return {};

    /// --- 预处理：毕达哥拉斯恒等式化简 ---
    /// 将 sin²(f) + cos²(f) 子模式替换为 1，简化后续分解
    auto pyth_simplified = tf_simplify_pythagorean(expr, var);

    /// --- 快速路径：乘积结构检测 ---
    /// 若表达式已为独立子表达式的乘积，直接返回各因子
    auto mult_factors = tf_detect_multiplicative_structure(pyth_simplified);
    if (!mult_factors.empty()) {
        return mult_factors;
    }

    /// --- 快速路径：指数分离 ---
    /// 若表达式为加法且所有项共享公共 exp(f(x)) 因子，提取之
    auto exp_factors = tf_detect_exponential_separation(pyth_simplified, var);
    if (!exp_factors.empty()) {
        return exp_factors;
    }

    /// 若表达式不含超越函数，直接返回原表达式
    if (!tf_contains_transcendental(lamina::detail::node(pyth_simplified), var)) {
        return {pyth_simplified};
    }

    /// --- Phase 1: 换元检测 ---
    TransSubstitutionResult sub_result = detect_trans_substitutions(pyth_simplified, var);

    /// 无换元映射：表达式不含可处理的超越子表达式
    if (sub_result.mappings.empty()) {
        return {pyth_simplified};
    }

    /// --- 快速路径：线性不可约检测 ---
    /// 若换元后表达式对所有不定元和原始变量均为线性，则不可约
    if (tf_is_linear_irreducible(sub_result, var)) {
        return {pyth_simplified};
    }

    /// --- Phase 2: 多项式构造 ---
    std::vector<std::string> indeterminates;
    indeterminates.reserve(sub_result.mappings.size());
    for (const auto& m : sub_result.mappings) {
        indeterminates.push_back(m.indeterminate);
    }

    TfPolyBuildResult poly_result = tf_build_polynomial(
        sub_result.poly_expr, indeterminates, var);

    if (!poly_result.success) {
        /// 多项式构造失败：表达式不可多项式化
        return {pyth_simplified};
    }

    /// 常数或线性多项式：不可约
    if (poly_result.poly.degree() <= 1) {
        return {pyth_simplified};
    }

    /// --- Phase 2.3: 无平方因子预处理 ---
    TfSquareFreeResult sqf_result = tf_square_free(poly_result.poly);

    /// 使用 square-free 部分进行后续分解
    Polynomial<Rational> work_poly = sqf_result.square_free;

    /// 若 square-free 部分为常数或线性，不可约
    if (work_poly.degree() <= 1) {
        return {pyth_simplified};
    }

    /// --- Phase 3: Berlekamp 模分解 ---
    /// 尝试多个素数以找到合适的分解
    static const int64_t TRIAL_PRIMES[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
    BerlekampResult berl_result;
    bool berlekamp_success = false;

    for (int64_t p : TRIAL_PRIMES) {
        try {
            berl_result = berlekamp_factor(work_poly, p);
            if (berl_result.prime > 0 && !berl_result.factors.empty()) {
                berlekamp_success = true;
                break;
            }
        } catch (...) {
            continue;
        }
    }

    if (!berlekamp_success || berl_result.factors.size() <= 1) {
        /// Berlekamp 分解失败或多项式在所有尝试素数下不可约
        return {pyth_simplified};
    }

    /// --- Phase 4: Hensel 提升 ---
    /// 将有理系数多项式转为整系数（乘以分母 LCM）
    BigInt lcm_denom(1);
    for (const auto& c : work_poly.coeffs) {
        BigInt d = c.get_denominator();
        lcm_denom = BigInt::lcm(lcm_denom, d);
    }

    std::vector<BigInt> int_coeffs;
    int_coeffs.reserve(work_poly.coeffs.size());
    for (const auto& c : work_poly.coeffs) {
        BigInt num = c.get_numerator() * (lcm_denom / c.get_denominator());
        int_coeffs.push_back(num);
    }
    Polynomial<BigInt> int_poly(int_coeffs, work_poly.variable_name);

    /// 计算提升界：使用 Mignotte bound 确定需要的精度
    int n = work_poly.degree();
    BigInt max_coeff(0);
    for (const auto& c : int_coeffs) {
        BigInt ac = c.Abs();
        if (ac > max_coeff) max_coeff = ac;
    }
    /// 保守 Mignotte 界：2^n * max_coeff
    BigInt mignotte = max_coeff;
    for (int i = 0; i < n; ++i) {
        mignotte = mignotte * BigInt(2);
    }

    /// 确定提升次数 k 使得 p^k > 2 * mignotte
    BigInt target = mignotte * BigInt(2);
    int64_t prime = berl_result.prime;
    int lift_bound = 1;
    BigInt pk(static_cast<long long>(prime));
    while (pk <= target && lift_bound < 100) {
        pk = pk * BigInt(static_cast<long long>(prime));
        lift_bound++;
    }

    std::vector<Polynomial<BigInt>> lifted_factors;
    try {
        lifted_factors = hensel_lift(int_poly, berl_result.factors, prime, lift_bound);
    } catch (...) {
        return {pyth_simplified};
    }

    if (lifted_factors.empty()) {
        return {pyth_simplified};
    }

    /// --- Phase 5: Zassenhaus 因子组合 ---
    int64_t prime_power = 1;
    for (int i = 0; i < lift_bound; ++i) {
        /// 防止溢出：若超出 int64_t 范围则截断
        if (prime_power > std::numeric_limits<int64_t>::max() / prime) {
            prime_power = std::numeric_limits<int64_t>::max();
            break;
        }
        prime_power *= prime;
    }

    std::vector<Polynomial<Rational>> true_factors;
    try {
        true_factors = zassenhaus_combine(work_poly, lifted_factors, prime_power);
    } catch (...) {
        return {pyth_simplified};
    }

    /// 若组合后仅得到一个因子（或无因子），表达式不可约
    if (true_factors.size() <= 1) {
        return {pyth_simplified};
    }

    /// --- Phase 6: 逆换元 ---
    std::vector<std::shared_ptr<SymbolicExpr>> symbolic_factors;
    symbolic_factors.reserve(true_factors.size());

    for (const auto& poly_factor : true_factors) {
        /// 将多项式因子转回符号表达式
        auto sym_factor = poly_to_symbolic(poly_factor);
        if (!sym_factor) continue;

        /// 逆换元：将不定元替换回原始超越表达式
        auto back_sub = tf_back_substitute(sym_factor, sub_result.mappings);
        if (back_sub) {
            symbolic_factors.push_back(back_sub);
        }
    }

    /// 若逆换元后因子数 ≤ 1，分解无效
    if (symbolic_factors.size() <= 1) {
        return {pyth_simplified};
    }

    /// --- Phase 6.2: 化简与常数提取 ---
    auto final_factors = tf_simplify_factors(symbolic_factors);

    /// 若化简后仅剩一个非常数因子，分解无效
    size_t non_const_count = 0;
    for (const auto& f : final_factors) {
        if (f && !f->is_number()) non_const_count++;
    }
    if (non_const_count <= 1 && final_factors.size() <= 1) {
        return {pyth_simplified};
    }

    /// 若存在重复因子（来自 square-free 预处理），需要恢复重数
    if (sqf_result.had_repeated_factors) {
        /// 将重复因子也进行逆换元并加入结果
        auto rep_sym = poly_to_symbolic(sqf_result.repeated_factor);
        if (rep_sym && !rep_sym->is_one()) {
            auto rep_back = tf_back_substitute(rep_sym, sub_result.mappings);
            if (rep_back && !rep_back->is_one()) {
                final_factors.push_back(rep_back);
            }
        }
    }

    return final_factors;
}

} // namespace lamina
