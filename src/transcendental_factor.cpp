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
#include "poly_utils_internal.hpp"

#include <string>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <limits>

namespace lamina {

// ============================================================
/// 文件局部辅助函数
// ============================================================

/**
 * @brief 判断 FunctionNode 的函数类型是否为需要换元的超越函数。
 * @param[in] ft 函数类型枚举
 * @return 属于 Sin/Cos/Exp/Ln/Tan 之一返回 true
 * @internal
 */
static bool tf_is_transcendental_type(FunctionNode::FuncType ft) {
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
            depends_on_var(func->arguments()[0], var)) {
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

// ============================================================
/// 公共 API 实现
// ============================================================

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

// ============================================================
/// Phase 2: 多项式构造
// ============================================================

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
static int tf_degree_in(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    if (!node) return 0;

    if (!depends_on_var(node, var)) return 0;

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
        if (depends_on_var(lamina::detail::node(poly_expr), var)) {
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
            if (depends_on_var(lamina::detail::node(poly_expr), ind)) {
                pre_vars.push_back(ind);
            }
        }
        if (depends_on_var(lamina::detail::node(poly_expr), original_var)) {
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
        if (depends_on_var(lamina::detail::node(expanded), ind)) {
            all_variables.push_back(ind);
        }
    }
    if (depends_on_var(lamina::detail::node(expanded), original_var)) {
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

// ============================================================
/// Phase 2.3: 无平方因子预处理
// ============================================================

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

// ============================================================
/// Phase 5: Zassenhaus 因子组合
// ============================================================

/**
 * @brief 计算整数的 popcount（二进制中 1 的个数）。
 *
 * @param[in] mask 待计算的无符号整数
 * @return 二进制表示中 1 的个数
 * @internal
 */
static int zc_popcount(uint64_t mask) {
    int count = 0;
    while (mask) {
        count += static_cast<int>(mask & 1);
        mask >>= 1;
    }
    return count;
}

/**
 * @brief 将 BigInt 系数归约到对称表示 [-m/2, m/2)。
 *
 * @param[in] c 待归约的系数
 * @param[in] m 模数（正整数）
 * @return 对称表示下的归约值
 * @internal
 */
static BigInt zc_symmetric_mod(const BigInt& c, const BigInt& m) {
    if (m.is_zero()) return c;

    BigInt r = c % m;
    if (r.IsNegative()) {
        r = r + m;
    }

    BigInt half_m = m / BigInt(2);
    if (r > half_m) {
        r = r - m;
    }
    return r;
}

/**
 * @brief 计算一组 Hensel 提升因子的乘积，系数模 prime_power 归约。
 *
 * 对给定的因子子集（由位掩码指定），逐个相乘并对每个系数取模归约到对称表示。
 *
 * @param[in] lifted_factors 所有提升后的因子
 * @param[in] mask           子集位掩码（第 i 位为 1 表示选取第 i 个因子）
 * @param[in] mod            模数 p^k
 * @return 子集因子乘积的系数向量（对称表示）
 * @internal
 */
static std::vector<BigInt> zc_subset_product(
    const std::vector<Polynomial<BigInt>>& lifted_factors,
    uint64_t mask,
    const BigInt& mod) {

    std::vector<BigInt> product = {BigInt(1)};

    for (size_t i = 0; i < lifted_factors.size(); ++i) {
        if (!((mask >> i) & 1)) continue;

        const auto& factor_coeffs = lifted_factors[i].coeffs;
        if (factor_coeffs.empty()) return {};

        /// 多项式乘法
        size_t new_size = product.size() + factor_coeffs.size() - 1;
        std::vector<BigInt> new_product(new_size, BigInt(0));

        for (size_t a = 0; a < product.size(); ++a) {
            if (product[a].is_zero()) continue;
            for (size_t b = 0; b < factor_coeffs.size(); ++b) {
                if (factor_coeffs[b].is_zero()) continue;
                new_product[a + b] = new_product[a + b] + product[a] * factor_coeffs[b];
            }
        }

        /// 对称归约
        for (auto& c : new_product) {
            c = zc_symmetric_mod(c, mod);
        }

        /// 去除高次零系数
        while (!new_product.empty() && new_product.back().is_zero()) {
            new_product.pop_back();
        }

        product = std::move(new_product);
    }

    return product;
}

/**
 * @brief 将 BigInt 系数向量转换为有理系数多项式。
 *
 * 每个 BigInt 系数直接转换为 Rational（分母为 1）。
 *
 * @param[in] coeffs BigInt 系数向量
 * @param[in] var    变量名
 * @return 有理系数多项式
 * @internal
 */
static Polynomial<Rational> zc_bigint_to_rational_poly(
    const std::vector<BigInt>& coeffs,
    const std::string& var) {

    std::vector<Rational> rat_coeffs;
    rat_coeffs.reserve(coeffs.size());
    for (const auto& c : coeffs) {
        rat_coeffs.emplace_back(c);
    }
    return Polynomial<Rational>(rat_coeffs, var);
}

/**
 * @brief 对 BigInt 系数执行有理数重构。
 *
 * 给定整数 a 和模数 m，使用扩展欧几里得算法寻找有理数 p/q 满足：
 * - a ≡ p/q (mod m)
 * - |p| ≤ floor(sqrt(m/2)), |q| ≤ floor(sqrt(m/2))
 *
 * 当模数超出 int64_t 范围时使用 BigInt 算术；否则委托给
 * modular_arithmetic.hpp 中的 int64_t 版本以获得更好性能。
 *
 * @param[in] a 待重构的整数系数（已归约到对称表示）
 * @param[in] m 模数（正整数，通常为 p^k）
 * @param[out] num 输出分子
 * @param[out] den 输出分母
 * @return 重构成功返回 true；失败（无法满足界约束或 gcd≠1）返回 false
 * @internal
 */
static bool zc_rational_reconstruction(
    const BigInt& a,
    const BigInt& m,
    BigInt& num,
    BigInt& den) {

    if (m.is_zero() || m.IsNegative()) return false;

    /// 对于可精确表示的小模数，委托给 int64_t 版本
    const auto a_small = a.try_to_int64();
    const auto m_small = m.try_to_int64();
    if (a_small && m_small) {
        const int64_t a_val = *a_small;
        const int64_t m_val = *m_small;

        auto [p, q] = lamina::rational_reconstruction(a_val, m_val);
        if (q == 0) return false;

        num = BigInt(static_cast<long long>(p));
        den = BigInt(static_cast<long long>(q));
        return true;
    }

    /// BigInt 版本的有理重构
    /// 将 a 归约到 [0, m)
    BigInt a_mod = a % m;
    if (a_mod.IsNegative()) {
        a_mod = a_mod + m;
    }

    /// 计算界 bound = floor(sqrt(m/2))
    BigInt half_m = m / BigInt(2);
    BigInt bound = half_m.sqrt();
    if (bound.is_zero()) bound = BigInt(1);

    /// 扩展欧几里得算法
    BigInt r0 = m, r1 = a_mod;
    BigInt s0 = BigInt(0), s1 = BigInt(1);

    while (r1 > bound) {
        BigInt q = r0 / r1;
        BigInt r_new = r0 - q * r1;
        BigInt s_new = s0 - q * s1;

        r0 = r1;
        r1 = r_new;
        s0 = s1;
        s1 = s_new;
    }

    BigInt p = r1;
    BigInt q_val = s1;

    /// 确保分母为正
    if (q_val.IsNegative()) {
        p = -p;
        q_val = -q_val;
    }

    /// 验证界约束
    if (p.Abs() > bound || q_val.is_zero() || q_val > bound) {
        return false;
    }

    /// 验证 gcd(|p|, q) == 1
    BigInt g = BigInt::gcd(p.Abs(), q_val);
    if (g != BigInt(1)) {
        return false;
    }

    num = p;
    den = q_val;
    return true;
}

/**
 * @brief 对子集乘积的所有系数执行有理重构，构造有理系数多项式。
 *
 * 对 product_coeffs 中的每个 BigInt 系数调用 zc_rational_reconstruction，
 * 若所有系数均成功重构，则返回对应的有理系数多项式。
 * 若任一系数重构失败，则回退到整数系数直接转换。
 *
 * @param[in] product_coeffs 子集乘积的 BigInt 系数向量（对称表示）
 * @param[in] mod            模数 p^k
 * @param[in] var            变量名
 * @return 有理系数多项式（通过有理重构或整数直接转换）
 * @internal
 */
static Polynomial<Rational> zc_reconstruct_candidate(
    const std::vector<BigInt>& product_coeffs,
    const BigInt& mod,
    const std::string& var) {

    std::vector<Rational> rat_coeffs;
    rat_coeffs.reserve(product_coeffs.size());

    /// 尝试对每个系数执行有理重构
    for (const auto& c : product_coeffs) {
        BigInt num, den;
        if (zc_rational_reconstruction(c, mod, num, den)) {
            rat_coeffs.emplace_back(num, den);
        } else {
            /// 任一系数重构失败，回退到整数系数
            return zc_bigint_to_rational_poly(product_coeffs, var);
        }
    }

    return Polynomial<Rational>(rat_coeffs, var);
}

/**
 * @brief 检验候选因子是否整除原多项式（精确有理除法）。
 *
 * 执行多项式带余除法 f / candidate，若余数为零则整除。
 *
 * @param[in] f         原多项式
 * @param[in] candidate 候选因子
 * @return 整除返回 true
 * @internal
 */
static bool zc_divides_exactly(
    const Polynomial<Rational>& f,
    const Polynomial<Rational>& candidate) {

    if (candidate.is_zero()) return false;
    if (candidate.degree() > f.degree()) return false;
    if (candidate.degree() == 0) return true;

    auto [quotient, remainder] = f.div_mod(candidate);
    return remainder.is_zero();
}

/**
 * @brief 使候选多项式首一化（首项系数归一）。
 *
 * @param[in] poly 输入多项式
 * @return 首一多项式
 * @internal
 */
static Polynomial<Rational> zc_make_primitive(const Polynomial<Rational>& poly) {
    if (poly.is_zero()) return poly;
    return poly.make_monic();
}

/**
 * @brief 计算多项式系数的 L1 范数（系数绝对值之和）。
 *
 * @param[in] coeffs BigInt 系数向量
 * @return L1 范数
 * @internal
 */
static BigInt zc_l1_norm(const std::vector<BigInt>& coeffs) {
    BigInt norm(0);
    for (const auto& c : coeffs) {
        norm = norm + c.Abs();
    }
    return norm;
}

/**
 * @brief 计算 Mignotte 界：用于判断候选因子系数的合理范围。
 *
 * Mignotte 界为 C(n, floor(n/2)) * ||f||_2 的近似上界。
 * 此处使用简化估计：2^n * max|coeff(f)| 作为保守上界。
 *
 * @param[in] poly 原多项式
 * @return Mignotte 界的近似值
 * @internal
 */
static BigInt zc_mignotte_bound(const Polynomial<Rational>& poly) {
    BigInt max_coeff(0);
    for (const auto& c : poly.coeffs) {
        BigInt abs_num = c.get_numerator().Abs();
        /// 取 |numerator| 作为系数大小的近似上界
        if (abs_num > max_coeff) {
            max_coeff = abs_num;
        }
    }

    int n = poly.degree();
    if (n <= 0) return max_coeff;

    /// 2^n * max_coeff 作为保守 Mignotte 界
    BigInt bound = max_coeff;
    for (int i = 0; i < n; ++i) {
        bound = bound * BigInt(2);
    }
    return bound;
}

/**
 * @brief 基于度数和范数剪枝的因子组合（用于因子数 > 15 的情形）。
 *
 * 当模因子数超过 15 时，穷举所有 2^r 个子集不再实际可行。
 * 本函数使用以下启发式剪枝策略加速搜索：
 *
 * 1. 度数剪枝：预计算每个提升因子的次数，仅枚举总度数为原多项式
 *    次数的合理因子度数（≤ deg(f)/2）的子集。
 * 2. 范数剪枝：对候选子集乘积计算 L1 范数，若超过 2 倍 Mignotte 界
 *    则立即拒绝，避免昂贵的有理重构和整除性检验。
 *
 * @note 完整的 LLL 格基约化算法可进一步优化此步骤，将搜索复杂度
 *       从指数级降至多项式级。对于典型用例（15 < r ≤ 30），
 *       度数+范数剪枝已足够高效。
 *
 * @param[in] poly           有理系数原多项式
 * @param[in] lifted_factors Hensel 提升后的整系数因子
 * @param[in] mod            模数 p^k
 * @return 有理系数真因子列表
 *
 * @see van Hoeij, M. "Factoring polynomials and the knapsack problem."
 *      Journal of Number Theory, 95(2), 2002.
 * @internal
 */
static std::vector<Polynomial<Rational>> zc_lll_pruned_combine(
    const Polynomial<Rational>& poly,
    const std::vector<Polynomial<BigInt>>& lifted_factors,
    const BigInt& mod) {

    std::vector<Polynomial<Rational>> true_factors;
    size_t r = lifted_factors.size();
    std::string var = poly.variable_name;

    /// 预计算每个提升因子的次数
    std::vector<int> factor_degrees(r);
    for (size_t i = 0; i < r; ++i) {
        factor_degrees[i] = static_cast<int>(lifted_factors[i].coeffs.size()) - 1;
        if (factor_degrees[i] < 0) factor_degrees[i] = 0;
    }

    /// 计算 Mignotte 界用于范数剪枝
    BigInt mignotte = zc_mignotte_bound(poly);
    BigInt norm_threshold = mignotte * BigInt(2);

    /// 工作副本
    Polynomial<Rational> remaining = poly;
    uint64_t active_mask = (1ULL << r) - 1;

    /// TODO: 完整的 LLL 格基约化实现可将此搜索从指数级优化为多项式级。
    /// 当前使用度数+范数启发式剪枝，对典型用例（15 < r ≤ 30）已足够高效。
    /// 参考：van Hoeij (2002) 的 LLL-based 因子组合算法。

    bool found_factor = true;
    while (found_factor) {
        found_factor = false;

        int active_count = zc_popcount(active_mask);
        if (active_count <= 1) break;
        if (remaining.degree() <= 1) break;

        int remaining_deg = remaining.degree();
        int max_subset_size = active_count / 2;

        for (int subset_size = 1; subset_size <= max_subset_size; ++subset_size) {
            uint64_t subset = (1ULL << subset_size) - 1;

            while (subset <= active_mask) {
                if ((subset & active_mask) == subset && zc_popcount(subset) == subset_size) {

                    /// --- 度数剪枝 ---
                    /// 计算子集因子的总度数
                    int subset_degree = 0;
                    for (size_t i = 0; i < r; ++i) {
                        if ((subset >> i) & 1) {
                            subset_degree += factor_degrees[i];
                        }
                    }

                    /// 跳过总度数超过剩余多项式度数一半的子集
                    if (subset_degree > remaining_deg / 2) {
                        goto next_pruned_subset;
                    }

                    /// 跳过总度数为 0 的子集（不可能产生有意义的因子）
                    if (subset_degree <= 0) {
                        goto next_pruned_subset;
                    }

                    {
                        /// 计算子集因子乘积 mod p^k
                        std::vector<BigInt> product_coeffs = zc_subset_product(
                            lifted_factors, subset, mod);

                        if (product_coeffs.empty()) {
                            goto next_pruned_subset;
                        }

                        /// --- 范数剪枝 ---
                        /// 若乘积的 L1 范数超过 2 倍 Mignotte 界，拒绝此候选
                        BigInt candidate_norm = zc_l1_norm(product_coeffs);
                        if (candidate_norm > norm_threshold) {
                            goto next_pruned_subset;
                        }

                        {
                            /// 有理重构
                            Polynomial<Rational> candidate = zc_reconstruct_candidate(
                                product_coeffs, mod, var);

                            candidate = zc_make_primitive(candidate);

                            /// 整除性检验
                            if (!candidate.is_zero() && candidate.degree() > 0 &&
                                candidate.degree() < remaining.degree() &&
                                zc_divides_exactly(remaining, candidate)) {

                                true_factors.push_back(candidate);

                                auto [quotient, rem] = remaining.div_mod(candidate);
                                remaining = quotient;

                                active_mask &= ~subset;
                                found_factor = true;
                                goto restart_pruned_enumeration;
                            }
                        }
                    }
                }

                next_pruned_subset:
                if (subset == 0) break;
                uint64_t c = subset & (-static_cast<int64_t>(subset));
                uint64_t rr = subset + c;
                subset = (((rr ^ subset) >> 2) / c) | rr;
                if (subset > active_mask || subset == 0) break;
            }
        }

        restart_pruned_enumeration:;
    }

    /// 剩余因子形成最后一个真因子
    if (!remaining.is_zero() && remaining.degree() > 0) {
        true_factors.push_back(remaining.make_monic());
    } else if (!remaining.is_zero() && remaining.degree() == 0) {
        if (remaining.coeffs[0] != Rational(1)) {
            true_factors.push_back(remaining);
        }
    }

    return true_factors;
}

/**
 * @brief Zassenhaus 因子组合：从 Hensel 提升后的因子中筛选真因子。
 *
 * 算法：
 * 1. 将原多项式转为首一形式
 * 2. 按子集大小从 1 到 floor(r/2) 枚举 lifted_factors 的子集
 * 3. 对每个子集，计算其乘积 mod p^k（对称表示）
 * 4. 对乘积系数执行有理重构（rational_reconstruction），恢复有理系数
 * 5. 首一化后检验候选因子是否整除当前剩余多项式
 * 6. 若整除，记录为真因子，从因子池中移除已用因子，更新剩余多项式
 * 7. 剩余因子形成最后一个真因子
 *
 * 当因子数超过 15 时，使用度数和范数启发式剪枝加速搜索。
 *
 * @param[in] poly           有理系数原多项式
 * @param[in] lifted_factors Hensel 提升后的整系数因子
 * @param[in] prime_power    素数幂 p^k
 * @return 有理系数真因子列表
 *
 * @see Zassenhaus, H. "On Hensel factorization, I."
 *      Journal of Number Theory, 1(3), 1969.
 */
std::vector<Polynomial<Rational>> zassenhaus_combine(
    const Polynomial<Rational>& poly,
    const std::vector<Polynomial<BigInt>>& lifted_factors,
    int64_t prime_power) {

    std::vector<Polynomial<Rational>> true_factors;

    /// 边界情形
    if (poly.is_zero() || lifted_factors.empty()) {
        if (!poly.is_zero()) true_factors.push_back(poly);
        return true_factors;
    }

    /// 单因子：原多项式本身不可约
    if (lifted_factors.size() == 1) {
        true_factors.push_back(poly);
        return true_factors;
    }

    /// 因子数上限
    size_t r = lifted_factors.size();
    if (r > 30) {
        /// 因子数过多（> 30），即使剪枝也不实际，返回原多项式
        true_factors.push_back(poly);
        return true_factors;
    }

    /// 当因子数 > 15 时，使用度数+范数剪枝的组合策略
    if (r > 15) {
        BigInt mod(prime_power);
        return zc_lll_pruned_combine(poly, lifted_factors, mod);
    }

    BigInt mod(prime_power);
    std::string var = poly.variable_name;

    /// 工作副本：当前剩余多项式和因子池
    Polynomial<Rational> remaining = poly;

    /// 活跃因子索引集合（用位掩码表示）
    uint64_t active_mask = (1ULL << r) - 1;  // 所有因子初始活跃

    /// 按子集大小从 1 到 floor(active_count/2) 枚举
    bool found_factor = true;
    while (found_factor) {
        found_factor = false;

        int active_count = zc_popcount(active_mask);

        /// 早期终止：仅剩 1 个活跃因子，剩余多项式本身即为不可约因子
        if (active_count <= 1) break;

        /// 早期终止：剩余多项式为线性或常数，必然不可约
        if (remaining.degree() <= 1) break;

        int max_subset_size = active_count / 2;

        for (int subset_size = 1; subset_size <= max_subset_size; ++subset_size) {
            /// 枚举所有大小为 subset_size 的活跃因子子集
            /// 使用位掩码枚举：遍历 active_mask 的所有子集中 popcount == subset_size 的
            uint64_t subset = 0;

            /// Gosper's hack 初始化：最小的 subset_size 位子集
            subset = (1ULL << subset_size) - 1;

            while (subset <= active_mask) {
                /// 检查 subset 是否为 active_mask 的子集
                if ((subset & active_mask) == subset && zc_popcount(subset) == subset_size) {
                    /// 计算子集因子的乘积 mod p^k
                    std::vector<BigInt> product_coeffs = zc_subset_product(
                        lifted_factors, subset, mod);

                    if (product_coeffs.empty()) {
                        /// 跳过空乘积
                        goto next_subset;
                    }

                    {
                        /// 通过有理重构构造候选因子
                        Polynomial<Rational> candidate = zc_reconstruct_candidate(
                            product_coeffs, mod, var);

                        /// 首一化
                        candidate = zc_make_primitive(candidate);

                        /// 检验整除性
                        if (!candidate.is_zero() && candidate.degree() > 0 &&
                            candidate.degree() < remaining.degree() &&
                            zc_divides_exactly(remaining, candidate)) {

                            /// 找到真因子
                            true_factors.push_back(candidate);

                            /// 更新剩余多项式
                            auto [quotient, rem] = remaining.div_mod(candidate);
                            remaining = quotient;

                            /// 从活跃集中移除已用因子
                            active_mask &= ~subset;

                            found_factor = true;
                            goto restart_enumeration;
                        }
                    }
                }

                next_subset:
                /// 下一个子集（Gosper's hack）
                if (subset == 0) break;
                uint64_t c = subset & (-static_cast<int64_t>(subset));
                uint64_t rr = subset + c;
                subset = (((rr ^ subset) >> 2) / c) | rr;

                /// 防止溢出
                if (subset > active_mask || subset == 0) break;
            }
        }

        restart_enumeration:;
    }

    /// 剩余因子形成最后一个真因子
    if (!remaining.is_zero() && remaining.degree() > 0) {
        true_factors.push_back(remaining.make_monic());
    } else if (!remaining.is_zero() && remaining.degree() == 0) {
        /// 常数因子：若非 1 则记录
        if (remaining.coeffs[0] != Rational(1)) {
            true_factors.push_back(remaining);
        }
    }

    return true_factors;
}

// ============================================================
/// Phase 6: 逆换元（Back-substitution）
// ============================================================

/**
 * @brief 将因子表达式中的不定元变量替换回原始超越子表达式。
 *
 * 对换元后的因子执行逆操作：遍历映射列表，将每个不定元（u0, u1, ...）
 * 替换为其对应的原始超越表达式（如 sin(x)、cos(x) 等）。
 * 利用 SymbolicExpr::substitute() 逐一执行变量替换。
 *
 * @param[in] factor_expr 以不定元表示的因子表达式
 * @param[in] mappings    换元映射列表（indeterminate → trans_expr）
 * @return 替换后的符号表达式，以原始变量和超越函数表示
 * @internal
 */
static std::shared_ptr<SymbolicExpr> tf_back_substitute(
    const std::shared_ptr<SymbolicExpr>& factor_expr,
    const std::vector<TransSubstitution>& mappings) {

    if (!factor_expr || !lamina::detail::node(factor_expr)) return factor_expr;
    if (mappings.empty()) return factor_expr;

    auto result = factor_expr;

    for (const auto& m : mappings) {
        if (!m.trans_expr || m.indeterminate.empty()) continue;

        /// 仅当表达式依赖该不定元时才执行替换
        if (depends_on_var(lamina::detail::node(result), m.indeterminate)) {
            result = result->substitute(m.indeterminate, m.trans_expr);
        }
    }

    return result;
}

// ============================================================
/// Phase 6.2: 因子化简与常数提取
// ============================================================

/**
 * @brief 从 NumberNode 中提取有理数值。
 *
 * 将 BigInt、Rational、lmmc_real_t 统一转换为 Rational 表示。
 * 对于浮点数，仅当其为精确整数时才转换；否则返回失败。
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
 * @brief 对逆换元后的因子列表执行化简与常数乘子提取。
 *
 * 算法：
 * 1. 对每个因子调用 simplify() 进行规范化
 * 2. 若化简后的因子为纯数值（NumberNode），将其累乘到常数积中
 * 3. 若化简后的因子为乘积形式（MultiplyNode）且含数值前导系数，
 *    提取该系数并保留非常数部分
 * 4. 若累积常数 ≠ 1，将其作为独立数值因子插入结果列表首位
 *
 * @param[in,out] factors 因子列表，就地修改为化简后的结果
 * @return 化简并提取常数后的因子列表
 */
std::vector<std::shared_ptr<SymbolicExpr>> tf_simplify_factors(
    std::vector<std::shared_ptr<SymbolicExpr>>& factors) {

    std::vector<std::shared_ptr<SymbolicExpr>> result;
    Rational constant_product(1);

    for (auto& factor : factors) {
        if (!factor || !lamina::detail::node(factor)) continue;

        /// 调用 simplify() 规范化因子
        auto simplified = factor->simplify();
        if (!simplified || !lamina::detail::node(simplified)) {
            simplified = factor;
        }

        /// 情形 1：因子为纯数值
        if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(simplified))) {
            Rational val;
            if (tf_extract_rational(num, val)) {
                if (val != Rational(0)) {
                    constant_product = constant_product * val;
                }
                /// 零因子不累乘，但保留（整个乘积为零）
                else {
                    result.clear();
                    result.push_back(SymbolicExpr::number(0));
                    return result;
                }
            } else {
                /// 无法提取为有理数的数值（如非整数浮点），保留原样
                result.push_back(simplified);
            }
            continue;
        }

        /// 情形 2：因子为乘积形式，检查是否含数值前导系数
        if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(simplified))) {
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
                        /// 无法提取的数值操作数保留在非常数部分
                        non_numeric_ops.push_back(nop);
                    }
                }

                /// 重建非常数部分
                if (non_numeric_ops.size() == 1) {
                    result.push_back(lamina::detail::make_expression_ptr(non_numeric_ops[0]));
                } else {
                    result.push_back(lamina::detail::make_expression_ptr(
                        lamina::detail::make_node<MultiplyNode>(std::move(non_numeric_ops))));
                }
            } else if (numeric_ops.empty()) {
                /// 无数值操作数，保留原因子
                result.push_back(simplified);
            } else {
                /// 全为数值操作数：整个因子为常数
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

        /// 情形 3：非数值、非乘积形式，直接保留
        result.push_back(simplified);
    }

    /// 若累积常数 ≠ 1，插入为首个因子
    if (constant_product != Rational(1)) {
        auto const_expr = SymbolicExpr::number(constant_product);
        result.insert(result.begin(), const_expr);
    }

    return result;
}

// ============================================================
/// 特殊情形快速路径
// ============================================================

/**
 * @brief 判断换元后的表达式是否对所有不定元和原始变量均为线性。
 *
 * 若换元后的表达式对每个不定元（u0, u1, ...）以及原始变量（如 x）的次数
 * 均不超过 1，则该表达式在超越多项式环中不可约——无法进一步因式分解。
 * 典型例子：a*sin(x) + b*x + c 换元后为 a*u0 + b*x + c，对 u0 和 x 均为线性。
 *
 * @param[in] sub_result 换元结果（含 poly_expr 和 mappings）
 * @param[in] var        原始目标变量名
 * @return 表达式对所有变量均为线性返回 true
 * @internal
 */
static bool tf_is_linear_irreducible(
    const TransSubstitutionResult& sub_result,
    const std::string& var) {

    if (!sub_result.poly_expr || !lamina::detail::node(sub_result.poly_expr)) return false;
    if (sub_result.mappings.empty()) return false;

    const auto& root = lamina::detail::node(sub_result.poly_expr);

    /// 检查每个不定元的次数是否 ≤ 1
    for (const auto& m : sub_result.mappings) {
        int deg = tf_degree_in(root, m.indeterminate);
        if (deg < 0 || deg > 1) return false;
    }

    /// 检查原始变量的次数是否 ≤ 1
    if (depends_on_var(root, var)) {
        int deg = tf_degree_in(root, var);
        if (deg < 0 || deg > 1) return false;
    }

    return true;
}

/**
 * @brief 检测表达式是否已为独立子表达式的乘积形式。
 *
 * 若表达式根节点为 MultiplyNode，则将各操作数视为独立因子直接返回，
 * 避免进入完整的 Berlekamp/Zassenhaus 分解流程。
 * 数值常数被单独累积，仅当非 1 时作为独立因子返回。
 *
 * @param[in] expr 待检测的符号表达式
 * @return 因子列表；若表达式非乘积形式则返回空向量（表示无快速路径）
 * @internal
 */
static std::vector<std::shared_ptr<SymbolicExpr>> tf_detect_multiplicative_structure(
    const std::shared_ptr<SymbolicExpr>& expr) {

    if (!expr || !lamina::detail::node(expr)) return {};

    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr));
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
                    /// 浮点数值：作为独立因子保留
                    factors.push_back(lamina::detail::make_expression_ptr(op));
                }
            }
            continue;
        }

        /// 非数值操作数作为独立因子
        factors.push_back(lamina::detail::make_expression_ptr(op));
    }

    /// 仅当存在至少两个非常数因子（或一个非常数因子加一个非 1 常数）时才视为有效乘积分解
    if (factors.size() < 2 && (factors.empty() || constant_acc == Rational(1))) {
        return {};
    }

    /// 插入累积常数因子（若非 1）
    if (constant_acc != Rational(1)) {
        auto const_expr = SymbolicExpr::number(constant_acc);
        factors.insert(factors.begin(), const_expr);
    }

    return factors;
}

/**
 * @brief 从乘积项中提取指数函数因子。
 *
 * 若节点本身为 exp(f(x)) 形式，直接返回该节点。
 * 若节点为 MultiplyNode，遍历其操作数寻找 exp(f(x)) 因子。
 * 仅提取第一个匹配的指数函数因子。
 *
 * @param[in] node 待检测的 AST 节点
 * @param[in] var  目标变量名
 * @return 找到的 exp 因子节点；未找到返回 nullptr
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
            depends_on_var(func->arguments()[0], var)) {
            return node;
        }
    }

    /// 乘积形式：遍历操作数寻找 exp 因子
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            auto func = std::dynamic_pointer_cast<const FunctionNode>(op);
            if (func &&
                func->type() == FunctionNode::FuncType::Exp &&
                func->arguments().size() == 1 &&
                depends_on_var(func->arguments()[0], var)) {
                return op;
            }
        }
    }

    return nullptr;
}

/**
 * @brief 从乘积项中移除指定的指数函数因子，返回剩余部分。
 *
 * 若节点本身即为该 exp 因子，返回数值 1。
 * 若节点为 MultiplyNode，移除匹配的 exp 操作数后重建乘积。
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
        return lamina::detail::make_node<NumberNode>(BigInt(1));
    }

    /// 乘积形式：移除匹配的操作数
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
            return lamina::detail::make_node<NumberNode>(BigInt(1));
        }
        if (remaining_ops.size() == 1) {
            return remaining_ops[0];
        }
        return lamina::detail::make_node<MultiplyNode>(std::move(remaining_ops));
    }

    return node;
}

/**
 * @brief 检测加法表达式中的公共指数因子并执行分离。
 *
 * 对于 AddNode 形式的表达式，检查所有加法项是否共享相同的 exp(f(x)) 因子。
 * 若是，则提取公因子：expr = exp(f(x)) * (t1' + t2' + ... + tn')，
 * 其中 ti' = ti / exp(f(x))。
 *
 * 典型用例：
 * - exp(x)*x + exp(x) → [exp(x), x+1]
 * - exp(x)*x² + 2*exp(x)*x + exp(x) → [exp(x), x²+2x+1]
 *
 * @param[in] expr 待检测的符号表达式
 * @param[in] var  目标变量名
 * @return 因子列表 [exp(f(x)), remaining_sum]；若无公共 exp 因子则返回空向量
 * @internal
 */
static std::vector<std::shared_ptr<SymbolicExpr>> tf_detect_exponential_separation(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    if (!expr || !lamina::detail::node(expr)) return {};

    auto add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(expr));
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

    /// 所有项共享相同的 exp(f(x))，执行分离
    /// 构造剩余和：对每个项移除 exp 因子
    std::vector<std::shared_ptr<const SymbolicNode>> remainder_terms;
    remainder_terms.reserve(add->operands().size());

    for (const auto& op : add->operands()) {
        auto remainder = tf_remove_exp_factor(op, common_exp);
        remainder_terms.push_back(remainder);
    }

    /// 构造结果
    auto exp_factor_expr = lamina::detail::make_expression_ptr(common_exp);

    std::shared_ptr<SymbolicExpr> sum_expr;
    if (remainder_terms.size() == 1) {
        sum_expr = lamina::detail::make_expression_ptr(remainder_terms[0]);
    } else {
        sum_expr = lamina::detail::make_expression_ptr(
            lamina::detail::make_node<AddNode>(std::move(remainder_terms)));
    }

    /// 化简剩余和
    auto simplified_sum = sum_expr->simplify();
    if (simplified_sum && lamina::detail::node(simplified_sum)) {
        sum_expr = simplified_sum;
    }

    return {exp_factor_expr, sum_expr};
}

// ============================================================
/// 公共 API 实现
// ============================================================

/**
 * @brief 判断表达式 AST 中是否包含依赖指定变量的超越函数。
 *
 * 递归遍历 AST，若发现任何 FunctionNode 类型为 Sin/Cos/Exp/Ln/Tan
 * 且其参数依赖 var，则返回 true。
 *
 * @param[in] node 当前 AST 节点
 * @param[in] var  目标变量名
 * @return 包含超越函数返回 true
 * @internal
 */
static bool tf_contains_transcendental(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var) {

    if (!node) return false;

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (func->arguments().size() == 1 &&
            tf_is_transcendental_type(func->type()) &&
            depends_on_var(func->arguments()[0], var)) {
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

// ============================================================
/// 特殊情形快速路径：毕达哥拉斯恒等式化简
// ============================================================

/**
 * @brief 判断节点是否为 sin²(f) 或 cos²(f) 形式，并提取函数类型和参数。
 *
 * 匹配模式：PowerNode(FunctionNode(Sin/Cos, [f]), NumberNode(2))
 *
 * @param[in]  node      待检测的 AST 节点
 * @param[out] func_type 输出函数类型（Sin 或 Cos）
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
 * @brief 从乘积项中提取系数和 sin²/cos² 核心部分。
 *
 * 对于形如 a*sin²(f) 的项，提取系数 a 和 sin²(f) 部分。
 * 若项本身即为 sin²(f)，系数为 1。
 *
 * @param[in]  node       待分析的加法操作数节点
 * @param[out] coeff      输出系数节点（nullptr 表示系数为 1）
 * @param[out] func_type  输出函数类型（Sin 或 Cos）
 * @param[out] argument   输出函数参数节点
 * @return 匹配成功返回 true
 * @internal
 */
static bool tf_extract_coeff_trig_squared(
    const std::shared_ptr<const SymbolicNode>& node,
    std::shared_ptr<const SymbolicNode>& coeff,
    FunctionNode::FuncType& func_type,
    std::shared_ptr<const SymbolicNode>& argument) {

    /// 直接为 sin²(f) 或 cos²(f)
    if (tf_is_trig_squared(node, func_type, argument)) {
        coeff = nullptr;  // 系数为 1
        return true;
    }

    /// 乘积形式：a * sin²(f) 或 sin²(f) * a
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!mul || mul->operands().size() < 2) return false;

    /// 在操作数中寻找 sin²(f) 或 cos²(f) 部分
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
                coeff = lamina::detail::make_node<MultiplyNode>(std::move(coeff_ops));
            }
            return true;
        }
    }

    return false;
}

/**
 * @brief 对表达式执行毕达哥拉斯恒等式化简：sin²(f) + cos²(f) → 1。
 *
 * 递归遍历 AST，在每个 AddNode 中扫描操作数对，寻找满足以下模式的配对：
 * - sin²(f) + cos²(f) → 替换为 1
 * - a*sin²(f) + a*cos²(f) → 替换为 a（公共系数）
 *
 * 要求 sin² 和 cos² 的参数 f 结构相等，且公共系数结构相等。
 *
 * @param[in] node 当前 AST 节点
 * @param[in] var  目标变量名（用于限定化简范围）
 * @return 化简后的节点；若无可化简的模式则返回原节点
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

        /// 在化简后的操作数中寻找 sin²(f) + cos²(f) 配对
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

                /// 找到配对：sin²(f) + cos²(f) → 1，或 a*sin²(f) + a*cos²(f) → a
                consumed[i] = true;
                consumed[j] = true;
                found_pair = true;

                if (!coeff_i) {
                    /// 系数为 1：替换为 NumberNode(1)
                    result_ops.push_back(lamina::detail::make_node<NumberNode>(BigInt(1)));
                } else {
                    /// 有公共系数：替换为系数本身
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
            return lamina::detail::make_node<NumberNode>(BigInt(0));
        }
        if (result_ops.size() == 1) {
            return result_ops[0];
        }
        return lamina::detail::make_node<AddNode>(std::move(result_ops));
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        new_ops.reserve(mul->operands().size());
        for (const auto& op : mul->operands()) {
            new_ops.push_back(tf_simplify_pythagorean_node(op, var));
        }
        return lamina::detail::make_node<MultiplyNode>(std::move(new_ops));
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto new_base = tf_simplify_pythagorean_node(pow->base(), var);
        auto new_exp = tf_simplify_pythagorean_node(pow->exponent(), var);
        return lamina::detail::make_node<PowerNode>(std::move(new_base), std::move(new_exp));
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        new_args.reserve(func->arguments().size());
        for (const auto& arg : func->arguments()) {
            new_args.push_back(tf_simplify_pythagorean_node(arg, var));
        }
        return lamina::detail::make_node<FunctionNode>(func->type(), std::move(new_args));
    }

    /// 叶节点（NumberNode、VariableNode）直接返回
    return node;
}

/**
 * @brief 对符号表达式执行毕达哥拉斯恒等式化简预处理。
 *
 * 在因式分解主流程之前调用，将表达式中的 sin²(f) + cos²(f) 子模式
 * 替换为 1，减少后续换元和分解的复杂度。
 *
 * @param[in] expr 待化简的符号表达式
 * @param[in] var  目标变量名
 * @return 化简后的表达式；若无可化简模式则返回原表达式
 * @internal
 */
static std::shared_ptr<SymbolicExpr> tf_simplify_pythagorean(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    if (!expr || !lamina::detail::node(expr)) return expr;

    auto simplified_root = tf_simplify_pythagorean_node(lamina::detail::node(expr), var);

    if (!simplified_root) return expr;

    /// 若化简后与原表达式结构相同，返回原表达式避免不必要的重建
    if (simplified_root->equals(*lamina::detail::node(expr))) {
        return expr;
    }

    return lamina::detail::make_expression_ptr(simplified_root);
}

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
