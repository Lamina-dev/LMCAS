/**
 * @file transcendental_rewrite.cpp
 * @brief 超越因子的逆换元、常数提取与乘法结构识别。
 */

#include "internal/expression_analysis.hpp"

#include <string>
#include <vector>
#include <cmath>

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

    if (!factor_expr || !LMCAS::detail::node(factor_expr)) {
        return factor_expr;
    }
    if (mappings.empty()) {
        return factor_expr;
    }

    auto result = factor_expr;

    for (const auto& m : mappings) {
        if (!m.trans_expr || m.indeterminate.empty()) {
            continue;
        }

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
    if (!num_node) {
        return false;
    }

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
            out = Rational::from_double(v);
            return true;
        }
        return false;
    }
    return false;
}

namespace {

class FactorSimplification {
public:
    std::vector<std::shared_ptr<SymbolicExpr>> run(
        const std::vector<std::shared_ptr<SymbolicExpr>>& factors) {
        for (const auto& factor : factors) {
            if (!factor || !LMCAS::detail::node(factor)) {
                continue;
            }
            auto simplified = factor->simplify();
            if (!simplified || !LMCAS::detail::node(simplified)) {
                simplified = factor;
            }
            if (!append(simplified)) {
                return std::move(result_);
            }
        }
        if (constant_ != Rational(1)) {
            result_.insert(result_.begin(), SymbolicExpr::number(constant_));
        }
        return std::move(result_);
    }

private:
    bool append(const std::shared_ptr<SymbolicExpr>& factor) {
        const auto& node = LMCAS::detail::node(factor);
        if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
            return append_number(number, factor);
        }
        if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
            append_product(*product, factor);
        } else {
            result_.push_back(factor);
        }
        return true;
    }

    bool append_number(const std::shared_ptr<const NumberNode>& number,
                       const std::shared_ptr<SymbolicExpr>& factor) {
        Rational value;
        if (!tf_extract_rational(number, value)) {
            result_.push_back(factor);
            return true;
        }
        if (value == Rational(0)) {
            result_.clear();
            result_.push_back(SymbolicExpr::number(0));
            return false;
        }
        constant_ = constant_ * value;
        return true;
    }

    void append_mixed_product(
        const std::vector<std::shared_ptr<const SymbolicNode>>& numeric,
        std::vector<std::shared_ptr<const SymbolicNode>>& symbolic) {
        for (const auto& operand : numeric) {
            auto number = std::dynamic_pointer_cast<const NumberNode>(operand);
            Rational value;
            if (tf_extract_rational(number, value)) {
                constant_ = constant_ * value;
            } else {
                symbolic.push_back(operand);
            }
        }
        if (symbolic.size() == 1) {
            result_.push_back(LMCAS::detail::make_expression_ptr(symbolic[0]));
        } else {
            result_.push_back(LMCAS::detail::make_expression_ptr(
                LMCAS::detail::make_node<MultiplyNode>(std::move(symbolic))));
        }
    }

    void accumulate_numeric_product(
        const std::vector<std::shared_ptr<const SymbolicNode>>& numeric) {
        Rational product(1);
        for (const auto& operand : numeric) {
            auto number = std::dynamic_pointer_cast<const NumberNode>(operand);
            Rational value;
            if (tf_extract_rational(number, value)) {
                product = product * value;
            }
        }
        constant_ = constant_ * product;
    }

    void append_product(const MultiplyNode& product,
                        const std::shared_ptr<SymbolicExpr>& factor) {
        std::vector<std::shared_ptr<const SymbolicNode>> numeric;
        std::vector<std::shared_ptr<const SymbolicNode>> symbolic;
        for (const auto& operand : product.operands()) {
            if (std::dynamic_pointer_cast<const NumberNode>(operand)) {
                numeric.push_back(operand);
            } else {
                symbolic.push_back(operand);
            }
        }
        if (!numeric.empty() && !symbolic.empty()) {
            append_mixed_product(numeric, symbolic);
        } else if (numeric.empty()) {
            result_.push_back(factor);
        } else {
            accumulate_numeric_product(numeric);
        }
    }

    std::vector<std::shared_ptr<SymbolicExpr>> result_;
    Rational constant_{1};
};

}

/**
 * @brief 化简逆换元后的因子，合并数值因子及乘积的数值前导系数。
 * 保留非常数部分，将非 1 的常数积置于结果首位。
 * @param[in,out] factors 就地化简的因子列表
 * @return 化简并提取常数后的因子列表
 */
std::vector<std::shared_ptr<SymbolicExpr>> tf_simplify_factors(
    std::vector<std::shared_ptr<SymbolicExpr>>& factors) {
    FactorSimplification simplification;
    return simplification.run(factors);
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

    if (!sub_result.poly_expr || !LMCAS::detail::node(sub_result.poly_expr)) {
        return false;
    }
    if (sub_result.mappings.empty()) {
        return false;
    }

    const auto& root = LMCAS::detail::node(sub_result.poly_expr);

    /// 检查每个不定元的次数是否 <= 1
    for (const auto& m : sub_result.mappings) {
        int deg = tf_degree_in(root, m.indeterminate);
        if (deg < 0 || deg > 1) {
            return false;
        }
    }

    /// 检查原始变量的次数是否 <= 1
    if (expression_depends_on_variable(root, var)) {
        int deg = tf_degree_in(root, var);
        if (deg < 0 || deg > 1) {
            return false;
        }
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

    if (!expr || !LMCAS::detail::node(expr)) {
        return {};
    }

    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr));
    if (!mul || mul->operands().size() < 2) {
        return {};
    }

    std::vector<std::shared_ptr<SymbolicExpr>> factors;
    Rational constant_acc(1);

    for (const auto& op : mul->operands()) {
        if (!op) {
            continue;
        }

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

static bool tf_is_variable_exponential(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    auto function = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!function) {
        return false;
    }
    if (function->type() != FunctionNode::FuncType::Exp ||
        function->arguments().size() != 1) {
        return false;
    }
    return expression_depends_on_variable(function->arguments()[0], var);
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

    if (!node) {
        return nullptr;
    }

    /// 直接为 exp(f(x)) 形式
    if (tf_is_variable_exponential(node, var)) {
        return node;
    }

    /// 乘积形式:遍历操作数寻找 exp 因子
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            if (tf_is_variable_exponential(op, var)) {
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

    if (!node || !exp_factor) {
        return node;
    }

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

        if (!removed) {
            return node;
        }

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

    if (!expr || !LMCAS::detail::node(expr)) {
        return {};
    }

    auto add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(expr));
    if (!add || add->operands().size() < 2) {
        return {};
    }

    /// 从第一个加法项中提取 exp 因子作为候选公因子
    std::shared_ptr<const SymbolicNode> common_exp = tf_extract_exp_factor(add->operands()[0], var);
    if (!common_exp) {
        return {};
    }

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
static bool any_transcendental(
    const std::vector<std::shared_ptr<const SymbolicNode>>& nodes,
    const std::string& var) {
    for (const auto& node : nodes) {
        if (tf_contains_transcendental(node, var)) {
            return true;
        }
    }
    return false;
}

bool tf_contains_transcendental(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var) {

    if (!node) {
        return false;
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (func->arguments().size() == 1 &&
            tf_is_transcendental_type(func->type()) &&
            expression_depends_on_variable(func->arguments()[0], var)) {
            return true;
        }
        return any_transcendental(func->arguments(), var);
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return any_transcendental(add->operands(), var);
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return any_transcendental(mul->operands(), var);
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return tf_contains_transcendental(pow->base(), var) ||
               tf_contains_transcendental(pow->exponent(), var);
    }

    return false;
}

}
