/**
 * @file limit_visitor_taylor.cpp
 * @brief LimitVisitor 的 Taylor 展开回退实现.
 *
 * L'Hôpital 规则达到最大迭代深度后,
 * Taylor 级数展开分子和分母,并以首项系数比继续求极限.
 *
 * 算法来源:标准 CAS Taylor 级数极限技术
 */

#include "symbolic.hpp"
#include "internal/visitors/limit_visitor.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "internal/expression_analysis.hpp"

namespace LMCAS {

namespace {
std::string substitute_infinity_series(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den, const std::string& var,
    std::shared_ptr<SymbolicExpr>& num_expr, std::shared_ptr<SymbolicExpr>& den_expr) {
    /// x -> infinity: 代换 x = 1/t,在 t -> 0 处展开
    std::string t_var = "__lim_t__";
    auto names = detail::all_variable_names(num);
    auto denominator_names = detail::all_variable_names(den);
    names.insert(denominator_names.begin(), denominator_names.end());
    while (names.count(t_var)) { t_var += "_"; }
    auto t_expr = SymbolicExpr::variable(t_var);
    auto one_over_t = SymbolicExpr::power(t_expr, SymbolicExpr::number(-1));

    num_expr = num_expr->substitute(var, one_over_t);
    den_expr = den_expr->substitute(var, one_over_t);

    if (num_expr) { num_expr = num_expr->simplify(); }
    if (den_expr) { den_expr = den_expr->simplify(); }

    return t_var;
}
} // namespace

std::shared_ptr<const SymbolicNode> LimitVisitor::divergent_leading_ratio(
    const std::shared_ptr<SymbolicExpr>& ratio, int difference) {
    auto sign = ratio ? get_node_sign(LMCAS::detail::node(ratio)) : std::nullopt;
    if (!sign || !*sign) { return nullptr; }
    if (difference % 2 && direction.empty() && !is_inf(point)) { return nullptr; }
    if (difference % 2 && direction == "-" && !is_inf(point)) { *sign = -*sign; }
    return signed_infinity(*sign);
}



/**
 * @brief 对 L'Hôpital 法则产生的导数比进行代数化简后求极限.
 *
 * 当 dN/dD 含有公因子(如 x^-2)时,直接分别求极限会导致
 * 无限 0/0 循环.此方法先通过 simplify() 约去公因子,
 * 再对化简后的表达式求极限.
 *
 * @param[in] ratio_node 导数比 dN * dD^(-1) 的 AST 节点
 * @return 极限结果,化简无效时返回 nullptr
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::simplify_and_eval_ratio(
    const std::shared_ptr<const SymbolicNode>& ratio_node) {

    if (!ratio_node) {
        return nullptr;
    }

    auto ratio_expr = LMCAS::detail::make_expression_ptr(ratio_node->clone());
    auto simplified = ratio_expr->simplify();
    if (!simplified || !LMCAS::detail::node(simplified)) {
        return nullptr;
    }

    if (detail::node(ratio_expr)->compare(*detail::node(simplified)) == 0) return nullptr;

    /// Reject simplified forms that are AddNodes (sums). When simplify() expands
    /// a fraction like (1-cos(x))/(sin(x)+x*cos(x)) into a sum of terms with
    /// negative powers, evaluating the limit of that sum can trigger infinity-infinity detection,
    /// which calls resolve_inf_minus_inf -> apply_lhopital -> simplify_and_eval_ratio
    /// again, creating an infinite loop.
    if (std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(simplified))) {
        return nullptr;
    }

    return derivative_limit(LMCAS::detail::node(simplified), true);
}

/**
 * @brief Taylor 展开回退策略实现.
 *
 * 当极限点为有限值时,直接在该点展开 Taylor 级数.
 * 当极限点为无穷时,先做 x = 1/t 代换,再在 t = 0 处展开.
 * 从 order=4 开始逐步增加到 max_order,直到找到非零首项.
 *
 * @param[in] num 分子 AST 节点
 * @param[in] den 分母 AST 节点
 * @param[in] max_order 最大展开阶数
 * @return 极限结果节点;nullptr 表示当前规则保持结果未知
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::taylor_fallback(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den,
    int max_order) {

    if (!num || !den) {
        return nullptr;
    }
    if (!is_inf(point) && (!analytic_at_point(num) || !analytic_at_point(den))) { return nullptr; }

    auto num_expr = LMCAS::detail::make_expression_ptr(num->clone());
    auto den_expr = LMCAS::detail::make_expression_ptr(den->clone());
    auto point_expr = LMCAS::detail::make_expression_ptr(point->clone());

    bool at_infinity = is_inf(point);

    std::string expand_var = var;
    std::shared_ptr<SymbolicExpr> expand_point;

    if (at_infinity) {
        expand_var = substitute_infinity_series(num, den, var, num_expr, den_expr);
        expand_point = SymbolicExpr::number(0);
    } else {
        expand_point = point_expr;
    }

    if (!num_expr || !den_expr || !expand_point) {
        return nullptr;
    }
    if (at_infinity) {
        LimitVisitor analytic(expand_var, detail::node(expand_point), "+", assumption_ctx_, context_);
        if (!analytic.analytic_at_point(detail::node(num_expr)) ||
            !analytic.analytic_at_point(detail::node(den_expr))) { return nullptr; }
    }

    /// 从 order=4 开始,逐步增加到 max_order
    for (int order = 4; order <= max_order; order += 2) {
        auto value = taylor_series_limit(num_expr, den_expr, expand_var, expand_point, order);
        if (value) {
            return value;
        }
    }

    return nullptr;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::taylor_series_limit(
    const std::shared_ptr<SymbolicExpr>& numerator,
    const std::shared_ptr<SymbolicExpr>& denominator,
    const std::string& expand_var,
    const std::shared_ptr<SymbolicExpr>& expand_point, int order) {
    auto numerator_series = numerator->series(expand_var, expand_point, order);
    auto denominator_series = denominator->series(expand_var, expand_point, order);
    if (!numerator_series || !denominator_series) {
        return nullptr;
    }
    numerator_series = numerator_series->simplify();
    denominator_series = denominator_series->simplify();
    auto numerator_leading = find_leading_term(numerator_series, expand_var, expand_point, order);
    auto denominator_leading = find_leading_term(denominator_series, expand_var, expand_point, order);
    return leading_terms_limit(numerator_leading, denominator_leading);
}

std::shared_ptr<const SymbolicNode> LimitVisitor::leading_terms_limit(
    const std::pair<std::shared_ptr<const SymbolicNode>, int>& numerator,
    const std::pair<std::shared_ptr<const SymbolicNode>, int>& denominator) {
    if (!numerator.first || !denominator.first) {
        return nullptr;
    }
    if (denominator.first->is_zero()) {
        return nullptr;
    }
    if (numerator.first->is_zero() || numerator.second > denominator.second) {
        return LMCAS::detail::make_node<NumberNode>(BigInt(0));
    }
    auto ratio = SymbolicExpr::multiply(
        LMCAS::detail::make_expression_ptr(numerator.first),
        SymbolicExpr::power(LMCAS::detail::make_expression_ptr(denominator.first),
                            SymbolicExpr::number(-1)));
    ratio = ratio->simplify();
    if (numerator.second < denominator.second) {
        return divergent_leading_ratio(ratio, denominator.second - numerator.second);
    }
    return ratio ? LMCAS::detail::node(ratio) : nullptr;
}

/**
 * @internal
 * @brief 从 Taylor 级数中提取关于 (x - a) 的首个非零项.
 *
 * 通过逐阶求导并在展开点求值来确定首个非零系数及其阶数.
 *
 * @param[in] series_expr 级数表达式
 * @param[in] expand_var 展开变量名
 * @param[in] expand_point 展开中心点
 * @param[in] max_order 最大检查阶数
 * @return (首项系数节点, 首项阶数) 对
 */
std::pair<std::shared_ptr<const SymbolicNode>, int> LimitVisitor::find_leading_term(
    const std::shared_ptr<SymbolicExpr>& series_expr,
    const std::string& expand_var,
    const std::shared_ptr<SymbolicExpr>& expand_point,
    int max_order) {

    if (!series_expr) {
        return {nullptr, 0};
    }

    auto current = series_expr;

    for (int n = 0; n <= max_order; ++n) {
        /// 在展开点求值得到第 n 阶系数(乘以 n!)
        auto val = current->substitute(expand_var, expand_point);
        if (!val) {
            return {nullptr, 0};
        }
        val = val->simplify();

        auto sign = val ? get_node_sign(LMCAS::detail::node(val)) : std::nullopt;
        if (!sign) return {nullptr, 0};
        if (*sign) return {LMCAS::detail::node(val), n};

        /// 对当前表达式求导以获取下一阶系数
        current = current->differentiate(expand_var);
        if (!current) {
            return {nullptr, 0};
        }
        current = current->simplify();
    }

    return {nullptr, max_order + 1};
}

/**
 * @internal
 * @brief 判断数值节点的符号.
 * @param[in] node AST 节点
 * @return 正数返回 1,负数返回 -1,零或未知符号返回 0
 */
int LimitVisitor::get_sign(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) {
        return 0;
    }
    if (auto num = std::dynamic_pointer_cast<const NumberNode>(node)) {
        if (std::holds_alternative<double>(num->value())) {
            double v = std::get<double>(num->value());
            if (v > 0) {
                return 1;
            }
            if (v < 0) {
                return -1;
            }
            return 0;
        }
        if (std::holds_alternative<BigInt>(num->value())) {
            BigInt v = std::get<BigInt>(num->value());
            if (v > BigInt(0)) {
                return 1;
            }
            if (v < BigInt(0)) {
                return -1;
            }
            return 0;
        }
        if (std::holds_alternative<Rational>(num->value())) {
            Rational v = std::get<Rational>(num->value());
            if (v > Rational(0)) {
                return 1;
            }
            if (v < Rational(0)) {
                return -1;
            }
            return 0;
        }
    }
    /// 对于乘法节点,符号为各因子符号之积
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        int sign = 1;
        for (auto& op : mul->operands()) {
            int s = get_sign(op);
            if (s == 0) {
                return 0;
            }
            sign *= s;
        }
        return sign;
    }
    return 0;
}


} // namespace LMCAS
