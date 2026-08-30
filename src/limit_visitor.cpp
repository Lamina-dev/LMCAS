/**
 * @file limit_visitor.cpp
 * @brief LimitVisitor 的 Taylor 展开回退实现。
 *
 * L'Hôpital 规则达到最大迭代深度后，
 * Taylor 级数展开分子和分母，并以首项系数比继续求极限。
 *
 * 算法来源：标准 CAS Taylor 级数极限技术
 */

#include "../include/symbolic.hpp"
#include "../include/visitors/limit_visitor.hpp"
#include "../include/visitors/normalization_visitor.hpp"

/**
 * @brief 对 L'Hôpital 法则产生的导数比进行代数化简后求极限。
 *
 * 当 dN/dD 含有公因子（如 x^-2）时，直接分别求极限会导致
 * 无限 0/0 循环。此方法先通过 simplify() 约去公因子，
 * 再对化简后的表达式求极限。
 *
 * @param[in] ratio_node 导数比 dN * dD^(-1) 的 AST 节点
 * @return 极限结果，化简无效时返回 nullptr
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::simplify_and_eval_ratio(
    const std::shared_ptr<const SymbolicNode>& ratio_node) {

    if (!ratio_node) return nullptr;

    auto ratio_expr = lamina::detail::make_expression_ptr(ratio_node->clone());
    auto simplified = ratio_expr->simplify();
    if (!simplified || !lamina::detail::node(simplified)) return nullptr;

    /// Check if simplification actually changed the expression (avoid infinite loops)
    auto orig_str = ratio_expr->to_string();
    auto simp_str = simplified->to_string();
    if (orig_str == simp_str) return nullptr;

    /// Reject simplified forms that are AddNodes (sums). When simplify() expands
    /// a fraction like (1-cos(x))/(sin(x)+x*cos(x)) into a sum of terms with
    /// negative powers, evaluating the limit of that sum can trigger ∞−∞ detection,
    /// which calls resolve_inf_minus_inf → apply_lhopital → simplify_and_eval_ratio
    /// again, creating an infinite loop.
    if (std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(simplified))) return nullptr;

    LimitVisitor sub_simp(var, point, direction, assumption_ctx_);
    sub_simp.lhopital_depth_ = this->lhopital_depth_ + 1;
    lamina::detail::node(simplified)->accept(sub_simp);
    auto simp_result = sub_simp.get_result();
    if (!simp_result) return nullptr;

    NormalizationVisitor n2;
    simp_result->accept(n2);
    simp_result = n2.get_result();
    if (!simp_result) return nullptr;

    /// Accept the result if it's a determinate finite value
    if (!is_inf(simp_result)) {
        return simp_result;
    }
    /// Also accept infinity results
    return simp_result;
}

/**
 * @brief Taylor 展开回退策略实现。
 *
 * 当极限点为有限值时，直接在该点展开 Taylor 级数。
 * 当极限点为无穷时，先做 x = 1/t 代换，再在 t = 0 处展开。
 * 从 order=4 开始逐步增加到 max_order，直到找到非零首项。
 *
 * @param[in] num 分子 AST 节点
 * @param[in] den 分母 AST 节点
 * @param[in] max_order 最大展开阶数
 * @return 极限结果节点；nullptr 表示当前规则保持结果未知
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::taylor_fallback(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den,
    int max_order) {

    if (!num || !den) return nullptr;

    auto num_expr = lamina::detail::make_expression_ptr(num->clone());
    auto den_expr = lamina::detail::make_expression_ptr(den->clone());
    auto point_expr = lamina::detail::make_expression_ptr(point->clone());

    bool at_infinity = is_inf(point);

    std::string expand_var = var;
    std::shared_ptr<SymbolicExpr> expand_point;

    if (at_infinity) {
        /// x → ∞: 代换 x = 1/t，在 t → 0 处展开
        std::string t_var = "__lim_t__";
        auto t_expr = SymbolicExpr::variable(t_var);
        auto one_over_t = SymbolicExpr::power(t_expr, SymbolicExpr::number(-1));

        num_expr = num_expr->substitute(var, one_over_t);
        den_expr = den_expr->substitute(var, one_over_t);

        if (num_expr) num_expr = num_expr->simplify();
        if (den_expr) den_expr = den_expr->simplify();

        expand_var = t_var;
        expand_point = SymbolicExpr::number(0);
    } else {
        expand_point = point_expr;
    }

    if (!num_expr || !den_expr || !expand_point) return nullptr;

    /// 从 order=4 开始，逐步增加到 max_order
    for (int order = 4; order <= max_order; order += 2) {
        auto num_series = num_expr->series(expand_var, expand_point, order);
        auto den_series = den_expr->series(expand_var, expand_point, order);

        if (!num_series || !den_series) continue;

        num_series = num_series->simplify();
        den_series = den_series->simplify();

        /// 提取首项：在展开点处求值得到常数项，
        /// 若为零则对 (x - a) 的各阶系数逐一检查
        auto num_leading = find_leading_term(num_series, expand_var, expand_point, order);
        auto den_leading = find_leading_term(den_series, expand_var, expand_point, order);

        if (!num_leading.first || !den_leading.first) continue;

        /// 两个首项都非零 → 可以确定极限
        if (!num_leading.first->is_zero() && !den_leading.first->is_zero()) {
            int power_diff = num_leading.second - den_leading.second;

            if (power_diff > 0) {
                /// 分子阶数更高 → 极限为 0
                return lamina::detail::make_node<NumberNode>(BigInt(0));
            } else if (power_diff < 0) {
                /// 分母阶数更高 → 极限为 ±∞
                /// 确定符号
                auto ratio = SymbolicExpr::multiply(
                    lamina::detail::make_expression_ptr(num_leading.first),
                    SymbolicExpr::power(
                        lamina::detail::make_expression_ptr(den_leading.first),
                        SymbolicExpr::number(-1)));
                ratio = ratio->simplify();
                if (!ratio) {
                    std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
                    return lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
                }
                int sign = get_sign(lamina::detail::node(ratio));
                std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
                auto inf_node = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
                if (sign < 0) {
                    std::vector<std::shared_ptr<const SymbolicNode>> m = {
                        lamina::detail::make_node<NumberNode>(BigInt(-1)), inf_node};
                    return lamina::detail::make_node<MultiplyNode>(m);
                }
                return inf_node;
            } else {
                /// 同阶 → 极限为系数之比
                auto ratio = SymbolicExpr::multiply(
                    lamina::detail::make_expression_ptr(num_leading.first),
                    SymbolicExpr::power(
                        lamina::detail::make_expression_ptr(den_leading.first),
                        SymbolicExpr::number(-1)));
                ratio = ratio->simplify();
                if (ratio) return lamina::detail::node(ratio);
            }
        }

        /// 如果分子首项为零但分母不为零，极限为 0
        if (num_leading.first->is_zero() && !den_leading.first->is_zero()) {
            return lamina::detail::make_node<NumberNode>(BigInt(0));
        }

        /// 两者都为零 → 需要更高阶展开
    }

    return nullptr;
}

/**
 * @internal
 * @brief 从 Taylor 级数中提取关于 (x - a) 的首个非零项。
 *
 * 通过逐阶求导并在展开点求值来确定首个非零系数及其阶数。
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

    if (!series_expr) return {nullptr, 0};

    auto current = series_expr;

    for (int n = 0; n <= max_order; ++n) {
        /// 在展开点求值得到第 n 阶系数（乘以 n!）
        auto val = current->substitute(expand_var, expand_point);
        if (!val) return {nullptr, 0};
        val = val->simplify();

        if (val && lamina::detail::node(val) && !lamina::detail::node(val)->is_zero()) {
            /// 系数为 val / n!（但对于比较比值，n! 会约掉，所以直接返回 val）
            return {lamina::detail::node(val), n};
        }

        /// 对当前表达式求导以获取下一阶系数
        current = current->differentiate(expand_var);
        if (!current) return {nullptr, 0};
        current = current->simplify();
    }

    return {lamina::detail::make_node<NumberNode>(BigInt(0)), max_order + 1};
}

/**
 * @internal
 * @brief 判断数值节点的符号。
 * @param[in] node AST 节点
 * @return 正数返回 1，负数返回 -1，零或未知符号返回 0
 */
int LimitVisitor::get_sign(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) return 0;
    if (auto num = std::dynamic_pointer_cast<const NumberNode>(node)) {
        if (std::holds_alternative<double>(num->value())) {
            double v = std::get<double>(num->value());
            if (v > 0) return 1;
            if (v < 0) return -1;
            return 0;
        }
        if (std::holds_alternative<BigInt>(num->value())) {
            BigInt v = std::get<BigInt>(num->value());
            if (v > BigInt(0)) return 1;
            if (v < BigInt(0)) return -1;
            return 0;
        }
        if (std::holds_alternative<Rational>(num->value())) {
            Rational v = std::get<Rational>(num->value());
            if (v > Rational(0)) return 1;
            if (v < Rational(0)) return -1;
            return 0;
        }
    }
    /// 对于乘法节点，符号为各因子符号之积
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        int sign = 1;
        for (auto& op : mul->operands()) {
            int s = get_sign(op);
            if (s == 0) return 0;
            sign *= s;
        }
        return sign;
    }
    return 0;
}


/**
 * @brief 处理分段函数节点的极限。
 *
 * 根据趋近方向选择满足条件的分支并计算极限；
 * 双侧结果相等时返回该值，差异时以 nullptr 表示极限未定义。
 */
void LimitVisitor::visit(const PiecewiseNode& node) {
    if (direction == "+" || direction == "-") {
        auto branch_expr = select_branch_by_direction(node, direction);
        if (branch_expr) {
            LimitVisitor sv(var, point, direction, assumption_ctx_);
            branch_expr->accept(sv);
            result = sv.get_result();
        } else if (node.default_expr()) {
            LimitVisitor sv(var, point, direction, assumption_ctx_);
            node.default_expr()->accept(sv);
            result = sv.get_result();
        } else {
            result = nullptr;
        }
    } else {
        /// 双侧极限：分别计算左右极限
        auto lb = select_branch_by_direction(node, "-");
        auto rb = select_branch_by_direction(node, "+");
        std::shared_ptr<const SymbolicNode> lr = nullptr, rr = nullptr;
        if (lb) {
            LimitVisitor lv(var, point, "-", assumption_ctx_);
            lb->accept(lv);
            lr = lv.get_result();
        } else if (node.default_expr()) {
            LimitVisitor lv(var, point, "-", assumption_ctx_);
            node.default_expr()->accept(lv);
            lr = lv.get_result();
        }
        if (rb) {
            LimitVisitor rv(var, point, "+", assumption_ctx_);
            rb->accept(rv);
            rr = rv.get_result();
        } else if (node.default_expr()) {
            LimitVisitor rv(var, point, "+", assumption_ctx_);
            node.default_expr()->accept(rv);
            rr = rv.get_result();
        }
        if (!lr || !rr) { result = nullptr; return; }
        result = (lr->compare(*rr) == 0) ? lr : nullptr;
    }
}

/**
 * @brief 根据趋近方向选择分段函数中满足条件的分支表达式。
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::select_branch_by_direction(
    const PiecewiseNode& node, const std::string& dir) {
    for (const auto& branch : node.branches()) {
        if (condition_satisfied_by_direction(branch.condition, dir))
            return branch.expression;
    }
    return nullptr;
}

/**
 * @brief 判断条件在给定趋近方向下是否满足。
 *
 * 右极限（"+"）意味着 var > point，左极限（"-"）意味着 var < point。
 */
bool LimitVisitor::condition_satisfied_by_direction(
    const std::shared_ptr<const SymbolicNode>& condition, const std::string& dir) {
    auto rel = std::dynamic_pointer_cast<const RelationalNode>(condition);
    if (!rel) {
        auto logical = std::dynamic_pointer_cast<const LogicalNode>(condition);
        if (logical) {
            if (logical->op() == LogicalNode::Op::And)
                return condition_satisfied_by_direction(logical->left(), dir) &&
                       condition_satisfied_by_direction(logical->right(), dir);
            if (logical->op() == LogicalNode::Op::Or)
                return condition_satisfied_by_direction(logical->left(), dir) ||
                       condition_satisfied_by_direction(logical->right(), dir);
            if (logical->op() == LogicalNode::Op::Not && logical->left())
                return !condition_satisfied_by_direction(logical->left(), dir);
        }
        return false;
    }
    auto sign = evaluate_relational_sign(rel, dir);
    if (!sign) return false;
    int s = *sign;
    switch (rel->op()) {
        case RelationalNode::Op::GT:  return s > 0;
        case RelationalNode::Op::GEQ: return s >= 0;
        case RelationalNode::Op::LT:  return s < 0;
        case RelationalNode::Op::LEQ: return s <= 0;
        case RelationalNode::Op::EQ:  return s == 0;
        case RelationalNode::Op::NEQ: return s != 0;
        default: return false;
    }
}

/**
 * @brief 评估关系表达式 left - right 在趋近方向下的符号。
 */
std::optional<int> LimitVisitor::evaluate_relational_sign(
    const std::shared_ptr<const RelationalNode>& rel, const std::string& dir) {
    /// 模式：var op number
    auto left_var = std::dynamic_pointer_cast<const VariableNode>(rel->left());
    auto right_num = std::dynamic_pointer_cast<const NumberNode>(rel->right());
    if (left_var && left_var->name() == var && right_num) {
        double rv = get_numeric_value(right_num);
        double pv = get_point_value();
        if (std::isnan(pv) || std::isnan(rv)) return std::nullopt;
        double av = (dir == "+") ? pv + 1e-10 : pv - 1e-10;
        double diff = av - rv;
        if (diff > 1e-15) return 1;
        if (diff < -1e-15) return -1;
        return 0;
    }
    /// 模式：number op var
    auto left_num = std::dynamic_pointer_cast<const NumberNode>(rel->left());
    auto right_var = std::dynamic_pointer_cast<const VariableNode>(rel->right());
    if (left_num && right_var && right_var->name() == var) {
        double lv = get_numeric_value(left_num);
        double pv = get_point_value();
        if (std::isnan(pv) || std::isnan(lv)) return std::nullopt;
        double av = (dir == "+") ? pv + 1e-10 : pv - 1e-10;
        double diff = lv - av;
        if (diff > 1e-15) return 1;
        if (diff < -1e-15) return -1;
        return 0;
    }
    /// 通用：计算 (left - right) 的极限符号
    auto neg_one = lamina::detail::make_node<NumberNode>(BigInt(-1));
    std::vector<std::shared_ptr<const SymbolicNode>> neg_ops = {neg_one, rel->right()->clone()};
    std::vector<std::shared_ptr<const SymbolicNode>> add_ops = {rel->left()->clone(), lamina::detail::make_node<MultiplyNode>(neg_ops)};
    auto diff_expr = lamina::detail::make_node<AddNode>(add_ops);
    LimitVisitor sv(var, point, dir, assumption_ctx_);
    diff_expr->accept(sv);
    auto val = sv.get_result();
    if (!val) return std::nullopt;
    NormalizationVisitor norm;
    val->accept(norm);
    val = norm.get_result();
    return get_node_sign(val);
}

/**
 * @brief 方向感知的 sgn 函数极限计算。
 *
 * 当 sgn 的参数在趋近点为零时，根据趋近方向确定符号。
 */
std::optional<std::shared_ptr<const SymbolicNode>> LimitVisitor::evaluate_sgn_limit(
    const std::shared_ptr<const SymbolicNode>& arg) {
    LimitVisitor sv(var, point, direction, assumption_ctx_);
    arg->accept(sv);
    auto al = sv.get_result();
    if (!al) return std::nullopt;
    NormalizationVisitor norm;
    al->accept(norm);
    al = norm.get_result();
    if (!al->is_zero()) {
        auto s = get_node_sign(al);
        if (s) return lamina::detail::make_node<NumberNode>(BigInt(*s));
        return std::nullopt;
    }
    /// 参数极限为零，根据方向确定符号
    int sign = determine_sign_near_point(arg, direction);
    if (sign != 0) return lamina::detail::make_node<NumberNode>(BigInt(sign));
    return lamina::detail::make_node<NumberNode>(BigInt(0));
}

/**
 * @brief 方向感知的绝对值函数极限计算。
 */
std::optional<std::shared_ptr<const SymbolicNode>> LimitVisitor::evaluate_abs_limit(
    const std::shared_ptr<const SymbolicNode>& arg) {
    LimitVisitor sv(var, point, direction, assumption_ctx_);
    arg->accept(sv);
    auto al = sv.get_result();
    if (!al) return std::nullopt;
    NormalizationVisitor norm;
    al->accept(norm);
    al = norm.get_result();
    if (is_inf(al)) {
        std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
        return lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
    }
    auto s = get_node_sign(al);
    if (s) {
        if (*s >= 0) return al;
        std::vector<std::shared_ptr<const SymbolicNode>> neg_ops = {lamina::detail::make_node<NumberNode>(BigInt(-1)), al};
        auto neg_result = lamina::detail::make_node<MultiplyNode>(neg_ops);
        neg_result->accept(norm);
        return norm.get_result();
    }
    return std::nullopt;
}



/**
 * @brief 判断当前极限点是否为正无穷。
 * @return 极限点为 +∞ 时返回 true
 */
bool LimitVisitor::is_limit_at_infinity() const {
    if (auto f = std::dynamic_pointer_cast<const FunctionNode>(point)) {
        return f->type() == FunctionNode::FuncType::Infinity;
    }
    return false;
}

/**
 * @brief 判断当前极限点是否为负无穷。
 *
 * 负无穷表示为 -1 * Infinity 的 MultiplyNode。
 * @return 极限点为 -∞ 时返回 true
 */
bool LimitVisitor::is_limit_at_neg_infinity() const {
    return is_neg_inf(point);
}

/**
 * @brief 获取多项式表达式关于 var 的次数。
 *
 * 仅处理简单的多项式结构：
 * - VariableNode → 1
 * - NumberNode → 0
 * - PowerNode(var, n) → n
 * - MultiplyNode → 各因子次数之和
 * - AddNode → 各项次数的最大值
 *
 * @param[in] node AST 节点
 * @return 多项式次数，非多项式返回 -1
 */
int LimitVisitor::get_polynomial_degree(const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) return -1;

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return 0;
    }

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        return (v->name() == var) ? 1 : 0;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base_var = std::dynamic_pointer_cast<const VariableNode>(pow->base());
        if (base_var && base_var->name() == var) {
            if (auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
                double e = 0;
                if (std::holds_alternative<double>(exp_num->value())) e = std::get<double>(exp_num->value());
                else if (std::holds_alternative<BigInt>(exp_num->value())) e = std::get<BigInt>(exp_num->value()).to_double();
                else if (std::holds_alternative<Rational>(exp_num->value())) e = std::get<Rational>(exp_num->value()).to_double();
                if (e == static_cast<int>(e) && e >= 0) return static_cast<int>(e);
            }
        }
        /// c^x or similar — not a polynomial
        int base_deg = get_polynomial_degree(pow->base());
        if (base_deg < 0) return -1;
        if (auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
            double e = 0;
            if (std::holds_alternative<double>(exp_num->value())) e = std::get<double>(exp_num->value());
            else if (std::holds_alternative<BigInt>(exp_num->value())) e = std::get<BigInt>(exp_num->value()).to_double();
            else if (std::holds_alternative<Rational>(exp_num->value())) e = std::get<Rational>(exp_num->value()).to_double();
            if (e == static_cast<int>(e) && e >= 0) return base_deg * static_cast<int>(e);
        }
        return -1;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        int total = 0;
        for (auto& op : mul->operands()) {
            int d = get_polynomial_degree(op);
            if (d < 0) return -1;
            total += d;
        }
        return total;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        int max_deg = 0;
        for (auto& op : add->operands()) {
            int d = get_polynomial_degree(op);
            if (d < 0) return -1;
            if (d > max_deg) max_deg = d;
        }
        return max_deg;
    }

    /// FunctionNode (sin, cos, exp, ln, etc.) — not a polynomial
    return -1;
}

/**
 * @brief 获取多项式的首项系数。
 *
 * 对于多项式 P(x) = a_n * x^n + ... + a_0，返回 a_n。
 *
 * @param[in] node AST 节点
 * @return 首项系数节点；nullptr 表示当前结构保持未知
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::get_leading_coefficient(const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) return nullptr;

    int deg = get_polynomial_degree(node);
    if (deg < 0) return nullptr;

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return node->clone();
    }

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (v->name() == var) return lamina::detail::make_node<NumberNode>(BigInt(1));
        return node->clone();
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto base_var = std::dynamic_pointer_cast<const VariableNode>(pow->base());
        if (base_var && base_var->name() == var) {
            return lamina::detail::make_node<NumberNode>(BigInt(1));
        }
        /// constant^n
        auto base_lc = get_leading_coefficient(pow->base());
        if (!base_lc) return nullptr;
        if (auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
            return lamina::detail::make_node<PowerNode>(base_lc, pow->exponent()->clone());
        }
        return nullptr;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> lc_parts;
        for (auto& op : mul->operands()) {
            auto lc = get_leading_coefficient(op);
            if (!lc) return nullptr;
            lc_parts.push_back(lc);
        }
        if (lc_parts.size() == 1) return lc_parts[0];
        auto prod = lamina::detail::make_node<MultiplyNode>(lc_parts);
        NormalizationVisitor norm;
        prod->accept(norm);
        return norm.get_result();
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        /// Find terms with the highest degree and sum their leading coefficients
        std::vector<std::shared_ptr<const SymbolicNode>> leading_terms;
        for (auto& op : add->operands()) {
            int d = get_polynomial_degree(op);
            if (d == deg) {
                auto lc = get_leading_coefficient(op);
                if (lc) leading_terms.push_back(lc);
            }
        }
        if (leading_terms.empty()) return nullptr;
        if (leading_terms.size() == 1) return leading_terms[0];
        auto sum = lamina::detail::make_node<AddNode>(leading_terms);
        NormalizationVisitor norm;
        sum->accept(norm);
        return norm.get_result();
    }

    return nullptr;
}

/**
 * @brief 通过多项式次数比较计算有理函数在无穷处的极限。
 *
 * 对于 P(x)/Q(x)：
 * - deg(P) < deg(Q) → 0
 * - deg(P) = deg(Q) → 首项系数之比
 * - deg(P) > deg(Q) → ±∞（符号由首项系数决定）
 *
 * @param[in] num 分子 AST 节点
 * @param[in] den 分母 AST 节点
 * @return 极限结果，非有理函数时返回 nullptr
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::limit_rational_at_infinity(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den) {

    int deg_num = get_polynomial_degree(num);
    int deg_den = get_polynomial_degree(den);

    if (deg_num < 0 || deg_den < 0) return nullptr;

    if (deg_num < deg_den) {
        return lamina::detail::make_node<NumberNode>(BigInt(0));
    }

    if (deg_num == deg_den) {
        auto lc_num = get_leading_coefficient(num);
        auto lc_den = get_leading_coefficient(den);
        if (!lc_num || !lc_den) return nullptr;

        auto ratio = lamina::detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>{
            lc_num, lamina::detail::make_node<PowerNode>(lc_den, lamina::detail::make_node<NumberNode>(BigInt(-1)))
        });
        NormalizationVisitor norm;
        ratio->accept(norm);
        return norm.get_result();
    }

    /// deg_num > deg_den → ±∞
    auto lc_num = get_leading_coefficient(num);
    auto lc_den = get_leading_coefficient(den);
    int sign_num = lc_num ? get_sign(lc_num) : 1;
    int sign_den = lc_den ? get_sign(lc_den) : 1;
    int final_sign = sign_num * sign_den;

    std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
    auto inf_node = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
    if (final_sign < 0) {
        std::vector<std::shared_ptr<const SymbolicNode>> m = {lamina::detail::make_node<NumberNode>(BigInt(-1)), inf_node};
        return lamina::detail::make_node<MultiplyNode>(m);
    }
    return inf_node;
}

/**
 * @brief 分类表达式的增长速率。
 *
 * 增长速率层次：Exponential > Polynomial > Logarithmic > Constant
 *
 * @param[in] node AST 节点
 * @return 增长速率分类
 */
LimitVisitor::GrowthClass LimitVisitor::classify_growth(const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) return GrowthClass::Unknown;

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return GrowthClass::Constant;
    }

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        return (v->name() == var) ? GrowthClass::Polynomial : GrowthClass::Constant;
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (func->type() == FunctionNode::FuncType::Exp) {
            /// Check if argument contains var with positive growth
            if (func->arguments().size() == 1) {
                int arg_deg = get_polynomial_degree(func->arguments()[0]);
                if (arg_deg > 0) return GrowthClass::Exponential;
            }
            return GrowthClass::Exponential;
        }
        if (func->type() == FunctionNode::FuncType::Ln || func->type() == FunctionNode::FuncType::Log) {
            return GrowthClass::Logarithmic;
        }
        if (func->type() == FunctionNode::FuncType::Infinity) {
            return GrowthClass::Unknown;
        }
        /// sin, cos, etc. are bounded
        if (func->type() == FunctionNode::FuncType::Sin || func->type() == FunctionNode::FuncType::Cos) {
            return GrowthClass::Constant;
        }
        return GrowthClass::Unknown;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        /// x^n is polynomial, e^x is exponential
        auto base_var = std::dynamic_pointer_cast<const VariableNode>(pow->base());
        if (base_var && base_var->name() == var) {
            return GrowthClass::Polynomial;
        }
        /// Check if base is exp-like: a^x where a > 1
        if (auto base_num = std::dynamic_pointer_cast<const NumberNode>(pow->base())) {
            int exp_deg = get_polynomial_degree(pow->exponent());
            if (exp_deg > 0) return GrowthClass::Exponential;
        }
        /// Check if base contains exp
        auto base_growth = classify_growth(pow->base());
        if (base_growth == GrowthClass::Exponential) return GrowthClass::Exponential;
        if (base_growth == GrowthClass::Polynomial) return GrowthClass::Polynomial;
        return GrowthClass::Unknown;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        GrowthClass max_class = GrowthClass::Constant;
        for (auto& op : mul->operands()) {
            GrowthClass g = classify_growth(op);
            if (g == GrowthClass::Unknown) return GrowthClass::Unknown;
            if (static_cast<int>(g) > static_cast<int>(max_class)) max_class = g;
        }
        return max_class;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        GrowthClass max_class = GrowthClass::Constant;
        for (auto& op : add->operands()) {
            GrowthClass g = classify_growth(op);
            if (g == GrowthClass::Unknown) return GrowthClass::Unknown;
            if (static_cast<int>(g) > static_cast<int>(max_class)) max_class = g;
        }
        return max_class;
    }

    return GrowthClass::Unknown;
}

/**
 * @brief 获取多项式增长的有效次数（用于增长速率比较）。
 *
 * 对于 x^n 返回 n，对于 x^n * ln(x)^m 返回 n（多项式部分主导）。
 *
 * @param[in] node AST 节点
 * @return 有效多项式次数
 */
int LimitVisitor::get_growth_polynomial_degree(const std::shared_ptr<const SymbolicNode>& node) const {
    int deg = get_polynomial_degree(node);
    if (deg >= 0) return deg;

    /// For multiply nodes with mixed polynomial and logarithmic factors
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        int total_poly_deg = 0;
        for (auto& op : mul->operands()) {
            GrowthClass g = classify_growth(op);
            if (g == GrowthClass::Polynomial || g == GrowthClass::Constant) {
                int d = get_polynomial_degree(op);
                if (d >= 0) total_poly_deg += d;
            }
            /// Logarithmic factors don't contribute to polynomial degree
        }
        return total_poly_deg;
    }

    return 0;
}

/**
 * @brief 通过增长速率比较计算极限。
 *
 * 当分子和分母都趋向无穷时，比较增长速率：
 * - 指数 > 多项式 > 对数
 * - 同类增长时比较具体次数
 *
 * @param[in] num 分子 AST 节点
 * @param[in] den 分母 AST 节点
 * @return 极限结果；nullptr 表示增长率比较保持未知
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::limit_by_growth_comparison(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den) {

    GrowthClass num_growth = classify_growth(num);
    GrowthClass den_growth = classify_growth(den);

    if (num_growth == GrowthClass::Unknown || den_growth == GrowthClass::Unknown) {
        return nullptr;
    }

    /// Different growth classes
    if (static_cast<int>(num_growth) > static_cast<int>(den_growth)) {
        /// Numerator grows faster → ±∞
        std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
        return lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
    }

    if (static_cast<int>(num_growth) < static_cast<int>(den_growth)) {
        /// Denominator grows faster → 0
        return lamina::detail::make_node<NumberNode>(BigInt(0));
    }

    /// Same growth class — compare within class
    if (num_growth == GrowthClass::Polynomial) {
        int num_deg = get_growth_polynomial_degree(num);
        int den_deg = get_growth_polynomial_degree(den);
        if (num_deg < den_deg) return lamina::detail::make_node<NumberNode>(BigInt(0));
        if (num_deg > den_deg) {
            std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
            return lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
        }
        /// Same degree — fall through to L'Hôpital or other methods
    }

    return nullptr;
}

/**
 * @brief 处理 x→-∞ 的极限，通过代换 x = -t 转化为 t→+∞。
 *
 * @param[in] expr 原始表达式
 * @return 极限结果；nullptr 表示当前规则集之外
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::handle_neg_infinity_limit(
    const std::shared_ptr<const SymbolicNode>& expr) {

    if (!expr) return nullptr;

    std::string t_var = "__neg_inf_t__";

    /// Substitute x = -t in the expression
    auto substituted = substitute_neg_t(expr, t_var);
    if (!substituted) return nullptr;

    /// Normalize the substituted expression
    NormalizationVisitor norm;
    substituted->accept(norm);
    substituted = norm.get_result();
    if (!substituted) return nullptr;

    /// Evaluate lim(t→+∞) of the substituted expression
    std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
    auto pos_inf = lamina::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);

    LimitVisitor sub_vis(t_var, pos_inf, "", assumption_ctx_);
    sub_vis.lhopital_depth_ = this->lhopital_depth_;
    substituted->accept(sub_vis);
    return sub_vis.get_result();
}

/**
 * @brief 在表达式中将 var 替换为 -t_var。
 *
 * @param[in] node AST 节点
 * @param[in] t_var 替换变量名
 * @return 替换后的节点
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::substitute_neg_t(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& t_var) const {

    if (!node) return nullptr;

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return node->clone();
    }

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        if (v->name() == var) {
            /// x = -t → multiply(-1, t)
            std::vector<std::shared_ptr<const SymbolicNode>> ops = {
                lamina::detail::make_node<NumberNode>(BigInt(-1)),
                lamina::detail::make_node<VariableNode>(t_var)
            };
            return lamina::detail::make_node<MultiplyNode>(ops);
        }
        return node->clone();
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        for (auto& op : add->operands()) {
            auto sub = substitute_neg_t(op, t_var);
            if (!sub) return nullptr;
            new_ops.push_back(sub);
        }
        return lamina::detail::make_node<AddNode>(new_ops);
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        for (auto& op : mul->operands()) {
            auto sub = substitute_neg_t(op, t_var);
            if (!sub) return nullptr;
            new_ops.push_back(sub);
        }
        return lamina::detail::make_node<MultiplyNode>(new_ops);
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto new_base = substitute_neg_t(pow->base(), t_var);
        auto new_exp = substitute_neg_t(pow->exponent(), t_var);
        if (!new_base || !new_exp) return nullptr;
        return lamina::detail::make_node<PowerNode>(new_base, new_exp);
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (func->type() == FunctionNode::FuncType::Infinity) {
            return node->clone();
        }
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        for (auto& arg : func->arguments()) {
            auto sub = substitute_neg_t(arg, t_var);
            if (!sub) return nullptr;
            new_args.push_back(sub);
        }
        return lamina::detail::make_node<FunctionNode>(func->type(), new_args);
    }

    /// For other node types, clone as-is
    return node->clone();
}
