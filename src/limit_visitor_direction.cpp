#include "internal/visitors/limit_visitor.hpp"
#include "internal/exact_constant_bounds.hpp"
#include "symbolic.hpp"
#include "internal/assumption_simplification.hpp"
#include "polynomial_conversion.hpp"
#include "internal/facts_query.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/normalization_utils.hpp"

namespace LMCAS {

std::optional<std::pair<Rational, int>> LimitVisitor::exact_polynomial_leading(
    const Polynomial<Rational>& polynomial, const Rational& point_value,
    detail::ExactBoundArithmetic& arithmetic) const {
    if (polynomial.is_zero()) { return std::make_pair(Rational(0), 0); }
    auto coefficients = polynomial.coeffs;
    for (int order = 0; !coefficients.empty(); ++order) {
        charge();
        Rational value(0);
        for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it)
            { value = arithmetic.add(arithmetic.mul(value, point_value), *it); }
        if (value != Rational(0)) { return std::make_pair(value, order); }
        for (std::size_t i = 1; i < coefficients.size(); ++i)
            { coefficients[i - 1] = arithmetic.mul(coefficients[i], Rational(BigInt(i))); }
        coefficients.pop_back();
    }
    return std::nullopt;
}

bool LimitVisitor::analytic_defined_at_point(
    const std::shared_ptr<const SymbolicNode>& expr,
    const std::shared_ptr<const SymbolicNode>& positive_base) const {
    auto substituted = detail::make_expression_ptr(expr)->substitute(var, detail::make_expression_ptr(point));
    if (!substituted) { return false; }
    std::optional<detail::AssumptionFacts> assumed;
    if (assumption_ctx_) { assumed.emplace(*assumption_ctx_); }
    const auto& facts = assumed ? static_cast<const FactsQuery&>(*assumed) : detail::no_facts();
    if (positive_base) {
        auto at_point = detail::make_expression_ptr(positive_base)->substitute(var, detail::make_expression_ptr(point));
        if (!at_point) { return false; }
        auto positive = detail::query_positive_value(detail::node(at_point), facts, *context_);
        if (!positive) { throw positive.error(); }
        if (positive.value() != Tribool::True) { return false; }
    }
    auto defined = detail::query_definedness(detail::node(substituted), facts, Domain::Real, *context_);
    if (!defined) { throw defined.error(); }
    return defined && defined.value() == Tribool::True;
}

bool LimitVisitor::analytic_power_at_point(const PowerNode& power,
    const std::shared_ptr<const SymbolicNode>& expr) const {
    std::shared_ptr<const SymbolicNode> positive_base;
    if (!analytic_at_point(power.base()) || !analytic_at_point(power.exponent())) { return false; }
    BigInt integer;
    auto number = std::dynamic_pointer_cast<const NumberNode>(power.exponent());
    const bool integral = number && try_get_integer_value(number, integer);
    if (integral && integer > BigInt(0)) { return true; }
    if (!integral) { positive_base = power.base(); }
    return analytic_defined_at_point(expr, positive_base);
}

bool LimitVisitor::analytic_function_at_point(const FunctionNode& function,
    const std::shared_ptr<const SymbolicNode>& expr) const {
    if (function.arguments().size() != 1 || !analytic_at_point(function.arguments()[0])) { return false; }
    switch (function.type()) {
        case FunctionNode::FuncType::Sin: { break; }
        case FunctionNode::FuncType::Cos: { break; }
        case FunctionNode::FuncType::Exp: { break; }
        case FunctionNode::FuncType::Ln: { break; }
        case FunctionNode::FuncType::Tan: { break; }
        case FunctionNode::FuncType::ArcTan: { break; }
        default: { return false; }
    }
    return analytic_defined_at_point(expr, nullptr);
}

bool LimitVisitor::analytic_operands_at_point(
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands) const {
    for (const auto& op : operands) {
        if (!analytic_at_point(op)) { return false; }
    }
    return true;
}

std::optional<int> LimitVisitor::rational_sign_at_infinity(
    const std::shared_ptr<const SymbolicNode>& expr) {
    auto fraction = extract_num_den(expr);
    auto numerator = symbolic_to_poly<Rational>(detail::make_expression_ptr(fraction.first), var);
    auto denominator = symbolic_to_poly<Rational>(detail::make_expression_ptr(fraction.second), var);
    if (!numerator && numerator.error().code != CasErrc::UnsupportedExpression) { throw numerator.error(); }
    if (!denominator && denominator.error().code != CasErrc::UnsupportedExpression) { throw denominator.error(); }
    if (numerator && denominator && !denominator.value().is_zero()) {
        if (numerator.value().is_zero()) { return 0; }
        auto coefficient = numerator.value().lead_coeff() / denominator.value().lead_coeff();
        int sign = coefficient > Rational(0) ? 1 : -1;
        if (is_limit_at_neg_infinity() && (numerator.value().degree() - denominator.value().degree()) % 2) { sign = -sign; }
        return sign;
    }
    return std::nullopt;
}

std::optional<int> LimitVisitor::product_sign_near_point(
    const MultiplyNode& product, const std::string& dir) {
    int sign = 1;
    for (const auto& op : product.operands()) {
        auto value = determine_sign_near_point(op, dir);
        if (!value) { return std::nullopt; }
        sign *= *value;
    }
    return sign;
}

std::optional<int> LimitVisitor::power_sign_near_point(
    const PowerNode& power, const std::string& dir) {
    BigInt integer;
    auto number = std::dynamic_pointer_cast<const NumberNode>(power.exponent());
    auto sign = determine_sign_near_point(power.base(), dir);
    if (number && sign && *sign && try_get_integer_value(number, integer))
        { return integer.is_odd() ? *sign : 1; }
    if (number && sign && *sign > 0) {
        auto real_exponent = extract_coeff_value<Rational>(detail::make_expression_ptr(power.exponent()));
        if (real_exponent) { return 1; }
    }
    return std::nullopt;
}

std::optional<int> LimitVisitor::analytic_sign_near_point(
    const std::shared_ptr<const SymbolicNode>& expr, const std::string& dir) {
    if (!analytic_at_point(expr)) { return std::nullopt; }
    DifferentiationVisitor differentiation(var);
    auto current = expr;
    for (int order = 0; order <= 8; ++order) {
        auto substituted = detail::make_expression_ptr(current)->substitute(var, detail::make_expression_ptr(point));
        if (!substituted) { return std::nullopt; }
        auto value = substituted->simplify();
        auto sign = value ? get_node_sign(detail::node(value)) : std::nullopt;
        if (!sign) { return std::nullopt; }
        if (*sign) {
            if (order % 2 && dir.empty()) { return std::nullopt; }
            return order % 2 && dir == "-" ? -*sign : *sign;
        }
        current->accept(differentiation);
        current = differentiation.get_result();
        if (!current) { return std::nullopt; }
    }
    return std::nullopt;
}

namespace {
std::optional<int> directional_leading_sign(
    const std::pair<Rational, int>& leading, const std::string& dir) {
    int sign = leading.first > Rational(0) ? 1 : leading.first < Rational(0) ? -1 : 0;
    if (leading.second % 2 && dir == "-") { sign = -sign; }
    if (leading.second % 2 && dir.empty()) { return std::nullopt; }
    return sign;
}

std::optional<bool> conjunctive_near_point_truth(
    std::optional<bool> left, std::optional<bool> right) {
    if ((left && !*left) || (right && !*right)) { return false; }
    if (left && right) { return *left && *right; }
    return std::nullopt;
}

std::optional<bool> disjunctive_near_point_truth(
    std::optional<bool> left, std::optional<bool> right) {
    if ((left && *left) || (right && *right)) { return true; }
    if (left && right) { return false; }
    return std::nullopt;
}
}

std::optional<bool> LimitVisitor::relation_truth_near_point(const RelationalNode& relation) {
    auto difference = detail::make_node<AddNode>(std::vector<std::shared_ptr<const SymbolicNode>>{
        relation.left(), detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>{
            detail::make_node<NumberNode>(BigInt(-1)), relation.right()})});
    auto normalized = detail::simplify_expression(detail::make_expression_ptr(difference), *context_);
    if (!normalized) { throw normalized.error(); }
    auto sign = determine_sign_near_point(detail::node(normalized.value()), direction);
    if (!sign) { return std::nullopt; }
    switch (relation.op()) {
        case RelationalNode::Op::GT: { return *sign > 0; }
        case RelationalNode::Op::GEQ: { return *sign >= 0; }
        case RelationalNode::Op::LT: { return *sign < 0; }
        case RelationalNode::Op::LEQ: { return *sign <= 0; }
        case RelationalNode::Op::EQ: { return *sign == 0; }
        case RelationalNode::Op::NEQ: { return *sign != 0; }
    }
    return std::nullopt;
}

std::optional<bool> LimitVisitor::logical_truth_near_point(const LogicalNode& logical) {
    auto left = condition_truth_near_point(logical.left());
    if (logical.op() == LogicalNode::Op::Not) {
        return left ? std::optional<bool>(!*left) : std::nullopt;
    }
    auto right = condition_truth_near_point(logical.right());
    if (logical.op() == LogicalNode::Op::And) {
        return conjunctive_near_point_truth(left, right);
    }
    if (logical.op() == LogicalNode::Op::Or) {
        return disjunctive_near_point_truth(left, right);
    }
    return std::nullopt;
}

std::optional<bool> LimitVisitor::condition_truth_near_point(
    const std::shared_ptr<const SymbolicNode>& node) {
    if (auto logical = std::dynamic_pointer_cast<const LogicalNode>(node)) {
        return logical_truth_near_point(*logical);
    }
    auto relation = std::dynamic_pointer_cast<const RelationalNode>(node);
    if (!relation) { return std::nullopt; }
    return relation_truth_near_point(*relation);
}



/**
 * @brief 处理分段函数节点的极限.
 *
 * 根据趋近方向选择满足条件的分支并计算极限;
 * 双侧结果相等时返回该值,差异时以 nullptr 表示极限未定义.
 */
void LimitVisitor::visit(const PiecewiseNode& node) {
    if (direction.empty()) {
        LimitVisitor left(var, point, "-", assumption_ctx_, context_);
        LimitVisitor right(var, point, "+", assumption_ctx_, context_);
        node.accept(left);
        node.accept(right);
        auto lv = left.get_result(), rv = right.get_result();
        result = lv && rv && lv->compare(*rv) == 0 ? lv : nullptr;
        return;
    }
    auto branch = select_branch_by_direction(node, direction);
    result = branch ? eval_limit(branch) : nullptr;
}

/**
 * @brief 根据趋近方向选择分段函数中满足条件的分支表达式.
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::select_branch_by_direction(
    const PiecewiseNode& node, const std::string& dir) {
    LimitVisitor witness(var, point, dir, assumption_ctx_, context_);
    for (const auto& branch : node.branches()) {
        if (witness.condition_holds_near_point(branch.condition)) return branch.expression;
        auto inverse = detail::make_node<LogicalNode>(branch.condition, nullptr, LogicalNode::Op::Not);
        if (!witness.condition_holds_near_point(inverse)) return nullptr;
    }
    return node.default_expr();
}


/**
 * @brief 方向感知的 sgn 函数极限计算.
 *
 * 当 sgn 的参数在趋近点为零时,根据趋近方向确定符号.
 */
std::optional<std::shared_ptr<const SymbolicNode>> LimitVisitor::evaluate_sgn_limit(
    const std::shared_ptr<const SymbolicNode>& arg) {
    LimitVisitor sv(var, point, direction, assumption_ctx_, context_);
    arg->accept(sv);
    auto al = sv.get_result();
    if (!al) {
        return std::nullopt;
    }
    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
    al->accept(norm);
    al = norm.get_result();
    if (!al->is_zero()) {
        auto s = get_node_sign(al);
        if (s) {
            return LMCAS::detail::make_node<NumberNode>(BigInt(*s));
        }
        return std::nullopt;
    }
    auto sign = determine_sign_near_point(arg, direction);
    if (!sign) return std::nullopt;
    return LMCAS::detail::make_node<NumberNode>(BigInt(*sign));
}

/**
 * @brief 方向感知的绝对值函数极限计算.
 */
std::optional<std::shared_ptr<const SymbolicNode>> LimitVisitor::evaluate_abs_limit(
    const std::shared_ptr<const SymbolicNode>& arg) {
    LimitVisitor sv(var, point, direction, assumption_ctx_, context_);
    arg->accept(sv);
    auto al = sv.get_result();
    if (!al) {
        return std::nullopt;
    }
    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
    al->accept(norm);
    al = norm.get_result();
    if (is_inf(al)) {
        std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
        return LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);
    }
    auto s = get_node_sign(al);
    if (s) {
        if (*s >= 0) {
            return al;
        }
        std::vector<std::shared_ptr<const SymbolicNode>> neg_ops = {LMCAS::detail::make_node<NumberNode>(BigInt(-1)), al};
        auto neg_result = LMCAS::detail::make_node<MultiplyNode>(neg_ops);
        neg_result->accept(norm);
        return norm.get_result();
    }
    return std::nullopt;
}


std::optional<std::pair<Rational, int>> LimitVisitor::exact_leading_at_point(
    const std::shared_ptr<const SymbolicNode>& expr) const {
    charge();
    detail::ExactBoundArithmetic arithmetic(*context_, "limit.leading_term");
    auto p = symbolic_to_poly<Rational>(detail::make_expression_ptr(expr), var);
    auto a = extract_coeff_value<Rational>(detail::make_expression_ptr(point));
    if (!p && p.error().code != CasErrc::UnsupportedExpression) { throw p.error(); }
    if (!a && a.error().code != CasErrc::UnsupportedExpression) { throw a.error(); }
    if (!a && !is_inf(point)) {
        auto names = detail::all_variable_names(expr);
        auto point_names = detail::all_variable_names(point);
        names.insert(point_names.begin(), point_names.end());
        std::string shifted_variable = "__limit_shift";
        while (names.count(shifted_variable)) { shifted_variable += "_"; }
        auto shifted = detail::make_expression_ptr(expr)->substitute(var,
            SymbolicExpr::add(detail::make_expression_ptr(point), SymbolicExpr::variable(shifted_variable)));
        if (shifted) { shifted = shifted->expand(); }
        if (!shifted) { return std::nullopt; }
        p = symbolic_to_poly<Rational>(shifted, shifted_variable);
        if (!p && p.error().code != CasErrc::UnsupportedExpression) { throw p.error(); }
        a = Result<Rational>::success(Rational(0));
    }
    if (!p || !a) { return std::nullopt; }
    return exact_polynomial_leading(p.value(), a.value(), arithmetic);
}

bool LimitVisitor::analytic_at_point(const std::shared_ptr<const SymbolicNode>& expr) const {
    if (!expr) { return false; }
    if (std::dynamic_pointer_cast<const NumberNode>(expr) ||
        std::dynamic_pointer_cast<const VariableNode>(expr)) { return true; }
    if (auto sum = std::dynamic_pointer_cast<const AddNode>(expr)) {
        return analytic_operands_at_point(sum->operands());
    }
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(expr)) {
        return analytic_operands_at_point(product->operands());
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(expr)) {
        return analytic_power_at_point(*power, expr);
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(expr)) {
        return analytic_function_at_point(*function, expr);
    }
    return false;
}

std::optional<int> LimitVisitor::determine_sign_near_point(
    const std::shared_ptr<const SymbolicNode>& expr, const std::string& dir) {
    if (!expr) { return std::nullopt; }
    if (is_limit_at_infinity() || is_limit_at_neg_infinity()) {
        auto sign = rational_sign_at_infinity(expr);
        if (sign) { return sign; }
    }
    if (auto leading = exact_leading_at_point(expr)) {
        return directional_leading_sign(*leading, dir);
    }
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(expr)) {
        return product_sign_near_point(*product, dir);
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(expr)) {
        auto sign = power_sign_near_point(*power, dir);
        if (sign) { return sign; }
    }
    if (is_limit_at_infinity() || is_limit_at_neg_infinity()) {
        auto sign = get_node_sign(eval_limit(expr));
        return sign && *sign != 0 ? sign : std::nullopt;
    }
    auto limiting_sign = get_node_sign(eval_limit(expr));
    if (limiting_sign && *limiting_sign != 0) { return limiting_sign; }
    return analytic_sign_near_point(expr, dir);
}
bool LimitVisitor::condition_holds_near_point(
    const std::shared_ptr<const SymbolicNode>& condition) {
    auto result = condition_truth_near_point(condition);
    return result && *result;
}


}
