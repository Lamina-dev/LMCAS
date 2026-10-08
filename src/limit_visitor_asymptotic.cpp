#include "symbolic.hpp"
#include "polynomial_conversion.hpp"
#include "internal/normalization_utils.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/assumption_simplification.hpp"
#include "internal/visitors/limit_visitor.hpp"

namespace LMCAS {

namespace {
bool compatible_leading_sign(const std::optional<int>& sign, const std::optional<int>& leading_sign) {
    if (!sign || !leading_sign) { return false; }
    return *sign != 0 && *sign == *leading_sign;
}
}



bool LimitVisitor::rational_denominator_supported(
    const std::shared_ptr<const SymbolicNode>& leading_denominator,
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den) const {
    std::optional<detail::AssumptionFacts> assumed;
    if (assumption_ctx_) { assumed.emplace(*assumption_ctx_); }
    const auto& facts = assumed ? static_cast<const FactsQuery&>(*assumed) : detail::no_facts();
    auto denominator_nonzero = detail::query_nonzero_value(leading_denominator, facts, Domain::Real, *context_);
    if (!denominator_nonzero) { throw denominator_nonzero.error(); }
    if (denominator_nonzero.value() != Tribool::True) { return false; }
    if (!defined_near_point(num) || !defined_near_point(den)) { return false; }
    return true;
}

std::shared_ptr<const SymbolicNode> LimitVisitor::rational_leading_limit(
    const std::shared_ptr<const SymbolicNode>& leading_numerator,
    const std::shared_ptr<const SymbolicNode>& leading_denominator, int degree_difference) {
    if (degree_difference == 0) {
        auto ratio = detail::make_node<MultiplyNode>(std::vector<std::shared_ptr<const SymbolicNode>>{
            leading_numerator, detail::make_node<PowerNode>(leading_denominator, detail::make_node<NumberNode>(BigInt(-1)))});
        NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
        ratio->accept(norm);
        return norm.get_result();
    }
    auto numerator_sign = get_node_sign(leading_numerator);
    auto denominator_sign = get_node_sign(leading_denominator);
    if (!denominator_sign || !*denominator_sign) { return nullptr; }
    if (!numerator_sign || !*numerator_sign) { return nullptr; }
    int sign = *numerator_sign * *denominator_sign;
    if (is_limit_at_neg_infinity() && degree_difference % 2) { sign = -sign; }
    return signed_infinity(sign);
}

std::optional<Rational> LimitVisitor::algebraic_sum_growth(const AddNode& sum) const {
    std::optional<Rational> maximum;
    std::optional<int> leading_sign;
    std::size_t leading_count = 0;
    bool common_sign = true;
    LimitVisitor witness(var, point, direction, assumption_ctx_, context_);
    for (const auto& term : sum.operands()) {
        if (term->is_zero()) { continue; }
        auto degree = algebraic_growth_degree(term);
        if (!degree) { return std::nullopt; }
        if (!maximum || *degree > *maximum) {
            maximum = *degree;
            leading_count = 1;
            leading_sign = witness.determine_sign_near_point(term, direction);
            common_sign = true;
        } else if (*degree == *maximum) {
            ++leading_count;
            auto sign = witness.determine_sign_near_point(term, direction);
            common_sign = common_sign && compatible_leading_sign(sign, leading_sign);
        }
    }
    if (maximum && (leading_count == 1 || common_sign)) { return maximum; }
    return std::nullopt;
}

std::optional<Rational> LimitVisitor::algebraic_product_growth(const MultiplyNode& product) const {
    Rational degree(0);
    for (const auto& operand : product.operands()) {
        auto d = algebraic_growth_degree(operand);
        if (!d) { return std::nullopt; }
        degree = degree + *d;
    }
    return degree;
}

std::optional<Rational> LimitVisitor::algebraic_power_growth(const PowerNode& power) const {
    auto exponent = extract_coeff_value<Rational>(detail::make_expression_ptr(power.exponent()));
    auto degree = algebraic_growth_degree(power.base());
    if (!exponent || !degree) { return std::nullopt; }
    auto number = std::dynamic_pointer_cast<const NumberNode>(power.exponent());
    BigInt integer;
    if (!number || !try_get_integer_value(number, integer)) {
        LimitVisitor witness(var, point, direction, assumption_ctx_, context_);
        auto sign = witness.determine_sign_near_point(power.base(), direction);
        if (!sign || *sign <= 0 || !defined_near_point(power.base())) { return std::nullopt; }
    }
    return *degree * exponent.value();
}

bool LimitVisitor::logarithmic_algebraic_bound(
    const std::shared_ptr<const SymbolicNode>& expr) const {
    auto function = std::dynamic_pointer_cast<const FunctionNode>(expr);
    if (!function || function->type() != FunctionNode::FuncType::Ln ||
        function->arguments().size() != 1) { return false; }
    const auto& argument = function->arguments()[0];
    if (!algebraic_growth_degree(argument) || !defined_near_point(argument)) { return false; }
    LimitVisitor witness(var, point, direction, assumption_ctx_, context_);
    auto sign = witness.determine_sign_near_point(argument, direction);
    return sign && *sign > 0;
}

bool LimitVisitor::polynomial_logarithmic_operands(
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands) const {
    for (const auto& op : operands) {
        if (!polynomial_logarithmic_bound(op)) { return false; }
    }
    return true;
}

bool LimitVisitor::exponential_dominates(
    const std::shared_ptr<const SymbolicNode>& exponential,
    const std::shared_ptr<const SymbolicNode>& other) {
    auto f = std::dynamic_pointer_cast<const FunctionNode>(exponential);
    if (!f || f->type() != FunctionNode::FuncType::Exp || f->arguments().size() != 1) {
        return false;
    }
    if (!proves_superlogarithmic_magnitude(f->arguments()[0]) ||
        !polynomial_logarithmic_bound(other)) { return false; }
    auto exponent_limit = eval_limit(f->arguments()[0]);
    return is_inf(exponent_limit) && !is_neg_inf(exponent_limit);
}



/**
 * @brief 判断当前极限点是否为正无穷.
 * @return 极限点为 +infinity 时返回 true
 */
bool LimitVisitor::is_limit_at_infinity() const {
    if (auto f = std::dynamic_pointer_cast<const FunctionNode>(point)) {
        return f->type() == FunctionNode::FuncType::Infinity;
    }
    return false;
}

/**
 * @brief 判断当前极限点是否为负无穷.
 *
 * 负无穷表示为 -1 * Infinity 的 MultiplyNode.
 * @return 极限点为 -infinity 时返回 true
 */
bool LimitVisitor::is_limit_at_neg_infinity() const {
    return is_neg_inf(point);
}

/**
 * @brief 转为关于 var 的有理系数多项式并求次数。
 * @param[in] node AST 节点
 * @return 多项式次数；空节点或不支持的表达式返回 -1。
 * @note 其他转换错误原样抛出。
 */
int LimitVisitor::get_polynomial_degree(const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) return -1;
    auto polynomial = symbolic_to_poly<Rational>(detail::make_expression_ptr(node), var);
    if (!polynomial && polynomial.error().code != CasErrc::UnsupportedExpression) throw polynomial.error();
    return polynomial ? polynomial.value().degree() : -1;
}

/**
 * @brief 获取多项式的首项系数.
 *
 * 对于多项式 P(x) = a_n * x^n + ... + a_0,返回 a_n.
 *
 * @param[in] node AST 节点
 * @return 首项系数节点；空节点、零多项式或不支持的表达式返回 nullptr。
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::get_leading_coefficient(
    const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) return nullptr;
    auto polynomial = symbolic_to_poly<Rational>(detail::make_expression_ptr(node), var);
    if (!polynomial && polynomial.error().code != CasErrc::UnsupportedExpression) throw polynomial.error();
    return polynomial && !polynomial.value().is_zero()
        ? detail::make_node<NumberNode>(polynomial.value().lead_coeff()) : nullptr;
}

/**
 * @brief 通过多项式次数比较计算有理函数在无穷处的极限.
 *
 * 对于 P(x)/Q(x):
 * - deg(P) < deg(Q) -> 0
 * - deg(P) = deg(Q) -> 首项系数之比
 * - deg(P) > deg(Q) -> +/-infinity，符号由首项系数及负无穷处次数差的奇偶性决定
 *
 * @param[in] num 分子 AST 节点
 * @param[in] den 分母 AST 节点
 * @return 极限结果,非有理函数时返回 nullptr
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::limit_rational_at_infinity(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den) {
    charge();
    auto numerator = symbolic_to_poly<SymbolicPolyCoeff>(detail::make_expression_ptr(num), var);
    auto denominator = symbolic_to_poly<SymbolicPolyCoeff>(detail::make_expression_ptr(den), var);
    if (!numerator && numerator.error().code != CasErrc::UnsupportedExpression) { throw numerator.error(); }
    if (!denominator && denominator.error().code != CasErrc::UnsupportedExpression) { throw denominator.error(); }
    if (!numerator || !denominator || denominator.value().is_zero()) { return nullptr; }
    auto normalized_denominator = detail::simplify_expression(
        denominator.value().lead_coeff().val, *context_);
    if (!normalized_denominator) { throw normalized_denominator.error(); }
    auto leading_denominator = detail::node(normalized_denominator.value());
    if (!rational_denominator_supported(leading_denominator, num, den)) { return nullptr; }
    if (numerator.value().is_zero() || numerator.value().degree() < denominator.value().degree())
        { return detail::make_node<NumberNode>(BigInt(0)); }
    auto normalized_numerator = detail::simplify_expression(
        numerator.value().lead_coeff().val, *context_);
    if (!normalized_numerator) { throw normalized_numerator.error(); }
    auto leading_numerator = detail::node(normalized_numerator.value());
    return rational_leading_limit(leading_numerator, leading_denominator,
        numerator.value().degree() - denominator.value().degree());
}

/**
 * @brief 分类表达式的增长速率.
 *
 * 增长速率层次:Exponential > Polynomial > Logarithmic > Constant
 *
 * @param[in] node AST 节点
 * @return 增长速率分类
 */
LimitVisitor::GrowthClass LimitVisitor::classify_growth(const std::shared_ptr<const SymbolicNode>& node) const {
    if (!node) {
        return GrowthClass::Unknown;
    }

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return GrowthClass::Constant;
    }

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        return (!v->is_constant() && v->name() == var)
            ? GrowthClass::Polynomial : GrowthClass::Constant;
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return function_growth(*func);
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return power_growth(*pow);
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return operands_growth(mul->operands());
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return operands_growth(add->operands());
    }

    return GrowthClass::Unknown;
}

LimitVisitor::GrowthClass LimitVisitor::function_growth(const FunctionNode& node) const {
    if (node.type() == FunctionNode::FuncType::Exp) {
        return GrowthClass::Exponential;
    }
    if (node.type() == FunctionNode::FuncType::Ln ||
        node.type() == FunctionNode::FuncType::Log) {
        return GrowthClass::Logarithmic;
    }
    if (node.type() == FunctionNode::FuncType::Sin ||
        node.type() == FunctionNode::FuncType::Cos) {
        return GrowthClass::Constant;
    }
    return GrowthClass::Unknown;
}

LimitVisitor::GrowthClass LimitVisitor::power_growth(const PowerNode& node) const {
    if (is_limit_variable(node.base())) {
        return GrowthClass::Polynomial;
    }
    if (std::dynamic_pointer_cast<const NumberNode>(node.base())) {
        if (get_polynomial_degree(node.exponent()) > 0) {
            return GrowthClass::Exponential;
        }
    }
    auto growth = classify_growth(node.base());
    if (growth == GrowthClass::Exponential) {
        return growth;
    }
    if (growth == GrowthClass::Polynomial) {
        return growth;
    }
    return GrowthClass::Unknown;
}

LimitVisitor::GrowthClass LimitVisitor::operands_growth(
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands) const {
    GrowthClass maximum = GrowthClass::Constant;
    for (const auto& operand : operands) {
        auto growth = classify_growth(operand);
        if (growth == GrowthClass::Unknown) {
            return growth;
        }
        if (static_cast<int>(growth) > static_cast<int>(maximum)) maximum = growth;
    }
    return maximum;
}


/**
 * @brief 求代数增长次数，供指数、多项式、对数增长率比较使用。
 * 乘积次数相加，幂次数相乘；和式须排除最高次项相消。
 * 非整数幂要求底数在邻域内有定义且为正。
 * @param[in] expr 表达式节点
 * @return 有理增长次数；无法证明时返回 std::nullopt。
 */
std::optional<Rational> LimitVisitor::algebraic_growth_degree(
    const std::shared_ptr<const SymbolicNode>& expr) const {
    charge();
    if (!expr) { return std::nullopt; }
    auto polynomial = symbolic_to_poly<Rational>(detail::make_expression_ptr(expr), var);
    if (!polynomial && polynomial.error().code != CasErrc::UnsupportedExpression) { throw polynomial.error(); }
    if (polynomial && !polynomial.value().is_zero())
        { return Rational(polynomial.value().degree()); }
    if (auto sum = std::dynamic_pointer_cast<const AddNode>(expr)) {
        auto cancelled = detail::make_expression_ptr(expr)->cancel();
        if (cancelled && detail::node(cancelled)->compare(*expr) != 0 && defined_near_point(expr))
            { return algebraic_growth_degree(detail::node(cancelled)); }
        return algebraic_sum_growth(*sum);
    }
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(expr)) {
        return algebraic_product_growth(*product);
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(expr)) {
        return algebraic_power_growth(*power);
    }
    return std::nullopt;
}
bool LimitVisitor::polynomial_logarithmic_bound(
    const std::shared_ptr<const SymbolicNode>& expr) const {
    auto polynomial = symbolic_to_poly<SymbolicPolyCoeff>(detail::make_expression_ptr(expr), var);
    if (!polynomial && polynomial.error().code != CasErrc::UnsupportedExpression) { throw polynomial.error(); }
    if (polynomial && defined_near_point(expr)) { return true; }
    if (algebraic_growth_degree(expr)) { return true; }
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(expr)) {
        return polynomial_logarithmic_operands(product->operands());
    }
    if (auto sum = std::dynamic_pointer_cast<const AddNode>(expr)) {
        return polynomial_logarithmic_operands(sum->operands());
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(expr)) {
        BigInt integer;
        auto n = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
        return n && try_get_integer_value(n, integer) && integer > BigInt(0) &&
            polynomial_logarithmic_bound(power->base());
    }
    return logarithmic_algebraic_bound(expr);
}
std::shared_ptr<const SymbolicNode> LimitVisitor::limit_by_growth_comparison(
    const std::shared_ptr<const SymbolicNode>& num,
    const std::shared_ptr<const SymbolicNode>& den) {
    if (exponential_dominates(den, num)) { return detail::make_node<NumberNode>(BigInt(0)); }
    if (exponential_dominates(num, den)) {
        auto denominator = eval_limit(den);
        auto sign = get_node_sign(denominator);
        if (sign && *sign) { return signed_infinity(*sign); }
    }
    auto log = std::dynamic_pointer_cast<const FunctionNode>(num);
    auto denominator_degree = algebraic_growth_degree(den);
    if (log && log->type() == FunctionNode::FuncType::Ln &&
        polynomial_logarithmic_bound(num) && denominator_degree && *denominator_degree > Rational(0))
        { return detail::make_node<NumberNode>(BigInt(0)); }
    return nullptr;
}

bool LimitVisitor::proves_superlogarithmic_magnitude(
    const std::shared_ptr<const SymbolicNode>& expr) const {
    auto degree = algebraic_growth_degree(expr);
    return degree && *degree > Rational(0);
}


/**
 * @brief 处理 x->-infinity 的极限,通过代换 x = -t 转化为 t->+infinity.
 *
 * @param[in] expr 原始表达式
 * @return 极限结果；不支持的表达式返回 nullptr。
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::handle_neg_infinity_limit(
    const std::shared_ptr<const SymbolicNode>& expr) {

    if (!expr) {
        return nullptr;
    }

    std::string t_var = "__neg_inf_t__";

    auto substituted = substitute_neg_t(expr, t_var);
    if (!substituted) {
        return nullptr;
    }

    NormalizationVisitor norm(*context_, detail::no_facts(), Domain::Real, rewrite_budget());
    substituted->accept(norm);
    substituted = norm.get_result();
    if (!substituted) {
        return nullptr;
    }

    std::vector<std::shared_ptr<const SymbolicNode>> inf_args;
    auto pos_inf = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Infinity, inf_args);

    LimitVisitor sub_vis(t_var, pos_inf, "", assumption_ctx_, context_);
    sub_vis.lhopital_depth_ = this->lhopital_depth_;
    substituted->accept(sub_vis);
    return sub_vis.get_result();
}

static std::shared_ptr<const SymbolicNode> substitute_negative_variable(
    const VariableNode& node, const std::string& var, const std::string& t_var) {
    if (!node.is_constant() && node.name() == var) {
        std::vector<std::shared_ptr<const SymbolicNode>> ops = {
            LMCAS::detail::make_node<NumberNode>(BigInt(-1)),
            LMCAS::detail::make_node<VariableNode>(t_var)
        };
        return LMCAS::detail::make_node<MultiplyNode>(ops);
    }
    return node.clone();
}

/**
 * @brief 在表达式中将 var 替换为 -t_var.
 *
 * @param[in] node AST 节点
 * @param[in] t_var 替换变量名
 * @return 替换后的节点
 */
std::shared_ptr<const SymbolicNode> LimitVisitor::substitute_neg_t(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& t_var) const {

    if (!node) {
        return nullptr;
    }

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return node->clone();
    }

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        return substitute_negative_variable(*v, var, t_var);
    }

    auto add = std::dynamic_pointer_cast<const AddNode>(node);
    auto mul = add ? std::shared_ptr<const MultiplyNode>{}
                   : std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (add || mul) {
        const auto& operands = add ? add->operands() : mul->operands();
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        for (const auto& op : operands) {
            auto sub = substitute_neg_t(op, t_var);
            if (!sub) {
                return nullptr;
            }
            new_ops.push_back(sub);
        }
        if (add) {
            return LMCAS::detail::make_node<AddNode>(new_ops);
        }
        return LMCAS::detail::make_node<MultiplyNode>(new_ops);
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto new_base = substitute_neg_t(pow->base(), t_var);
        auto new_exp = substitute_neg_t(pow->exponent(), t_var);
        if (!new_base || !new_exp) {
            return nullptr;
        }
        return LMCAS::detail::make_node<PowerNode>(new_base, new_exp);
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return substitute_function_neg_t(*func, t_var);
    }
    return node->clone();
}

std::shared_ptr<const SymbolicNode> LimitVisitor::substitute_function_neg_t(
    const FunctionNode& node, const std::string& t_var) const {
    if (node.type() == FunctionNode::FuncType::Infinity) {
        return node.clone();
    }
    std::vector<std::shared_ptr<const SymbolicNode>> arguments;
    for (const auto& argument : node.arguments()) {
        auto substituted = substitute_neg_t(argument, t_var);
        if (!substituted) {
            return nullptr;
        }
        arguments.push_back(substituted);
    }
    return LMCAS::detail::make_node<FunctionNode>(node.type(), arguments);
}

}
