#include "internal/expression_transform.hpp"
#include "internal/integration_support.hpp"
#include "internal/assumption_facts.hpp"

namespace LMCAS {

namespace {

bool weier_trig_function(const FunctionNode& fn) {
    using FT = FunctionNode::FuncType;
    if (fn.arguments().size() != 1) {
        return false;
    }
    return fn.type() == FT::Sin || fn.type() == FT::Cos ||
           fn.type() == FT::Tan || fn.type() == FT::Sec ||
           fn.type() == FT::Csc || fn.type() == FT::Cot;
}

bool weier_rational_function(const FunctionNode& fn, const std::string& var) {
    if (weier_trig_function(fn)) {
        auto variable = std::dynamic_pointer_cast<const VariableNode>(fn.arguments()[0]);
        if (variable && !variable->is_constant() && variable->name() == var) { return true; }
        return !expression_depends_on_variable(fn.arguments()[0], var);
    }
    for (const auto& argument : fn.arguments()) {
        if (expression_depends_on_variable(argument, var)) { return false; }
    }
    return true;
}

/// 判断表达式是否仅由 var 通过 sin(var)/cos(var)/tan(var) 以及常数,四则,整数幂构成,
/// 即关于 sin/cos 的有理函数.含有其它依赖 var 的函数(exp/ln/sqrt 等)时返回 false.
bool weier_is_rational_trig(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    if (!node) { return true; }
    if (auto vn = std::dynamic_pointer_cast<const VariableNode>(node)) {
        /// 裸 var 表示输入超出 sin/cos 有理函数域.
        return vn->is_constant() || vn->name() != var;
    }
    if (std::dynamic_pointer_cast<const NumberNode>(node)) { return true; }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (auto& op : add->operands()) if (!weier_is_rational_trig(op, var)) { return false; }
        return true;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (auto& op : mul->operands()) if (!weier_is_rational_trig(op, var)) { return false; }
        return true;
    }
    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
        /// 指数必须是不依赖 var 的整数常数
        auto en = std::dynamic_pointer_cast<const NumberNode>(pw->exponent());
        if (!en) { return false; }
        return weier_is_rational_trig(pw->base(), var);
    }
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return weier_rational_function(*fn, var);
    }
    return false;
}

bool weier_trig_of_variable(const FunctionNode& function, const std::string& var) {
    if (!weier_trig_function(function)) {
        return false;
    }
    auto argument = std::dynamic_pointer_cast<const VariableNode>(function.arguments()[0]);
    return argument && !argument->is_constant() && argument->name() == var;
}

/// 是否至少包含一个 sin(var)/cos(var)/tan(var)... 形式(确保确实是三角有理函数)
bool weier_has_trig_of_var(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    if (!node) { return false; }
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (weier_trig_of_variable(*fn, var)) {
            return true;
        }
        for (auto& a : fn->arguments()) if (weier_has_trig_of_var(a, var)) { return true; }
        return false;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (auto& op : add->operands()) if (weier_has_trig_of_var(op, var)) { return true; }
        return false;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (auto& op : mul->operands()) if (weier_has_trig_of_var(op, var)) { return true; }
        return false;
    }
    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return weier_has_trig_of_var(pw->base(), var) || weier_has_trig_of_var(pw->exponent(), var);
    }
    return false;
}

std::shared_ptr<const SymbolicNode> weier_replace_function(
    const FunctionNode& fn, const std::string& var, const std::string& tvar) {
    using FT = FunctionNode::FuncType;
    if (fn.arguments().size() != 1) {
        return nullptr;
    }
    auto av = std::dynamic_pointer_cast<const VariableNode>(fn.arguments()[0]);
    if (!av || av->is_constant() || av->name() != var) {
        return nullptr;
    }
    auto t = SymbolicExpr::variable(tvar);
    auto one = SymbolicExpr::number(1);
    auto t2 = SymbolicExpr::power(t, SymbolicExpr::number(2));
    auto onep = SymbolicExpr::add(one, t2);                          /**< 多项式 1+t^2。 */
    auto onem = SymbolicExpr::add(one, SymbolicExpr::multiply(SymbolicExpr::number(-1), t2)); /**< 多项式 1-t^2。 */
    auto two_t = SymbolicExpr::multiply(SymbolicExpr::number(2), t);
    switch (fn.type()) {
        case FT::Sin: {
            return LMCAS::detail::node(SymbolicExpr::divide(two_t, onep));
        }
        case FT::Cos: {
            return LMCAS::detail::node(SymbolicExpr::divide(onem, onep));
        }
        case FT::Tan: {
            return LMCAS::detail::node(SymbolicExpr::divide(two_t, onem));
        }
        case FT::Csc: {
            return LMCAS::detail::node(SymbolicExpr::divide(onep, two_t));
        }
        case FT::Sec: {
            return LMCAS::detail::node(SymbolicExpr::divide(onep, onem));
        }
        case FT::Cot: {
            return LMCAS::detail::node(SymbolicExpr::divide(onem, two_t));
        }
        default: break;
    }
    return nullptr;
}

std::shared_ptr<const SymbolicNode> weier_replace(const std::shared_ptr<const SymbolicNode>& node,
                                            const std::string& var, const std::string& tvar) {
    if (!node) { return node; }
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        auto replaced = weier_replace_function(*fn, var, tvar);
        if (replaced) { return replaced; }
        /// 其它函数:递归替换参数
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        for (auto& a : fn->arguments()) new_args.push_back(weier_replace(a, var, tvar));
        return LMCAS::detail::make_node<FunctionNode>(fn->type(), new_args);
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> ops;
        for (auto& op : add->operands()) ops.push_back(weier_replace(op, var, tvar));
        return LMCAS::detail::make_node<AddNode>(ops);
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> ops;
        for (auto& op : mul->operands()) ops.push_back(weier_replace(op, var, tvar));
        return LMCAS::detail::make_node<MultiplyNode>(ops);
    }
    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return LMCAS::detail::make_node<PowerNode>(weier_replace(pw->base(), var, tvar),
                                           weier_replace(pw->exponent(), var, tvar));
    }
    return node->clone();
}

/// 递归求值已执行 sin/cos->t 代换的表达式,生成关于 t 的多项式分子与分母;
/// nullopt 表示当前节点位于有理函数支持域之外.
typedef std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>> RatPair;
std::optional<RatPair> weier_to_rational(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& tvar, ComputationContext& context);

std::optional<RatPair> weier_rational_power(const PowerNode& pw, const std::string& tvar,
    const std::shared_ptr<SymbolicExpr>& one, ComputationContext& context) {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
        auto en = std::dynamic_pointer_cast<const NumberNode>(pw.exponent());
        if (!en) { return std::nullopt; }
        const auto bounded_exponent = integration_bounded_exponent(*en, -32, 32);
        if (!bounded_exponent) { return std::nullopt; }

        auto base = weier_to_rational(pw.base(), tvar, context);
        if (!base) { return std::nullopt; }
        auto [bn, bd] = *base;
        const bool neg = *bounded_exponent < 0;
        const int power_count = neg ? -*bounded_exponent : *bounded_exponent;
        std::shared_ptr<SymbolicExpr> num = one, den = one;
        for (int i = 0; i < power_count; ++i) {
            num = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(num, bn)), transform_budget));
            den = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(den, bd)), transform_budget));
        }
        if (neg) std::swap(num, den);
        return RatPair{num, den};
}

std::optional<RatPair> weier_to_rational(const std::shared_ptr<const SymbolicNode>& node,
                                         const std::string& tvar, ComputationContext& context) {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    auto one = SymbolicExpr::number(1);
    if (!node) { return RatPair{SymbolicExpr::number(0), one}; }

    if (std::dynamic_pointer_cast<const NumberNode>(node) ||
        std::dynamic_pointer_cast<const VariableNode>(node)) {
        return RatPair{LMCAS::detail::make_expression_ptr(node->clone()), one};
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        /// 累加:a/b + c/d = (a*d + c*b)/(b*d)
        std::shared_ptr<SymbolicExpr> num = SymbolicExpr::number(0);
        std::shared_ptr<SymbolicExpr> den = one;
        for (auto& op : add->operands()) {
            auto r = weier_to_rational(op, tvar, context);
            if (!r) { return std::nullopt; }
            auto [n2, d2] = *r;
            auto new_num = SymbolicExpr::add(
                SymbolicExpr::multiply(num, d2), SymbolicExpr::multiply(n2, den));
            den = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(den, d2)), transform_budget));
            num = detail::make_expression_ptr(detail::normalize_node(detail::node(new_num), transform_budget));
        }
        return RatPair{num, den};
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::shared_ptr<SymbolicExpr> num = one;
        std::shared_ptr<SymbolicExpr> den = one;
        for (auto& op : mul->operands()) {
            auto r = weier_to_rational(op, tvar, context);
            if (!r) { return std::nullopt; }
            auto [n2, d2] = *r;
            num = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(num, n2)), transform_budget));
            den = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(den, d2)), transform_budget));
        }
        return RatPair{num, den};
    }
    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return weier_rational_power(*pw, tvar, one, context);
    }
    /// 其他节点结构由 nullopt 标记为 t 有理函数域之外.
    return std::nullopt;
}

Result<std::shared_ptr<SymbolicExpr>> weier_reduce_rational(
    const std::shared_ptr<SymbolicExpr>& num_poly,
    const std::shared_ptr<SymbolicExpr>& den_poly, const std::string& tvar, ComputationContext& context) try {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    auto numerator = symbolic_to_poly<Rational>(detail::make_expression_ptr(detail::expand_node(detail::node(num_poly), transform_budget)), tvar);
    if (!numerator && numerator.error().code != CasErrc::UnsupportedExpression) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(numerator.error());
    }
    auto denominator = symbolic_to_poly<Rational>(detail::make_expression_ptr(detail::expand_node(detail::node(den_poly), transform_budget)), tvar);
    if (!denominator && denominator.error().code != CasErrc::UnsupportedExpression) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(denominator.error());
    }
    if (!numerator || !denominator) { return SymbolicExpr::divide(num_poly, den_poly); }
    auto Np = std::move(numerator.value());
    auto Dp = std::move(denominator.value());
    if (Dp.is_zero()) { return SymbolicExpr::divide(num_poly, den_poly); }
    auto g = Polynomial<Rational>::gcd(Np, Dp);
    if (g.degree() >= 1) {
        auto [q1, r1] = Np.div_mod(g);
        auto [q2, r2] = Dp.div_mod(g);
        if (r1.is_zero() && r2.is_zero()) { Np = q1; Dp = q2; }
    }
    auto np_expr = poly_to_symbolic(Np);
    auto dp_expr = poly_to_symbolic(Dp);
    return SymbolicExpr::divide(np_expr, dp_expr);
} catch (const CasError& error) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(error);
}


}

Result<std::optional<detail::IntegralParts>> Integrator::try_weierstrass(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& computation, int depth) try {
    detail::RewriteBudget transform_budget(computation, (computation).limits().max_recursion_depth,
                                            (computation).limits().max_ast_nodes, "transform");
    if (!LMCAS::detail::node(expr)) { return std::optional<detail::IntegralParts>{}; }

    /// 必须确实含有 sin/cos(var) 且整体为其有理函数
    if (!weier_has_trig_of_var(LMCAS::detail::node(expr), var)) {
        return std::optional<detail::IntegralParts>{};
    }
    if (!weier_is_rational_trig(LMCAS::detail::node(expr), var)) {
        return std::optional<detail::IntegralParts>{};
    }

    std::string tvar = "__weier_t";
    std::size_t suffix = 0;
    while (cycle_state_.occupied_names.count(tvar)) {
        transform_budget.consume();
        tvar = "__weier_t_" + std::to_string(suffix++);
    }
    cycle_state_.occupied_names.insert(tvar);

    /// 替换 sin/cos -> t 的有理式,并乘以 dx = 2/(1+t^2) dt
    auto replaced = LMCAS::detail::make_expression_ptr(weier_replace(LMCAS::detail::node(expr), var, tvar));
    auto t = SymbolicExpr::variable(tvar);
    auto t2 = SymbolicExpr::power(t, SymbolicExpr::number(2));
    auto onep = SymbolicExpr::add(SymbolicExpr::number(1), t2);
    auto dx = SymbolicExpr::divide(SymbolicExpr::number(2), onep); // 2/(1+t^2)
    auto integrand_raw = SymbolicExpr::multiply(replaced, dx);

    /// sin/cos 有理函数经代换后仍为 t 的有理函数.
    /// 通过多项式对递归求值整棵表达式树,直接得到 N(t)/D(t).
    auto rat = weier_to_rational(LMCAS::detail::node(integrand_raw), tvar, computation);
    if (!rat) { return std::optional<detail::IntegralParts>{}; }
    auto [num_poly, den_poly] = *rat;
    if (!den_poly || LMCAS::detail::node(den_poly)->is_zero()) {
        return std::optional<detail::IntegralParts>{};
    }

    /// 使用多项式 GCD 将 num/den 化为最简有理函数,并保持单一分式结构.
    auto reduced = weier_reduce_rational(num_poly, den_poly, tvar, computation);
    if (!reduced) {
        return Result<std::optional<detail::IntegralParts>>::failure(reduced.error());
    }
    auto integrand_t = std::move(reduced.value());
    Integrator inner;
    inner.cycle_state_.occupied_names = detail::all_variable_names(
        detail::node(integrand_t), &transform_budget);
    inner.cycle_state_.occupied_names.insert(tvar);
    if (computation.assumptions()) {
        detail::AssumptionFacts::collect_names(
            *computation.assumptions(), inner.cycle_state_.occupied_names, transform_budget);
    }
    auto attempted = inner.try_rational_decomposition(
        *integrand_t, tvar, computation, depth + 1);
    if (!attempted) {
        return Result<std::optional<detail::IntegralParts>>::failure(attempted.error());
    }
    std::optional<detail::IntegralParts> integrated;
    if (attempted.value() && detail::node(attempted.value()->remainder)->is_zero() &&
        !inner.has_cycle_placeholder(attempted.value()->primitive) &&
        !inner.has_cycle_placeholder(attempted.value()->remainder)) {
        integrated = std::move(*attempted.value());
    }
    if (!integrated) {
        auto recursive = inner.integrate_recursive(*integrand_t, tvar, computation, 0);
        if (!recursive) {
            return Result<std::optional<detail::IntegralParts>>::failure(recursive.error());
        }
        integrated = std::move(recursive.value());
    }
    if (!detail::node(integrated->remainder)->is_zero() ||
        inner.has_cycle_placeholder(integrated->primitive) ||
        inner.has_cycle_placeholder(integrated->remainder) ||
        contains_unevaluated_integral(detail::node(integrated->primitive))) {
        return std::optional<detail::IntegralParts>{};
    }

    /// 回代 t = tan(x/2)
    auto half_x = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)),
                                         SymbolicExpr::variable(var));
    auto tan_half = LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Tan, std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(half_x)}));
    auto result = detail::substitute_node(
        detail::node(integrated->primitive), tvar, detail::node(tan_half), transform_budget);
    auto primitive = detail::expression_from_node(detail::normalize_node(result, transform_budget));
    if (depends_on_integration_variable(primitive, tvar)) {
        return std::optional<detail::IntegralParts>{};
    }
    return std::optional<detail::IntegralParts>{{std::move(primitive), zero_}};
} catch (const CasError& error) {
    return Result<std::optional<detail::IntegralParts>>::failure(error);
}

}
