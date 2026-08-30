#include "internal/integration_support.hpp"

namespace lamina {

namespace {

/// 判断表达式是否仅由 var 通过 sin(var)/cos(var)/tan(var) 以及常数、四则、整数幂构成，
/// 即关于 sin/cos 的有理函数。含有其它依赖 var 的函数（exp/ln/sqrt 等）时返回 false。
bool weier_is_rational_trig(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    if (!node) return true;
    if (auto vn = std::dynamic_pointer_cast<const VariableNode>(node)) {
        /// 裸 var 不允许（如 x*sin(x) 不是 sin/cos 的有理函数）
        return vn->name() != var;
    }
    if (std::dynamic_pointer_cast<const NumberNode>(node)) return true;
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (auto& op : add->operands()) if (!weier_is_rational_trig(op, var)) return false;
        return true;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (auto& op : mul->operands()) if (!weier_is_rational_trig(op, var)) return false;
        return true;
    }
    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
        /// 指数必须是不依赖 var 的整数常数
        auto en = std::dynamic_pointer_cast<const NumberNode>(pw->exponent());
        if (!en) return false;
        return weier_is_rational_trig(pw->base(), var);
    }
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        using FT = FunctionNode::FuncType;
        if ((fn->type() == FT::Sin || fn->type() == FT::Cos || fn->type() == FT::Tan ||
             fn->type() == FT::Sec || fn->type() == FT::Csc || fn->type() == FT::Cot) &&
            fn->arguments().size() == 1) {
            /// 参数必须恰为 var
            auto av = std::dynamic_pointer_cast<const VariableNode>(fn->arguments()[0]);
            if (av && av->name() == var) return true;
            /// 参数不依赖 var 时也算常数
            return !expression_depends_on_variable(fn->arguments()[0], var);
        }
        /// 其它函数：仅当不依赖 var 才允许
        for (auto& a : fn->arguments()) if (expression_depends_on_variable(a, var)) return false;
        return true;
    }
    return false;
}

/// 是否至少包含一个 sin(var)/cos(var)/tan(var)... 形式（确保确实是三角有理函数）
bool weier_has_trig_of_var(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    if (!node) return false;
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        using FT = FunctionNode::FuncType;
        if ((fn->type() == FT::Sin || fn->type() == FT::Cos || fn->type() == FT::Tan ||
             fn->type() == FT::Sec || fn->type() == FT::Csc || fn->type() == FT::Cot) &&
            fn->arguments().size() == 1) {
            auto av = std::dynamic_pointer_cast<const VariableNode>(fn->arguments()[0]);
            if (av && av->name() == var) return true;
        }
        for (auto& a : fn->arguments()) if (weier_has_trig_of_var(a, var)) return true;
        return false;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (auto& op : add->operands()) if (weier_has_trig_of_var(op, var)) return true;
        return false;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (auto& op : mul->operands()) if (weier_has_trig_of_var(op, var)) return true;
        return false;
    }
    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return weier_has_trig_of_var(pw->base(), var) || weier_has_trig_of_var(pw->exponent(), var);
    }
    return false;
}

/// 递归替换：将 sin(var)/cos(var)/tan(var)... 替换为关于 t 的有理表达式。
///   sin = 2t/(1+t²), cos = (1-t²)/(1+t²), tan = 2t/(1-t²)
std::shared_ptr<const SymbolicNode> weier_replace(const std::shared_ptr<const SymbolicNode>& node,
                                            const std::string& var, const std::string& tvar) {
    if (!node) return node;
    using FT = FunctionNode::FuncType;
    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (fn->arguments().size() == 1) {
            auto av = std::dynamic_pointer_cast<const VariableNode>(fn->arguments()[0]);
            bool is_var = av && av->name() == var;
            if (is_var) {
                auto t = SymbolicExpr::variable(tvar);
                auto one = SymbolicExpr::number(1);
                auto t2 = SymbolicExpr::power(t, SymbolicExpr::number(2));
                auto onep = SymbolicExpr::add(one, t2);                          // 1+t²
                auto onem = SymbolicExpr::add(one, SymbolicExpr::multiply(SymbolicExpr::number(-1), t2)); // 1-t²
                auto two_t = SymbolicExpr::multiply(SymbolicExpr::number(2), t);
                switch (fn->type()) {
                    case FT::Sin: return lamina::detail::node(SymbolicExpr::divide(two_t, onep));
                    case FT::Cos: return lamina::detail::node(SymbolicExpr::divide(onem, onep));
                    case FT::Tan: return lamina::detail::node(SymbolicExpr::divide(two_t, onem));
                    case FT::Csc: return lamina::detail::node(SymbolicExpr::divide(onep, two_t));
                    case FT::Sec: return lamina::detail::node(SymbolicExpr::divide(onep, onem));
                    case FT::Cot: return lamina::detail::node(SymbolicExpr::divide(onem, two_t));
                    default: break;
                }
            }
        }
        /// 其它函数：递归替换参数
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        for (auto& a : fn->arguments()) new_args.push_back(weier_replace(a, var, tvar));
        return lamina::detail::make_node<FunctionNode>(fn->type(), new_args);
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> ops;
        for (auto& op : add->operands()) ops.push_back(weier_replace(op, var, tvar));
        return lamina::detail::make_node<AddNode>(ops);
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> ops;
        for (auto& op : mul->operands()) ops.push_back(weier_replace(op, var, tvar));
        return lamina::detail::make_node<MultiplyNode>(ops);
    }
    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return lamina::detail::make_node<PowerNode>(weier_replace(pw->base(), var, tvar),
                                           weier_replace(pw->exponent(), var, tvar));
    }
    return node->clone();
}

/// 将（已做 sin/cos→t 代换的）表达式树递归求值为有理数对 (分子, 分母)，
/// 分子分母均为关于 t 的多项式表达式。返回 nullopt 表示遇到无法处理的结构。
typedef std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>> RatPair;

std::optional<RatPair> weier_to_rational(const std::shared_ptr<const SymbolicNode>& node,
                                         const std::string& tvar) {
    auto one = SymbolicExpr::number(1);
    if (!node) return RatPair{SymbolicExpr::number(0), one};

    if (std::dynamic_pointer_cast<const NumberNode>(node) ||
        std::dynamic_pointer_cast<const VariableNode>(node)) {
        return RatPair{lamina::detail::make_expression_ptr(node->clone()), one};
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        /// 累加：a/b + c/d = (a*d + c*b)/(b*d)
        std::shared_ptr<SymbolicExpr> num = SymbolicExpr::number(0);
        std::shared_ptr<SymbolicExpr> den = one;
        for (auto& op : add->operands()) {
            auto r = weier_to_rational(op, tvar);
            if (!r) return std::nullopt;
            auto [n2, d2] = *r;
            auto new_num = SymbolicExpr::add(
                SymbolicExpr::multiply(num, d2), SymbolicExpr::multiply(n2, den));
            den = SymbolicExpr::multiply(den, d2)->simplify();
            num = new_num->simplify();
        }
        return RatPair{num, den};
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::shared_ptr<SymbolicExpr> num = one;
        std::shared_ptr<SymbolicExpr> den = one;
        for (auto& op : mul->operands()) {
            auto r = weier_to_rational(op, tvar);
            if (!r) return std::nullopt;
            auto [n2, d2] = *r;
            num = SymbolicExpr::multiply(num, n2)->simplify();
            den = SymbolicExpr::multiply(den, d2)->simplify();
        }
        return RatPair{num, den};
    }
    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto en = std::dynamic_pointer_cast<const NumberNode>(pw->exponent());
        if (!en) return std::nullopt;
        long long e;
        if (std::holds_alternative<BigInt>(en->value())) e = (long long)std::get<BigInt>(en->value()).to_int();
        else if (std::holds_alternative<lmmc_real_t>(en->value())) {
            double d = std::get<lmmc_real_t>(en->value());
            if (d != (long long)d) return std::nullopt;
            e = (long long)d;
        } else if (std::holds_alternative<Rational>(en->value())) {
            double d = std::get<Rational>(en->value()).to_double();
            if (d != (long long)d) return std::nullopt;
            e = (long long)d;
        } else return std::nullopt;

        auto base = weier_to_rational(pw->base(), tvar);
        if (!base) return std::nullopt;
        auto [bn, bd] = *base;
        bool neg = e < 0;
        long long k = neg ? -e : e;
        if (k > 32) return std::nullopt; // 防止指数爆炸
        std::shared_ptr<SymbolicExpr> num = one, den = one;
        for (long long i = 0; i < k; ++i) {
            num = SymbolicExpr::multiply(num, bn)->simplify();
            den = SymbolicExpr::multiply(den, bd)->simplify();
        }
        if (neg) std::swap(num, den);
        return RatPair{num, den};
    }
    /// 其它结构（含未代换的函数）无法表示为 t 的有理函数
    return std::nullopt;
}

} // anonymous namespace

std::shared_ptr<SymbolicExpr> WeierstrassStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator&,
    ComputationContext& computation, int depth) {
    if (!lamina::detail::node(expr)) return nullptr;

    /// 必须确实含有 sin/cos(var) 且整体为其有理函数
    if (!weier_has_trig_of_var(lamina::detail::node(expr), var)) return nullptr;
    if (!weier_is_rational_trig(lamina::detail::node(expr), var)) return nullptr;

    const std::string tvar = "__weier_t";

    /// 替换 sin/cos -> t 的有理式，并乘以 dx = 2/(1+t²) dt
    auto replaced = lamina::detail::make_expression_ptr(weier_replace(lamina::detail::node(expr), var, tvar));
    auto t = SymbolicExpr::variable(tvar);
    auto t2 = SymbolicExpr::power(t, SymbolicExpr::number(2));
    auto onep = SymbolicExpr::add(SymbolicExpr::number(1), t2);
    auto dx = SymbolicExpr::divide(SymbolicExpr::number(2), onep); // 2/(1+t²)
    auto integrand_raw = SymbolicExpr::multiply(replaced, dx);

    /// 关键：被积函数是 sin/cos 的有理函数，代换后仍是 t 的有理函数，但
    /// 嵌套分式 simplify() 无法约化。用「有理数对 (分子多项式, 分母多项式)」
    /// 递归求值整棵表达式树，得到干净的 N(t)/D(t)，再交给有理函数积分。
    auto rat = weier_to_rational(lamina::detail::node(integrand_raw), tvar);
    if (!rat) return nullptr;  // 出现非多项式结构，放弃
    auto [num_poly, den_poly] = *rat;
    if (!den_poly || lamina::detail::node(den_poly)->is_zero()) return nullptr;

    /// 用多项式 GCD 约简 num/den，得到最简有理函数，避免 simplify() 把
    /// 单一分式重新展开成分式之和（那样 RationalDecomposition 无法识别）。
    std::shared_ptr<SymbolicExpr> integrand_t;
    try {
        Polynomial<Rational> Np = symbolic_to_poly<Rational>(num_poly->expand(), tvar);
        Polynomial<Rational> Dp = symbolic_to_poly<Rational>(den_poly->expand(), tvar);
        if (Dp.degree() < 0) { integrand_t = SymbolicExpr::divide(num_poly, den_poly); }
        else {
            auto g = Polynomial<Rational>::gcd(Np, Dp);
            if (g.degree() >= 1) {
                auto [q1, r1] = Np.div_mod(g);
                auto [q2, r2] = Dp.div_mod(g);
                if (r1.is_zero() && r2.is_zero()) { Np = q1; Dp = q2; }
            }
            auto np_expr = poly_to_symbolic(Np);
            auto dp_expr = poly_to_symbolic(Dp);
            integrand_t = SymbolicExpr::divide(np_expr, dp_expr);
        }
    } catch (...) {
        integrand_t = SymbolicExpr::divide(num_poly, den_poly)->simplify();
    }

    /// 用有理分解策略直接积分（用独立的 Integrator 实例，避免污染外层
    /// 积分器的循环检测状态导致结果被错误改写）。
    Integrator inner;
    RationalDecompositionStrategy rds;
    auto rational_attempt = rds.try_integrate(
        *integrand_t, tvar, inner, computation, depth + 1);
    if (!rational_attempt) {
        throw detail::ResultPropagation(rational_attempt.error());
    }
    std::shared_ptr<SymbolicExpr> integrated;
    if (auto* candidate =
            std::get_if<IntegrationCandidate>(&rational_attempt.value())) {
        integrated = candidate->expression;
    }
    if (!integrated) {
        integrated = detail::propagate_result(
            inner.integrate_recursive(*integrand_t, tvar, computation, 0));
    }
    if (!integrated) return nullptr;

    /// 若结果仍含未求值积分节点，视为失败
    if (expression_depends_on_variable(lamina::detail::node(integrated), tvar) &&
        lamina::detail::contains_node_type<IntegralNode>(
            lamina::detail::node(integrated))) {
        return nullptr;
    }

    /// 回代 t = tan(x/2)
    auto half_x = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)),
                                         SymbolicExpr::variable(var));
    auto tan_half = lamina::detail::make_expression_ptr(lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::Tan, std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(half_x)}));
    auto result = integrated->substitute(tvar, tan_half);
    if (!result) return nullptr;
    return result->simplify();
}

} // namespace lamina
