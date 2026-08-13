#include "internal/integration_support.hpp"

namespace lamina {

namespace {

/// 检测形如 (c0 + c2*x²)^p 的二次根式幂，其中 p = ±1/2。
/// 返回 a²=|c0/c2| 信息与符号模式。pattern:
///   1 => a²-x²  (c2<0, c0>0, 归一化后 a²-x²)
///   2 => a²+x²  (c2>0, c0>0)
///   3 => x²-a²  (c2>0, c0<0)
/// 同时输出 a_sq=|c0|/|c2| 以及 c2 的绝对值（要求 c2=±1 以保持简单）。
struct QuadRadical {
    int pattern = 0;       // 1,2,3
    double exponent = 0;   // +0.5 or -0.5
    double a_sq = 0;       // a²
};

bool trigsub_match_radical(const std::shared_ptr<const SymbolicNode>& node, const std::string& var,
                           QuadRadical& out) {
    auto pw = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!pw) return false;
    auto en = std::dynamic_pointer_cast<const NumberNode>(pw->exponent());
    if (!en) return false;
    double e;
    if (std::holds_alternative<lmmc_real_t>(en->value())) e = std::get<lmmc_real_t>(en->value());
    else if (std::holds_alternative<Rational>(en->value())) e = std::get<Rational>(en->value()).to_double();
    else if (std::holds_alternative<BigInt>(en->value())) e = std::get<BigInt>(en->value()).to_double();
    else return false;
    if (std::abs(e - 0.5) > 1e-9 && std::abs(e + 0.5) > 1e-9) return false;

    /// base 必须是 c0 + c2*x²（关于 var 的二次、无一次项）
    auto base = lamina::detail::expression_from_node(pw->base());
    auto b = base.expand();
    if (!b) b = lamina::detail::make_expression_ptr(pw->base());
    /// 提取关于 var 的系数：c0（常数）、c1（一次）、c2（二次）
    /// 用求导法：c2 = (1/2) d²/dx² ; c1 = d/dx |_{x=0} ; c0 = base|_{x=0}
    auto d1 = b->differentiate(var);
    auto d2 = d1->differentiate(var);
    auto zero = SymbolicExpr::number(0);
    auto c0e = b->substitute(var, zero)->simplify();
    auto c1e = d1->substitute(var, zero)->simplify();
    auto c2e = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1,2)), d2)->simplify();
    /// 必须 c1=0，且 c2 为非零常数，c0 常数，且 d2 不依赖 var（纯二次）
    if (!lamina::detail::node(c1e) || !lamina::detail::node(c1e)->is_zero()) return false;
    if (expression_depends_on_variable(lamina::detail::node(c2e), var)) return false;
    if (expression_depends_on_variable(lamina::detail::node(c0e), var)) return false;
    if (!lamina::detail::node(c2e)->is_number() || !lamina::detail::node(c0e)->is_number()) return false;
    auto c0_checked = try_checked_numeric_constant(*c0e);
    auto c2_checked = try_checked_numeric_constant(*c2e);
    if (!c0_checked || !c2_checked) return false;
    double c0 = *c0_checked;
    double c2 = *c2_checked;
    if (std::abs(c2) < 1e-12) return false;
    /// 仅支持 c2 = ±1（标准型 a²±x² / x²-a²）
    if (std::abs(std::abs(c2) - 1.0) > 1e-9) return false;

    out.exponent = e;
    if (c2 < 0 && c0 > 0) { out.pattern = 1; out.a_sq = c0; }          // a² - x²
    else if (c2 > 0 && c0 > 0) { out.pattern = 2; out.a_sq = c0; }     // a² + x²
    else if (c2 > 0 && c0 < 0) { out.pattern = 3; out.a_sq = -c0; }    // x² - a²
    else return false;
    return true;
}

} // anonymous namespace

std::shared_ptr<SymbolicExpr> TrigSubstitutionStrategy::try_integrate(
    const SymbolicExpr& expr, const std::string& var, Integrator&,
    ComputationContext&, int) {
    if (!lamina::detail::node(expr)) return nullptr;

    auto x = SymbolicExpr::variable(var);

    /// 仅处理被积函数恰为单个二次根式幂的常见闭式情形：
    ///   ∫ (a²-x²)^(-1/2) dx = arcsin(x/a)
    ///   ∫ (a²+x²)^(-1/2) dx = arcsinh(x/a) = ln(x + √(x²+a²))
    ///   ∫ (x²-a²)^(-1/2) dx = arccosh(x/a) = ln(x + √(x²-a²))
    ///   ∫ (a²-x²)^( 1/2) dx = (x/2)√(a²-x²) + (a²/2)arcsin(x/a)
    QuadRadical qr;
    if (!trigsub_match_radical(lamina::detail::node(expr), var, qr)) return nullptr;

    double a_sq = qr.a_sq;
    auto a_sq_expr = SymbolicExpr::number(a_sq);
    auto a_expr = SymbolicExpr::sqrt(a_sq_expr); // a = √(a²)

    /// 构造 √(模式) 表达式
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto make_radicand = [&]() -> std::shared_ptr<SymbolicExpr> {
        if (qr.pattern == 1) // a²-x²
            return SymbolicExpr::add(a_sq_expr, SymbolicExpr::multiply(SymbolicExpr::number(-1), x_sq));
        if (qr.pattern == 2) // a²+x²
            return SymbolicExpr::add(a_sq_expr, x_sq);
        return SymbolicExpr::add(x_sq, SymbolicExpr::multiply(SymbolicExpr::number(-1), a_sq_expr)); // x²-a²
    };
    auto x_over_a = SymbolicExpr::divide(x, a_expr);

    auto arcsin = [&](const std::shared_ptr<SymbolicExpr>& u) {
        return lamina::detail::make_expression_ptr(lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::ArcSin, std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(u)}));
    };
    auto ln = [&](const std::shared_ptr<SymbolicExpr>& u) { return SymbolicExpr::ln(u); };

    if (qr.exponent < 0) {
        /// (...)^(-1/2)
        if (qr.pattern == 1) {
            /// arcsin(x/a)
            return arcsin(x_over_a)->simplify();
        }
        if (qr.pattern == 2) {
            /// ln(x + √(x²+a²))
            auto rad = SymbolicExpr::sqrt(SymbolicExpr::add(x_sq, a_sq_expr));
            return ln(SymbolicExpr::add(x, rad))->simplify();
        }
        /// pattern 3: ln(x + √(x²-a²))
        auto rad = SymbolicExpr::sqrt(SymbolicExpr::add(x_sq,
            SymbolicExpr::multiply(SymbolicExpr::number(-1), a_sq_expr)));
        return ln(SymbolicExpr::add(x, rad))->simplify();
    } else {
        /// (...)^(+1/2)
        auto rad = SymbolicExpr::sqrt(make_radicand());
        auto half = SymbolicExpr::number(Rational(1, 2));
        if (qr.pattern == 1) {
            /// (x/2)√(a²-x²) + (a²/2)arcsin(x/a)
            auto term1 = SymbolicExpr::multiply(SymbolicExpr::multiply(half, x), rad);
            auto term2 = SymbolicExpr::multiply(SymbolicExpr::multiply(half, a_sq_expr), arcsin(x_over_a));
            return SymbolicExpr::add(term1, term2)->simplify();
        }
        if (qr.pattern == 2) {
            /// (x/2)√(a²+x²) + (a²/2)ln(x+√(x²+a²))
            auto lrad = SymbolicExpr::sqrt(SymbolicExpr::add(x_sq, a_sq_expr));
            auto term1 = SymbolicExpr::multiply(SymbolicExpr::multiply(half, x), rad);
            auto term2 = SymbolicExpr::multiply(SymbolicExpr::multiply(half, a_sq_expr),
                ln(SymbolicExpr::add(x, lrad)));
            return SymbolicExpr::add(term1, term2)->simplify();
        }
        /// pattern 3: (x/2)√(x²-a²) - (a²/2)ln(x+√(x²-a²))
        auto lrad = SymbolicExpr::sqrt(SymbolicExpr::add(x_sq,
            SymbolicExpr::multiply(SymbolicExpr::number(-1), a_sq_expr)));
        auto term1 = SymbolicExpr::multiply(SymbolicExpr::multiply(half, x), rad);
        auto term2 = SymbolicExpr::multiply(
            SymbolicExpr::multiply(SymbolicExpr::multiply(SymbolicExpr::number(-1), half), a_sq_expr),
            ln(SymbolicExpr::add(x, lrad)));
        return SymbolicExpr::add(term1, term2)->simplify();
    }
}

} // namespace lamina
