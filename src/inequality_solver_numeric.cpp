#include "inequality_solver.hpp"
#include "numeric_evaluation.hpp"
#include "poly_utils.hpp"
#include "solve_polynomial.hpp"
#include "solve_strategies.hpp"
#include "symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/numeric_probe.hpp"
#include "internal/inequality_solver_support.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace lamina::detail::inequality_support {

std::optional<double> try_checked_numeric_constant(const SymbolicExpr& expr) {
    return lamina::detail::try_finite_numeric(expr);
}

int exact_numeric_sign(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !lamina::detail::node(expr)) return 0;
    auto simplified = expr->simplify();
    if (!simplified || !lamina::detail::node(simplified)) return 0;
    auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(simplified));
    if (!num) return 0;
    if (std::holds_alternative<BigInt>(num->value())) {
        const auto& value = std::get<BigInt>(num->value());
        if (value.is_zero()) return 0;
        return value.IsNegative() ? -1 : 1;
    }
    if (std::holds_alternative<Rational>(num->value())) {
        const auto& value = std::get<Rational>(num->value());
        if (value.get_numerator().is_zero()) return 0;
        return value.get_numerator().IsNegative() ? -1 : 1;
    }
    const auto value = std::get<lmmc_real_t>(num->value());
    if (!std::isfinite(value) || value == 0.0) return 0;
    return value < 0 ? -1 : 1;
}

int determine_leading_sign(const Polynomial<SymbolicPolyCoeff>& poly) {
    if (poly.is_zero()) return 0;
    auto lc = poly.lead_coeff().val;
    if (!lc) return 1;
    auto simplified = lc->simplify();
    if (!simplified) return 1;

    if (auto val = try_checked_numeric_constant(*simplified)) {
        if (*val > 0) return 1;
        if (*val < 0) return -1;
    }

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(simplified))) {
        if (std::holds_alternative<BigInt>(num->value())) {
            return std::get<BigInt>(num->value()).IsNegative() ? -1 : 1;
        }
        if (std::holds_alternative<Rational>(num->value())) {
            return std::get<Rational>(num->value()).get_numerator().IsNegative() ? -1 : 1;
        }
        if (std::holds_alternative<lmmc_real_t>(num->value())) {
            return std::get<lmmc_real_t>(num->value()) < 0 ? -1 : 1;
        }
    }
    return 1;
}

static std::shared_ptr<SymbolicExpr> snap_verified_integer_root(
    const Polynomial<Rational>& poly,
    const std::shared_ptr<SymbolicExpr>& root) {
    if (!root || poly.is_zero()) return root;
    auto numeric = try_checked_numeric_constant(*root);
    if (!numeric || !std::isfinite(*numeric)) return root;

    double rounded = std::round(*numeric);
    if (std::abs(*numeric - rounded) > 1e-8) return root;
    if (rounded < static_cast<double>(std::numeric_limits<long long>::min()) ||
        rounded > static_cast<double>(std::numeric_limits<long long>::max())) {
        return root;
    }

    Rational candidate(BigInt(static_cast<long long>(rounded)));
    if (poly.eval(candidate) == Rational(0)) {
        return SymbolicExpr::number(candidate);
    }
    return root;
}

std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>> find_roots_with_multiplicity(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& variable) {

    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>> result;

    auto poly_rat = symbolic_to_poly<Rational>(expr, variable);
    if (poly_rat.is_zero() || poly_rat.degree() <= 0) {
        return result;
    }

    auto factors = square_free_factorization(poly_rat);
    if (factors.empty()) {

        auto poly_spc = symbolic_to_poly<SymbolicPolyCoeff>(expr, variable);
        if (!poly_spc.is_zero() && poly_spc.degree() >= 1) {
            auto roots = solve_by_factoring(poly_spc, variable);
            for (const auto& root : roots) {
                if (!root) continue;
                auto verified_root = snap_verified_integer_root(poly_rat, root);
                if (try_checked_numeric_constant(*verified_root)) {
                    result.push_back({verified_root, 1});
                }
            }
        }
        return result;
    }

    for (const auto& [factor, mult] : factors) {
        if (factor.degree() <= 0) continue;

        if (factor.degree() == 1) {

            Rational a = factor.coeffs[1];
            Rational b = factor.coeffs[0];
            if (a != Rational(0)) {
                Rational root_val = Rational(0) - b / a;
                auto root_expr = SymbolicExpr::number(root_val);
                if (try_checked_numeric_constant(*root_expr)) {
                    result.push_back({root_expr, mult});
                }
            }
            continue;
        }

        if (factor.degree() == 2) {
            Rational a = factor.coeffs[2];
            Rational b = factor.coeffs[1];
            Rational c = factor.coeffs[0];

            Rational disc = b * b - Rational(4) * a * c;
            double disc_val = disc.to_double();
            if (disc_val < -1e-10) continue;
            if (disc_val < 0) disc_val = 0;

            double a_val = a.to_double();
            double b_val = b.to_double();
            double sqrt_disc = std::sqrt(disc_val);

            double r1 = (-b_val + sqrt_disc) / (2.0 * a_val);
            double r2 = (-b_val - sqrt_disc) / (2.0 * a_val);

            if (std::isfinite(r1)) {
                auto root = snap_verified_integer_root(factor, SymbolicExpr::number(r1));
                result.push_back({root, mult});
            }
            if (std::isfinite(r2) && std::abs(r1 - r2) > 1e-10) {
                auto root = snap_verified_integer_root(factor, SymbolicExpr::number(r2));
                result.push_back({root, mult});
            }
            continue;
        }

        std::vector<SymbolicPolyCoeff> spc_coeffs;
        for (int i = 0; i <= factor.degree(); ++i) {
            spc_coeffs.push_back(SymbolicPolyCoeff(SymbolicExpr::number(factor.coeffs[i])));
        }
        Polynomial<SymbolicPolyCoeff> factor_spc(spc_coeffs, variable);

        auto factor_roots = solve_by_factoring(factor_spc, variable);
        for (const auto& root : factor_roots) {
            if (!root) continue;
            auto verified_root = snap_verified_integer_root(factor, root);
            if (try_checked_numeric_constant(*verified_root)) {
                result.push_back({verified_root, mult});
            }
        }
    }

    return result;
}

bool root_less_than(const std::shared_ptr<SymbolicExpr>& a,
                           const std::shared_ptr<SymbolicExpr>& b) {
    if (!a || !b) return false;
    auto va = try_checked_numeric_constant(*a);
    auto vb = try_checked_numeric_constant(*b);
    if (va && vb) {
        return *va < *vb;
    }

    /// 参数根的顺序通过 (a - b) 的符号推导;
    /// 二次公式的根差可化简为 +/-sqrt(disc)/a.
    auto diff = SymbolicExpr::add(a, SymbolicExpr::multiply(b, SymbolicExpr::number(-1)))->simplify();
    if (auto vd = try_checked_numeric_constant(*diff)) {
        return *vd < 0;
    }

    /// 尝试判断差值表达式的符号结构:
    /// 如果差值形如 k * sqrt(...) / denom,判断各因子的符号.
    /// 这覆盖了二次公式根差 = sqrt(delta)/a 的情形.
    auto try_sign_of_node = [](const std::shared_ptr<const SymbolicNode>& node) -> int {
        if (!node) return 0;

        if (auto num_node = std::dynamic_pointer_cast<const NumberNode>(node)) {
            if (std::holds_alternative<BigInt>(num_node->value())) {
                auto& v = std::get<BigInt>(num_node->value());
                if (v.is_zero()) return 0;
                return v.IsNegative() ? -1 : 1;
            }
            if (std::holds_alternative<Rational>(num_node->value())) {
                auto& v = std::get<Rational>(num_node->value());
                if (v.get_numerator().is_zero()) return 0;
                return v.get_numerator().IsNegative() ? -1 : 1;
            }
            if (std::holds_alternative<lmmc_real_t>(num_node->value())) {
                auto v = std::get<lmmc_real_t>(num_node->value());
                if (v == 0.0) return 0;
                return v < 0 ? -1 : 1;
            }
        }

        /// sqrt(...) 非负(假设参数使判别式非负)
        if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
            if (fn->type() == FunctionNode::FuncType::Sqrt) return 1;
        }

        /// x^(1/2) 或 x^0.5 也是平方根,非负
        if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
            auto exp_expr = lamina::detail::make_expression_ptr(pw->exponent());
            if (auto ev = try_checked_numeric_constant(*exp_expr)) {
                if (*ev > 0 && *ev < 1.0) {
                    /// base^(正分数) >= 0(假设 base 为判别式等非负量)
                    return 1;
                }
            }

            /// 对于整数指数,判断底数符号
            auto base_expr = lamina::detail::make_expression_ptr(pw->base());
            auto bv_checked = try_checked_numeric_constant(*base_expr);
            auto ev_checked = try_checked_numeric_constant(*exp_expr);
            if (bv_checked && ev_checked) {
                double bv = *bv_checked;
                double ev = *ev_checked;
                int ei = static_cast<int>(ev);
                if (std::abs(ev - ei) < 1e-10) {
                    if (bv > 0) return 1;
                    if (bv < 0) return (ei % 2 == 0) ? 1 : -1;
                }
            }
        }

        return 0;
    };

    /// 对乘积节点,各因子符号之积
    auto try_sign_of_expr = [&try_sign_of_node](const std::shared_ptr<SymbolicExpr>& expr) -> int {
        if (!expr || !lamina::detail::node(expr)) return 0;

        /// 直接节点
        int s = try_sign_of_node(lamina::detail::node(expr));
        if (s != 0) return s;

        /// 乘积:各因子符号之积
        if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
            int sign = 1;
            for (const auto& op : mul->operands()) {
                int os = try_sign_of_node(op);
                if (os == 0) return 0;
                sign *= os;
            }
            return sign;
        }

        return 0;
    };

    int diff_sign = try_sign_of_expr(diff);
    if (diff_sign < 0) return true;
    if (diff_sign > 0) return false;

    return false;
}

bool roots_equal(const std::shared_ptr<SymbolicExpr>& a,
                        const std::shared_ptr<SymbolicExpr>& b) {
    if (!a || !b) return false;
    auto va = try_checked_numeric_constant(*a);
    auto vb = try_checked_numeric_constant(*b);
    return va && vb && std::abs(*va - *vb) < 1e-10;
}
bool depends_on_any_param(const std::shared_ptr<SymbolicExpr>& expr,
                                  const std::vector<std::string>& parameters) {
    if (!expr || !lamina::detail::node(expr)) return false;
    for (const auto& param : parameters) {
        if (expression_depends_on_variable(lamina::detail::node(expr), param)) return true;
    }
    return false;
}

std::vector<std::shared_ptr<SymbolicExpr>> solve_symbolic_poly(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& variable) {

    if (poly.is_zero() || poly.degree() < 1) return {};

    int deg = poly.degree();
    auto get_coeff = [&](int d) -> std::shared_ptr<SymbolicExpr> {
        if (d < 0 || d > deg) return SymbolicExpr::number(0);
        return poly.coeffs[d].val ? poly.coeffs[d].val : SymbolicExpr::number(0);
    };

    if (deg == 1) {

        auto a = get_coeff(1);
        auto b = get_coeff(0);
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        auto root = SymbolicExpr::divide(neg_b, a)->simplify();
        return { root };
    }

    auto results = solve_by_factoring(poly, variable);
    return results;
}
} // namespace lamina::detail::inequality_support
