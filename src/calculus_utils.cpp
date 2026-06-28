/**
 * @file calculus_utils.cpp
 * @brief 微积分工具函数实现：连续性判定、渐近线分析、对数微分、微分、全微分、反函数导数、反函数求解。
 */

#include "calculus_utils.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "integration.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <memory>

namespace lamina {

// ============================================================
// 内部辅助函数（连续性・渐近线）
// ============================================================

/**
 * @internal
 * @brief 判断表达式是否包含无穷大节点。
 */
static bool calculus_utils_is_infinity(const std::shared_ptr<SymbolicExpr>& expr)
{
    if (!expr || !expr->root) return false;

    if (auto f = std::dynamic_pointer_cast<FunctionNode>(expr->root)) {
        return f->type == FunctionNode::FuncType::Infinity;
    }
    if (auto m = std::dynamic_pointer_cast<MultiplyNode>(expr->root)) {
        for (auto& op : m->operands) {
            if (auto f = std::dynamic_pointer_cast<FunctionNode>(op)) {
                if (f->type == FunctionNode::FuncType::Infinity) return true;
            }
        }
    }
    return false;
}

/**
 * @internal
 * @brief 判断两个表达式是否结构相等。
 */
static bool calculus_utils_expr_equal(const std::shared_ptr<SymbolicExpr>& a,
                                      const std::shared_ptr<SymbolicExpr>& b)
{
    if (!a || !a->root || !b || !b->root) {
        return (!a || !a->root) && (!b || !b->root);
    }
    return a->root->equals(*b->root);
}

/**
 * @internal
 * @brief 从表达式中提取分母部分。
 *
 * 对于 PowerNode 且指数为负，base 即为分母。
 * 对于 MultiplyNode，查找指数为负的因子作为分母。
 */
static std::shared_ptr<SymbolicExpr> calculus_utils_extract_denominator(
    const std::shared_ptr<SymbolicExpr>& expr)
{
    if (!expr || !expr->root) return nullptr;

    // 检查是否为 PowerNode 且指数为负（如 x^(-1) 表示 1/x）
    if (auto pow = std::dynamic_pointer_cast<PowerNode>(expr->root)) {
        if (auto num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
            double e = 0;
            if (std::holds_alternative<lmmc_real_t>(num->value))
                e = std::get<lmmc_real_t>(num->value);
            else if (std::holds_alternative<BigInt>(num->value))
                e = std::get<BigInt>(num->value).to_double();
            else if (std::holds_alternative<Rational>(num->value))
                e = std::get<Rational>(num->value).to_double();

            if (e < 0) {
                if (e == -1.0) {
                    return std::make_shared<SymbolicExpr>(pow->base);
                }
                auto pos_exp = std::make_shared<NumberNode>(BigInt(static_cast<int>(-e)));
                auto den_node = std::make_shared<PowerNode>(pow->base, pos_exp);
                return std::make_shared<SymbolicExpr>(den_node);
            }
        }
    }

    // 检查 MultiplyNode 中是否有负指数因子
    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(expr->root)) {
        std::vector<std::shared_ptr<SymbolicNode>> den_factors;

        for (auto& op : mul->operands) {
            if (auto pow = std::dynamic_pointer_cast<PowerNode>(op)) {
                if (auto num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
                    double e = 0;
                    if (std::holds_alternative<lmmc_real_t>(num->value))
                        e = std::get<lmmc_real_t>(num->value);
                    else if (std::holds_alternative<BigInt>(num->value))
                        e = std::get<BigInt>(num->value).to_double();
                    else if (std::holds_alternative<Rational>(num->value))
                        e = std::get<Rational>(num->value).to_double();

                    if (e < 0) {
                        if (e == -1.0) {
                            den_factors.push_back(pow->base);
                        } else {
                            auto pos_exp = std::make_shared<NumberNode>(BigInt(static_cast<int>(-e)));
                            den_factors.push_back(std::make_shared<PowerNode>(pow->base, pos_exp));
                        }
                    }
                }
            }
        }

        if (!den_factors.empty()) {
            std::shared_ptr<SymbolicNode> den_node;
            if (den_factors.size() == 1) {
                den_node = den_factors[0];
            } else {
                den_node = std::make_shared<MultiplyNode>(den_factors);
            }
            return std::make_shared<SymbolicExpr>(den_node);
        }
    }

    return nullptr;
}

// ============================================================
// 连续性判定 (Requirement 7)
// ============================================================

ContinuityType continuity_at(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point)
{
    if (!f || !point) return ContinuityType::Essential;

    // 计算左极限
    auto left_lim = f->limit(var, point, "-");
    // 计算右极限
    auto right_lim = f->limit(var, point, "+");
    // 计算函数值（直接代入）
    auto func_val = f->substitute(var, point);

    // 简化结果
    if (left_lim) left_lim = left_lim->simplify();
    if (right_lim) right_lim = right_lim->simplify();
    if (func_val) func_val = func_val->simplify();

    // 检查极限是否存在（不为无穷）
    bool left_exists = left_lim && !calculus_utils_is_infinity(left_lim);
    bool right_exists = right_lim && !calculus_utils_is_infinity(right_lim);

    // 如果任一侧极限不存在或为无穷 → 本性间断点
    if (!left_exists || !right_exists) {
        return ContinuityType::Essential;
    }

    // 检查左极限是否等于右极限
    bool limits_equal = calculus_utils_expr_equal(left_lim, right_lim);

    if (!limits_equal) {
        // 左极限 ≠ 右极限 → 跳跃间断点
        return ContinuityType::Jump;
    }

    // 左极限 = 右极限，检查是否等于函数值
    bool val_exists = func_val && !calculus_utils_is_infinity(func_val);

    if (!val_exists || !calculus_utils_expr_equal(left_lim, func_val)) {
        // 极限存在但不等于函数值（或函数无定义） → 可去间断点
        return ContinuityType::Removable;
    }

    // 左极限 = 右极限 = 函数值 → 连续
    return ContinuityType::Continuous;
}

// ============================================================
// 渐近线分析 (Requirement 15)
// ============================================================

AsymptoteResult asymptotes(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    AsymptoteResult result;
    if (!f) return result;

    auto pos_inf = SymbolicExpr::infinity(1);
    auto neg_inf = SymbolicExpr::infinity(-1);

    // ========== 垂直渐近线 ==========
    // 提取分母，求解分母 = 0 的点
    auto denominator = calculus_utils_extract_denominator(f);
    if (denominator) {
        auto zeros = SymbolicExpr::solve(denominator, var);
        for (auto& z : zeros) {
            if (!z || calculus_utils_is_infinity(z)) continue;

            // 验证该点处极限为 ±∞
            auto lim_at_z = f->limit(var, z);
            if (calculus_utils_is_infinity(lim_at_z)) {
                result.vertical.push_back(z);
            } else {
                // 尝试单侧极限
                auto lim_right = f->limit(var, z, "+");
                auto lim_left = f->limit(var, z, "-");
                if (calculus_utils_is_infinity(lim_right) || calculus_utils_is_infinity(lim_left)) {
                    result.vertical.push_back(z);
                }
            }
        }
    }

    // ========== 水平渐近线 ==========
    // 计算 x→+∞ 和 x→-∞ 的极限
    auto lim_pos = f->limit(var, pos_inf);
    auto lim_neg = f->limit(var, neg_inf);

    if (lim_pos) lim_pos = lim_pos->simplify();
    if (lim_neg) lim_neg = lim_neg->simplify();

    bool has_horiz_pos = lim_pos && !calculus_utils_is_infinity(lim_pos);
    bool has_horiz_neg = lim_neg && !calculus_utils_is_infinity(lim_neg);

    if (has_horiz_pos) {
        result.horizontal.push_back(lim_pos);
    }
    if (has_horiz_neg) {
        // 避免重复添加相同的水平渐近线
        bool duplicate = has_horiz_pos && calculus_utils_expr_equal(lim_pos, lim_neg);
        if (!duplicate) {
            result.horizontal.push_back(lim_neg);
        }
    }

    // ========== 斜渐近线 ==========
    // 仅当对应方向无水平渐近线时才检查斜渐近线
    auto x_expr = SymbolicExpr::variable(var);

    if (!has_horiz_pos) {
        // 计算 slope = lim(f/x) as x→+∞
        auto f_over_x = SymbolicExpr::multiply(f, SymbolicExpr::power(x_expr, SymbolicExpr::number(-1)));
        auto slope_pos = f_over_x->limit(var, pos_inf);
        if (slope_pos) slope_pos = slope_pos->simplify();

        if (slope_pos && !calculus_utils_is_infinity(slope_pos) && !slope_pos->is_zero()) {
            // 计算 intercept = lim(f - slope*x) as x→+∞
            auto slope_times_x = SymbolicExpr::multiply(slope_pos, x_expr);
            auto f_minus_mx = SymbolicExpr::add(f, SymbolicExpr::multiply(SymbolicExpr::number(-1), slope_times_x));
            auto intercept_pos = f_minus_mx->limit(var, pos_inf);
            if (intercept_pos) intercept_pos = intercept_pos->simplify();

            if (intercept_pos && !calculus_utils_is_infinity(intercept_pos)) {
                result.oblique.emplace_back(slope_pos, intercept_pos);
            }
        }
    }

    if (!has_horiz_neg) {
        // 计算 slope = lim(f/x) as x→-∞
        auto f_over_x = SymbolicExpr::multiply(f, SymbolicExpr::power(x_expr, SymbolicExpr::number(-1)));
        auto slope_neg = f_over_x->limit(var, neg_inf);
        if (slope_neg) slope_neg = slope_neg->simplify();

        if (slope_neg && !calculus_utils_is_infinity(slope_neg) && !slope_neg->is_zero()) {
            // 计算 intercept = lim(f - slope*x) as x→-∞
            auto slope_times_x = SymbolicExpr::multiply(slope_neg, x_expr);
            auto f_minus_mx = SymbolicExpr::add(f, SymbolicExpr::multiply(SymbolicExpr::number(-1), slope_times_x));
            auto intercept_neg = f_minus_mx->limit(var, neg_inf);
            if (intercept_neg) intercept_neg = intercept_neg->simplify();

            if (intercept_neg && !calculus_utils_is_infinity(intercept_neg)) {
                // 避免重复：如果 +∞ 方向已有相同斜渐近线则跳过
                bool duplicate = false;
                for (auto& [s, i] : result.oblique) {
                    if (calculus_utils_expr_equal(s, slope_neg) &&
                        calculus_utils_expr_equal(i, intercept_neg)) {
                        duplicate = true;
                        break;
                    }
                }
                if (!duplicate) {
                    result.oblique.emplace_back(slope_neg, intercept_neg);
                }
            }
        }
    }

    return result;
}

// ============================================================
// 对数微分、微分、全微分、反函数 (Requirements 12, 13, 14, 77)
// ============================================================

std::shared_ptr<SymbolicExpr> log_differentiate(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    if (!f || !f->root) return nullptr;

    /// 计算 ln(f)
    auto ln_f = SymbolicExpr::ln(f);

    /// 对 ln(f) 求导 — 化简器会自动应用对数规则
    /// （ln(a*b) = ln(a)+ln(b), ln(a^n) = n*ln(a)）
    auto ln_f_simplified = ln_f->simplify();
    auto d_ln_f = ln_f_simplified->differentiate(var);

    /// 结果 = f * d/dx[ln(f)]
    auto result = SymbolicExpr::multiply(f, d_ln_f);
    return result->simplify();
}

std::shared_ptr<SymbolicExpr> differential(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    if (!f || !f->root) return nullptr;

    /// 微分 df = f'(var) * dx，返回系数 f'(var)
    return f->differentiate(var);
}

std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::string>> total_differential(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars)
{
    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::string>> result;
    if (!f || !f->root) return result;

    result.reserve(vars.size());
    for (const auto& v : vars) {
        auto partial = f->differentiate(v);
        result.emplace_back(partial, v);
    }
    return result;
}

std::shared_ptr<SymbolicExpr> inverse_derivative(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point)
{
    if (!f || !f->root || !point || !point->root) return nullptr;

    /// 求解 f(x) = point，即 f(x) - point = 0
    auto eq = SymbolicExpr::add(f,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), point));
    auto solutions = SymbolicExpr::solve(eq, var);

    if (solutions.empty()) return nullptr;

    /// 取第一个解作为 f⁻¹(point)
    auto x0 = solutions[0];

    /// 计算 f'(var)
    auto f_prime = f->differentiate(var);

    /// 在 x = x0 处求值 f'(x0)
    auto f_prime_at_x0 = f_prime->substitute(var, x0)->simplify();

    /// 若 f'(x0) = 0，反函数导数不存在
    if (f_prime_at_x0->is_zero()) return nullptr;

    /// 返回 1 / f'(x0)
    auto one = SymbolicExpr::number(1);
    auto result = SymbolicExpr::divide(one, f_prime_at_x0);
    return result->simplify();
}

std::vector<std::shared_ptr<SymbolicExpr>> inverse_function(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& y)
{
    if (!f || !f->root || !y || !y->root) return {};

    /// 求解 f(var) = y，即 f(var) - y = 0
    auto eq = SymbolicExpr::add(f,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), y));
    return SymbolicExpr::solve(eq, var);
}

// ============================================================
// 内部辅助函数（曲率・旋转体表面积）
// ============================================================

/// 构造绝对值表达式 |expr|
static std::shared_ptr<SymbolicExpr> calculus_utils_make_abs(
    const std::shared_ptr<SymbolicExpr>& expr)
{
    return std::make_shared<SymbolicExpr>(
        std::make_shared<FunctionNode>(
            FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<SymbolicNode>>{expr->root}));
}

/// 尝试符号定积分，若结果仍含未求值积分节点则返回 nullptr
static std::shared_ptr<SymbolicExpr> calculus_utils_try_symbolic_definite(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    Integrator integrator;
    SymbolicExpr result = integrator.integrate_def(*integrand, var, *a, *b);

    // 检查结果是否仍含未求值的积分节点
    if (auto func = std::dynamic_pointer_cast<FunctionNode>(result.root)) {
        if (func->type == FunctionNode::FuncType::Calculus_Integral) {
            return nullptr;
        }
    }
    auto res = std::make_shared<SymbolicExpr>(result);
    auto simplified = res->simplify();
    return simplified ? simplified : res;
}

/// 数值定积分回退（复合 Simpson 法）
static std::shared_ptr<SymbolicExpr> calculus_utils_numerical_definite(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    double a_val = a->to_numeric();
    double b_val = b->to_numeric();

    if (std::isnan(a_val) || std::isnan(b_val) ||
        std::isinf(a_val) || std::isinf(b_val)) {
        return nullptr;
    }

    int n = 1000;
    double h = (b_val - a_val) / n;
    double sum = 0.0;

    for (int i = 0; i <= n; ++i) {
        double xi = a_val + i * h;
        auto xi_expr = SymbolicExpr::number(xi);
        auto fi = integrand->substitute(var, xi_expr);
        double fi_val = fi->to_numeric();
        if (std::isnan(fi_val) || std::isinf(fi_val)) {
            return nullptr;
        }

        double weight = 1.0;
        if (i == 0 || i == n) {
            weight = 1.0;
        } else if (i % 2 == 1) {
            weight = 4.0;
        } else {
            weight = 2.0;
        }
        sum += weight * fi_val;
    }
    sum *= h / 3.0;

    return SymbolicExpr::number(sum);
}

/// 符号积分优先，失败时回退到数值积分
static std::shared_ptr<SymbolicExpr> calculus_utils_integrate_with_fallback(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    auto symbolic_result = calculus_utils_try_symbolic_definite(integrand, var, a, b);
    if (symbolic_result) {
        return symbolic_result;
    }
    return calculus_utils_numerical_definite(integrand, var, a, b);
}

// ============================================================
// 曲率 (Requirement 16)
// ============================================================

std::shared_ptr<SymbolicExpr> curvature(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    // f' = df/dvar
    auto f_prime = f->differentiate(var);
    // f'' = d²f/dvar²
    auto f_double_prime = f_prime->differentiate(var);

    // |f''|
    auto abs_f_pp = calculus_utils_make_abs(f_double_prime);

    // 1 + f'²
    auto f_prime_sq = SymbolicExpr::power(f_prime, SymbolicExpr::number(2));
    auto one_plus_fp_sq = SymbolicExpr::add(SymbolicExpr::number(1), f_prime_sq);

    // (1 + f'²)^(3/2)
    auto three_half = SymbolicExpr::number(Rational(3, 2));
    auto denom = SymbolicExpr::power(one_plus_fp_sq, three_half);

    // κ = |f''| / (1 + f'²)^(3/2)
    auto result = SymbolicExpr::divide(abs_f_pp, denom);
    auto simplified = result->simplify();
    return simplified ? simplified : result;
}

std::shared_ptr<SymbolicExpr> curvature_parametric(
    const std::shared_ptr<SymbolicExpr>& x_t,
    const std::shared_ptr<SymbolicExpr>& y_t, const std::string& t)
{
    // x' = dx/dt, x'' = d²x/dt²
    auto x_prime = x_t->differentiate(t);
    auto x_double_prime = x_prime->differentiate(t);

    // y' = dy/dt, y'' = d²y/dt²
    auto y_prime = y_t->differentiate(t);
    auto y_double_prime = y_prime->differentiate(t);

    // x'y'' - y'x''
    auto cross_term = SymbolicExpr::add(
        SymbolicExpr::multiply(x_prime, y_double_prime),
        SymbolicExpr::multiply(
            SymbolicExpr::number(-1),
            SymbolicExpr::multiply(y_prime, x_double_prime)));

    // |x'y'' - y'x''|
    auto abs_cross = calculus_utils_make_abs(cross_term);

    // x'² + y'²
    auto x_prime_sq = SymbolicExpr::power(x_prime, SymbolicExpr::number(2));
    auto y_prime_sq = SymbolicExpr::power(y_prime, SymbolicExpr::number(2));
    auto sum_sq = SymbolicExpr::add(x_prime_sq, y_prime_sq);

    // (x'² + y'²)^(3/2)
    auto three_half = SymbolicExpr::number(Rational(3, 2));
    auto denom = SymbolicExpr::power(sum_sq, three_half);

    // κ = |x'y'' - y'x''| / (x'² + y'²)^(3/2)
    auto result = SymbolicExpr::divide(abs_cross, denom);
    auto simplified = result->simplify();
    return simplified ? simplified : result;
}

// ============================================================
// 拐点 (Requirement 34)
// ============================================================

std::vector<std::shared_ptr<SymbolicExpr>> inflection_points(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    if (!f || !f->root) return {};
    
    // f' = df/dvar
    auto f_prime = f->differentiate(var);
    if (!f_prime) return {};
    
    // f'' = d²f/dvar²
    auto f_double_prime = f_prime->differentiate(var);
    if (!f_double_prime) return {};
    
    // 寻找二阶导数为0的点
    return SymbolicExpr::solve(f_double_prime, var);
}

// ============================================================
// 旋转体表面积 (Requirement 17)
// ============================================================

std::shared_ptr<SymbolicExpr> surface_area_revolution_x(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b)
{
    // f' = df/dvar
    auto f_prime = f->differentiate(var);

    // √(1 + f'²)
    auto f_prime_sq = SymbolicExpr::power(f_prime, SymbolicExpr::number(2));
    auto one_plus_fp_sq = SymbolicExpr::add(SymbolicExpr::number(1), f_prime_sq);
    auto arc_factor = SymbolicExpr::sqrt(one_plus_fp_sq);

    // |f(x)| · √(1 + f'²)
    auto abs_f = calculus_utils_make_abs(f);
    auto integrand = SymbolicExpr::multiply(abs_f, arc_factor);

    // 2π
    auto two_pi = SymbolicExpr::number(2.0 * M_PI);

    // S = 2π ∫ₐᵇ |f(x)| · √(1 + f'²) dx
    auto integral = calculus_utils_integrate_with_fallback(integrand, var, a, b);
    if (!integral) {
        return nullptr;
    }
    auto result = SymbolicExpr::multiply(two_pi, integral);
    auto simplified = result->simplify();
    return simplified ? simplified : result;
}

std::shared_ptr<SymbolicExpr> surface_area_revolution_y(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b)
{
    // f' = df/dvar
    auto f_prime = f->differentiate(var);

    // √(1 + f'²)
    auto f_prime_sq = SymbolicExpr::power(f_prime, SymbolicExpr::number(2));
    auto one_plus_fp_sq = SymbolicExpr::add(SymbolicExpr::number(1), f_prime_sq);
    auto arc_factor = SymbolicExpr::sqrt(one_plus_fp_sq);

    // |var| · √(1 + f'²)
    auto var_expr = SymbolicExpr::variable(var);
    auto abs_var = calculus_utils_make_abs(var_expr);
    auto integrand = SymbolicExpr::multiply(abs_var, arc_factor);

    // 2π
    auto two_pi = SymbolicExpr::number(2.0 * M_PI);

    // S = 2π ∫ₐᵇ |x| · √(1 + f'²) dx
    auto integral = calculus_utils_integrate_with_fallback(integrand, var, a, b);
    if (!integral) {
        return nullptr;
    }
    auto result = SymbolicExpr::multiply(two_pi, integral);
    auto simplified = result->simplify();
    return simplified ? simplified : result;
}

} // namespace lamina
