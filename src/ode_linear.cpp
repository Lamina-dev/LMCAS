/**
 * @file symbolic_ode_engine.cpp
 * @brief 统一 ODE 求解引擎实现：类型检测与分类。
 */
#include "../include/symbolic_ode_engine.hpp"
#include "symbolic_ast.hpp"
#include "../include/symbolic.hpp"
#include "../include/poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "internal/ode_support.hpp"
#include <cmath>
#include <memory>
#include <string>

namespace lamina {

static ODESolution solve_higher_order_ode_impl(
    const std::vector<double>&,
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&, const std::string&);
static ODESolution solve_euler_ode_impl(
    const std::vector<double>&,
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&, const std::string&);

/**
 * @internal
 * @brief 表示特征多项式的一个根及其重数。
 */
struct CharRoot {
    double real_part;   ///< 实部
    double imag_part;   ///< 虚部（为零表示实根）
    int multiplicity;   ///< 重数
    bool is_complex;    ///< 是否为复根
};

/**
 * @internal
 * @brief 使用数值方法求解特征多项式的所有根。
 *
 * 对于度数 ≤ 4 的多项式使用解析公式，
 * 对于度数 5-6 使用 Durand-Kerner 迭代法。
 */
[[maybe_unused]] static std::vector<CharRoot> find_characteristic_roots(
    const std::vector<double>& coeffs)
{
    int n = static_cast<int>(coeffs.size()) - 1;
    if (n <= 0) return {};

    std::vector<CharRoot> roots;

    /// 归一化系数（使最高次项系数为 1）
    double leading = coeffs[0];
    if (std::abs(leading) < 1e-15) return {};

    std::vector<double> norm_coeffs(coeffs.size());
    for (size_t i = 0; i < coeffs.size(); ++i) {
        norm_coeffs[i] = coeffs[i] / leading;
    }

    /// 对于低阶多项式，使用解析公式
    if (n == 1) {
        /// r + norm_coeffs[1] = 0
        double r = -norm_coeffs[1];
        roots.push_back({r, 0.0, 1, false});
        return roots;
    }

    if (n == 2) {
        double b = norm_coeffs[1];
        double c = norm_coeffs[2];
        double D = b * b - 4.0 * c;
        int eq;
        lmmc_double_nearly_equal_tol(D, 0.0, 1e-10, 1e-10, &eq);
        if (eq) {
            roots.push_back({-b / 2.0, 0.0, 2, false});
        } else if (D > 0) {
            double r1 = (-b + std::sqrt(D)) / 2.0;
            double r2 = (-b - std::sqrt(D)) / 2.0;
            roots.push_back({r1, 0.0, 1, false});
            roots.push_back({r2, 0.0, 1, false});
        } else {
            double re = -b / 2.0;
            double im = std::sqrt(-D) / 2.0;
            roots.push_back({re, im, 1, true});
        }
        return roots;
    }

    /// 对于 n >= 3，使用 Durand-Kerner 方法求所有根
    /// 初始化：在单位圆上均匀分布初始猜测
    struct Complex {
        double re, im;
        Complex(double r = 0, double i = 0) : re(r), im(i) {}
        Complex operator*(const Complex& o) const {
            return {re * o.re - im * o.im, re * o.im + im * o.re};
        }
        Complex operator+(const Complex& o) const {
            return {re + o.re, im + o.im};
        }
        Complex operator-(const Complex& o) const {
            return {re - o.re, im - o.im};
        }
        Complex operator/(const Complex& o) const {
            double d = o.re * o.re + o.im * o.im;
            if (d < 1e-30) return {0, 0};
            return {(re * o.re + im * o.im) / d,
                    (im * o.re - re * o.im) / d};
        }
        double mag() const { return std::sqrt(re * re + im * im); }
    };

    std::vector<Complex> z(n);
    /// 初始猜测：使用不同半径的点避免对称性问题
    for (int i = 0; i < n; ++i) {
        double angle = 2.0 * 3.14159265358979323846 * i / n + 0.1;
        double radius = 1.0 + 0.3 * i;
        z[i] = {radius * std::cos(angle), radius * std::sin(angle)};
    }

    /// 求值多项式 p(z)
    auto eval_poly = [&](const Complex& val) -> Complex {
        Complex result = {1.0, 0.0};
        for (int i = 1; i <= n; ++i) {
            result = result * val + Complex{norm_coeffs[i], 0.0};
        }
        return result;
    };

    /// Durand-Kerner 迭代
    for (int iter = 0; iter < 1000; ++iter) {
        double max_change = 0.0;
        for (int i = 0; i < n; ++i) {
            Complex num = eval_poly(z[i]);
            Complex denom = {1.0, 0.0};
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    denom = denom * (z[i] - z[j]);
                }
            }
            Complex delta = num / denom;
            z[i] = z[i] - delta;
            double change = delta.mag();
            if (change > max_change) max_change = change;
        }
        if (max_change < 1e-12) break;
    }

    /// 将数值根分类为实根和复共轭对，并检测重根
    std::vector<bool> used(n, false);
    for (int i = 0; i < n; ++i) {
        if (used[i]) continue;

        /// 检查是否为实根（虚部接近零）
        int is_real;
        lmmc_double_nearly_equal_tol(z[i].im, 0.0, 1e-8, 1e-8, &is_real);

        if (is_real) {
            double r = z[i].re;
            int mult = 1;
            used[i] = true;
            /// 检查重根
            for (int j = i + 1; j < n; ++j) {
                if (used[j]) continue;
                int j_real;
                lmmc_double_nearly_equal_tol(z[j].im, 0.0, 1e-8, 1e-8, &j_real);
                if (!j_real) continue;
                int same;
                lmmc_double_nearly_equal_tol(z[j].re, r, 1e-8, 1e-8, &same);
                if (same) {
                    mult++;
                    used[j] = true;
                }
            }
            roots.push_back({r, 0.0, mult, false});
        } else {
            /// 复根：找共轭对
            double re = z[i].re;
            double im = std::abs(z[i].im);
            used[i] = true;
            int mult = 1;

            /// 找到共轭根并标记
            for (int j = i + 1; j < n; ++j) {
                if (used[j]) continue;
                int same_re, conj_im;
                lmmc_double_nearly_equal_tol(z[j].re, re, 1e-8, 1e-8, &same_re);
                lmmc_double_nearly_equal_tol(z[j].im, -z[i].im, 1e-8, 1e-8, &conj_im);
                if (same_re && conj_im) {
                    used[j] = true;
                    break;
                }
            }
            /// 检查重复的复共轭对
            for (int j = i + 1; j < n; ++j) {
                if (used[j]) continue;
                int same_re, same_im;
                lmmc_double_nearly_equal_tol(z[j].re, re, 1e-8, 1e-8, &same_re);
                lmmc_double_nearly_equal_tol(std::abs(z[j].im), im, 1e-8, 1e-8, &same_im);
                if (same_re && same_im) {
                    mult++;
                    used[j] = true;
                    /// 也标记其共轭
                    for (int k = j + 1; k < n; ++k) {
                        if (used[k]) continue;
                        int k_re, k_im;
                        lmmc_double_nearly_equal_tol(z[k].re, re, 1e-8, 1e-8, &k_re);
                        lmmc_double_nearly_equal_tol(z[k].im, -z[j].im, 1e-8, 1e-8, &k_im);
                        if (k_re && k_im) {
                            used[k] = true;
                            break;
                        }
                    }
                }
            }
            roots.push_back({re, im, mult, true});
        }
    }

    return roots;
}

/**
 * @internal
 * @brief 将 double 值转为"干净"的数值表达式。
 *
 * 若值接近整数或简单分数，使用精确表示。
 */
static std::shared_ptr<SymbolicExpr> clean_number(double val) {
    /// 检查是否接近整数
    double rounded = std::round(val);
    int eq;
    lmmc_double_nearly_equal_tol(val, rounded, 1e-10, 1e-10, &eq);
    if (eq) {
        return SymbolicExpr::number(static_cast<int>(rounded));
    }

    /// 检查是否接近简单分数 p/q (q <= 12)
    for (int q = 2; q <= 12; ++q) {
        double p = val * q;
        double p_rounded = std::round(p);
        lmmc_double_nearly_equal_tol(p, p_rounded, 1e-10, 1e-10, &eq);
        if (eq) {
            return SymbolicExpr::divide(
                SymbolicExpr::number(static_cast<int>(p_rounded)),
                SymbolicExpr::number(q));
        }
    }

    return SymbolicExpr::number(val);
}

static bool has_nonzero_forcing(const std::shared_ptr<SymbolicExpr>& forcing) {
    return forcing && lamina::detail::node(forcing) && !forcing->is_zero();
}

/**
 * @internal
 * @brief 根据特征根构造齐次通解。
 *
 * - 实根 r（重数 m）：C_k * x^k * e^(rx)，k = 0, ..., m-1
 * - 复根 α±βi（重数 m）：x^k * e^(αx) * (C_a*cos(βx) + C_b*sin(βx))
 */
[[maybe_unused]] static std::shared_ptr<SymbolicExpr> build_homogeneous_solution(
    const std::vector<CharRoot>& roots,
    const std::string& x,
    std::vector<std::string>& constants)
{
    auto x_var = SymbolicExpr::variable(x);
    std::shared_ptr<SymbolicExpr> solution = nullptr;
    int const_idx = 1;

    for (const auto& root : roots) {
        for (int k = 0; k < root.multiplicity; ++k) {
            if (!root.is_complex) {
                /// 实根: C_i * x^k * e^(r*x)
                std::string c_name = "C" + std::to_string(const_idx++);
                constants.push_back(c_name);
                auto C = SymbolicExpr::variable(c_name);

                std::shared_ptr<SymbolicExpr> term = C;

                /// x^k 因子
                if (k > 0) {
                    auto x_pow = SymbolicExpr::power(x_var, SymbolicExpr::number(k));
                    term = SymbolicExpr::multiply(term, x_pow);
                }

                /// e^(r*x) 因子
                int r_zero;
                lmmc_double_nearly_equal_tol(root.real_part, 0.0, 1e-10, 1e-10, &r_zero);
                if (!r_zero) {
                    auto r_expr = clean_number(root.real_part);
                    auto exp_arg = SymbolicExpr::multiply(r_expr, x_var);
                    auto exp_term = SymbolicExpr::exp(exp_arg);
                    term = SymbolicExpr::multiply(term, exp_term);
                }

                solution = solution ? SymbolicExpr::add(solution, term) : term;
            } else {
                /// 复根 α±βi: 产生两个基本解
                /// C_a * x^k * e^(αx) * cos(βx)
                /// C_b * x^k * e^(αx) * sin(βx)
                std::string ca_name = "C" + std::to_string(const_idx++);
                std::string cb_name = "C" + std::to_string(const_idx++);
                constants.push_back(ca_name);
                constants.push_back(cb_name);
                auto Ca = SymbolicExpr::variable(ca_name);
                auto Cb = SymbolicExpr::variable(cb_name);

                auto beta_expr = clean_number(root.imag_part);
                auto beta_x = SymbolicExpr::multiply(beta_expr, x_var);
                auto cos_term = SymbolicExpr::cos(beta_x);
                auto sin_term = SymbolicExpr::sin(beta_x);

                /// e^(αx) 因子
                std::shared_ptr<SymbolicExpr> exp_factor = nullptr;
                int alpha_zero;
                lmmc_double_nearly_equal_tol(root.real_part, 0.0, 1e-10, 1e-10, &alpha_zero);
                if (!alpha_zero) {
                    auto alpha_expr = clean_number(root.real_part);
                    auto exp_arg = SymbolicExpr::multiply(alpha_expr, x_var);
                    exp_factor = SymbolicExpr::exp(exp_arg);
                }

                /// x^k 因子
                std::shared_ptr<SymbolicExpr> x_pow_factor = nullptr;
                if (k > 0) {
                    x_pow_factor = SymbolicExpr::power(x_var, SymbolicExpr::number(k));
                }

                /// 构造 cos 项: Ca * x^k * e^(αx) * cos(βx)
                auto term_cos = Ca;
                if (x_pow_factor) term_cos = SymbolicExpr::multiply(term_cos, x_pow_factor);
                if (exp_factor) term_cos = SymbolicExpr::multiply(term_cos, exp_factor);
                term_cos = SymbolicExpr::multiply(term_cos, cos_term);

                /// 构造 sin 项: Cb * x^k * e^(αx) * sin(βx)
                auto term_sin = Cb;
                if (x_pow_factor) term_sin = SymbolicExpr::multiply(term_sin, x_pow_factor);
                if (exp_factor) term_sin = SymbolicExpr::multiply(term_sin, exp_factor);
                term_sin = SymbolicExpr::multiply(term_sin, sin_term);

                solution = solution ? SymbolicExpr::add(solution, term_cos) : term_cos;
                solution = SymbolicExpr::add(solution, term_sin);
            }
        }
    }

    return solution;
}

static std::shared_ptr<SymbolicExpr> build_euler_solution(
    const std::vector<CharRoot>& roots,
    const std::string& x,
    std::vector<std::string>& constants)
{
    auto x_var = SymbolicExpr::variable(x);
    auto ln_x = SymbolicExpr::ln(x_var);
    std::shared_ptr<SymbolicExpr> solution = nullptr;
    int const_idx = 1;

    for (const auto& root : roots) {
        for (int k = 0; k < root.multiplicity; ++k) {
            std::shared_ptr<SymbolicExpr> log_factor = nullptr;
            if (k > 0) {
                log_factor = SymbolicExpr::power(ln_x, SymbolicExpr::number(k));
            }

            if (!root.is_complex) {
                std::string c_name = "C" + std::to_string(const_idx++);
                constants.push_back(c_name);
                auto term = SymbolicExpr::variable(c_name);

                auto exponent = clean_number(root.real_part);
                term = SymbolicExpr::multiply(
                    term,
                    SymbolicExpr::power(x_var, exponent));
                if (log_factor) {
                    term = SymbolicExpr::multiply(term, log_factor);
                }
                solution = solution ? SymbolicExpr::add(solution, term) : term;
            } else {
                std::string ca_name = "C" + std::to_string(const_idx++);
                std::string cb_name = "C" + std::to_string(const_idx++);
                constants.push_back(ca_name);
                constants.push_back(cb_name);

                auto alpha_expr = clean_number(root.real_part);
                auto beta_expr = clean_number(root.imag_part);
                auto beta_ln_x = SymbolicExpr::multiply(beta_expr, ln_x);
                auto trig = SymbolicExpr::add(
                    SymbolicExpr::multiply(SymbolicExpr::variable(ca_name),
                                           SymbolicExpr::cos(beta_ln_x)),
                    SymbolicExpr::multiply(SymbolicExpr::variable(cb_name),
                                           SymbolicExpr::sin(beta_ln_x)));
                auto term = SymbolicExpr::multiply(
                    SymbolicExpr::power(x_var, alpha_expr),
                    trig);
                if (log_factor) {
                    term = SymbolicExpr::multiply(term, log_factor);
                }
                solution = solution ? SymbolicExpr::add(solution, term) : term;
            }
        }
    }

    return solution;
}

ODESolutionResult solve_higher_order_ode_checked(
    const std::vector<double>& coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string& y,
    ComputationContext& context)
{
    const std::string operation = "solve_higher_order_ode";
    auto variables = validate_ode_variables(x, y, context, operation);
    if (!variables) return ODESolutionResult::failure(variables.error());

    auto coeff_check = validate_numeric_ode_coefficients(coeffs, 2, 7, operation);
    if (!coeff_check) return ODESolutionResult::failure(coeff_check.error());

    auto budget = context.consume_steps(coeffs.size() * 20 + 20, operation);
    if (!budget) return ODESolutionResult::failure(budget.error());

    const bool nonzero_forcing = has_nonzero_forcing(forcing);
    if (nonzero_forcing &&
        (!forcing->is_number() || coeffs.back() == 0.0)) {
        return ODESolutionResult::failure(
            CasErrc::Inconclusive,
            "non-homogeneous constant-coefficient ODE currently requires a constant forcing and nonzero y coefficient",
            operation);
    }

    try {
        auto solution = solve_higher_order_ode_impl(
            coeffs, nullptr, x, y);
        if (nonzero_forcing && solution.general_solution) {
            auto particular = SymbolicExpr::divide(
                forcing, SymbolicExpr::number(coeffs.back()))->simplify();
            solution.general_solution = SymbolicExpr::add(
                solution.general_solution, particular)->simplify();
        }
        return wrap_ode_solution(
            std::move(solution),
            ODEType::HigherOrder_ConstCoeff,
            operation);
    } catch (const std::bad_alloc&) {
        return ODESolutionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while solving higher-order ODE",
            operation);
    } catch (const std::exception& ex) {
        return ODESolutionResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            operation);
    }
}

ODESolutionResult solve_higher_order_ode_checked(
    const std::vector<double>& coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string& y)
{
    ComputationContext context;
    return solve_higher_order_ode_checked(coeffs, forcing, x, y, context);
}

static ODESolution solve_higher_order_ode_impl(
    const std::vector<double>& coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string&)
{
    ODESolution result;
    result.method_used = ODEType::HigherOrder_ConstCoeff;

    if (coeffs.size() < 2 || has_nonzero_forcing(forcing)) {
        result.general_solution = nullptr;
        return result;
    }

    auto roots = find_characteristic_roots(coeffs);
    if (roots.empty()) {
        result.general_solution = nullptr;
        return result;
    }

    result.general_solution = build_homogeneous_solution(roots, x, result.constants);
    if (result.general_solution) {
        result.general_solution = result.general_solution->simplify();
    }
    return result;
}

ODESolutionResult solve_euler_ode_checked(
    const std::vector<double>& euler_coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string& y,
    ComputationContext& context)
{
    const std::string operation = "solve_euler_ode";
    auto variables = validate_ode_variables(x, y, context, operation);
    if (!variables) return ODESolutionResult::failure(variables.error());

    auto coeff_check = validate_numeric_ode_coefficients(euler_coeffs, 3, 4, operation);
    if (!coeff_check) return ODESolutionResult::failure(coeff_check.error());

    auto budget = context.consume_steps(euler_coeffs.size() * 20 + 20, operation);
    if (!budget) return ODESolutionResult::failure(budget.error());

    const bool nonzero_forcing = has_nonzero_forcing(forcing);
    if (nonzero_forcing &&
        (!forcing->is_number() || euler_coeffs.back() == 0.0)) {
        return ODESolutionResult::failure(
            CasErrc::Inconclusive,
            "non-homogeneous Euler ODE currently requires a constant forcing and nonzero y coefficient",
            operation);
    }

    try {
        auto solution = solve_euler_ode_impl(
            euler_coeffs, nullptr, x, y);
        if (nonzero_forcing && solution.general_solution) {
            auto particular = SymbolicExpr::divide(
                forcing,
                SymbolicExpr::number(euler_coeffs.back()))->simplify();
            solution.general_solution = SymbolicExpr::add(
                solution.general_solution, particular)->simplify();
        }
        return wrap_ode_solution(
            std::move(solution),
            ODEType::Euler,
            operation);
    } catch (const std::bad_alloc&) {
        return ODESolutionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while solving Euler ODE",
            operation);
    } catch (const std::exception& ex) {
        return ODESolutionResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            operation);
    }
}

ODESolutionResult solve_euler_ode_checked(
    const std::vector<double>& euler_coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string& y)
{
    ComputationContext context;
    return solve_euler_ode_checked(euler_coeffs, forcing, x, y, context);
}

static ODESolution solve_euler_ode_impl(
    const std::vector<double>& euler_coeffs,
    const std::shared_ptr<SymbolicExpr>& forcing,
    const std::string& x,
    const std::string&)
{
    ODESolution result;
    result.method_used = ODEType::Euler;

    if ((euler_coeffs.size() != 3 && euler_coeffs.size() != 4) ||
        has_nonzero_forcing(forcing)) {
        result.general_solution = nullptr;
        return result;
    }

    std::vector<double> characteristic;
    if (euler_coeffs.size() == 3) {
        const double a = euler_coeffs[0];
        const double b = euler_coeffs[1];
        const double c = euler_coeffs[2];
        characteristic = {a, b - a, c};
    } else {
        const double a = euler_coeffs[0];
        const double b = euler_coeffs[1];
        const double c = euler_coeffs[2];
        const double d = euler_coeffs[3];
        characteristic = {a, b - 3.0 * a, 2.0 * a - b + c, d};
    }

    auto roots = find_characteristic_roots(characteristic);
    if (roots.empty()) {
        result.general_solution = nullptr;
        return result;
    }

    result.general_solution = build_euler_solution(roots, x, result.constants);
    if (result.general_solution) {
        result.general_solution = result.general_solution->simplify();
    }
    return result;
}

/**
 * @internal
 * @brief 检测非齐次项的类型，用于待定系数法。
 *
 * 支持的类型：
 * - 多项式: x^n
 * - 指数: e^(ax)
 * - 三角: cos(bx), sin(bx)
 */

} // namespace lamina
