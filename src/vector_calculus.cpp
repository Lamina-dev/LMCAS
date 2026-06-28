/**
 * @file vector_calculus.cpp
 * @brief 向量微积分模块实现：梯度、散度、旋度、拉普拉斯算子、方向导数、曲线积分、曲面积分、极值分析。
 */

#include "vector_calculus.hpp"
#include "integration.hpp"
#include "solver.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <cmath>
#include <map>
#include <set>

namespace lamina {

// ============================================================
/// 梯度 ∇f (Requirement 9)
// ============================================================

VectorField gradient(const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars)
{
    if (!f) {
        throw std::invalid_argument("gradient: f must not be null");
    }
    VectorField result;
    result.reserve(vars.size());
    for (const auto& var : vars) {
        auto partial = f->differentiate(var);
        if (partial) {
            partial = partial->simplify();
        }
        result.push_back(partial);
    }
    return result;
}

// ============================================================
/// 散度 ∇·F (Requirement 45)
// ============================================================

std::shared_ptr<SymbolicExpr> divergence(const VectorField& F,
    const std::vector<std::string>& vars)
{
    if (F.size() != vars.size()) {
        throw std::invalid_argument(
            "divergence: F and vars must have the same dimension");
    }
    if (F.empty()) {
        return SymbolicExpr::number(0);
    }

    std::shared_ptr<SymbolicExpr> sum = nullptr;
    for (size_t i = 0; i < F.size(); ++i) {
        if (!F[i]) continue;
        auto partial = F[i]->differentiate(vars[i]);
        if (!partial) continue;
        partial = partial->simplify();
        if (!sum) {
            sum = partial;
        } else {
            sum = SymbolicExpr::add(sum, partial);
            sum = sum->simplify();
        }
    }

    if (!sum) {
        return SymbolicExpr::number(0);
    }
    return sum->simplify();
}

// ============================================================
/// 旋度 ∇×F (Requirement 46)
// ============================================================

VectorField curl(const VectorField& F,
    const std::vector<std::string>& vars)
{
    if (F.size() != vars.size()) {
        throw std::invalid_argument(
            "curl: F and vars must have the same dimension");
    }

    /// 二维标量旋度: ∂F₂/∂x₁ - ∂F₁/∂x₂
    if (F.size() == 2 && vars.size() == 2) {
        auto dF2_dx1 = F[1]->differentiate(vars[0]);
        auto dF1_dx2 = F[0]->differentiate(vars[1]);

        auto neg_dF1_dx2 = SymbolicExpr::multiply(
            SymbolicExpr::number(-1), dF1_dx2);
        auto scalar_curl = SymbolicExpr::add(dF2_dx1, neg_dF1_dx2);
        scalar_curl = scalar_curl->simplify();

        return VectorField{scalar_curl};
    }

    /// 三维旋度
    if (F.size() != 3 || vars.size() != 3) {
        throw std::invalid_argument(
            "curl: requires 2D or 3D vector field");
    }

    /// curl_x = ∂F₃/∂x₂ - ∂F₂/∂x₃
    auto dF3_dx2 = F[2]->differentiate(vars[1]);
    auto dF2_dx3 = F[1]->differentiate(vars[2]);
    auto curl_x = SymbolicExpr::add(
        dF3_dx2,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dF2_dx3));
    curl_x = curl_x->simplify();

    /// curl_y = ∂F₁/∂x₃ - ∂F₃/∂x₁
    auto dF1_dx3 = F[0]->differentiate(vars[2]);
    auto dF3_dx1 = F[2]->differentiate(vars[0]);
    auto curl_y = SymbolicExpr::add(
        dF1_dx3,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dF3_dx1));
    curl_y = curl_y->simplify();

    /// curl_z = ∂F₂/∂x₁ - ∂F₁/∂x₂
    auto dF2_dx1 = F[1]->differentiate(vars[0]);
    auto dF1_dx2 = F[0]->differentiate(vars[1]);
    auto curl_z = SymbolicExpr::add(
        dF2_dx1,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dF1_dx2));
    curl_z = curl_z->simplify();

    return VectorField{curl_x, curl_y, curl_z};
}

// ============================================================
/// 拉普拉斯算子 ∇²f (Requirement 47)
// ============================================================

std::shared_ptr<SymbolicExpr> laplacian(const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars)
{
    if (!f) {
        throw std::invalid_argument("laplacian: f must not be null");
    }
    if (vars.empty()) {
        return SymbolicExpr::number(0);
    }

    std::shared_ptr<SymbolicExpr> sum = nullptr;
    for (const auto& var : vars) {
        auto first = f->differentiate(var);
        if (!first) continue;
        auto second = first->differentiate(var);
        if (!second) continue;
        second = second->simplify();
        if (!sum) {
            sum = second;
        } else {
            sum = SymbolicExpr::add(sum, second);
            sum = sum->simplify();
        }
    }

    if (!sum) {
        return SymbolicExpr::number(0);
    }
    return sum->simplify();
}

// ============================================================
/// 方向导数 (Requirement 8, 87, 88)
// ============================================================

/**
 * @internal
 * @brief 计算向量的模长平方（符号表达式），并化简。
 */
static std::shared_ptr<SymbolicExpr> vector_calculus_magnitude_squared(
    const VectorField& v)
{
    std::shared_ptr<SymbolicExpr> sum = nullptr;
    for (const auto& comp : v) {
        if (!comp) continue;
        auto sq = SymbolicExpr::power(comp, SymbolicExpr::number(2));
        sq = sq->simplify();
        if (!sum) {
            sum = sq;
        } else {
            sum = SymbolicExpr::add(sum, sq);
            sum = sum->simplify();
        }
    }
    if (!sum) return SymbolicExpr::number(0);
    return sum->simplify();
}

/**
 * @internal
 * @brief 计算单阶方向导数 D_u f = ∑(∂f/∂xᵢ · uᵢ)。
 *
 * 此处 u 已经是单位向量分量。
 */
static std::shared_ptr<SymbolicExpr> vector_calculus_single_dir_deriv(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    const VectorField& unit_dir)
{
    std::shared_ptr<SymbolicExpr> sum = nullptr;
    for (size_t i = 0; i < vars.size(); ++i) {
        auto partial = f->differentiate(vars[i]);
        if (!partial) continue;
        partial = partial->simplify();
        auto term = SymbolicExpr::multiply(partial, unit_dir[i]);
        term = term->simplify();
        if (!sum) {
            sum = term;
        } else {
            sum = SymbolicExpr::add(sum, term);
            sum = sum->simplify();
        }
    }
    if (!sum) return SymbolicExpr::number(0);
    return sum->simplify();
}

std::shared_ptr<SymbolicExpr> directional_derivative(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars,
    const VectorField& direction, int order)
{
    if (!f) {
        throw std::invalid_argument("directional_derivative: f must not be null");
    }
    if (direction.size() != vars.size()) {
        throw std::invalid_argument(
            "directional_derivative: direction and vars must have the same dimension");
    }
    if (order < 1) {
        throw std::invalid_argument(
            "directional_derivative: order must be >= 1");
    }

    /// 计算方向向量的模长平方
    auto mag_sq = vector_calculus_magnitude_squared(direction);
    if (!mag_sq || mag_sq->is_zero()) {
        /// 零向量，返回 nullptr 表示错误
        return nullptr;
    }

    /// 构造单位向量分量: uᵢ = dirᵢ / |dir|
    auto magnitude = SymbolicExpr::sqrt(mag_sq);
    magnitude = magnitude->simplify();
    VectorField unit_dir;
    unit_dir.reserve(direction.size());
    for (const auto& comp : direction) {
        if (!comp || comp->is_zero()) {
            unit_dir.push_back(SymbolicExpr::number(0));
        } else {
            auto u_comp = SymbolicExpr::multiply(
                comp,
                SymbolicExpr::power(magnitude, SymbolicExpr::number(-1)));
            u_comp = u_comp->simplify();
            unit_dir.push_back(u_comp);
        }
    }

    /// 递归应用方向导数 order 次
    auto result = f;
    for (int k = 0; k < order; ++k) {
        result = vector_calculus_single_dir_deriv(result, vars, unit_dir);
        if (!result) return SymbolicExpr::number(0);
    }

    return result;
}

// ============================================================
/// 雅可比矩阵 (Requirement 10)
// ============================================================

std::shared_ptr<SymbolicExpr> jacobian(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars)
{
    size_t m = functions.size();
    size_t n = vars.size();

    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid;
    grid.reserve(m);

    for (size_t i = 0; i < m; ++i) {
        std::vector<std::shared_ptr<SymbolicExpr>> row;
        row.reserve(n);
        for (size_t j = 0; j < n; ++j) {
            auto partial = functions[i]->differentiate(vars[j]);
            partial = partial->simplify();
            row.push_back(partial);
        }
        grid.push_back(std::move(row));
    }

    return SymbolicExpr::matrix(grid);
}

// ============================================================
/// 海森矩阵 (Requirement 11)
// ============================================================

std::shared_ptr<SymbolicExpr> hessian(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars)
{
    size_t n = vars.size();

    /// 先计算一阶偏导数
    std::vector<std::shared_ptr<SymbolicExpr>> first_partials;
    first_partials.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        first_partials.push_back(f->differentiate(vars[i]));
    }

    /// 构建 n×n 对称矩阵
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            auto second_partial = first_partials[i]->differentiate(vars[j]);
            second_partial = second_partial->simplify();
            grid[i][j] = second_partial;
            if (i != j) {
                grid[j][i] = second_partial;
            }
        }
    }

    return SymbolicExpr::matrix(grid);
}

// ============================================================
/// 曲线积分与曲面积分辅助函数
// ============================================================

/**
 * @internal
 * @brief 获取参数化曲线对应的坐标变量名列表。
 *
 * 根据维度返回 {"x", "y"} 或 {"x", "y", "z"}。
 */
static std::vector<std::string> vector_calculus_coord_vars(size_t dim)
{
    if (dim == 2) return {"x", "y"};
    return {"x", "y", "z"};
}

/**
 * @internal
 * @brief 尝试符号定积分，若结果仍含未求值积分节点则返回 nullptr。
 */
static std::shared_ptr<SymbolicExpr> vector_calculus_try_definite(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    Integrator integrator;
    SymbolicExpr result = integrator.integrate_def(*integrand, var, *a, *b);

    if (auto func = std::dynamic_pointer_cast<FunctionNode>(result.root)) {
        if (func->type == FunctionNode::FuncType::Calculus_Integral) {
            return nullptr;
        }
    }
    auto res = std::make_shared<SymbolicExpr>(result);
    auto simplified = res->simplify();
    return simplified ? simplified : res;
}

/**
 * @internal
 * @brief 数值定积分回退（复合 Simpson 法）。
 */
static std::shared_ptr<SymbolicExpr> vector_calculus_numerical_definite(
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

/**
 * @internal
 * @brief 符号积分优先，失败时回退到数值积分。
 */
static std::shared_ptr<SymbolicExpr> vector_calculus_integrate_with_fallback(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    auto symbolic_result = vector_calculus_try_definite(integrand, var, a, b);
    if (symbolic_result) {
        return symbolic_result;
    }
    return vector_calculus_numerical_definite(integrand, var, a, b);
}

// ============================================================
/// 第一类曲线积分 (Requirement 48.1, 48.3, 48.4)
// ============================================================

std::shared_ptr<SymbolicExpr> curve_integral_scalar(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    if (!f) {
        throw std::invalid_argument("curve_integral_scalar: f must not be null");
    }
    if (parametrization.empty()) {
        throw std::invalid_argument("curve_integral_scalar: parametrization must not be empty");
    }
    if (!a || !b) {
        throw std::invalid_argument("curve_integral_scalar: bounds must not be null");
    }

    size_t dim = parametrization.size();
    auto coord_vars = vector_calculus_coord_vars(dim);

    /// 将 f 中的坐标变量替换为参数化表达式: f(r(t))
    auto f_composed = f;
    for (size_t i = 0; i < dim; ++i) {
        f_composed = f_composed->substitute(coord_vars[i], parametrization[i]);
    }
    f_composed = f_composed->simplify();

    /// 计算 |r'(t)| = √(∑(r_i'(t))²)
    std::shared_ptr<SymbolicExpr> speed_sq = nullptr;
    for (size_t i = 0; i < dim; ++i) {
        auto deriv = parametrization[i]->differentiate(t);
        if (!deriv) continue;
        deriv = deriv->simplify();
        auto sq = SymbolicExpr::power(deriv, SymbolicExpr::number(2));
        sq = sq->simplify();
        if (!speed_sq) {
            speed_sq = sq;
        } else {
            speed_sq = SymbolicExpr::add(speed_sq, sq);
            speed_sq = speed_sq->simplify();
        }
    }

    if (!speed_sq) {
        return SymbolicExpr::number(0);
    }

    auto speed = SymbolicExpr::sqrt(speed_sq);
    speed = speed->simplify();

    /// 被积函数: f(r(t)) · |r'(t)|
    auto integrand = SymbolicExpr::multiply(f_composed, speed);
    integrand = integrand->simplify();

    /// 计算定积分 ∫ₐᵇ f(r(t))·|r'(t)| dt
    return vector_calculus_integrate_with_fallback(integrand, t, a, b);
}

// ============================================================
/// 第二类曲线积分 (Requirement 48.2)
// ============================================================

std::shared_ptr<SymbolicExpr> curve_integral_vector(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    if (F.empty()) {
        throw std::invalid_argument("curve_integral_vector: F must not be empty");
    }
    if (parametrization.empty()) {
        throw std::invalid_argument("curve_integral_vector: parametrization must not be empty");
    }
    if (F.size() != parametrization.size()) {
        throw std::invalid_argument(
            "curve_integral_vector: F and parametrization must have the same dimension");
    }
    if (!a || !b) {
        throw std::invalid_argument("curve_integral_vector: bounds must not be null");
    }

    size_t dim = parametrization.size();
    auto coord_vars = vector_calculus_coord_vars(dim);

    /// 计算 F(r(t))·r'(t) = ∑ Fᵢ(r(t)) · rᵢ'(t)
    std::shared_ptr<SymbolicExpr> dot_product = nullptr;
    for (size_t i = 0; i < dim; ++i) {
        if (!F[i]) continue;

        /// 将 Fᵢ 中的坐标变量替换为参数化表达式
        auto Fi_composed = F[i];
        for (size_t j = 0; j < dim; ++j) {
            Fi_composed = Fi_composed->substitute(coord_vars[j], parametrization[j]);
        }
        Fi_composed = Fi_composed->simplify();

        /// rᵢ'(t)
        auto ri_prime = parametrization[i]->differentiate(t);
        if (!ri_prime) continue;
        ri_prime = ri_prime->simplify();

        /// Fᵢ(r(t)) · rᵢ'(t)
        auto term = SymbolicExpr::multiply(Fi_composed, ri_prime);
        term = term->simplify();

        if (!dot_product) {
            dot_product = term;
        } else {
            dot_product = SymbolicExpr::add(dot_product, term);
            dot_product = dot_product->simplify();
        }
    }

    if (!dot_product) {
        return SymbolicExpr::number(0);
    }

    /// 计算定积分 ∫ₐᵇ F(r(t))·r'(t) dt
    return vector_calculus_integrate_with_fallback(dot_product, t, a, b);
}

// ============================================================
/// 曲面积分辅助：计算 r_u × r_v (Requirement 49.3)
// ============================================================

/**
 * @internal
 * @brief 计算参数化曲面的法向量 r_u × r_v。
 *
 * 参数化曲面 r(u,v) = [x(u,v), y(u,v), z(u,v)]，
 * 返回叉积 r_u × r_v 的三个分量。
 */
static VectorField vector_calculus_cross_product_partials(
    const VectorField& parametrization,
    const std::string& u, const std::string& v)
{
    /// r_u = [∂x/∂u, ∂y/∂u, ∂z/∂u]
    VectorField r_u;
    r_u.reserve(3);
    for (size_t i = 0; i < 3; ++i) {
        auto deriv = parametrization[i]->differentiate(u);
        r_u.push_back(deriv ? deriv->simplify() : SymbolicExpr::number(0));
    }

    /// r_v = [∂x/∂v, ∂y/∂v, ∂z/∂v]
    VectorField r_v;
    r_v.reserve(3);
    for (size_t i = 0; i < 3; ++i) {
        auto deriv = parametrization[i]->differentiate(v);
        r_v.push_back(deriv ? deriv->simplify() : SymbolicExpr::number(0));
    }

    /// 叉积: r_u × r_v = (r_u[1]*r_v[2] - r_u[2]*r_v[1],
    ///                     r_u[2]*r_v[0] - r_u[0]*r_v[2],
    ///                     r_u[0]*r_v[1] - r_u[1]*r_v[0])
    auto cross_x = SymbolicExpr::add(
        SymbolicExpr::multiply(r_u[1], r_v[2]),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(r_u[2], r_v[1])));
    cross_x = cross_x->simplify();

    auto cross_y = SymbolicExpr::add(
        SymbolicExpr::multiply(r_u[2], r_v[0]),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(r_u[0], r_v[2])));
    cross_y = cross_y->simplify();

    auto cross_z = SymbolicExpr::add(
        SymbolicExpr::multiply(r_u[0], r_v[1]),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(r_u[1], r_v[0])));
    cross_z = cross_z->simplify();

    return VectorField{cross_x, cross_y, cross_z};
}

// ============================================================
/// 第一类曲面积分 (Requirement 49.1, 49.3)
// ============================================================

std::shared_ptr<SymbolicExpr> surface_integral_scalar(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper)
{
    if (!f) {
        throw std::invalid_argument("surface_integral_scalar: f must not be null");
    }
    if (parametrization.size() != 3) {
        throw std::invalid_argument(
            "surface_integral_scalar: parametrization must have 3 components");
    }
    if (!u_lower || !u_upper || !v_lower || !v_upper) {
        throw std::invalid_argument("surface_integral_scalar: bounds must not be null");
    }

    std::vector<std::string> coord_vars = {"x", "y", "z"};

    /// 将 f 中的坐标变量替换为参数化表达式: f(r(u,v))
    auto f_composed = f;
    for (size_t i = 0; i < 3; ++i) {
        f_composed = f_composed->substitute(coord_vars[i], parametrization[i]);
    }
    f_composed = f_composed->simplify();

    /// 计算 r_u × r_v
    auto cross = vector_calculus_cross_product_partials(parametrization, u, v);

    /// |r_u × r_v| = √(cross_x² + cross_y² + cross_z²)
    std::shared_ptr<SymbolicExpr> mag_sq = nullptr;
    for (size_t i = 0; i < 3; ++i) {
        auto sq = SymbolicExpr::power(cross[i], SymbolicExpr::number(2));
        sq = sq->simplify();
        if (!mag_sq) {
            mag_sq = sq;
        } else {
            mag_sq = SymbolicExpr::add(mag_sq, sq);
            mag_sq = mag_sq->simplify();
        }
    }

    auto magnitude = SymbolicExpr::sqrt(mag_sq);
    magnitude = magnitude->simplify();

    /// 被积函数: f(r(u,v)) · |r_u × r_v|
    auto integrand = SymbolicExpr::multiply(f_composed, magnitude);
    integrand = integrand->simplify();

    /// 使用 MultipleIntegralEngine 计算二重积分
    MultipleIntegralEngine engine;
    Integrator integrator;

    std::vector<MultipleIntegralEngine::IntegrationStep> steps;
    steps.push_back({v, v_lower, v_upper});  // 内层积分
    steps.push_back({u, u_lower, u_upper});  // 外层积分

    auto result = engine.evaluate(*integrand, steps, integrator);
    if (result) {
        auto simplified = result->simplify();
        return simplified ? simplified : result;
    }
    return nullptr;
}

// ============================================================
/// 第二类曲面积分 (Requirement 49.2, 49.3)
// ============================================================

std::shared_ptr<SymbolicExpr> surface_integral_vector(
    const VectorField& F, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper)
{
    if (F.size() != 3) {
        throw std::invalid_argument(
            "surface_integral_vector: F must have 3 components");
    }
    if (parametrization.size() != 3) {
        throw std::invalid_argument(
            "surface_integral_vector: parametrization must have 3 components");
    }
    if (!u_lower || !u_upper || !v_lower || !v_upper) {
        throw std::invalid_argument("surface_integral_vector: bounds must not be null");
    }

    std::vector<std::string> coord_vars = {"x", "y", "z"};

    /// 计算 r_u × r_v
    auto cross = vector_calculus_cross_product_partials(parametrization, u, v);

    /// 将 F 中的坐标变量替换为参数化表达式，然后计算 F·(r_u × r_v)
    std::shared_ptr<SymbolicExpr> dot_product = nullptr;
    for (size_t i = 0; i < 3; ++i) {
        if (!F[i]) continue;

        /// 将 Fᵢ 中的坐标变量替换为参数化表达式
        auto Fi_composed = F[i];
        for (size_t j = 0; j < 3; ++j) {
            Fi_composed = Fi_composed->substitute(coord_vars[j], parametrization[j]);
        }
        Fi_composed = Fi_composed->simplify();

        /// Fᵢ(r(u,v)) · (r_u × r_v)ᵢ
        auto term = SymbolicExpr::multiply(Fi_composed, cross[i]);
        term = term->simplify();

        if (!dot_product) {
            dot_product = term;
        } else {
            dot_product = SymbolicExpr::add(dot_product, term);
            dot_product = dot_product->simplify();
        }
    }

    if (!dot_product) {
        return SymbolicExpr::number(0);
    }

    /// 使用 MultipleIntegralEngine 计算二重积分
    MultipleIntegralEngine engine;
    Integrator integrator;

    std::vector<MultipleIntegralEngine::IntegrationStep> steps;
    steps.push_back({v, v_lower, v_upper});  // 内层积分
    steps.push_back({u, u_lower, u_upper});  // 外层积分

    auto result = engine.evaluate(*dot_product, steps, integrator);
    if (result) {
        auto simplified = result->simplify();
        return simplified ? simplified : result;
    }
    return nullptr;
}

// ============================================================
/// 格林定理 (Requirements 50.1, 89.1, 89.2, 89.3)
// ============================================================

std::shared_ptr<SymbolicExpr> greens_theorem(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds)
{
    if (!P || !Q) {
        throw std::invalid_argument("greens_theorem: P and Q must not be null");
    }
    if (vars.size() != 2) {
        throw std::invalid_argument("greens_theorem: vars must have exactly 2 elements");
    }
    if (!x_bounds.first || !x_bounds.second ||
        !y_bounds.first || !y_bounds.second) {
        throw std::invalid_argument("greens_theorem: bounds must not be null");
    }

    const std::string& x_var = vars[0];
    const std::string& y_var = vars[1];

    /// 计算被积函数: ∂Q/∂x - ∂P/∂y
    auto dQ_dx = Q->differentiate(x_var);
    auto dP_dy = P->differentiate(y_var);

    if (!dQ_dx) dQ_dx = SymbolicExpr::number(0);
    if (!dP_dy) dP_dy = SymbolicExpr::number(0);

    dQ_dx = dQ_dx->simplify();
    dP_dy = dP_dy->simplify();

    auto integrand = SymbolicExpr::add(
        dQ_dx,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dP_dy));
    integrand = integrand->simplify();

    /// 使用 MultipleIntegralEngine 计算二重积分 ∬(∂Q/∂x - ∂P/∂y) dA
    MultipleIntegralEngine engine;
    Integrator integrator;

    std::vector<MultipleIntegralEngine::IntegrationStep> steps;
    steps.push_back({y_var, y_bounds.first, y_bounds.second});  // 内层积分 (y)
    steps.push_back({x_var, x_bounds.first, x_bounds.second});  // 外层积分 (x)

    auto result = engine.evaluate(*integrand, steps, integrator);
    if (result) {
        auto simplified = result->simplify();
        return simplified ? simplified : result;
    }
    return nullptr;
}

// ============================================================
/// 格林定理面积公式 (Requirement 50.4, 89.2)
// ============================================================

std::shared_ptr<SymbolicExpr> greens_theorem_area(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    if (parametrization.size() != 2) {
        throw std::invalid_argument(
            "greens_theorem_area: parametrization must have 2 components (x(t), y(t))");
    }
    if (!a || !b) {
        throw std::invalid_argument("greens_theorem_area: bounds must not be null");
    }

    auto x_t = parametrization[0];
    auto y_t = parametrization[1];

    if (!x_t || !y_t) {
        throw std::invalid_argument("greens_theorem_area: parametrization components must not be null");
    }

    /// A = (1/2) ∮ (x dy - y dx)
    // = (1/2) ∫ₐᵇ (x(t)·y'(t) - y(t)·x'(t)) dt
    auto dx_dt = x_t->differentiate(t);
    auto dy_dt = y_t->differentiate(t);

    if (!dx_dt) dx_dt = SymbolicExpr::number(0);
    if (!dy_dt) dy_dt = SymbolicExpr::number(0);

    dx_dt = dx_dt->simplify();
    dy_dt = dy_dt->simplify();

    /// x(t)·y'(t) - y(t)·x'(t)
    auto term1 = SymbolicExpr::multiply(x_t, dy_dt);
    auto term2 = SymbolicExpr::multiply(y_t, dx_dt);
    auto integrand = SymbolicExpr::add(
        term1,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), term2));
    integrand = integrand->simplify();

    /// 计算定积分
    auto integral_result = vector_calculus_integrate_with_fallback(integrand, t, a, b);
    if (!integral_result) {
        return nullptr;
    }

    /// 乘以 1/2
    auto half = std::make_shared<SymbolicExpr>(
        std::make_shared<NumberNode>(Rational(1, 2)));
    auto area = SymbolicExpr::multiply(half, integral_result);
    area = area->simplify();

    return area;
}

// ============================================================
/// 散度定理 / 高斯定理 (Requirements 50.2, 92.1, 92.2, 92.3)
// ============================================================

std::shared_ptr<SymbolicExpr> divergence_theorem(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& z_bounds)
{
    if (F.size() != 3) {
        throw std::invalid_argument("divergence_theorem: F must have 3 components");
    }
    if (vars.size() != 3) {
        throw std::invalid_argument("divergence_theorem: vars must have exactly 3 elements");
    }
    if (!x_bounds.first || !x_bounds.second ||
        !y_bounds.first || !y_bounds.second ||
        !z_bounds.first || !z_bounds.second) {
        throw std::invalid_argument("divergence_theorem: bounds must not be null");
    }

    /// 计算散度 ∇·F = ∂F₁/∂x + ∂F₂/∂y + ∂F₃/∂z
    auto div_F = divergence(F, vars);
    if (!div_F) {
        return SymbolicExpr::number(0);
    }

    /// 使用 MultipleIntegralEngine 计算三重积分 ∭_V ∇·F dV
    MultipleIntegralEngine engine;
    Integrator integrator;

    std::vector<MultipleIntegralEngine::IntegrationStep> steps;
    steps.push_back({vars[2], z_bounds.first, z_bounds.second});  // 内层积分 (z)
    steps.push_back({vars[1], y_bounds.first, y_bounds.second});  // 中层积分 (y)
    steps.push_back({vars[0], x_bounds.first, x_bounds.second});  // 外层积分 (x)

    auto result = engine.evaluate(*div_F, steps, integrator);
    if (result) {
        auto simplified = result->simplify();
        return simplified ? simplified : result;
    }
    return nullptr;
}

// ============================================================
/// 斯托克斯定理 (Requirements 50.3, 93.1, 93.2, 93.3)
// ============================================================

std::shared_ptr<SymbolicExpr> stokes_theorem(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& v_bounds)
{
    if (F.size() != 3) {
        throw std::invalid_argument("stokes_theorem: F must have 3 components");
    }
    if (vars.size() != 3) {
        throw std::invalid_argument("stokes_theorem: vars must have exactly 3 elements");
    }
    if (parametrization.size() != 3) {
        throw std::invalid_argument(
            "stokes_theorem: parametrization must have 3 components");
    }
    if (!u_bounds.first || !u_bounds.second ||
        !v_bounds.first || !v_bounds.second) {
        throw std::invalid_argument("stokes_theorem: bounds must not be null");
    }

    /// 计算旋度 ∇×F
    auto curl_F = curl(F, vars);

    /// 计算 r_u × r_v（曲面法向量）
    auto cross = vector_calculus_cross_product_partials(parametrization, u, v);

    /// 将 (∇×F) 中的坐标变量替换为参数化表达式，然后计算 (∇×F)·(r_u × r_v)
    std::vector<std::string> coord_vars = {"x", "y", "z"};
    std::shared_ptr<SymbolicExpr> dot_product = nullptr;

    for (size_t i = 0; i < 3; ++i) {
        if (!curl_F[i]) continue;

        /// 将 (∇×F)ᵢ 中的坐标变量替换为参数化表达式
        auto curl_i_composed = curl_F[i];
        for (size_t j = 0; j < 3; ++j) {
            curl_i_composed = curl_i_composed->substitute(coord_vars[j], parametrization[j]);
        }
        curl_i_composed = curl_i_composed->simplify();

        /// (∇×F)ᵢ(r(u,v)) · (r_u × r_v)ᵢ
        auto term = SymbolicExpr::multiply(curl_i_composed, cross[i]);
        term = term->simplify();

        if (!dot_product) {
            dot_product = term;
        } else {
            dot_product = SymbolicExpr::add(dot_product, term);
            dot_product = dot_product->simplify();
        }
    }

    if (!dot_product) {
        return SymbolicExpr::number(0);
    }

    /// 使用 MultipleIntegralEngine 计算二重积分
    MultipleIntegralEngine engine;
    Integrator integrator;

    std::vector<MultipleIntegralEngine::IntegrationStep> steps;
    steps.push_back({v, v_bounds.first, v_bounds.second});  // 内层积分 (v)
    steps.push_back({u, u_bounds.first, u_bounds.second});  // 外层积分 (u)

    auto result = engine.evaluate(*dot_product, steps, integrator);
    if (result) {
        auto simplified = result->simplify();
        return simplified ? simplified : result;
    }
    return nullptr;
}

// ============================================================
/// 多元极值 (Requirements 44, 51, 91)
// ============================================================

/**
 * @internal
 * @brief 将海森矩阵在临界点处求值，返回数值矩阵元素。
 *
 * @param[in] H    海森矩阵表达式
 * @param[in] vars 变量名列表
 * @param[in] pt   临界点坐标映射
 * @param[in] n    矩阵维度
 * @param[out] numeric_H 数值化的海森矩阵（行优先）
 * @return 是否所有元素都成功数值化
 */
static bool vector_calculus_evaluate_hessian_at_point(
    const std::shared_ptr<SymbolicExpr>& H,
    const std::vector<std::string>& vars,
    const std::map<std::string, std::shared_ptr<SymbolicExpr>>& pt,
    size_t n,
    std::vector<double>& numeric_H)
{
    auto mat_node = std::dynamic_pointer_cast<MatrixNode>(H->root);
    if (!mat_node || mat_node->rows != n || mat_node->cols != n) {
        return false;
    }

    numeric_H.resize(n * n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            auto elem_node = mat_node->get(i, j);
            if (!elem_node) {
                numeric_H[i * n + j] = 0.0;
                continue;
            }
            auto elem = std::make_shared<SymbolicExpr>(elem_node);
            /// 代入临界点坐标
            for (const auto& [var_name, val] : pt) {
                elem = elem->substitute(var_name, val);
                if (!elem) return false;
            }
            elem = elem->simplify();
            if (!elem) return false;

            if (elem->is_number()) {
                try {
                    numeric_H[i * n + j] = elem->to_numeric();
                } catch (...) {
                    return false;
                }
            } else {
                /// 尝试直接数值求值
                try {
                    numeric_H[i * n + j] = elem->to_numeric();
                } catch (...) {
                    return false;
                }
            }
        }
    }
    return true;
}

/**
 * @internal
 * @brief 使用特征值分析对临界点进行分类。
 *
 * 对 n×n 实对称矩阵，直接通过数值分析特征值符号判断正定性。
 * 对于小矩阵（n ≤ 3），使用解析公式；对于大矩阵，使用特征多项式求解。
 *
 * @param[in] numeric_H 数值化的海森矩阵（行优先，n×n）
 * @param[in] n         矩阵维度
 * @return 分类字符串: "minimum", "maximum", "saddle", "degenerate"
 */
static std::string vector_calculus_classify_critical_point(
    const std::vector<double>& numeric_H, size_t n)
{
    const double tol = 1e-10;

    /// 1×1 情况：直接判断
    if (n == 1) {
        double val = numeric_H[0];
        if (std::abs(val) < tol) return "degenerate";
        if (val > 0) return "minimum";
        return "maximum";
    }

    /// 2×2 情况：使用行列式和迹直接判断
    if (n == 2) {
        double a = numeric_H[0], b = numeric_H[1];
        double c = numeric_H[2], d = numeric_H[3];
        double det = a * d - b * c;
        double trace = a + d;

        if (std::abs(det) < tol) return "degenerate";
        if (det > 0 && trace > 0) return "minimum";
        if (det > 0 && trace < 0) return "maximum";
        return "saddle";
    }

    /// 一般情况：构建符号矩阵并求特征值
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            /// 使用 Rational 避免浮点精度问题
            double val = numeric_H[i * n + j];
            int int_val = static_cast<int>(std::round(val));
            if (std::abs(val - int_val) < tol) {
                grid[i][j] = SymbolicExpr::number(int_val);
            } else {
                grid[i][j] = SymbolicExpr::number(val);
            }
        }
    }
    auto H_mat = SymbolicExpr::matrix(grid);

    /// 计算特征多项式并求解特征值
    auto cp = SymbolicExpr::charpoly(H_mat, "lambda");
    if (!cp) return "degenerate";

    auto eigenvals = SymbolicExpr::solve(cp, "lambda");
    if (eigenvals.empty()) return "degenerate";

    bool all_positive = true;
    bool all_negative = true;
    bool has_zero = false;

    for (const auto& ev : eigenvals) {
        if (!ev) {
            return "degenerate";
        }
        auto ev_simplified = ev->simplify();
        double val = 0.0;
        try {
            val = ev_simplified->to_numeric();
        } catch (...) {
            return "degenerate";
        }

        if (std::abs(val) < tol) {
            has_zero = true;
            all_positive = false;
            all_negative = false;
        } else if (val > 0) {
            all_negative = false;
        } else {
            all_positive = false;
        }
    }

    if (has_zero) return "degenerate";
    if (all_positive) return "minimum";
    if (all_negative) return "maximum";
    return "saddle";
}

std::vector<CriticalPoint> find_extrema(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars)
{
    if (!f || vars.empty()) {
        return {};
    }

    size_t n = vars.size();

    /// 计算梯度 ∇f
    std::vector<std::shared_ptr<SymbolicExpr>> grad_eqs;
    grad_eqs.reserve(n);
    for (const auto& var : vars) {
        auto partial = f->differentiate(var);
        if (partial) {
            partial = partial->simplify();
        }
        grad_eqs.push_back(partial);
    }

    /// 求解 ∇f = 0 系统（使用多项式系统求解器）
    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> solutions;

    std::vector<SymbolicExpr> poly_eqs;
    poly_eqs.reserve(n);
    for (const auto& eq : grad_eqs) {
        if (eq) poly_eqs.push_back(*eq);
    }
    if (poly_eqs.size() == n) {
        auto poly_solutions = Solver::solve_polynomial_system(poly_eqs, vars);
        for (const auto& sol : poly_solutions) {
            std::map<std::string, std::shared_ptr<SymbolicExpr>> pt;
            for (const auto& [name, val] : sol) {
                pt[name] = std::make_shared<SymbolicExpr>(val);
            }
            solutions.push_back(pt);
        }
    }

    /// 如果多项式求解器失败，尝试线性求解器
    if (solutions.empty()) {
        solutions = SymbolicExpr::solve_system(grad_eqs, vars);
    }

    if (solutions.empty()) {
        return {};
    }

    /// 计算海森矩阵
    auto H = hessian(f, vars);
    if (!H) return {};

    /// 对每个临界点进行分类
    std::vector<CriticalPoint> result;
    result.reserve(solutions.size());

    for (const auto& sol : solutions) {
        CriticalPoint cp;
        cp.point = sol;

        /// 在临界点处求值海森矩阵
        std::vector<double> numeric_H;
        if (vector_calculus_evaluate_hessian_at_point(H, vars, sol, n, numeric_H)) {
            cp.classification = vector_calculus_classify_critical_point(numeric_H, n);
        } else {
            cp.classification = "degenerate";
        }

        result.push_back(std::move(cp));
    }

    return result;
}

// ============================================================
/// 拉格朗日乘数法 (Requirement 44)
// ============================================================

std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> lagrange_multipliers(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars)
{
    if (!f || constraints.empty() || vars.empty()) {
        return {};
    }

    size_t n = vars.size();
    size_t m = constraints.size();

    /// 计算目标函数的梯度
    VectorField grad_f;
    grad_f.reserve(n);
    for (const auto& var : vars) {
        auto partial = f->differentiate(var);
        if (partial) partial = partial->simplify();
        grad_f.push_back(partial);
    }

    /// 计算每个约束的梯度
    std::vector<VectorField> grad_constraints;
    grad_constraints.reserve(m);
    for (const auto& g : constraints) {
        VectorField grad_g;
        grad_g.reserve(n);
        for (const auto& var : vars) {
            auto partial = g->differentiate(var);
            if (partial) partial = partial->simplify();
            grad_g.push_back(partial);
        }
        grad_constraints.push_back(std::move(grad_g));
    }

    /// 构造乘数变量名
    std::vector<std::string> lambda_names;
    lambda_names.reserve(m);
    for (size_t k = 0; k < m; ++k) {
        lambda_names.push_back("lambda_" + std::to_string(k + 1));
    }

    /// 构造方程组: ∂f/∂xᵢ - ∑ λₖ · ∂gₖ/∂xᵢ = 0 (对每个 i)
    /// 加上约束方程: gₖ = 0 (对每个 k)
    std::vector<std::shared_ptr<SymbolicExpr>> equations;
    equations.reserve(n + m);

    for (size_t i = 0; i < n; ++i) {
        /// ∂f/∂xᵢ - λ₁·∂g₁/∂xᵢ - λ₂·∂g₂/∂xᵢ - ...
        auto eq = grad_f[i];
        for (size_t k = 0; k < m; ++k) {
            auto lambda_var = SymbolicExpr::variable(lambda_names[k]);
            auto term = SymbolicExpr::multiply(lambda_var, grad_constraints[k][i]);
            eq = SymbolicExpr::add(eq, SymbolicExpr::multiply(SymbolicExpr::number(-1), term));
        }
        if (eq) eq = eq->simplify();
        equations.push_back(eq);
    }

    /// 添加约束方程
    for (const auto& g : constraints) {
        equations.push_back(g);
    }

    /// 构造所有未知数列表: 原始变量 + 乘数
    std::vector<std::string> all_vars;
    all_vars.reserve(n + m);
    for (const auto& v : vars) {
        all_vars.push_back(v);
    }
    for (const auto& lam : lambda_names) {
        all_vars.push_back(lam);
    }

    /// 使用多项式系统求解器求解增广系统
    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> solutions;

    std::vector<SymbolicExpr> poly_eqs;
    poly_eqs.reserve(equations.size());
    for (const auto& eq : equations) {
        if (eq) poly_eqs.push_back(*eq);
    }

    if (poly_eqs.size() == equations.size()) {
        auto poly_solutions = Solver::solve_polynomial_system(poly_eqs, all_vars);
        for (const auto& sol : poly_solutions) {
            std::map<std::string, std::shared_ptr<SymbolicExpr>> pt;
            for (const auto& [name, val] : sol) {
                pt[name] = std::make_shared<SymbolicExpr>(val);
            }
            solutions.push_back(pt);
        }
    }

    /// 如果多项式求解器失败，尝试线性求解器
    if (solutions.empty()) {
        solutions = SymbolicExpr::solve_system(equations, all_vars);
    }

    /// 从解中提取仅原始变量（去除乘数）
    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> result;
    std::set<std::string> var_set(vars.begin(), vars.end());

    for (const auto& sol : solutions) {
        std::map<std::string, std::shared_ptr<SymbolicExpr>> filtered;
        for (const auto& [name, val] : sol) {
            if (var_set.count(name)) {
                filtered[name] = val;
            }
        }
        if (!filtered.empty()) {
            result.push_back(std::move(filtered));
        }
    }

    return result;
}

// ============================================================
/// 向量代数运算 (Requirements 40, 41, 42)
// ============================================================

std::shared_ptr<SymbolicExpr> dot_product(const VectorField& a, const VectorField& b)
{
    if (a.size() != b.size()) {
        throw std::invalid_argument("dot_product: vectors must have the same dimension");
    }
    std::shared_ptr<SymbolicExpr> sum = SymbolicExpr::number(0);
    for (size_t i = 0; i < a.size(); ++i) {
        if (!a[i] || !b[i]) continue;
        sum = SymbolicExpr::add(sum, SymbolicExpr::multiply(a[i], b[i]));
    }
    return sum->simplify();
}

VectorField cross_product(const VectorField& a, const VectorField& b)
{
    if (a.size() != 3 || b.size() != 3) {
        throw std::invalid_argument("cross_product: vectors must be 3-dimensional");
    }
    auto sub = [](const std::shared_ptr<SymbolicExpr>& p,
                  const std::shared_ptr<SymbolicExpr>& q) {
        return SymbolicExpr::add(p, SymbolicExpr::multiply(SymbolicExpr::number(-1), q));
    };
    VectorField result(3);
    result[0] = sub(SymbolicExpr::multiply(a[1], b[2]), SymbolicExpr::multiply(a[2], b[1]))->simplify();
    result[1] = sub(SymbolicExpr::multiply(a[2], b[0]), SymbolicExpr::multiply(a[0], b[2]))->simplify();
    result[2] = sub(SymbolicExpr::multiply(a[0], b[1]), SymbolicExpr::multiply(a[1], b[0]))->simplify();
    return result;
}

VectorField vector_project(const VectorField& a, const VectorField& b)
{
    if (a.size() != b.size()) {
        throw std::invalid_argument("vector_project: vectors must have the same dimension");
    }
    auto bb = dot_product(b, b);
    if (bb->root && bb->root->is_zero()) {
        return VectorField(a.size(), SymbolicExpr::number(0));
    }
    auto ab = dot_product(a, b);
    auto coeff = SymbolicExpr::divide(ab, bb);
    VectorField result;
    result.reserve(b.size());
    for (const auto& comp : b) {
        result.push_back(SymbolicExpr::multiply(coeff, comp)->simplify());
    }
    return result;
}

std::shared_ptr<SymbolicExpr> scalar_project(const VectorField& a, const VectorField& b)
{
    auto bb = dot_product(b, b);
    if (bb->root && bb->root->is_zero()) {
        return nullptr;
    }
    auto ab = dot_product(a, b);
    return SymbolicExpr::divide(ab, SymbolicExpr::sqrt(bb))->simplify();
}

std::shared_ptr<SymbolicExpr> vector_angle_symbolic(const VectorField& a, const VectorField& b)
{
    auto aa = dot_product(a, a);
    auto bb = dot_product(b, b);
    if ((aa->root && aa->root->is_zero()) || (bb->root && bb->root->is_zero())) {
        return nullptr;
    }
    auto ab = dot_product(a, b);
    auto denom = SymbolicExpr::multiply(SymbolicExpr::sqrt(aa), SymbolicExpr::sqrt(bb));
    auto cos_theta = SymbolicExpr::divide(ab, denom);
    auto arccos_node = std::make_shared<FunctionNode>(
        FunctionNode::FuncType::ArcCos,
        std::vector<std::shared_ptr<SymbolicNode>>{cos_theta->root});
    return std::make_shared<SymbolicExpr>(arccos_node)->simplify();
}

std::shared_ptr<SymbolicExpr> mixed_product(const VectorField& a, const VectorField& b,
    const VectorField& c)
{
    return dot_product(a, cross_product(b, c));
}

} // namespace lamina
