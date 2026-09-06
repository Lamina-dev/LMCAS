#include "internal/vector_calculus_support.hpp"
#include "integration.hpp"
#include "numeric_evaluation.hpp"
#include "solver.hpp"
#include "symbolic_ast.hpp"

#include <cmath>
#include <exception>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace LMCAS {

using namespace vector_calculus_detail;

VectorCalculusFieldResult gradient_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    const std::string operation = "gradient";
    auto valid = vector_calculus_validate_expr_vars(f, vars, context, operation);
    if (!valid) return VectorCalculusFieldResult::failure(valid.error());
    auto step = context.consume_steps(vars.size(), operation);
    if (!step) return VectorCalculusFieldResult::failure(step.error());
    try {
        return vector_calculus_wrap_field(gradient(f, vars), operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusFieldResult::failure(CasErrc::ResourceLimit,
                                                  "vector calculus allocation failed",
                                                  operation);
    } catch (const std::exception& e) {
        return VectorCalculusFieldResult::failure(CasErrc::InternalInvariant,
                                                  e.what(), operation);
    }
}

VectorCalculusFieldResult gradient_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return gradient_checked(f, vars, context);
}

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


VectorCalculusExprResult divergence_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    const std::string operation = "divergence";
    auto valid = vector_calculus_validate_field_vars(F, vars, context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());
    auto step = context.consume_steps(F.size(), operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());
    try {
        return vector_calculus_wrap_expr(divergence(F, vars), operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "vector calculus allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult divergence_checked(
    const VectorField& F,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return divergence_checked(F, vars, context);
}

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


VectorCalculusFieldResult curl_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    const std::string operation = "curl";
    auto valid = vector_calculus_validate_field_vars(F, vars, context, operation);
    if (!valid) return VectorCalculusFieldResult::failure(valid.error());
    if (F.size() != 2 && F.size() != 3) {
        return VectorCalculusFieldResult::failure(
            CasErrc::InvalidArgument,
            "curl requires a 2D or 3D vector field",
            operation);
    }
    auto step = context.consume_steps(F.size(), operation);
    if (!step) return VectorCalculusFieldResult::failure(step.error());
    try {
        return vector_calculus_wrap_field(curl(F, vars), operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusFieldResult::failure(CasErrc::ResourceLimit,
                                                  "vector calculus allocation failed",
                                                  operation);
    } catch (const std::exception& e) {
        return VectorCalculusFieldResult::failure(CasErrc::InternalInvariant,
                                                  e.what(), operation);
    }
}

VectorCalculusFieldResult curl_checked(
    const VectorField& F,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return curl_checked(F, vars, context);
}

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


VectorCalculusExprResult laplacian_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    const std::string operation = "laplacian";
    auto valid = vector_calculus_validate_expr_vars(f, vars, context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());
    auto step = context.consume_steps(vars.size() * 2, operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());
    try {
        return vector_calculus_wrap_expr(laplacian(f, vars), operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "vector calculus allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult laplacian_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return laplacian_checked(f, vars, context);
}

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

VectorCalculusExprResult directional_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    const VectorField& direction,
    int order,
    ComputationContext& context)
{
    const std::string operation = "directional_derivative";
    auto valid = vector_calculus_validate_expr_vars(f, vars, context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());
    auto dir_valid = vector_calculus_validate_field_vars(direction, vars, context, operation);
    if (!dir_valid) return VectorCalculusExprResult::failure(dir_valid.error());
    if (order < 1) {
        return VectorCalculusExprResult::failure(CasErrc::InvalidArgument,
                                                 "directional derivative order must be positive",
                                                 operation);
    }
    auto mag_sq = vector_calculus_magnitude_squared(direction);
    if (!mag_sq || !LMCAS::detail::node(mag_sq)) {
        return VectorCalculusExprResult::failure(
            CasErrc::InternalInvariant,
            "direction magnitude construction failed",
            operation);
    }
    if (mag_sq->is_zero()) {
        return VectorCalculusExprResult::failure(CasErrc::DomainError,
                                                 "direction vector cannot be zero",
                                                 operation);
    }
    auto step = context.consume_steps(vars.size() * static_cast<std::size_t>(order),
                                      operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());
    try {
        return vector_calculus_wrap_expr(
            directional_derivative(f, vars, direction, order), operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "vector calculus allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult directional_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    const VectorField& direction,
    int order)
{
    ComputationContext context;
    return directional_derivative_checked(f, vars, direction, order, context);
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


VectorCalculusExprResult jacobian_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    const std::string operation = "jacobian";
    auto valid = vector_calculus_validate_functions_vars(functions, vars,
                                                        context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());
    auto step = context.consume_steps(functions.size() * vars.size(), operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());
    try {
        return vector_calculus_wrap_expr(jacobian(functions, vars), operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "vector calculus allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult jacobian_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return jacobian_checked(functions, vars, context);
}

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


VectorCalculusExprResult hessian_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    const std::string operation = "hessian";
    auto valid = vector_calculus_validate_expr_vars(f, vars, context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());
    auto step = context.consume_steps(vars.size() * vars.size(), operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());
    try {
        return vector_calculus_wrap_expr(hessian(f, vars), operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "vector calculus allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult hessian_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return hessian_checked(f, vars, context);
}

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
} // namespace LMCAS
