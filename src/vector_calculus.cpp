/**
 * @file vector_calculus.cpp
 * @brief 向量微积分模块实现：梯度、散度、旋度、拉普拉斯算子、方向导数、曲线积分、曲面积分、极值分析。
 */

#include "vector_calculus.hpp"
#include "integration.hpp"
#include "numeric_evaluation.hpp"
#include "solver.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <exception>
#include <cmath>
#include <map>
#include <set>

namespace lamina {

namespace {

bool vector_calculus_checked_finite_numeric(
    const std::shared_ptr<SymbolicExpr>& expr,
    double& value,
    ComputationContext* context = nullptr)
{
    if (!expr || !lamina::detail::node(expr)) return false;
    ComputationContext local_context;
    ComputationContext& eval_context = context ? *context : local_context;
    auto evaluated = evaluate_numeric(*expr, NumericBindings{}, eval_context);
    if (!evaluated || !evaluated.value().is_finite()) return false;
    value = evaluated.value().value;
    return std::isfinite(value);
}

bool vector_calculus_contains_unevaluated_integral(
    const std::shared_ptr<const SymbolicNode>& node,
    std::size_t depth = 0)
{
    if (!node || depth > 200) return false;
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (func->type() == FunctionNode::FuncType::Calculus_Integral) {
            return true;
        }
        for (const auto& arg : func->arguments()) {
            if (vector_calculus_contains_unevaluated_integral(arg, depth + 1)) {
                return true;
            }
        }
        return false;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& op : add->operands()) {
            if (vector_calculus_contains_unevaluated_integral(op, depth + 1)) {
                return true;
            }
        }
        return false;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            if (vector_calculus_contains_unevaluated_integral(op, depth + 1)) {
                return true;
            }
        }
        return false;
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return vector_calculus_contains_unevaluated_integral(pow->base(), depth + 1) ||
               vector_calculus_contains_unevaluated_integral(pow->exponent(), depth + 1);
    }
    if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        return vector_calculus_contains_unevaluated_integral(complex->real(), depth + 1) ||
               vector_calculus_contains_unevaluated_integral(complex->imag(), depth + 1);
    }
    if (auto matrix = std::dynamic_pointer_cast<const MatrixNode>(node)) {
        for (size_t r = 0; r < matrix->rows(); ++r) {
            for (size_t c = 0; c < matrix->cols(); ++c) {
                if (vector_calculus_contains_unevaluated_integral(
                        matrix->get(r, c), depth + 1)) {
                    return true;
                }
            }
        }
        return false;
    }
    if (auto rel = std::dynamic_pointer_cast<const RelationalNode>(node)) {
        return vector_calculus_contains_unevaluated_integral(rel->left(), depth + 1) ||
               vector_calculus_contains_unevaluated_integral(rel->right(), depth + 1);
    }
    if (auto logical = std::dynamic_pointer_cast<const LogicalNode>(node)) {
        return vector_calculus_contains_unevaluated_integral(logical->left(), depth + 1) ||
               vector_calculus_contains_unevaluated_integral(logical->right(), depth + 1);
    }
    if (auto piecewise = std::dynamic_pointer_cast<const PiecewiseNode>(node)) {
        for (const auto& branch : piecewise->branches()) {
            if (vector_calculus_contains_unevaluated_integral(branch.expression, depth + 1) ||
                vector_calculus_contains_unevaluated_integral(branch.condition, depth + 1)) {
                return true;
            }
        }
        return vector_calculus_contains_unevaluated_integral(
            piecewise->default_expr(), depth + 1);
    }
    if (auto sum = std::dynamic_pointer_cast<const SummationNode>(node)) {
        return vector_calculus_contains_unevaluated_integral(sum->body(), depth + 1) ||
               vector_calculus_contains_unevaluated_integral(sum->lower_bound(), depth + 1) ||
               vector_calculus_contains_unevaluated_integral(sum->upper_bound(), depth + 1);
    }
    if (auto product = std::dynamic_pointer_cast<const ProductNode_Op>(node)) {
        return vector_calculus_contains_unevaluated_integral(product->body(), depth + 1) ||
               vector_calculus_contains_unevaluated_integral(product->lower_bound(), depth + 1) ||
               vector_calculus_contains_unevaluated_integral(product->upper_bound(), depth + 1);
    }
    if (auto transform = std::dynamic_pointer_cast<const TransformNode>(node)) {
        return vector_calculus_contains_unevaluated_integral(transform->body(), depth + 1);
    }
    if (auto quantifier = std::dynamic_pointer_cast<const QuantifierNode>(node)) {
        return vector_calculus_contains_unevaluated_integral(quantifier->domain(), depth + 1) ||
               vector_calculus_contains_unevaluated_integral(quantifier->predicate(), depth + 1);
    }
    if (auto set_builder = std::dynamic_pointer_cast<const SetBuilderNode>(node)) {
        return vector_calculus_contains_unevaluated_integral(set_builder->domain(), depth + 1) ||
               vector_calculus_contains_unevaluated_integral(set_builder->predicate(), depth + 1);
    }
    return false;
}

Result<void> vector_calculus_validate_expr_vars(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!f || !lamina::detail::node(f)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "expression cannot be null", operation);
    }
    if (vars.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "variable list cannot be empty", operation);
    }
    for (const auto& var : vars) {
        if (var.empty()) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "variable name cannot be empty",
                                         operation);
        }
    }
    return Result<void>::success();
}

Result<void> vector_calculus_validate_field_vars(
    const VectorField& field,
    const std::vector<std::string>& vars,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (field.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "vector field cannot be empty", operation);
    }
    if (field.size() != vars.size()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "vector field and variables must have the same dimension",
                                     operation);
    }
    for (const auto& var : vars) {
        if (var.empty()) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "variable name cannot be empty",
                                         operation);
        }
    }
    for (const auto& component : field) {
        if (!component || !lamina::detail::node(component)) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "vector field components cannot be null",
                                         operation);
        }
    }
    return Result<void>::success();
}

Result<void> vector_calculus_validate_functions_vars(
    const std::vector<std::shared_ptr<SymbolicExpr>>& functions,
    const std::vector<std::string>& vars,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (functions.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "function list cannot be empty",
                                     operation);
    }
    if (vars.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "variable list cannot be empty",
                                     operation);
    }
    for (const auto& var : vars) {
        if (var.empty()) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "variable name cannot be empty",
                                         operation);
        }
    }
    for (const auto& function : functions) {
        if (!function || !lamina::detail::node(function)) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "function expressions cannot be null",
                                         operation);
        }
    }
    return Result<void>::success();
}

Result<void> vector_calculus_validate_curve_parametrization(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (parametrization.size() != 2 && parametrization.size() != 3) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "curve parametrization must be two- or three-dimensional",
                                     operation);
    }
    if (t.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "parameter variable name cannot be empty",
                                     operation);
    }
    if (!a || !lamina::detail::node(a) || !b || !lamina::detail::node(b)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "curve integral bounds cannot be null",
                                     operation);
    }
    for (const auto& component : parametrization) {
        if (!component || !lamina::detail::node(component)) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "curve parametrization components cannot be null",
                                         operation);
        }
    }
    return Result<void>::success();
}

Result<void> vector_calculus_validate_curve_scalar_inputs(
    const std::shared_ptr<SymbolicExpr>& f,
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!f || !lamina::detail::node(f)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "scalar field cannot be null",
                                     operation);
    }
    return vector_calculus_validate_curve_parametrization(
        parametrization, t, a, b, context, operation);
}

Result<void> vector_calculus_validate_curve_vector_inputs(
    const VectorField& F,
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (F.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "vector field cannot be empty",
                                     operation);
    }
    if (F.size() != parametrization.size()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "vector field and parametrization must have the same dimension",
                                     operation);
    }
    for (const auto& component : F) {
        if (!component || !lamina::detail::node(component)) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "vector field components cannot be null",
                                         operation);
        }
    }
    return vector_calculus_validate_curve_parametrization(
        parametrization, t, a, b, context, operation);
}

Result<void> vector_calculus_validate_surface_parametrization(
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower,
    const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower,
    const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (parametrization.size() != 3) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "surface parametrization must be three-dimensional",
                                     operation);
    }
    if (u.empty() || v.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "surface parameter variable names cannot be empty",
                                     operation);
    }
    if (u == v) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "surface parameter variables must be distinct",
                                     operation);
    }
    if (!u_lower || !lamina::detail::node(u_lower) || !u_upper || !lamina::detail::node(u_upper) ||
        !v_lower || !lamina::detail::node(v_lower) || !v_upper || !lamina::detail::node(v_upper)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "surface integral bounds cannot be null",
                                     operation);
    }
    for (const auto& component : parametrization) {
        if (!component || !lamina::detail::node(component)) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "surface parametrization components cannot be null",
                                         operation);
        }
    }
    return Result<void>::success();
}

Result<void> vector_calculus_validate_surface_scalar_inputs(
    const std::shared_ptr<SymbolicExpr>& f,
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower,
    const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower,
    const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!f || !lamina::detail::node(f)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "scalar field cannot be null",
                                     operation);
    }
    return vector_calculus_validate_surface_parametrization(
        parametrization, u, v, u_lower, u_upper, v_lower, v_upper,
        context, operation);
}

Result<void> vector_calculus_validate_surface_vector_inputs(
    const VectorField& F,
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower,
    const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower,
    const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (F.size() != 3) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "surface vector field must have three components",
                                     operation);
    }
    for (const auto& component : F) {
        if (!component || !lamina::detail::node(component)) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "surface vector field components cannot be null",
                                         operation);
        }
    }
    return vector_calculus_validate_surface_parametrization(
        parametrization, u, v, u_lower, u_upper, v_lower, v_upper,
        context, operation);
}

Result<void> vector_calculus_validate_bound_pair(
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& bounds,
    ComputationContext& context,
    const std::string& operation,
    const std::string& label)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!bounds.first || !lamina::detail::node(bounds.first) || !bounds.second || !lamina::detail::node(bounds.second)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     label + " bounds cannot be null",
                                     operation);
    }
    return Result<void>::success();
}

Result<void> vector_calculus_validate_distinct_vars(
    const std::vector<std::string>& vars,
    std::size_t expected,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (vars.size() != expected) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "variable list has the wrong dimension",
                                     operation);
    }
    std::set<std::string> seen;
    for (const auto& var : vars) {
        if (var.empty()) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "variable name cannot be empty",
                                         operation);
        }
        if (!seen.insert(var).second) {
            return Result<void>::failure(CasErrc::InvalidArgument,
                                         "variable names must be distinct",
                                         operation);
        }
    }
    return Result<void>::success();
}

bool vector_calculus_expr_zero_after_substitution(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::map<std::string, std::shared_ptr<SymbolicExpr>>& point)
{
    if (!expr || !lamina::detail::node(expr)) return false;
    auto substituted = expr;
    for (const auto& [var, value] : point) {
        if (!value || !lamina::detail::node(value)) return false;
        substituted = substituted->substitute(var, value);
        if (!substituted || !lamina::detail::node(substituted)) return false;
    }
    substituted = substituted->simplify();
    if (!substituted || !lamina::detail::node(substituted)) return false;
    if (substituted->is_zero()) return true;
    double value = 0.0;
    return vector_calculus_checked_finite_numeric(substituted, value) &&
           std::abs(value) < 1e-8;
}

bool vector_calculus_point_has_vars(
    const std::map<std::string, std::shared_ptr<SymbolicExpr>>& point,
    const std::vector<std::string>& vars)
{
    for (const auto& var : vars) {
        auto it = point.find(var);
        if (it == point.end() || !it->second || !lamina::detail::node(it->second)) {
            return false;
        }
    }
    return true;
}

VectorCalculusExprResult vector_calculus_wrap_expr(
    std::shared_ptr<SymbolicExpr> expr,
    const std::string& operation)
{
    if (!expr || !lamina::detail::node(expr)) {
        return VectorCalculusExprResult::failure(
            CasErrc::InternalInvariant,
            "vector calculus construction produced a null expression",
            operation);
    }
    return VectorCalculusExprResult::success(std::move(expr));
}

VectorCalculusFieldResult vector_calculus_wrap_field(
    VectorField field,
    const std::string& operation)
{
    for (const auto& component : field) {
        if (!component || !lamina::detail::node(component)) {
            return VectorCalculusFieldResult::failure(
                CasErrc::InternalInvariant,
                "vector calculus construction produced a null component",
                operation);
        }
    }
    return VectorCalculusFieldResult::success(std::move(field));
}

} // namespace


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
    if (!mag_sq || !lamina::detail::node(mag_sq)) {
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

    if (vector_calculus_contains_unevaluated_integral(lamina::detail::node(result))) {
        return nullptr;
    }
    auto res = lamina::detail::make_expression_ptr(result);
    auto simplified = res->simplify();
    if (simplified && vector_calculus_contains_unevaluated_integral(lamina::detail::node(simplified))) {
        return nullptr;
    }
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
    ComputationContext context;
    double a_val = 0.0;
    double b_val = 0.0;
    if (!vector_calculus_checked_finite_numeric(a, a_val, &context) ||
        !vector_calculus_checked_finite_numeric(b, b_val, &context)) {
        return nullptr;
    }

    int n = 1000;
    double h = (b_val - a_val) / n;
    if (!std::isfinite(h)) return nullptr;
    double sum = 0.0;

    for (int i = 0; i <= n; ++i) {
        double xi = a_val + i * h;
        if (!std::isfinite(xi)) return nullptr;
        auto xi_expr = SymbolicExpr::number(xi);
        auto fi = integrand->substitute(var, xi_expr);
        double fi_val = 0.0;
        if (!vector_calculus_checked_finite_numeric(fi, fi_val, &context)) {
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
        if (!std::isfinite(sum)) return nullptr;
    }
    sum *= h / 3.0;
    if (!std::isfinite(sum)) return nullptr;

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

static VectorCalculusExprResult vector_calculus_inconclusive(
    const std::string& operation,
    const std::string& message)
{
    return VectorCalculusExprResult::failure(CasErrc::Inconclusive,
                                             message, operation);
}

static VectorCalculusExprResult vector_calculus_simplify_strict(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& operation,
    const std::string& message)
{
    if (!expr || !lamina::detail::node(expr)) {
        return vector_calculus_inconclusive(operation, message);
    }
    auto simplified = expr->simplify();
    if (!simplified || !lamina::detail::node(simplified)) {
        return vector_calculus_inconclusive(operation, message);
    }
    return VectorCalculusExprResult::success(std::move(simplified));
}

static VectorCalculusExprResult vector_calculus_differentiate_strict(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const std::string& operation)
{
    if (!expr || !lamina::detail::node(expr)) {
        return vector_calculus_inconclusive(
            operation, "expression derivative is outside the supported domain");
    }
    try {
        auto derivative = expr->differentiate(var);
        return vector_calculus_simplify_strict(
            derivative, operation,
            "expression derivative is outside the supported domain");
    } catch (const std::exception&) {
        return vector_calculus_inconclusive(
            operation, "expression derivative is outside the supported domain");
    }
}

static VectorCalculusExprResult vector_calculus_substitute_coords_strict(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::vector<std::string>& coord_vars,
    const VectorField& values,
    const std::string& operation)
{
    if (!expr || !lamina::detail::node(expr) || coord_vars.size() != values.size()) {
        return vector_calculus_inconclusive(
            operation, "coordinate substitution is outside the supported domain");
    }
    auto substituted = expr;
    for (size_t i = 0; i < coord_vars.size(); ++i) {
        if (!values[i] || !lamina::detail::node(values[i])) {
            return vector_calculus_inconclusive(
                operation, "coordinate substitution is outside the supported domain");
        }
        substituted = substituted->substitute(coord_vars[i], values[i]);
        if (!substituted || !lamina::detail::node(substituted)) {
            return vector_calculus_inconclusive(
                operation, "coordinate substitution is outside the supported domain");
        }
    }
    return vector_calculus_simplify_strict(
        substituted, operation,
        "coordinate substitution is outside the supported domain");
}

static VectorCalculusExprResult vector_calculus_definite_integral_strict(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::string& operation)
{
    auto result = vector_calculus_try_definite(integrand, var, a, b);
    if (!result || !lamina::detail::node(result)) {
        return vector_calculus_inconclusive(
            operation, "integral could not be evaluated exactly in the supported domain");
    }
    return VectorCalculusExprResult::success(std::move(result));
}

static VectorCalculusExprResult vector_calculus_multiple_integral_strict(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::vector<MultipleIntegralEngine::IntegrationStep>& steps,
    const std::string& operation)
{
    if (!integrand || !lamina::detail::node(integrand)) {
        return vector_calculus_inconclusive(
            operation, "integrand construction is outside the supported domain");
    }
    MultipleIntegralEngine engine;
    Integrator integrator;
    auto result = engine.evaluate(*integrand, steps, integrator);
    if (!result || !lamina::detail::node(result)) {
        return vector_calculus_inconclusive(
            operation, "integral could not be evaluated in the supported domain");
    }
    return vector_calculus_simplify_strict(
        result, operation,
        "integral result simplification is outside the supported domain");
}

static Result<VectorField> vector_calculus_cross_product_partials_strict(
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::string& operation)
{
    VectorField r_u;
    VectorField r_v;
    r_u.reserve(3);
    r_v.reserve(3);
    for (size_t i = 0; i < 3; ++i) {
        auto du = vector_calculus_differentiate_strict(
            parametrization[i], u, operation);
        if (!du) return Result<VectorField>::failure(du.error());
        auto dv = vector_calculus_differentiate_strict(
            parametrization[i], v, operation);
        if (!dv) return Result<VectorField>::failure(dv.error());
        r_u.push_back(std::move(du.value()));
        r_v.push_back(std::move(dv.value()));
    }

    auto cross_x = SymbolicExpr::add(
        SymbolicExpr::multiply(r_u[1], r_v[2]),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(r_u[2], r_v[1])));
    auto cross_x_checked = vector_calculus_simplify_strict(
        cross_x, operation, "surface normal construction is outside the supported domain");
    if (!cross_x_checked) return Result<VectorField>::failure(cross_x_checked.error());

    auto cross_y = SymbolicExpr::add(
        SymbolicExpr::multiply(r_u[2], r_v[0]),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(r_u[0], r_v[2])));
    auto cross_y_checked = vector_calculus_simplify_strict(
        cross_y, operation, "surface normal construction is outside the supported domain");
    if (!cross_y_checked) return Result<VectorField>::failure(cross_y_checked.error());

    auto cross_z = SymbolicExpr::add(
        SymbolicExpr::multiply(r_u[0], r_v[1]),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(r_u[1], r_v[0])));
    auto cross_z_checked = vector_calculus_simplify_strict(
        cross_z, operation, "surface normal construction is outside the supported domain");
    if (!cross_z_checked) return Result<VectorField>::failure(cross_z_checked.error());

    return Result<VectorField>::success(
        VectorField{std::move(cross_x_checked.value()),
                    std::move(cross_y_checked.value()),
                    std::move(cross_z_checked.value())});
}

static VectorCalculusExprResult curve_integral_scalar_strict(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::string& operation)
{
    const auto coord_vars = vector_calculus_coord_vars(parametrization.size());
    auto f_composed = vector_calculus_substitute_coords_strict(
        f, coord_vars, parametrization, operation);
    if (!f_composed) return f_composed;

    std::shared_ptr<SymbolicExpr> speed_sq;
    for (const auto& component : parametrization) {
        auto derivative = vector_calculus_differentiate_strict(
            component, t, operation);
        if (!derivative) return derivative;
        auto sq = SymbolicExpr::power(derivative.value(), SymbolicExpr::number(2));
        auto sq_checked = vector_calculus_simplify_strict(
            sq, operation, "curve speed construction is outside the supported domain");
        if (!sq_checked) return sq_checked;
        speed_sq = speed_sq ? SymbolicExpr::add(speed_sq, sq_checked.value())
                            : sq_checked.value();
        auto speed_sq_checked = vector_calculus_simplify_strict(
            speed_sq, operation, "curve speed construction is outside the supported domain");
        if (!speed_sq_checked) return speed_sq_checked;
        speed_sq = std::move(speed_sq_checked.value());
    }

    auto speed = SymbolicExpr::sqrt(speed_sq);
    auto speed_checked = vector_calculus_simplify_strict(
        speed, operation, "curve speed construction is outside the supported domain");
    if (!speed_checked) return speed_checked;

    auto integrand = SymbolicExpr::multiply(f_composed.value(), speed_checked.value());
    auto integrand_checked = vector_calculus_simplify_strict(
        integrand, operation, "curve integrand construction is outside the supported domain");
    if (!integrand_checked) return integrand_checked;

    return vector_calculus_definite_integral_strict(
        integrand_checked.value(), t, a, b, operation);
}

static VectorCalculusExprResult curve_integral_vector_strict(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::string& operation)
{
    const auto coord_vars = vector_calculus_coord_vars(parametrization.size());
    std::shared_ptr<SymbolicExpr> dot_product;
    for (size_t i = 0; i < F.size(); ++i) {
        auto Fi_composed = vector_calculus_substitute_coords_strict(
            F[i], coord_vars, parametrization, operation);
        if (!Fi_composed) return Fi_composed;

        auto ri_prime = vector_calculus_differentiate_strict(
            parametrization[i], t, operation);
        if (!ri_prime) return ri_prime;

        auto term = SymbolicExpr::multiply(Fi_composed.value(), ri_prime.value());
        auto term_checked = vector_calculus_simplify_strict(
            term, operation, "curve integrand construction is outside the supported domain");
        if (!term_checked) return term_checked;

        dot_product = dot_product ? SymbolicExpr::add(dot_product, term_checked.value())
                                  : term_checked.value();
        auto dot_checked = vector_calculus_simplify_strict(
            dot_product, operation, "curve integrand construction is outside the supported domain");
        if (!dot_checked) return dot_checked;
        dot_product = std::move(dot_checked.value());
    }

    return vector_calculus_definite_integral_strict(
        dot_product, t, a, b, operation);
}

static VectorCalculusExprResult surface_integral_scalar_strict(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper,
    const std::string& operation)
{
    const std::vector<std::string> coord_vars = {"x", "y", "z"};
    auto f_composed = vector_calculus_substitute_coords_strict(
        f, coord_vars, parametrization, operation);
    if (!f_composed) return f_composed;

    auto cross = vector_calculus_cross_product_partials_strict(
        parametrization, u, v, operation);
    if (!cross) return VectorCalculusExprResult::failure(cross.error());

    std::shared_ptr<SymbolicExpr> mag_sq;
    for (const auto& component : cross.value()) {
        auto sq = SymbolicExpr::power(component, SymbolicExpr::number(2));
        auto sq_checked = vector_calculus_simplify_strict(
            sq, operation, "surface magnitude construction is outside the supported domain");
        if (!sq_checked) return sq_checked;
        mag_sq = mag_sq ? SymbolicExpr::add(mag_sq, sq_checked.value())
                        : sq_checked.value();
        auto mag_sq_checked = vector_calculus_simplify_strict(
            mag_sq, operation, "surface magnitude construction is outside the supported domain");
        if (!mag_sq_checked) return mag_sq_checked;
        mag_sq = std::move(mag_sq_checked.value());
    }

    auto magnitude = SymbolicExpr::sqrt(mag_sq);
    auto magnitude_checked = vector_calculus_simplify_strict(
        magnitude, operation, "surface magnitude construction is outside the supported domain");
    if (!magnitude_checked) return magnitude_checked;

    auto integrand = SymbolicExpr::multiply(f_composed.value(), magnitude_checked.value());
    auto integrand_checked = vector_calculus_simplify_strict(
        integrand, operation, "surface integrand construction is outside the supported domain");
    if (!integrand_checked) return integrand_checked;

    std::vector<MultipleIntegralEngine::IntegrationStep> steps;
    steps.push_back({v, v_lower, v_upper});
    steps.push_back({u, u_lower, u_upper});
    return vector_calculus_multiple_integral_strict(
        integrand_checked.value(), steps, operation);
}

static VectorCalculusExprResult surface_integral_vector_strict(
    const VectorField& F, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper,
    const std::string& operation)
{
    const std::vector<std::string> coord_vars = {"x", "y", "z"};
    auto cross = vector_calculus_cross_product_partials_strict(
        parametrization, u, v, operation);
    if (!cross) return VectorCalculusExprResult::failure(cross.error());

    std::shared_ptr<SymbolicExpr> dot_product;
    for (size_t i = 0; i < 3; ++i) {
        auto Fi_composed = vector_calculus_substitute_coords_strict(
            F[i], coord_vars, parametrization, operation);
        if (!Fi_composed) return Fi_composed;

        auto term = SymbolicExpr::multiply(Fi_composed.value(), cross.value()[i]);
        auto term_checked = vector_calculus_simplify_strict(
            term, operation, "surface integrand construction is outside the supported domain");
        if (!term_checked) return term_checked;

        dot_product = dot_product ? SymbolicExpr::add(dot_product, term_checked.value())
                                  : term_checked.value();
        auto dot_checked = vector_calculus_simplify_strict(
            dot_product, operation, "surface integrand construction is outside the supported domain");
        if (!dot_checked) return dot_checked;
        dot_product = std::move(dot_checked.value());
    }

    std::vector<MultipleIntegralEngine::IntegrationStep> steps;
    steps.push_back({v, v_lower, v_upper});
    steps.push_back({u, u_lower, u_upper});
    return vector_calculus_multiple_integral_strict(dot_product, steps, operation);
}

static VectorCalculusExprResult greens_theorem_strict(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::string& operation)
{
    auto dQ_dx = vector_calculus_differentiate_strict(Q, vars[0], operation);
    if (!dQ_dx) return dQ_dx;
    auto dP_dy = vector_calculus_differentiate_strict(P, vars[1], operation);
    if (!dP_dy) return dP_dy;

    auto integrand = SymbolicExpr::add(
        dQ_dx.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dP_dy.value()));
    auto integrand_checked = vector_calculus_simplify_strict(
        integrand, operation, "Green's theorem integrand is outside the supported domain");
    if (!integrand_checked) return integrand_checked;

    std::vector<MultipleIntegralEngine::IntegrationStep> steps;
    steps.push_back({vars[1], y_bounds.first, y_bounds.second});
    steps.push_back({vars[0], x_bounds.first, x_bounds.second});
    return vector_calculus_multiple_integral_strict(
        integrand_checked.value(), steps, operation);
}

static VectorCalculusExprResult greens_theorem_area_strict(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    const std::string& operation)
{
    auto dx_dt = vector_calculus_differentiate_strict(
        parametrization[0], t, operation);
    if (!dx_dt) return dx_dt;
    auto dy_dt = vector_calculus_differentiate_strict(
        parametrization[1], t, operation);
    if (!dy_dt) return dy_dt;

    auto term1 = SymbolicExpr::multiply(parametrization[0], dy_dt.value());
    auto term2 = SymbolicExpr::multiply(parametrization[1], dx_dt.value());
    auto integrand = SymbolicExpr::add(
        term1,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), term2));
    auto integrand_checked = vector_calculus_simplify_strict(
        integrand, operation, "Green's area integrand is outside the supported domain");
    if (!integrand_checked) return integrand_checked;

    auto integral = vector_calculus_definite_integral_strict(
        integrand_checked.value(), t, a, b, operation);
    if (!integral) return integral;

    auto half = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(Rational(1, 2)));
    auto area = SymbolicExpr::multiply(half, integral.value());
    return vector_calculus_simplify_strict(
        area, operation, "Green's area result is outside the supported domain");
}

static VectorCalculusExprResult divergence_theorem_strict(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& z_bounds,
    const std::string& operation)
{
    std::shared_ptr<SymbolicExpr> div_F;
    for (size_t i = 0; i < 3; ++i) {
        auto partial = vector_calculus_differentiate_strict(F[i], vars[i], operation);
        if (!partial) return partial;
        div_F = div_F ? SymbolicExpr::add(div_F, partial.value()) : partial.value();
        auto div_checked = vector_calculus_simplify_strict(
            div_F, operation, "divergence theorem integrand is outside the supported domain");
        if (!div_checked) return div_checked;
        div_F = std::move(div_checked.value());
    }

    std::vector<MultipleIntegralEngine::IntegrationStep> steps;
    steps.push_back({vars[2], z_bounds.first, z_bounds.second});
    steps.push_back({vars[1], y_bounds.first, y_bounds.second});
    steps.push_back({vars[0], x_bounds.first, x_bounds.second});
    return vector_calculus_multiple_integral_strict(div_F, steps, operation);
}

static Result<VectorField> vector_calculus_curl_strict(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::string& operation)
{
    auto dFz_dy = vector_calculus_differentiate_strict(F[2], vars[1], operation);
    if (!dFz_dy) return Result<VectorField>::failure(dFz_dy.error());
    auto dFy_dz = vector_calculus_differentiate_strict(F[1], vars[2], operation);
    if (!dFy_dz) return Result<VectorField>::failure(dFy_dz.error());
    auto dFx_dz = vector_calculus_differentiate_strict(F[0], vars[2], operation);
    if (!dFx_dz) return Result<VectorField>::failure(dFx_dz.error());
    auto dFz_dx = vector_calculus_differentiate_strict(F[2], vars[0], operation);
    if (!dFz_dx) return Result<VectorField>::failure(dFz_dx.error());
    auto dFy_dx = vector_calculus_differentiate_strict(F[1], vars[0], operation);
    if (!dFy_dx) return Result<VectorField>::failure(dFy_dx.error());
    auto dFx_dy = vector_calculus_differentiate_strict(F[0], vars[1], operation);
    if (!dFx_dy) return Result<VectorField>::failure(dFx_dy.error());

    auto cx = SymbolicExpr::add(
        dFz_dy.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dFy_dz.value()));
    auto cx_checked = vector_calculus_simplify_strict(
        cx, operation, "curl construction is outside the supported domain");
    if (!cx_checked) return Result<VectorField>::failure(cx_checked.error());

    auto cy = SymbolicExpr::add(
        dFx_dz.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dFz_dx.value()));
    auto cy_checked = vector_calculus_simplify_strict(
        cy, operation, "curl construction is outside the supported domain");
    if (!cy_checked) return Result<VectorField>::failure(cy_checked.error());

    auto cz = SymbolicExpr::add(
        dFy_dx.value(),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dFx_dy.value()));
    auto cz_checked = vector_calculus_simplify_strict(
        cz, operation, "curl construction is outside the supported domain");
    if (!cz_checked) return Result<VectorField>::failure(cz_checked.error());

    return Result<VectorField>::success(
        VectorField{std::move(cx_checked.value()),
                    std::move(cy_checked.value()),
                    std::move(cz_checked.value())});
}

static VectorCalculusExprResult stokes_theorem_strict(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const VectorField& parametrization,
    const std::string& u,
    const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& v_bounds,
    const std::string& operation)
{
    auto curl_F = vector_calculus_curl_strict(F, vars, operation);
    if (!curl_F) return VectorCalculusExprResult::failure(curl_F.error());

    auto cross = vector_calculus_cross_product_partials_strict(
        parametrization, u, v, operation);
    if (!cross) return VectorCalculusExprResult::failure(cross.error());

    const std::vector<std::string> coord_vars = {"x", "y", "z"};
    std::shared_ptr<SymbolicExpr> dot_product;
    for (size_t i = 0; i < 3; ++i) {
        auto curl_i_composed = vector_calculus_substitute_coords_strict(
            curl_F.value()[i], coord_vars, parametrization, operation);
        if (!curl_i_composed) return curl_i_composed;

        auto term = SymbolicExpr::multiply(curl_i_composed.value(), cross.value()[i]);
        auto term_checked = vector_calculus_simplify_strict(
            term, operation, "Stokes theorem integrand is outside the supported domain");
        if (!term_checked) return term_checked;

        dot_product = dot_product ? SymbolicExpr::add(dot_product, term_checked.value())
                                  : term_checked.value();
        auto dot_checked = vector_calculus_simplify_strict(
            dot_product, operation, "Stokes theorem integrand is outside the supported domain");
        if (!dot_checked) return dot_checked;
        dot_product = std::move(dot_checked.value());
    }

    std::vector<MultipleIntegralEngine::IntegrationStep> steps;
    steps.push_back({v, v_bounds.first, v_bounds.second});
    steps.push_back({u, u_bounds.first, u_bounds.second});
    return vector_calculus_multiple_integral_strict(dot_product, steps, operation);
}


VectorCalculusExprResult curve_integral_scalar_checked(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context)
{
    const std::string operation = "curve_integral_scalar";
    auto valid = vector_calculus_validate_curve_scalar_inputs(
        f, parametrization, t, a, b, context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());

    auto step = context.consume_steps(parametrization.size() * 4 + 4, operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return curve_integral_scalar_strict(
            f, parametrization, t, a, b, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "curve integral allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult curve_integral_scalar_checked(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    ComputationContext context;
    return curve_integral_scalar_checked(f, parametrization, t, a, b, context);
}

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


VectorCalculusExprResult curve_integral_vector_checked(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context)
{
    const std::string operation = "curve_integral_vector";
    auto valid = vector_calculus_validate_curve_vector_inputs(
        F, parametrization, t, a, b, context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());

    auto step = context.consume_steps(F.size() * parametrization.size() +
                                      parametrization.size() * 3 + 4,
                                      operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return curve_integral_vector_strict(
            F, parametrization, t, a, b, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "curve integral allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult curve_integral_vector_checked(
    const VectorField& F, const VectorField& parametrization,
    const std::string& t, const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    ComputationContext context;
    return curve_integral_vector_checked(F, parametrization, t, a, b, context);
}

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


VectorCalculusExprResult surface_integral_scalar_checked(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context)
{
    const std::string operation = "surface_integral_scalar";
    auto valid = vector_calculus_validate_surface_scalar_inputs(
        f, parametrization, u, v, u_lower, u_upper, v_lower, v_upper,
        context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());

    auto step = context.consume_steps(parametrization.size() * 6 + 8, operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return surface_integral_scalar_strict(
            f, parametrization, u, v, u_lower, u_upper, v_lower, v_upper,
            operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "surface integral allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult surface_integral_scalar_checked(
    const std::shared_ptr<SymbolicExpr>& f, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper)
{
    ComputationContext context;
    return surface_integral_scalar_checked(
        f, parametrization, u, v, u_lower, u_upper, v_lower, v_upper, context);
}

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


VectorCalculusExprResult surface_integral_vector_checked(
    const VectorField& F, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper,
    ComputationContext& context)
{
    const std::string operation = "surface_integral_vector";
    auto valid = vector_calculus_validate_surface_vector_inputs(
        F, parametrization, u, v, u_lower, u_upper, v_lower, v_upper,
        context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());

    auto step = context.consume_steps(F.size() * parametrization.size() +
                                      parametrization.size() * 6 + 8,
                                      operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return surface_integral_vector_strict(
            F, parametrization, u, v, u_lower, u_upper, v_lower, v_upper,
            operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "surface integral allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult surface_integral_vector_checked(
    const VectorField& F, const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::shared_ptr<SymbolicExpr>& u_lower, const std::shared_ptr<SymbolicExpr>& u_upper,
    const std::shared_ptr<SymbolicExpr>& v_lower, const std::shared_ptr<SymbolicExpr>& v_upper)
{
    ComputationContext context;
    return surface_integral_vector_checked(
        F, parametrization, u, v, u_lower, u_upper, v_lower, v_upper, context);
}

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


VectorCalculusExprResult greens_theorem_checked(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    ComputationContext& context)
{
    const std::string operation = "greens_theorem";
    auto vars_valid = vector_calculus_validate_distinct_vars(vars, 2, context, operation);
    if (!vars_valid) return VectorCalculusExprResult::failure(vars_valid.error());
    if (!P || !lamina::detail::node(P) || !Q || !lamina::detail::node(Q)) {
        return VectorCalculusExprResult::failure(
            CasErrc::InvalidArgument,
            "Green's theorem vector-field components cannot be null",
            operation);
    }
    auto x_valid = vector_calculus_validate_bound_pair(x_bounds, context, operation, "x");
    if (!x_valid) return VectorCalculusExprResult::failure(x_valid.error());
    auto y_valid = vector_calculus_validate_bound_pair(y_bounds, context, operation, "y");
    if (!y_valid) return VectorCalculusExprResult::failure(y_valid.error());
    auto step = context.consume_steps(10, operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return greens_theorem_strict(P, Q, vars, x_bounds, y_bounds, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "Green's theorem allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult greens_theorem_checked(
    const std::shared_ptr<SymbolicExpr>& P,
    const std::shared_ptr<SymbolicExpr>& Q,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds)
{
    ComputationContext context;
    return greens_theorem_checked(P, Q, vars, x_bounds, y_bounds, context);
}

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


VectorCalculusExprResult greens_theorem_area_checked(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context)
{
    const std::string operation = "greens_theorem_area";
    auto valid = vector_calculus_validate_curve_parametrization(
        parametrization, t, a, b, context, operation);
    if (!valid) return VectorCalculusExprResult::failure(valid.error());
    if (parametrization.size() != 2) {
        return VectorCalculusExprResult::failure(
            CasErrc::InvalidArgument,
            "Green's area parametrization must be two-dimensional",
            operation);
    }
    auto step = context.consume_steps(10, operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return greens_theorem_area_strict(parametrization, t, a, b, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "Green's area allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult greens_theorem_area_checked(
    const VectorField& parametrization,
    const std::string& t,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    ComputationContext context;
    return greens_theorem_area_checked(parametrization, t, a, b, context);
}

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
    auto half = lamina::detail::make_expression_ptr(
        lamina::detail::make_node<NumberNode>(Rational(1, 2)));
    auto area = SymbolicExpr::multiply(half, integral_result);
    area = area->simplify();

    return area;
}


VectorCalculusExprResult divergence_theorem_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& z_bounds,
    ComputationContext& context)
{
    const std::string operation = "divergence_theorem";
    auto field_valid = vector_calculus_validate_field_vars(F, vars, context, operation);
    if (!field_valid) return VectorCalculusExprResult::failure(field_valid.error());
    auto vars_valid = vector_calculus_validate_distinct_vars(vars, 3, context, operation);
    if (!vars_valid) return VectorCalculusExprResult::failure(vars_valid.error());
    auto x_valid = vector_calculus_validate_bound_pair(x_bounds, context, operation, "x");
    if (!x_valid) return VectorCalculusExprResult::failure(x_valid.error());
    auto y_valid = vector_calculus_validate_bound_pair(y_bounds, context, operation, "y");
    if (!y_valid) return VectorCalculusExprResult::failure(y_valid.error());
    auto z_valid = vector_calculus_validate_bound_pair(z_bounds, context, operation, "z");
    if (!z_valid) return VectorCalculusExprResult::failure(z_valid.error());
    auto step = context.consume_steps(14, operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return divergence_theorem_strict(
            F, vars, x_bounds, y_bounds, z_bounds, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "divergence theorem allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult divergence_theorem_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& x_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& y_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& z_bounds)
{
    ComputationContext context;
    return divergence_theorem_checked(F, vars, x_bounds, y_bounds, z_bounds, context);
}

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


VectorCalculusExprResult stokes_theorem_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& v_bounds,
    ComputationContext& context)
{
    const std::string operation = "stokes_theorem";
    auto field_valid = vector_calculus_validate_field_vars(F, vars, context, operation);
    if (!field_valid) return VectorCalculusExprResult::failure(field_valid.error());
    auto vars_valid = vector_calculus_validate_distinct_vars(vars, 3, context, operation);
    if (!vars_valid) return VectorCalculusExprResult::failure(vars_valid.error());
    auto param_valid = vector_calculus_validate_surface_parametrization(
        parametrization, u, v, u_bounds.first, u_bounds.second,
        v_bounds.first, v_bounds.second, context, operation);
    if (!param_valid) return VectorCalculusExprResult::failure(param_valid.error());
    auto step = context.consume_steps(18, operation);
    if (!step) return VectorCalculusExprResult::failure(step.error());

    try {
        return stokes_theorem_strict(
            F, vars, parametrization, u, v, u_bounds, v_bounds, operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "Stokes theorem allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

VectorCalculusExprResult stokes_theorem_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    const VectorField& parametrization,
    const std::string& u, const std::string& v,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& u_bounds,
    const std::pair<std::shared_ptr<SymbolicExpr>, std::shared_ptr<SymbolicExpr>>& v_bounds)
{
    ComputationContext context;
    return stokes_theorem_checked(F, vars, parametrization, u, v,
                                  u_bounds, v_bounds, context);
}

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
    const std::vector<std::string>&,
    const std::map<std::string, std::shared_ptr<SymbolicExpr>>& pt,
    size_t n,
    std::vector<double>& numeric_H)
{
    auto mat_node = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(H));
    if (!mat_node || mat_node->rows() != n || mat_node->cols() != n) {
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
            auto elem = lamina::detail::make_expression_ptr(elem_node);
            /// 代入临界点坐标
            for (const auto& [var_name, val] : pt) {
                elem = elem->substitute(var_name, val);
                if (!elem) return false;
            }
            elem = elem->simplify();
            if (!elem) return false;

            if (!vector_calculus_checked_finite_numeric(
                    elem, numeric_H[i * n + j])) {
                return false;
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
        if (!vector_calculus_checked_finite_numeric(ev_simplified, val)) {
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

ExtremaResult find_extrema_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    const std::string operation = "find_extrema";
    auto valid = vector_calculus_validate_expr_vars(f, vars, context, operation);
    if (!valid) return ExtremaResult::failure(valid.error());
    auto distinct = vector_calculus_validate_distinct_vars(vars, vars.size(), context, operation);
    if (!distinct) return ExtremaResult::failure(distinct.error());
    auto step = context.consume_steps(vars.size() * vars.size() + vars.size() * 8 + 8,
                                      operation);
    if (!step) return ExtremaResult::failure(step.error());

    try {
        std::vector<std::shared_ptr<SymbolicExpr>> gradient;
        gradient.reserve(vars.size());
        for (const auto& var : vars) {
            auto partial = vector_calculus_differentiate_strict(f, var, operation);
            if (!partial) return ExtremaResult::failure(partial.error());
            gradient.push_back(std::move(partial.value()));
        }

        auto extrema = find_extrema(f, vars);
        if (extrema.empty()) {
            return ExtremaResult::failure(
                CasErrc::Inconclusive,
                "extrema solver produced no verifiable critical points in the supported domain",
                operation);
        }

        for (const auto& cp : extrema) {
            if (!vector_calculus_point_has_vars(cp.point, vars)) {
                return ExtremaResult::failure(
                    CasErrc::Inconclusive,
                    "extrema candidate omits one or more variables",
                    operation);
            }
            for (const auto& partial : gradient) {
                if (!vector_calculus_expr_zero_after_substitution(partial, cp.point)) {
                    return ExtremaResult::failure(
                        CasErrc::Inconclusive,
                        "extrema candidate does not verify against the gradient",
                        operation);
                }
            }
            if (cp.classification.empty()) {
                return ExtremaResult::failure(
                    CasErrc::Inconclusive,
                    "extrema candidate has no verified classification",
                    operation);
            }
        }
        return ExtremaResult::success(std::move(extrema));
    } catch (const std::bad_alloc&) {
        return ExtremaResult::failure(CasErrc::ResourceLimit,
                                      "extrema allocation failed",
                                      operation);
    } catch (const std::exception& e) {
        return ExtremaResult::failure(CasErrc::InternalInvariant,
                                      e.what(), operation);
    }
}

ExtremaResult find_extrema_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return find_extrema_checked(f, vars, context);
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
                pt[name] = lamina::detail::make_expression_ptr(val);
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


static LagrangeResult lagrange_multipliers_strict(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars,
    const std::string& operation)
{
    const size_t n = vars.size();
    const size_t m = constraints.size();

    VectorField grad_f;
    grad_f.reserve(n);
    for (const auto& var : vars) {
        auto partial = vector_calculus_differentiate_strict(f, var, operation);
        if (!partial) return LagrangeResult::failure(partial.error());
        grad_f.push_back(std::move(partial.value()));
    }

    std::vector<VectorField> grad_constraints;
    grad_constraints.reserve(m);
    for (const auto& constraint : constraints) {
        VectorField grad_g;
        grad_g.reserve(n);
        for (const auto& var : vars) {
            auto partial = vector_calculus_differentiate_strict(
                constraint, var, operation);
            if (!partial) return LagrangeResult::failure(partial.error());
            grad_g.push_back(std::move(partial.value()));
        }
        grad_constraints.push_back(std::move(grad_g));
    }

    std::vector<std::string> lambda_names;
    lambda_names.reserve(m);
    for (size_t k = 0; k < m; ++k) {
        lambda_names.push_back("lambda_" + std::to_string(k + 1));
    }

    std::vector<std::shared_ptr<SymbolicExpr>> equations;
    equations.reserve(n + m);
    for (size_t i = 0; i < n; ++i) {
        auto eq = grad_f[i];
        for (size_t k = 0; k < m; ++k) {
            auto lambda_var = SymbolicExpr::variable(lambda_names[k]);
            auto term = SymbolicExpr::multiply(lambda_var, grad_constraints[k][i]);
            eq = SymbolicExpr::add(
                eq, SymbolicExpr::multiply(SymbolicExpr::number(-1), term));
        }
        auto eq_checked = vector_calculus_simplify_strict(
            eq, operation, "Lagrange stationarity equation is outside the supported domain");
        if (!eq_checked) return LagrangeResult::failure(eq_checked.error());
        equations.push_back(std::move(eq_checked.value()));
    }
    for (const auto& constraint : constraints) {
        equations.push_back(constraint);
    }

    std::vector<std::string> all_vars;
    all_vars.reserve(n + m);
    all_vars.insert(all_vars.end(), vars.begin(), vars.end());
    all_vars.insert(all_vars.end(), lambda_names.begin(), lambda_names.end());

    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> full_solutions;
    std::vector<SymbolicExpr> poly_eqs;
    poly_eqs.reserve(equations.size());
    for (const auto& eq : equations) {
        if (!eq || !lamina::detail::node(eq)) {
            return LagrangeResult::failure(
                CasErrc::Inconclusive,
                "Lagrange equation construction failed in the supported domain",
                operation);
        }
        poly_eqs.push_back(*eq);
    }

    auto poly_solutions = Solver::solve_polynomial_system(poly_eqs, all_vars);
    for (const auto& sol : poly_solutions) {
        std::map<std::string, std::shared_ptr<SymbolicExpr>> pt;
        for (const auto& [name, val] : sol) {
            pt[name] = lamina::detail::make_expression_ptr(val);
        }
        full_solutions.push_back(std::move(pt));
    }

    if (full_solutions.empty()) {
        full_solutions = SymbolicExpr::solve_system(equations, all_vars);
    }
    if (full_solutions.empty()) {
        return LagrangeResult::failure(
            CasErrc::Inconclusive,
            "Lagrange solver produced no verifiable candidates in the supported domain",
            operation);
    }

    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> verified;
    std::set<std::string> var_set(vars.begin(), vars.end());
    for (const auto& full_solution : full_solutions) {
        if (!vector_calculus_point_has_vars(full_solution, all_vars)) {
            return LagrangeResult::failure(
                CasErrc::Inconclusive,
                "Lagrange candidate omits one or more variables or multipliers",
                operation);
        }
        for (const auto& eq : equations) {
            if (!vector_calculus_expr_zero_after_substitution(eq, full_solution)) {
                return LagrangeResult::failure(
                    CasErrc::Inconclusive,
                    "Lagrange candidate does not verify against the full stationarity system",
                    operation);
            }
        }

        std::map<std::string, std::shared_ptr<SymbolicExpr>> filtered;
        for (const auto& [name, val] : full_solution) {
            if (var_set.count(name)) {
                filtered[name] = val;
            }
        }
        if (!vector_calculus_point_has_vars(filtered, vars)) {
            return LagrangeResult::failure(
                CasErrc::Inconclusive,
                "Lagrange candidate omits one or more variables",
                operation);
        }
        verified.push_back(std::move(filtered));
    }

    if (verified.empty()) {
        return LagrangeResult::failure(
            CasErrc::Inconclusive,
            "Lagrange solver produced no verifiable candidates in the supported domain",
            operation);
    }
    return LagrangeResult::success(std::move(verified));
}

LagrangeResult lagrange_multipliers_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    const std::string operation = "lagrange_multipliers";
    auto valid = vector_calculus_validate_expr_vars(f, vars, context, operation);
    if (!valid) return LagrangeResult::failure(valid.error());
    auto distinct = vector_calculus_validate_distinct_vars(vars, vars.size(), context, operation);
    if (!distinct) return LagrangeResult::failure(distinct.error());
    if (constraints.empty()) {
        return LagrangeResult::failure(CasErrc::InvalidArgument,
                                       "constraint list cannot be empty",
                                       operation);
    }
    for (const auto& constraint : constraints) {
        if (!constraint || !lamina::detail::node(constraint)) {
            return LagrangeResult::failure(CasErrc::InvalidArgument,
                                           "constraint expressions cannot be null",
                                           operation);
        }
    }
    auto step = context.consume_steps(vars.size() * constraints.size() * 6 +
                                      vars.size() * 6 + constraints.size() * 4 + 8,
                                      operation);
    if (!step) return LagrangeResult::failure(step.error());

    try {
        return lagrange_multipliers_strict(f, constraints, vars, operation);
    } catch (const std::bad_alloc&) {
        return LagrangeResult::failure(CasErrc::ResourceLimit,
                                       "Lagrange multiplier allocation failed",
                                       operation);
    } catch (const std::exception& e) {
        return LagrangeResult::failure(CasErrc::InternalInvariant,
                                       e.what(), operation);
    }
}

LagrangeResult lagrange_multipliers_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return lagrange_multipliers_checked(f, constraints, vars, context);
}

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
                pt[name] = lamina::detail::make_expression_ptr(val);
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
    if (lamina::detail::node(bb) && lamina::detail::node(bb)->is_zero()) {
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
    if (lamina::detail::node(bb) && lamina::detail::node(bb)->is_zero()) {
        return nullptr;
    }
    auto ab = dot_product(a, b);
    return SymbolicExpr::divide(ab, SymbolicExpr::sqrt(bb))->simplify();
}

std::shared_ptr<SymbolicExpr> vector_angle_symbolic(const VectorField& a, const VectorField& b)
{
    auto aa = dot_product(a, a);
    auto bb = dot_product(b, b);
    if ((lamina::detail::node(aa) && lamina::detail::node(aa)->is_zero()) || (lamina::detail::node(bb) && lamina::detail::node(bb)->is_zero())) {
        return nullptr;
    }
    auto ab = dot_product(a, b);
    auto denom = SymbolicExpr::multiply(SymbolicExpr::sqrt(aa), SymbolicExpr::sqrt(bb));
    auto cos_theta = SymbolicExpr::divide(ab, denom);
    auto arccos_node = lamina::detail::make_node<FunctionNode>(
        FunctionNode::FuncType::ArcCos,
        std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(cos_theta)});
    return lamina::detail::make_expression_ptr(arccos_node)->simplify();
}

std::shared_ptr<SymbolicExpr> mixed_product(const VectorField& a, const VectorField& b,
    const VectorField& c)
{
    return dot_product(a, cross_product(b, c));
}

} // namespace lamina
