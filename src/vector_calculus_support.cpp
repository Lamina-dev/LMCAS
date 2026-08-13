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

namespace lamina::vector_calculus_detail {

bool vector_calculus_checked_finite_numeric(
    const std::shared_ptr<SymbolicExpr>& expr,
    double& value,
    ComputationContext* context)
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
    std::size_t depth)
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
    if (auto product = std::dynamic_pointer_cast<const ProductNode>(node)) {
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

} // namespace lamina::vector_calculus_detail
