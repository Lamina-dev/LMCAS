#include "internal/expression_transform.hpp"
#include "internal/vector_calculus_support.hpp"
#include "vector_calculus_differential.hpp"
#include <stdexcept>

namespace LMCAS {

using namespace vector_calculus_detail;

VectorCalculusFieldResult gradient_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "gradient");
    const std::string operation = "gradient";
    auto valid = vector_calculus_validate_expr_vars(f, vars, context, operation);
    if (!valid) {
        return VectorCalculusFieldResult::failure(valid.error());
    }
    auto step = context.consume_steps(vars.size(), operation);
    if (!step) {
        return VectorCalculusFieldResult::failure(step.error());
    }
    try {
        return vector_calculus_wrap_field(gradient_impl(f, vars, transform_budget), operation);
    } catch (const CasError& error) {
        return VectorCalculusFieldResult::failure(error);
    } catch (const std::bad_alloc&) {
        return VectorCalculusFieldResult::failure(CasErrc::ResourceLimit,
                                                  "vector calculus allocation failed",
                                                  operation);
    } catch (const std::length_error& error) {
        return VectorCalculusFieldResult::failure(CasErrc::ResourceLimit, error.what(), operation);
    } catch (const std::exception& e) {
        return VectorCalculusFieldResult::failure(CasErrc::InternalInvariant, e.what(), operation);
    }
}

VectorCalculusFieldResult gradient_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return gradient_checked(f, vars, context);
}

VectorField vector_calculus_detail::gradient_impl(const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    detail::RewriteBudget& transform_budget)
{
    if (!f) {
        throw std::invalid_argument("gradient: f must not be null");
    }
    VectorField result;
    transform_budget.require_nodes(vars.size());
    std::size_t nodes = 0;
    result.reserve(vars.size());
    for (const auto& var : vars) {
        auto partial = detail::make_expression_ptr(detail::differentiate_node(detail::node(f), var, transform_budget));
        if (partial) {
            partial = detail::make_expression_ptr(detail::normalize_node(detail::node(partial), transform_budget));
        }
        nodes = transform_budget.append_size(
            nodes, transform_budget.measure(detail::node(partial)));
        result.push_back(partial);
    }
    return result;
}


VectorCalculusExprResult divergence_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "divergence");
    const std::string operation = "divergence";
    auto valid = vector_calculus_validate_field_vars(F, vars, context, operation);
    if (!valid) {
        return VectorCalculusExprResult::failure(valid.error());
    }
    auto step = context.consume_steps(F.size(), operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }
    try {
        return vector_calculus_wrap_expr(divergence_impl(F, vars, transform_budget), operation);
    } catch (const CasError& error) {
        return VectorCalculusExprResult::failure(error);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "vector calculus allocation failed",
                                                 operation);
    } catch (const std::length_error& error) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit, error.what(), operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant, e.what(), operation);
    }
}

VectorCalculusExprResult divergence_checked(
    const VectorField& F,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return divergence_checked(F, vars, context);
}

std::shared_ptr<SymbolicExpr> vector_calculus_detail::divergence_impl(const VectorField& F,
    const std::vector<std::string>& vars,
    detail::RewriteBudget& transform_budget)
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
        if (!F[i]) {
            continue;
        }
        auto partial = detail::make_expression_ptr(detail::differentiate_node(detail::node(F[i]), vars[i], transform_budget));
        if (!partial) {
            continue;
        }
        partial = detail::make_expression_ptr(detail::normalize_node(detail::node(partial), transform_budget));
        if (!sum) {
            sum = partial;
        } else {
            sum = SymbolicExpr::add(sum, partial);
            sum = detail::make_expression_ptr(detail::normalize_node(detail::node(sum), transform_budget));
        }
    }

    if (!sum) {
        return SymbolicExpr::number(0);
    }
    return detail::make_expression_ptr(detail::normalize_node(detail::node(sum), transform_budget));
}


VectorCalculusFieldResult curl_checked(
    const VectorField& F,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "curl");
    const std::string operation = "curl";
    auto valid = vector_calculus_validate_field_vars(F, vars, context, operation);
    if (!valid) {
        return VectorCalculusFieldResult::failure(valid.error());
    }
    if (F.size() != 2 && F.size() != 3) {
        return VectorCalculusFieldResult::failure(
            CasErrc::InvalidArgument,
            "curl requires a 2D or 3D vector field",
            operation);
    }
    auto step = context.consume_steps(F.size(), operation);
    if (!step) {
        return VectorCalculusFieldResult::failure(step.error());
    }
    try {
        return vector_calculus_wrap_field(curl_impl(F, vars, transform_budget), operation);
    } catch (const CasError& error) {
        return VectorCalculusFieldResult::failure(error);
    } catch (const std::bad_alloc&) {
        return VectorCalculusFieldResult::failure(CasErrc::ResourceLimit,
                                                  "vector calculus allocation failed",
                                                  operation);
    } catch (const std::length_error& error) {
        return VectorCalculusFieldResult::failure(CasErrc::ResourceLimit, error.what(), operation);
    } catch (const std::exception& e) {
        return VectorCalculusFieldResult::failure(CasErrc::InternalInvariant, e.what(), operation);
    }
}

VectorCalculusFieldResult curl_checked(
    const VectorField& F,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return curl_checked(F, vars, context);
}

VectorField vector_calculus_detail::curl_impl(const VectorField& F,
    const std::vector<std::string>& vars,
    detail::RewriteBudget& transform_budget)
{
    if (F.size() != vars.size()) {
        throw std::invalid_argument(
            "curl: F and vars must have the same dimension");
    }

    /// 二维标量旋度: ∂F₂/∂x₁ - ∂F₁/∂x₂
    if (F.size() == 2 && vars.size() == 2) {
        auto dF2_dx1 = detail::make_expression_ptr(detail::differentiate_node(detail::node(F[1]), vars[0], transform_budget));
        auto dF1_dx2 = detail::make_expression_ptr(detail::differentiate_node(detail::node(F[0]), vars[1], transform_budget));

        auto neg_dF1_dx2 = SymbolicExpr::multiply(
            SymbolicExpr::number(-1), dF1_dx2);
        auto scalar_curl = SymbolicExpr::add(dF2_dx1, neg_dF1_dx2);
        scalar_curl = detail::make_expression_ptr(detail::normalize_node(detail::node(scalar_curl), transform_budget));

        return VectorField{scalar_curl};
    }

    /// 三维旋度
    if (F.size() != 3 || vars.size() != 3) {
        throw std::invalid_argument(
            "curl: requires 2D or 3D vector field");
    }

    /// curl_x = ∂F₃/∂x₂ - ∂F₂/∂x₃
    auto dF3_dx2 = detail::make_expression_ptr(detail::differentiate_node(detail::node(F[2]), vars[1], transform_budget));
    auto dF2_dx3 = detail::make_expression_ptr(detail::differentiate_node(detail::node(F[1]), vars[2], transform_budget));
    auto curl_x = SymbolicExpr::add(
        dF3_dx2,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dF2_dx3));
    curl_x = detail::make_expression_ptr(detail::normalize_node(detail::node(curl_x), transform_budget));

    /// curl_y = ∂F₁/∂x₃ - ∂F₃/∂x₁
    auto dF1_dx3 = detail::make_expression_ptr(detail::differentiate_node(detail::node(F[0]), vars[2], transform_budget));
    auto dF3_dx1 = detail::make_expression_ptr(detail::differentiate_node(detail::node(F[2]), vars[0], transform_budget));
    auto curl_y = SymbolicExpr::add(
        dF1_dx3,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dF3_dx1));
    curl_y = detail::make_expression_ptr(detail::normalize_node(detail::node(curl_y), transform_budget));

    /// curl_z = ∂F₂/∂x₁ - ∂F₁/∂x₂
    auto dF2_dx1 = detail::make_expression_ptr(detail::differentiate_node(detail::node(F[1]), vars[0], transform_budget));
    auto dF1_dx2 = detail::make_expression_ptr(detail::differentiate_node(detail::node(F[0]), vars[1], transform_budget));
    auto curl_z = SymbolicExpr::add(
        dF2_dx1,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), dF1_dx2));
    curl_z = detail::make_expression_ptr(detail::normalize_node(detail::node(curl_z), transform_budget));

    VectorField result{curl_x, curl_y, curl_z};
    std::size_t nodes = 0;
    for (const auto& component : result) {
        nodes = transform_budget.append_size(
            nodes, transform_budget.measure(detail::node(component)));
    }
    return result;
}


VectorCalculusExprResult laplacian_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "laplacian");
    const std::string operation = "laplacian";
    auto valid = vector_calculus_validate_expr_vars(f, vars, context, operation);
    if (!valid) {
        return VectorCalculusExprResult::failure(valid.error());
    }
    auto step = context.consume_steps(vars.size() * 2, operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }
    try {
        return vector_calculus_wrap_expr(laplacian_impl(f, vars, transform_budget), operation);
    } catch (const CasError& error) {
        return VectorCalculusExprResult::failure(error);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "vector calculus allocation failed",
                                                 operation);
    } catch (const std::length_error& error) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit, error.what(), operation);
    } catch (const std::exception& e) {
        return VectorCalculusExprResult::failure(CasErrc::InternalInvariant, e.what(), operation);
    }
}

VectorCalculusExprResult laplacian_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return laplacian_checked(f, vars, context);
}

std::shared_ptr<SymbolicExpr> vector_calculus_detail::laplacian_impl(const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    detail::RewriteBudget& transform_budget)
{
    if (!f) {
        throw std::invalid_argument("laplacian: f must not be null");
    }
    if (vars.empty()) {
        return SymbolicExpr::number(0);
    }

    std::shared_ptr<SymbolicExpr> sum = nullptr;
    for (const auto& var : vars) {
        auto first = detail::make_expression_ptr(detail::differentiate_node(detail::node(f), var, transform_budget));
        if (!first) {
            continue;
        }
        auto second = detail::make_expression_ptr(detail::differentiate_node(detail::node(first), var, transform_budget));
        if (!second) {
            continue;
        }
        second = detail::make_expression_ptr(detail::normalize_node(detail::node(second), transform_budget));
        if (!sum) {
            sum = second;
        } else {
            sum = SymbolicExpr::add(sum, second);
            sum = detail::make_expression_ptr(detail::normalize_node(detail::node(sum), transform_budget));
        }
    }

    if (!sum) {
        return SymbolicExpr::number(0);
    }
    return detail::make_expression_ptr(detail::normalize_node(detail::node(sum), transform_budget));
}


VectorField gradient(const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "gradient");
    return gradient_impl(f, vars, transform_budget);
}

std::shared_ptr<SymbolicExpr> divergence(const VectorField& F,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "divergence");
    return divergence_impl(F, vars, transform_budget);
}

VectorField curl(const VectorField& F,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "curl");
    return curl_impl(F, vars, transform_budget);
}

std::shared_ptr<SymbolicExpr> laplacian(const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                          context.limits().max_ast_nodes, "laplacian");
    return laplacian_impl(f, vars, transform_budget);
}

}
