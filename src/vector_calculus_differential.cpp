#include "internal/vector_calculus_support.hpp"
#include "vector_calculus_differential.hpp"
#include "internal/visitors/differentiation_visitor.hpp"
#include <stdexcept>

namespace LMCAS {

using namespace vector_calculus_detail;

VectorCalculusFieldResult gradient_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
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
        return vector_calculus_wrap_field(gradient(f, vars), operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusFieldResult::failure(CasErrc::ResourceLimit,
                                                  "vector calculus allocation failed",
                                                  operation);
    } catch (const detail::UnsupportedDifferentiation& error) {
        return VectorCalculusFieldResult::failure(
            CasErrc::UnsupportedExpression, error.what(), operation);
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
    if (!valid) {
        return VectorCalculusExprResult::failure(valid.error());
    }
    auto step = context.consume_steps(F.size(), operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }
    try {
        return vector_calculus_wrap_expr(divergence(F, vars), operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "vector calculus allocation failed",
                                                 operation);
    } catch (const detail::UnsupportedDifferentiation& error) {
        return VectorCalculusExprResult::failure(
            CasErrc::UnsupportedExpression, error.what(), operation);
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
        if (!F[i]) {
            continue;
        }
        auto partial = F[i]->differentiate(vars[i]);
        if (!partial) {
            continue;
        }
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
        return vector_calculus_wrap_field(curl(F, vars), operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusFieldResult::failure(CasErrc::ResourceLimit,
                                                  "vector calculus allocation failed",
                                                  operation);
    } catch (const detail::UnsupportedDifferentiation& error) {
        return VectorCalculusFieldResult::failure(
            CasErrc::UnsupportedExpression, error.what(), operation);
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
    if (!valid) {
        return VectorCalculusExprResult::failure(valid.error());
    }
    auto step = context.consume_steps(vars.size() * 2, operation);
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }
    try {
        return vector_calculus_wrap_expr(laplacian(f, vars), operation);
    } catch (const std::bad_alloc&) {
        return VectorCalculusExprResult::failure(CasErrc::ResourceLimit,
                                                 "vector calculus allocation failed",
                                                 operation);
    } catch (const detail::UnsupportedDifferentiation& error) {
        return VectorCalculusExprResult::failure(
            CasErrc::UnsupportedExpression, error.what(), operation);
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
        if (!first) {
            continue;
        }
        auto second = first->differentiate(var);
        if (!second) {
            continue;
        }
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

}
