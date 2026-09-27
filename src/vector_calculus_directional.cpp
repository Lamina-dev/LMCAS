#include "internal/vector_calculus_support.hpp"
#include "vector_calculus_differential.hpp"
#include "internal/visitors/differentiation_visitor.hpp"
#include <stdexcept>

namespace LMCAS {

using namespace vector_calculus_detail;

static std::shared_ptr<SymbolicExpr> vector_calculus_magnitude_squared(
    const VectorField& v)
{
    std::shared_ptr<SymbolicExpr> sum = nullptr;
    for (const auto& comp : v) {
        if (!comp) {
            continue;
        }
        auto sq = SymbolicExpr::power(comp, SymbolicExpr::number(2));
        sq = sq->simplify();
        if (!sum) {
            sum = sq;
        } else {
            sum = SymbolicExpr::add(sum, sq);
            sum = sum->simplify();
        }
    }
    if (!sum) {
        return SymbolicExpr::number(0);
    }
    return sum->simplify();
}

static std::shared_ptr<SymbolicExpr> vector_calculus_single_dir_deriv(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    const VectorField& unit_dir)
{
    std::shared_ptr<SymbolicExpr> sum = nullptr;
    for (size_t i = 0; i < vars.size(); ++i) {
        auto partial = f->differentiate(vars[i]);
        if (!partial) {
            continue;
        }
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
    if (!sum) {
        return SymbolicExpr::number(0);
    }
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
    if (!valid) {
        return VectorCalculusExprResult::failure(valid.error());
    }
    auto dir_valid = vector_calculus_validate_field_vars(direction, vars, context, operation);
    if (!dir_valid) {
        return VectorCalculusExprResult::failure(dir_valid.error());
    }
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
    if (!step) {
        return VectorCalculusExprResult::failure(step.error());
    }
    try {
        return vector_calculus_wrap_expr(
            directional_derivative(f, vars, direction, order), operation);
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

    auto mag_sq = vector_calculus_magnitude_squared(direction);
    if (!mag_sq || mag_sq->is_zero()) {
        return nullptr;
    }

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

    auto result = f;
    for (int k = 0; k < order; ++k) {
        result = vector_calculus_single_dir_deriv(result, vars, unit_dir);
        if (!result) {
            return SymbolicExpr::number(0);
        }
    }

    return result;
}

}
