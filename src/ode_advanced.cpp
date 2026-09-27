#include "limit_result.hpp"
/**
 * @file ode_advanced.cpp
 * @brief 通过受检接口用常数变易法构造解。
 */
#include "symbolic_ode_engine.hpp"
#include "internal/symbolic_ast.hpp"
#include "symbolic.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "internal/ode_support.hpp"
#include "integrator.hpp"
#include "residual_verification.hpp"
#include <cmath>
#include <memory>
#include <string>

namespace LMCAS {

static ODESolutionResult solve_variation_of_parameters_core(
    const std::shared_ptr<SymbolicExpr>& y1,
    const std::shared_ptr<SymbolicExpr>& y2,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& x,
    ComputationContext& context);

ODESolutionResult solve_variation_of_parameters_checked(
    const std::shared_ptr<SymbolicExpr>& y1,
    const std::shared_ptr<SymbolicExpr>& y2,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& x,
    ComputationContext& context)
{
    const std::string operation = "solve_variation_of_parameters";
    auto valid = validate_ode_three_expr_one_var(y1, y2, g, x, context, operation);
    if (!valid) {
        return ODESolutionResult::failure(valid.error());
    }

    auto budget = context.consume_steps(24, operation);
    if (!budget) {
        return ODESolutionResult::failure(budget.error());
    }

    try {
        auto solved = solve_variation_of_parameters_core(
            y1, y2, g, x, context);
        if (!solved) {
            return solved;
        }
        return wrap_ode_solution(
            std::move(solved.value()),
            ODEType::HigherOrder_ConstCoeff,
            operation);
    } catch (const std::bad_alloc&) {
        return ODESolutionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while applying variation of parameters",
            operation);
    } catch (const std::exception& ex) {
        return ODESolutionResult::failure(
            CasErrc::InternalInvariant,
            ex.what(),
            operation);
    }
}

ODESolutionResult solve_variation_of_parameters_checked(
    const std::shared_ptr<SymbolicExpr>& y1,
    const std::shared_ptr<SymbolicExpr>& y2,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& x)
{
    ComputationContext context;
    return solve_variation_of_parameters_checked(y1, y2, g, x, context);
}

static ODESolutionResult solve_variation_of_parameters_core(
    const std::shared_ptr<SymbolicExpr>& y1,
    const std::shared_ptr<SymbolicExpr>& y2,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& x,
    ComputationContext& context)
{
    ODESolution result;
    result.method_used = ODEType::HigherOrder_ConstCoeff;
    result.constants = {};

    if (!y1 || !y2 || !g) {
        result.general_solution = nullptr;
        return result;
    }

    /// 计算 y₁' 和 y₂'
    auto y1_prime = y1->differentiate(x);
    auto y2_prime = y2->differentiate(x);

    if (!y1_prime || !y2_prime) {
        result.general_solution = nullptr;
        return result;
    }

    /// 计算 Wronskian: W = y₁·y₂' - y₂·y₁'
    auto term1 = SymbolicExpr::multiply(y1, y2_prime);
    auto term2 = SymbolicExpr::multiply(y2, y1_prime);
    auto W = SymbolicExpr::add(term1,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), term2));
    W = W->simplify();
    EqvOptions trig_options;
    trig_options.profile = EqvProfile::TrigBasic;
    auto unit_wronskian = check_equivalent(
        W, SymbolicExpr::number(1), context, trig_options);
    if (!unit_wronskian) {
        return ODESolutionResult::failure(unit_wronskian.error());
    }
    if (std::holds_alternative<ProvedZeroResidual>(
            unit_wronskian.value())) {
        W = SymbolicExpr::number(1);
    }

    if (!W || W->is_zero()) {
        result.general_solution = nullptr;
        return result;
    }

    /// 计算 u₁' = -y₂·g(x)/W
    auto u1_prime = SymbolicExpr::divide(
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
            SymbolicExpr::multiply(y2, g)),
        W);
    u1_prime = u1_prime->simplify();

    /// 计算 u₂' = y₁·g(x)/W
    auto u2_prime = SymbolicExpr::divide(
        SymbolicExpr::multiply(y1, g),
        W);
    u2_prime = u2_prime->simplify();

    Integrator integrator;
    auto u1_value =
        integrator.integrate_checked(*u1_prime, x, context);
    if (!u1_value) {
        return ODESolutionResult::failure(u1_value.error());
    }
    auto u2_value =
        integrator.integrate_checked(*u2_prime, x, context);
    if (!u2_value) {
        return ODESolutionResult::failure(u2_value.error());
    }
    auto u1 = std::make_shared<SymbolicExpr>(
        std::move(u1_value.value()));
    auto u2 = std::make_shared<SymbolicExpr>(
        std::move(u2_value.value()));

    /// 特解: y_p = u₁·y₁ + u₂·y₂
    auto y_p = SymbolicExpr::add(
        SymbolicExpr::multiply(u1, y1),
        SymbolicExpr::multiply(u2, y2));
    y_p = y_p->simplify();

    result.general_solution = y_p;
    return result;
}

}
