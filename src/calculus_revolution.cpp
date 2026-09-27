#include "integration.hpp"
#include "internal/calculus_utils_support.hpp"
#include "internal/integration_support.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS {

using namespace calculus_utils_detail;
static ExpressionResult calculus_utils_try_symbolic_definite(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context)
{
    Integrator integrator;
    auto integrated = integrator.integrate_def_checked(
        *integrand, var, *a, *b, context);
    if (!integrated) {
        return ExpressionResult::failure(integrated.error());
    }
    SymbolicExpr result = std::move(integrated.value());

    if (contains_unevaluated_integral(LMCAS::detail::node(result))) {
        return std::shared_ptr<SymbolicExpr>{};
    }
    auto res = LMCAS::detail::make_expression_ptr(result);
    auto simplified = res->simplify();
    if (simplified &&
        contains_unevaluated_integral(LMCAS::detail::node(simplified))) {
        return std::shared_ptr<SymbolicExpr>{};
    }
    return simplified ? simplified : res;
}
ExpressionResult surface_area_revolution_x_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context)
{
    const std::string operation = "surface_area_revolution_x";
    auto input = calculus_utils_validate_expr_bounds(f, var, a, b, context, operation);
    if (!input) {
        return ExpressionResult::failure(input.error());
    }
    auto step = context.consume_steps(8, operation);
    if (!step) {
        return ExpressionResult::failure(step.error());
    }

    try {
        auto f_prime = f->differentiate(var);
        if (!f_prime || !LMCAS::detail::node(f_prime)) {
            return ExpressionResult::failure(
                CasErrc::Inconclusive,
                "surface derivative is outside the supported domain",
                operation);
        }

        auto f_prime_sq = SymbolicExpr::power(f_prime, SymbolicExpr::number(2));
        auto one_plus_fp_sq = SymbolicExpr::add(SymbolicExpr::number(1), f_prime_sq);
        auto arc_factor = SymbolicExpr::sqrt(one_plus_fp_sq);
        auto abs_f = calculus_utils_make_abs(f);
        if (!abs_f || !LMCAS::detail::node(abs_f)) {
            return ExpressionResult::failure(
                CasErrc::InternalInvariant,
                "absolute value node construction failed",
                operation);
        }
        auto integrand = SymbolicExpr::multiply(abs_f, arc_factor);

        auto integral = calculus_utils_try_symbolic_definite(
            integrand, var, a, b, context);
        if (!integral) {
            return integral;
        }
        if (!integral.value() ||
            !LMCAS::detail::node(integral.value())) {
            return ExpressionResult::failure(
                CasErrc::Inconclusive,
                "surface area integral could not be evaluated exactly",
                operation);
        }

        auto two_pi = SymbolicExpr::multiply(
            SymbolicExpr::number(2), SymbolicExpr::variable("pi"));
        auto result = SymbolicExpr::multiply(
            two_pi, std::move(integral.value()));
        auto simplified = result->simplify();
        return ExpressionResult::success(simplified ? simplified : result);
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                           "surface area allocation failed",
                                           operation);
    } catch (const std::exception& e) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                           e.what(), operation);
    }
}

ExpressionResult surface_area_revolution_x_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b)
{
    ComputationContext context;
    return surface_area_revolution_x_checked(f, var, a, b, context);
}


ExpressionResult surface_area_revolution_y_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context)
{
    const std::string operation = "surface_area_revolution_y";
    auto input = calculus_utils_validate_expr_bounds(f, var, a, b, context, operation);
    if (!input) {
        return ExpressionResult::failure(input.error());
    }
    auto step = context.consume_steps(8, operation);
    if (!step) {
        return ExpressionResult::failure(step.error());
    }

    try {
        auto f_prime = f->differentiate(var);
        if (!f_prime || !LMCAS::detail::node(f_prime)) {
            return ExpressionResult::failure(
                CasErrc::Inconclusive,
                "surface derivative is outside the supported domain",
                operation);
        }

        auto f_prime_sq = SymbolicExpr::power(f_prime, SymbolicExpr::number(2));
        auto one_plus_fp_sq = SymbolicExpr::add(SymbolicExpr::number(1), f_prime_sq);
        auto arc_factor = SymbolicExpr::sqrt(one_plus_fp_sq);
        auto var_expr = SymbolicExpr::variable(var);
        auto abs_var = calculus_utils_make_abs(var_expr);
        if (!abs_var || !LMCAS::detail::node(abs_var)) {
            return ExpressionResult::failure(
                CasErrc::InternalInvariant,
                "absolute value node construction failed",
                operation);
        }
        auto integrand = SymbolicExpr::multiply(abs_var, arc_factor);

        auto integral = calculus_utils_try_symbolic_definite(
            integrand, var, a, b, context);
        if (!integral) {
            return integral;
        }
        if (!integral.value() ||
            !LMCAS::detail::node(integral.value())) {
            return ExpressionResult::failure(
                CasErrc::Inconclusive,
                "surface area integral could not be evaluated exactly",
                operation);
        }

        auto two_pi = SymbolicExpr::multiply(
            SymbolicExpr::number(2), SymbolicExpr::variable("pi"));
        auto result = SymbolicExpr::multiply(
            two_pi, std::move(integral.value()));
        auto simplified = result->simplify();
        return ExpressionResult::success(simplified ? simplified : result);
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                           "surface area allocation failed",
                                           operation);
    } catch (const std::exception& e) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                           e.what(), operation);
    }
}

ExpressionResult surface_area_revolution_y_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b)
{
    ComputationContext context;
    return surface_area_revolution_y_checked(f, var, a, b, context);
}

}
