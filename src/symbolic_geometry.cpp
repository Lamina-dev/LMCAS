#include <stdexcept>
#include "internal/expression_transform.hpp"
#include "symbolic_geometry.hpp"
#include "internal/integration_support.hpp"
#include "internal/symbolic_ast.hpp"
#include "symbolic.hpp"
#include <cmath>
#include <exception>

namespace LMCAS {

namespace {

Result<void> validate_geometry_inputs(const std::shared_ptr<SymbolicExpr>& function,
                                      const std::shared_ptr<SymbolicExpr>& lower,
                                      const std::shared_ptr<SymbolicExpr>& upper,
                                      ComputationContext& context,
                                      const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!function || !LMCAS::detail::node(function)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "profile expression cannot be null",
                                     operation);
    }
    if (!lower || !LMCAS::detail::node(lower) || !upper || !LMCAS::detail::node(upper)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "integration bounds cannot be null",
                                     operation);
    }
    return Result<void>::success();
}


ExpressionResult simplify_geometry_checked(std::shared_ptr<SymbolicExpr> expr,
                                             const std::string& operation,
                                             const std::string& message, ComputationContext& context)
try {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    if (!expr || !LMCAS::detail::node(expr)) {
        return ExpressionResult::failure(CasErrc::Inconclusive,
                                           message, operation);
    }
    auto simplified = detail::make_expression_ptr(detail::normalize_node(detail::node(expr), transform_budget));
    if (!simplified || !LMCAS::detail::node(simplified)) {
        return ExpressionResult::failure(CasErrc::Inconclusive,
                                           message, operation);
    }
    if (contains_unevaluated_integral(LMCAS::detail::node(simplified))) {
        return ExpressionResult::failure(CasErrc::Inconclusive,
                                           "geometry integral is outside the supported domain",
                                           operation);
    }
    return ExpressionResult::success(std::move(simplified));
} catch (const CasError& error) {
    return ExpressionResult::failure(error);
} catch (const std::bad_alloc&) {
    return ExpressionResult::failure(CasErrc::ResourceLimit, "allocation failed", operation);
} catch (const std::length_error& error) {
    return ExpressionResult::failure(CasErrc::ResourceLimit, error.what(), operation);
} catch (const std::exception& error) {
    return ExpressionResult::failure(CasErrc::InternalInvariant, error.what(), operation);
}

ExpressionResult differentiate_geometry_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& variable,
    const std::string& operation, ComputationContext& context)
try {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    auto derivative = detail::make_expression_ptr(
        detail::differentiate_node(detail::node(expr), variable, transform_budget));
    return simplify_geometry_checked(
        derivative, operation,
        "geometry derivative is outside the supported domain", context);
} catch (const CasError& error) {
    if (error.code == CasErrc::UnsupportedExpression) {
        return ExpressionResult::failure(CasErrc::Inconclusive,
            "geometry derivative is outside the supported domain", operation);
    }
    return ExpressionResult::failure(error);
} catch (const std::bad_alloc&) {
    return ExpressionResult::failure(CasErrc::ResourceLimit, "allocation failed", operation);
} catch (const std::length_error& error) {
    return ExpressionResult::failure(CasErrc::ResourceLimit, error.what(), operation);
} catch (const std::exception& error) {
    return ExpressionResult::failure(CasErrc::InternalInvariant, error.what(), operation);
}

ExpressionResult definite_integral_geometry_checked(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& lower_bound,
    const std::shared_ptr<SymbolicExpr>& upper_bound,
    const std::string& operation, ComputationContext& context)
try {

    if (!integrand || !LMCAS::detail::node(integrand)) {
        return ExpressionResult::failure(
            CasErrc::Inconclusive,
            "geometry integrand is outside the supported domain",
            operation);
    }
    Integrator integrator;
    auto integrated = integrator.integrate_def_checked(
        *integrand, variable, *lower_bound, *upper_bound, context);
    if (!integrated) { return ExpressionResult::failure(integrated.error()); }
    auto result = detail::make_expression_ptr(std::move(integrated.value()));
    if (contains_unevaluated_integral(detail::node(result))) {
        return ExpressionResult::failure(
            CasErrc::Inconclusive,
            "geometry integral is outside the supported domain",
            operation);
    }
    return simplify_geometry_checked(
        result, operation,
        "geometry definite integral result is outside the supported domain", context);
} catch (const CasError& error) {
    return ExpressionResult::failure(error);
} catch (const std::bad_alloc&) {
    return ExpressionResult::failure(CasErrc::ResourceLimit, "allocation failed", operation);
} catch (const std::length_error& error) {
    return ExpressionResult::failure(CasErrc::ResourceLimit, error.what(), operation);
} catch (const std::exception& error) {
    return ExpressionResult::failure(CasErrc::InternalInvariant, error.what(), operation);
}

std::shared_ptr<SymbolicExpr> squared_profile(const std::shared_ptr<SymbolicExpr>& f, ComputationContext& context) {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    if (f) {
        if (auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(f))) {
            if (func->type() == FunctionNode::FuncType::Sqrt && func->arguments().size() == 1) {
                return detail::make_expression_ptr(detail::normalize_node(func->arguments()[0], transform_budget));
            }
        }
        if (auto power = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(f))) {
            auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
            if (exponent && std::holds_alternative<Rational>(exponent->value())) {
                const auto& rational = std::get<Rational>(exponent->value());
                if (rational.get_numerator() == BigInt(1) &&
                    rational.get_denominator() == BigInt(2)) {
                    return detail::make_expression_ptr(detail::normalize_node(power->base(), transform_budget));
                }
            }
        }
    }
    return detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::power(f, SymbolicExpr::number(2))), transform_budget));
}

ExpressionResult volume_of_revolution_checked_impl(
    std::shared_ptr<SymbolicExpr> f,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    const std::string& variable,
    ComputationContext& context,
    const std::string& operation)
try {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    auto valid = validate_geometry_inputs(f, a, b, context, operation);
    if (!valid) return ExpressionResult::failure(valid.error());

    auto step = context.consume_steps(8, operation);
    if (!step) return ExpressionResult::failure(step.error());

    try {
        auto pi = LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<VariableNode>("pi", true));
        auto integrand = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::multiply(pi, squared_profile(f, context))), transform_budget));
        return definite_integral_geometry_checked(
            integrand, variable, a, b, operation, context);
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                           "symbolic geometry allocation failed",
                                           operation);
    } catch (const std::length_error& error) { return ExpressionResult::failure(CasErrc::ResourceLimit, error.what(), operation); } catch (const std::exception& e) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                           e.what(), operation);
    }
} catch (const CasError& error) {
    return ExpressionResult::failure(error);
} catch (const std::bad_alloc&) {
    return ExpressionResult::failure(CasErrc::ResourceLimit, "allocation failed", operation);
} catch (const std::length_error& error) {
    return ExpressionResult::failure(CasErrc::ResourceLimit, error.what(), operation);
} catch (const std::exception& error) {
    return ExpressionResult::failure(CasErrc::InternalInvariant, error.what(), operation);
}

ExpressionResult arc_length_checked_impl(
    std::shared_ptr<SymbolicExpr> f,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    const std::string& variable,
    ComputationContext& context,
    const std::string& operation)
{
    auto valid = validate_geometry_inputs(f, a, b, context, operation);
    if (!valid) return ExpressionResult::failure(valid.error());

    auto step = context.consume_steps(10, operation);
    if (!step) return ExpressionResult::failure(step.error());

    try {
        auto derivative_result = differentiate_geometry_checked(f, variable, operation, context);
        if (!derivative_result) return derivative_result;
        auto integrand = SymbolicExpr::sqrt(SymbolicExpr::add(
            SymbolicExpr::number(1),
            SymbolicExpr::power(derivative_result.value(), SymbolicExpr::number(2))));
        return definite_integral_geometry_checked(
            integrand, variable, a, b, operation, context);
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                           "symbolic geometry allocation failed",
                                           operation);
    } catch (const std::exception& e) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                           e.what(), operation);
    }
}

}

ExpressionResult volume_of_revolution_x_checked(
    std::shared_ptr<SymbolicExpr> fx,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    ComputationContext& context
) {
    return volume_of_revolution_checked_impl(
        std::move(fx), std::move(a), std::move(b), "x",
        context, "volume_of_revolution_x");
}

ExpressionResult volume_of_revolution_x_checked(
    std::shared_ptr<SymbolicExpr> fx,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
) {
    ComputationContext context;
    return volume_of_revolution_x_checked(std::move(fx), std::move(a), std::move(b),
                                          context);
}


ExpressionResult arc_length_x_checked(
    std::shared_ptr<SymbolicExpr> fx,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    ComputationContext& context
) {
    return arc_length_checked_impl(
        std::move(fx), std::move(a), std::move(b), "x",
        context, "arc_length_x");
}

ExpressionResult arc_length_x_checked(
    std::shared_ptr<SymbolicExpr> fx,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
) {
    ComputationContext context;
    return arc_length_x_checked(std::move(fx), std::move(a), std::move(b), context);
}


ExpressionResult volume_of_revolution_y_checked(
    std::shared_ptr<SymbolicExpr> fy,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    ComputationContext& context
) {
    return volume_of_revolution_checked_impl(
        std::move(fy), std::move(a), std::move(b), "y",
        context, "volume_of_revolution_y");
}

ExpressionResult volume_of_revolution_y_checked(
    std::shared_ptr<SymbolicExpr> fy,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
) {
    ComputationContext context;
    return volume_of_revolution_y_checked(std::move(fy), std::move(a), std::move(b),
                                          context);
}


ExpressionResult arc_length_y_checked(
    std::shared_ptr<SymbolicExpr> fy,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    ComputationContext& context
) {
    return arc_length_checked_impl(
        std::move(fy), std::move(a), std::move(b), "y",
        context, "arc_length_y");
}

ExpressionResult arc_length_y_checked(
    std::shared_ptr<SymbolicExpr> fy,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b
) {
    ComputationContext context;
    return arc_length_y_checked(std::move(fy), std::move(a), std::move(b), context);
}


}
