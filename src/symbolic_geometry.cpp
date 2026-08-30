#include "../include/symbolic_geometry.hpp"
#include "symbolic_ast.hpp"
#include "../include/symbolic.hpp"
#include <cmath>
#include <exception>

namespace lamina {

namespace {

Result<void> validate_geometry_inputs(const std::shared_ptr<SymbolicExpr>& function,
                                      const std::shared_ptr<SymbolicExpr>& lower,
                                      const std::shared_ptr<SymbolicExpr>& upper,
                                      ComputationContext& context,
                                      const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!function || !lamina::detail::node(function)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "profile expression cannot be null",
                                     operation);
    }
    if (!lower || !lamina::detail::node(lower) || !upper || !lamina::detail::node(upper)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "integration bounds cannot be null",
                                     operation);
    }
    return Result<void>::success();
}

bool contains_unevaluated_integral(
    const std::shared_ptr<const SymbolicNode>& node,
    std::size_t = 0) {
    return lamina::detail::contains_node_type<IntegralNode>(node);
}

ExpressionResult simplify_geometry_checked(std::shared_ptr<SymbolicExpr> expr,
                                             const std::string& operation,
                                             const std::string& message)
{
    if (!expr || !lamina::detail::node(expr)) {
        return ExpressionResult::failure(CasErrc::Inconclusive,
                                           message, operation);
    }
    auto simplified = expr->simplify();
    if (!simplified || !lamina::detail::node(simplified)) {
        return ExpressionResult::failure(CasErrc::Inconclusive,
                                           message, operation);
    }
    if (contains_unevaluated_integral(lamina::detail::node(simplified))) {
        return ExpressionResult::failure(CasErrc::Inconclusive,
                                           "geometry integral is outside the supported domain",
                                           operation);
    }
    return ExpressionResult::success(std::move(simplified));
}

ExpressionResult differentiate_geometry_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& variable,
    const std::string& operation)
{
    try {
        auto derivative = expr->differentiate(variable);
        return simplify_geometry_checked(
            derivative, operation,
            "geometry derivative is outside the supported domain");
    } catch (const std::exception&) {
        return ExpressionResult::failure(
            CasErrc::Inconclusive,
            "geometry derivative is outside the supported domain",
            operation);
    }
}

ExpressionResult definite_integral_geometry_checked(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& lower_bound,
    const std::shared_ptr<SymbolicExpr>& upper_bound,
    const std::string& operation)
{
    if (!integrand || !lamina::detail::node(integrand)) {
        return ExpressionResult::failure(
            CasErrc::Inconclusive,
            "geometry integrand is outside the supported domain",
            operation);
    }
    auto integral = integrand->integrate(variable);
    if (!integral || !lamina::detail::node(integral) ||
        contains_unevaluated_integral(lamina::detail::node(integral))) {
        return ExpressionResult::failure(
            CasErrc::Inconclusive,
            "geometry integral is outside the supported domain",
            operation);
    }
    auto upper = integral->substitute(variable, upper_bound);
    auto lower = integral->substitute(variable, lower_bound);
    auto result = SymbolicExpr::add(
        upper,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), lower));
    return simplify_geometry_checked(
        result, operation,
        "geometry definite integral result is outside the supported domain");
}

std::shared_ptr<SymbolicExpr> squared_profile(const std::shared_ptr<SymbolicExpr>& f) {
    if (f) {
        if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(f))) {
            if (func->type() == FunctionNode::FuncType::Sqrt && func->arguments().size() == 1) {
                return lamina::detail::make_expression_ptr(func->arguments()[0])->simplify();
            }
        }
        if (auto power = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(f))) {
            auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
            if (exponent && std::holds_alternative<Rational>(exponent->value())) {
                const auto& rational = std::get<Rational>(exponent->value());
                if (rational.get_numerator() == BigInt(1) &&
                    rational.get_denominator() == BigInt(2)) {
                    return lamina::detail::make_expression_ptr(power->base())->simplify();
                }
            }
        }
    }
    return SymbolicExpr::power(f, SymbolicExpr::number(2))->simplify();
}

ExpressionResult volume_of_revolution_checked_impl(
    std::shared_ptr<SymbolicExpr> f,
    std::shared_ptr<SymbolicExpr> a,
    std::shared_ptr<SymbolicExpr> b,
    const std::string& variable,
    ComputationContext& context,
    const std::string& operation)
{
    auto valid = validate_geometry_inputs(f, a, b, context, operation);
    if (!valid) return ExpressionResult::failure(valid.error());

    auto step = context.consume_steps(8, operation);
    if (!step) return ExpressionResult::failure(step.error());

    try {
        auto pi = SymbolicExpr::variable("pi");
        auto integrand = SymbolicExpr::multiply(pi, squared_profile(f))->simplify();
        return definite_integral_geometry_checked(
            integrand, variable, a, b, operation);
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                           "symbolic geometry allocation failed",
                                           operation);
    } catch (const std::exception& e) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                           e.what(), operation);
    }
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
        auto derivative_result = differentiate_geometry_checked(f, variable, operation);
        if (!derivative_result) return derivative_result;
        auto integrand = SymbolicExpr::sqrt(SymbolicExpr::add(
            SymbolicExpr::number(1),
            SymbolicExpr::power(derivative_result.value(), SymbolicExpr::number(2))));
        return definite_integral_geometry_checked(
            integrand, variable, a, b, operation);
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
