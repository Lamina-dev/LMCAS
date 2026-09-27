#include "internal/transform_support.hpp"
#include "internal/normalization_utils.hpp"
#include <cmath>
#include <exception>

namespace LMCAS {
using namespace transform_detail;

static Result<void> te_validate_convolution_inputs(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& var,
    ComputationContext& context,
    const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!f || !LMCAS::detail::node(f) || !g || !LMCAS::detail::node(g)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "convolution expressions cannot be null",
                                     operation);
    }
    if (var.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "convolution variable name cannot be empty",
                                     operation);
    }
    return Result<void>::success();
}

static TransformEngineResult te_convolve_core(
    const std::shared_ptr<SymbolicExpr>& f, const std::shared_ptr<SymbolicExpr>& g,
    const std::string& variable, ComputationContext& context);

static TransformEngineResult te_convolve_sum(
    const AddNode& addition, const std::shared_ptr<SymbolicExpr>& other,
    bool sum_is_left, const std::string& variable, ComputationContext& context) {
    std::shared_ptr<SymbolicExpr> accumulated;
    for (const auto& operand : addition.operands()) {
        auto term = detail::make_expression_ptr(operand);
        auto transformed = sum_is_left
            ? te_convolve_core(term, other, variable, context)
            : te_convolve_core(other, term, variable, context);
        if (!transformed) return transformed;
        auto expression = transformed.value().value.expression;
        accumulated = accumulated ? SymbolicExpr::add(accumulated, expression) : expression;
    }
    return te_wrap_transform_result(
        accumulated ? accumulated->simplify() : SymbolicExpr::number(0), "convolve");
}

static TransformEngineResult te_convolve_core(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& variable,
    ComputationContext& context)
{
    constexpr const char* operation = "convolve";
    auto step = context.consume_steps(1, operation);
    if (!step) return TransformEngineResult::failure(step.error());

    if (te_is_zero(f) || te_is_zero(g)) {
        return te_wrap_transform_result(SymbolicExpr::number(0), operation);
    }

    auto left_sum = std::dynamic_pointer_cast<const AddNode>(detail::node(f));
    if (left_sum) {
        return te_convolve_sum(*left_sum, g, true, variable, context);
    }
    auto right_sum = std::dynamic_pointer_cast<const AddNode>(detail::node(g));
    if (right_sum) {
        return te_convolve_sum(*right_sum, f, false, variable, context);
    }

    auto left_parts = te_split_coeff(f, variable);
    if (te_depends_on(left_parts.second, variable) &&
        left_parts.second->to_string() != f->to_string()) {
        auto inner =
            te_convolve_core(left_parts.second, g, variable, context);
        if (!inner) return inner;
        auto scaled = SymbolicExpr::multiply(
            left_parts.first, inner.value().value.expression)->simplify();
        return te_wrap_transform_result(scaled, operation);
    }
    auto right_parts = te_split_coeff(g, variable);
    if (te_depends_on(right_parts.second, variable) &&
        right_parts.second->to_string() != g->to_string()) {
        auto inner =
            te_convolve_core(f, right_parts.second, variable, context);
        if (!inner) return inner;
        auto scaled = SymbolicExpr::multiply(
            right_parts.first, inner.value().value.expression)->simplify();
        return te_wrap_transform_result(scaled, operation);
    }

    const auto left = te_match_scaled_gaussian(f, variable);
    const auto right = te_match_scaled_gaussian(g, variable);
    if (!left || !right) {
        return TransformEngineResult::failure(
            CasErrc::Inconclusive,
            "当前双边卷积规则不支持该表达式对", operation);
    }

    const double rate_sum = left->rate + right->rate;
    const double output_rate = left->rate * right->rate / rate_sum;
    auto scale = SymbolicExpr::multiply(
        SymbolicExpr::multiply(left->coefficient, right->coefficient),
        SymbolicExpr::sqrt(SymbolicExpr::divide(
            SymbolicExpr::variable("pi"),
            SymbolicExpr::number(rate_sum))));
    auto square = SymbolicExpr::power(
        SymbolicExpr::variable(variable), SymbolicExpr::number(2));
    auto exponent = SymbolicExpr::multiply(
        SymbolicExpr::number(-output_rate), square);
    auto gaussian = SymbolicExpr::exp(exponent);
    return te_wrap_transform_result(
        SymbolicExpr::multiply(scale, gaussian)->simplify(), operation);
}

TransformEngineResult convolve_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& var,
    ComputationContext& context) {
    const std::string operation = "convolve";
    auto valid =
        te_validate_convolution_inputs(f, g, var, context, operation);
    if (!valid) return TransformEngineResult::failure(valid.error());
    try {
        return te_convolve_core(f, g, var, context);
    } catch (const std::bad_alloc&) {
        return TransformEngineResult::failure(
            CasErrc::ResourceLimit, "双边卷积分配失败", operation);
    } catch (const std::exception& error) {
        return TransformEngineResult::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

TransformEngineResult convolve_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& g,
    const std::string& var) {
    ComputationContext context;
    return convolve_checked(f, g, var, context);
}

}
