#include "internal/transform_support.hpp"
#include "internal/normalization_utils.hpp"
#include "assumption_context.hpp"
#include "inference_engine.hpp"
#include "residual_verification.hpp"
#include <cmath>
#include <exception>

namespace LMCAS {
using namespace transform_detail;


static std::shared_ptr<SymbolicExpr> te_fourier_unevaluated(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& t,
    const std::string& omega) {
    return detail::make_expression_ptr(detail::make_node<TransformNode>(
        TransformNode::TransformType::Fourier, detail::node(f)->clone(), t,
        SymbolicFactory::create_variable(omega)));
}

static bool te_is_abs_variable(const std::shared_ptr<const SymbolicNode>& node,
                               const std::string& variable) {
    auto function = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!function || function->type() != FunctionNode::FuncType::Abs ||
        function->arguments().size() != 1) return false;
    auto argument = std::dynamic_pointer_cast<const VariableNode>(function->arguments()[0]);
    return argument && !argument->is_constant() && argument->name() == variable;
}
static std::shared_ptr<SymbolicExpr> te_abs_decay(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& t) {
    auto function = std::dynamic_pointer_cast<const FunctionNode>(detail::node(f));
    if (!function || function->type() != FunctionNode::FuncType::Exp ||
        function->arguments().size() != 1) return nullptr;
    auto argument = detail::make_expression_ptr(function->arguments()[0]);
    auto parts = te_split_coeff(argument, t);
    if (!te_is_abs_variable(detail::node(parts.second), t)) return nullptr;
    return detail::make_expression_ptr(SymbolicFactory::create_multiply({
        SymbolicFactory::create_number(BigInt(-1)), detail::node(parts.first)}))->simplify();
}

static std::shared_ptr<SymbolicExpr> te_decay_condition(
    const std::shared_ptr<SymbolicExpr>& rate) {
    auto product = std::dynamic_pointer_cast<const MultiplyNode>(detail::node(rate));
    if (!product) return te_gt_condition(rate, SymbolicExpr::number(0));
    std::vector<std::shared_ptr<const SymbolicNode>> factors;
    for (const auto& factor : product->operands()) {
        auto number = std::dynamic_pointer_cast<const NumberNode>(factor);
        if (!number || !number->is_positive()) factors.push_back(factor);
    }
    auto lhs = detail::make_expression_ptr(SymbolicFactory::create_multiply(std::move(factors)));
    return te_gt_condition(lhs, SymbolicExpr::number(0));
}

static bool te_accept_decay(const std::shared_ptr<SymbolicExpr>& rate,
                            const AssumptionContext& assumptions,
                            std::vector<std::shared_ptr<SymbolicExpr>>& conditions,
                            ComputationContext& context) {
    InferenceEngine inference(assumptions);
    auto real = inference.query_real_checked(*rate, context);
    if (!real) {
        if (real.error().code == CasErrc::DomainError) return false;
        throw real.error();
    }
    if (real.value() != Tribool::True) {
        return false;
    }
    auto positive = inference.query_positive_checked(*rate, context);
    if (!positive) {
        throw positive.error();
    }
    if (positive.value() == Tribool::True) {
        return true;
    }
    auto nonpositive = inference.query_nonpositive_checked(*rate, context);
    if (!nonpositive) {
        throw nonpositive.error();
    }
    if (positive.value() == Tribool::False || nonpositive.value() == Tribool::True) {
        return false;
    }
    conditions.push_back(te_decay_condition(rate));
    return true;
}

static std::shared_ptr<SymbolicExpr> te_abs_fourier(
    const std::shared_ptr<SymbolicExpr>& rate, const std::string& omega) {
    auto two = SymbolicFactory::create_number(BigInt(2));
    auto denominator = SymbolicFactory::create_add({
        SymbolicFactory::create_power(detail::node(rate), two),
        SymbolicFactory::create_power(SymbolicFactory::create_variable(omega), two)});
    auto reciprocal = SymbolicFactory::create_power(
        denominator, SymbolicFactory::create_number(BigInt(-1)));
    return detail::make_expression_ptr(SymbolicFactory::create_multiply({
        two, detail::node(rate), reciprocal}))->simplify();
}

static std::shared_ptr<SymbolicExpr> te_gaussian_fourier(
    double rate, const std::string& omega) {
    auto coefficient = SymbolicExpr::number(std::sqrt(3.14159265358979323846 / rate));
    auto exponent = SymbolicExpr::multiply(SymbolicExpr::number(-1.0 / (4.0 * rate)),
        SymbolicExpr::power(SymbolicExpr::variable(omega), SymbolicExpr::number(2)));
    return SymbolicExpr::multiply(coefficient, SymbolicExpr::exp(exponent));
}


static std::shared_ptr<SymbolicExpr> fourier_transform_core(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& t,
    const std::string& omega, ComputationContext& context,
    const AssumptionContext& assumptions,
    std::vector<std::shared_ptr<SymbolicExpr>>& conditions);

static std::shared_ptr<SymbolicExpr> te_fourier_sum(
    const AddNode& addition, const std::string& t, const std::string& omega,
    ComputationContext& context, const AssumptionContext& assumptions,
    std::vector<std::shared_ptr<SymbolicExpr>>& conditions) {
    std::shared_ptr<SymbolicExpr> result;
    for (const auto& operand : addition.operands()) {
        auto transformed = fourier_transform_core(detail::make_expression_ptr(operand),
            t, omega, context, assumptions, conditions);
        if (!transformed) return nullptr;
        result = result ? SymbolicExpr::add(result, transformed) : transformed;
    }
    return result;
}

static std::shared_ptr<SymbolicExpr> fourier_transform_core(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& t,
    const std::string& omega, ComputationContext& context,
    const AssumptionContext& assumptions,
    std::vector<std::shared_ptr<SymbolicExpr>>& conditions) {
    auto step = context.consume_steps(1, "fourier_transform.recursive");
    if (!step) {
        return nullptr;
    }
    if (!f || !detail::node(f)) {
        return nullptr;
    }
    if (!te_depends_on(f, t)) {
        return te_fourier_unevaluated(f, t, omega);
    }
    auto addition = std::dynamic_pointer_cast<const AddNode>(detail::node(f));
    if (addition) {
        return te_fourier_sum(*addition, t, omega, context, assumptions, conditions);
    }
    auto [coefficient, body] = te_split_coeff(f, t);
    if (!coefficient->is_one()) {
        auto transformed = fourier_transform_core(body, t, omega, context, assumptions, conditions);
        if (transformed) {
            return SymbolicExpr::multiply(coefficient, transformed);
        }
    }
    if (const auto gaussian_rate = te_match_pure_gaussian(f, t)) {
        return te_gaussian_fourier(*gaussian_rate, omega);
    }
    auto decay = te_abs_decay(f, t);
    if (decay) {
        if (te_accept_decay(decay, assumptions, conditions, context)) {
            return te_abs_fourier(decay, omega);
        }
        return te_fourier_unevaluated(f, t, omega);
    }
    return te_fourier_unevaluated(f, t, omega);
}

TransformEngineResult fourier_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& t,
    const std::string& omega,
    ComputationContext& context) {
    const std::string operation = "fourier_transform";
    auto valid = te_validate_expr_vars(f, t, omega, context, operation);
    if (!valid) return TransformEngineResult::failure(valid.error());
    auto step = context.consume_steps(8, operation);
    if (!step) return TransformEngineResult::failure(step.error());
    try {
        AssumptionContext empty_assumptions;
        const auto& assumptions = context.assumptions() ? *context.assumptions() : empty_assumptions;
        std::vector<std::shared_ptr<SymbolicExpr>> conditions;
        auto expression = fourier_transform_core(f, t, omega, context, assumptions, conditions);
        auto final_access = context.consume_steps(0, operation);
        if (!final_access) {
            return TransformEngineResult::failure(final_access.error());
        }
        auto result = te_wrap_transform_result(std::move(expression), operation);
        if (result) result.value().value.conditions = std::move(conditions);
        return result;
    } catch (const CasError& error) {
        return TransformEngineResult::failure(error);
    } catch (const std::bad_alloc&) {
        return TransformEngineResult::failure(CasErrc::ResourceLimit,
                                              "transform allocation failed",
                                              operation);
    } catch (const std::exception& e) {
        return TransformEngineResult::failure(CasErrc::InternalInvariant,
                                              e.what(), operation);
    }
}

TransformEngineResult fourier_transform_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& t,
    const std::string& omega) {
    ComputationContext context;
    return fourier_transform_checked(f, t, omega, context);
}

static std::shared_ptr<SymbolicExpr> inverse_fourier_transform_core(
    const std::shared_ptr<SymbolicExpr>& F, const std::string& omega,
    const std::string& t, ComputationContext& context) {
    auto step = context.consume_steps(1, "inverse_fourier_transform.recursive");
    if (!step) return nullptr;
    if (!F || !LMCAS::detail::node(F)) return nullptr;
    auto tv = SymbolicExpr::variable(t);
    if (!te_depends_on(F, omega)) {
        return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
            TransformNode::TransformType::InverseFourier, LMCAS::detail::node(F)->clone(), omega,
            SymbolicFactory::create_variable(t)));
    }
    auto add_node = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(F));
    if (add_node) {
        std::shared_ptr<SymbolicExpr> result;
        for (auto& op : add_node->operands()) {
            auto ift = inverse_fourier_transform_core(
                LMCAS::detail::make_expression_ptr(op), omega, t, context);
            if (!ift) {
                return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
                    TransformNode::TransformType::InverseFourier, LMCAS::detail::node(F)->clone(), omega,
                    SymbolicFactory::create_variable(t)));
            }
            result = result ? SymbolicExpr::add(result, ift) : ift;
        }
        return result;
    }

    if (!contains_inexact_number(detail::node(F))) {
        auto canonical = te_abs_fourier(SymbolicExpr::number(1), omega);
        if (detail::node(F->simplify())->equals(*detail::node(canonical))) {
            auto abs_t = detail::make_expression_ptr(detail::make_node<FunctionNode>(
                FunctionNode::FuncType::Abs,
                std::vector<std::shared_ptr<const SymbolicNode>>{detail::node(tv)}));
            return SymbolicExpr::exp(SymbolicExpr::multiply(
                SymbolicExpr::number(-1), abs_t));
        }
    }

    if (const auto rate = te_match_pure_gaussian(F, omega)) {
        const double coefficient =
            1.0 / std::sqrt(4.0 * 3.14159265358979323846 * *rate);
        auto exponent = SymbolicExpr::multiply(
            SymbolicExpr::number(-1.0 / (4.0 * *rate)),
            SymbolicExpr::power(tv, SymbolicExpr::number(2)));
        return SymbolicExpr::multiply(
            SymbolicExpr::number(coefficient), SymbolicExpr::exp(exponent));
    }
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<TransformNode>(
        TransformNode::TransformType::InverseFourier, LMCAS::detail::node(F)->clone(), omega,
        SymbolicFactory::create_variable(t)));
}

TransformEngineResult inverse_fourier_transform_checked(
    const std::shared_ptr<SymbolicExpr>& F,
    const std::string& omega,
    const std::string& t,
    ComputationContext& context) {
    const std::string operation = "inverse_fourier_transform";
    auto valid = te_validate_expr_vars(F, omega, t, context, operation);
    if (!valid) return TransformEngineResult::failure(valid.error());
    auto step = context.consume_steps(8, operation);
    if (!step) return TransformEngineResult::failure(step.error());
    try {
        auto expression = inverse_fourier_transform_core(
            F, omega, t, context);
        auto final_access = context.consume_steps(0, operation);
        if (!final_access) {
            return TransformEngineResult::failure(final_access.error());
        }
        if (contains_inexact_number(detail::node(expression))) {
            return TransformEngineResult::failure(
                CasErrc::Inconclusive,
                "inverse Fourier round trip requires an exact candidate",
                operation);
        }
        auto result = te_wrap_transform_result(expression, operation);
        if (!result) return result;
        AssumptionContext empty_assumptions;
        const auto& assumptions = context.assumptions() ? *context.assumptions() : empty_assumptions;
        std::vector<std::shared_ptr<SymbolicExpr>> conditions;
        auto round_trip = fourier_transform_core(
            expression, t, omega, context, assumptions, conditions);
        if (!round_trip || !conditions.empty() ||
            te_contains_transform(LMCAS::detail::node(round_trip))) {
            return TransformEngineResult::failure(
                CasErrc::Inconclusive,
                "inverse Fourier round trip is not proved", operation);
        }
        auto verified = check_equivalent(round_trip, F, context);
        if (!verified) return TransformEngineResult::failure(verified.error());
        if (!std::holds_alternative<ProvedZeroResidual>(verified.value())) {
            return TransformEngineResult::failure(
                CasErrc::Inconclusive,
                "inverse Fourier round trip is not proved", operation);
        }
        result.value().certificate =
            ExactRoundTripProof{SymbolicExpr::number(0)};
        return result;
    } catch (const CasError& error) {
        return TransformEngineResult::failure(error);
    } catch (const std::bad_alloc&) {
        return TransformEngineResult::failure(CasErrc::ResourceLimit,
                                              "transform allocation failed",
                                              operation);
    } catch (const std::exception& e) {
        return TransformEngineResult::failure(CasErrc::InternalInvariant,
                                              e.what(), operation);
    }
}

TransformEngineResult inverse_fourier_transform_checked(
    const std::shared_ptr<SymbolicExpr>& F,
    const std::string& omega,
    const std::string& t) {
    ComputationContext context;
    return inverse_fourier_transform_checked(F, omega, t, context);
}

}
