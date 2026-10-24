#include "internal/expression_analysis.hpp"
#include <stdexcept>
#include "internal/expression_transform.hpp"
#include "internal/integration_support.hpp"
#include "internal/assumption_facts.hpp"

namespace LMCAS {

namespace {
std::optional<SymbolicExpr> substitution_argument(const SymbolicExpr& term) {
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(detail::node(term))) {
        return detail::expression_from_node(power->base());
    }
    auto function = std::dynamic_pointer_cast<const FunctionNode>(detail::node(term));
    if (function && !function->arguments().empty()) {
        return detail::expression_from_node(function->arguments()[0]);
    }
    return std::nullopt;
}

Result<std::shared_ptr<SymbolicExpr>> substitution_primitive(
    const SymbolicExpr& candidate_term, const SymbolicExpr& u,
    const std::string& var, ComputationContext& context) {
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(detail::node(candidate_term))) {
        auto primitive = detail::integrate_power_primitive(u,
            detail::expression_from_node(power->exponent()), var, context);
        if (!primitive) {
            return Result<std::shared_ptr<SymbolicExpr>>::failure(primitive.error());
        }
        return primitive.value() ? detail::make_expression_ptr(*primitive.value()) : nullptr;
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(detail::node(candidate_term))) {
        if (function->type() == FunctionNode::FuncType::Cos) {
            return SymbolicExpr::sin(detail::make_expression_ptr(u));
        }
        if (function->type() == FunctionNode::FuncType::Sin) {
            return SymbolicExpr::multiply(SymbolicExpr::number(-1),
                SymbolicExpr::cos(detail::make_expression_ptr(u)));
        }
        if (function->type() == FunctionNode::FuncType::Exp) {
            return SymbolicExpr::exp(detail::make_expression_ptr(u));
        }
    }
    return nullptr;
}

Result<std::shared_ptr<SymbolicExpr>> substitute_candidate(
    const SymbolicExpr& expr, const SymbolicExpr& candidate_term,
    const std::string& var, ComputationContext& computation) try {
    detail::RewriteBudget transform_budget(computation, (computation).limits().max_recursion_depth,
                                            (computation).limits().max_ast_nodes, "transform");
    auto u = substitution_argument(candidate_term);
    if (!u || !depends_on_integration_variable(*u, var)) { return nullptr; }
            auto d_ptr = detail::make_expression_ptr(detail::differentiate_node(detail::node(*u), var, transform_budget));
            if (!d_ptr) { return nullptr; }
            auto du = detail::make_expression_ptr(detail::normalize_node(detail::node(d_ptr), transform_budget));
            if (du->is_zero()) { return nullptr; }

            auto f_u = candidate_term;
            auto term_times_du = SymbolicExpr::multiply(detail::make_expression_ptr(f_u), detail::make_expression_ptr(*du));
            auto ratio = detail::make_expression_ptr(detail::normalize_node(detail::node(SymbolicExpr::divide(detail::make_expression_ptr(expr), term_times_du)), transform_budget));
    if (depends_on_integration_variable(*ratio, var)) { return nullptr; }
    auto primitive = substitution_primitive(candidate_term, *u, var, computation);
    if (!primitive) { return primitive; }
    if (!primitive.value()) { return nullptr; }
    return SymbolicExpr::multiply(ratio, primitive.value());
} catch (const CasError& error) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(error);
} catch (const std::bad_alloc&) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(CasErrc::ResourceLimit, "allocation failed", "substitute_candidate");
} catch (const std::length_error& error) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(CasErrc::ResourceLimit, error.what(), "substitute_candidate");
} catch (const std::exception& error) {
    return Result<std::shared_ptr<SymbolicExpr>>::failure(CasErrc::InternalInvariant, error.what(), "substitute_candidate");
}
}

Result<std::optional<detail::IntegralParts>> Integrator::try_substitution(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& computation, int) {
    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        ops = mul->operands();
    } else {
        ops.push_back(LMCAS::detail::node(expr));
    }
    for (const auto& operand : ops) {
        auto candidate = detail::expression_from_node(operand);
        auto result = substitute_candidate(expr, candidate, var, computation);
        if (!result) {
            return Result<std::optional<detail::IntegralParts>>::failure(result.error());
        }
        if (result.value()) {
            return std::optional<detail::IntegralParts>{{*result.value(), zero_}};
        }
    }
    return std::optional<detail::IntegralParts>{};
}

namespace {
enum class LinearWrapperKind { None, Function, PowerOfFunction, PowerOfLinear };
struct LinearWrapper {
    LinearWrapperKind kind = LinearWrapperKind::None;
    std::shared_ptr<const FunctionNode> function;
    std::shared_ptr<const PowerNode> power;
    std::optional<SymbolicExpr> argument;
};

bool linear_substitutable_function(const FunctionNode& function) {
    using FT = FunctionNode::FuncType;
    switch (function.type()) {
        case FT::Infinity:
        case FT::Atan2:
        case FT::Log: return false;
        default: return function.arguments().size() == 1;
    }
}

LinearWrapper classify_linear_wrapper(const SymbolicExpr& expr, const std::string& var) {
    LinearWrapper wrapper;
    wrapper.function = std::dynamic_pointer_cast<const FunctionNode>(detail::node(expr));
    if (wrapper.function) {
        if (!linear_substitutable_function(*wrapper.function)) { return wrapper; }
        wrapper.argument = detail::expression_from_node(wrapper.function->arguments()[0]);
        wrapper.kind = LinearWrapperKind::Function;
        return wrapper;
    }
    wrapper.power = std::dynamic_pointer_cast<const PowerNode>(detail::node(expr));
    if (!wrapper.power) { return wrapper; }
    if (expression_depends_on_variable(wrapper.power->exponent(), var)) { return wrapper; }
    wrapper.function = std::dynamic_pointer_cast<const FunctionNode>(wrapper.power->base());
    if (wrapper.function) {
        if (!linear_substitutable_function(*wrapper.function)) { return wrapper; }
        wrapper.argument = detail::expression_from_node(wrapper.function->arguments()[0]);
        wrapper.kind = LinearWrapperKind::PowerOfFunction;
    } else {
        wrapper.argument = detail::expression_from_node(wrapper.power->base());
        wrapper.kind = LinearWrapperKind::PowerOfLinear;
    }
    return wrapper;
}

std::shared_ptr<SymbolicExpr> linear_reduced_expression(
    const LinearWrapper& wrapper, const std::shared_ptr<SymbolicExpr>& variable) {
    if (wrapper.kind == LinearWrapperKind::PowerOfLinear) {
        return detail::make_expression_ptr(detail::make_node<PowerNode>(
            detail::node(variable), wrapper.power->exponent()));
    }
    std::vector<std::shared_ptr<const SymbolicNode>> args{detail::node(variable)};
    auto function = detail::make_node<FunctionNode>(wrapper.function->type(), args);
    if (wrapper.kind == LinearWrapperKind::Function) {
        return detail::make_expression_ptr(function);
    }
    return detail::make_expression_ptr(
        detail::make_node<PowerNode>(function, wrapper.power->exponent()));
}
}


Result<std::optional<detail::IntegralParts>> Integrator::try_linear_substitution(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& computation, int) {
    detail::RewriteBudget budget(computation, computation.limits().max_recursion_depth,
                                 computation.limits().max_ast_nodes, "integrate.linear_substitution");
    auto wrapper = classify_linear_wrapper(expr, var);
    if (!wrapper.argument) { return std::optional<detail::IntegralParts>{}; }
    auto affine = detail::recognize_affine(*wrapper.argument, var, computation);
    if (!affine) {
        return Result<std::optional<detail::IntegralParts>>::failure(affine.error());
    }
    if (!affine.value()) { return std::optional<detail::IntegralParts>{}; }
    auto slope = detail::make_expression_ptr(detail::normalize_node(
        detail::node(affine.value()->slope), budget));
    auto offset = detail::make_expression_ptr(detail::normalize_node(
        detail::node(affine.value()->offset), budget));
    if (slope->is_one() && offset->is_zero()) {
        return std::optional<detail::IntegralParts>{};
    }
    std::optional<detail::AssumptionFacts> facts;
    if (computation.assumptions()) { facts.emplace(*computation.assumptions()); }
    auto nonzero = detail::query_nonzero_value(detail::node(slope),
        facts ? static_cast<const FactsQuery&>(*facts) : detail::no_facts(), Domain::Real, computation);
    if (!nonzero) {
        return Result<std::optional<detail::IntegralParts>>::failure(nonzero.error());
    }
    if (nonzero.value() != Tribool::True) { return std::optional<detail::IntegralParts>{}; }
    if (wrapper.kind == LinearWrapperKind::Function &&
        wrapper.function->type() == FunctionNode::FuncType::Abs) {
        auto real = detail::query_real_value(detail::node(*wrapper.argument),
            facts ? static_cast<const FactsQuery&>(*facts) : detail::no_facts(), computation);
        if (!real) {
            return Result<std::optional<detail::IntegralParts>>::failure(real.error());
        }
        if (real.value() != Tribool::True) { return std::optional<detail::IntegralParts>{}; }
        auto primitive = SymbolicExpr::divide(SymbolicExpr::multiply(
            detail::make_expression_ptr(*wrapper.argument), detail::make_expression_ptr(expr)),
            SymbolicExpr::multiply(SymbolicExpr::number(2), slope));
        return std::optional<detail::IntegralParts>{{detail::expression_from_node(
            detail::normalize_node(detail::node(primitive), budget)), zero_}};
    }

    auto variable = SymbolicExpr::variable(var);
    auto test = linear_reduced_expression(wrapper, variable);
    std::shared_ptr<SymbolicExpr> primitive;
    if (wrapper.kind == LinearWrapperKind::PowerOfLinear) {
        auto power = detail::integrate_power_primitive(*variable,
            detail::expression_from_node(wrapper.power->exponent()), var, computation);
        if (!power) {
            return Result<std::optional<detail::IntegralParts>>::failure(power.error());
        }
        if (power.value()) { primitive = detail::make_expression_ptr(*power.value()); }
    } else {
        const auto& query = facts ? static_cast<const FactsQuery&>(*facts) : detail::no_facts();
        auto real_variable = detail::query_real_value(detail::node(variable), query, computation);
        if (!real_variable) {
            return Result<std::optional<detail::IntegralParts>>::failure(real_variable.error());
        }
        if (real_variable.value() == Tribool::True) {
            auto real_argument = detail::query_real_value(
                detail::node(*wrapper.argument), query, computation);
            if (!real_argument) {
                return Result<std::optional<detail::IntegralParts>>::failure(real_argument.error());
            }
            if (real_argument.value() != Tribool::True) {
                return std::optional<detail::IntegralParts>{};
            }
        }
        primitive = lookup_integral(entries_, *test, var, budget,
            [](std::shared_ptr<SymbolicExpr> candidate) { return candidate; });
    }
    if (!primitive) { return std::optional<detail::IntegralParts>{}; }
    auto substituted = detail::make_expression_ptr(detail::substitute_node(
        detail::node(primitive), var, detail::node(*wrapper.argument), budget));
    auto result = detail::expression_from_node(detail::normalize_node(detail::node(
        SymbolicExpr::divide(substituted, slope)), budget));
    return std::optional<detail::IntegralParts>{{std::move(result), zero_}};
}

}
