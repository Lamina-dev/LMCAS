#include "internal/expression_transform.hpp"
#include "internal/integration_support.hpp"

namespace LMCAS {

namespace {
int integration_by_parts_score(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
        auto e = LMCAS::detail::expression_from_node(node);
        if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
            if (fn->type() == FunctionNode::FuncType::Ln || fn->type() == FunctionNode::FuncType::Log) { return 1; }
            if (fn->type() == FunctionNode::FuncType::ArcSin || fn->type() == FunctionNode::FuncType::ArcTan) { return 2; }
            if (fn->type() == FunctionNode::FuncType::Sin || fn->type() == FunctionNode::FuncType::Cos) { return 4; }
            if (fn->type() == FunctionNode::FuncType::Exp) { return 5; }
        }
        if (!depends_on_integration_variable(e, var)) { return 10; }
        if (std::dynamic_pointer_cast<const VariableNode>(node)) { return 3; }
        if (const auto* power = dynamic_cast<const PowerNode*>(node.get())) {
            const auto exponent = detail::exact_rational_value(power->exponent());
            const auto* base = dynamic_cast<const VariableNode*>(power->base().get());
            if (base && base->name() == var && exponent &&
                exponent->is_integer() && *exponent > Rational(0)) { return 3; }
        }
        return 10;
}

int integration_by_parts_choice(
    const std::vector<std::shared_ptr<const SymbolicNode>>& ops, const std::string& var) {
    int best_u_idx = -1;
    int best_score = 10;
    if (ops.size() == 1) {
        int s = integration_by_parts_score(ops[0], var);
        if (s <= 2) best_u_idx = 0;
    } else {
        for (size_t i = 0; i < ops.size(); ++i) {
            if (!depends_on_integration_variable(LMCAS::detail::expression_from_node(ops[i]), var)) { continue; }
            int s = integration_by_parts_score(ops[i], var);
            if (s < best_score) {
                best_score = s;
                best_u_idx = (int)i;
            }
        }
    }
    return best_u_idx;
}

std::shared_ptr<SymbolicExpr> integration_differential_factor(
    const std::vector<std::shared_ptr<const SymbolicNode>>& ops, int best_u_idx) {
    std::vector<std::shared_ptr<const SymbolicNode>> dv_ops;
    for (size_t i = 0; i < ops.size(); ++i) {
        if ((int)i != best_u_idx) {
            dv_ops.push_back(ops[i]);
        }
    }
    if (dv_ops.empty()) {
        return SymbolicExpr::number(1);
    }
    if (dv_ops.size() == 1) {
        return detail::make_expression_ptr(detail::expression_from_node(dv_ops[0]));
    }
    return detail::make_expression_ptr(detail::make_node<MultiplyNode>(dv_ops));
}
}

Result<std::optional<detail::IntegralParts>> Integrator::try_integration_by_parts(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& computation, int depth) try {
    detail::RewriteBudget transform_budget(computation, (computation).limits().max_recursion_depth,
                                            (computation).limits().max_ast_nodes, "transform");

    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        ops = mul->operands();
    } else {
        ops.push_back(LMCAS::detail::node(expr));
    }
    int best_u_idx = integration_by_parts_choice(ops, var);
    if (best_u_idx == -1) { return std::optional<detail::IntegralParts>{}; }
    if (integration_by_parts_score(ops[best_u_idx], var) >= 4) {
        for (std::size_t i = 0; i < ops.size(); ++i) {
            if (i == static_cast<std::size_t>(best_u_idx) ||
                !expression_depends_on_variable(ops[i], var)) { continue; }
            const auto* function = dynamic_cast<const FunctionNode*>(ops[i].get());
            if (!function || (function->type() != FunctionNode::FuncType::Sin &&
                function->type() != FunctionNode::FuncType::Cos &&
                function->type() != FunctionNode::FuncType::Exp)) {
                return std::optional<detail::IntegralParts>{};
            }
        }
    }

    auto u = LMCAS::detail::expression_from_node(ops[best_u_idx]);
    auto dv = integration_differential_factor(ops, best_u_idx);

    auto integrated_v = integrate_recursive(*dv, var, computation, depth + 1);
    if (!integrated_v) {
        return Result<std::optional<detail::IntegralParts>>::failure(integrated_v.error());
    }
    if (!detail::node(integrated_v.value().remainder)->is_zero()) {
        return std::optional<detail::IntegralParts>{};
    }
    auto v = detail::make_expression_ptr(std::move(integrated_v.value().primitive));

    auto du = detail::make_expression_ptr(detail::normalize_node(
        detail::differentiate_node(detail::node(u), var, transform_budget), transform_budget));

    auto uv = SymbolicExpr::multiply(detail::make_expression_ptr(u), v);
    auto vdu = SymbolicExpr::multiply(v, du);
    vdu = detail::make_expression_ptr(detail::normalize_node(detail::node(vdu->cancel()), transform_budget));
    auto integrated_vdu = integrate_recursive(*vdu, var, computation, depth + 1);
    if (!integrated_vdu) {
        return Result<std::optional<detail::IntegralParts>>::failure(integrated_vdu.error());
    }
    auto parts = std::move(integrated_vdu.value());
    auto primitive = detail::expression_from_node(detail::normalize_node(
        detail::node(sym_sub(*uv, parts.primitive)), transform_budget));
    auto remainder = detail::expression_from_node(detail::normalize_node(
        detail::node(SymbolicExpr::multiply(
            SymbolicExpr::number(-1), detail::make_expression_ptr(parts.remainder))),
        transform_budget));
    return std::optional<detail::IntegralParts>{{std::move(primitive), std::move(remainder)}};
} catch (const CasError& error) {
    return Result<std::optional<detail::IntegralParts>>::failure(error);
}

}
