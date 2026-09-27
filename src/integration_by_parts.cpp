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
        if (std::dynamic_pointer_cast<const PowerNode>(node)) { return 3; }
        return 10;
}

int integration_by_parts_choice(
    const std::vector<std::shared_ptr<const SymbolicNode>>& ops, const std::string& var) {
    int best_u_idx = -1;
    int best_score = 100;
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

bool integration_primitive_complete(const std::shared_ptr<SymbolicExpr>& primitive) {
    return primitive && detail::node(primitive) &&
        !contains_unevaluated_integral(detail::node(primitive));
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

Result<std::shared_ptr<SymbolicExpr>> IBPStrategy::try_integrate_raw(
    const SymbolicExpr& expr, const std::string& var, Integrator& ctx,
    ComputationContext& computation, int depth) {

    std::vector<std::shared_ptr<const SymbolicNode>> ops;
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        ops = mul->operands();
    } else {
        ops.push_back(LMCAS::detail::node(expr));
    }
    int best_u_idx = integration_by_parts_choice(ops, var);
    if (best_u_idx == -1) { return nullptr; }

    auto u = LMCAS::detail::expression_from_node(ops[best_u_idx]);
    auto dv = integration_differential_factor(ops, best_u_idx);

    auto integrated_v =
        ctx.integrate_recursive(*dv, var, computation, depth + 1);
    if (!integrated_v) { return integrated_v; }
    auto v = std::move(integrated_v.value());
    if (!integration_primitive_complete(v)) {
        return nullptr;
    }

    auto du_ptr = u.differentiate(var);
    if (!du_ptr) { return nullptr; }
    auto simplified_du = du_ptr->simplify();
    if (!simplified_du || !LMCAS::detail::node(simplified_du)) { return nullptr; }
    auto du = detail::make_expression_ptr(*simplified_du);

    auto uv = SymbolicExpr::multiply(detail::make_expression_ptr(u), v);
    auto vdu = SymbolicExpr::multiply(v, du);
    if (!uv || !vdu) { return nullptr; }
    vdu = vdu->cancel()->simplify();
    auto integrated_vdu =
        ctx.integrate_recursive(*vdu, var, computation, depth + 1);
    if (!integrated_vdu) { return integrated_vdu; }
    auto int_vdu = std::move(integrated_vdu.value());
    if (!integration_primitive_complete(int_vdu)) {
        return nullptr;
    }

    return sym_sub(*uv, *int_vdu);
}

}
