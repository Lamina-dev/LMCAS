#include "internal/visitors/differentiation_visitor.hpp"

namespace LMCAS {

void DifferentiationVisitor::visit(const NumberNode&) {
    result = SymbolicFactory::create_number(BigInt(0));
}

void DifferentiationVisitor::visit(const VariableNode& node) {
    if (node.name() == var) {
        result = SymbolicFactory::create_number(BigInt(1));
    } else {
        result = SymbolicFactory::create_number(BigInt(0));
    }
}

void DifferentiationVisitor::visit(const AddNode& node) {
    std::vector<std::shared_ptr<const SymbolicNode>> diff_ops;
    for (const auto& op : node.operands()) {
        op->accept(*this);
        diff_ops.push_back(result);
    }
    result = SymbolicFactory::create_add(diff_ops);
}

void DifferentiationVisitor::visit(const MultiplyNode& node) {
    if (node.operands().empty()) {
         result = SymbolicFactory::create_number(BigInt(0));
         return;
    }

    std::vector<std::shared_ptr<const SymbolicNode>> sum_terms;

    for (size_t i = 0; i < node.operands().size(); ++i) {

        node.operands()[i]->accept(*this);
        auto d_term = result;

        if (d_term->is_zero()) continue;

        std::vector<std::shared_ptr<const SymbolicNode>> prod_terms;

        for (size_t j = 0; j < node.operands().size(); ++j) {
            if (i == j) {
                prod_terms.push_back(d_term);
            } else {
                prod_terms.push_back(node.operands()[j]);
            }
        }
        sum_terms.push_back(SymbolicFactory::create_multiply(prod_terms));
    }

    if (sum_terms.empty()) {
         result = SymbolicFactory::create_number(BigInt(0));
    } else {
         result = SymbolicFactory::create_add(sum_terms);
    }
}

void DifferentiationVisitor::visit(const PowerNode& node) {

    node.base()->accept(*this);
    auto du = result;
    node.exponent()->accept(*this);
    auto dv = result;

    auto u = node.base();
    auto v = node.exponent();

    if (dv->is_zero()) {

        auto n = v;

        auto n_minus_1 = SymbolicFactory::create_add({n, SymbolicFactory::create_number(BigInt(-1))});

        auto u_pow = LMCAS::detail::make_node<PowerNode>(u, n_minus_1);

        result = SymbolicFactory::create_multiply({n, u_pow, du});
    } else {

        auto ln_u = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln, std::vector<std::shared_ptr<const SymbolicNode>>{u});

        auto t1 = SymbolicFactory::create_multiply({dv, ln_u});

        auto u_inv = LMCAS::detail::make_node<PowerNode>(u, SymbolicFactory::create_number(BigInt(-1)));
        auto t2 = SymbolicFactory::create_multiply({v, du, u_inv});

        auto sum = SymbolicFactory::create_add({t1, t2});

        auto u_pow_v = LMCAS::detail::make_node<PowerNode>(u, v);
        result = SymbolicFactory::create_multiply({u_pow_v, sum});
    }
}

}
