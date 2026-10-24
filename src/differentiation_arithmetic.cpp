#include "internal/visitors/differentiation_visitor.hpp"

namespace LMCAS {

void DifferentiationVisitor::visit(const NumberNode&) {
    set_result(SymbolicFactory::create_number(BigInt(0)));
}

void DifferentiationVisitor::visit(const VariableNode& node) {
    if (!node.is_constant() && node.name() == var) {
        set_result(SymbolicFactory::create_number(BigInt(1)));
    } else {
        set_result(SymbolicFactory::create_number(BigInt(0)));
    }
}

void DifferentiationVisitor::visit(const AddNode& node) {
    std::vector<std::shared_ptr<const SymbolicNode>> diff_ops;
    std::size_t size = 1;
    rewrite_budget()->append_size(size, node.operands().size());
    diff_ops.reserve(node.operands().size());
    for (const auto& op : node.operands()) {
        op->accept(*this);
        append_result(diff_ops, size, result);
    }
    set_result(SymbolicFactory::create_add(std::move(diff_ops), rewrite_budget()));
}

void DifferentiationVisitor::visit(const MultiplyNode& node) {
    if (node.operands().empty()) {
         set_result(SymbolicFactory::create_number(BigInt(0)));
         return;
    }

    std::vector<std::shared_ptr<const SymbolicNode>> sum_terms;
    std::size_t sum_size = 1;

    for (size_t i = 0; i < node.operands().size(); ++i) {
        rewrite_budget()->consume();

        node.operands()[i]->accept(*this);
        auto d_term = result;

        if (d_term->is_zero()) continue;

        std::vector<std::shared_ptr<const SymbolicNode>> prod_terms;
        std::size_t product_size = 1;
        rewrite_budget()->append_size(product_size, node.operands().size());
        prod_terms.reserve(node.operands().size());

        for (size_t j = 0; j < node.operands().size(); ++j) {
            rewrite_budget()->consume();
            append_result(prod_terms, product_size,
                          i == j ? d_term : node.operands()[j]);
        }
        append_result(sum_terms, sum_size,
                      SymbolicFactory::create_multiply(std::move(prod_terms), rewrite_budget()));
    }

    if (sum_terms.empty()) {
         set_result(SymbolicFactory::create_number(BigInt(0)));
    } else {
         set_result(SymbolicFactory::create_add(std::move(sum_terms), rewrite_budget()));
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

        auto n_minus_1 = SymbolicFactory::create_add(
            {n, SymbolicFactory::create_number(BigInt(-1))}, rewrite_budget());

        auto u_pow = LMCAS::detail::make_node<PowerNode>(u, n_minus_1);

        set_result(SymbolicFactory::create_multiply({n, u_pow, du}, rewrite_budget()));
    } else {

        auto ln_u = LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Ln, std::vector<std::shared_ptr<const SymbolicNode>>{u});

        auto t1 = SymbolicFactory::create_multiply({dv, ln_u}, rewrite_budget());

        auto u_inv = LMCAS::detail::make_node<PowerNode>(u, SymbolicFactory::create_number(BigInt(-1)));
        auto t2 = SymbolicFactory::create_multiply({v, du, u_inv}, rewrite_budget());

        auto sum = SymbolicFactory::create_add({t1, t2}, rewrite_budget());

        auto u_pow_v = LMCAS::detail::make_node<PowerNode>(u, v);
        set_result(SymbolicFactory::create_multiply({u_pow_v, sum}, rewrite_budget()));
    }
}

}
