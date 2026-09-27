#include "internal/symbolic_ast.hpp"
#include "internal/normalization_utils.hpp"
#include <optional>

namespace LMCAS {

std::shared_ptr<const SymbolicNode> SymbolicFactory::create_number(const ::LMCAS::BigInt& v) { return LMCAS::detail::make_node<NumberNode>(v); }
std::shared_ptr<const SymbolicNode> SymbolicFactory::create_number(const ::LMCAS::Rational& v) { return LMCAS::detail::make_node<NumberNode>(v); }
std::shared_ptr<const SymbolicNode> SymbolicFactory::create_number(lmmc_real_t v) { return LMCAS::detail::make_node<NumberNode>(v); }
std::shared_ptr<const SymbolicNode> SymbolicFactory::create_variable(const std::string& name) { return LMCAS::detail::make_node<VariableNode>(name); }

std::shared_ptr<const SymbolicNode> SymbolicFactory::create_add(
    std::vector<std::shared_ptr<const SymbolicNode>> ops,
    detail::RewriteBudget* budget) {
    if (ops.empty()) {
        normalization_check_children(budget, 1);
        return create_number(::LMCAS::BigInt(0));
    }
    std::size_t nodes = 0;
    for (const auto& op : ops) {
        if (!op) {
            throw std::invalid_argument("create_add operand cannot be null");
        }
        if (budget) {
            nodes = budget->append_size(nodes, budget->measure(op));
        }
    }
    if (ops.size() == 1) {
        return ops[0];
    }

    std::vector<std::shared_ptr<const SymbolicNode>> flat_ops;
    nodes = 0;
    flat_ops.reserve(ops.size());
    for (const auto& op : ops) {
        if (op->is_zero()) {
            continue;
        }
        if (auto add = std::dynamic_pointer_cast<const AddNode>(op)) {
            if (budget) {
                for (const auto& child : add->operands()) {
                    nodes = budget->append_size(nodes, budget->measure(child));
                }
            }
            flat_ops.insert(flat_ops.end(), add->operands().begin(), add->operands().end());
        } else {
            normalization_append(budget, nodes, flat_ops, op);
        }
    }

    if (flat_ops.empty()) {
        normalization_check_children(budget, 1);
        return create_number(::LMCAS::BigInt(0));
    }
    if (flat_ops.size() == 1) {
        return flat_ops[0];
    }
    normalization_check_arithmetic<AddNode>(budget, flat_ops, nodes);
    return LMCAS::detail::make_node<AddNode>(std::move(flat_ops));
}

namespace {

detail::SymbolicNodePtr discardable_zero_factor(
    const std::vector<detail::SymbolicNodePtr>& operands, detail::RewriteBudget* budget) {
    const auto zero = std::find_if(operands.begin(), operands.end(), [](const auto& operand) {
        return operand->is_zero();
    });
    if (zero == operands.end()) { return nullptr; }
    std::optional<ComputationContext> local_context;
    if (!budget) { local_context.emplace(); }
    auto& context = budget ? budget->context() : *local_context;
    if (std::all_of(operands.begin(), operands.end(), [&](const auto& operand) {
            return normalization_can_discard(operand, detail::no_facts(), Domain::Real, context);
        })) {
        return *zero;
    }
    return nullptr;
}

}

std::shared_ptr<const SymbolicNode> SymbolicFactory::create_multiply(
    std::vector<std::shared_ptr<const SymbolicNode>> ops,
    detail::RewriteBudget* budget) {
    if (ops.empty()) {
        normalization_check_children(budget, 1);
        return create_number(::LMCAS::BigInt(1));
    }

    std::size_t nodes = 0;
    for (const auto& op : ops) {
        if (!op) {
            throw std::invalid_argument("create_multiply operand cannot be null");
        }
        if (budget) {
            nodes = budget->append_size(nodes, budget->measure(op));
        }
    }
    if (auto zero = discardable_zero_factor(ops, budget)) {
        return zero;
    }

    if (ops.size() == 1) {
        return ops[0];
    }

    std::vector<std::shared_ptr<const SymbolicNode>> flat_ops;
    nodes = 0;
    flat_ops.reserve(ops.size());
    for (const auto& op : ops) {
        if (op->is_one()) {
            continue;
        }
        if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(op)) {
            if (budget) {
                for (const auto& child : mul->operands()) {
                    nodes = budget->append_size(nodes, budget->measure(child));
                }
            }
            flat_ops.insert(flat_ops.end(), mul->operands().begin(), mul->operands().end());
        } else {
            normalization_append(budget, nodes, flat_ops, op);
        }
    }

    if (flat_ops.empty()) {
        normalization_check_children(budget, 1);
        return create_number(::LMCAS::BigInt(1));
    }
    if (flat_ops.size() == 1) {
        return flat_ops[0];
    }
    normalization_check_arithmetic<MultiplyNode>(budget, flat_ops, nodes);
    return LMCAS::detail::make_node<MultiplyNode>(std::move(flat_ops));
}


std::shared_ptr<const SymbolicNode> SymbolicFactory::create_power(std::shared_ptr<const SymbolicNode> base, std::shared_ptr<const SymbolicNode> exponent) {
    if (!base || !exponent) {
        throw std::invalid_argument("create_power operands cannot be null");
    }
    ComputationContext context;
    if (exponent->is_zero()) {
        if (normalization_can_discard(base, detail::no_facts(), Domain::Real, context) &&
            normalization_proved(detail::query_nonzero_value(base, detail::no_facts(), Domain::Complex, context))) {
            return create_number(::LMCAS::BigInt(1));
        }
        return LMCAS::detail::make_node<PowerNode>(std::move(base), std::move(exponent));
    }
    if (exponent->is_one()) {
        return base;
    }
    if (base->is_zero()) {
        if (normalization_can_discard(exponent, detail::no_facts(), Domain::Real, context) &&
            normalization_proved(detail::query_positive_value(exponent, detail::no_facts(), context))) {
            return create_number(::LMCAS::BigInt(0));
        }
        return LMCAS::detail::make_node<PowerNode>(std::move(base), std::move(exponent));
    }
    if (base->is_one() && normalization_can_discard(exponent, detail::no_facts(), Domain::Real, context)) {
        return create_number(::LMCAS::BigInt(1));
    }
    return LMCAS::detail::make_node<PowerNode>(std::move(base), std::move(exponent));
}

std::shared_ptr<const SymbolicNode> SymbolicFactory::create_complex(std::shared_ptr<const SymbolicNode> real, std::shared_ptr<const SymbolicNode> imag) {
    if (!real || !imag) {
        throw std::invalid_argument("create_complex operands cannot be null");
    }
    if (imag->is_zero()) {
        return real;
    }
    return LMCAS::detail::make_node<ComplexNode>(std::move(real), std::move(imag));
}

}
