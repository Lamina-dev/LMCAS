#include "internal/visitors/normalization_visitor.hpp"
#include "internal/normalization_utils.hpp"

namespace LMCAS {
namespace {
struct FactorAccum {
    std::shared_ptr<const NumberNode> exponent;
    std::vector<std::shared_ptr<const SymbolicNode>> originals;
    std::size_t nodes;
};

std::shared_ptr<const NumberNode> collectable_power_exponent(
    const std::shared_ptr<const SymbolicNode>& factor,
    std::shared_ptr<const SymbolicNode>& base, const FactsQuery& facts, Domain domain,
    ComputationContext& context) {
    auto exponent = detail::make_node<NumberNode>(BigInt(1));
    const auto power = std::dynamic_pointer_cast<const PowerNode>(factor);
    if (!power) { return exponent; }
    const auto numeric = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
    BigInt integer;
    if (!try_get_integer_value(numeric, integer)) { return exponent; }
    if (integer <= BigInt(0) &&
        !(normalization_can_discard(power->base(), facts, domain, context) &&
          normalization_proved(detail::query_nonzero_value(power->base(), facts, Domain::Complex, context)))) {
        return exponent;
    }
    base = power->base();
    return detail::make_node<NumberNode>(integer);
}

class ProductAccumulator {
public:
    bool add_number(const std::shared_ptr<const SymbolicNode>& factor) {
        const auto number = std::dynamic_pointer_cast<const NumberNode>(factor);
        if (!number) {
            return false;
        }
        constant_ = multiply_numbers(constant_, number);
        return true;
    }

    void add_base(const std::shared_ptr<const SymbolicNode>& factor,
                  const std::shared_ptr<const SymbolicNode>& base,
                  const std::shared_ptr<const NumberNode>& exponent, detail::RewriteBudget* budget) {
        auto it = bases_.find(base);
        if (it == bases_.end()) {
            const auto nodes = budget ? budget->measure(factor) : 0;
            bases_.emplace(base, FactorAccum{exponent, {factor}, nodes});
        } else {
            it->second.exponent = add_numbers(it->second.exponent, exponent);
            normalization_append(budget, it->second.nodes, it->second.originals, factor);
        }
    }

    std::shared_ptr<const SymbolicNode> finish(
        const FactsQuery& facts, Domain domain, ComputationContext& context, bool expanding,
        detail::RewriteBudget* budget) const {
        if (budget) { budget->require_nodes(1); }
        if (constant_->is_zero() && factors_defined(facts, domain, context)) { return constant_; }
        std::vector<std::shared_ptr<const SymbolicNode>> operands;
        std::size_t nodes = 0;
        if (!constant_->is_one()) {
            normalization_append(budget, nodes, operands, constant_);
        }
        std::vector<std::shared_ptr<const SymbolicNode>> variables;
        std::size_t variable_nodes = 0;
        collect_variables(facts, domain, budget, variables, variable_nodes, context);
        std::sort(variables.begin(), variables.end(), normalization_factor_less);
        if (budget) {
            nodes = budget->append_size(nodes, variable_nodes);
        }
        operands.insert(operands.end(), variables.begin(), variables.end());
        if (expanding && operands.size() > 1 &&
            std::dynamic_pointer_cast<const NumberNode>(operands.back())) {
            std::rotate(operands.begin(), operands.end() - 1, operands.end());
        }
        return make_normalized_multiply_node(operands, budget);
    }

private:
    void collect_variables(const FactsQuery& facts, Domain domain,
        detail::RewriteBudget* budget,
        std::vector<std::shared_ptr<const SymbolicNode>>& variables,
        std::size_t& nodes, ComputationContext& context) const {
        for (const auto& [base, accumulated] : bases_) {
            const auto& exponent = accumulated.exponent;
            if (exponent->is_zero()) {
                if (normalization_can_discard(base, facts, domain, context) &&
                    normalization_proved(detail::query_nonzero_value(base, facts, Domain::Complex, context))) {
                    continue;
                }
                for (const auto& original : accumulated.originals) {
                    normalization_append(budget, nodes, variables, original);
                }
            } else if (exponent->is_one()) {
                normalization_append(budget, nodes, variables, base);
            } else {
                normalization_check_children(budget, 1, base, exponent);
                normalization_append(budget, nodes, variables,
                                     detail::make_node<PowerNode>(base, exponent));
            }
        }
    }

    bool factors_defined(const FactsQuery& facts, Domain domain, ComputationContext& context) const {
        for (const auto& [base, accumulated] : bases_) {
            for (const auto& original : accumulated.originals) {
                if (!normalization_can_discard(original, facts, domain, context)) { return false; }
            }
        }
        return true;
    }

    std::shared_ptr<const NumberNode> constant_ = detail::make_node<NumberNode>(BigInt(1));
    std::map<std::shared_ptr<const SymbolicNode>, FactorAccum, NodeCompare> bases_;
};

bool has_matrix_factor(const std::vector<std::shared_ptr<const SymbolicNode>>& operands) {
    for (const auto& operand : operands) {
        if (std::dynamic_pointer_cast<const MatrixNode>(operand)) {
            return true;
        }
        const auto power = std::dynamic_pointer_cast<const PowerNode>(operand);
        if (power && std::dynamic_pointer_cast<const MatrixNode>(power->base())) {
            return true;
        }
    }
    return false;
}
}

std::shared_ptr<const SymbolicNode> NormalizationVisitor::collect_product(
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands) {
    ProductAccumulator accumulated;
    for (const auto& operand : operands) {
        if (accumulated.add_number(operand)) {
            continue;
        }
        auto base = operand;
        auto exponent = collectable_power_exponent(operand, base, facts_, domain_, context_);
        if (base == operand) {
            const auto power = std::dynamic_pointer_cast<const PowerNode>(operand);
            if (power && can_combine_nested_power(*power, exponent)) {
                base = power->base();
                exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
            }
        }
        accumulated.add_base(operand, base, exponent, rewrite_budget());
    }
    return accumulated.finish(facts_, domain_, context_, false, rewrite_budget());
}

void NormalizationVisitor::visit(const MultiplyNode& node) {
    std::vector<std::shared_ptr<const SymbolicNode>> operands;
    std::size_t nodes = 0;
    for (const auto& operand : node.operands()) {
        operand->accept(*this);
        if (const auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(result)) {
            for (const auto& child : multiply->operands()) {
                normalization_append(rewrite_budget(), nodes, operands, child);
            }
        } else {
            normalization_append(rewrite_budget(), nodes, operands, result);
        }
    }
    if (has_matrix_factor(operands)) {
        set_result(normalize_matrix_product(operands));
        return;
    }
    set_result(collect_product(operands));
}

std::shared_ptr<const SymbolicNode> NormalizationVisitor::expand_complex_product(
    const std::shared_ptr<const SymbolicNode>& lhs,
    const std::shared_ptr<const SymbolicNode>& rhs) {
    const auto left = std::dynamic_pointer_cast<const ComplexNode>(lhs);
    const auto right = std::dynamic_pointer_cast<const ComplexNode>(rhs);
    if (!left && !right) {
        return nullptr;
    }
    std::shared_ptr<const SymbolicNode> real;
    std::shared_ptr<const SymbolicNode> imag;
    if (left && right) {
        auto ac = expand_product(left->real(), right->real());
        auto bd = expand_product(left->imag(), right->imag());
        auto ad = expand_product(left->real(), right->imag());
        auto bc = expand_product(left->imag(), right->real());
        auto negative_bd = expand_product(detail::make_node<NumberNode>(BigInt(-1)), bd);
        real = SymbolicFactory::create_add({ac, negative_bd}, rewrite_budget());
        imag = SymbolicFactory::create_add({ad, bc}, rewrite_budget());
    } else if (left) {
        real = expand_product(left->real(), rhs);
        imag = expand_product(left->imag(), rhs);
    } else {
        real = expand_product(lhs, right->real());
        imag = expand_product(lhs, right->imag());
    }
    NormalizationVisitor visitor(context_, facts_, domain_, rewrite_budget());
    real->accept(visitor);
    auto normalized_real = visitor.get_result();
    imag->accept(visitor);
    normalization_check_children(rewrite_budget(), 1, normalized_real, visitor.get_result());
    return SymbolicFactory::create_complex(normalized_real, visitor.get_result());
}

std::shared_ptr<const SymbolicNode> NormalizationVisitor::expand_product(
    const std::shared_ptr<const SymbolicNode>& lhs,
    const std::shared_ptr<const SymbolicNode>& rhs) {
    detail::RewriteScope scope(rewrite_budget());
    if (auto complex = expand_complex_product(lhs, rhs)) {
        return complex;
    }
    const auto left = std::dynamic_pointer_cast<const AddNode>(lhs);
    const auto right = std::dynamic_pointer_cast<const AddNode>(rhs);
    if (left || right) {
        return distribute_product(lhs, rhs);
    }
    return collect_expanded_product(lhs, rhs);
}

std::shared_ptr<const SymbolicNode> NormalizationVisitor::distribute_product(
    const std::shared_ptr<const SymbolicNode>& lhs,
    const std::shared_ptr<const SymbolicNode>& rhs) {
    const auto left = std::dynamic_pointer_cast<const AddNode>(lhs);
    const auto right = std::dynamic_pointer_cast<const AddNode>(rhs);
    std::size_t left_count = 1, right_count = 1;
    if (left) {
        left_count = left->operands().size();
    }
    if (right) {
        right_count = right->operands().size();
    }
    normalization_check_count(rewrite_budget(), left_count, right_count, false);
    std::vector<std::shared_ptr<const SymbolicNode>> terms;
    std::size_t nodes = 0;
    for (std::size_t i = 0; i < left_count; ++i) {
        const auto& a = left ? left->operands()[i] : lhs;
        for (std::size_t j = 0; j < right_count; ++j) {
            const auto& b = right ? right->operands()[j] : rhs;
            normalization_append(rewrite_budget(), nodes, terms, expand_product(a, b));
        }
    }
    normalization_check_arithmetic<AddNode>(rewrite_budget(), terms, nodes);
    return detail::make_node<AddNode>(terms);
}

std::shared_ptr<const SymbolicNode> NormalizationVisitor::collect_expanded_product(
    const std::shared_ptr<const SymbolicNode>& lhs,
    const std::shared_ptr<const SymbolicNode>& rhs) {
    ProductAccumulator accumulated;
    const auto process = [&](const std::shared_ptr<const SymbolicNode>& factor) {
        if (accumulated.add_number(factor)) {
            return;
        }
        auto base = factor;
        auto exponent = collectable_power_exponent(factor, base, facts_, domain_, context_);
        accumulated.add_base(factor, base, exponent, rewrite_budget());
    };
    for (const auto* node : {&lhs, &rhs}) {
        const auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(*node);
        if (!multiply) {
            process(*node);
            continue;
        }
        for (const auto& operand : multiply->operands()) {
            process(operand);
        }
    }
    return accumulated.finish(facts_, domain_, context_, true, rewrite_budget());
}
}
