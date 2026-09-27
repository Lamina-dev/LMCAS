#include "internal/equivalence_support.hpp"

#include "internal/expression_analysis.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

namespace LMCAS::equivalence_detail {

ExprPtr rewritten_complex(const ExprPtr& real, const ExprPtr& imag,
                           const detail::SymbolicNodePtr& original,
                           detail::RewriteBudget& budget) {
    budget.append_size(budget.append_size(1, budget.measure(detail::node(real))),
                       budget.measure(detail::node(imag)));
    ExprPtr value;
    try {
        value = detail::make_expression_ptr(SymbolicFactory::create_complex(
            detail::node(real), detail::node(imag)));
    } catch (const std::invalid_argument&) {
        return detail::make_expression_ptr(original);
    }
    return normalize_equivalence(value, budget);
}

namespace {

ExprPtr canonicalize_variable(const VariableNode& variable,
                               const detail::SymbolicNodePtr& original,
                               detail::RewriteBudget& budget) {
    if (!detail::is_imaginary_unit_name(variable.name())) {
        return detail::make_expression_ptr(original);
    }
    budget.require_nodes(3);
    try {
        const auto zero = SymbolicExpr::number(0);
        const auto one = SymbolicExpr::number(1);
        return detail::make_expression_ptr(SymbolicFactory::create_complex(
            detail::node(zero), detail::node(one)));
    } catch (const std::invalid_argument&) {
        return detail::make_expression_ptr(original);
    }
}

ExprPtr canonicalize_sum(const AddNode& add,
                         const detail::SymbolicNodePtr& original,
                         detail::RewriteBudget& budget) {
    std::vector<detail::SymbolicNodePtr> operands;
    bool changed = false;
    operands.reserve(add.operands().size());
    for (const auto& operand : add.operands()) {
        auto canonical = canonicalize_complex_product(
            *detail::make_expression_ptr(operand), budget);
        changed = changed || detail::node(canonical) != operand;
        operands.push_back(detail::node(canonical));
    }
    if (!changed) {
        budget.measure(original);
        return detail::make_expression_ptr(original);
    }
    std::size_t size = 0;
    for (const auto& operand : operands) {
        size = budget.append_size(size, budget.measure(operand));
    }
    return normalize_equivalence(detail::make_expression_ptr(
        SymbolicFactory::create_add(std::move(operands), &budget)), budget);
}

ExprPtr canonicalize_power(const PowerNode& power,
                            const detail::SymbolicNodePtr& original,
                            detail::RewriteBudget& budget) {
    auto exponent = exact_small_integer_node(power.exponent(), 0, 16);
    if (exponent) {
        auto base = canonicalize_complex_product(
            *detail::make_expression_ptr(power.base()), budget);
        if (std::dynamic_pointer_cast<const ComplexNode>(detail::node(base))) {
            budget.require_nodes(1);
            auto result = SymbolicExpr::number(1);
            for (int index = 0; index < *exponent; ++index) {
                normalization_check_arithmetic<MultiplyNode>(
                    &budget, 0, detail::node(result), detail::node(base));
                result = canonicalize_complex_product(
                    *SymbolicExpr::multiply(result, base), budget);
            }
            return normalize_equivalence(result, budget);
        }
    }
    budget.measure(original);
    return detail::make_expression_ptr(original);
}

struct ComplexProduct {
    ExprPtr real = SymbolicExpr::number(1);
    ExprPtr imag = SymbolicExpr::number(0);
    bool saw_complex = false;

    void accumulate(const ExprPtr& factor, detail::RewriteBudget& budget) {
        ExprPtr factor_real;
        ExprPtr factor_imag;
        if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(
                detail::node(factor))) {
            saw_complex = true;
            factor_real = detail::make_expression_ptr(complex->real());
            factor_imag = detail::make_expression_ptr(complex->imag());
        } else {
            factor_real = factor;
            budget.require_nodes(1);
            factor_imag = SymbolicExpr::number(0);
        }
        const auto ac_size = normalization_check_arithmetic<MultiplyNode>(
            &budget, 0, detail::node(real), detail::node(factor_real));
        auto ac = SymbolicExpr::multiply(real, factor_real);
        const auto bd_size = normalization_check_arithmetic<MultiplyNode>(
            &budget, 0, detail::node(imag), detail::node(factor_imag));
        auto bd = SymbolicExpr::multiply(imag, factor_imag);
        const auto ad_size = normalization_check_arithmetic<MultiplyNode>(
            &budget, 0, detail::node(real), detail::node(factor_imag));
        auto ad = SymbolicExpr::multiply(real, factor_imag);
        const auto bc_size = normalization_check_arithmetic<MultiplyNode>(
            &budget, 0, detail::node(imag), detail::node(factor_real));
        auto bc = SymbolicExpr::multiply(imag, factor_real);
        const auto negative_size = budget.append_size(1, bd_size);
        budget.append_size(budget.append_size(1, ac_size), negative_size);
        auto next_real = normalize_equivalence(SymbolicExpr::add(
            ac, SymbolicExpr::multiply(SymbolicExpr::number(-1), bd)), budget);
        budget.append_size(budget.append_size(1, ad_size), bc_size);
        auto next_imag = normalize_equivalence(SymbolicExpr::add(ad, bc), budget);
        real = std::move(next_real);
        imag = std::move(next_imag);
    }
};

ExprPtr canonicalize_product(const MultiplyNode& multiply,
                              const detail::SymbolicNodePtr& original,
                              detail::RewriteBudget& budget) {
    std::vector<ExprPtr> factors;
    bool changed = false;
    bool saw_complex = false;
    factors.reserve(multiply.operands().size());
    for (const auto& operand : multiply.operands()) {
        auto canonical = canonicalize_complex_product(
            *detail::make_expression_ptr(operand), budget);
        changed = changed || detail::node(canonical) != operand;
        saw_complex = saw_complex ||
            static_cast<bool>(std::dynamic_pointer_cast<const ComplexNode>(
                detail::node(canonical)));
        factors.push_back(std::move(canonical));
    }
    if (!changed && !saw_complex) {
        budget.measure(original);
        return detail::make_expression_ptr(original);
    }
    if (!saw_complex) {
        std::vector<detail::SymbolicNodePtr> operands;
        std::size_t size = 0;
        operands.reserve(factors.size());
        for (const auto& factor : factors) {
            size = budget.append_size(size, budget.measure(detail::node(factor)));
            operands.push_back(detail::node(factor));
        }
        return normalize_equivalence(detail::make_expression_ptr(
            SymbolicFactory::create_multiply(std::move(operands), &budget)), budget);
    }
    ComplexProduct product;
    for (const auto& factor : factors) {
        product.accumulate(factor, budget);
    }
    return rewritten_complex(product.real, product.imag, original, budget);
}

}

ExprPtr canonicalize_complex_product(const SymbolicExpr& expression,
                                     detail::RewriteBudget& budget) {
    detail::RewriteScope scope(&budget);
    const auto& node = detail::node(expression);
    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        budget.require_nodes(1);
        return canonicalize_variable(*variable, node, budget);
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return canonicalize_sum(*add, node, budget);
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return canonicalize_power(*power, node, budget);
    }
    if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        auto real = canonicalize_complex_product(
            *detail::make_expression_ptr(complex->real()), budget);
        auto imag = canonicalize_complex_product(
            *detail::make_expression_ptr(complex->imag()), budget);
        return rewritten_complex(real, imag, node, budget);
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return canonicalize_product(*multiply, node, budget);
    }
    budget.measure(node);
    return detail::make_expression_ptr(node);
}

}
