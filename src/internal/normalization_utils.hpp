/** @file internal/normalization_utils.hpp */
#pragma once
#include "lmmc/config.h"
#include "lmmc/numeric.h"
#include "internal/symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/facts_query.hpp"
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

namespace LMCAS {

inline bool normalization_proved(const Result<Tribool>& result) {
    if (!result) { throw result.error(); }
    return result.value() == Tribool::True;
}
inline bool normalization_can_discard(
    const std::shared_ptr<const SymbolicNode>& node, const FactsQuery& facts,
    Domain domain, ComputationContext& context) {
    return (domain == Domain::Complex ||
            normalization_proved(detail::query_definedness(node, facts, Domain::Real, context))) &&
           normalization_proved(detail::query_definedness(node, facts, Domain::Complex, context));
}
void normalization_append(detail::RewriteBudget* budget, std::size_t& nodes,
    std::vector<std::shared_ptr<const SymbolicNode>>& children,
    std::shared_ptr<const SymbolicNode> child);
template <class... Nodes>
void normalization_check_children(detail::RewriteBudget* budget,
    std::size_t roots, const std::shared_ptr<const Nodes>&... children) {
    if (!budget) {
        return;
    }
    budget->require_nodes(roots);
    ((roots = budget->append_size(roots, budget->measure(children))), ...);
}

inline void normalization_check_children(detail::RewriteBudget* budget,
    std::size_t roots,
    const std::vector<std::shared_ptr<const SymbolicNode>>& children) {
    if (!budget) {
        return;
    }
    budget->require_nodes(roots);
    for (const auto& child : children) {
        roots = budget->append_size(roots, budget->measure(child));
    }
}
template <class Parent, class... Nodes>
std::size_t normalization_check_arithmetic(detail::RewriteBudget* budget,
    std::size_t extra_nodes, const std::shared_ptr<const Nodes>&... children) {
    if (!budget) {
        return 0;
    }
    std::size_t nodes = extra_nodes, flattened = 0;
    const auto append = [&](const auto& child) {
        nodes = budget->append_size(nodes, budget->measure(child));
        if (dynamic_cast<const Parent*>(child.get())) {
            ++flattened;
        }
    };
    (append(children), ...);
    return budget->append_size(1, nodes - flattened);
}

template <class Parent>
std::size_t normalization_check_arithmetic(detail::RewriteBudget* budget,
    const std::vector<std::shared_ptr<const SymbolicNode>>& children,
    std::size_t nodes, bool wrapped = true) {
    if (!budget) {
        return 0;
    }
    budget->require_nodes(nodes);
    if (!wrapped) {
        return nodes;
    }
    for (const auto& child : children) {
        if (dynamic_cast<const Parent*>(child.get())) {
            --nodes;
        }
    }
    return budget->append_size(1, nodes);
}
std::size_t normalization_check_count(detail::RewriteBudget* budget,
    std::size_t rows, std::size_t cols = 1, bool wrapped = true);


bool try_get_integer_value(const std::shared_ptr<const NumberNode>& node, BigInt& value);
BigInt get_node_degree_helper(const std::shared_ptr<const SymbolicNode>& node);
std::shared_ptr<const SymbolicNode> make_normalized_multiply_node(
    const std::vector<std::shared_ptr<const SymbolicNode>>& operands,
    detail::RewriteBudget* budget = nullptr);
std::shared_ptr<const SymbolicNode> norm_subst_index(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& index_var,
    const std::shared_ptr<const SymbolicNode>& value, detail::RewriteBudget* budget = nullptr);
struct NodeCompare {
    bool operator()(const std::shared_ptr<const SymbolicNode>& left,
                    const std::shared_ptr<const SymbolicNode>& right) const;
};
bool normalization_factor_less(const std::shared_ptr<const SymbolicNode>& left,
                               const std::shared_ptr<const SymbolicNode>& right);
lmmc_real_t normalization_number_value(const NumberNode& node);
bool is_approximate_number(const NumberNode& node);
bool contains_inexact_number(
    const std::shared_ptr<const SymbolicNode>& node);
Rational exact_number_as_rational(const NumberNode& node);
std::shared_ptr<const NumberNode> add_numbers(
    const std::shared_ptr<const NumberNode>& left, const std::shared_ptr<const NumberNode>& right);
std::shared_ptr<const NumberNode> multiply_numbers(
    const std::shared_ptr<const NumberNode>& left, const std::shared_ptr<const NumberNode>& right);
bool check_negative_arg(const std::shared_ptr<const SymbolicNode>& argument,
                        std::shared_ptr<const SymbolicNode>& positive,
                        detail::RewriteBudget* budget = nullptr);
bool get_pi_coeff(const std::shared_ptr<const SymbolicNode>& node, Rational& coefficient);
std::shared_ptr<const SymbolicNode> normalization_imaginary_power(
    const BigInt& exponent, detail::RewriteBudget* budget = nullptr);
std::shared_ptr<const NumberNode> normalization_exact_square_root(const NumberNode& number);
std::shared_ptr<const SymbolicNode> normalization_numeric_power(
    const std::shared_ptr<const NumberNode>& base,
    const std::shared_ptr<const NumberNode>& exponent, detail::RewriteBudget* budget = nullptr);
std::shared_ptr<const SymbolicNode> normalization_numeric_function(
    FunctionNode::FuncType type, const std::shared_ptr<const NumberNode>& number,
    detail::RewriteBudget* budget = nullptr);
std::shared_ptr<const SymbolicNode> normalization_pi_function(
    FunctionNode::FuncType type, const Rational& coefficient,
    detail::RewriteBudget* budget = nullptr);

}
