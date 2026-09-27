/**
 * @file normalization_visitor.hpp
 * @brief 符号规范化访问器接口。
 */
#pragma once

#include "internal/symbolic_ast.hpp"
#include "internal/facts_query.hpp"

#include <memory>
#include <optional>

namespace LMCAS {

/**
 * @brief 统一规范化数值幂。
 * 精确操作数保留整数和有理数运算；近似求根要求指数精确等于 1/2。
 * 有限浮点幂直接折叠，省去倒数中间值；非有限结果保留幂表达式。
 */
class LMCAS_API NormalizationVisitor final : public LMCAS::detail::SymbolicVisitor {
public:
    NormalizationVisitor() : NormalizationVisitor(nullptr) {}
    explicit NormalizationVisitor(detail::RewriteBudget* budget)
        : SymbolicVisitor(budget), local_context_(budget ? std::nullopt
              : std::optional<ComputationContext>(std::in_place)),
          context_(budget ? budget->context() : *local_context_) {}
    explicit NormalizationVisitor(ComputationContext& context,
                                  const FactsQuery& facts = detail::no_facts(),
                                  Domain domain = Domain::Real,
                                  detail::RewriteBudget* budget = nullptr)
        : SymbolicVisitor(budget), context_(context), facts_(facts), domain_(domain) {}
    explicit NormalizationVisitor(const FactsQuery& facts, Domain domain = Domain::Real,
                                  detail::RewriteBudget* budget = nullptr)
        : SymbolicVisitor(budget), local_context_(budget ? std::nullopt
              : std::optional<ComputationContext>(std::in_place)),
          context_(budget ? budget->context() : *local_context_), facts_(facts), domain_(domain) {}
    NormalizationVisitor(ComputationContext&, FactsQuery&&, Domain = Domain::Real,
                         detail::RewriteBudget* = nullptr) = delete;
    NormalizationVisitor(FactsQuery&&, Domain = Domain::Real,
                         detail::RewriteBudget* = nullptr) = delete;

    [[nodiscard]] std::shared_ptr<const SymbolicNode> get_result() const;

    std::shared_ptr<const SymbolicNode> expand_product(
        const std::shared_ptr<const SymbolicNode>& lhs,
        const std::shared_ptr<const SymbolicNode>& rhs);

    void visit(const NumberNode& node) override;
    void visit(const VariableNode& node) override;
    void visit(const AddNode& node) override;
    void visit(const MultiplyNode& node) override;
    void visit(const PowerNode& node) override;
    void visit(const FunctionNode& node) override;
    void visit(const UninterpretedFunctionNode& node) override;
    void visit(const MatrixNode& node) override;
    void visit(const RelationalNode& node) override;
    void visit(const LogicalNode& node) override;
    void visit(const PiecewiseNode& node) override;
    void visit(const SummationNode& node) override;
    void visit(const ProductNode& node) override;
    void visit(const TransformNode& node) override;
    void visit(const QuantifierNode& node) override;
    void visit(const SetBuilderNode& node) override;
    void visit(const FiniteSetNode& node) override;
    void visit(const IntervalNode& node) override;
    void visit(const MembershipNode& node) override;
    void visit(const QuantityNode& node) override;
    void visit(const ComplexNode& node) override;
    void visit(const IntegralNode& node) override;
    void visit(const LimitNode& node) override;
    void visit(const RootOfNode& node) override;

private:
    std::optional<ComputationContext> local_context_;
    ComputationContext& context_;
    std::shared_ptr<const SymbolicNode> result;
    const FactsQuery& facts_ = detail::no_facts();
    Domain domain_ = Domain::Real;

    void set_result(std::shared_ptr<const SymbolicNode> candidate);

    std::shared_ptr<const SymbolicNode> normalize_bound_body(const SymbolicNode& node);

    void merge_complex_terms(std::vector<std::shared_ptr<const SymbolicNode>>& operands);
    std::shared_ptr<const NumberNode> extract_term_coefficient(
        std::shared_ptr<const SymbolicNode>& term);
    std::shared_ptr<const SymbolicNode> collect_sum(
        const std::vector<std::shared_ptr<const SymbolicNode>>& operands);
    std::shared_ptr<const SymbolicNode> collect_product(
        const std::vector<std::shared_ptr<const SymbolicNode>>& operands);
    std::shared_ptr<const SymbolicNode> collect_expanded_product(
        const std::shared_ptr<const SymbolicNode>& lhs,
        const std::shared_ptr<const SymbolicNode>& rhs);
    std::shared_ptr<const SymbolicNode> expand_complex_product(
        const std::shared_ptr<const SymbolicNode>& lhs,
        const std::shared_ptr<const SymbolicNode>& rhs);
    std::shared_ptr<const SymbolicNode> distribute_product(
        const std::shared_ptr<const SymbolicNode>& lhs,
        const std::shared_ptr<const SymbolicNode>& rhs);
    bool normalize_matrix_sum(
        const std::vector<std::shared_ptr<const SymbolicNode>>& operands);
    std::shared_ptr<const SymbolicNode> normalize_matrix_product(
        const std::vector<std::shared_ptr<const SymbolicNode>>& operands);
    std::shared_ptr<const SymbolicNode> fuse_matrix_product(
        const std::shared_ptr<const SymbolicNode>& left,
        const std::shared_ptr<const SymbolicNode>& right);
    std::shared_ptr<const SymbolicNode> matrix_product_element(
        const MatrixNode& left, const MatrixNode& right, std::size_t row, std::size_t col);
    std::shared_ptr<const SymbolicNode> scale_matrix(
        const MatrixNode& matrix, const std::shared_ptr<const NumberNode>& scalar);
    bool normalize_power_identity(
        const std::shared_ptr<const SymbolicNode>& base,
        const std::shared_ptr<const SymbolicNode>& exponent);
    bool normalize_product_power(
        const MultiplyNode& base, const std::shared_ptr<const SymbolicNode>& exponent);
    void distribute_power(
        const MultiplyNode& base, const std::shared_ptr<const SymbolicNode>& exponent);
    bool normalize_nested_imaginary_power(
        const PowerNode& base, const std::shared_ptr<const NumberNode>& exponent);
    bool normalize_nested_power(
        const PowerNode& base, const std::shared_ptr<const SymbolicNode>& exponent);
    bool finish_squared_norm(const SymbolicNode& node,
        std::shared_ptr<const SymbolicNode>& argument, bool normalize_argument);
    bool normalize_logarithm(FunctionNode::FuncType type,
        const std::vector<std::shared_ptr<const SymbolicNode>>& arguments);
    bool normalize_natural_logarithm(const std::shared_ptr<const SymbolicNode>& argument);
    bool normalize_function_parity(FunctionNode::FuncType type,
        const std::shared_ptr<const SymbolicNode>& argument);
    bool normalize_unary_function(FunctionNode::FuncType type,
        const std::shared_ptr<const SymbolicNode>& argument);

    bool try_normalize_squared_norm(
        const SymbolicNode& node, std::shared_ptr<const SymbolicNode>& argument);
    static bool is_positive_integer_number(
        const std::shared_ptr<const NumberNode>& node);
    bool can_combine_nested_power(
        const PowerNode& inner,
        const std::shared_ptr<const NumberNode>& outer_exponent) const;
    std::shared_ptr<const SymbolicNode> try_assumption_simplify(
        const std::shared_ptr<const SymbolicNode>& node);
    std::shared_ptr<const SymbolicNode> simplify_assumed_square_root(
        const std::shared_ptr<const SymbolicNode>& argument) const;
    std::shared_ptr<const SymbolicNode> simplify_assumed_absolute(
        const std::shared_ptr<const SymbolicNode>& argument) const;
};

}
