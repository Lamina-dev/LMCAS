/**
 * @file normalization_visitor.hpp
 * @brief Symbolic normalization visitor interface.
 */
#pragma once

#include "../symbolic_ast.hpp"

#include <memory>

namespace lamina {
class AssumptionContext;
}

class LAMINA_API NormalizationVisitor final : public lamina::detail::SymbolicVisitor {
public:
    NormalizationVisitor() = default;
    explicit NormalizationVisitor(const lamina::AssumptionContext* assumptions)
        : assumptions_(assumptions) {}

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
    std::shared_ptr<const SymbolicNode> result;
    const lamina::AssumptionContext* assumptions_ = nullptr;

    static bool try_get_integer_value(
        const std::shared_ptr<const NumberNode>& node, long long& value);
    static bool is_positive_integer_number(
        const std::shared_ptr<const NumberNode>& node);
    bool is_provably_nonzero(
        const std::shared_ptr<const SymbolicNode>& node) const;
    std::shared_ptr<const SymbolicNode> try_assumption_simplify(
        const std::shared_ptr<const SymbolicNode>& node);
};
