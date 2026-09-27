/**
 * @file differentiation_visitor.hpp
 * @brief 微分访问器，对 AST 执行符号求导。
 */
#pragma once
#include "internal/symbolic_ast.hpp"
#include <stdexcept>
#include <string>

namespace LMCAS {

namespace detail {
class UnsupportedDifferentiation : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
}
class DifferentiationVisitor : public LMCAS::detail::SymbolicVisitor {
    std::string var;

    [[noreturn]] void unsupported(const char* node_type) const {
        throw detail::UnsupportedDifferentiation(std::string(node_type) +
                                 " differentiation is outside the current supported domain");
    }

public:
    std::shared_ptr<const SymbolicNode> result;  /**< 求导结果节点。 */


    /**
     * @brief 构造普通求导访问器。
     * @param v 求导变量名
     */
    DifferentiationVisitor(const std::string& v) : var(v), result(nullptr) {}


    /**
     * @brief 获取求导结果。
     * @return 求导后的 AST 节点
     */
    std::shared_ptr<const SymbolicNode> get_result() const {
        return result;
    }
    void visit(const NumberNode&) override;
    void visit(const VariableNode& node) override;
    void visit(const AddNode& node) override;
    void visit(const MultiplyNode& node) override;
    void visit(const PowerNode& node) override;
    void visit(const FunctionNode& node) override;
    void visit(const MatrixNode& node) override;
    void visit(const RelationalNode&) override;
    void visit(const LogicalNode&) override;
    void visit(const ComplexNode& node) override;
    void visit(const PiecewiseNode& node) override;
    void visit(const UninterpretedFunctionNode&) override;
    void visit(const SummationNode&) override;
    void visit(const ProductNode&) override;
    void visit(const TransformNode&) override;
    void visit(const QuantifierNode&) override;
    void visit(const SetBuilderNode&) override;
    void visit(const FiniteSetNode&) override;
    void visit(const IntervalNode&) override;
    void visit(const MembershipNode&) override;
    void visit(const QuantityNode& node) override;
    void visit(const IntegralNode& node) override;
    void visit(const LimitNode&) override;
    void visit(const RootOfNode&) override;

private:
    void differentiate_atan2(const FunctionNode& node);
};

}
