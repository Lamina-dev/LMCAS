/**
 * @file print_visitor.hpp
 * @brief 打印访问器，将 AST 转换为可读字符串。
 */
#pragma once
#include <sstream>
#include <string>
#include "../lamina_export.hpp"
#include "../symbolic_ast.hpp"

/** @brief 打印访问器，遍历 AST 并生成可读的数学表达式字符串 */
class LAMINA_API PrintVisitor : public lamina::detail::SymbolicVisitor {
public:
    PrintVisitor() = default;

    /**
     * @brief 获取打印结果
     * @return 格式化后的表达式字符串
     */
    std::string get_result() const {
        return buffer.str();
    }

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

private:
    std::stringstream buffer;
};
