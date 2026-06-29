/**
 * @file print_visitor.hpp
 * @brief 打印访问器，将 AST 转换为可读字符串。
 */
#pragma once
#include <sstream>
#include <string>
#include "../symbolic_ast.hpp"

#ifdef _WIN32
#ifdef LAMINA_CORE_EXPORTS
#define LAMINA_API __declspec(dllexport)
#else
#define LAMINA_API __declspec(dllimport)
#endif
#else
#define LAMINA_API
#endif

/** @brief 打印访问器，遍历 AST 并生成可读的数学表达式字符串 */
class LAMINA_API PrintVisitor : public SymbolicVisitor {
public:
    PrintVisitor() = default;

    /**
     * @brief 获取打印结果
     * @return 格式化后的表达式字符串
     */
    std::string get_result() const {
        return buffer.str();
    }

    void visit(NumberNode& node) override;
    void visit(VariableNode& node) override;
    void visit(AddNode& node) override;
    void visit(MultiplyNode& node) override;
    void visit(PowerNode& node) override;
    void visit(FunctionNode& node) override;
    void visit(MatrixNode& node) override;
    void visit(RelationalNode& node) override;
    void visit(LogicalNode& node) override;
    void visit(PiecewiseNode& node) override;
    void visit(SummationNode& node) override;
    void visit(ProductNode_Op& node) override;
    void visit(TransformNode& node) override;
    void visit(QuantifierNode& node) override;
    void visit(SetBuilderNode& node) override;
    void visit(ComplexNode& node) override;

private:
    std::stringstream buffer;
};
