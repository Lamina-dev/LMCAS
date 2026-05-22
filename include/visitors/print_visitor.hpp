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

class LAMINA_API PrintVisitor : public SymbolicVisitor {
public:
    PrintVisitor() = default;
    
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

private:
    std::stringstream buffer;
};
