#pragma once
#include <sstream>
#include <string>
#include "../symbolic_ast.hpp"

class PrintVisitor : public SymbolicVisitor {
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
    
private:
    std::stringstream buffer;
};
