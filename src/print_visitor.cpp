#include "../include/visitors/print_visitor.hpp"
#include <iomanip>

void PrintVisitor::visit(NumberNode& node) {
    if (std::holds_alternative<BigInt>(node.value)) {
        buffer << std::get<BigInt>(node.value).to_string();
    } else if (std::holds_alternative<Rational>(node.value)) {
        buffer << std::get<Rational>(node.value).to_string();
    } else {
        buffer << std::get<double>(node.value);
    }
}

void PrintVisitor::visit(VariableNode& node) {
    buffer << node.name;
}

void PrintVisitor::visit(AddNode& node) {
    buffer << "(";
    for (size_t i = 0; i < node.operands.size(); ++i) {
        node.operands[i]->accept(*this);
        if (i < node.operands.size() - 1) {
            buffer << " + ";
        }
    }
    buffer << ")";
}

void PrintVisitor::visit(MultiplyNode& node) {
    buffer << "(";
    for (size_t i = 0; i < node.operands.size(); ++i) {
        node.operands[i]->accept(*this);
        if (i < node.operands.size() - 1) {
            buffer << " * ";
        }
    }
    buffer << ")";
}

void PrintVisitor::visit(PowerNode& node) {
    buffer << "(";
    node.base->accept(*this);
    buffer << " ^ ";
    node.exponent->accept(*this);
    buffer << ")";
}

void PrintVisitor::visit(FunctionNode& node) {
    switch (node.type) {
        case FunctionNode::FuncType::Sin: buffer << "sin"; break;
        case FunctionNode::FuncType::Cos: buffer << "cos"; break;
        case FunctionNode::FuncType::Tan: buffer << "tan"; break;
        case FunctionNode::FuncType::Cot: buffer << "cot"; break;
        case FunctionNode::FuncType::Sec: buffer << "sec"; break;
        case FunctionNode::FuncType::Csc: buffer << "csc"; break;
        case FunctionNode::FuncType::ArcSin: buffer << "asin"; break;
        case FunctionNode::FuncType::ArcCos: buffer << "acos"; break;
        case FunctionNode::FuncType::ArcTan: buffer << "atan"; break;
        case FunctionNode::FuncType::Sinh: buffer << "sinh"; break;
        case FunctionNode::FuncType::Cosh: buffer << "cosh"; break;
        case FunctionNode::FuncType::Tanh: buffer << "tanh"; break;
        case FunctionNode::FuncType::Ln: buffer << "ln"; break;
        case FunctionNode::FuncType::Log: buffer << "log"; break;
        case FunctionNode::FuncType::Abs: buffer << "abs"; break;
        case FunctionNode::FuncType::Sqrt: buffer << "sqrt"; break;
        case FunctionNode::FuncType::Exp: buffer << "exp"; break;
    }
    buffer << "(";
    for (size_t i = 0; i < node.arguments.size(); ++i) {
        node.arguments[i]->accept(*this);
        if (i < node.arguments.size() - 1) {
            buffer << ", ";
        }
    }
    buffer << ")";
}

void PrintVisitor::visit(MatrixNode& node) {
    buffer << "[";
    for (size_t i = 0; i < node.rows; ++i) {
        buffer << "[";
        for (size_t j = 0; j < node.cols; ++j) {
            auto element = node.get(i, j);
            if (element) {
                element->accept(*this);
            } else {
                buffer << "0"; // Fallback if no element (shouldn't happen with get())
            }
            if (j < node.cols - 1) {
                buffer << ", ";
            }
        }
        buffer << "]";
        if (i < node.rows - 1) {
            buffer << ", ";
        }
    }
    buffer << "]";
}

void PrintVisitor::visit(RelationalNode& node) {
    node.left->accept(*this);
    buffer << " " << RelationalNode::op_to_string(node.op) << " ";
    node.right->accept(*this);
}
