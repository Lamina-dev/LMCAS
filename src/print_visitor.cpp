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
    if (node.operands.empty()) {
        buffer << "0";
        return;
    }
    for (size_t i = 0; i < node.operands.size(); ++i) {
        if (i > 0) {
            std::ostringstream sub;
            PrintVisitor sub_v;
            node.operands[i]->accept(sub_v);
            std::string s = sub_v.get_result();
            if (!s.empty() && s[0] == '-') {
                buffer << " - " << s.substr(1);
            } else {
                buffer << " + " << s;
            }
        } else {
            node.operands[i]->accept(*this);
        }
    }
}

void PrintVisitor::visit(MultiplyNode& node) {
    if (node.operands.empty()) {
        buffer << "1";
        return;
    }
    for (size_t i = 0; i < node.operands.size(); ++i) {
        if (i > 0) buffer << "*";

        bool needs_parens = std::dynamic_pointer_cast<AddNode>(node.operands[i]) != nullptr ||
                           std::dynamic_pointer_cast<PowerNode>(node.operands[i]) != nullptr;

        if (!needs_parens) {
            if (auto num = std::dynamic_pointer_cast<NumberNode>(node.operands[i])) {
                if (std::holds_alternative<Rational>(num->value)) {
                    needs_parens = true;
                }
            }
        }

        if (needs_parens) buffer << "(";
        node.operands[i]->accept(*this);
        if (needs_parens) buffer << ")";
    }
}

void PrintVisitor::visit(PowerNode& node) {
    bool base_parens = std::dynamic_pointer_cast<AddNode>(node.base) ||
                      std::dynamic_pointer_cast<MultiplyNode>(node.base);
    if (base_parens) buffer << "(";
    node.base->accept(*this);
    if (base_parens) buffer << ")";
    buffer << "^";
    bool exp_parens = std::dynamic_pointer_cast<AddNode>(node.exponent) ||
                     std::dynamic_pointer_cast<MultiplyNode>(node.exponent) ||
                     std::dynamic_pointer_cast<PowerNode>(node.exponent);
    if (exp_parens) buffer << "(";
    node.exponent->accept(*this);
    if (exp_parens) buffer << ")";
}

void PrintVisitor::visit(FunctionNode& node) {
    if (node.type == FunctionNode::FuncType::Infinity) {
        buffer << "inf";
        return;
    }

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
        case FunctionNode::FuncType::Atan2: buffer << "atan2"; break;
        case FunctionNode::FuncType::Sinh: buffer << "sinh"; break;
        case FunctionNode::FuncType::Cosh: buffer << "cosh"; break;
        case FunctionNode::FuncType::Tanh: buffer << "tanh"; break;
        case FunctionNode::FuncType::Ln: buffer << "ln"; break;
        case FunctionNode::FuncType::Log: buffer << "log"; break;
        case FunctionNode::FuncType::Abs: buffer << "abs"; break;
        case FunctionNode::FuncType::Sqrt: buffer << "sqrt"; break;
        case FunctionNode::FuncType::Exp: buffer << "exp"; break;
        case FunctionNode::FuncType::LambertW: buffer << "lambertw"; break;
        case FunctionNode::FuncType::RootOf: buffer << "rootof"; break;
        case FunctionNode::FuncType::Calculus_Integral: buffer << "integral"; break;
        case FunctionNode::FuncType::Limit: buffer << "limit"; break;
        case FunctionNode::FuncType::Infinity: buffer << "inf"; break;
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
                buffer << "0";
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

void PrintVisitor::visit(LogicalNode& node) {
    buffer << "(";
    node.left->accept(*this);
    buffer << " " << LogicalNode::op_to_string(node.op) << " ";
    node.right->accept(*this);
    buffer << ")";
}
