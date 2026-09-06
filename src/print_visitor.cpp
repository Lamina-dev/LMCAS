#include "../include/visitors/print_visitor.hpp"
#include <cmath>
#include <iomanip>

namespace LMCAS {

void PrintVisitor::visit(const NumberNode& node) {
    if (std::holds_alternative<BigInt>(node.value())) {
        buffer << std::get<BigInt>(node.value()).to_string();
    } else if (std::holds_alternative<Rational>(node.value())) {
        buffer << std::get<Rational>(node.value()).to_string();
    } else {
        buffer << std::get<double>(node.value());
    }
}

void PrintVisitor::visit(const VariableNode& node) {
    buffer << node.name();
}

void PrintVisitor::visit(const AddNode& node) {
    if (node.operands().empty()) {
        buffer << "0";
        return;
    }
    for (size_t i = 0; i < node.operands().size(); ++i) {
        if (i > 0) {
            std::ostringstream sub;
            PrintVisitor sub_v;
            node.operands()[i]->accept(sub_v);
            std::string s = sub_v.get_result();
            if (!s.empty() && s[0] == '-') {
                buffer << " - " << s.substr(1);
            } else {
                buffer << " + " << s;
            }
        } else {
            node.operands()[i]->accept(*this);
        }
    }
}

void PrintVisitor::visit(const MultiplyNode& node) {
    if (node.operands().empty()) {
        buffer << "1";
        return;
    }
    for (size_t i = 0; i < node.operands().size(); ++i) {
        if (i > 0) buffer << "*";

        bool needs_parens = std::dynamic_pointer_cast<const AddNode>(node.operands()[i]) != nullptr ||
                           std::dynamic_pointer_cast<const PowerNode>(node.operands()[i]) != nullptr;

        if (!needs_parens) {
            if (auto num = std::dynamic_pointer_cast<const NumberNode>(node.operands()[i])) {
                if (std::holds_alternative<Rational>(num->value())) {
                    needs_parens = true;
                }
            }
        }

        if (needs_parens) buffer << "(";
        node.operands()[i]->accept(*this);
        if (needs_parens) buffer << ")";
    }
}

void PrintVisitor::visit(const PowerNode& node) {
    bool base_parens = dynamic_cast<const AddNode*>(node.base().get()) ||
                       dynamic_cast<const MultiplyNode*>(node.base().get()) ||
                       dynamic_cast<const PowerNode*>(node.base().get());
    if (const auto* number = dynamic_cast<const NumberNode*>(node.base().get())) {
        const auto& value = number->value();
        base_parens = std::holds_alternative<Rational>(value) ||
                      (std::holds_alternative<BigInt>(value) &&
                       std::get<BigInt>(value).IsNegative()) ||
                      (std::holds_alternative<double>(value) &&
                       std::signbit(std::get<double>(value)));
    }
    if (base_parens) buffer << "(";
    node.base()->accept(*this);
    if (base_parens) buffer << ")";
    buffer << "^";
    bool exp_parens = dynamic_cast<const AddNode*>(node.exponent().get()) ||
                      dynamic_cast<const MultiplyNode*>(node.exponent().get()) ||
                      dynamic_cast<const PowerNode*>(node.exponent().get());
    if (const auto* number = dynamic_cast<const NumberNode*>(node.exponent().get())) {
        const auto* rational = std::get_if<Rational>(&number->value());
        exp_parens = rational && !rational->is_integer();
    }
    if (exp_parens) buffer << "(";
    node.exponent()->accept(*this);
    if (exp_parens) buffer << ")";
}

void PrintVisitor::visit(const FunctionNode& node) {
    if (node.type() == FunctionNode::FuncType::Infinity) {
        buffer << "inf";
        return;
    }

    switch (node.type()) {
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
        case FunctionNode::FuncType::Infinity: buffer << "inf"; break;
        case FunctionNode::FuncType::Erf: buffer << "erf"; break;
        case FunctionNode::FuncType::Ei: buffer << "Ei"; break;
        case FunctionNode::FuncType::Si: buffer << "Si"; break;
        case FunctionNode::FuncType::Ci: buffer << "Ci"; break;
        case FunctionNode::FuncType::Li: buffer << "Li"; break;
        case FunctionNode::FuncType::Max: buffer << "max"; break;
        case FunctionNode::FuncType::Min: buffer << "min"; break;
        case FunctionNode::FuncType::Sgn: buffer << "sgn"; break;
        case FunctionNode::FuncType::Floor: buffer << "floor"; break;
        case FunctionNode::FuncType::Ceil: buffer << "ceil"; break;
        case FunctionNode::FuncType::Round: buffer << "round"; break;
        case FunctionNode::FuncType::RealPart: buffer << "re"; break;
        case FunctionNode::FuncType::ImagPart: buffer << "im"; break;
        case FunctionNode::FuncType::Conjugate: buffer << "conj"; break;
        case FunctionNode::FuncType::ComplexAbs: buffer << "cabs"; break;
        case FunctionNode::FuncType::ComplexArg: buffer << "carg"; break;
    }
    buffer << "(";
    for (size_t i = 0; i < node.arguments().size(); ++i) {
        node.arguments()[i]->accept(*this);
        if (i < node.arguments().size() - 1) {
            buffer << ", ";
        }
    }
    buffer << ")";
}

void PrintVisitor::visit(const UninterpretedFunctionNode& node) {
    buffer << node.name() << "(";
    for (std::size_t index = 0; index < node.arguments().size(); ++index) {
        if (index != 0) buffer << ", ";
        node.arguments()[index]->accept(*this);
    }
    buffer << ")";
}

void PrintVisitor::visit(const MatrixNode& node) {
    buffer << "[";
    for (size_t i = 0; i < node.rows(); ++i) {
        buffer << "[";
        for (size_t j = 0; j < node.cols(); ++j) {
            auto element = node.get(i, j);
            if (element) {
                element->accept(*this);
            } else {
                buffer << "0";
            }
            if (j < node.cols() - 1) {
                buffer << ", ";
            }
        }
        buffer << "]";
        if (i < node.rows() - 1) {
            buffer << ", ";
        }
    }
    buffer << "]";
}

void PrintVisitor::visit(const RelationalNode& node) {
    node.left()->accept(*this);
    buffer << " " << RelationalNode::op_to_string(node.op()) << " ";
    node.right()->accept(*this);
}

void PrintVisitor::visit(const LogicalNode& node) {
    if (node.op() == LogicalNode::Op::Not) {
        buffer << "(not ";
        node.left()->accept(*this);
        buffer << ")";
    } else if (node.op() == LogicalNode::Op::Implies) {
        buffer << "(";
        node.left()->accept(*this);
        buffer << " implies ";
        node.right()->accept(*this);
        buffer << ")";
    } else {
        buffer << "(";
        node.left()->accept(*this);
        buffer << " " << (node.op() == LogicalNode::Op::And ? "and" : "or") << " ";
        node.right()->accept(*this);
        buffer << ")";
    }
}

void PrintVisitor::visit(const PiecewiseNode& node) {
    buffer << "piecewise(";
    for (size_t i = 0; i < node.branches().size(); ++i) {
        if (i > 0) buffer << ", ";
        node.branches()[i].expression->accept(*this);
        buffer << " if ";
        node.branches()[i].condition->accept(*this);
    }
    if (node.default_expr()) {
        buffer << ", default: ";
        node.default_expr()->accept(*this);
    }
    buffer << ")";
}

void PrintVisitor::visit(const SummationNode& node) {
    buffer << "\xce\xa3(";
    node.body()->accept(*this);
    buffer << ", " << node.index_var() << "=";
    node.lower_bound()->accept(*this);
    buffer << "..";
    node.upper_bound()->accept(*this);
    buffer << ")";
}

void PrintVisitor::visit(const ProductNode& node) {
    buffer << "\xce\xa0(";
    node.body()->accept(*this);
    buffer << ", " << node.index_var() << "=";
    node.lower_bound()->accept(*this);
    buffer << "..";
    node.upper_bound()->accept(*this);
    buffer << ")";
}

void PrintVisitor::visit(const TransformNode& node) {
    switch (node.transform_type()) {
        case TransformNode::TransformType::Laplace:
            buffer << "L{";
            node.body()->accept(*this);
            buffer << "}(";
            node.target()->accept(*this);
            buffer << ")";
            break;
        case TransformNode::TransformType::InverseLaplace:
            buffer << "L\xe2\x81\xbb\xc2\xb9{";
            node.body()->accept(*this);
            buffer << "}(";
            node.target()->accept(*this);
            buffer << ")";
            break;
        case TransformNode::TransformType::Fourier:
            buffer << "F{";
            node.body()->accept(*this);
            buffer << "}(";
            node.target()->accept(*this);
            buffer << ")";
            break;
        case TransformNode::TransformType::InverseFourier:
            buffer << "F\xe2\x81\xbb\xc2\xb9{";
            node.body()->accept(*this);
            buffer << "}(";
            node.target()->accept(*this);
            buffer << ")";
            break;
        case TransformNode::TransformType::ZTransform:
            buffer << "Z{";
            node.body()->accept(*this);
            buffer << "}(";
            node.target()->accept(*this);
            buffer << ")";
            break;
    }
}

void PrintVisitor::visit(const QuantifierNode& node) {
    if (node.quantifier_type() == QuantifierNode::Type::ForAll) {
        buffer << "\xe2\x88\x80";
    } else {
        buffer << "\xe2\x88\x83";
    }
    buffer << node.bound_var() << "\xe2\x88\x88";
    node.domain()->accept(*this);
    buffer << ": ";
    node.predicate()->accept(*this);
}

void PrintVisitor::visit(const SetBuilderNode& node) {
    buffer << "{" << node.element_var() << " \xe2\x88\x88 ";
    node.domain()->accept(*this);
    buffer << " | ";
    node.predicate()->accept(*this);
    buffer << "}";
}

void PrintVisitor::visit(const FiniteSetNode& node) {
    buffer << "{";
    for (std::size_t index = 0; index < node.elements().size(); ++index) {
        if (index != 0) buffer << ", ";
        node.elements()[index]->accept(*this);
    }
    buffer << "}";
}

void PrintVisitor::visit(const IntervalNode& node) {
    buffer << (node.lower_closed() ? "[" : "(");
    node.lower()->accept(*this);
    buffer << ", ";
    node.upper()->accept(*this);
    buffer << (node.upper_closed() ? "]" : ")");
}

void PrintVisitor::visit(const MembershipNode& node) {
    node.element()->accept(*this);
    buffer << " in ";
    node.set()->accept(*this);
}

void PrintVisitor::visit(const QuantityNode& node) {
    node.value()->accept(*this);
    buffer << "<";
    buffer << (node.display_unit().empty()
                   ? node.dimension().to_string()
                   : node.display_unit());
    buffer << ">";
}

void PrintVisitor::visit(const ComplexNode& node) {
    buffer << "(";
    node.real()->accept(*this);
    buffer << " + ";
    node.imag()->accept(*this);
    buffer << "*I)";
}

void PrintVisitor::visit(const IntegralNode& node) {
    buffer << "Integral(";
    node.body()->accept(*this);
    buffer << ", " << node.variable();
    if (node.is_definite()) {
        buffer << ", ";
        node.lower()->accept(*this);
        buffer << ", ";
        node.upper()->accept(*this);
    }
    buffer << ")";
}

void PrintVisitor::visit(const LimitNode& node) {
    buffer << "limit(";
    node.body()->accept(*this);
    buffer << ", " << node.variable() << ", ";
    node.point()->accept(*this);
    switch (node.direction()) {
        case LimitDirection::Both: buffer << ", both)"; break;
        case LimitDirection::FromBelow: buffer << ", below)"; break;
        case LimitDirection::FromAbove: buffer << ", above)"; break;
    }
}

void PrintVisitor::visit(const RootOfNode& node) {
    buffer << "rootof(";
    node.polynomial()->accept(*this);
    buffer << ", " << node.variable() << ", " << node.index() << ")";
}

} // namespace LMCAS
