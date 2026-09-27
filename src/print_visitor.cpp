#include "internal/visitors/print_visitor.hpp"
#include <cmath>
#include <charconv>
#include <limits>
#include <stdexcept>
#include <array>

namespace LMCAS {

PrintVisitor::Precedence PrintVisitor::precedence(const SymbolicNode& node) {
    if (dynamic_cast<const AddNode*>(&node)) { return Precedence::Add; }
    if (dynamic_cast<const MultiplyNode*>(&node)) { return Precedence::Multiply; }
    if (dynamic_cast<const PowerNode*>(&node)) { return Precedence::Power; }
    if (dynamic_cast<const MembershipNode*>(&node)) { return Precedence::Membership; }
    if (const auto* relation = dynamic_cast<const RelationalNode*>(&node)) {
        return relation->op() == RelationOp::EQ || relation->op() == RelationOp::NEQ
            ? Precedence::Equality : Precedence::Relation;
    }
    if (const auto* number = dynamic_cast<const NumberNode*>(&node)) {
        if (const auto* rational = std::get_if<Rational>(&number->value())) {
            if (!rational->is_integer()) { return Precedence::Multiply; }
            if (*rational < Rational(0)) { return Precedence::Unary; }
        }
        if (const auto* integer = std::get_if<BigInt>(&number->value())) {
            if (integer->is_negative()) { return Precedence::Unary; }
        }
    }
    return Precedence::Atom;
}

void PrintVisitor::print_operand(const SymbolicNode& node, Precedence parent,
                                bool parenthesize_equal) {
    const auto child = precedence(node);
    const bool parentheses = child < parent ||
        (parenthesize_equal && child == parent) ||
        (parent == Precedence::Power && parenthesize_equal && child == Precedence::Unary);
    if (parentheses) { buffer << "("; }
    node.accept(*this);
    if (parentheses) { buffer << ")"; }
}

void PrintVisitor::visit(const NumberNode& node) {
    if (std::holds_alternative<BigInt>(node.value())) {
        buffer << std::get<BigInt>(node.value()).to_string();
    } else if (std::holds_alternative<Rational>(node.value())) {
        buffer << std::get<Rational>(node.value()).to_string();
    } else {
        const double value = std::get<double>(node.value());
        if (!std::isfinite(value)) {
            buffer << value;
            return;
        }
        char token[64];
        const auto result = std::to_chars(token, token + sizeof(token), value,
            std::chars_format::general, std::numeric_limits<double>::max_digits10);
        if (result.ec != std::errc{}) {
            throw std::logic_error("finite double formatting failed");
        }
        buffer << "approx(";
        if (value == 0 && std::signbit(value)) buffer << "-0";
        else buffer.write(token, result.ptr - token);
        buffer << ")";
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
        if (i > 0) buffer << " + ";
        print_operand(*node.operands()[i], Precedence::Add, true);
    }
}

void PrintVisitor::visit(const MultiplyNode& node) {
    if (node.operands().empty()) {
        buffer << "1";
        return;
    }
    for (size_t i = 0; i < node.operands().size(); ++i) {
        if (i > 0) buffer << "*";

        print_operand(*node.operands()[i], Precedence::Multiply, true);
    }
}

namespace {

const char* function_name(FunctionNode::FuncType type) {
    using Type = FunctionNode::FuncType;
    static constexpr auto names = [] {
        std::array<const char*, static_cast<std::size_t>(Type::ComplexArg) + 1> values{};
        values[static_cast<std::size_t>(Type::Sin)] = "sin";
        values[static_cast<std::size_t>(Type::Cos)] = "cos";
        values[static_cast<std::size_t>(Type::Tan)] = "tan";
        values[static_cast<std::size_t>(Type::Cot)] = "cot";
        values[static_cast<std::size_t>(Type::Sec)] = "sec";
        values[static_cast<std::size_t>(Type::Csc)] = "csc";
        values[static_cast<std::size_t>(Type::ArcSin)] = "asin";
        values[static_cast<std::size_t>(Type::ArcCos)] = "acos";
        values[static_cast<std::size_t>(Type::ArcTan)] = "atan";
        values[static_cast<std::size_t>(Type::Atan2)] = "atan2";
        values[static_cast<std::size_t>(Type::Sinh)] = "sinh";
        values[static_cast<std::size_t>(Type::Cosh)] = "cosh";
        values[static_cast<std::size_t>(Type::Tanh)] = "tanh";
        values[static_cast<std::size_t>(Type::Ln)] = "ln";
        values[static_cast<std::size_t>(Type::Log)] = "log";
        values[static_cast<std::size_t>(Type::Abs)] = "abs";
        values[static_cast<std::size_t>(Type::Sqrt)] = "sqrt";
        values[static_cast<std::size_t>(Type::Exp)] = "exp";
        values[static_cast<std::size_t>(Type::LambertW)] = "lambertw";
        values[static_cast<std::size_t>(Type::Infinity)] = "inf";
        values[static_cast<std::size_t>(Type::Erf)] = "erf";
        values[static_cast<std::size_t>(Type::Ei)] = "Ei";
        values[static_cast<std::size_t>(Type::Si)] = "Si";
        values[static_cast<std::size_t>(Type::Ci)] = "Ci";
        values[static_cast<std::size_t>(Type::Li)] = "Li";
        values[static_cast<std::size_t>(Type::Max)] = "max";
        values[static_cast<std::size_t>(Type::Min)] = "min";
        values[static_cast<std::size_t>(Type::Sgn)] = "sgn";
        values[static_cast<std::size_t>(Type::Floor)] = "floor";
        values[static_cast<std::size_t>(Type::Ceil)] = "ceil";
        values[static_cast<std::size_t>(Type::Round)] = "round";
        values[static_cast<std::size_t>(Type::RealPart)] = "re";
        values[static_cast<std::size_t>(Type::ImagPart)] = "im";
        values[static_cast<std::size_t>(Type::Conjugate)] = "conj";
        values[static_cast<std::size_t>(Type::ComplexAbs)] = "cabs";
        values[static_cast<std::size_t>(Type::ComplexArg)] = "carg";
        return values;
    }();
    const auto index = static_cast<std::size_t>(type);
    return index < names.size() && names[index] ? names[index] : "";
}

}

void PrintVisitor::visit(const PowerNode& node) {
    print_operand(*node.base(), Precedence::Power, true);
    buffer << "^";
    print_operand(*node.exponent(), Precedence::Power);
}

void PrintVisitor::visit(const FunctionNode& node) {
    if (node.type() == FunctionNode::FuncType::Infinity) {
        buffer << "inf";
        return;
    }

    buffer << function_name(node.type());
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
    print_operand(*node.left(), precedence(node), true);
    buffer << " " << (node.op() == RelationOp::EQ ? "==" :
        RelationalNode::op_to_string(node.op())) << " ";
    print_operand(*node.right(), precedence(node), true);
}

void PrintVisitor::visit(const LogicalNode& node) {
    if (node.op() == LogicalNode::Op::Not) {
        buffer << "(not ";
        print_operand(*node.left(), Precedence::Not, true);
        buffer << ")";
    } else if (node.op() == LogicalNode::Op::Implies) {
        buffer << "(";
        node.left()->accept(*this);
        buffer << " implies ";
        node.right()->accept(*this);
        buffer << ")";
    } else {
        buffer << "(";
        print_operand(*node.left(), node.op() == LogicalNode::Op::And
            ? Precedence::And : Precedence::Or, true);
        buffer << " " << (node.op() == LogicalNode::Op::And ? "and" : "or") << " ";
        print_operand(*node.right(), node.op() == LogicalNode::Op::And
            ? Precedence::And : Precedence::Or, true);
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
    print_operand(*node.element(), Precedence::Membership, true);
    buffer << " in ";
    print_operand(*node.set(), Precedence::Membership, true);
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
    print_operand(*node.real(), Precedence::Add, true);
    buffer << " + ";
    print_operand(*node.imag(), Precedence::Multiply, true);
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
