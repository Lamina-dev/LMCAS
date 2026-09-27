#include "internal/solver_support.hpp"
#include <cmath>

namespace LMCAS::solver_detail {

bool solver_number_value(const NumberNode& number, double& out_value) {
    if (std::holds_alternative<BigInt>(number.value())) {
        out_value = std::get<BigInt>(number.value()).to_double();
        return true;
    }
    if (std::holds_alternative<Rational>(number.value())) {
        out_value = std::get<Rational>(number.value()).to_double();
        return true;
    }
    if (std::holds_alternative<lmmc_real_t>(number.value())) {
        out_value = std::get<lmmc_real_t>(number.value());
        return true;
    }
    return false;
}

namespace {

bool square_root_is_imaginary(const FunctionNode& function) {
    if (function.type() != FunctionNode::FuncType::Sqrt || function.arguments().size() != 1) {
        return false;
    }
    auto number = std::dynamic_pointer_cast<const NumberNode>(function.arguments()[0]);
    if (!number) {
        return false;
    }
    const auto& value = number->value();
    if (const auto* integer = std::get_if<BigInt>(&value)) {
        return *integer == BigInt(-1) || integer->is_negative();
    }
    if (const auto* rational = std::get_if<Rational>(&value)) {
        return *rational == Rational(-1) || *rational < Rational(0);
    }
    if (const auto* real = std::get_if<lmmc_real_t>(&value)) {
        int equal;
        lmmc_double_nearly_equal(*real, -1.0, &equal);
        return equal || (std::isfinite(*real) && *real < 0.0);
    }
    return false;
}

bool power_is_imaginary(const PowerNode& power) {
    auto base = std::dynamic_pointer_cast<const NumberNode>(power.base());
    auto exponent = std::dynamic_pointer_cast<const NumberNode>(power.exponent());
    if (!base || !exponent) {
        return false;
    }
    double base_value = 0.0;
    bool base_numeric = solver_number_value(*base, base_value);
    double exponent_value = 0.0;
    bool exponent_numeric = solver_number_value(*exponent, exponent_value);
    if (!base_numeric || !exponent_numeric || base_value >= 0.0) {
        return false;
    }
    return base_value < 0.0 && std::abs(exponent_value - std::round(exponent_value)) > 1e-12;
}

/**
 * @brief 检查表达式树中的虚数分量。
 * @note 虚数单位表示为 sqrt(-1)，即 Sqrt 类型且唯一参数为数值 -1 的 FunctionNode。
 */
class ContainsImaginaryVisitor : public LMCAS::detail::RecursiveSymbolicVisitor {
public:
    bool found = false;

    void visit(const NumberNode&) override {}
    void visit(const VariableNode&) override {}

    void visit(const AddNode& node) override {
        for (auto& op : node.operands()) {
            if (found) {
                return;
            }
            op->accept(*this);
        }
    }

    void visit(const MultiplyNode& node) override {
        for (auto& op : node.operands()) {
            if (found) {
                return;
            }
            op->accept(*this);
        }
    }

    void visit(const PowerNode& node) override {
        if (found) {
            return;
        }

        if (power_is_imaginary(node)) {
            found = true;
            return;
        }
        node.base()->accept(*this);
        if (found) {
            return;
        }
        node.exponent()->accept(*this);
    }

    void visit(const FunctionNode& node) override {
        if (square_root_is_imaginary(node)) {
            found = true;
            return;
        }
        for (auto& arg : node.arguments()) {
            if (found) {
                return;
            }
            arg->accept(*this);
        }
    }

    void visit(const MatrixNode& node) override {
        if (std::holds_alternative<MatrixNode::DenseStorage>(node.storage())) {
            for (auto& item : std::get<MatrixNode::DenseStorage>(node.storage())) {
                if (found) {
                    return;
                }
                if (item) {
                    item->accept(*this);
                }
            }
        } else {
            for (auto& [idx, item] : std::get<MatrixNode::SparseStorage>(node.storage())) {
                (void)idx;
                if (found) {
                    return;
                }
                if (item) {
                    item->accept(*this);
                }
            }
        }
    }

    void visit(const RelationalNode& node) override {
        node.left()->accept(*this);
        if (!found) {
            node.right()->accept(*this);
        }
    }

    void visit(const LogicalNode& node) override {
        node.left()->accept(*this);
        if (!found && node.right()) {
            node.right()->accept(*this);
        }
    }

    void visit(const PiecewiseNode& node) override {
        for (const auto& branch : node.branches()) {
            if (found) {
                return;
            }
            branch.expression->accept(*this);
            if (!found) {
                branch.condition->accept(*this);
            }
        }
        if (!found && node.default_expr()) {
            node.default_expr()->accept(*this);
        }
    }

    void visit(const SummationNode& node) override {
        node.body()->accept(*this);
        if (!found) {
            node.lower_bound()->accept(*this);
        }
        if (!found) {
            node.upper_bound()->accept(*this);
        }
    }

    void visit(const ProductNode& node) override {
        node.body()->accept(*this);
        if (!found) {
            node.lower_bound()->accept(*this);
        }
        if (!found) {
            node.upper_bound()->accept(*this);
        }
    }

    void visit(const TransformNode& node) override {
        node.body()->accept(*this);
    }

    void visit(const QuantifierNode& node) override {
        node.domain()->accept(*this);
        if (!found) {
            node.predicate()->accept(*this);
        }
    }

    void visit(const SetBuilderNode& node) override {
        node.domain()->accept(*this);
        if (!found) {
            node.predicate()->accept(*this);
        }
    }

    void visit(const ComplexNode&) override {
        found = true;
    }
    void visit(const FiniteSetNode& node) override {
        for (const auto& element : node.elements()) {
            if (found) {
                return;
            }
            element->accept(*this);
        }
    }
    void visit(const IntervalNode& node) override {
        node.lower()->accept(*this);
        if (!found) {
            node.upper()->accept(*this);
        }
    }
    void visit(const MembershipNode& node) override {
        node.element()->accept(*this);
        if (!found) {
            node.set()->accept(*this);
        }
    }
    void visit(const QuantityNode& node) override { node.value()->accept(*this); }
};

}
bool contains_imaginary(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return false;
    }
    ContainsImaginaryVisitor visitor;
    LMCAS::detail::node(expr)->accept(visitor);
    return visitor.found;
}

}
