#include "solve_mixed_transcendental.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/normalization_utils.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS {
namespace {
/**
 * @internal
 * @brief 判断函数类型是否为 sin、cos、tan、exp 或 ln。
 */
static bool is_transcendental_func(FunctionNode::FuncType t) {
    switch (t) {
    case FunctionNode::FuncType::Sin:
    case FunctionNode::FuncType::Cos:
    case FunctionNode::FuncType::Tan:
    case FunctionNode::FuncType::Exp:
    case FunctionNode::FuncType::Ln:
        return true;
    default:
        return false;
    }
}

    struct TranscendentalDetector : public LMCAS::detail::RecursiveSymbolicVisitor {
        const std::string& target_var;
        bool found = false;

        explicit TranscendentalDetector(const std::string& v) : target_var(v) {}

        void visit(const NumberNode&) override {}
        void visit(const VariableNode&) override {}
        void visit(const MatrixNode&) override {}
        void visit(const RelationalNode& n) override {
            if (found) {
                return;
            }
            if (n.left()) {
                n.left()->accept(*this);
            }
            if (!found && n.right()) {
                n.right()->accept(*this);
            }
        }
        void visit(const LogicalNode& n) override {
            if (found) {
                return;
            }
            if (n.left()) {
                n.left()->accept(*this);
            }
            if (!found && n.right()) {
                n.right()->accept(*this);
            }
        }
        void visit(const PiecewiseNode& n) override {
            for (const auto& branch : n.branches()) {
                if (found) {
                    return;
                }
                branch.expression->accept(*this);
                if (!found) {
                    branch.condition->accept(*this);
                }
            }
            if (!found && n.default_expr()) {
                n.default_expr()->accept(*this);
            }
        }
        void visit(const SummationNode& n) override {
            n.body()->accept(*this);
            if (!found) {
                n.lower_bound()->accept(*this);
            }
            if (!found) {
                n.upper_bound()->accept(*this);
            }
        }
        void visit(const ProductNode& n) override {
            n.body()->accept(*this);
            if (!found) {
                n.lower_bound()->accept(*this);
            }
            if (!found) {
                n.upper_bound()->accept(*this);
            }
        }
        void visit(const TransformNode& n) override {
            n.body()->accept(*this);
        }
        void visit(const QuantifierNode& n) override {
            n.domain()->accept(*this);
            if (!found) {
                n.predicate()->accept(*this);
            }
        }
        void visit(const SetBuilderNode& n) override {
            n.domain()->accept(*this);
            if (!found) {
                n.predicate()->accept(*this);
            }
        }
        void visit(const ComplexNode& n) override {
            n.real()->accept(*this);
            if (!found) {
                n.imag()->accept(*this);
            }
        }
        void visit(const FiniteSetNode& n) override { for (const auto& e : n.elements()) { if (found) {
            return;
        } e->accept(*this); } }
        void visit(const IntervalNode& n) override { n.lower()->accept(*this); if (!found) {
            n.upper()->accept(*this);
        } }
        void visit(const MembershipNode& n) override { n.element()->accept(*this); if (!found) {
            n.set()->accept(*this);
        } }
        void visit(const QuantityNode& n) override { n.value()->accept(*this); }

        void visit(const AddNode& n) override {
            for (auto& op : n.operands()) {
                if (found) {
                    return;
                }
                op->accept(*this);
            }
        }

        void visit(const MultiplyNode& n) override {
            for (auto& op : n.operands()) {
                if (found) {
                    return;
                }
                op->accept(*this);
            }
        }

        void visit(const PowerNode& n) override {
            if (found) {
                return;
            }
            n.base()->accept(*this);
            if (!found) {
                n.exponent()->accept(*this);
            }
        }

        void visit(const FunctionNode& n) override {
            if (found) {
                return;
            }

            if (is_transcendental_func(n.type())) {
                for (auto& arg : n.arguments()) {
                    if (expression_depends_on_variable(arg, target_var)) {
                        found = true;
                        return;
                    }
                }
            }
            for (auto& arg : n.arguments()) {
                if (found) {
                    return;
                }
                arg->accept(*this);
            }
        }
};

static BigInt degree_in_var(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& var);

static BigInt sum_degree(const AddNode& sum, const std::string& var) {
    BigInt maximum;
    for (const auto& operand : sum.operands()) {
        auto degree = degree_in_var(operand, var);
        if (degree < 0) {
            return -1;
        }
        if (degree > maximum) {
            maximum = std::move(degree);
        }
    }
    return maximum;
}

static BigInt product_degree(const MultiplyNode& product, const std::string& var) {
    BigInt total;
    for (const auto& operand : product.operands()) {
        auto degree = degree_in_var(operand, var);
        if (degree < 0) {
            return -1;
        }
        total += degree;
    }
    return total;
}

/**
 * @internal
 * @brief 计算节点关于指定变量的多项式次数。
 * 与 var 无关的子表达式次数为 0。
 * @param[in] node 符号节点。
 * @param[in] var 变量名。
 * @return 关于 var 的次数；超出多项式次数支持域时返回 -1。
 */
static BigInt degree_in_var(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    if (!node) {
        return 0;
    }

    if (!expression_depends_on_variable(node, var)) {
        return 0;
    }

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(node)) {
        return (v->name() == var) ? 1 : 0;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return sum_degree(*add, var);
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return product_degree(*mul, var);
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent());
        BigInt e_val;
        if (!try_get_integer_value(exp_num, e_val) || e_val.is_negative()) {
            return -1;
        }

        BigInt base_deg = degree_in_var(pow->base(), var);
        if (base_deg < 0) {
            return -1;
        }
        return base_deg * e_val;
    }
    if (std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return -1;
    }

    return 0;
}

}

bool contains_transcendental_of_var(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var)
{
    if (!expr || !LMCAS::detail::node(expr)) {
        return false;
    }
    TranscendentalDetector detector(var);
    LMCAS::detail::node(expr)->accept(detector);
    return detector.found;
}

bool is_polynomial_after_substitution(
    const TransSubstitutionResult& sub_result)
{
    if (sub_result.mappings.empty()) {
        return false;
    }
    if (!sub_result.poly_expr || !LMCAS::detail::node(sub_result.poly_expr)) {
        return false;
    }

    const auto& root = LMCAS::detail::node(sub_result.poly_expr);

    for (const auto& m : sub_result.mappings) {
        if (expression_depends_on_variable(root, m.indeterminate) &&
            degree_in_var(root, m.indeterminate).is_negative()) {
            return false;
        }
    }

    return true;
}

}
