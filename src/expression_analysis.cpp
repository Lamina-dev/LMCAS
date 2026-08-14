#include "internal/expression_analysis.hpp"

namespace lamina {
namespace {

class FreeVariableQuery final : public detail::SymbolicVisitor {
public:
    explicit FreeVariableQuery(std::string variable)
        : variable_(std::move(variable)) {}

    bool found() const noexcept { return found_; }

    void visit(const NumberNode&) override {}
    void visit(const VariableNode& node) override {
        found_ = node.name() == variable_;
    }
    void visit(const AddNode& node) override { visit_all(node.operands()); }
    void visit(const MultiplyNode& node) override { visit_all(node.operands()); }
    void visit(const PowerNode& node) override {
        visit_one(node.base());
        visit_one(node.exponent());
    }
    void visit(const FunctionNode& node) override { visit_all(node.arguments()); }
    void visit(const MatrixNode& node) override {
        if (std::holds_alternative<MatrixNode::DenseStorage>(node.storage())) {
            visit_all(std::get<MatrixNode::DenseStorage>(node.storage()));
            return;
        }
        for (const auto& [index, value] :
             std::get<MatrixNode::SparseStorage>(node.storage())) {
            (void)index;
            visit_one(value);
        }
    }
    void visit(const RelationalNode& node) override {
        visit_one(node.left());
        visit_one(node.right());
    }
    void visit(const LogicalNode& node) override {
        visit_one(node.left());
        visit_one(node.right());
    }
    void visit(const PiecewiseNode& node) override {
        for (const auto& branch : node.branches()) {
            visit_one(branch.expression);
            visit_one(branch.condition);
        }
        visit_one(node.default_expr());
    }
    void visit(const SummationNode& node) override {
        visit_one(node.lower_bound());
        visit_one(node.upper_bound());
        if (node.index_var() != variable_) visit_one(node.body());
    }
    void visit(const ProductNode& node) override {
        visit_one(node.lower_bound());
        visit_one(node.upper_bound());
        if (node.index_var() != variable_) visit_one(node.body());
    }
    void visit(const TransformNode& node) override {
        if (node.source_var() != variable_) visit_one(node.body());
        found_ = found_ || node.target_var() == variable_;
    }
    void visit(const QuantifierNode& node) override {
        visit_one(node.domain());
        if (node.bound_var() != variable_) visit_one(node.predicate());
    }
    void visit(const SetBuilderNode& node) override {
        visit_one(node.domain());
        if (node.element_var() != variable_) visit_one(node.predicate());
    }
    void visit(const FiniteSetNode& node) override { visit_all(node.elements()); }
    void visit(const IntervalNode& node) override {
        visit_one(node.lower());
        visit_one(node.upper());
    }
    void visit(const MembershipNode& node) override {
        visit_one(node.element());
        visit_one(node.set());
    }
    void visit(const QuantityNode& node) override { visit_one(node.value()); }
    void visit(const ComplexNode& node) override {
        visit_one(node.real());
        visit_one(node.imag());
    }

private:
    void visit_one(const std::shared_ptr<const SymbolicNode>& node) {
        if (!found_ && node) node->accept(*this);
    }

    void visit_all(
        const std::vector<std::shared_ptr<const SymbolicNode>>& nodes) {
        for (const auto& node : nodes) visit_one(node);
    }

    std::string variable_;
    bool found_ = false;
};

} // namespace

bool expression_depends_on_variable(
    const std::shared_ptr<const SymbolicNode>& expression,
    const std::string& variable) {
    if (!expression) return false;
    FreeVariableQuery query(variable);
    expression->accept(query);
    return query.found();
}

} // namespace lamina
