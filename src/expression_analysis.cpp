#include "internal/expression_analysis.hpp"

#include <unordered_map>
#include <utility>

namespace lamina {
namespace {

class FreeVariableCollector final : public detail::RecursiveSymbolicVisitor {
public:
    std::set<std::string> variables;

    void visit(const VariableNode& node) override {
        if (!is_bound(node.name())) variables.insert(node.name());
    }

    void visit(const SummationNode& node) override {
        visit_child(node.lower_bound());
        visit_child(node.upper_bound());
        visit_scoped(node.index_var(), node.body());
    }

    void visit(const ProductNode& node) override {
        visit_child(node.lower_bound());
        visit_child(node.upper_bound());
        visit_scoped(node.index_var(), node.body());
    }

    void visit(const IntegralNode& node) override {
        visit_child(node.lower());
        visit_child(node.upper());
        visit_scoped(node.variable(), node.body());
    }

    void visit(const TransformNode& node) override {
        visit_child(node.target());
        visit_scoped(node.source_var(), node.body());
    }

    void visit(const QuantifierNode& node) override {
        visit_child(node.domain());
        visit_scoped(node.bound_var(), node.predicate());
    }

    void visit(const SetBuilderNode& node) override {
        visit_child(node.domain());
        visit_scoped(node.element_var(), node.predicate());
    }

    void visit(const LimitNode& node) override {
        visit_child(node.point());
        visit_scoped(node.variable(), node.body());
    }

    void visit(const RootOfNode&) override {}

private:
    bool is_bound(const std::string& name) const {
        const auto found = bound_.find(name);
        return found != bound_.end() && found->second != 0;
    }

    void visit_scoped(const std::string& name, const detail::SymbolicNodePtr& body) {
        ++bound_[name];
        visit_child(body);
        auto found = bound_.find(name);
        if (--found->second == 0) bound_.erase(found);
    }

    std::unordered_map<std::string, std::size_t> bound_;
};

class AllNameCollector final : public detail::RecursiveSymbolicVisitor {
public:
    std::set<std::string> names;

    void visit(const VariableNode& node) override { names.insert(node.name()); }
    void visit(const SummationNode& node) override {
        names.insert(node.index_var());
        detail::RecursiveSymbolicVisitor::visit(node);
    }
    void visit(const ProductNode& node) override {
        names.insert(node.index_var());
        detail::RecursiveSymbolicVisitor::visit(node);
    }
    void visit(const IntegralNode& node) override {
        names.insert(node.variable());
        detail::RecursiveSymbolicVisitor::visit(node);
    }
    void visit(const TransformNode& node) override {
        names.insert(node.source_var());
        detail::RecursiveSymbolicVisitor::visit(node);
    }
    void visit(const QuantifierNode& node) override {
        names.insert(node.bound_var());
        detail::RecursiveSymbolicVisitor::visit(node);
    }
    void visit(const SetBuilderNode& node) override {
        names.insert(node.element_var());
        detail::RecursiveSymbolicVisitor::visit(node);
    }
    void visit(const LimitNode& node) override {
        names.insert(node.variable());
        detail::RecursiveSymbolicVisitor::visit(node);
    }
    void visit(const RootOfNode&) override {}
};

std::set<std::string> all_names(const detail::SymbolicNodePtr& expression) {
    AllNameCollector collector;
    if (expression) expression->accept(collector);
    return collector.names;
}

std::string fresh_name(const std::string& base, std::set<std::string>& occupied) {
    for (std::size_t suffix = 1;; ++suffix) {
        auto candidate = base + "_" + std::to_string(suffix);
        if (occupied.insert(candidate).second) return candidate;
    }
}

class FreeSubstitution final : public detail::SymbolicRewriter {
public:
    FreeSubstitution(std::string target, detail::SymbolicNodePtr replacement,
                     std::set<std::string> replacement_free,
                     std::set<std::string> occupied)
        : target_(std::move(target)), replacement_(std::move(replacement)),
          replacement_free_(std::move(replacement_free)),
          occupied_(std::move(occupied)) {}

    void visit(const VariableNode& node) override {
        set_result(node.name() == target_ ? replacement_ : current());
    }

    void visit(const SummationNode& node) override {
        bool changed = false;
        auto lower = rewrite_child(node.lower_bound(), changed);
        auto upper = rewrite_child(node.upper_bound(), changed);
        auto binder = node.index_var();
        auto body = rewrite_scoped(node.body(), binder, changed);
        set_result(changed ? detail::make_node<SummationNode>(body, binder, lower, upper)
                           : current());
    }

    void visit(const ProductNode& node) override {
        bool changed = false;
        auto lower = rewrite_child(node.lower_bound(), changed);
        auto upper = rewrite_child(node.upper_bound(), changed);
        auto binder = node.index_var();
        auto body = rewrite_scoped(node.body(), binder, changed);
        set_result(changed ? detail::make_node<ProductNode>(body, binder, lower, upper)
                           : current());
    }

    void visit(const IntegralNode& node) override {
        bool changed = false;
        auto lower = node.lower() ? rewrite_child(node.lower(), changed) : nullptr;
        auto upper = node.upper() ? rewrite_child(node.upper(), changed) : nullptr;
        auto binder = node.variable();
        auto body = rewrite_scoped(node.body(), binder, changed);
        set_result(changed ? detail::make_node<IntegralNode>(body, binder, lower, upper)
                           : current());
    }

    void visit(const TransformNode& node) override {
        bool changed = false;
        auto target = rewrite_child(node.target(), changed);
        auto binder = node.source_var();
        auto body = rewrite_scoped(node.body(), binder, changed);
        set_result(changed ? detail::make_node<TransformNode>(node.transform_type(), body,
                                                              binder, target)
                           : current());
    }

    void visit(const QuantifierNode& node) override {
        bool changed = false;
        auto domain = rewrite_child(node.domain(), changed);
        auto binder = node.bound_var();
        auto predicate = rewrite_scoped(node.predicate(), binder, changed);
        set_result(changed ? detail::make_node<QuantifierNode>(node.quantifier_type(), binder,
                                                               domain, predicate)
                           : current());
    }

    void visit(const SetBuilderNode& node) override {
        bool changed = false;
        auto domain = rewrite_child(node.domain(), changed);
        auto binder = node.element_var();
        auto predicate = rewrite_scoped(node.predicate(), binder, changed);
        set_result(changed ? detail::make_node<SetBuilderNode>(binder, domain, predicate)
                           : current());
    }

    void visit(const LimitNode& node) override {
        bool changed = false;
        auto point = rewrite_child(node.point(), changed);
        auto binder = node.variable();
        auto body = rewrite_scoped(node.body(), binder, changed);
        set_result(changed ? detail::make_node<LimitNode>(
                                 body, binder, point, node.direction())
                           : current());
    }

    void visit(const RootOfNode&) override { set_result(current()); }

private:
    detail::SymbolicNodePtr rewrite_scoped(
        const detail::SymbolicNodePtr& original_body,
        std::string& binder,
        bool& changed) {
        if (binder == target_) return original_body;

        auto body = original_body;
        if (expression_depends_on_variable(body, target_) &&
            replacement_free_.find(binder) != replacement_free_.end()) {
            const auto renamed = fresh_name(binder, occupied_);
            body = substitute_free(body, binder,
                                   SymbolicFactory::create_variable(renamed));
            binder = renamed;
            changed = true;
        }
        return rewrite_child(body, changed);
    }

    std::string target_;
    detail::SymbolicNodePtr replacement_;
    const std::set<std::string> replacement_free_;
    std::set<std::string> occupied_;
};

} // namespace

std::set<std::string> free_variables(const detail::SymbolicNodePtr& expression) {
    FreeVariableCollector collector;
    if (expression) expression->accept(collector);
    return collector.variables;
}

bool expression_depends_on_variable(
    const detail::SymbolicNodePtr& expression,
    const std::string& variable) {
    if (!expression || variable.empty()) return false;
    const auto variables = free_variables(expression);
    return variables.find(variable) != variables.end();
}

detail::SymbolicNodePtr substitute_free(
    const detail::SymbolicNodePtr& expression,
    const std::string& variable,
    const detail::SymbolicNodePtr& replacement) {
    if (!expression) return nullptr;
    if (variable.empty()) throw std::invalid_argument("substitution variable cannot be empty");
    if (!replacement) throw std::invalid_argument("substitution replacement cannot be null");

    auto replacement_free = free_variables(replacement);
    auto occupied = all_names(expression);
    auto replacement_names = all_names(replacement);
    occupied.insert(replacement_names.begin(), replacement_names.end());
    occupied.insert(variable);

    FreeSubstitution substitution(variable, replacement,
                                  std::move(replacement_free),
                                  std::move(occupied));
    return substitution.rewrite(expression);
}

} // namespace lamina
