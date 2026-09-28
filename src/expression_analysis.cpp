#include "internal/expression_analysis.hpp"

#include <unordered_map>
#include <utility>

namespace LMCAS {
namespace {

class FreeVariableCollector final : public detail::RecursiveSymbolicVisitor {
public:
    using RecursiveSymbolicVisitor::RecursiveSymbolicVisitor;
    std::set<std::string> variables;

    void visit(const VariableNode& node) override {
        if (!node.is_constant() && !is_bound(node.name())) variables.insert(node.name());
    }

    void visit(const SummationNode& node) override { visit_binder(node); }
    void visit(const ProductNode& node) override { visit_binder(node); }
    void visit(const IntegralNode& node) override { visit_binder(node); }
    void visit(const TransformNode& node) override { visit_binder(node); }
    void visit(const QuantifierNode& node) override { visit_binder(node); }
    void visit(const SetBuilderNode& node) override { visit_binder(node); }
    void visit(const LimitNode& node) override { visit_binder(node); }

    void visit(const RootOfNode&) override {}

private:
    bool is_bound(const std::string& name) const {
        const auto found = bound_.find(name);
        return found != bound_.end() && found->second != 0;
    }

    void visit_binder(const SymbolicNode& node) {
        const auto binder = detail::binder_view(node);
        visit_children(binder->outside_scope);
        ++bound_[binder->bound_name];
        visit_child(binder->scoped_body);
        auto found = bound_.find(binder->bound_name);
        if (--found->second == 0) bound_.erase(found);
    }

    std::unordered_map<std::string_view, std::size_t> bound_;
};
class FreeVariableFinder final : public detail::RecursiveSymbolicVisitor {
public:
    explicit FreeVariableFinder(std::string_view name, detail::RewriteBudget* budget = nullptr)
        : RecursiveSymbolicVisitor(budget), name_(name) {}
    bool found = false;

    void visit(const VariableNode& node) override {
        found = found || (!node.is_constant() && node.name() == name_);
    }
    void visit(const SummationNode& node) override { visit_binder(node); }
    void visit(const ProductNode& node) override { visit_binder(node); }
    void visit(const IntegralNode& node) override { visit_binder(node); }
    void visit(const TransformNode& node) override { visit_binder(node); }
    void visit(const QuantifierNode& node) override { visit_binder(node); }
    void visit(const SetBuilderNode& node) override { visit_binder(node); }
    void visit(const LimitNode& node) override { visit_binder(node); }

private:
    void visit_binder(const SymbolicNode& node) {
        if (found) return;
        const auto binder = detail::binder_view(node);
        visit_children(binder->outside_scope);
        if (!found && binder->bound_name != name_) visit_child(binder->scoped_body);
    }

    const std::string_view name_;
};

class AllNameCollector final : public detail::RecursiveSymbolicVisitor {
public:
    using RecursiveSymbolicVisitor::RecursiveSymbolicVisitor;
    std::set<std::string> names;

    void visit(const VariableNode& node) override {
        if (!node.is_constant()) names.insert(node.name());
    }
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

std::set<std::string> all_names(const detail::SymbolicNodePtr& expression,
                                detail::RewriteBudget* budget = nullptr) {
    AllNameCollector collector(budget);
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
                     std::set<std::string> occupied, detail::RewriteBudget* budget)
        : SymbolicRewriter(budget), target_(std::move(target)),
          replacement_(std::move(replacement)),
          replacement_free_(std::move(replacement_free)),
          occupied_(std::move(occupied)) {}

    void visit(const VariableNode& node) override {
        set_result(!node.is_constant() && node.name() == target_ ? replacement_ : current());
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
        if (expression_depends_on_variable(body, target_, rewrite_budget()) &&
            replacement_free_.find(binder) != replacement_free_.end()) {
            const auto renamed = fresh_name(binder, occupied_);
            body = substitute_free(body, binder,
                                   SymbolicFactory::create_variable(renamed), rewrite_budget());
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
std::optional<Rational> detail::exact_rational_value(
    const NumberNode& number) {
    if (const auto* integer = std::get_if<BigInt>(&number.value())) {
        return Rational(*integer);
    }
    if (const auto* rational = std::get_if<Rational>(&number.value())) {
        return *rational;
    }
    return std::nullopt;
}

std::optional<Rational> detail::exact_rational_value(
    const detail::SymbolicNodePtr& node) {
    const auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    return number ? exact_rational_value(*number) : std::nullopt;
}


Result<std::optional<detail::AffineForm>> detail::recognize_affine(
    const SymbolicExpr& expression, const std::string& variable,
    ComputationContext& context) {
    using AffineResult = Result<std::optional<AffineForm>>;
    auto converted = symbolic_to_poly_checked(expression, variable, context);
    if (!converted) {
        if (converted.error().code == CasErrc::UnsupportedExpression) {
            return std::optional<AffineForm>{};
        }
        return AffineResult::failure(converted.error());
    }
    const auto& polynomial = converted.value();
    if (polynomial.degree() > 1) return std::optional<AffineForm>{};
    return std::optional<AffineForm>{AffineForm{
        polynomial.coeffs.size() > 1 ? polynomial.coeffs[1].val : SymbolicExpr::number(0),
        polynomial.coeffs.empty() ? SymbolicExpr::number(0) : polynomial.coeffs[0].val}};
}

std::optional<int> exact_small_integer_node(
    const detail::SymbolicNodePtr& node, int min_value, int max_value) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return std::nullopt;

    BigInt value;
    if (std::holds_alternative<BigInt>(number->value())) {
        value = std::get<BigInt>(number->value());
    } else if (std::holds_alternative<Rational>(number->value())) {
        const Rational& rational = std::get<Rational>(number->value());
        if (!rational.is_integer()) return std::nullopt;
        value = rational.to_bigint();
    } else {
        return std::nullopt;
    }

    const auto exact = value.try_to_int64();
    if (!exact || *exact < min_value || *exact > max_value) {
        return std::nullopt;
    }
    return static_cast<int>(*exact);
}

bool detail::is_imaginary_unit_name(const std::string& name) {
    return name == "I";
}

std::set<std::string> detail::all_variable_names(const detail::SymbolicNodePtr& expression) {
    return all_names(expression);
}

std::set<std::string> free_variables(const detail::SymbolicNodePtr& expression,
                                     detail::RewriteBudget* budget) {
    FreeVariableCollector collector(budget);
    if (expression) expression->accept(collector);
    return collector.variables;
}

bool expression_depends_on_variable(
    const detail::SymbolicNodePtr& expression,
    const std::string& variable, detail::RewriteBudget* budget) {
    if (!expression || variable.empty()) return false;
    FreeVariableFinder finder(variable, budget);
    expression->accept(finder);
    return finder.found;
}

detail::SymbolicNodePtr substitute_free(
    const detail::SymbolicNodePtr& expression,
    const std::string& variable,
    const detail::SymbolicNodePtr& replacement, detail::RewriteBudget* budget) {
    if (!expression) return nullptr;
    if (variable.empty()) throw std::invalid_argument("substitution variable cannot be empty");
    if (!replacement) throw std::invalid_argument("substitution replacement cannot be null");

    auto replacement_free = free_variables(replacement, budget);
    auto occupied = all_names(expression, budget);
    auto replacement_names = all_names(replacement, budget);
    occupied.insert(replacement_names.begin(), replacement_names.end());
    occupied.insert(variable);

    FreeSubstitution substitution(variable, replacement,
                                  std::move(replacement_free),
                                  std::move(occupied), budget);
    return substitution.rewrite(expression);
}

} // namespace LMCAS
