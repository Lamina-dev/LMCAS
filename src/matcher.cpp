#include "matcher.hpp"
#include "symbolic_ast.hpp"
#include "assumption_context.hpp"
#include <algorithm>
#include <iostream>

namespace lamina {

SymbolicExpr wildcard(const std::string& name) {
    return lamina::detail::expression_from_node(SymbolicFactory::create_variable(name));
}

static bool match_recursive(const std::shared_ptr<const SymbolicNode>& p_node,
                            const std::shared_ptr<const SymbolicNode>& t_node,
                            const std::unordered_set<std::string>& wildcards,
                            MatchMap& results);

static bool match_commutative_recursive(const std::vector<std::shared_ptr<const SymbolicNode>>& p_ops,
                                        const std::vector<std::shared_ptr<const SymbolicNode>>& t_ops,
                                        std::vector<bool>& used_t,
                                        size_t p_index,
                                        const std::unordered_set<std::string>& wildcards,
                                        MatchMap& results) {

    if (p_index == p_ops.size()) {
        return true;
    }

    for (size_t j = 0; j < t_ops.size(); ++j) {
        if (!used_t[j]) {

            MatchMap saved_results = results;

            if (match_recursive(p_ops[p_index], t_ops[j], wildcards, results)) {
                used_t[j] = true;
                if (match_commutative_recursive(p_ops, t_ops, used_t, p_index + 1, wildcards, results)) {
                    return true;
                }

                used_t[j] = false;
                results = saved_results;
            } else {

                results = saved_results;
            }
        }
    }

    return false;
}

static bool is_wildcard(const std::shared_ptr<const SymbolicNode>& node,
                       const std::unordered_set<std::string>& wildcards,
                       std::string& name_out) {
    if (!node) return false;

    auto var = std::dynamic_pointer_cast<const VariableNode>(node);
    if (var) {
        if (wildcards.count(var->name())) {
            name_out = var->name();
            return true;
        }
    }
    return false;
}

static bool match_recursive(const std::shared_ptr<const SymbolicNode>& p_node,
                            const std::shared_ptr<const SymbolicNode>& t_node,
                            const std::unordered_set<std::string>& wildcards,
                            MatchMap& results) {
    if (!p_node) return !t_node;
    if (!t_node) return false;

    std::string wname;
    if (is_wildcard(p_node, wildcards, wname)) {

        if (const auto existing = results.find(wname); existing != results.end()) {
            return lamina::detail::node(existing->second)->equals(*t_node);
        } else {
            results.emplace(wname, lamina::detail::expression_from_node(t_node));
            return true;
        }
    }

    if (p_node->type_priority() != t_node->type_priority()) return false;

    if (p_node->is_number() || std::dynamic_pointer_cast<const VariableNode>(p_node)) {
        return p_node->equals(*t_node);
    }

    if (auto p_add = std::dynamic_pointer_cast<const AddNode>(p_node)) {
        auto t_add = std::dynamic_pointer_cast<const AddNode>(t_node);

        if (p_add->operands().size() > t_add->operands().size()) return false;

        std::vector<bool> used(t_add->operands().size(), false);
        if (match_commutative_recursive(p_add->operands(), t_add->operands(), used, 0, wildcards, results)) {

            std::vector<std::shared_ptr<const SymbolicNode>> rest;
            for (size_t i = 0; i < used.size(); ++i) {
                if (!used[i]) {
                    rest.push_back(t_add->operands()[i]);
                }
            }

            if (!rest.empty()) {

                std::shared_ptr<const SymbolicNode> rest_node;
                if (rest.size() == 1) rest_node = rest[0];
                else rest_node = SymbolicFactory::create_add(rest);

                if (const auto found = results.find("__Add_REST__"); found != results.end()) {
                    auto existing = found->second;

                    std::vector<std::shared_ptr<const SymbolicNode>> combined_ops;
                    if (auto existing_add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(existing))) {
                        combined_ops.insert(combined_ops.end(), existing_add->operands().begin(), existing_add->operands().end());
                    } else {
                        combined_ops.push_back(lamina::detail::node(existing));
                    }
                    combined_ops.push_back(rest_node);
                    found->second = lamina::detail::expression_from_node(
                        SymbolicFactory::create_add(combined_ops));
                } else {
                    results.emplace("__Add_REST__", lamina::detail::expression_from_node(rest_node));
                }
            }
            return true;
        }
        return false;
    }

    if (auto p_mul = std::dynamic_pointer_cast<const MultiplyNode>(p_node)) {
        auto t_mul = std::dynamic_pointer_cast<const MultiplyNode>(t_node);

        if (p_mul->operands().size() > t_mul->operands().size()) return false;

        std::vector<bool> used(t_mul->operands().size(), false);
        if (match_commutative_recursive(p_mul->operands(), t_mul->operands(), used, 0, wildcards, results)) {

            std::vector<std::shared_ptr<const SymbolicNode>> rest;
            for (size_t i = 0; i < used.size(); ++i) {
                if (!used[i]) {
                    rest.push_back(t_mul->operands()[i]);
                }
            }
            if (!rest.empty()) {
                std::shared_ptr<const SymbolicNode> rest_node;
                if (rest.size() == 1) rest_node = rest[0];
                else rest_node = SymbolicFactory::create_multiply(rest);

                if (const auto found = results.find("__Mul_REST__"); found != results.end()) {
                    auto existing = found->second;
                    std::vector<std::shared_ptr<const SymbolicNode>> combined_ops;
                    if (auto existing_mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(existing))) {
                        combined_ops.insert(combined_ops.end(), existing_mul->operands().begin(), existing_mul->operands().end());
                    } else {
                        combined_ops.push_back(lamina::detail::node(existing));
                    }
                    combined_ops.push_back(rest_node);
                    found->second = lamina::detail::expression_from_node(
                        SymbolicFactory::create_multiply(combined_ops));
                } else {
                    results.emplace("__Mul_REST__", lamina::detail::expression_from_node(rest_node));
                }
            }
            return true;
        }
        return false;
    }

    if (auto p_pow = std::dynamic_pointer_cast<const PowerNode>(p_node)) {
        auto t_pow = std::dynamic_pointer_cast<const PowerNode>(t_node);
        return match_recursive(p_pow->base(), t_pow->base(), wildcards, results) &&
               match_recursive(p_pow->exponent(), t_pow->exponent(), wildcards, results);
    }

    if (auto p_func = std::dynamic_pointer_cast<const FunctionNode>(p_node)) {
        auto t_func = std::dynamic_pointer_cast<const FunctionNode>(t_node);
        if (p_func->type() != t_func->type()) return false;
        if (p_func->arguments().size() != t_func->arguments().size()) return false;
        for (size_t i = 0; i < p_func->arguments().size(); ++i) {
            if (!match_recursive(p_func->arguments()[i], t_func->arguments()[i], wildcards, results)) {
                return false;
            }
        }
        return true;
    }

    return p_node->equals(*t_node);
}

bool Matcher::match(const SymbolicExpr& pattern, const SymbolicExpr& target,
                  const std::unordered_set<std::string>& wildcards,
                  MatchMap& results) {
    if (!lamina::detail::node(pattern)) return !lamina::detail::node(target);
    return match_recursive(lamina::detail::node(pattern), lamina::detail::node(target), wildcards, results);
}

class ReplacementVisitor : public lamina::detail::SymbolicVisitor {
    const MatchMap& bindings;
    std::shared_ptr<const SymbolicNode> result;

public:
    explicit ReplacementVisitor(const MatchMap& b) : bindings(b) {}

    std::shared_ptr<const SymbolicNode> get_result() const { return result; }

    void visit(const NumberNode& node) override {
        result = node.clone();
    }

    void visit(const VariableNode& node) override {

        if (bindings.count(node.name())) {

            result = lamina::detail::node(bindings.at(node.name()))->clone();
        } else {
            result = node.clone();
        }
    }

    void visit(const AddNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        for (auto& op : node.operands()) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        result = SymbolicFactory::create_add(new_ops);
    }

    void visit(const MultiplyNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        for (auto& op : node.operands()) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        result = SymbolicFactory::create_multiply(new_ops);
    }

    void visit(const PowerNode& node) override {
        node.base()->accept(*this);
        auto new_base = result;
        node.exponent()->accept(*this);
        auto new_exp = result;

        result = lamina::detail::make_node<PowerNode>(new_base, new_exp);
    }

    void visit(const FunctionNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        for (auto& arg : node.arguments()) {
            arg->accept(*this);
            new_args.push_back(result);
        }
        result = lamina::detail::make_node<FunctionNode>(node.type(), new_args);
    }

    void visit(const MatrixNode& node) override {

        result = node.clone();
    }

    void visit(const RelationalNode& node) override {
        node.left()->accept(*this);
        auto left = result;
        node.right()->accept(*this);
        auto right = result;
        result = lamina::detail::make_node<RelationalNode>(left, right, node.op());
    }

    void visit(const LogicalNode& node) override {
        node.left()->accept(*this);
        auto left = result;
        std::shared_ptr<const SymbolicNode> right = nullptr;
        if (node.right()) {
            node.right()->accept(*this);
            right = result;
        }
        result = lamina::detail::make_node<LogicalNode>(left, right, node.op());
    }

    void visit(const PiecewiseNode& node) override {
        std::vector<PiecewiseNode::Branch> branches;
        branches.reserve(node.branches().size());
        for (const auto& branch : node.branches()) {
            branch.expression->accept(*this);
            auto expression = result;
            branch.condition->accept(*this);
            auto condition = result;
            branches.push_back({expression, condition});
        }
        std::shared_ptr<const SymbolicNode> default_expr = nullptr;
        if (node.default_expr()) {
            node.default_expr()->accept(*this);
            default_expr = result;
        }
        result = lamina::detail::make_node<PiecewiseNode>(std::move(branches), default_expr);
    }

    void visit(const SummationNode& node) override {
        node.lower_bound()->accept(*this);
        auto lower = result;
        node.upper_bound()->accept(*this);
        auto upper = result;
        std::shared_ptr<const SymbolicNode> body;
        if (bindings.count(node.index_var())) {
            body = node.body()->clone();
        } else {
            node.body()->accept(*this);
            body = result;
        }
        result = lamina::detail::make_node<SummationNode>(body, node.index_var(), lower, upper);
    }

    void visit(const ProductNode& node) override {
        node.lower_bound()->accept(*this);
        auto lower = result;
        node.upper_bound()->accept(*this);
        auto upper = result;
        std::shared_ptr<const SymbolicNode> body;
        if (bindings.count(node.index_var())) {
            body = node.body()->clone();
        } else {
            node.body()->accept(*this);
            body = result;
        }
        result = lamina::detail::make_node<ProductNode>(body, node.index_var(), lower, upper);
    }

    void visit(const TransformNode& node) override {
        std::shared_ptr<const SymbolicNode> body;
        if (bindings.count(node.source_var())) {
            body = node.body()->clone();
        } else {
            node.body()->accept(*this);
            body = result;
        }
        result = lamina::detail::make_node<TransformNode>(
            node.transform_type(), body, node.source_var(), node.target_var());
    }

    void visit(const QuantifierNode& node) override {
        node.domain()->accept(*this);
        auto domain = result;
        std::shared_ptr<const SymbolicNode> predicate;
        if (bindings.count(node.bound_var())) {
            predicate = node.predicate()->clone();
        } else {
            node.predicate()->accept(*this);
            predicate = result;
        }
        result = lamina::detail::make_node<QuantifierNode>(
            node.quantifier_type(), node.bound_var(), domain, predicate);
    }

    void visit(const SetBuilderNode& node) override {
        node.domain()->accept(*this);
        auto domain = result;
        std::shared_ptr<const SymbolicNode> predicate;
        if (bindings.count(node.element_var())) {
            predicate = node.predicate()->clone();
        } else {
            node.predicate()->accept(*this);
            predicate = result;
        }
        result = lamina::detail::make_node<SetBuilderNode>(node.element_var(), domain, predicate);
    }

    void visit(const ComplexNode& node) override {
        node.real()->accept(*this);
        auto real = result;
        node.imag()->accept(*this);
        auto imag = result;
        result = SymbolicFactory::create_complex(real, imag);
    }
    void visit(const FiniteSetNode& node) override {
        std::vector<std::shared_ptr<const SymbolicNode>> elements;
        for (const auto& element : node.elements()) {
            element->accept(*this);
            elements.push_back(result);
        }
        result = lamina::detail::make_node<FiniteSetNode>(std::move(elements));
    }
    void visit(const IntervalNode& node) override {
        node.lower()->accept(*this);
        auto lower = result;
        node.upper()->accept(*this);
        result = lamina::detail::make_node<IntervalNode>(
            lower, result, node.lower_closed(), node.upper_closed());
    }
    void visit(const MembershipNode& node) override {
        node.element()->accept(*this);
        auto element = result;
        node.set()->accept(*this);
        result = lamina::detail::make_node<MembershipNode>(element, result);
    }
    void visit(const QuantityNode& node) override {
        node.value()->accept(*this);
        result = lamina::detail::make_node<QuantityNode>(
            result, node.dimension(), node.scale_to_base(), node.display_unit());
    }
};

SymbolicExpr Matcher::replace(const SymbolicExpr& template_expr, const MatchMap& bindings, bool use_rest) {
    if (!lamina::detail::node(template_expr)) return template_expr;

    ReplacementVisitor visitor(bindings);
    lamina::detail::node(template_expr)->accept(visitor);

    auto res = lamina::detail::expression_from_node(visitor.get_result());
    if (use_rest) {

        if (bindings.find("__Add_REST__") != bindings.end()) {
            auto rest = bindings.at("__Add_REST__");

            std::vector<std::shared_ptr<const SymbolicNode>> ops;
            if (auto add_node = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(res))) {
                ops = add_node->operands();
            } else {
                ops.push_back(lamina::detail::node(res));
            }

            if (auto rest_add = std::dynamic_pointer_cast<const AddNode>(lamina::detail::node(rest))) {
                ops.insert(ops.end(), rest_add->operands().begin(), rest_add->operands().end());
            } else {
                ops.push_back(lamina::detail::node(rest));
            }
            res = lamina::detail::expression_from_node(SymbolicFactory::create_add(ops));
        }

        if (bindings.find("__Mul_REST__") != bindings.end()) {
            auto rest = bindings.at("__Mul_REST__");

            std::vector<std::shared_ptr<const SymbolicNode>> ops;
            if (auto mul_node = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(res))) {
                ops = mul_node->operands();
            } else {
                ops.push_back(lamina::detail::node(res));
            }

            if (auto rest_mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(rest))) {
                ops.insert(ops.end(), rest_mul->operands().begin(), rest_mul->operands().end());
            } else {
                ops.push_back(lamina::detail::node(rest));
            }
            res = lamina::detail::expression_from_node(SymbolicFactory::create_multiply(ops));
        }
    }

    return res;
}

void RewriteEngine::add_rule(const Rule& rule) {
    rules.push_back(rule);
}

class RewriteVisitor : public lamina::detail::SymbolicVisitor {
public:
    const RewriteEngine& engine;
    ComputationContext& context;
    std::shared_ptr<const SymbolicNode> result;
    bool changed = false;

    RewriteVisitor(const RewriteEngine& e, ComputationContext& ctx)
        : engine(e), context(ctx) {}

    std::shared_ptr<const SymbolicNode> get_result() const { return result; }

    std::shared_ptr<const SymbolicNode> try_match(std::shared_ptr<const SymbolicNode> node) {
        auto current_expr = lamina::detail::expression_from_node(node);
        const auto& rules = engine.get_rules();
        const AssumptionContext* ctx = context.assumptions().get();
        for (const auto& rule : rules) {
            MatchMap bindings;

            if (Matcher::match(rule.pattern, current_expr, rule.wildcards, bindings)) {

                // Evaluate condition: prefer assumption_condition when context is available
                if (rule.assumption_condition && ctx) {
                    if (!rule.assumption_condition(bindings, ctx)) {
                        continue;
                    }
                } else if (rule.assumption_condition && !ctx) {
                    // No context available; fall back to plain condition if present
                    if (rule.condition && !rule.condition(bindings)) {
                        continue;
                    }
                } else if (rule.condition) {
                    if (!rule.condition(bindings)) {
                        continue;
                    }
                }

                SymbolicExpr new_expr = Matcher::replace(rule.replacement, bindings, true);

                changed = true;
                return lamina::detail::node(new_expr);
            }
        }
        return node;
    }

    std::shared_ptr<const SymbolicNode> visit_child(const std::shared_ptr<const SymbolicNode>& child,
                                             bool& child_changed) {
        bool saved_changed = changed;
        changed = false;
        child->accept(*this);
        auto rewritten = result;
        if (changed || !rewritten->equals(*child)) {
            child_changed = true;
        }
        changed = saved_changed;
        return rewritten;
    }

    void finish_rewrite(std::shared_ptr<const SymbolicNode> node, bool child_changed,
                        bool original_changed) {
        changed = child_changed;
        auto matched = try_match(node);
        if (matched != node && !matched->equals(*node)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(const NumberNode& node) override {

        result = try_match(node.clone());
    }

    void visit(const VariableNode& node) override {
        result = try_match(node.clone());
    }

    void visit(const AddNode& node) override {

        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        bool child_changed = false;

        bool original_changed = changed;
        changed = false;

        for (auto& op : node.operands()) {
            bool current_changed = changed;
            changed = false;
            op->accept(*this);
            new_ops.push_back(result);
            if (changed) child_changed = true;
            else changed = current_changed;
        }

        std::shared_ptr<const SymbolicNode> node_to_match;
        if (child_changed) {
            node_to_match = SymbolicFactory::create_add(new_ops);
        } else {
            node_to_match = lamina::detail::make_node<AddNode>(node.operands());
        }

        changed = child_changed;
        auto matched = try_match(node_to_match);
        if (matched != node_to_match && !matched->equals(*node_to_match)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(const MultiplyNode& node) override {
         std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        bool child_changed = false;
        bool original_changed = changed;
        changed = false;

        for (auto& op : node.operands()) {
            bool current_changed = changed;
            changed = false;
            op->accept(*this);
            new_ops.push_back(result);
            if (changed) child_changed = true;
            else changed = current_changed;
        }

        std::shared_ptr<const SymbolicNode> node_to_match;
        if (child_changed) {
            node_to_match = SymbolicFactory::create_multiply(new_ops);
        } else {
             node_to_match = lamina::detail::make_node<MultiplyNode>(node.operands());
        }

        changed = child_changed;
        auto matched = try_match(node_to_match);
        if (matched != node_to_match && !matched->equals(*node_to_match)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(const PowerNode& node) override {
        bool original_changed = changed;
        changed = false;

        node.base()->accept(*this);
        auto new_base = result;
        bool base_changed = changed;

        changed = false;
        node.exponent()->accept(*this);
        auto new_exp = result;
        bool exp_changed = changed;

        bool child_changed = base_changed || exp_changed;
        auto node_to_match = lamina::detail::make_node<PowerNode>(new_base, new_exp);

        changed = child_changed;
        auto matched = try_match(node_to_match);
        if (matched != node_to_match && !matched->equals(*node_to_match)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(const FunctionNode& node) override {
        bool original_changed = changed;
        changed = false;
        bool child_changed = false;

        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        for (auto& arg : node.arguments()) {
            bool current_changed = changed;
            changed = false;
            arg->accept(*this);
            new_args.push_back(result);
            if (changed) child_changed = true;
            else changed = current_changed;
        }

        auto node_to_match = lamina::detail::make_node<FunctionNode>(node.type(), new_args);

        changed = child_changed;
        auto matched = try_match(node_to_match);
        if (matched != node_to_match && !matched->equals(*node_to_match)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(const MatrixNode& node) override {
        result = try_match(node.clone());
    }

    void visit(const RelationalNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto left = visit_child(node.left(), child_changed);
        auto right = visit_child(node.right(), child_changed);
        finish_rewrite(lamina::detail::make_node<RelationalNode>(left, right, node.op()),
                       child_changed, original_changed);
    }

    void visit(const LogicalNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto left = visit_child(node.left(), child_changed);
        std::shared_ptr<const SymbolicNode> right = nullptr;
        if (node.right()) right = visit_child(node.right(), child_changed);
        finish_rewrite(lamina::detail::make_node<LogicalNode>(left, right, node.op()),
                       child_changed, original_changed);
    }

    void visit(const PiecewiseNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        std::vector<PiecewiseNode::Branch> branches;
        branches.reserve(node.branches().size());
        for (const auto& branch : node.branches()) {
            auto expression = visit_child(branch.expression, child_changed);
            auto condition = visit_child(branch.condition, child_changed);
            branches.push_back({expression, condition});
        }
        std::shared_ptr<const SymbolicNode> default_expr = nullptr;
        if (node.default_expr()) default_expr = visit_child(node.default_expr(), child_changed);
        finish_rewrite(lamina::detail::make_node<PiecewiseNode>(std::move(branches), default_expr),
                       child_changed, original_changed);
    }

    void visit(const SummationNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto body = visit_child(node.body(), child_changed);
        auto lower = visit_child(node.lower_bound(), child_changed);
        auto upper = visit_child(node.upper_bound(), child_changed);
        finish_rewrite(lamina::detail::make_node<SummationNode>(body, node.index_var(), lower, upper),
                       child_changed, original_changed);
    }

    void visit(const ProductNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto body = visit_child(node.body(), child_changed);
        auto lower = visit_child(node.lower_bound(), child_changed);
        auto upper = visit_child(node.upper_bound(), child_changed);
        finish_rewrite(lamina::detail::make_node<ProductNode>(body, node.index_var(), lower, upper),
                       child_changed, original_changed);
    }

    void visit(const TransformNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto body = visit_child(node.body(), child_changed);
        finish_rewrite(lamina::detail::make_node<TransformNode>(
                           node.transform_type(), body, node.source_var(), node.target_var()),
                       child_changed, original_changed);
    }

    void visit(const QuantifierNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto domain = visit_child(node.domain(), child_changed);
        auto predicate = visit_child(node.predicate(), child_changed);
        finish_rewrite(lamina::detail::make_node<QuantifierNode>(
                           node.quantifier_type(), node.bound_var(), domain, predicate),
                       child_changed, original_changed);
    }

    void visit(const SetBuilderNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto domain = visit_child(node.domain(), child_changed);
        auto predicate = visit_child(node.predicate(), child_changed);
        finish_rewrite(lamina::detail::make_node<SetBuilderNode>(
                           node.element_var(), domain, predicate),
                       child_changed, original_changed);
    }

    void visit(const ComplexNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto real = visit_child(node.real(), child_changed);
        auto imag = visit_child(node.imag(), child_changed);
        finish_rewrite(SymbolicFactory::create_complex(real, imag),
                       child_changed, original_changed);
    }
    void visit(const FiniteSetNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        std::vector<std::shared_ptr<const SymbolicNode>> elements;
        for (const auto& element : node.elements()) {
            elements.push_back(visit_child(element, child_changed));
        }
        finish_rewrite(lamina::detail::make_node<FiniteSetNode>(std::move(elements)),
                       child_changed, original_changed);
    }
    void visit(const IntervalNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto lower = visit_child(node.lower(), child_changed);
        auto upper = visit_child(node.upper(), child_changed);
        finish_rewrite(lamina::detail::make_node<IntervalNode>(
                           lower, upper, node.lower_closed(), node.upper_closed()),
                       child_changed, original_changed);
    }
    void visit(const MembershipNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto element = visit_child(node.element(), child_changed);
        auto set = visit_child(node.set(), child_changed);
        finish_rewrite(lamina::detail::make_node<MembershipNode>(element, set),
                       child_changed, original_changed);
    }
    void visit(const QuantityNode& node) override {
        bool original_changed = changed;
        bool child_changed = false;
        auto value = visit_child(node.value(), child_changed);
        finish_rewrite(lamina::detail::make_node<QuantityNode>(
                           value, node.dimension(), node.scale_to_base(), node.display_unit()),
                       child_changed, original_changed);
    }
};

Result<SymbolicExpr> RewriteEngine::apply_step_checked(
    const SymbolicExpr& expr,
    ComputationContext& context) const {
    auto budget = context.consume_steps(1, "rewrite.apply_step");
    if (!budget) return Result<SymbolicExpr>::failure(budget.error());
    if (!lamina::detail::node(expr)) {
        return Result<SymbolicExpr>::failure(
            CasErrc::InvalidArgument, "rewrite expression must not be null",
            "rewrite.apply_step");
    }
    RewriteVisitor v(*this, context);
    lamina::detail::node(expr)->accept(v);
    return Result<SymbolicExpr>::success(
        lamina::detail::expression_from_node(v.get_result()));
}

Result<SymbolicExpr> RewriteEngine::apply_checked(
    const SymbolicExpr& expr,
    ComputationContext& context,
    int max_iterations) const {
    if (max_iterations < 0) {
        return Result<SymbolicExpr>::failure(
            CasErrc::InvalidArgument, "rewrite iteration count must be non-negative",
            "rewrite.apply");
    }
    SymbolicExpr current = expr;
    for (int i = 0; i < max_iterations; ++i) {
        auto next_result = apply_step_checked(current, context);
        if (!next_result) return next_result;
        SymbolicExpr next = std::move(next_result.value());
        if (lamina::detail::node(next)->equals(*lamina::detail::node(current))) {
            return Result<SymbolicExpr>::success(std::move(current));
        }
        current = next;
    }
    return Result<SymbolicExpr>::success(std::move(current));
}

SymbolicExpr RewriteEngine::apply_step(const SymbolicExpr& expr) const {
    ComputationContext context;
    auto result = apply_step_checked(expr, context);
    if (!result) throw std::runtime_error(result.error().message);
    return std::move(result.value());
}

SymbolicExpr RewriteEngine::apply(const SymbolicExpr& expr, int max_iterations) const {
    ComputationContext context;
    auto result = apply_checked(expr, context, max_iterations);
    if (!result) throw std::runtime_error(result.error().message);
    return std::move(result.value());
}

}
