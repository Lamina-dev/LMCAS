#include "matcher.hpp"
#include "assumption_context.hpp"
#include <algorithm>
#include <iostream>

namespace lamina {

SymbolicExpr wildcard(const std::string& name) {
    return SymbolicExpr(SymbolicFactory::create_variable(name));
}

static bool match_recursive(const std::shared_ptr<SymbolicNode>& p_node,
                            const std::shared_ptr<SymbolicNode>& t_node,
                            const std::unordered_set<std::string>& wildcards,
                            MatchMap& results);

static bool match_commutative_recursive(const std::vector<std::shared_ptr<SymbolicNode>>& p_ops,
                                        const std::vector<std::shared_ptr<SymbolicNode>>& t_ops,
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

static bool is_wildcard(const std::shared_ptr<SymbolicNode>& node,
                       const std::unordered_set<std::string>& wildcards,
                       std::string& name_out) {
    if (!node) return false;

    auto var = std::dynamic_pointer_cast<VariableNode>(node);
    if (var) {
        if (wildcards.count(var->name)) {
            name_out = var->name;
            return true;
        }
    }
    return false;
}

static bool match_recursive(const std::shared_ptr<SymbolicNode>& p_node,
                            const std::shared_ptr<SymbolicNode>& t_node,
                            const std::unordered_set<std::string>& wildcards,
                            MatchMap& results) {
    if (!p_node) return !t_node;
    if (!t_node) return false;

    std::string wname;
    if (is_wildcard(p_node, wildcards, wname)) {

        if (results.count(wname)) {

            return results[wname].root->equals(*t_node);
        } else {

            results[wname] = SymbolicExpr(t_node);
            return true;
        }
    }

    if (p_node->type_priority() != t_node->type_priority()) return false;

    if (p_node->is_number() || std::dynamic_pointer_cast<VariableNode>(p_node)) {
        return p_node->equals(*t_node);
    }

    if (auto p_add = std::dynamic_pointer_cast<AddNode>(p_node)) {
        auto t_add = std::dynamic_pointer_cast<AddNode>(t_node);

        if (p_add->operands.size() > t_add->operands.size()) return false;

        std::vector<bool> used(t_add->operands.size(), false);
        if (match_commutative_recursive(p_add->operands, t_add->operands, used, 0, wildcards, results)) {

            std::vector<std::shared_ptr<SymbolicNode>> rest;
            for (size_t i = 0; i < used.size(); ++i) {
                if (!used[i]) {
                    rest.push_back(t_add->operands[i]);
                }
            }

            if (!rest.empty()) {

                std::shared_ptr<SymbolicNode> rest_node;
                if (rest.size() == 1) rest_node = rest[0];
                else rest_node = SymbolicFactory::create_add(rest);

                if (results.find("__Add_REST__") != results.end()) {
                    auto existing = results["__Add_REST__"];

                    std::vector<std::shared_ptr<SymbolicNode>> combined_ops;
                    if (auto existing_add = std::dynamic_pointer_cast<AddNode>(existing.root)) {
                        combined_ops.insert(combined_ops.end(), existing_add->operands.begin(), existing_add->operands.end());
                    } else {
                        combined_ops.push_back(existing.root);
                    }
                    combined_ops.push_back(rest_node);
                    results["__Add_REST__"] = SymbolicExpr(SymbolicFactory::create_add(combined_ops));
                } else {
                    results["__Add_REST__"] = SymbolicExpr(rest_node);
                }
            }
            return true;
        }
        return false;
    }

    if (auto p_mul = std::dynamic_pointer_cast<MultiplyNode>(p_node)) {
        auto t_mul = std::dynamic_pointer_cast<MultiplyNode>(t_node);

        if (p_mul->operands.size() > t_mul->operands.size()) return false;

        std::vector<bool> used(t_mul->operands.size(), false);
        if (match_commutative_recursive(p_mul->operands, t_mul->operands, used, 0, wildcards, results)) {

            std::vector<std::shared_ptr<SymbolicNode>> rest;
            for (size_t i = 0; i < used.size(); ++i) {
                if (!used[i]) {
                    rest.push_back(t_mul->operands[i]);
                }
            }
            if (!rest.empty()) {
                std::shared_ptr<SymbolicNode> rest_node;
                if (rest.size() == 1) rest_node = rest[0];
                else rest_node = SymbolicFactory::create_multiply(rest);

                if (results.find("__Mul_REST__") != results.end()) {
                    auto existing = results["__Mul_REST__"];
                    std::vector<std::shared_ptr<SymbolicNode>> combined_ops;
                    if (auto existing_mul = std::dynamic_pointer_cast<MultiplyNode>(existing.root)) {
                        combined_ops.insert(combined_ops.end(), existing_mul->operands.begin(), existing_mul->operands.end());
                    } else {
                        combined_ops.push_back(existing.root);
                    }
                    combined_ops.push_back(rest_node);
                    results["__Mul_REST__"] = SymbolicExpr(SymbolicFactory::create_multiply(combined_ops));
                } else {
                    results["__Mul_REST__"] = SymbolicExpr(rest_node);
                }
            }
            return true;
        }
        return false;
    }

    if (auto p_pow = std::dynamic_pointer_cast<PowerNode>(p_node)) {
        auto t_pow = std::dynamic_pointer_cast<PowerNode>(t_node);
        return match_recursive(p_pow->base, t_pow->base, wildcards, results) &&
               match_recursive(p_pow->exponent, t_pow->exponent, wildcards, results);
    }

    if (auto p_func = std::dynamic_pointer_cast<FunctionNode>(p_node)) {
        auto t_func = std::dynamic_pointer_cast<FunctionNode>(t_node);
        if (p_func->type != t_func->type) return false;
        if (p_func->arguments.size() != t_func->arguments.size()) return false;
        for (size_t i = 0; i < p_func->arguments.size(); ++i) {
            if (!match_recursive(p_func->arguments[i], t_func->arguments[i], wildcards, results)) {
                return false;
            }
        }
        return true;
    }

    return p_node->equals(*t_node);
}

bool Matcher::match(const SymbolicExpr& pattern, const SymbolicExpr& target,
                  const std::unordered_set<std::string>& wildcards,
                  MatchMap& results,
                  const AssumptionContext* /*ctx*/) {
    if (!pattern.root) return !target.root;
    return match_recursive(pattern.root, target.root, wildcards, results);
}

class ReplacementVisitor : public SymbolicVisitor {
    const MatchMap& bindings;
    std::shared_ptr<SymbolicNode> result;
    bool use_rest;

public:
    ReplacementVisitor(const MatchMap& b, bool use_rest = false) : bindings(b), use_rest(use_rest) {}

    std::shared_ptr<SymbolicNode> get_result() const { return result; }

    void visit(NumberNode& node) override {
        result = node.clone();
    }

    void visit(VariableNode& node) override {

        if (bindings.count(node.name)) {

            result = bindings.at(node.name).root->clone();
        } else {
            result = node.clone();
        }
    }

    void visit(AddNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        for (auto& op : node.operands) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        result = SymbolicFactory::create_add(new_ops);
    }

    void visit(MultiplyNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        for (auto& op : node.operands) {
            op->accept(*this);
            new_ops.push_back(result);
        }
        result = SymbolicFactory::create_multiply(new_ops);
    }

    void visit(PowerNode& node) override {
        node.base->accept(*this);
        auto new_base = result;
        node.exponent->accept(*this);
        auto new_exp = result;

        result = std::make_shared<PowerNode>(new_base, new_exp);
    }

    void visit(FunctionNode& node) override {
        std::vector<std::shared_ptr<SymbolicNode>> new_args;
        for (auto& arg : node.arguments) {
            arg->accept(*this);
            new_args.push_back(result);
        }
        result = std::make_shared<FunctionNode>(node.type, new_args);
    }

    void visit(MatrixNode& node) override {

        result = node.clone();
    }

    void visit(RelationalNode& node) override {
        result = node.clone();
    }

    void visit(LogicalNode& node) override {
        result = node.clone();
    }
};

SymbolicExpr Matcher::replace(const SymbolicExpr& template_expr, const MatchMap& bindings, bool use_rest) {
    if (!template_expr.root) return template_expr;

    ReplacementVisitor visitor(bindings, false);
    template_expr.root->accept(visitor);

    SymbolicExpr res(visitor.get_result());

    if (use_rest) {

        if (bindings.find("__Add_REST__") != bindings.end()) {
            auto rest = bindings.at("__Add_REST__");

            std::vector<std::shared_ptr<SymbolicNode>> ops;
            if (auto add_node = std::dynamic_pointer_cast<AddNode>(res.root)) {
                ops = add_node->operands;
            } else {
                ops.push_back(res.root);
            }

            if (auto rest_add = std::dynamic_pointer_cast<AddNode>(rest.root)) {
                ops.insert(ops.end(), rest_add->operands.begin(), rest_add->operands.end());
            } else {
                ops.push_back(rest.root);
            }
            res = SymbolicExpr(SymbolicFactory::create_add(ops));
        }

        if (bindings.find("__Mul_REST__") != bindings.end()) {
            auto rest = bindings.at("__Mul_REST__");

            std::vector<std::shared_ptr<SymbolicNode>> ops;
            if (auto mul_node = std::dynamic_pointer_cast<MultiplyNode>(res.root)) {
                ops = mul_node->operands;
            } else {
                ops.push_back(res.root);
            }

            if (auto rest_mul = std::dynamic_pointer_cast<MultiplyNode>(rest.root)) {
                ops.insert(ops.end(), rest_mul->operands.begin(), rest_mul->operands.end());
            } else {
                ops.push_back(rest.root);
            }
            res = SymbolicExpr(SymbolicFactory::create_multiply(ops));
        }
    }

    return res;
}

void RewriteEngine::add_rule(const Rule& rule) {
    rules.push_back(rule);
}

void RewriteEngine::set_assumption_context(const AssumptionContext* ctx) {
    assumption_ctx_ = ctx;
}

class RewriteVisitor : public SymbolicVisitor {
public:
    RewriteEngine& engine;
    std::shared_ptr<SymbolicNode> result;
    bool changed = false;

    RewriteVisitor(RewriteEngine& e) : engine(e) {}

    std::shared_ptr<SymbolicNode> get_result() const { return result; }

    std::shared_ptr<SymbolicNode> try_match(std::shared_ptr<SymbolicNode> node) {
        SymbolicExpr current_expr(node);
        const auto& rules = engine.get_rules();
        const AssumptionContext* ctx = engine.get_assumption_context();
        for (const auto& rule : rules) {
            MatchMap bindings;

            if (Matcher::match(rule.pattern, current_expr, rule.wildcards, bindings, ctx)) {

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
                return new_expr.root;
            }
        }
        return node;
    }

    void visit(NumberNode& node) override {

        result = try_match(node.clone());
    }

    void visit(VariableNode& node) override {
        result = try_match(node.clone());
    }

    void visit(AddNode& node) override {

        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        bool child_changed = false;

        bool original_changed = changed;
        changed = false;

        for (auto& op : node.operands) {
            bool current_changed = changed;
            changed = false;
            op->accept(*this);
            new_ops.push_back(result);
            if (changed) child_changed = true;
            else changed = current_changed;
        }

        std::shared_ptr<SymbolicNode> node_to_match;
        if (child_changed) {
            node_to_match = SymbolicFactory::create_add(new_ops);
        } else {
            node_to_match = std::make_shared<AddNode>(node.operands);
        }

        changed = child_changed;
        auto matched = try_match(node_to_match);
        if (matched != node_to_match && !matched->equals(*node_to_match)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(MultiplyNode& node) override {
         std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        bool child_changed = false;
        bool original_changed = changed;
        changed = false;

        for (auto& op : node.operands) {
            bool current_changed = changed;
            changed = false;
            op->accept(*this);
            new_ops.push_back(result);
            if (changed) child_changed = true;
            else changed = current_changed;
        }

        std::shared_ptr<SymbolicNode> node_to_match;
        if (child_changed) {
            node_to_match = SymbolicFactory::create_multiply(new_ops);
        } else {
             node_to_match = std::make_shared<MultiplyNode>(node.operands);
        }

        changed = child_changed;
        auto matched = try_match(node_to_match);
        if (matched != node_to_match && !matched->equals(*node_to_match)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(PowerNode& node) override {
        bool original_changed = changed;
        changed = false;

        node.base->accept(*this);
        auto new_base = result;
        bool base_changed = changed;

        changed = false;
        node.exponent->accept(*this);
        auto new_exp = result;
        bool exp_changed = changed;

        bool child_changed = base_changed || exp_changed;
        auto node_to_match = std::make_shared<PowerNode>(new_base, new_exp);

        changed = child_changed;
        auto matched = try_match(node_to_match);
        if (matched != node_to_match && !matched->equals(*node_to_match)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(FunctionNode& node) override {
        bool original_changed = changed;
        changed = false;
        bool child_changed = false;

        std::vector<std::shared_ptr<SymbolicNode>> new_args;
        for (auto& arg : node.arguments) {
            bool current_changed = changed;
            changed = false;
            arg->accept(*this);
            new_args.push_back(result);
            if (changed) child_changed = true;
            else changed = current_changed;
        }

        auto node_to_match = std::make_shared<FunctionNode>(node.type, new_args);

        changed = child_changed;
        auto matched = try_match(node_to_match);
        if (matched != node_to_match && !matched->equals(*node_to_match)) changed = true;
        result = matched;

        if (changed) original_changed = true;
        changed = original_changed;
    }

    void visit(MatrixNode& node) override {
        result = try_match(node.clone());
    }

    void visit(RelationalNode& node) override {
        result = try_match(node.clone());
    }

    void visit(LogicalNode& node) override {
        result = try_match(node.clone());
    }
};

SymbolicExpr RewriteEngine::apply_step(const SymbolicExpr& expr) {
    if (!expr.root) return expr;
    RewriteVisitor v(*this);
    expr.root->accept(v);
    return SymbolicExpr(v.get_result());
}

SymbolicExpr RewriteEngine::apply(const SymbolicExpr& expr, int max_iterations) {
    SymbolicExpr current = expr;
    for (int i = 0; i < max_iterations; ++i) {
        SymbolicExpr next = apply_step(current);
        if (next.root->equals(*current.root)) {
            return current;
        }
        current = next;
    }
    return current;
}

}
