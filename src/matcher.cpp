#include "matcher.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/expression_analysis.hpp"
#include <algorithm>
#include <type_traits>

namespace LMCAS {

SymbolicExpr wildcard(const std::string& name) {
    return LMCAS::detail::expression_from_node(SymbolicFactory::create_variable(name));
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
    if (var && !var->is_constant()) {
        if (wildcards.count(var->name())) {
            name_out = var->name();
            return true;
        }
    }
    return false;
}

template <typename Operation>
static std::shared_ptr<const SymbolicNode> combine_operands(
    std::vector<std::shared_ptr<const SymbolicNode>> operands) {
    if constexpr (std::is_same_v<Operation, AddNode>) {
        return SymbolicFactory::create_add(std::move(operands));
    } else {
        return SymbolicFactory::create_multiply(std::move(operands));
    }
}

template <typename Operation>
static void record_remainder(const Operation& target, const std::vector<bool>& used,
                             const char* key, MatchMap& results) {
    std::vector<std::shared_ptr<const SymbolicNode>> rest;
    for (size_t index = 0; index < used.size(); ++index) {
        if (!used[index]) {
            rest.push_back(target.operands()[index]);
        }
    }
    if (rest.empty()) {
        return;
    }
    auto remainder = rest.size() == 1 ? rest.front() : combine_operands<Operation>(std::move(rest));
    const auto found = results.find(key);
    if (found == results.end()) {
        results.emplace(key, detail::expression_from_node(remainder));
        return;
    }
    std::vector<std::shared_ptr<const SymbolicNode>> operands;
    if (const auto* existing = dynamic_cast<const Operation*>(detail::node(found->second).get())) {
        operands.insert(operands.end(), existing->operands().begin(), existing->operands().end());
    } else {
        operands.push_back(detail::node(found->second));
    }
    operands.push_back(remainder);
    found->second = detail::expression_from_node(combine_operands<Operation>(std::move(operands)));
}

template <typename Operation>
static bool match_operation(const Operation& pattern, const Operation& target,
                            const std::unordered_set<std::string>& wildcards,
                            const char* rest_key, MatchMap& results) {
    if (pattern.operands().size() > target.operands().size()) {
        return false;
    }
    std::vector<bool> used(target.operands().size(), false);
    if (!match_commutative_recursive(pattern.operands(), target.operands(), used, 0, wildcards, results)) {
        return false;
    }
    record_remainder(target, used, rest_key, results);
    return true;
}

static bool match_function(const FunctionNode& pattern, const FunctionNode& target,
                           const std::unordered_set<std::string>& wildcards, MatchMap& results) {
    if (pattern.type() != target.type() || pattern.arguments().size() != target.arguments().size()) {
        return false;
    }
    for (size_t index = 0; index < pattern.arguments().size(); ++index) {
        if (!match_recursive(pattern.arguments()[index], target.arguments()[index], wildcards, results)) {
            return false;
        }
    }
    return true;
}

static bool match_recursive(const std::shared_ptr<const SymbolicNode>& pattern,
                            const std::shared_ptr<const SymbolicNode>& target,
                            const std::unordered_set<std::string>& wildcards, MatchMap& results) {
    if (!pattern) {
        return !target;
    }
    if (!target) {
        return false;
    }
    std::string name;
    if (is_wildcard(pattern, wildcards, name)) {
        const auto existing = results.find(name);
        if (existing != results.end()) {
            return detail::node(existing->second)->equals(*target);
        }
        results.emplace(name, detail::expression_from_node(target));
        return true;
    }
    if (pattern->type_priority() != target->type_priority()) {
        return false;
    }
    if (pattern->is_number() || dynamic_cast<const VariableNode*>(pattern.get())) {
        return pattern->equals(*target);
    }
    if (const auto* add = dynamic_cast<const AddNode*>(pattern.get())) {
        return match_operation(*add, static_cast<const AddNode&>(*target), wildcards, "__Add_REST__", results);
    }
    if (const auto* multiply = dynamic_cast<const MultiplyNode*>(pattern.get())) {
        return match_operation(*multiply, static_cast<const MultiplyNode&>(*target), wildcards, "__Mul_REST__", results);
    }
    if (const auto* power = dynamic_cast<const PowerNode*>(pattern.get())) {
        const auto& other = static_cast<const PowerNode&>(*target);
        return match_recursive(power->base(), other.base(), wildcards, results) &&
               match_recursive(power->exponent(), other.exponent(), wildcards, results);
    }
    if (const auto* function = dynamic_cast<const FunctionNode*>(pattern.get())) {
        return match_function(*function, static_cast<const FunctionNode&>(*target), wildcards, results);
    }
    return pattern->equals(*target);
}

bool Matcher::match(const SymbolicExpr& pattern, const SymbolicExpr& target,
                  const std::unordered_set<std::string>& wildcards,
                  MatchMap& results) {
    if (!LMCAS::detail::node(pattern)) return !LMCAS::detail::node(target);
    return match_recursive(LMCAS::detail::node(pattern), LMCAS::detail::node(target), wildcards, results);
}


SymbolicExpr Matcher::replace(const SymbolicExpr& template_expr, const MatchMap& bindings, bool use_rest) {
    if (!LMCAS::detail::node(template_expr)) return template_expr;

    auto replaced = LMCAS::detail::node(template_expr);
    std::set<std::string> occupied = LMCAS::free_variables(replaced);
    for (const auto& [name, expression] : bindings) {
        const auto names = LMCAS::free_variables(LMCAS::detail::node(expression));
        occupied.insert(names.begin(), names.end());
    }

    std::vector<std::pair<std::string, std::string>> placeholders;
    placeholders.reserve(bindings.size());
    std::size_t suffix = 0;
    for (const auto& [name, expression] : bindings) {
        (void)expression;
        std::string placeholder;
        do {
            placeholder = "__lmcas_match_binding_" + std::to_string(suffix++);
        } while (!occupied.insert(placeholder).second);
        replaced = LMCAS::substitute_free(
            replaced, name, SymbolicFactory::create_variable(placeholder));
        placeholders.emplace_back(name, std::move(placeholder));
    }
    for (const auto& [name, placeholder] : placeholders) {
        replaced = LMCAS::substitute_free(
            replaced, placeholder, LMCAS::detail::node(bindings.at(name)));
    }

    auto res = LMCAS::detail::expression_from_node(std::move(replaced));
    if (use_rest) {

        if (bindings.find("__Add_REST__") != bindings.end()) {
            auto rest = bindings.at("__Add_REST__");

            std::vector<std::shared_ptr<const SymbolicNode>> ops;
            if (auto add_node = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(res))) {
                ops = add_node->operands();
            } else {
                ops.push_back(LMCAS::detail::node(res));
            }

            if (auto rest_add = std::dynamic_pointer_cast<const AddNode>(LMCAS::detail::node(rest))) {
                ops.insert(ops.end(), rest_add->operands().begin(), rest_add->operands().end());
            } else {
                ops.push_back(LMCAS::detail::node(rest));
            }
            res = LMCAS::detail::expression_from_node(SymbolicFactory::create_add(ops));
        }

        if (bindings.find("__Mul_REST__") != bindings.end()) {
            auto rest = bindings.at("__Mul_REST__");

            std::vector<std::shared_ptr<const SymbolicNode>> ops;
            if (auto mul_node = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(res))) {
                ops = mul_node->operands();
            } else {
                ops.push_back(LMCAS::detail::node(res));
            }

            if (auto rest_mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(rest))) {
                ops.insert(ops.end(), rest_mul->operands().begin(), rest_mul->operands().end());
            } else {
                ops.push_back(LMCAS::detail::node(rest));
            }
            res = LMCAS::detail::expression_from_node(SymbolicFactory::create_multiply(ops));
        }
    }

    return res;
}

}
