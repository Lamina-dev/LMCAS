#include "solve_transcendental.hpp"
#include "internal/symbolic_ast.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/normalization_utils.hpp"
#include <vector>

namespace LMCAS {
namespace {

static bool function_candidate(const FunctionNode& function, const std::string& var) {
    if (function.arguments().size() != 1 ||
        !expression_depends_on_variable(function.arguments()[0], var)) {
        return false;
    }
    using Type = FunctionNode::FuncType;
    return function.type() == Type::Exp || function.type() == Type::Sin ||
           function.type() == Type::Cos || function.type() == Type::Tan;
}

static bool power_candidate(const PowerNode& power, const std::string& var) {
    if (!expression_depends_on_variable(power.base(), var) ||
        expression_depends_on_variable(power.exponent(), var)) {
        return false;
    }
    auto exponent = std::dynamic_pointer_cast<const NumberNode>(power.exponent());
    BigInt value;
    return try_get_integer_value(exponent, value) && value >= 2;
}


static void collect_transcendental_subexprs(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var,
    std::vector<std::shared_ptr<SymbolicExpr>>& candidates) {

    if (!node) {
        return;
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (function_candidate(*func, var)) {
            candidates.push_back(LMCAS::detail::make_expression_ptr(node));
        }

        for (auto& arg : func->arguments()) {
            collect_transcendental_subexprs(arg, var, candidates);
        }
        return;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        if (power_candidate(*pow, var)) {
            candidates.push_back(LMCAS::detail::make_expression_ptr(pow->base()));
        }

        collect_transcendental_subexprs(pow->base(), var, candidates);
        collect_transcendental_subexprs(pow->exponent(), var, candidates);
        return;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (auto& op : add->operands()) {
            collect_transcendental_subexprs(op, var, candidates);
        }
        return;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (auto& op : mul->operands()) {
            collect_transcendental_subexprs(op, var, candidates);
        }
        return;
    }
}

static void deduplicate_candidates(std::vector<std::shared_ptr<SymbolicExpr>>& candidates) {
    std::vector<std::shared_ptr<SymbolicExpr>> unique;
    std::vector<std::string> seen;
    for (auto& c : candidates) {
        std::string s = c->to_string();
        bool found = false;
        for (auto& existing : seen) {
            if (existing == s) { found = true; break; }
        }
        if (!found) {
            seen.push_back(s);
            unique.push_back(c);
        }
    }
    candidates = std::move(unique);
}

static BigInt product_exp_multiplier(const MultiplyNode& product, const std::string& var) {
    std::shared_ptr<const NumberNode> number;
    bool has_variable = false;
    bool has_other = false;
    for (const auto& operand : product.operands()) {
        if (auto numeric = std::dynamic_pointer_cast<const NumberNode>(operand)) {
            if (!number) {
                number = numeric;
            } else {
                has_other = true;
            }
        } else if (auto variable = std::dynamic_pointer_cast<const VariableNode>(operand)) {
            if (!variable->is_constant() && variable->name() == var) {
                has_variable = true;
            } else {
                has_other = true;
            }
        } else {
            has_other = true;
        }
    }
    BigInt value;
    if (has_variable && number && !has_other &&
        try_get_integer_value(number, value) && value > 0) {
        return value;
    }
    return 0;
}

static BigInt extract_exp_multiplier(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    auto func = std::dynamic_pointer_cast<const FunctionNode>(node);
    if (!func || func->type() != FunctionNode::FuncType::Exp || func->arguments().size() != 1) {
        return 0;
    }

    auto arg = func->arguments()[0];

    if (auto v = std::dynamic_pointer_cast<const VariableNode>(arg)) {
        if (!v->is_constant() && v->name() == var) {
            return 1;
        }
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(arg)) {
        return product_exp_multiplier(*mul, var);
    }

    return 0;
}

static std::shared_ptr<const SymbolicNode> rewrite_exp_as_u_power(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var,
    const std::string& u_var) {

    if (!node) {
        return node;
    }

    const BigInt k = extract_exp_multiplier(node, var);
    if (k > 0) {
        auto u_node = LMCAS::detail::make_node<VariableNode>(u_var);
        if (k == 1) {
            return u_node;
        }
        return SymbolicFactory::create_power(u_node, LMCAS::detail::make_node<NumberNode>(k));
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        new_ops.reserve(add->operands().size());
        for (auto& op : add->operands()) {
            new_ops.push_back(rewrite_exp_as_u_power(op, var, u_var));
        }
        return SymbolicFactory::create_add(std::move(new_ops));
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        new_ops.reserve(mul->operands().size());
        for (auto& op : mul->operands()) {
            new_ops.push_back(rewrite_exp_as_u_power(op, var, u_var));
        }
        return SymbolicFactory::create_multiply(std::move(new_ops));
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto new_base = rewrite_exp_as_u_power(pow->base(), var, u_var);
        auto new_exp = rewrite_exp_as_u_power(pow->exponent(), var, u_var);
        return SymbolicFactory::create_power(new_base, new_exp);
    }

    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        new_args.reserve(func->arguments().size());
        for (auto& arg : func->arguments()) {
            new_args.push_back(rewrite_exp_as_u_power(arg, var, u_var));
        }
        return LMCAS::detail::make_node<FunctionNode>(func->type(), std::move(new_args));
    }

    return node;
}

static bool has_exp_k_var_terms(const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    if (!node) {
        return false;
    }

    if (extract_exp_multiplier(node, var) > 0) {
        return true;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (auto& op : add->operands()) {
            if (has_exp_k_var_terms(op, var)) {
                return true;
            }
        }
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (auto& op : mul->operands()) {
            if (has_exp_k_var_terms(op, var)) {
                return true;
            }
        }
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        if (has_exp_k_var_terms(pow->base(), var)) {
            return true;
        }
        if (has_exp_k_var_terms(pow->exponent(), var)) {
            return true;
        }
    }
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        for (auto& arg : func->arguments()) {
            if (has_exp_k_var_terms(arg, var)) {
                return true;
            }
        }
    }

    return false;
}

        struct SubstNodeVisitor {
            std::shared_ptr<const SymbolicNode> target;
            std::string u_name;

            std::shared_ptr<const SymbolicNode> replace(const std::shared_ptr<const SymbolicNode>& node) {
                if (!node) {
                    return node;
                }

                if (node->equals(*target)) {
                    return LMCAS::detail::make_node<VariableNode>(u_name);
                }

                if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
                    std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
                    new_ops.reserve(add->operands().size());
                    for (auto& op : add->operands()) {
                        new_ops.push_back(replace(op));
                    }
                    return SymbolicFactory::create_add(std::move(new_ops));
                }
                if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
                    std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
                    new_ops.reserve(mul->operands().size());
                    for (auto& op : mul->operands()) {
                        new_ops.push_back(replace(op));
                    }
                    return SymbolicFactory::create_multiply(std::move(new_ops));
                }
                if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
                    auto new_base = replace(pow->base());
                    auto new_exp = replace(pow->exponent());
                    return SymbolicFactory::create_power(new_base, new_exp);
                }
                if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
                    std::vector<std::shared_ptr<const SymbolicNode>> new_args;
                    new_args.reserve(func->arguments().size());
                    for (auto& arg : func->arguments()) {
                        new_args.push_back(replace(arg));
                    }
                    return LMCAS::detail::make_node<FunctionNode>(func->type(), std::move(new_args));
                }

                return node;
            }
        };

static Result<bool> redundant_variable_candidate(
    const std::shared_ptr<SymbolicExpr>& candidate,
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    auto variable = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(candidate));
    if (!variable || variable->is_constant() || variable->name() != var) {
        return false;
    }
    auto polynomial = symbolic_to_poly<SymbolicPolyCoeff>(expr, var);
    if (!polynomial) {
        if (polynomial.error().code == CasErrc::UnsupportedExpression) {
            return false;
        }
        return Result<bool>::failure(polynomial.error());
    }
    return polynomial.value().degree() >= 2;
}

static Result<bool> is_substitution_polynomial(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var, const std::string& u_var) {
    if (expression_depends_on_variable(LMCAS::detail::node(expr), var)) {
        return false;
    }
    auto polynomial = symbolic_to_poly<SymbolicPolyCoeff>(expr, u_var);
    if (!polynomial) {
        if (polynomial.error().code == CasErrc::UnsupportedExpression) {
            return false;
        }
        return Result<bool>::failure(polynomial.error());
    }
    return polynomial.value().degree() >= 2;
}

}

Result<std::optional<SubstitutionResult>> detect_substitution(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var) {

    if (!expr || !LMCAS::detail::node(expr)) {
        return std::optional<SubstitutionResult>{};
    }

    if (!expression_depends_on_variable(LMCAS::detail::node(expr), var)) {
        return std::optional<SubstitutionResult>{};
    }

    const std::string u_var = "_u_subst";

    std::vector<std::shared_ptr<SymbolicExpr>> candidates;
    collect_transcendental_subexprs(LMCAS::detail::node(expr), var, candidates);
    deduplicate_candidates(candidates);

    for (auto& h : candidates) {

        auto redundant = redundant_variable_candidate(h, expr, var);
        if (!redundant) {
            return Result<std::optional<SubstitutionResult>>::failure(redundant.error());
        }
        if (redundant.value()) {
            continue;
        }



        SubstNodeVisitor visitor;
        visitor.target = LMCAS::detail::node(h);
        visitor.u_name = u_var;

        auto substituted_node = visitor.replace(LMCAS::detail::node(expr));
        auto substituted_expr = LMCAS::detail::make_expression_ptr(substituted_node);
        substituted_expr = substituted_expr->simplify();

        auto polynomial = is_substitution_polynomial(substituted_expr, var, u_var);
        if (!polynomial) {
            return Result<std::optional<SubstitutionResult>>::failure(polynomial.error());
        }
        if (polynomial.value()) {
            return std::optional<SubstitutionResult>{{h, substituted_expr, u_var}};
        }
    }

    if (has_exp_k_var_terms(LMCAS::detail::node(expr), var)) {
        auto rewritten_node = rewrite_exp_as_u_power(LMCAS::detail::node(expr), var, u_var);
        auto rewritten_expr = LMCAS::detail::make_expression_ptr(rewritten_node);
        rewritten_expr = rewritten_expr->simplify();

        auto polynomial = is_substitution_polynomial(rewritten_expr, var, u_var);
        if (!polynomial) {
            return Result<std::optional<SubstitutionResult>>::failure(polynomial.error());
        }
        if (polynomial.value()) {
            auto exp_x = SymbolicExpr::exp(SymbolicExpr::variable(var));
            return std::optional<SubstitutionResult>{{exp_x, rewritten_expr, u_var}};
        }
    }

    return std::optional<SubstitutionResult>{};
}

}
