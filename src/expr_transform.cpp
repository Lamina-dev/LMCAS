#include "expr.hpp"
#include "matcher.hpp"
#include "symbolic_ast.hpp"
#include "internal/expr_common.hpp"
#include <exception>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace LMCAS {

using namespace expr_detail::expr_common;

ExprResult simplify(const ExprPtr& expression, ComputationContext& context) {
    return checked_transform_expr(
        expression, context, kSimplifyOperation, "simplify",
        [](const SymbolicExpr& value) { return value.simplify(); });
}

ExprResult simplify(const ExprPtr& expression) {
    ComputationContext context;
    return simplify(expression, context);
}

ExprResult expand(const ExprPtr& expression, ComputationContext& context) {
    return checked_transform_expr(
        expression, context, kExpandOperation, "expand",
        [](const SymbolicExpr& value) { return value.expand(); });
}

ExprResult expand(const ExprPtr& expression) {
    ComputationContext context;
    return expand(expression, context);
}

ExprResult differentiate(const ExprPtr& expression,
                         const std::string& variable,
                         ComputationContext& context) {
    if (variable.empty()) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "differentiate variable cannot be empty",
                                  kDifferentiateOperation);
    }
    return checked_transform_expr(
        expression, context, kDifferentiateOperation, "differentiate",
        [&variable](const SymbolicExpr& value) {
            return value.differentiate(variable);
        });
}

ExprResult differentiate(const ExprPtr& expression,
                         const std::string& variable) {
    ComputationContext context;
    return differentiate(expression, variable, context);
}

ExprResult substitute(const ExprPtr& expression,
                      const std::string& variable,
                      const ExprPtr& value,
                      ComputationContext& context) {
    auto step = context.consume_steps(1, kSubstituteOperation);
    if (!step) return ExprResult::failure(step.error());
    if (!expression) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "expression cannot be null",
                                  kSubstituteOperation);
    }
    if (variable.empty()) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "substitution variable cannot be empty",
                                  kSubstituteOperation);
    }
    if (!value) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "substitution value cannot be null",
                                  kSubstituteOperation);
    }
    try {
        auto result = expression->substitute(variable, value);
        if (!result || !LMCAS::detail::node(result)) {
            return expression_failure(CasErrc::InternalInvariant,
                                      "substitution returned an empty expression",
                                      kSubstituteOperation);
        }
        return ExprResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return expression_failure(CasErrc::ResourceLimit,
                                  "substitution allocation failed",
                                  kSubstituteOperation);
    } catch (const std::exception& error) {
        return expression_failure(CasErrc::InvalidArgument, error.what(),
                                  kSubstituteOperation);
    }
}

ExprResult substitute(const ExprPtr& expression,
                      const std::string& variable,
                      const ExprPtr& value) {
    ComputationContext context;
    return substitute(expression, variable, value, context);
}

BindingResult binding(ExprPtr symbol, ExprPtr value) {
    if (!symbol || !LMCAS::detail::node(symbol)) {
        return BindingResult::failure(
            CasErrc::InvalidArgument, "binding symbol cannot be null",
            kSubstituteOperation);
    }
    if (!value || !LMCAS::detail::node(value)) {
        return BindingResult::failure(
            CasErrc::InvalidArgument, "binding value cannot be null",
            kSubstituteOperation);
    }
    if (!std::dynamic_pointer_cast<const VariableNode>(
            LMCAS::detail::node(symbol))) {
        return BindingResult::failure(
            CasErrc::InvalidArgument,
            "binding left-hand side must be a symbol",
            kSubstituteOperation);
    }
    return BindingResult::success(Binding{std::move(symbol), std::move(value)});
}

ExprResult substitute(const ExprPtr& expression,
                      const Binding& replacement,
                      ComputationContext& context) {
    auto checked = binding(replacement.symbol, replacement.value);
    if (!checked) return ExprResult::failure(checked.error());
    const auto variable = std::dynamic_pointer_cast<const VariableNode>(
        LMCAS::detail::node(checked.value().symbol));
    return substitute(expression, variable->name(), checked.value().value, context);
}

ExprResult substitute(const ExprPtr& expression,
                      const Binding& replacement) {
    ComputationContext context;
    return substitute(expression, replacement, context);
}

ExprResult substitute(const ExprPtr& expression,
                      const std::vector<Binding>& replacements,
                      ComputationContext& context) {
    if (!expression || !LMCAS::detail::node(expression)) {
        return expression_failure(CasErrc::InvalidArgument,
                                  "expression cannot be null",
                                  kSubstituteOperation);
    }
    auto result = expression;
    for (const auto& replacement : replacements) {
        auto next = substitute(result, replacement, context);
        if (!next) return next;
        result = std::move(next.value());
    }
    return ExprResult::success(std::move(result));
}

ExprResult substitute(const ExprPtr& expression,
                      const std::vector<Binding>& replacements) {
    ComputationContext context;
    return substitute(expression, replacements, context);
}

ExprMatchResult expr_match(const ExprPtr& pattern,
                           const ExprPtr& target,
                           const std::vector<std::string>& wildcards,
                           ComputationContext& context) {
    auto step = context.consume_steps(1, kExprMatchOperation);
    if (!step) return ExprMatchResult::failure(step.error());
    if (!pattern || !LMCAS::detail::node(pattern)) {
        return expr_match_failure(CasErrc::InvalidArgument,
                                  "match pattern cannot be null",
                                  kExprMatchOperation);
    }
    if (!target || !LMCAS::detail::node(target)) {
        return expr_match_failure(CasErrc::InvalidArgument,
                                  "match target cannot be null",
                                  kExprMatchOperation);
    }

    std::unordered_set<std::string> wildcard_set;
    wildcard_set.reserve(wildcards.size());
    for (const auto& wildcard : wildcards) {
        if (wildcard.empty()) {
            return expr_match_failure(CasErrc::InvalidArgument,
                                      "wildcard names cannot be empty",
                                      kExprMatchOperation);
        }
        if (!wildcard_set.insert(wildcard).second) {
            return expr_match_failure(CasErrc::InvalidArgument,
                                      "wildcard names must be unique",
                                      kExprMatchOperation);
        }
    }

    try {
        MatchMap raw_bindings;
        const bool matched = Matcher::match(
            *pattern, *target, wildcard_set, raw_bindings);
        if (!matched) {
            return ExprMatchResult::success(ExprMatch{false, {}});
        }

        std::vector<std::string> names;
        names.reserve(raw_bindings.size());
        for (const auto& binding : raw_bindings) {
            names.push_back(binding.first);
        }
        std::sort(names.begin(), names.end());

        ExprMatch result;
        result.matched = true;
        result.bindings.reserve(names.size());
        for (const auto& name : names) {
            const auto& value = raw_bindings.at(name);
            auto node = LMCAS::detail::node(value);
            if (!node) {
                return expr_match_failure(CasErrc::InternalInvariant,
                                          "matcher produced a null binding",
                                          kExprMatchOperation);
            }
            result.bindings.push_back(
                ExprMatchBinding{name,
                                 LMCAS::detail::make_expression_ptr(node)});
        }
        return ExprMatchResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return expr_match_failure(CasErrc::ResourceLimit,
                                  "expression matching allocation failed",
                                  kExprMatchOperation);
    } catch (const std::exception& error) {
        return expr_match_failure(CasErrc::InternalInvariant, error.what(),
                                  kExprMatchOperation);
    }
}

ExprMatchResult expr_match(const ExprPtr& pattern,
                           const ExprPtr& target,
                           const std::vector<std::string>& wildcards) {
    ComputationContext context;
    return expr_match(pattern, target, wildcards, context);
}
} // namespace LMCAS
