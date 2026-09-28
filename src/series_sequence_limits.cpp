#include "limit_result.hpp"
#include "series_engine.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "internal/series_support.hpp"
#include <memory>
#include <string>

namespace LMCAS {

using detail::series_support::series_is_infinity;

static bool is_alternating_factor(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& n) {
    auto power = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!power) {
        return false;
    }
    auto base = std::dynamic_pointer_cast<const NumberNode>(power->base());
    auto exponent = std::dynamic_pointer_cast<const VariableNode>(power->exponent());
    if (!base || !exponent || exponent->is_constant() || exponent->name() != n) {
        return false;
    }
    return base->is_negative_one();
}

static std::shared_ptr<SymbolicExpr> product_without_factor(
    const MultiplyNode& product, size_t omitted) {
    std::vector<std::shared_ptr<const SymbolicNode>> rest;
    for (size_t j = 0; j < product.operands().size(); ++j) {
        if (j != omitted) {
            rest.push_back(product.operands()[j]);
        }
    }
    if (rest.empty()) {
        return SymbolicExpr::number(1);
    }
    if (rest.size() == 1) {
        return detail::make_expression_ptr(rest[0]);
    }
    return detail::make_expression_ptr(detail::make_node<MultiplyNode>(rest));
}

static bool series_extract_alternating(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& n, std::shared_ptr<SymbolicExpr>& remainder) {
    if (is_alternating_factor(node, n)) {
        remainder = SymbolicExpr::number(1);
        return true;
    }
    auto product = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!product) {
        return false;
    }
    for (size_t i = 0; i < product->operands().size(); ++i) {
        if (is_alternating_factor(product->operands()[i], n)) {
            remainder = product_without_factor(*product, i);
            return true;
        }
    }
    return false;
}


static std::shared_ptr<SymbolicExpr> series_abs(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr) {
        return nullptr;
    }
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<FunctionNode>(FunctionNode::FuncType::Abs,
        std::vector<std::shared_ptr<const SymbolicNode>>{LMCAS::detail::node(expr)}))->simplify();
}

static std::shared_ptr<SymbolicExpr> series_negate(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr) {
        return nullptr;
    }
    return SymbolicExpr::multiply(SymbolicExpr::number(-1), expr)->simplify();
}

static ExpressionResult alternating_sequence_bound(
    const std::shared_ptr<SymbolicExpr>& remainder, const std::string& n,
    bool upper, ComputationContext& context) {
    auto limited = limit_expression_checked(
        remainder, n, SymbolicExpr::infinity(), LimitDirection::Both, context);
    if (!limited) {
        return limited;
    }
    auto limit = std::move(limited.value());
    auto simplified = limit ? limit->simplify() : nullptr;
    if (simplified && !series_is_infinity(simplified)) {
        auto magnitude = series_abs(simplified);
        return upper ? magnitude : series_negate(magnitude);
    }
    if (series_is_infinity(simplified)) {
        return upper ? SymbolicExpr::infinity() : SymbolicExpr::infinity(-1);
    }
    auto magnitude = series_abs(limit);
    return upper ? magnitude : series_negate(magnitude);
}


static ExpressionResult sequence_extreme_limit(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n,
    bool upper, ComputationContext& context) {
    std::shared_ptr<SymbolicExpr> remainder;
    if (series_extract_alternating(detail::node(a_n), n, remainder)) {
        return alternating_sequence_bound(remainder, n, upper, context);
    }
    auto limited = limit_expression_checked(
        a_n, n, SymbolicExpr::infinity(),
        LimitDirection::Both, context);
    if (!limited) {
        return limited;
    }
    auto limit = std::move(limited.value());
    auto simplified = limit ? limit->simplify() : nullptr;
    return ExpressionResult::success(
        simplified ? std::move(simplified) : std::move(limit));
}

ExpressionResult lim_sup_checked(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n,
    ComputationContext& context) {
    constexpr const char* operation = "lim_sup";
    if (!a_n || !detail::node(a_n) || n.empty()) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument,
            "upper limit requires a sequence term and index variable",
            operation);
    }
    auto budget = context.consume_steps(1, operation);
    if (!budget) {
        return ExpressionResult::failure(budget.error());
    }
    try {
        return sequence_extreme_limit(a_n, n, true, context);
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while calculating upper limit", operation);
    } catch (const std::exception& ex) {
        return ExpressionResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

ExpressionResult lim_sup_checked(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n) {
    ComputationContext context;
    return lim_sup_checked(a_n, n, context);
}

ExpressionResult lim_inf_checked(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n,
    ComputationContext& context) {
    constexpr const char* operation = "lim_inf";
    if (!a_n || !detail::node(a_n) || n.empty()) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument,
            "lower limit requires a sequence term and index variable",
            operation);
    }
    auto budget = context.consume_steps(1, operation);
    if (!budget) {
        return ExpressionResult::failure(budget.error());
    }
    try {
        return sequence_extreme_limit(a_n, n, false, context);
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while calculating lower limit", operation);
    } catch (const std::exception& ex) {
        return ExpressionResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

ExpressionResult lim_inf_checked(
    const std::shared_ptr<SymbolicExpr>& a_n, const std::string& n) {
    ComputationContext context;
    return lim_inf_checked(a_n, n, context);
}
}
