#include "internal/integration_support.hpp"
#include "internal/assumption_facts.hpp"

namespace LMCAS {
namespace {
bool is_zero(const SymbolicExpr& expression) {
    const auto value = detail::exact_rational_value(detail::node(expression));
    return value && value->is_zero();
}
}

namespace detail {
SymbolicExpr materialize_integral(
    const IntegralParts& parts, const std::string& variable, RewriteBudget& budget) {
    const auto* primitive = dynamic_cast<const PiecewiseNode*>(node(parts.primitive).get());
    const auto* remainder = dynamic_cast<const PiecewiseNode*>(node(parts.remainder).get());
    if (primitive && remainder) {
        std::vector<PiecewiseNode::Branch> branches;
        branches.reserve(primitive->branches().size());
        for (std::size_t i = 0; i < primitive->branches().size(); ++i) {
            budget.consume();
            auto branch = materialize_integral({
                expression_from_node(primitive->branches()[i].expression),
                expression_from_node(remainder->branches()[i].expression)}, variable, budget);
            branches.push_back({node(branch), primitive->branches()[i].condition});
        }
        SymbolicNodePtr fallback;
        if (primitive->default_expr()) {
            fallback = node(materialize_integral({expression_from_node(primitive->default_expr()),
                expression_from_node(remainder->default_expr())}, variable, budget));
        }
        auto result = make_node<PiecewiseNode>(std::move(branches), fallback);
        budget.measure(result);
        return expression_from_node(result);
    }
    if (is_zero(parts.remainder)) { return parts.primitive; }
    auto integral = make_node<IntegralNode>(node(parts.remainder), variable);
    if (is_zero(parts.primitive)) { return expression_from_node(integral); }
    return expression_from_node(SymbolicFactory::create_add(
        {node(parts.primitive), integral}, &budget));
}
}

Integrator::Integrator() : entries_(detail::make_integration_entries()),
    zero_(*SymbolicExpr::number(0)) {}
Integrator::~Integrator() = default;
Integrator::Integrator(Integrator&&) noexcept = default;
Integrator& Integrator::operator=(Integrator&&) noexcept = default;


std::shared_ptr<SymbolicExpr> Integrator::check_cycle(
    const SymbolicExpr& expr, const std::string& var, ComputationContext& context) {
    detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
        context.limits().max_ast_nodes, "integrate.cycle");
    for (const auto& entry : cycle_state_.history) {
        budget.consume();
        if (detail::node(expr)->equals(*detail::node(entry.integrand))) {
            return SymbolicExpr::variable(entry.placeholder);
        }
        auto ratio = detail::make_expression_ptr(detail::normalize_node(
            detail::node(SymbolicExpr::divide(detail::make_expression_ptr(expr),
                detail::make_expression_ptr(entry.integrand))), budget));
        if (!depends_on_integration_variable(*ratio, var)) {
            return SymbolicExpr::multiply(ratio, SymbolicExpr::variable(entry.placeholder));
        }
    }
    return nullptr;
}

bool Integrator::has_cycle_placeholder(const SymbolicExpr& expression) const {
    for (const auto& entry : cycle_state_.history) {
        if (depends_on_integration_variable(expression, entry.placeholder)) { return true; }
    }
    return false;
}

bool Integrator::resolve_cycle(detail::IntegralParts& result,
    const std::string& placeholder, const std::string& var, ComputationContext& context) {
    if (depends_on_integration_variable(result.remainder, placeholder)) { return false; }
    if (!depends_on_integration_variable(result.primitive, placeholder)) { return true; }
    detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
        context.limits().max_ast_nodes, "integrate.cycle");
    auto B = detail::normalize_node(detail::differentiate_node(
        detail::node(result.primitive), placeholder, budget), budget);
    if (expression_depends_on_variable(B, placeholder) ||
        expression_depends_on_variable(B, var)) { return false; }
    auto A = detail::substitute_node(detail::node(result.primitive), placeholder,
        detail::node(zero_), budget);
    auto denominator = detail::normalize_node(detail::node(sym_sub(
        *SymbolicExpr::number(1), detail::expression_from_node(B))), budget);
    auto exact = detail::exact_rational_value(denominator);
    bool invertible = exact && !exact->is_zero();
    if (!exact && context.assumptions()) {
        detail::AssumptionFacts facts(*context.assumptions());
        auto nonzero = facts.is_nonzero(denominator, context);
        if (!nonzero) { throw nonzero.error(); }
        invertible = nonzero.value() == Tribool::True;
    }
    if (!invertible) { return false; }
    auto divisor = detail::make_expression_ptr(denominator);
    result.primitive = detail::expression_from_node(detail::normalize_node(detail::node(
        SymbolicExpr::divide(detail::make_expression_ptr(A), divisor)), budget));
    result.remainder = detail::expression_from_node(detail::normalize_node(detail::node(
        SymbolicExpr::divide(detail::make_expression_ptr(result.remainder), divisor)), budget));
    return !depends_on_integration_variable(result.primitive, placeholder) &&
        !depends_on_integration_variable(result.remainder, placeholder);
}

Result<std::optional<detail::IntegralParts>> Integrator::dispatch_builtins(
    const SymbolicExpr& expr, const std::string& var, ComputationContext& context, int depth) {
    using Attempt = Result<std::optional<detail::IntegralParts>>;
    using Method = Attempt (Integrator::*)(const SymbolicExpr&, const std::string&,
        ComputationContext&, int);
    constexpr Method methods[] = {
        &Integrator::try_power_rule, &Integrator::try_table_lookup,
        &Integrator::try_linear_substitution, &Integrator::try_substitution,
        &Integrator::try_trig_substitution, &Integrator::try_trig_combination,
        &Integrator::try_weierstrass, &Integrator::try_rational_decomposition,
        &Integrator::try_special_function, &Integrator::try_integration_by_parts};
    detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
        context.limits().max_ast_nodes, "integrate.builtins");
    const auto placeholder = cycle_state_.history.back().placeholder;
    std::optional<detail::IntegralParts> partial;
    for (auto method : methods) {
        auto access = context.consume_steps(1, "integrate.builtins");
        if (!access) { return Attempt::failure(access.error()); }
        auto attempted = (this->*method)(expr, var, context, depth);
        if (!attempted) { return attempted; }
        access = context.consume_steps(0, "integrate.builtins");
        if (!access) { return Attempt::failure(access.error()); }
        if (!attempted.value()) { continue; }
        auto candidate = std::move(*attempted.value());
        if (!resolve_cycle(candidate, placeholder, var, context)) { continue; }
        candidate.primitive = detail::expression_from_node(detail::normalize_node(
            detail::node(candidate.primitive), budget));
        candidate.remainder = detail::expression_from_node(detail::normalize_node(
            detail::node(candidate.remainder), budget));
        if (is_zero(candidate.remainder)) { return Attempt::success(std::move(candidate)); }
        if (!partial && (!is_zero(candidate.primitive) ||
            !detail::node(candidate.remainder)->equals(*detail::node(expr)))) {
            partial = std::move(candidate);
        }
    }
    return Attempt::success(std::move(partial));
}

Result<detail::IntegralParts> Integrator::integrate_recursive(
    const SymbolicExpr& source, const std::string& var,
    ComputationContext& context, int depth, const SymbolicExpr* original) {
    using PartsResult = Result<detail::IntegralParts>;
    auto entered = context.enter_recursion("integrate.recursive");
    if (!entered) { return PartsResult::failure(entered.error()); }
    struct RecursionExit {
        ComputationContext& context;
        ~RecursionExit() { context.leave_recursion(); }
    } recursion_exit{context};
    detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
        context.limits().max_ast_nodes, "integrate.recursive");
    budget.measure(detail::node(source));
    const auto& integrand = original ? *original : source;
    if (const auto* piecewise = dynamic_cast<const PiecewiseNode*>(detail::node(source).get())) {
        std::vector<PiecewiseNode::Branch> primitives, remainders;
        for (const auto& branch : piecewise->branches()) {
            auto term = detail::expression_from_node(branch.expression);
            auto parts = integrate_recursive(term, var, context, depth + 1);
            if (!parts) { return parts; }
            primitives.push_back({detail::node(parts.value().primitive), branch.condition});
            remainders.push_back({detail::node(parts.value().remainder), branch.condition});
        }
        detail::SymbolicNodePtr primitive_default, remainder_default;
        if (piecewise->default_expr()) {
            auto term = detail::expression_from_node(piecewise->default_expr());
            auto parts = integrate_recursive(term, var, context, depth + 1);
            if (!parts) { return parts; }
            primitive_default = detail::node(parts.value().primitive);
            remainder_default = detail::node(parts.value().remainder);
        }
        return PartsResult::success({detail::expression_from_node(detail::make_node<PiecewiseNode>(
            std::move(primitives), primitive_default)),
            detail::expression_from_node(detail::make_node<PiecewiseNode>(
                std::move(remainders), remainder_default))});
    }
    auto linear = apply_linearity(source, var, context, depth);
    if (!linear) { return PartsResult::failure(linear.error()); }
    if (linear.value()) { return PartsResult::success(std::move(*linear.value())); }
    const auto& working = detail::node(source);
    std::optional<detail::AssumptionFacts> facts;
    if (context.assumptions()) { facts.emplace(*context.assumptions()); }
    auto expr = detail::expression_from_node(detail::normalize_node(working, budget,
        facts ? static_cast<const FactsQuery&>(*facts) : detail::no_facts()));
    if (is_zero(expr)) { return PartsResult::success({zero_, zero_}); }
    if (!depends_on_integration_variable(expr, var)) {
        return PartsResult::success({*SymbolicExpr::multiply(
            detail::make_expression_ptr(expr), SymbolicExpr::variable(var)), zero_});
    }
    auto cycle = check_cycle(expr, var, context);
    if (cycle) { return PartsResult::success({*cycle, zero_}); }
    std::string placeholder;
    do {
        budget.consume();
        placeholder = "INT_CYCLE_" + std::to_string(cycle_state_.next_placeholder++);
    } while (cycle_state_.occupied_names.count(placeholder));
    cycle_state_.occupied_names.insert(placeholder);
    cycle_state_.history.push_back({expr, integrand, placeholder});
    struct CycleExit {
        CycleState& state;
        ~CycleExit() { state.history.pop_back(); }
    } cycle_exit{cycle_state_};
    auto result = dispatch_builtins(expr, var, context, depth);
    if (!result) { return PartsResult::failure(result.error()); }
    if (result.value()) {
        auto parts = std::move(*result.value());
        if (!depends_on_integration_variable(parts.primitive, placeholder) &&
            !depends_on_integration_variable(parts.remainder, placeholder)) {
            return PartsResult::success(std::move(parts));
        }
    }
    return PartsResult::success({zero_, integrand});
}

Result<SymbolicExpr> Integrator::integrate_checked(
    const SymbolicExpr& expr, const std::string& var_name, ComputationContext& context) try {
    if (var_name.empty()) {
        return Result<SymbolicExpr>::failure(CasErrc::InvalidArgument,
            "integration variable must not be empty", "integrate");
    }
    auto access = context.consume_steps(1, "integrate");
    if (!access) { return Result<SymbolicExpr>::failure(access.error()); }
    detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
        context.limits().max_ast_nodes, "integrate");
    budget.measure(detail::node(expr));
    const bool top_level = query_depth_ == 0;
    ++query_depth_;
    struct QueryExit {
        Integrator& integrator;
        ~QueryExit() {
            if (--integrator.query_depth_ == 0) {
                integrator.cycle_state_.history.clear();
                integrator.cycle_state_.occupied_names.clear();
            }
        }
    } query_exit{*this};
    if (top_level) {
        cycle_state_.history.clear();
        cycle_state_.next_placeholder = 0;
        cycle_state_.occupied_names = detail::all_variable_names(detail::node(expr), &budget);
        cycle_state_.occupied_names.insert(var_name);
        if (context.assumptions()) {
            detail::AssumptionFacts::collect_names(*context.assumptions(),
                cycle_state_.occupied_names, budget);
        }
    }
    auto recursive = integrate_recursive(expr, var_name, context, 0, &expr);
    if (!recursive) { return Result<SymbolicExpr>::failure(recursive.error()); }
    auto parts = std::move(recursive.value());
    if (has_cycle_placeholder(parts.primitive) || has_cycle_placeholder(parts.remainder)) {
        parts = {zero_, expr};
    }
    auto result = detail::materialize_integral(parts, var_name, budget);
    budget.measure(detail::node(result));
    access = context.consume_steps(0, "integrate");
    if (!access) { return Result<SymbolicExpr>::failure(access.error()); }
    return Result<SymbolicExpr>::success(std::move(result));
} catch (const CasError& error) {
    return Result<SymbolicExpr>::failure(error);
} catch (const std::bad_alloc&) {
    return Result<SymbolicExpr>::failure(CasErrc::ResourceLimit,
        "integration allocation failed", "integrate");
} catch (const std::length_error& error) {
    return Result<SymbolicExpr>::failure(CasErrc::ResourceLimit, error.what(), "integrate");
} catch (const std::exception& error) {
    return Result<SymbolicExpr>::failure(CasErrc::InternalInvariant, error.what(), "integrate");
}

}
