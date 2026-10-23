#include "internal/expression_transform.hpp"
#include "internal/integration_support.hpp"
#include "residual_verification.hpp"
#include "internal/assumption_facts.hpp"

namespace LMCAS {

IntegrationStrategyResult IntegrationStrategy::try_integrate(
    const SymbolicExpr& expr,
    const std::string& var,
    Integrator& integrator,
    ComputationContext& computation,
    int depth) {
    std::string operation;
    try {
        auto strategy_name = name();
        operation = "integrate.strategy." + strategy_name;
        auto access = computation.consume_steps(0, operation);
        if (!access) { return IntegrationStrategyResult::failure(access.error()); }
        auto raw_result = try_integrate_raw(
            expr, var, integrator, computation, depth);
        if (!raw_result) {
            return IntegrationStrategyResult::failure(raw_result.error());
        }
        auto expression = std::move(raw_result.value());
        auto final_access = computation.consume_steps(0, operation);
        if (!final_access) {
            return IntegrationStrategyResult::failure(final_access.error());
        }
        if (!expression) {
            return IntegrationStrategyResult::success(
                IntegrationNotApplicable{});
        }
        return IntegrationStrategyResult::success(IntegrationCandidate{
            std::move(expression), std::move(strategy_name)});
    } catch (const CasError& error) {
        return IntegrationStrategyResult::failure(error);
    } catch (const std::bad_alloc&) {
        return IntegrationStrategyResult::failure(
            CasErrc::ResourceLimit,
            "integration strategy allocation failed",
            operation.empty() ? std::string("integrate.strategy") : std::move(operation));
    } catch (const std::length_error& error) {
        return IntegrationStrategyResult::failure(
            CasErrc::ResourceLimit, error.what(),
            operation.empty() ? std::string("integrate.strategy") : std::move(operation));
    } catch (const std::exception& error) {
        return IntegrationStrategyResult::failure(
            CasErrc::InternalInvariant, error.what(),
            operation.empty() ? std::string("integrate.strategy") : std::move(operation));
    }
}

Integrator::Integrator() {

    strategies_.push_back(std::make_unique<PowerRuleStrategy>());
    strategies_.push_back(std::make_unique<TableLookupStrategy>());
    strategies_.push_back(std::make_unique<LinearSubstitutionStrategy>());
    strategies_.push_back(std::make_unique<SubstitutionStrategy>());
    strategies_.push_back(std::make_unique<TrigSubstitutionStrategy>());
    strategies_.push_back(std::make_unique<TrigCombinationStrategy>());
    strategies_.push_back(std::make_unique<WeierstrassStrategy>());
    strategies_.push_back(std::make_unique<RationalDecompositionStrategy>());
    strategies_.push_back(std::make_unique<SpecialFunctionStrategy>());
    strategies_.push_back(std::make_unique<IBPStrategy>());
}

Result<void> Integrator::add_strategy(
    std::unique_ptr<IntegrationStrategy> strategy,
    int position) {
    if (!strategy) {
        return Result<void>::failure(
            CasErrc::InvalidArgument, "strategy must not be null", "integrator.add_strategy");
    }
    try {
        if (position < 0 || position >= static_cast<int>(strategies_.size())) {
            strategies_.push_back(std::move(strategy));
        } else {
            strategies_.insert(strategies_.begin() + position, std::move(strategy));
        }
    } catch (const std::bad_alloc&) {
        return Result<void>::failure(
            CasErrc::ResourceLimit, "strategy allocation failed", "integrator.add_strategy");
    }
    return Result<void>::success();
}

bool Integrator::depends_on(const SymbolicExpr& expr, const std::string& var) {
    return expression_depends_on_variable(LMCAS::detail::node(expr), var);
}

std::shared_ptr<SymbolicExpr> Integrator::make_integral_node(
    const SymbolicExpr& expr, const std::string& var) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<IntegralNode>(
            LMCAS::detail::node(expr), var));
}

std::shared_ptr<SymbolicExpr> Integrator::check_cycle(
    const SymbolicExpr& expr, const std::string& var, ComputationContext& context) {
    detail::RewriteBudget transform_budget(context, context.limits().max_recursion_depth,
                                            context.limits().max_ast_nodes, "integrate.cycle");
    if (contains_inexact_number(detail::node(expr))) { return nullptr; }
    for (const auto& entry : cycle_state_.history) {
        transform_budget.consume();
        if (detail::node(expr)->equals(*detail::node(entry.integrand))) {
            return SymbolicExpr::variable(entry.placeholder);
        }
        auto ratio = detail::make_expression_ptr(detail::normalize_node(
            detail::node(SymbolicExpr::divide(detail::make_expression_ptr(expr),
                detail::make_expression_ptr(entry.integrand))), transform_budget));
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

const SymbolicExpr& Integrator::original_integrand(const SymbolicExpr& expression) const {
    if (!cycle_state_.history.empty()) {
        const auto& entry = cycle_state_.history.back();
        if (detail::node(expression) == detail::node(entry.integrand)) { return entry.original; }
    }
    return expression;
}

bool Integrator::resolve_cycle(std::shared_ptr<SymbolicExpr>& result,
    const std::string& placeholder, ComputationContext& context) {
    if (!depends_on_integration_variable(*result, placeholder)) { return true; }
    detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
                                context.limits().max_ast_nodes, "integrate.cycle");
    auto B = detail::make_expression_ptr(detail::differentiate_node(
        detail::node(result), placeholder, budget));
    if (depends_on_integration_variable(*B, placeholder)) { return false; }
    auto A = detail::make_expression_ptr(detail::substitute_node(
        detail::node(result), placeholder, detail::node(SymbolicExpr::number(0)), budget));
    auto denominator = detail::make_expression_ptr(detail::normalize_node(
        detail::node(sym_sub(*SymbolicExpr::number(1), *B)), budget));
    auto exact = detail::exact_rational_value(detail::node(denominator));
    bool invertible = exact && *exact != Rational(0);
    if (!exact && context.assumptions()) {
        detail::AssumptionFacts facts(*context.assumptions());
        auto nonzero = facts.is_nonzero(detail::node(denominator), context);
        if (!nonzero) { throw nonzero.error(); }
        invertible = nonzero.value() == Tribool::True;
    }
    if (!invertible) { return false; }
    result = detail::make_expression_ptr(detail::normalize_node(
        detail::node(SymbolicExpr::divide(A, denominator)), budget));
    return !depends_on_integration_variable(*result, placeholder);
}

Result<std::shared_ptr<SymbolicExpr>> Integrator::dispatch_strategies(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& context, int depth, const SymbolicExpr& original) {
    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    const std::string placeholder = cycle_state_.history.back().placeholder;
    for (auto& strategy : strategies_) {
        auto budget = context.consume_steps(1, "integrate.strategy");
        if (!budget) {
            return Result<std::shared_ptr<SymbolicExpr>>::failure(
                budget.error());
        }
        auto attempt = strategy->try_integrate(
            expr, var, *this, context, depth);
        if (!attempt) {
            return Result<std::shared_ptr<SymbolicExpr>>::failure(
                attempt.error());
        }
        if (auto* candidate =
                std::get_if<IntegrationCandidate>(&attempt.value())) {
            auto expression = candidate->expression;
            if (!expression || contains_unevaluated_integral(detail::node(expression))) { continue; }
            if (!resolve_cycle(expression, placeholder, context)) { continue; }
            // An ancestor's obligation is private generation state, not a closed primitive.
            if (has_cycle_placeholder(*expression) ||
                certify_integral(*expression, original, var, transform_budget)) {
                return Result<std::shared_ptr<SymbolicExpr>>::success(std::move(expression));
            }
        }
    }
    return Result<std::shared_ptr<SymbolicExpr>>::success(nullptr);
}

Result<std::shared_ptr<SymbolicExpr>> Integrator::integrate_recursive(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& context, int depth, const SymbolicExpr* original) {
    auto entered = context.enter_recursion("integrate.recursive");
    if (!entered) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(
            entered.error());
    }
    struct RecursionExit {
        ComputationContext& context;
        ~RecursionExit() { context.leave_recursion(); }
    } recursion_exit{context};

    detail::RewriteBudget budget(context, context.limits().max_recursion_depth,
                                context.limits().max_ast_nodes, "integrate.recursive");
    budget.measure(detail::node(expr));
    const SymbolicExpr& integrand = original ? *original : expr;
    auto linear_result = apply_linearity(expr, var, context, depth);
    if (!linear_result) { return linear_result; }
    if (linear_result.value() &&
        !contains_unevaluated_integral(detail::node(linear_result.value())) &&
        (has_cycle_placeholder(*linear_result.value()) ||
         certify_integral(*linear_result.value(), integrand, var, budget))) {
        return linear_result;
    }

    auto cycle_result = check_cycle(expr, var, context);
    if (cycle_result) { return cycle_result; }

    if (expr.is_number() || !depends_on_integration_variable(expr, var)) {
        auto candidate = SymbolicExpr::multiply(
            detail::make_expression_ptr(expr), SymbolicExpr::variable(var));
        if (certify_integral(*candidate, integrand, var, budget)) { return candidate; }
    }

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
    auto result = dispatch_strategies(expr, var, context, depth, integrand);
    if (!result) { return result; }
    return result.value() ? result
        : Result<std::shared_ptr<SymbolicExpr>>::success(make_integral_node(integrand, var));
}

Result<SymbolicExpr> Integrator::integrate_checked(
    const SymbolicExpr& expr,
    const std::string& var_name,
    ComputationContext& context) try {

    if (var_name.empty()) {
        return Result<SymbolicExpr>::failure(
            CasErrc::InvalidArgument, "integration variable must not be empty",
            "integrate");
    }
    auto access = context.consume_steps(1, "integrate");
    if (!access) { return Result<SymbolicExpr>::failure(access.error()); }
    detail::RewriteBudget transform_budget(
        context, context.limits().max_recursion_depth,
        context.limits().max_ast_nodes, "integrate.preprocess");
    transform_budget.measure(detail::node(expr));
    if (contains_inexact_number(detail::node(expr))) {
        auto obligation = make_integral_node(expr, var_name);
        transform_budget.measure(detail::node(obligation));
        return Result<SymbolicExpr>::success(*obligation);
    }

    const bool top_level_query = query_depth_ == 0;
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
    if (top_level_query) {
        cycle_state_.history.clear();
        cycle_state_.next_placeholder = 0;
        cycle_state_.occupied_names = detail::all_variable_names(detail::node(expr), &transform_budget);
        cycle_state_.occupied_names.insert(var_name);
        if (context.assumptions()) {
            detail::AssumptionFacts::collect_names(
                *context.assumptions(), cycle_state_.occupied_names, transform_budget);
        }
    }

    if (auto pw = std::dynamic_pointer_cast<const PiecewiseNode>(LMCAS::detail::node(expr))) {
        std::size_t nodes = 1;
        std::vector<PiecewiseNode::Branch> new_brs;
        for (const auto& br : pw->branches()) {
            PiecewiseNode::Branch nb;
            auto branch = integrate_checked(
                LMCAS::detail::expression_from_node(br.expression), var_name, context);
            if (!branch) { return branch; }
            nb.expression = LMCAS::detail::node(branch.value());
            nb.condition = br.condition;
            nodes = transform_budget.append_size(nodes, transform_budget.measure(nb.expression));
            nodes = transform_budget.append_size(nodes, transform_budget.measure(nb.condition));
            new_brs.push_back(nb);
        }
        std::shared_ptr<const SymbolicNode> new_def = nullptr;
        if (pw->default_expr()) {
            auto default_result = integrate_checked(
                LMCAS::detail::expression_from_node(pw->default_expr()), var_name, context);
            if (!default_result) { return default_result; }
            new_def = LMCAS::detail::node(default_result.value());
            nodes = transform_budget.append_size(nodes, transform_budget.measure(new_def));
        }
        transform_budget.require_nodes(nodes);
        auto candidate = detail::make_expression_ptr(
            LMCAS::detail::make_node<PiecewiseNode>(std::move(new_brs), new_def));
        if (certify_integral(*candidate, expr, var_name, transform_budget)) {
            transform_budget.measure(detail::node(candidate));
            return Result<SymbolicExpr>::success(*candidate);
        }
    }

    auto working_expr = [&] {
        std::optional<detail::AssumptionFacts> facts;
        if (context.assumptions()) { facts.emplace(*context.assumptions()); }
        return detail::expression_from_node(detail::normalize_node(
            detail::node(expr), transform_budget,
            facts ? static_cast<const FactsQuery&>(*facts) : detail::no_facts()));
    }();

    auto recursive = integrate_recursive(
        working_expr, var_name, context, 0, &expr);
    if (!recursive) {
        return Result<SymbolicExpr>::failure(recursive.error());
    }
    auto result = std::move(recursive.value());
    auto final_access = context.consume_steps(0, "integrate");
    if (!final_access) { return Result<SymbolicExpr>::failure(final_access.error()); }
    if (!result) {
        return Result<SymbolicExpr>::failure(
            CasErrc::InternalInvariant, "integration produced no result", "integrate");
    }
    if (contains_unevaluated_integral(detail::node(result)) ||
        has_cycle_placeholder(*result) ||
        !certify_integral(*result, expr, var_name, transform_budget)) {
        result = make_integral_node(expr, var_name);
    }
    transform_budget.measure(detail::node(result));
    return Result<SymbolicExpr>::success(*result);
}
catch (const CasError& error) {
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
