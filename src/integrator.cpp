#include "internal/integration_support.hpp"
#include "residual_verification.hpp"

namespace LMCAS {

IntegrationStrategyResult IntegrationStrategy::try_integrate(
    const SymbolicExpr& expr,
    const std::string& var,
    Integrator& integrator,
    ComputationContext& computation,
    int depth) {
    const std::string operation = "integrate.strategy." + name();
    auto access = computation.consume_steps(0, operation);
    if (!access) { return IntegrationStrategyResult::failure(access.error()); }
    try {
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
            std::move(expression), name()});
    } catch (const std::bad_alloc&) {
        return IntegrationStrategyResult::failure(
            CasErrc::ResourceLimit,
            "integration strategy allocation failed", operation);
    } catch (const std::exception& error) {
        return IntegrationStrategyResult::failure(
            CasErrc::InternalInvariant, error.what(), operation);
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
    strategies_.push_back(std::make_unique<PartialFractionStrategy>());
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
    const SymbolicExpr& expr, const std::string& var) {
    for (size_t i = 0; i < cycle_state_.history.size(); ++i) {
        auto ratio = SymbolicExpr::divide(detail::make_expression_ptr(expr), detail::make_expression_ptr(cycle_state_.history[i]))->simplify();
        if (!depends_on_integration_variable(*ratio, var)) {
            return SymbolicExpr::multiply(ratio, SymbolicExpr::variable("INT_CYCLE_" + std::to_string(i)));
        }
    }
    return nullptr;
}

void Integrator::resolve_cycle(std::shared_ptr<SymbolicExpr>& result, size_t cycle_idx) {
    std::string cycle_var = "INT_CYCLE_" + std::to_string(cycle_idx);
    if (depends_on_integration_variable(*result, cycle_var)) {
        auto B = result->differentiate(cycle_var)->simplify();
        auto A = result->substitute(cycle_var, SymbolicExpr::number(0))->simplify();
        auto one_minus_B = sym_sub(*SymbolicExpr::number(1), *B)->simplify();
        if (!one_minus_B->is_zero()) {
            result = SymbolicExpr::divide(A, one_minus_B);
        }
    }
}

Result<std::shared_ptr<SymbolicExpr>> Integrator::dispatch_strategies(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& context, int depth) {
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
            if (!strategy->requires_residual_verification()) {
                return Result<std::shared_ptr<SymbolicExpr>>::success(
                    candidate->expression);
            }
            auto derivative = candidate->expression
                ? candidate->expression->differentiate(var) : nullptr;
            if (!derivative) {
                if (depth > 0) {
                    return Result<std::shared_ptr<SymbolicExpr>>::success(
                        nullptr);
                }
                continue;
            }
            auto normalized_derivative = derivative->simplify();
            auto normalized_integrand = detail::make_expression_ptr(expr)->simplify();
            if (normalized_derivative && normalized_integrand &&
                LMCAS::detail::node(normalized_derivative)->equals(
                    *LMCAS::detail::node(normalized_integrand))) {
                return Result<std::shared_ptr<SymbolicExpr>>::success(
                    candidate->expression);
            }
            auto delta = SymbolicExpr::add(
                derivative,
                SymbolicExpr::multiply(
                    SymbolicExpr::number(-1), detail::make_expression_ptr(expr)));
            auto proof = check_zero_residual(delta, context);
            if (!proof) return Result<std::shared_ptr<SymbolicExpr>>::failure(proof.error());
            if (std::holds_alternative<ProvedZeroResidual>(proof.value())) {
                return Result<std::shared_ptr<SymbolicExpr>>::success(
                    candidate->expression);
            }
            if (depth > 0) {
                return Result<std::shared_ptr<SymbolicExpr>>::success(nullptr);
            }
        }
    }
    return Result<std::shared_ptr<SymbolicExpr>>::success(nullptr);
}

Result<std::shared_ptr<SymbolicExpr>> Integrator::integrate_recursive(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& context, int depth) {
    auto entered = context.enter_recursion("integrate.recursive");
    if (!entered) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(
            entered.error());
    }
    struct RecursionExit {
        ComputationContext& context;
        ~RecursionExit() { context.leave_recursion(); }
    } recursion_exit{context};

    auto linear_result = apply_linearity(expr, var, context, depth);
    if (!linear_result) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(
            linear_result.error());
    }
    if (linear_result.value()) { return linear_result; }

    auto cycle_result = check_cycle(expr, var);
    if (cycle_result) {
        return Result<std::shared_ptr<SymbolicExpr>>::success(
            std::move(cycle_result));
    }

    cycle_state_.history.push_back(expr);
    const size_t my_idx = cycle_state_.history.size() - 1;

    if (expr.is_number() || !depends_on_integration_variable(expr, var)) {
        cycle_state_.history.pop_back();
        return Result<std::shared_ptr<SymbolicExpr>>::success(
            SymbolicExpr::multiply(
                detail::make_expression_ptr(expr), SymbolicExpr::variable(var)));
    }

    auto result = dispatch_strategies(expr, var, context, depth);
    cycle_state_.history.pop_back();
    if (!result) {
        return Result<std::shared_ptr<SymbolicExpr>>::failure(result.error());
    }
    if (!result.value()) {
        return Result<std::shared_ptr<SymbolicExpr>>::success(
            make_integral_node(expr, var));
    }

    resolve_cycle(result.value(), my_idx);
    return result;
}

Result<SymbolicExpr> Integrator::integrate_checked(
    const SymbolicExpr& expr,
    const std::string& var_name,
    ComputationContext& context) {

    if (var_name.empty()) {
        return Result<SymbolicExpr>::failure(
            CasErrc::InvalidArgument, "integration variable must not be empty",
            "integrate");
    }
    auto access = context.consume_steps(1, "integrate");
    if (!access) { return Result<SymbolicExpr>::failure(access.error()); }

    const bool top_level_query = query_depth_ == 0;
    if (top_level_query) cycle_state_.history.clear();
    ++query_depth_;
    struct QueryExit {
        Integrator& integrator;
        ~QueryExit() { --integrator.query_depth_; }
    } query_exit{*this};

    if (auto pw = std::dynamic_pointer_cast<const PiecewiseNode>(LMCAS::detail::node(expr))) {
        std::vector<PiecewiseNode::Branch> new_brs;
        for (const auto& br : pw->branches()) {
            PiecewiseNode::Branch nb;
            auto branch = integrate_checked(
                LMCAS::detail::expression_from_node(br.expression), var_name, context);
            if (!branch) { return branch; }
            nb.expression = LMCAS::detail::node(branch.value());
            nb.condition = br.condition;
            new_brs.push_back(nb);
        }
        std::shared_ptr<const SymbolicNode> new_def = nullptr;
        if (pw->default_expr()) {
            auto default_result = integrate_checked(
                LMCAS::detail::expression_from_node(pw->default_expr()), var_name, context);
            if (!default_result) { return default_result; }
            new_def = LMCAS::detail::node(default_result.value());
        }
        return Result<SymbolicExpr>::success(LMCAS::detail::expression_from_node(
            LMCAS::detail::make_node<PiecewiseNode>(std::move(new_brs), new_def)));
    }

    SymbolicExpr working_expr = expr;
    try {
        auto assumed = apply_assumption_simplifications(
            expr, var_name, context.assumptions().get());
        if (!assumed) { return assumed; }
        working_expr = std::move(assumed.value());
        auto normalized = working_expr.simplify();
        if (!normalized) {
            return Result<SymbolicExpr>::failure(
                CasErrc::InternalInvariant,
                "integrand normalization produced no expression",
                "integrate.preprocess");
        }
        working_expr = *normalized;
    } catch (const std::exception& error) {
        return Result<SymbolicExpr>::failure(
            CasErrc::InternalInvariant, error.what(), "integrate.preprocess");
    }

    auto recursive = integrate_recursive(
        working_expr, var_name, context, 0);
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
    return Result<SymbolicExpr>::success(*result);
}

Result<SymbolicExpr> Integrator::integrate(
    const SymbolicExpr& expr,
    const std::string& var_name) {
    ComputationContext context;
    return integrate_checked(expr, var_name, context);
}


}
