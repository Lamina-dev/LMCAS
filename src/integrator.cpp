#include "internal/integration_support.hpp"

namespace LMCAS {

// Assumption-aware integrand simplification helpers

/**
 * @brief Recursively replace |var| with var in an AST when the variable is known Positive.
 *
 * Traverses the expression tree and replaces any FunctionNode(Abs, [arg]) where
 * arg is exactly the integration variable with just the variable itself.
 * This is valid when the AssumptionContext confirms the variable is Positive.
 */
static std::shared_ptr<const SymbolicNode> simplify_abs_positive(
    const std::shared_ptr<const SymbolicNode>& node, const std::string& var) {
    if (!node) return node;

    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (fn->type() == FunctionNode::FuncType::Abs && fn->arguments().size() == 1) {
            // Check if the argument is exactly the integration variable
            if (auto vn = std::dynamic_pointer_cast<const VariableNode>(fn->arguments()[0])) {
                if (vn->name() == var) {
                    // |x| → x when x is Positive
                    return fn->arguments()[0];
                }
            }
        }
        // Recurse into function arguments
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        bool changed = false;
        for (auto& op : fn->arguments()) {
            auto new_op = simplify_abs_positive(op, var);
            if (new_op != op) changed = true;
            new_ops.push_back(new_op);
        }
        if (changed) {
            return LMCAS::detail::make_node<FunctionNode>(fn->type(), new_ops);
        }
        return node;
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        bool changed = false;
        for (auto& op : add->operands()) {
            auto new_op = simplify_abs_positive(op, var);
            if (new_op != op) changed = true;
            new_ops.push_back(new_op);
        }
        if (changed) return LMCAS::detail::make_node<AddNode>(new_ops);
        return node;
    }

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        std::vector<std::shared_ptr<const SymbolicNode>> new_ops;
        bool changed = false;
        for (auto& op : mul->operands()) {
            auto new_op = simplify_abs_positive(op, var);
            if (new_op != op) changed = true;
            new_ops.push_back(new_op);
        }
        if (changed) return LMCAS::detail::make_node<MultiplyNode>(new_ops);
        return node;
    }

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto new_base = simplify_abs_positive(pow->base(), var);
        auto new_exp = simplify_abs_positive(pow->exponent(), var);
        if (new_base != pow->base() || new_exp != pow->exponent()) {
            return LMCAS::detail::make_node<PowerNode>(new_base, new_exp);
        }
        return node;
    }

    // Leaf nodes (NumberNode, VariableNode, etc.) — no change
    return node;
}

/**
 * @brief Apply assumption-aware simplifications to an integrand.
 *
 * When the AssumptionContext indicates the integration variable is Positive,
 * replaces |var| with var throughout the integrand. Returns the original
 * expression unchanged if no simplifications apply.
 */
static Result<SymbolicExpr> apply_assumption_simplifications(
    const SymbolicExpr& expr, const std::string& var,
    const AssumptionContext* ctx) {
    if (!ctx) return expr;

    SymbolicExpr var_expr = *SymbolicExpr::variable(var);
    auto positive = ctx->is_positive(var_expr);
    if (!positive) {
        return Result<SymbolicExpr>::failure(positive.error());
    }

    if (positive.value() == Tribool::True) {
        auto new_root = simplify_abs_positive(
            LMCAS::detail::node(expr), var);
        if (new_root != LMCAS::detail::node(expr)) {
            return LMCAS::detail::expression_from_node(new_root);
        }
    }

    return expr;
}

IntegrationStrategyResult IntegrationStrategy::try_integrate(
    const SymbolicExpr& expr,
    const std::string& var,
    Integrator& integrator,
    ComputationContext& computation,
    int depth) {
    const std::string operation = "integrate.strategy." + name();
    auto access = computation.consume_steps(0, operation);
    if (!access) return IntegrationStrategyResult::failure(access.error());
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
        auto ratio = SymbolicExpr::divide(make_expr_ptr(expr), make_expr_ptr(cycle_state_.history[i]))->simplify();
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

Result<std::shared_ptr<SymbolicExpr>> Integrator::apply_linearity(
    const SymbolicExpr& expr, const std::string& var,
    ComputationContext& context, int depth) {

    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        std::vector<std::shared_ptr<const SymbolicNode>> constants;
        std::vector<std::shared_ptr<const SymbolicNode>> dependents;

        for (auto& op : mul->operands()) {
            auto term = LMCAS::detail::expression_from_node(op);
            if (!depends_on_integration_variable(term, var)) {
                constants.push_back(op);
            } else {
                dependents.push_back(op);
            }
        }

        if (!constants.empty() && dependents.size() < mul->operands().size()) {
            auto const_part = (constants.size() == 1)
                ? LMCAS::detail::expression_from_node(constants[0])
                : LMCAS::detail::expression_from_node(
                      LMCAS::detail::make_node<MultiplyNode>(constants));
            SymbolicExpr dep_part = (dependents.empty()) ?
                *SymbolicExpr::number(1) :
                ((dependents.size() == 1)
                     ? LMCAS::detail::expression_from_node(dependents[0])
                     : LMCAS::detail::expression_from_node(
                           LMCAS::detail::make_node<MultiplyNode>(dependents)));
            if (dependents.size() > 1) {
                std::vector<std::shared_ptr<const SymbolicNode>> exponents;
                exponents.reserve(dependents.size());
                bool all_exponential = true;
                for (const auto& dependent : dependents) {
                    auto function =
                        std::dynamic_pointer_cast<const FunctionNode>(dependent);
                    if (!function ||
                        function->type() != FunctionNode::FuncType::Exp ||
                        function->arguments().size() != 1) {
                        all_exponential = false;
                        break;
                    }
                    exponents.push_back(function->arguments()[0]);
                }
                if (all_exponential) {
                    dep_part = LMCAS::detail::expression_from_node(
                        LMCAS::detail::make_node<FunctionNode>(
                            FunctionNode::FuncType::Exp,
                            std::vector<std::shared_ptr<const SymbolicNode>>{
                                SymbolicFactory::create_add(
                                    std::move(exponents))}));
                }
            }

            auto int_part = integrate_recursive(dep_part, var, context, depth + 1);
            if (!int_part) {
                return Result<std::shared_ptr<SymbolicExpr>>::failure(
                    int_part.error());
            }
            return Result<std::shared_ptr<SymbolicExpr>>::success(
                SymbolicExpr::multiply(
                    LMCAS::detail::make_expression_ptr(const_part),
                    int_part.value()));
        }
    }

    if (auto add = std::dynamic_pointer_cast<const AddNode>(
            LMCAS::detail::node(expr))) {
        std::vector<std::shared_ptr<const SymbolicNode>> results;
        for (auto& op : add->operands()) {
            auto term = LMCAS::detail::expression_from_node(op);
            auto int_term = integrate_recursive(
                term, var, context, depth + 1);
            if (!int_term) {
                return Result<std::shared_ptr<SymbolicExpr>>::failure(
                    int_term.error());
            }
            results.push_back(LMCAS::detail::node(int_term.value()));
        }
        return Result<std::shared_ptr<SymbolicExpr>>::success(
            LMCAS::detail::make_expression_ptr(
                LMCAS::detail::make_node<AddNode>(results)));
    }
    return Result<std::shared_ptr<SymbolicExpr>>::success(nullptr);
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
            auto normalized_integrand = make_expr_ptr(expr)->simplify();
            if (normalized_derivative && normalized_integrand &&
                LMCAS::detail::node(normalized_derivative)->equals(
                    *LMCAS::detail::node(normalized_integrand))) {
                return Result<std::shared_ptr<SymbolicExpr>>::success(
                    candidate->expression);
            }
            auto delta = SymbolicExpr::add(
                derivative,
                SymbolicExpr::multiply(
                    SymbolicExpr::number(-1), make_expr_ptr(expr)));
            delta = delta->cancel()->simplify();
            if (delta && delta->is_zero()) {
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
    if (linear_result.value()) return linear_result;

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
                make_expr_ptr(expr), SymbolicExpr::variable(var)));
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
    if (!access) return Result<SymbolicExpr>::failure(access.error());

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
            if (!branch) return branch;
            nb.expression = LMCAS::detail::node(branch.value());
            nb.condition = br.condition;
            new_brs.push_back(nb);
        }
        std::shared_ptr<const SymbolicNode> new_def = nullptr;
        if (pw->default_expr()) {
            auto default_result = integrate_checked(
                LMCAS::detail::expression_from_node(pw->default_expr()), var_name, context);
            if (!default_result) return default_result;
            new_def = LMCAS::detail::node(default_result.value());
        }
        return Result<SymbolicExpr>::success(LMCAS::detail::expression_from_node(
            LMCAS::detail::make_node<PiecewiseNode>(std::move(new_brs), new_def)));
    }

    SymbolicExpr working_expr = expr;
    try {
        auto assumed = apply_assumption_simplifications(
            expr, var_name, context.assumptions().get());
        if (!assumed) return assumed;
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
    if (!final_access) return Result<SymbolicExpr>::failure(final_access.error());
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

Result<SymbolicExpr> Integrator::integrate_def_checked(
    const SymbolicExpr& expr,
    const std::string& var_name,
    const SymbolicExpr& lower,
    const SymbolicExpr& upper,
    ComputationContext& context) {

    auto access = context.consume_steps(1, "integrate.definite");
    if (!access) return Result<SymbolicExpr>::failure(access.error());
    SymbolicExpr simp_expr_val = expr;
    try {
        auto assumed = apply_assumption_simplifications(
            expr, var_name, context.assumptions().get());
        if (!assumed) return assumed;
        auto simplified = assumed.value().simplify();
        if (!simplified) {
            return Result<SymbolicExpr>::failure(
                CasErrc::InternalInvariant,
                "definite integrand normalization produced no expression",
                "integrate.definite.preprocess");
        }
        simp_expr_val = *simplified;
    } catch (const std::exception& error) {
        return Result<SymbolicExpr>::failure(
            CasErrc::InternalInvariant, error.what(),
            "integrate.definite.preprocess");
    }
    bool is_inv_x = false;

    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(LMCAS::detail::node(simp_expr_val))) {
        if (auto v = std::dynamic_pointer_cast<const VariableNode>(pow->base())) {
            if (v->name() == var_name) {
                if (auto en = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
                    int eq_minus_one = 0;
                    if (std::holds_alternative<lmmc_real_t>(en->value())) {
                        lmmc_double_nearly_equal_tol(std::get<lmmc_real_t>(en->value()), -1.0, 1e-9, 1e-9, &eq_minus_one);
                    }
                    if ((std::holds_alternative<lmmc_real_t>(en->value()) && eq_minus_one != 0) ||
                        (std::holds_alternative<BigInt>(en->value()) && std::get<BigInt>(en->value()).to_int() == -1) ||
                        (std::holds_alternative<Rational>(en->value()) && std::get<Rational>(en->value()).to_double() == -1.0)) {
                        is_inv_x = true;
                    }
                }
            }
        }
    }

    bool numeric_bounds = (LMCAS::detail::node(lower) && std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(lower))) &&
                          (LMCAS::detail::node(upper) && std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(upper)));

    if (is_inv_x && numeric_bounds) {
        auto l_checked = try_checked_numeric_constant(lower, context);
        auto u_checked = try_checked_numeric_constant(upper, context);
        if (l_checked && u_checked) {
            double l_val = *l_checked;
            double u_val = *u_checked;
            int eq_l, eq_u;
            lmmc_double_nearly_equal_tol(l_val, 0.0, 1e-9, 1e-9, &eq_l);
            lmmc_double_nearly_equal_tol(u_val, 0.0, 1e-9, 1e-9, &eq_u);
            if (!eq_l && l_val < 0 && !eq_u && u_val > 0) {
                auto t = LMCAS::detail::make_expression_ptr(*SymbolicExpr::variable("t"));
                auto zero = LMCAS::detail::make_expression_ptr(*SymbolicExpr::number(0));
                auto int_left_result = integrate_def_checked(
                    expr, var_name, lower, *t, context);
                if (!int_left_result) return int_left_result;
                auto int_left = int_left_result.value();
                auto lim_left_result = limit_expression_checked(
                    std::make_shared<SymbolicExpr>(int_left), "t", zero,
                    LimitDirection::FromBelow, context);
                if (!lim_left_result) {
                    return Result<SymbolicExpr>::failure(
                        lim_left_result.error());
                }
                auto lim_left = lim_left_result.value();
                auto int_right_result = integrate_def_checked(
                    expr, var_name, *t, upper, context);
                if (!int_right_result) return int_right_result;
                auto int_right = int_right_result.value();
                auto lim_right_result = limit_expression_checked(
                    std::make_shared<SymbolicExpr>(int_right), "t", zero,
                    LimitDirection::FromAbove, context);
                if (!lim_right_result) {
                    return Result<SymbolicExpr>::failure(
                        lim_right_result.error());
                }
                auto lim_right = lim_right_result.value();
                if (lim_left && lim_right) {
                    return Result<SymbolicExpr>::success(
                        *SymbolicExpr::add(lim_left, lim_right));
                }
            }
        }
    }

    auto indefinite_result = integrate_checked(expr, var_name, context);
    if (!indefinite_result) return indefinite_result;
    SymbolicExpr indefinite = indefinite_result.value();

    if (std::dynamic_pointer_cast<const IntegralNode>(
            LMCAS::detail::node(indefinite))) {
        return Result<SymbolicExpr>::success(
            LMCAS::detail::expression_from_node(
                LMCAS::detail::make_node<IntegralNode>(
                    LMCAS::detail::node(expr), var_name,
                    LMCAS::detail::node(lower),
                    LMCAS::detail::node(upper))));
    }

    try {
        auto F_b = indefinite.substitute(var_name, make_expr_ptr(upper));
        auto F_a = indefinite.substitute(var_name, make_expr_ptr(lower));
        if (!F_b || !F_a) {
            return Result<SymbolicExpr>::failure(
                CasErrc::InternalInvariant,
                "definite integration substitution produced no expression",
                "integrate.definite.evaluate");
        }
        auto result = sym_sub(*F_b, *F_a);
        auto simplified = result ? result->simplify() : nullptr;
        if (!simplified) {
            return Result<SymbolicExpr>::failure(
                CasErrc::InternalInvariant,
                "definite integration simplification produced no expression",
                "integrate.definite.evaluate");
        }
        return Result<SymbolicExpr>::success(*simplified);
    } catch (const std::exception& error) {
        return Result<SymbolicExpr>::failure(
            CasErrc::InternalInvariant, error.what(),
            "integrate.definite.evaluate");
    }
}

Result<SymbolicExpr> Integrator::integrate_def(
    const SymbolicExpr& expr,
    const std::string& var_name,
    const SymbolicExpr& lower,
    const SymbolicExpr& upper) {
    ComputationContext context;
    return integrate_def_checked(expr, var_name, lower, upper, context);
}

} // namespace LMCAS
