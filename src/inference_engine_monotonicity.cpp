#define _USE_MATH_DEFINES
#include "internal/inference_engine_impl.hpp"

namespace lamina {

// Periodicity inference


InferenceTriboolResult InferenceEngine::query_periodic_checked(const SymbolicExpr& expr) const {
    return checked_inference_result<Tribool>(expr, "query_periodic_checked",
        [&]() {
            // FunctionNode: sin, cos, tan are periodic
            if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
                switch (func->type()) {
                    case FunctionNode::FuncType::Sin:
                    case FunctionNode::FuncType::Cos:
                    case FunctionNode::FuncType::Tan:
                        return Tribool::True;
                    case FunctionNode::FuncType::ArcTan: {
                        if (func->arguments().size() != 1) {
                            return Tribool::Unknown;
                        }
                        auto argument = lamina::detail::expression_from_node(
                            func->arguments()[0]);
                        if (auto variable =
                                std::dynamic_pointer_cast<const VariableNode>(
                                    func->arguments()[0])) {
                            const auto& properties =
                                impl_->ctx.current_properties();
                            return properties.is_periodic(variable->name())
                                ? Tribool::True
                                : Tribool::False;
                        }
                        auto periodic = query_periodic_checked(argument);
                        if (periodic &&
                            periodic.value() == Tribool::True) {
                            return Tribool::True;
                        }
                        return Tribool::Unknown;
                    }
                    default:
                        break;
                }
            }

            // VariableNode: check PropertyStore for declared periodicity
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                const auto& props = impl_->ctx.current_properties();
                if (props.is_periodic(var->name())) return Tribool::True;
            }
            return Tribool::Unknown;
        });
}


InferencePeriodResult InferenceEngine::infer_period_checked(const SymbolicExpr& expr) const {
    return checked_inference_result<std::optional<SymbolicExpr>>(
        expr, "infer_period_checked", [&]() -> std::optional<SymbolicExpr> {
            // FunctionNode: known periods for trig functions
            if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
                switch (func->type()) {
                    case FunctionNode::FuncType::Sin:
                    case FunctionNode::FuncType::Cos: {
                        // Period = 2*pi
                        auto two = lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(2.0));
                        auto pi_val = lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(LMMC_CONST_PI));
                        auto two_pi = lamina::detail::make_node<MultiplyNode>(
                            std::vector<std::shared_ptr<const SymbolicNode>>{two, pi_val});
                        auto period = lamina::detail::expression_from_node(two_pi);
                        return period;
                    }
                    case FunctionNode::FuncType::Tan: {
                        // Period = pi
                        auto pi_val = lamina::detail::make_node<NumberNode>(static_cast<lmmc_real_t>(LMMC_CONST_PI));
                        auto period = lamina::detail::expression_from_node(pi_val);
                        return period;
                    }
                    default:
                        break;
                }
            }

            // VariableNode: check PropertyStore for declared period
            if (auto var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
                const auto& props = impl_->ctx.current_properties();
                auto stored_period = props.get_period(var->name());
                if (stored_period.has_value() && *stored_period) {
                    return *(*stored_period); // dereference optional<shared_ptr<SymbolicExpr>>
                }
            }
            return std::nullopt;
        });
}

// Monotonicity inference

/**
 * @brief Helper: check if a MultiplyNode represents negation (multiplication by -1).
 * Returns the inner expression if it's a negation, nullptr otherwise.
 */
static std::shared_ptr<const SymbolicNode> detect_negation(const std::shared_ptr<const SymbolicNode>& node) {
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!mul || mul->operands().size() != 2) return nullptr;

    // Check if one operand is -1
    for (size_t i = 0; i < 2; ++i) {
        auto num = std::dynamic_pointer_cast<const NumberNode>(mul->operands()[i]);
        if (!num) continue;

        bool is_neg_one = false;
        if (std::holds_alternative<BigInt>(num->value())) {
            is_neg_one = (std::get<BigInt>(num->value()) == BigInt(-1));
        } else if (std::holds_alternative<Rational>(num->value())) {
            is_neg_one = (std::get<Rational>(num->value()) == Rational(-1));
        } else if (std::holds_alternative<lmmc_real_t>(num->value())) {
            is_neg_one = (std::get<lmmc_real_t>(num->value()) == -1.0);
        }

        if (is_neg_one) {
            return mul->operands()[1 - i]; // return the other operand
        }
    }
    return nullptr;
}

Monotonicity InferenceEngine::infer_monotonicity(const SymbolicExpr& expr,
                                                  const std::string& var,
                                                  const Interval& interval) const {
    if (!lamina::detail::node(expr)) return Monotonicity::Unknown;

    // FunctionNode: auto-infer for known functions
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
        // Check that the function argument is the queried variable
        if (func->arguments().empty()) return Monotonicity::Unknown;
        auto arg_var = std::dynamic_pointer_cast<const VariableNode>(func->arguments()[0]);
        if (!arg_var || arg_var->name() != var) return Monotonicity::Unknown;

        switch (func->type()) {
            case FunctionNode::FuncType::Exp:
                // exp is strictly increasing on all of ℝ
                return Monotonicity::Increasing;

            case FunctionNode::FuncType::ArcTan:
                // atan 在整个实数域严格递增.
                return Monotonicity::Increasing;

            case FunctionNode::FuncType::Ln: {
                // ln is strictly increasing on ℝ⁺ (positive reals)
                // Check if the interval is within positive reals
                // If the lower bound is > 0 (or open at 0), ln is increasing
                auto lo = endpoint_to_double(interval.lower);
                if (lo.has_value() && *lo > 0.0) {
                    return Monotonicity::Increasing;
                }
                // If lower is open at 0, also increasing
                if (lo.has_value() && *lo == 0.0 && interval.lower.is_open) {
                    return Monotonicity::Increasing;
                }
                // If interval is (0, +inf) or similar positive interval
                if (interval.lower.is_neg_infinity) {
                    return Monotonicity::Unknown; // ln not defined on negative reals
                }
                return Monotonicity::Unknown;
            }

            default:
                break;
        }
    }

    // Negation: multiply by -1 reverses monotonicity
    if (auto inner = detect_negation(lamina::detail::node(expr))) {
        auto inner_expr = lamina::detail::expression_from_node(inner);
        Monotonicity inner_mono = infer_monotonicity(inner_expr, var, interval);
        switch (inner_mono) {
            case Monotonicity::Increasing:    return Monotonicity::Decreasing;
            case Monotonicity::Decreasing:    return Monotonicity::Increasing;
            case Monotonicity::NonDecreasing:  return Monotonicity::NonIncreasing;
            case Monotonicity::NonIncreasing:  return Monotonicity::NonDecreasing;
            default:                           return Monotonicity::Unknown;
        }
    }

    // VariableNode: check PropertyStore for declared monotonicity
    if (auto var_node = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(expr))) {
        const auto& props = impl_->ctx.current_properties();
        return detail::propagate_result(
            props.get_monotonicity(var_node->name(), var, interval));
    }

    return Monotonicity::Unknown;
}

// Monotonicity deduction (apply_monotonicity_rules)

/**
 * @brief Helper: collect positive integer exponents that appear in PowerNode
 * expressions within the RelationStore involving the given variable names.
 */
static std::vector<int> collect_power_exponents(const RelationStore& store,
                                                const std::string& lhs_name,
                                                const std::string& rhs_name) {
    std::unordered_set<int> exponents;

    auto check_node = [&](const std::shared_ptr<const SymbolicNode>& node) {
        if (!node) return;
        auto pow_node = std::dynamic_pointer_cast<const PowerNode>(node);
        if (!pow_node) return;

        // Check if the base is one of our variables
        auto base_var = std::dynamic_pointer_cast<const VariableNode>(pow_node->base());
        if (!base_var) return;
        if (base_var->name() != lhs_name && base_var->name() != rhs_name) return;

        // Check if the exponent is a positive integer
        auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow_node->exponent());
        if (!exp_num) return;
        if (!is_positive_integer_number(*exp_num)) return;

        // Extract the integer value
        int n = 0;
        if (std::holds_alternative<BigInt>(exp_num->value())) {
            const auto& b = std::get<BigInt>(exp_num->value());
            // Only consider small exponents to avoid explosion
            if (b > BigInt(100)) return;
            n = static_cast<int>(b.IsNegative() ? 0 : b.is_zero() ? 0 : std::stoi(b.to_string()));
        } else if (std::holds_alternative<Rational>(exp_num->value())) {
            BigInt num_val = std::get<Rational>(exp_num->value()).get_numerator();
            if (num_val > BigInt(100)) return;
            n = static_cast<int>(std::stoi(num_val.to_string()));
        } else {
            double v = std::get<lmmc_real_t>(exp_num->value());
            if (v <= 0.0 || v > 100.0 || v != std::floor(v)) return;
            n = static_cast<int>(v);
        }

        if (n > 0) {
            exponents.insert(n);
        }
    };

    // Scan all relations in the store for PowerNode expressions
    for (const auto& rel : store.get_relations()) {
        check_node(lamina::detail::node(rel.lhs));
        check_node(lamina::detail::node(rel.rhs));
    }

    return std::vector<int>(exponents.begin(), exponents.end());
}

void InferenceEngine::apply_monotonicity_rules(const Relation& rel, RelationStore& store,
                                               PropertyStore& prop_store, int depth) {
    (void)apply_monotonicity_rules_checked(rel, store, prop_store, depth);
}

Result<void> InferenceEngine::apply_monotonicity_rules_checked(
    const Relation& rel,
    RelationStore& store,
    PropertyStore& prop_store,
    int depth) {
    // Stop recursion at maximum depth
    if (depth >= MAX_MONOTONICITY_DEPTH) return Result<void>::success();

    // Only apply to GT (greater-than) relations
    if (rel.op != RelationalNode::Op::GT) return Result<void>::success();

    // Extract LHS and RHS — both must be single VariableNodes
    auto lhs_var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(rel.lhs));
    auto rhs_var = std::dynamic_pointer_cast<const VariableNode>(lamina::detail::node(rel.rhs));
    if (!lhs_var || !rhs_var) return Result<void>::success();

    const std::string& x_name = lhs_var->name();
    const std::string& y_name = rhs_var->name();

    // Helper to add a deduced relation and recursively apply monotonicity rules
    auto add_deduced = [&](const SymbolicExpr& new_lhs,
                           const SymbolicExpr& new_rhs) -> Result<void> {
        // Don't add if already present
        if (store.has_relation(new_lhs, new_rhs, RelationalNode::Op::GT)) {
            return Result<void>::success();
        }

        auto inserted = store.add_relation_checked(
            new_lhs, new_rhs, RelationalNode::Op::GT, prop_store);
        if (!inserted.has_value()) {
            return Result<void>::failure(inserted.error());
        }

        // Recursively apply monotonicity rules to the newly added relation
        Relation new_rel{new_lhs, new_rhs, RelationalNode::Op::GT};
        return apply_monotonicity_rules_checked(new_rel, store, prop_store, depth + 1);
    };

    // Check domain assumptions using the AssumptionContext (read-through all scopes)
    bool both_positive = impl_->ctx.has_sign(x_name, Sign::Positive) &&
                         impl_->ctx.has_sign(y_name, Sign::Positive);
    bool both_real = impl_->ctx.has_domain(x_name, Domain::Real) &&
                     impl_->ctx.has_domain(y_name, Domain::Real);
    bool both_nonnegative = impl_->ctx.has_sign(x_name, Sign::NonNegative) &&
                            impl_->ctx.has_sign(y_name, Sign::NonNegative);

    if (both_positive) {
        auto ln_x_node = lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Ln,
            std::vector<std::shared_ptr<const SymbolicNode>>{lhs_var->clone()});
        auto ln_y_node = lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Ln,
            std::vector<std::shared_ptr<const SymbolicNode>>{rhs_var->clone()});

        auto ln_x = lamina::detail::expression_from_node(ln_x_node);
        auto ln_y = lamina::detail::expression_from_node(ln_y_node);
        auto deduced = add_deduced(ln_x, ln_y);
        if (!deduced.has_value()) return deduced;
    }

    if (both_positive) {
        auto sqrt_x_node = lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sqrt,
            std::vector<std::shared_ptr<const SymbolicNode>>{lhs_var->clone()});
        auto sqrt_y_node = lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Sqrt,
            std::vector<std::shared_ptr<const SymbolicNode>>{rhs_var->clone()});

        auto sqrt_x = lamina::detail::expression_from_node(sqrt_x_node);
        auto sqrt_y = lamina::detail::expression_from_node(sqrt_y_node);
        auto deduced = add_deduced(sqrt_x, sqrt_y);
        if (!deduced.has_value()) return deduced;
    }

    if (both_real) {
        auto exp_x_node = lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Exp,
            std::vector<std::shared_ptr<const SymbolicNode>>{lhs_var->clone()});
        auto exp_y_node = lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Exp,
            std::vector<std::shared_ptr<const SymbolicNode>>{rhs_var->clone()});

        auto exp_x = lamina::detail::expression_from_node(exp_x_node);
        auto exp_y = lamina::detail::expression_from_node(exp_y_node);
        auto deduced = add_deduced(exp_x, exp_y);
        if (!deduced.has_value()) return deduced;
    }

    if (both_nonnegative) {
        std::vector<int> exponents = collect_power_exponents(store, x_name, y_name);

        for (int n : exponents) {
            auto exp_node = lamina::detail::make_node<NumberNode>(BigInt(n));
            auto pow_x_node = lamina::detail::make_node<PowerNode>(lhs_var->clone(), exp_node->clone());
            auto pow_y_node = lamina::detail::make_node<PowerNode>(rhs_var->clone(), exp_node->clone());

            auto pow_x = lamina::detail::expression_from_node(pow_x_node);
            auto pow_y = lamina::detail::expression_from_node(pow_y_node);
            auto deduced = add_deduced(pow_x, pow_y);
            if (!deduced.has_value()) return deduced;
        }
    }

    return Result<void>::success();
}


} // namespace lamina
