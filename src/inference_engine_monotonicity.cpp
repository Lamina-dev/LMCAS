#define _USE_MATH_DEFINES
#include "internal/inference_engine_impl.hpp"
#include "internal/interval_endpoint.hpp"
#include "internal/normalization_utils.hpp"

namespace LMCAS {

Monotonicity InferenceEngine::Impl::logarithm_monotonicity(const Interval& interval) {
    ComputationContext context;
    auto lower = detail::comparable_endpoint(
        interval.lower, context, "inference.logarithm_monotonicity");
    if (!lower) return Monotonicity::Unknown;
    auto comparison = detail::compare_comparable(
        lower.value(), detail::ComparableEndpoint{}, context);
    if (!comparison) return Monotonicity::Unknown;
    if (comparison.value() > 0 ||
        (comparison.value() == 0 && interval.lower.is_open)) {
        return Monotonicity::Increasing;
    }
    return Monotonicity::Unknown;
}

Monotonicity InferenceEngine::Impl::function_monotonicity(
    const FunctionNode& func,
    const std::string& var,
    const Interval& interval) {
    if (func.arguments().empty()) {
        return Monotonicity::Unknown;
    }
    auto arg_var = std::dynamic_pointer_cast<const VariableNode>(func.arguments()[0]);
    if (!arg_var || arg_var->name() != var) {
        return Monotonicity::Unknown;
    }

    switch (func.type()) {
        case FunctionNode::FuncType::Exp:
            return Monotonicity::Increasing;

        case FunctionNode::FuncType::ArcTan:
            return Monotonicity::Increasing;

        case FunctionNode::FuncType::Ln: {
            return logarithm_monotonicity(interval);
        }

        default:
            break;
    }
    return Monotonicity::Unknown;
}

std::optional<int> InferenceEngine::Impl::small_positive_exponent(
    const std::shared_ptr<const NumberNode>& number) {
    BigInt exponent;
    if (!try_get_integer_value(number, exponent) ||
        exponent <= BigInt(0) || exponent > BigInt(100)) {
        return std::nullopt;
    }
    return static_cast<int>(*exponent.try_to_int64());
}

void InferenceEngine::Impl::collect_node_exponent(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& lhs_name,
    const std::string& rhs_name,
    std::unordered_set<int>& exponents) {
    if (!node) {
        return;
    }
    auto pow_node = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!pow_node) {
        return;
    }

    auto base_var = std::dynamic_pointer_cast<const VariableNode>(pow_node->base());
    if (!base_var) {
        return;
    }
    if (base_var->name() != lhs_name && base_var->name() != rhs_name) {
        return;
    }

    auto exp_num = std::dynamic_pointer_cast<const NumberNode>(pow_node->exponent());
    if (!exp_num) {
        return;
    }
    const auto exponent = small_positive_exponent(exp_num);
    if (exponent) {
        exponents.insert(*exponent);
    }
}

std::vector<int> InferenceEngine::Impl::collect_power_exponents(
    const RelationStore& store,
    const std::string& lhs_name,
    const std::string& rhs_name) {
    std::unordered_set<int> exponents;
    for (const auto& rel : store.get_relations()) {
        collect_node_exponent(detail::node(rel.lhs), lhs_name, rhs_name, exponents);
        collect_node_exponent(detail::node(rel.rhs), lhs_name, rhs_name, exponents);
    }
    return std::vector<int>(exponents.begin(), exponents.end());
}

Result<void> InferenceEngine::Impl::add_deduced(
    InferenceEngine& engine,
    const SymbolicExpr& new_lhs,
    const SymbolicExpr& new_rhs,
    RelationStore& store,
    PropertyStore& prop_store,
    int depth) const {
    if (store.has_relation(new_lhs, new_rhs, RelationalNode::Op::GT)) {
        return Result<void>::success();
    }

    auto inserted = store.add_relation_checked(
        new_lhs, new_rhs, RelationalNode::Op::GT, prop_store);
    if (!inserted.has_value()) {
        return Result<void>::failure(inserted.error());
    }

    Relation new_rel{new_lhs, new_rhs, RelationalNode::Op::GT};
    return engine.apply_monotonicity_rules_checked(new_rel, store, prop_store, depth + 1);
}

Result<void> InferenceEngine::Impl::deduce_function_relation(
    InferenceEngine& engine,
    const VariableNode& lhs,
    const VariableNode& rhs,
    FunctionNode::FuncType type,
    RelationStore& store,
    PropertyStore& prop_store,
    int depth) const {
    auto lhs_node = detail::make_node<FunctionNode>(
        type, std::vector<std::shared_ptr<const SymbolicNode>>{lhs.clone()});
    auto rhs_node = detail::make_node<FunctionNode>(
        type, std::vector<std::shared_ptr<const SymbolicNode>>{rhs.clone()});
    auto lhs_expr = detail::expression_from_node(lhs_node);
    auto rhs_expr = detail::expression_from_node(rhs_node);
    return add_deduced(engine, lhs_expr, rhs_expr, store, prop_store, depth);
}

Result<void> InferenceEngine::Impl::deduce_power_relations(
    InferenceEngine& engine,
    const VariableNode& lhs,
    const VariableNode& rhs,
    RelationStore& store,
    PropertyStore& prop_store,
    int depth) const {
    std::vector<int> exponents = collect_power_exponents(store, lhs.name(), rhs.name());

    for (int n : exponents) {
        auto exp_node = LMCAS::detail::make_node<NumberNode>(BigInt(n));
        auto pow_x_node = LMCAS::detail::make_node<PowerNode>(lhs.clone(), exp_node->clone());
        auto pow_y_node = LMCAS::detail::make_node<PowerNode>(rhs.clone(), exp_node->clone());

        auto pow_x = LMCAS::detail::expression_from_node(pow_x_node);
        auto pow_y = LMCAS::detail::expression_from_node(pow_y_node);
        auto deduced = add_deduced(engine, pow_x, pow_y, store, prop_store, depth);
        if (!deduced.has_value()) {
            return deduced;
        }
    }
    return Result<void>::success();
}

// Monotonicity inference

/**
 * @brief If the node is multiplication by exactly -1, returns the other
 * operand.
 */
static std::shared_ptr<const SymbolicNode> detect_negation(
    const std::shared_ptr<const SymbolicNode>& node) {
    const auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!multiply || multiply->operands().size() != 2) {
        return nullptr;
    }
    for (std::size_t index = 0; index < 2; ++index) {
        const auto number = std::dynamic_pointer_cast<const NumberNode>(
            multiply->operands()[index]);
        if (number && number->is_negative_one()) {
            return multiply->operands()[1 - index];
        }
    }
    return nullptr;
}

Monotonicity InferenceEngine::infer_monotonicity(const SymbolicExpr& expr,
                                                  const std::string& var,
                                                  const Interval& interval) const {
    if (!LMCAS::detail::node(expr)) {
        return Monotonicity::Unknown;
    }

    // FunctionNode: auto-infer for known functions
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(LMCAS::detail::node(expr))) {
        return Impl::function_monotonicity(*func, var, interval);
    }

    // Negation: multiply by -1 reverses monotonicity
    if (auto inner = detect_negation(LMCAS::detail::node(expr))) {
        auto inner_expr = LMCAS::detail::expression_from_node(inner);
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
    if (auto var_node = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(expr))) {
        const auto& props = impl_->ctx;
        auto monotonicity =
            props.get_monotonicity_checked(var_node->name(), var, interval);
        return monotonicity ? monotonicity.value() : Monotonicity::Unknown;
    }

    return Monotonicity::Unknown;
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
    if (depth >= MAX_MONOTONICITY_DEPTH) {
        return Result<void>::success();
    }

    // Only apply to GT (greater-than) relations
    if (rel.op != RelationalNode::Op::GT) {
        return Result<void>::success();
    }

    // Extract LHS and RHS — both must be single VariableNodes
    auto lhs_var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(rel.lhs));
    auto rhs_var = std::dynamic_pointer_cast<const VariableNode>(LMCAS::detail::node(rel.rhs));
    if (!lhs_var || !rhs_var) {
        return Result<void>::success();
    }

    const std::string& x_name = lhs_var->name();
    const std::string& y_name = rhs_var->name();


    // Check domain assumptions using the AssumptionContext (read-through all scopes)
    bool both_positive = impl_->ctx.has_sign(x_name, Sign::Positive) &&
                         impl_->ctx.has_sign(y_name, Sign::Positive);
    bool both_real = impl_->ctx.has_domain(x_name, Domain::Real) &&
                     impl_->ctx.has_domain(y_name, Domain::Real);
    bool both_nonnegative = impl_->ctx.has_sign(x_name, Sign::NonNegative) &&
                            impl_->ctx.has_sign(y_name, Sign::NonNegative);

    if (both_positive) {
        auto deduced = impl_->deduce_function_relation(
            *this, *lhs_var, *rhs_var, FunctionNode::FuncType::Ln,
            store, prop_store, depth);
        if (!deduced.has_value()) {
            return deduced;
        }
    }

    if (both_positive) {
        auto deduced = impl_->deduce_function_relation(
            *this, *lhs_var, *rhs_var, FunctionNode::FuncType::Sqrt,
            store, prop_store, depth);
        if (!deduced.has_value()) {
            return deduced;
        }
    }

    if (both_real) {
        auto deduced = impl_->deduce_function_relation(
            *this, *lhs_var, *rhs_var, FunctionNode::FuncType::Exp,
            store, prop_store, depth);
        if (!deduced.has_value()) {
            return deduced;
        }
    }

    if (both_nonnegative) {
        return impl_->deduce_power_relations(
            *this, *lhs_var, *rhs_var, store, prop_store, depth);
    }

    return Result<void>::success();
}

} // namespace LMCAS
