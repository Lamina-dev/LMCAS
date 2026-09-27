#pragma once

#include "assumption_context.hpp"
#include "computation_context.hpp"
#include "inference_engine.hpp"
#include "interval.hpp"
#include "numeric_evaluation.hpp"
#include "property_store.hpp"
#include "relation_store.hpp"
#include "internal/symbolic_ast.hpp"

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cmath>
#include <functional>
#include <optional>
#include <unordered_set>

namespace LMCAS {

bool is_integer_number(const NumberNode& number);
bool is_even_integer_number(const NumberNode& number);
bool is_positive_integer_number(const NumberNode& number);
bool is_zero_number(const NumberNode& number);

inline bool is_negative_infinity_coefficient(const NumberNode& number) {
    if (std::holds_alternative<BigInt>(number.value())) {
        return std::get<BigInt>(number.value()) < BigInt(0);
    }
    if (std::holds_alternative<Rational>(number.value())) {
        return std::get<Rational>(number.value()) < Rational(0);
    }
    if (std::holds_alternative<lmmc_real_t>(number.value())) {
        return std::get<lmmc_real_t>(number.value()) < 0;
    }
    return false;
}

inline int infinity_sign(const SymbolicExpr& expression, ComputationContext& context) {
    const auto& node = detail::node(expression);
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return function->type() == FunctionNode::FuncType::Infinity ? 1 : 0;
    }
    auto product = std::dynamic_pointer_cast<const MultiplyNode>(node);
    if (!product) {
        return 0;
    }

    bool has_infinity = false;
    int sign = 1;
    for (const auto& operand : product->operands()) {
        auto step = context.consume_steps(1, "inference.infinity");
        if (!step) throw step.error();
        if (auto function = std::dynamic_pointer_cast<const FunctionNode>(operand)) {
            has_infinity = has_infinity ||
                function->type() == FunctionNode::FuncType::Infinity;
        }
        if (auto number = std::dynamic_pointer_cast<const NumberNode>(operand)) {
            if (is_negative_infinity_coefficient(*number)) {
                sign = -sign;
            }
        }
    }
    return has_infinity ? sign : 0;
}

template <typename T>
Result<T> checked_inference_result(
    const SymbolicExpr& expr,
    const std::string& operation,
    ComputationContext& context,
    const std::function<Result<T>()>& query) {
    if (!detail::node(expr)) {
        return Result<T>::failure(
            CasErrc::InvalidArgument,
            "inference expression must not be null", operation);
    }
    try {
        auto entered = context.enter_recursion(operation);
        if (!entered) return Result<T>::failure(entered.error());
        struct Leave {
            ComputationContext& context;
            ~Leave() { context.leave_recursion(); }
        } leave{context};
        return query();
    } catch (const CasError& error) {
        return Result<T>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<T>::failure(
            CasErrc::ResourceLimit,
            "inference query allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<T>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}


struct InferenceEngine::Impl {
    explicit Impl(const AssumptionContext& context)
        : ctx(context), max_depth(context.get_max_query_depth()) {}
    InferenceTriboolResult query_zero(
        const InferenceEngine& engine, const SymbolicExpr& expr, ComputationContext& context) const;
    static InferenceTriboolResult real_definedness(
        const InferenceEngine& engine, const SymbolicNode& node, ComputationContext& context);
    static Tribool number_positive(const NumberNode& num);
    static Tribool number_negative(const NumberNode& num);
    static Tribool number_nonnegative(const NumberNode& num);
    static Tribool number_nonpositive(const NumberNode& num);
    static Tribool number_real(const NumberNode& num);
    static Tribool number_integer(const NumberNode& num);
    static Tribool number_nonzero(const NumberNode& num);
    static Tribool natural_number(const NumberNode& num);
    InferenceTriboolResult composite_positive(
        const InferenceEngine& engine, const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult composite_negative(
        const InferenceEngine& engine, const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult composite_nonnegative(
        const InferenceEngine& engine, const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult composite_nonpositive(
        const InferenceEngine& engine, const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult composite_real(
        const InferenceEngine& engine, const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult composite_integer(
        const InferenceEngine& engine, const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult composite_nonzero(
        const InferenceEngine& engine, const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult sum_nonzero(
        const InferenceEngine& engine, const AddNode& add, ComputationContext& context) const;
    InferenceTriboolResult query_rational(
        const InferenceEngine& engine, const SymbolicExpr& expr, ComputationContext& context) const;
    InferenceTriboolResult query_natural(const SymbolicExpr& expr, ComputationContext& context) const;
    static Tribool known_sign(int value, Sign target);
    static Tribool nonnegative_sign(Sign target);
    static std::shared_ptr<const SymbolicNode> subtrahend(
        const std::shared_ptr<const SymbolicNode>& operand, ComputationContext& context);
    static bool has_negated_operand(const AddNode& node, ComputationContext& context);
    InferenceTriboolResult sum_sign(
        const InferenceEngine& engine, const AddNode& node, Sign target, ComputationContext& context) const;
    InferenceTriboolResult subtraction_sign(
        const InferenceEngine& engine, const SymbolicExpr& minuend_expr,
        const SymbolicExpr& subtrahend_expr, Sign target, ComputationContext& context) const;
    InferenceTriboolResult weak_subtraction_sign(
        const InferenceEngine& engine, const SymbolicExpr& minuend_expr,
        const SymbolicExpr& subtrahend_expr, Sign target, ComputationContext& context) const;
    bool operand_sign(
        const std::shared_ptr<const SymbolicNode>& operand, Sign target, ComputationContext& context) const;
    Result<void> accumulate_sum_relation_sign(
        const std::shared_ptr<const SymbolicNode>& operand, const SymbolicExpr& zero_expr,
        bool& all_gt_zero, bool& all_geq_zero, ComputationContext& context) const;
    InferenceTriboolResult sum_relation_sign(
        const AddNode& add, const SymbolicExpr& zero_expr, Sign target, ComputationContext& context) const;
    InferenceTriboolResult product_relation_sign(
        const MultiplyNode& mul, const SymbolicExpr& zero_expr, Sign target, ComputationContext& context) const;
    Result<bool> relation_rhs_nonnegative(
        const Relation& rel, const SymbolicExpr& zero_expr, ComputationContext& context) const;
    InferenceTriboolResult positive_relation_chain(
        const SymbolicExpr& expr, const SymbolicExpr& zero_expr, ComputationContext& context) const;
    InferenceTriboolResult accumulate_product_sign(
        const InferenceEngine& engine,
        const std::shared_ptr<const SymbolicNode>& operand,
        std::uint8_t& possible, ComputationContext& context) const;
    static InferenceTriboolResult product_parity_sign(
        std::uint8_t possible, Sign target);
    InferenceTriboolResult product_sign(
        const InferenceEngine& engine, const MultiplyNode& node, Sign target, ComputationContext& context) const;
    InferenceTriboolResult quotient_sign(
        const InferenceEngine& engine, const SymbolicExpr& num_expr,
        const SymbolicExpr& den_expr, Sign target, ComputationContext& context) const;
    InferenceTriboolResult multiply_domain(
        const InferenceEngine& engine, const MultiplyNode& node, Domain target, ComputationContext& context) const;
    static InferenceTriboolResult integral_power_domain(
        const InferenceEngine& engine, const SymbolicExpr& base,
        const NumberNode& exponent, Domain target, ComputationContext& context);
    InferenceTriboolResult zero_power_sign(
        const InferenceEngine& engine, const SymbolicExpr& base_expr,
        const NumberNode* exp_num, Sign target, ComputationContext& context) const;
    InferenceTriboolResult positive_power_sign(
        const InferenceEngine& engine, const SymbolicExpr& base_expr,
        const SymbolicExpr& exp_expr, Sign target, ComputationContext& context) const;
    InferenceTriboolResult even_power_sign(
        const InferenceEngine& engine, const SymbolicExpr& base_expr,
        const NumberNode& exponent, Sign target, ComputationContext& context) const;
    InferenceTriboolResult positive_integer_power_sign(
        const InferenceEngine& engine, const SymbolicExpr& base_expr,
        const NumberNode& exponent, Sign target, ComputationContext& context) const;
    InferenceTriboolResult nonzero_power_sign(
        const InferenceEngine& engine, const SymbolicExpr& base_expr,
        const NumberNode& exponent, Sign target, ComputationContext& context) const;
    InferenceTriboolResult integer_power_sign(
        const InferenceEngine& engine, const SymbolicExpr& base_expr,
        const NumberNode& exponent, Sign target, ComputationContext& context) const;
    InferenceTriboolResult exponential_sign(
        const InferenceEngine& engine, const SymbolicExpr& arg_expr, Sign target, ComputationContext& context) const;
    InferenceTriboolResult trigonometric_sign(
        const InferenceEngine& engine, const SymbolicExpr& arg_expr, ComputationContext& context) const;
    InferenceTriboolResult absolute_sign(
        const InferenceEngine& engine, const SymbolicExpr& arg_expr, Sign target, ComputationContext& context) const;
    InferenceTriboolResult logarithm_sign(
        const InferenceEngine& engine, const SymbolicExpr& arg_expr, Sign target, ComputationContext& context) const;
    InferenceTriboolResult square_root_sign(
        const InferenceEngine& engine, const SymbolicExpr& arg_expr, Sign target, ComputationContext& context) const;
    InferenceTriboolResult real_or_integer(
        const InferenceEngine& engine, const SymbolicExpr& arg_expr, ComputationContext& context) const;
    InferenceTriboolResult exponential_domain(
        const InferenceEngine& engine, const SymbolicExpr& arg_expr, Domain target, ComputationContext& context) const;
    InferenceTriboolResult absolute_domain(
        const InferenceEngine& engine, const SymbolicExpr& arg_expr, Domain target, ComputationContext& context) const;
    InferenceTriboolResult logarithm_domain(
        const InferenceEngine& engine, const SymbolicExpr& arg_expr, Domain target, ComputationContext& context) const;
    InferenceTriboolResult square_root_domain(
        const InferenceEngine& engine, const SymbolicExpr& arg_expr, Domain target, ComputationContext& context) const;
    InferenceTriboolResult arctangent_domain(
        const InferenceEngine& engine, const SymbolicExpr& arg_expr, Domain target, ComputationContext& context) const;
    InferenceTriboolResult function_sign(
        const InferenceEngine& engine, FunctionNode::FuncType type,
        const SymbolicExpr& argument, Sign target, ComputationContext& context) const;
    InferenceTriboolResult function_domain(
        const InferenceEngine& engine, FunctionNode::FuncType type,
        const SymbolicExpr& argument, Domain target, ComputationContext& context) const;
    static Monotonicity logarithm_monotonicity(const Interval& interval);
    static Monotonicity function_monotonicity(
        const FunctionNode& func, const std::string& var, const Interval& interval);
    static std::optional<int> small_positive_exponent(
        const std::shared_ptr<const NumberNode>& number);
    static void collect_node_exponent(
        const std::shared_ptr<const SymbolicNode>& node,
        const std::string& lhs_name, const std::string& rhs_name,
        std::unordered_set<int>& exponents);
    static std::vector<int> collect_power_exponents(
        const RelationStore& store,
        const std::string& lhs_name, const std::string& rhs_name);
    Result<void> add_deduced(
        InferenceEngine& engine, const SymbolicExpr& new_lhs,
        const SymbolicExpr& new_rhs, RelationStore& store,
        PropertyStore& prop_store, int depth) const;
    Result<void> deduce_function_relation(
        InferenceEngine& engine, const VariableNode& lhs,
        const VariableNode& rhs, FunctionNode::FuncType type,
        RelationStore& store, PropertyStore& prop_store, int depth) const;
    Result<void> deduce_power_relations(
        InferenceEngine& engine, const VariableNode& lhs,
        const VariableNode& rhs, RelationStore& store,
        PropertyStore& prop_store, int depth) const;

    const AssumptionContext& ctx;
    int max_depth;
    mutable std::unordered_set<const SymbolicNode*> visited;
    mutable int current_depth = 0;
};

class InferenceEngine::DepthGuard {
public:
    DepthGuard(const InferenceEngine& engine, const SymbolicNode& node);
    ~DepthGuard();

    bool should_abort() const { return abort_; }

    DepthGuard(const DepthGuard&) = delete;
    DepthGuard& operator=(const DepthGuard&) = delete;

private:
    const InferenceEngine& engine_;
    const SymbolicNode* node_;
    bool abort_ = false;
    bool inserted_ = false;
};

} // namespace LMCAS
