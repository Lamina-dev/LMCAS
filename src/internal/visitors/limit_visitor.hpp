/**
 * @file limit_visitor.hpp
 * @brief 极限访问器，通过代入和 L'Hôpital 法则计算极限。
 */
#pragma once

#include "lmcas_export.hpp"
#include "computation_context.hpp"
#include "internal/symbolic_ast.hpp"
#include "assumption_context.hpp"
#include "internal/visitors/normalization_visitor.hpp"
#include "internal/visitors/differentiation_visitor.hpp"
#include <iostream>
#include <cmath>
#include <optional>
#include <limits>
#include <utility>

namespace LMCAS {

class SymbolicExpr;
template <typename CoeffType>
class Polynomial;
namespace detail {
class ExactBoundArithmetic;
}
class LMCAS_API LimitVisitor : public LMCAS::detail::SymbolicVisitor {
    std::string var;
    std::shared_ptr<const SymbolicNode> point;
    std::string direction;
    const LMCAS::AssumptionContext* assumption_ctx_ = nullptr;
    std::optional<ComputationContext> local_context_;
    ComputationContext* context_ = nullptr;
    void charge() const {
        if (!context_) return;
        auto step = context_->consume_steps(1, "limit");
        if (!step) throw step.error();
    }
    int lhopital_depth_ = 0;
    static constexpr int max_lhopital_depth = 5;
    static int& active_limit_depth();
    static constexpr int max_active_limit_depth = 128;

    static std::shared_ptr<const SymbolicNode> make_product_or_one(
        const std::vector<std::shared_ptr<const SymbolicNode>>& factors);

    struct EvaluationScope {
        bool entered = false;
        EvaluationScope();
        ~EvaluationScope();
        explicit operator bool() const noexcept { return entered; }
    };

    enum class IndeterminateForm { None, ZeroTimesInf, InfMinusInf, OnePowInf, ZeroPowZero, InfPowZero };

    bool is_inf(const std::shared_ptr<const SymbolicNode>& node) const;
    bool is_neg_inf(const std::shared_ptr<const SymbolicNode>& node) const;
    std::optional<int> get_node_sign(const std::shared_ptr<const SymbolicNode>& node) const;
    double get_numeric_value(const std::shared_ptr<const NumberNode>& num) const;
    bool is_bounded(const std::shared_ptr<const SymbolicNode>& node) const;
    bool is_bounded_expression(const std::shared_ptr<const SymbolicNode>& node) const;
    bool tends_to_zero(const std::shared_ptr<const SymbolicNode>& node) const;
    std::shared_ptr<const SymbolicNode> eval_limit(const std::shared_ptr<const SymbolicNode>& expr);

    IndeterminateForm classify_product_form(const std::vector<std::shared_ptr<const SymbolicNode>>& vals);
    IndeterminateForm classify_power_form(const std::shared_ptr<const SymbolicNode>& bv, const std::shared_ptr<const SymbolicNode>& ev);
    IndeterminateForm classify_add_form(const std::vector<std::shared_ptr<const SymbolicNode>>& vals);

    std::shared_ptr<const SymbolicNode> try_squeeze(const std::shared_ptr<const SymbolicNode>& expr);
    std::pair<std::shared_ptr<const SymbolicNode>, std::shared_ptr<const SymbolicNode>> extract_num_den(const std::shared_ptr<const SymbolicNode>& expr);
    std::shared_ptr<const SymbolicNode> resolve_zero_times_inf(const std::vector<std::shared_ptr<const SymbolicNode>>& factors, const std::vector<std::shared_ptr<const SymbolicNode>>& factor_vals);
    std::shared_ptr<const SymbolicNode> resolve_inf_minus_inf(const AddNode& node, const std::vector<std::shared_ptr<const SymbolicNode>>& operand_vals);
    std::shared_ptr<const SymbolicNode> resolve_exponential_form(const std::shared_ptr<const SymbolicNode>& base, const std::shared_ptr<const SymbolicNode>& exponent);
    std::shared_ptr<const SymbolicNode> apply_lhopital(const std::shared_ptr<const SymbolicNode>& num, const std::shared_ptr<const SymbolicNode>& den);

    std::shared_ptr<const SymbolicNode> taylor_fallback(const std::shared_ptr<const SymbolicNode>& num, const std::shared_ptr<const SymbolicNode>& den, int max_order = 8);
    std::shared_ptr<const SymbolicNode> simplify_and_eval_ratio(const std::shared_ptr<const SymbolicNode>& ratio_node);
    std::pair<std::shared_ptr<const SymbolicNode>, int> find_leading_term(const std::shared_ptr<SymbolicExpr>& series_expr, const std::string& expand_var, const std::shared_ptr<SymbolicExpr>& expand_point, int max_order);
    static int get_sign(const std::shared_ptr<const SymbolicNode>& node);

    bool is_limit_at_infinity() const;
    bool is_limit_at_neg_infinity() const;
    int get_polynomial_degree(const std::shared_ptr<const SymbolicNode>& node) const;
    std::shared_ptr<const SymbolicNode> get_leading_coefficient(const std::shared_ptr<const SymbolicNode>& node) const;
    std::shared_ptr<const SymbolicNode> limit_rational_at_infinity(const std::shared_ptr<const SymbolicNode>& num, const std::shared_ptr<const SymbolicNode>& den);
    enum class GrowthClass { Constant, Logarithmic, Polynomial, Exponential, Unknown };
    GrowthClass classify_growth(const std::shared_ptr<const SymbolicNode>& node) const;
    std::shared_ptr<const SymbolicNode> limit_by_growth_comparison(const std::shared_ptr<const SymbolicNode>& num, const std::shared_ptr<const SymbolicNode>& den);
    std::shared_ptr<const SymbolicNode> handle_neg_infinity_limit(const std::shared_ptr<const SymbolicNode>& expr);
    std::shared_ptr<const SymbolicNode> positive_tail_log_power(const FunctionNode& node) const;
    bool defined_near_point(const std::shared_ptr<const SymbolicNode>& node) const;
    std::shared_ptr<const SymbolicNode> substitute_neg_t(const std::shared_ptr<const SymbolicNode>& node, const std::string& t_var) const;

    std::optional<int> determine_sign_near_point(const std::shared_ptr<const SymbolicNode>& expr, const std::string& dir);
    bool analytic_at_point(const std::shared_ptr<const SymbolicNode>& expr) const;
    std::optional<Rational> algebraic_growth_degree(const std::shared_ptr<const SymbolicNode>& expr) const;
    bool proves_superlogarithmic_magnitude(const std::shared_ptr<const SymbolicNode>& expr) const;
    bool polynomial_logarithmic_bound(const std::shared_ptr<const SymbolicNode>& expr) const;
    std::optional<std::pair<Rational, int>> exact_leading_at_point(
        const std::shared_ptr<const SymbolicNode>& expr) const;
    std::shared_ptr<const SymbolicNode> select_branch_by_direction(const PiecewiseNode& node, const std::string& dir);
    std::optional<std::shared_ptr<const SymbolicNode>> evaluate_sgn_limit(const std::shared_ptr<const SymbolicNode>& arg);
    std::optional<std::shared_ptr<const SymbolicNode>> evaluate_abs_limit(const std::shared_ptr<const SymbolicNode>& arg);

    bool is_infinite_power(const PowerNode& node) const;
    bool is_bounded_product(const MultiplyNode& node) const;
    bool is_zero_limit(const std::shared_ptr<const SymbolicNode>& node) const;
    std::shared_ptr<const SymbolicNode> squeeze_product(const MultiplyNode& node);
    std::shared_ptr<const SymbolicNode> squeeze_sum(const AddNode& node);
    std::shared_ptr<const SymbolicNode> denominator_factor(const std::shared_ptr<const SymbolicNode>& node) const;
    std::shared_ptr<const SymbolicNode> resolve_zero_infinity_pair(
        const std::shared_ptr<const SymbolicNode>& zero,
        const std::shared_ptr<const SymbolicNode>& infinity);
    std::shared_ptr<const SymbolicNode> multiply_remaining_limits(
        const std::shared_ptr<const SymbolicNode>& value,
        const std::vector<std::shared_ptr<const SymbolicNode>>& factors);
    std::shared_ptr<const SymbolicNode> resolve_combined_ratio(
        const std::shared_ptr<const SymbolicNode>& num,
        const std::shared_ptr<const SymbolicNode>& den);
    std::shared_ptr<const SymbolicNode> derivative_limit(
        const std::shared_ptr<const SymbolicNode>& node, bool normalize);
    std::shared_ptr<const SymbolicNode> derivative_ratio_limit(
        const std::shared_ptr<const SymbolicNode>& num,
        const std::shared_ptr<const SymbolicNode>& den,
        const std::shared_ptr<const SymbolicNode>& num_value,
        const std::shared_ptr<const SymbolicNode>& den_value);
    std::shared_ptr<const SymbolicNode> signed_infinity(int sign) const;
    bool is_limit_variable(const std::shared_ptr<const SymbolicNode>& node) const;
    bool is_positive_variable_power(const std::shared_ptr<const SymbolicNode>& node) const;
    bool is_variable_reciprocal(const std::shared_ptr<const SymbolicNode>& node) const;
    bool is_log_variable(const std::shared_ptr<const SymbolicNode>& node) const;
    bool logarithmic_zero_product(const MultiplyNode& node) const;
    std::shared_ptr<const SymbolicNode> one_plus_term(const AddNode& node) const;
    std::shared_ptr<const SymbolicNode> variable_fraction_numerator(
        const std::shared_ptr<const SymbolicNode>& node) const;
    bool logarithmic_product_limit(const MultiplyNode& node);
    std::optional<std::shared_ptr<const SymbolicNode>> logarithmic_remaining_limit(
        const MultiplyNode& node, const SymbolicNode* logarithm,
        const std::shared_ptr<const SymbolicNode>& numerator);
    bool has_slower_factors(const MultiplyNode& node, const SymbolicNode* exponential) const;
    bool decaying_exponential_product(const MultiplyNode& node);
    bool trigonometric_product_limit(const MultiplyNode& node);
    std::shared_ptr<const SymbolicNode> matching_reciprocal(
        const MultiplyNode& node, const std::shared_ptr<const SymbolicNode>& argument) const;
    std::shared_ptr<const SymbolicNode> trigonometric_remaining_limit(
        const MultiplyNode& node, const SymbolicNode* function, const SymbolicNode* reciprocal);
    std::shared_ptr<const SymbolicNode> product_ratio_limit(
        const std::shared_ptr<const SymbolicNode>& num,
        const std::shared_ptr<const SymbolicNode>& den);
    std::shared_ptr<const SymbolicNode> linear_exponent_coefficient(
        const std::shared_ptr<const SymbolicNode>& node) const;
    bool exponential_power_limit(const PowerNode& node);
    std::shared_ptr<const SymbolicNode> singular_power_limit(
        const PowerNode& node, const std::shared_ptr<const SymbolicNode>& base,
        const std::shared_ptr<const SymbolicNode>& exponent);
    bool composed_function_limit(const FunctionNode& node);
    bool infinite_function_limit(FunctionNode::FuncType type, bool negative);
    GrowthClass function_growth(const FunctionNode& node) const;
    GrowthClass power_growth(const PowerNode& node) const;
    GrowthClass operands_growth(const std::vector<std::shared_ptr<const SymbolicNode>>& operands) const;
    std::shared_ptr<const SymbolicNode> substitute_function_neg_t(const FunctionNode& node, const std::string& t_var) const;
    std::shared_ptr<const SymbolicNode> leading_terms_limit(
        const std::pair<std::shared_ptr<const SymbolicNode>, int>& numerator,
        const std::pair<std::shared_ptr<const SymbolicNode>, int>& denominator);
    std::shared_ptr<const SymbolicNode> taylor_series_limit(
        const std::shared_ptr<SymbolicExpr>& numerator,
        const std::shared_ptr<SymbolicExpr>& denominator,
        const std::string& expand_var, const std::shared_ptr<SymbolicExpr>& expand_point, int order);

    std::shared_ptr<const SymbolicNode> divergent_leading_ratio(
        const std::shared_ptr<SymbolicExpr>& ratio, int difference);
    bool is_infinite_product(const MultiplyNode& product) const;
    bool is_infinite_sum(const AddNode& sum) const;
    std::optional<int> assumed_node_sign(
        const std::shared_ptr<const SymbolicNode>& node) const;
    std::shared_ptr<const SymbolicNode> infinite_base_power_limit(
        const std::shared_ptr<const SymbolicNode>& base,
        const std::shared_ptr<const NumberNode>& number, bool integral, const BigInt& integer);
    bool rational_denominator_supported(
        const std::shared_ptr<const SymbolicNode>& leading_denominator,
        const std::shared_ptr<const SymbolicNode>& num,
        const std::shared_ptr<const SymbolicNode>& den) const;
    std::shared_ptr<const SymbolicNode> rational_leading_limit(
        const std::shared_ptr<const SymbolicNode>& leading_numerator,
        const std::shared_ptr<const SymbolicNode>& leading_denominator, int degree_difference);
    std::optional<Rational> algebraic_sum_growth(const AddNode& sum) const;
    std::optional<Rational> algebraic_product_growth(const MultiplyNode& product) const;
    std::optional<Rational> algebraic_power_growth(const PowerNode& power) const;
    bool logarithmic_algebraic_bound(
        const std::shared_ptr<const SymbolicNode>& expr) const;
    bool polynomial_logarithmic_operands(
        const std::vector<std::shared_ptr<const SymbolicNode>>& operands) const;
    bool exponential_dominates(
        const std::shared_ptr<const SymbolicNode>& exponential,
        const std::shared_ptr<const SymbolicNode>& other);
    bool can_rationalize_sum(const AddNode& node,
        const std::vector<std::shared_ptr<const SymbolicNode>>& values) const;
    std::optional<Polynomial<Rational>> square_algebraic_term(
        const std::shared_ptr<const SymbolicNode>& term) const;
    std::shared_ptr<const SymbolicNode> rationalized_sum_limit(
        const AddNode& node, const std::vector<std::shared_ptr<const SymbolicNode>>& values);
    std::pair<std::shared_ptr<const SymbolicNode>, std::shared_ptr<const SymbolicNode>>
    combine_sum_ratio(
        const std::vector<std::shared_ptr<const SymbolicNode>>& numerators,
        const std::vector<std::shared_ptr<const SymbolicNode>>& denominators);
    std::optional<std::pair<Rational, int>> exact_polynomial_leading(
        const Polynomial<Rational>& polynomial, const Rational& point_value,
        detail::ExactBoundArithmetic& arithmetic) const;
    bool analytic_defined_at_point(
        const std::shared_ptr<const SymbolicNode>& expr,
        const std::shared_ptr<const SymbolicNode>& positive_base) const;
    bool analytic_power_at_point(const PowerNode& power,
        const std::shared_ptr<const SymbolicNode>& expr) const;
    bool analytic_function_at_point(const FunctionNode& function,
        const std::shared_ptr<const SymbolicNode>& expr) const;
    bool analytic_operands_at_point(
        const std::vector<std::shared_ptr<const SymbolicNode>>& operands) const;
    std::optional<int> rational_sign_at_infinity(
        const std::shared_ptr<const SymbolicNode>& expr);
    std::optional<int> product_sign_near_point(
        const MultiplyNode& product, const std::string& dir);
    std::optional<int> power_sign_near_point(
        const PowerNode& power, const std::string& dir);
    std::optional<int> analytic_sign_near_point(
        const std::shared_ptr<const SymbolicNode>& expr, const std::string& dir);
    std::optional<bool> relation_truth_near_point(const RelationalNode& relation);
    std::optional<bool> logical_truth_near_point(const LogicalNode& logical);
    std::optional<bool> condition_truth_near_point(
        const std::shared_ptr<const SymbolicNode>& node);
    std::shared_ptr<const SymbolicNode> singular_ratio_limit(
        const std::shared_ptr<const SymbolicNode>& num,
        const std::shared_ptr<const SymbolicNode>& den);
    std::optional<std::shared_ptr<const SymbolicNode>> exact_ratio_at_point(
        const std::shared_ptr<const SymbolicNode>& num,
        const std::shared_ptr<const SymbolicNode>& den);
    std::shared_ptr<const SymbolicNode> evaluated_product_ratio(
        const std::shared_ptr<const SymbolicNode>& num,
        const std::shared_ptr<const SymbolicNode>& den, bool at_infinity);
    bool positive_tail_product_limit(const MultiplyNode& node);
    std::shared_ptr<const SymbolicNode> resolve_product_factors(
        const MultiplyNode& node, const std::vector<std::shared_ptr<const SymbolicNode>>& values);

public:
    bool condition_holds_near_point(const std::shared_ptr<const SymbolicNode>& condition);
    std::optional<int> sign_near_point(const std::shared_ptr<const SymbolicNode>& node) {
        return determine_sign_near_point(node, direction);
    }
    std::shared_ptr<const SymbolicNode> select_piecewise_near_point(const PiecewiseNode& node) {
        return select_branch_by_direction(node, direction);
    }
    std::shared_ptr<const SymbolicNode> result;

    LimitVisitor(std::string v, std::shared_ptr<const SymbolicNode> p, std::string dir = "",
                 const LMCAS::AssumptionContext* ctx = nullptr, ComputationContext* context = nullptr)
        : var(std::move(v)), point(std::move(p)), direction(std::move(dir)), assumption_ctx_(ctx),
          local_context_(context ? std::nullopt
              : std::optional<ComputationContext>(std::in_place)),
          context_(context ? context : &*local_context_) {
        if (context_) {
            auto entered = context_->enter_recursion("limit");
            if (!entered) throw entered.error();
        }
    }
    ~LimitVisitor() override { if (context_) context_->leave_recursion(); }
    LimitVisitor(const LimitVisitor&) = delete;
    LimitVisitor& operator=(const LimitVisitor&) = delete;

    std::shared_ptr<const SymbolicNode> get_result() const;

    void visit(const NumberNode& node) override;
    void visit(const VariableNode& node) override;
    void visit(const AddNode& node) override;
    void visit(const MultiplyNode& node) override;
    void visit(const PowerNode& node) override;
    void visit(const FunctionNode& node) override;
    void visit(const UninterpretedFunctionNode& node) override;
    void visit(const MatrixNode& node) override;
    void visit(const RelationalNode& node) override;
    void visit(const LogicalNode& node) override;
    void visit(const PiecewiseNode& node) override;
    void visit(const SummationNode& node) override;
    void visit(const ProductNode& node) override;
    void visit(const TransformNode& node) override;
    void visit(const QuantifierNode& node) override;
    void visit(const SetBuilderNode& node) override;
    void visit(const FiniteSetNode& node) override;
    void visit(const IntervalNode& node) override;
    void visit(const MembershipNode& node) override;
    void visit(const QuantityNode& node) override;
    void visit(const ComplexNode& node) override;
    void visit(const IntegralNode& node) override;
    void visit(const LimitNode& node) override;
    void visit(const RootOfNode& node) override;
};

}
