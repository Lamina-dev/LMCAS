/**
 * @file limit_visitor.hpp
 * @brief 极限访问器，通过代入和 L'Hôpital 法则计算极限。
 */
#pragma once

#include "../lamina_export.hpp"
#include "../symbolic_ast.hpp"
#include "../assumption_context.hpp"
#include "normalization_visitor.hpp"
#include "differentiation_visitor.hpp"
#include <iostream>
#include <cmath>
#include <optional>
#include <limits>
#include <utility>

class SymbolicExpr;

/** @brief 极限访问器，通过代入求值和 L'Hôpital 法则计算符号表达式的极限 */
class LAMINA_API LimitVisitor : public lamina::detail::SymbolicVisitor {
    std::string var;
    std::shared_ptr<const SymbolicNode> point;
    std::string direction;
    const lamina::AssumptionContext* assumption_ctx_ = nullptr;
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
    double get_point_value() const;
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
    int get_growth_polynomial_degree(const std::shared_ptr<const SymbolicNode>& node) const;
    std::shared_ptr<const SymbolicNode> limit_by_growth_comparison(const std::shared_ptr<const SymbolicNode>& num, const std::shared_ptr<const SymbolicNode>& den);
    std::shared_ptr<const SymbolicNode> handle_neg_infinity_limit(const std::shared_ptr<const SymbolicNode>& expr);
    std::shared_ptr<const SymbolicNode> substitute_neg_t(const std::shared_ptr<const SymbolicNode>& node, const std::string& t_var) const;

    int determine_sign_near_point(const std::shared_ptr<const SymbolicNode>& expr, const std::string& dir);
    std::shared_ptr<const SymbolicNode> select_branch_by_direction(const PiecewiseNode& node, const std::string& dir);
    bool condition_satisfied_by_direction(const std::shared_ptr<const SymbolicNode>& condition, const std::string& dir);
    std::optional<int> evaluate_relational_sign(const std::shared_ptr<const RelationalNode>& rel, const std::string& dir);
    std::optional<std::shared_ptr<const SymbolicNode>> evaluate_sgn_limit(const std::shared_ptr<const SymbolicNode>& arg);
    std::optional<std::shared_ptr<const SymbolicNode>> evaluate_abs_limit(const std::shared_ptr<const SymbolicNode>& arg);

public:
    std::shared_ptr<const SymbolicNode> result;

    LimitVisitor(std::string v, std::shared_ptr<const SymbolicNode> p, std::string dir = "",
                 const lamina::AssumptionContext* ctx = nullptr)
        : var(std::move(v)), point(std::move(p)), direction(std::move(dir)), assumption_ctx_(ctx) {}

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
