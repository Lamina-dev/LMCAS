#include "internal/equivalence_support.hpp"

#include "quantity.hpp"
#include "internal/visitors/expand_visitor.hpp"
#include "internal/visitors/normalization_visitor.hpp"

#include <exception>
#include <limits>
#include <memory>
#include <new>
#include <utility>

namespace LMCAS {
namespace equivalence_detail {

ExprPtr normalize_equivalence(const ExprPtr& expression, detail::RewriteBudget& budget) {
    if (!expression || !detail::node(expression)) {
        throw CasError{CasErrc::InternalInvariant,
                       "equivalence normalization returned null", kEquivalentOperation};
    }
    budget.measure(detail::node(expression));
    NormalizationVisitor visitor(&budget);
    detail::node(expression)->accept(visitor);
    auto result = visitor.get_result();
    if (!result) {
        throw CasError{CasErrc::InternalInvariant,
                       "equivalence normalization returned null", kEquivalentOperation};
    }
    budget.measure(result);
    return detail::make_expression_ptr(std::move(result));
}

ExprPtr expand_equivalence(const ExprPtr& expression, detail::RewriteBudget& budget) {
    if (!expression || !detail::node(expression)) {
        throw CasError{CasErrc::InternalInvariant,
                       "equivalence expansion returned null", kEquivalentOperation};
    }
    budget.measure(detail::node(expression));
    ExpandVisitor visitor(&budget);
    detail::node(expression)->accept(visitor);
    auto result = visitor.get_result();
    if (!result) {
        throw CasError{CasErrc::InternalInvariant,
                       "equivalence expansion returned null", kEquivalentOperation};
    }
    return normalize_equivalence(detail::make_expression_ptr(std::move(result)), budget);
}

}
namespace {

Result<bool> equivalent_with_budget(const SymbolicExpr& lhs, const SymbolicExpr& rhs,
                                    ComputationContext& context, const EqvOptions& options,
                                    detail::RewriteBudget* budget);

Result<std::pair<ExprPtr, ExprPtr>> prepare_operands(
    const SymbolicExpr& lhs, const SymbolicExpr& rhs,
    ComputationContext& context, detail::RewriteBudget& budget) {
    using OperandsResult = Result<std::pair<ExprPtr, ExprPtr>>;
    auto left = std::make_shared<SymbolicExpr>(lhs);
    auto right = std::make_shared<SymbolicExpr>(rhs);
    if (auto quantity = std::dynamic_pointer_cast<const QuantityNode>(detail::node(lhs))) {
        normalization_check_arithmetic<MultiplyNode>(&budget, 1, quantity->value());
        auto stripped = strip_unit(left, UnitStripMode::BaseValue, context);
        if (!stripped) { return OperandsResult::failure(stripped.error()); }
        left = stripped.value();
    }
    if (auto quantity = std::dynamic_pointer_cast<const QuantityNode>(detail::node(rhs))) {
        normalization_check_arithmetic<MultiplyNode>(&budget, 1, quantity->value());
        auto stripped = strip_unit(right, UnitStripMode::BaseValue, context);
        if (!stripped) { return OperandsResult::failure(stripped.error()); }
        right = stripped.value();
    }
    return OperandsResult::success({std::move(left), std::move(right)});
}
Result<bool> apply_profile(const SymbolicExpr& lhs, const SymbolicExpr& rhs,
                           ComputationContext& context,
                           const EqvOptions& options, detail::RewriteBudget& budget) {
    if (options.profile != EqvProfile::TrigBasic &&
        options.profile != EqvProfile::ExpLogBasic) {
        return Result<bool>::success(false);
    }
    const bool trig = options.profile == EqvProfile::TrigBasic;
    budget.consume();
    ExprPtr rewritten_lhs;
    ExprPtr rewritten_rhs;
    if (trig) {
        rewritten_lhs = equivalence_detail::rewrite_trig_basic_identity(detail::node(lhs), budget);
        rewritten_rhs = equivalence_detail::rewrite_trig_basic_identity(detail::node(rhs), budget);
    } else {
        rewritten_lhs = equivalence_detail::rewrite_exp_log_basic_identity(
            detail::node(lhs), context.assumptions().get(), budget);
        rewritten_rhs = equivalence_detail::rewrite_exp_log_basic_identity(
            detail::node(rhs), context.assumptions().get(), budget);
    }
    rewritten_lhs = equivalence_detail::normalize_equivalence(rewritten_lhs, budget);
    rewritten_rhs = equivalence_detail::normalize_equivalence(rewritten_rhs, budget);
    EqvOptions core_options = options;
    core_options.profile = EqvProfile::Core;
    return equivalent_with_budget(*rewritten_lhs, *rewritten_rhs, context, core_options, &budget);
}

Result<bool> prove_operands(const SymbolicExpr& lhs, const SymbolicExpr& rhs,
                            ComputationContext& context,
                            const EqvOptions& options, detail::RewriteBudget& budget) {
    budget.consume();
    auto operands = prepare_operands(lhs, rhs, context, budget);
    if (!operands) { return Result<bool>::failure(operands.error()); }
    auto left = equivalence_detail::canonicalize_complex_product(
        *operands.value().first, budget);
    auto right = equivalence_detail::canonicalize_complex_product(
        *operands.value().second, budget);
    if (structurally_equal(*left, *right)) {
        return Result<bool>::success(true);
    }
    normalization_check_arithmetic<MultiplyNode>(&budget, 1, detail::node(right));
    auto negative = SymbolicExpr::multiply(SymbolicExpr::number(-1), right);
    normalization_check_arithmetic<AddNode>(&budget, 0, detail::node(left), detail::node(negative));
    auto difference = SymbolicExpr::add(left, negative);
    if (!difference) {
        return Result<bool>::failure(CasErrc::InternalInvariant,
                                     "equivalence difference construction failed",
                                     kEquivalentOperation);
    }
    difference = equivalence_detail::normalize_equivalence(difference, budget);
    if (difference->is_zero()) {
        return Result<bool>::success(true);
    }
    auto polynomial = equivalence_detail::prove_rational_polynomial_equivalence(
        difference, context, options, budget);
    if (!polynomial) { return Result<bool>::failure(polynomial.error()); }
    if (polynomial.value()) { return Result<bool>::success(*polynomial.value()); }
    return Result<bool>::success(false);
}

Result<bool> equivalent_with_budget(const SymbolicExpr& lhs, const SymbolicExpr& rhs,
                                    ComputationContext& context, const EqvOptions& options,
                                    detail::RewriteBudget* budget) {
    try {
        auto valid = equivalence_detail::validate_eqv_options(options);
        if (!valid) { return Result<bool>::failure(valid.error()); }
        auto step = context.consume_steps(1, kEquivalentOperation);
        if (!step) { return Result<bool>::failure(step.error()); }
        auto lhs_dimension = dimension_of(lhs);
        if (!lhs_dimension) { return Result<bool>::failure(lhs_dimension.error()); }
        auto rhs_dimension = dimension_of(rhs);
        if (!rhs_dimension) { return Result<bool>::failure(rhs_dimension.error()); }
        if (lhs_dimension.value() != rhs_dimension.value()) {
            return Result<bool>::success(false);
        }
        if (budget) {
            return options.profile == EqvProfile::Core
                ? prove_operands(lhs, rhs, context, options, *budget)
                : apply_profile(lhs, rhs, context, options, *budget);
        }
        const auto input_nodes = detail::RewriteBudget::input_nodes(
            detail::node(lhs), detail::node(rhs), context,
            options.budget.max_rewrite_depth, kEquivalentOperation);
        const auto maximum = std::numeric_limits<std::size_t>::max();
        const auto factor = options.budget.max_node_growth_factor;
        const auto max_nodes = input_nodes > maximum / factor
            ? maximum : input_nodes * factor;
        detail::RewriteBudget initial_budget(
            context, options.budget.max_rewrite_steps,
            options.budget.max_rewrite_depth, max_nodes, kEquivalentOperation);
        return options.profile == EqvProfile::Core
            ? prove_operands(lhs, rhs, context, options, initial_budget)
            : apply_profile(lhs, rhs, context, options, initial_budget);
    } catch (const CasError& error) {
        return Result<bool>::failure(error);
    } catch (const std::bad_alloc&) {
        return Result<bool>::failure(CasErrc::ResourceLimit,
                                     "equivalence check allocation failed",
                                     kEquivalentOperation);
    } catch (const std::exception& error) {
        return Result<bool>::failure(CasErrc::Inconclusive, error.what(),
                                     kEquivalentOperation);
    }
}

}

bool structurally_equal(const SymbolicExpr& lhs, const SymbolicExpr& rhs) {
    const auto& left = detail::node(lhs);
    const auto& right = detail::node(rhs);
    if (!left || !right) return left == right;
    return left->equals(*right);
}

Result<bool> equivalent_core(const SymbolicExpr& lhs,
                             const SymbolicExpr& rhs,
                             ComputationContext& context,
                             const EqvOptions& options) {
    return equivalent_with_budget(lhs, rhs, context, options, nullptr);
}

Result<bool> equivalent_core(const SymbolicExpr& lhs,
                             const SymbolicExpr& rhs,
                             ComputationContext& context) {
    return equivalent_core(lhs, rhs, context, EqvOptions{});
}

}
