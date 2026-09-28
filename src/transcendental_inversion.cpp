#include "solve_polynomial.hpp"
#include "assumption_context.hpp"
#include "internal/assumption_simplification.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/facts_query.hpp"
#include "internal/numeric_probe.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/normalization_utils.hpp"
#include "lmmc/config.h"
#include <cmath>
#include <vector>
#include "internal/transcendental_solver_support.hpp"

namespace LMCAS::detail::transcendental {


bool try_evaluate_numeric(
    const std::shared_ptr<SymbolicExpr>& expr,
    lmmc_real_t& out) {
    auto numeric = detail::try_finite_numeric(expr);
    if (!numeric) {
        return false;
    }
    out = static_cast<lmmc_real_t>(*numeric);
    return true;
}

std::vector<std::shared_ptr<SymbolicExpr>> invert_lambert_w(
    const std::shared_ptr<SymbolicExpr>& c) {

    auto sol = SymbolicExpr::lambertw(c);
    if (!sol) {
        return {};
    }
    return {sol->simplify()};
}

std::vector<std::shared_ptr<SymbolicExpr>> invert_exp_base(
    const std::shared_ptr<SymbolicExpr>& base,
    const std::shared_ptr<SymbolicExpr>& c) {

    lmmc_real_t c_val;
    if (try_evaluate_numeric(c, c_val)) {
        if (c_val <= 0.0) {
            return {};
        }
    }

    auto ln_c = SymbolicExpr::ln(c);
    auto ln_a = SymbolicExpr::ln(base);
    auto sol = SymbolicExpr::divide(ln_c, ln_a)->simplify();
    return {sol};
}


SolutionSet conditional_preimage(const ExprPtr& original, const std::string& var) {
    auto relation = std::dynamic_pointer_cast<const RelationalNode>(detail::node(original));
    auto predicate = relation ? original : SymbolicExpr::eq(original, SymbolicExpr::number(0));
    return SolutionSet{ConditionalSolutions{ConditionSet{var, std::move(predicate), {}}}};
}

namespace {

ExprPtr subtract(const ExprPtr& left, const ExprPtr& right) {
    return SymbolicExpr::add(left, SymbolicExpr::multiply(SymbolicExpr::number(-1), right));
}

ExprPtr inverse_trig(FunctionNode::FuncType type, const ExprPtr& target) {
    return make_expression_ptr(make_node<FunctionNode>(
        type, std::vector<SymbolicNodePtr>{detail::node(target)}));
}

Result<std::optional<detail::AffineForm>> prepare_trigonometric_affine(
    const ExprPtr& inner, const std::string& var, ComputationContext& context) {
    auto affine = detail::recognize_affine(*inner, var, context);
    if (!affine || !affine.value()) { return affine; }
    /**
     * @brief 查询值前先规范化尚未求值的多项式系数。
     * a 的实数性未知时，偏移量 0*a 仍为零；已声明为零的斜率须在除法前识别。
     */
    for (auto* coefficient : {&affine.value()->slope, &affine.value()->offset}) {
        auto normalized = detail::simplify_expression(*coefficient, context);
        if (!normalized) { return Result<std::optional<detail::AffineForm>>::failure(normalized.error()); }
        *coefficient = std::move(normalized.value());
    }
    return affine;
}

Result<bool> real_affine_coefficients(const detail::AffineForm& affine, Tribool nonzero,
    const FactsQuery& facts, ComputationContext& context) {
    for (const auto& coefficient : {affine.slope, affine.offset}) {
        if (coefficient == affine.slope && nonzero == Tribool::False) {
            continue; /**< 已证有定义的零系数属于实数。 */
        }
        auto real_coefficient = detail::query_real_value(detail::node(coefficient), facts, context);
        if (!real_coefficient) { return Result<bool>::failure(real_coefficient.error()); }
        auto defined = detail::query_definedness(detail::node(coefficient), facts, Domain::Real, context);
        if (!defined) { return Result<bool>::failure(defined.error()); }
        if (real_coefficient.value() != Tribool::True || defined.value() != Tribool::True) {
            return false;
        }
    }
    return true;
}

std::optional<SolveResult> trigonometric_target_range(FunctionNode::FuncType type,
    const ExprPtr& rhs, const ExprPtr& original, const std::string& var,
    const FactsQuery& facts, ComputationContext& context, std::vector<ExprPtr>& conditions) {
    auto real = detail::query_real_value(detail::node(rhs), facts, context);
    if (!real) { return SolveResult::failure(real.error()); }
    if (real.value() == Tribool::False) { return SolveResult{SolutionSet{EmptySolutions{}}}; }
    if (type != FunctionNode::FuncType::Tan) {
        /**
         * @brief 构造分支条件前先消除精确端点等式。
         * 统一 +/-1 分支的恒真条件，使其可合并为完整的半周期解族。
         */
        auto lower_difference = detail::simplify_expression(SymbolicExpr::add(rhs, SymbolicExpr::number(1)), context);
        if (!lower_difference) { return SolveResult::failure(lower_difference.error()); }
        auto lower = detail::query_nonnegative_value(detail::node(lower_difference.value()), facts, context);
        if (!lower) { return SolveResult::failure(lower.error()); }
        auto upper_difference = detail::simplify_expression(subtract(SymbolicExpr::number(1), rhs), context);
        if (!upper_difference) { return SolveResult::failure(upper_difference.error()); }
        auto upper = detail::query_nonnegative_value(detail::node(upper_difference.value()), facts, context);
        if (!upper) { return SolveResult::failure(upper.error()); }
        if (lower.value() == Tribool::False || upper.value() == Tribool::False) {
            return SolveResult{SolutionSet{EmptySolutions{}}};
        }
        if (lower.value() != Tribool::True) {
            conditions.push_back(make_expression_ptr(make_node<RelationalNode>(
                detail::node(rhs), detail::node(SymbolicExpr::number(-1)), RelationalNode::Op::GEQ)));
        }
        if (upper.value() != Tribool::True) {
            conditions.push_back(make_expression_ptr(make_node<RelationalNode>(
                detail::node(rhs), detail::node(SymbolicExpr::number(1)), RelationalNode::Op::LEQ)));
        }
    } else if (real.value() != Tribool::True) {
        return SolveResult{conditional_preimage(original, var)};
    }
    return std::nullopt;
}

SolveResult constant_trigonometric_preimage(FunctionNode::FuncType type,
    const ExprPtr& offset, const ExprPtr& rhs, const ExprPtr& original, const std::string& var,
    const FactsQuery& facts, ComputationContext& context) {
        auto constant = make_expression_ptr(make_node<FunctionNode>(
            type, std::vector<SymbolicNodePtr>{detail::node(offset)}));
        auto residual = detail::simplify_expression(subtract(constant, rhs), context);
        if (!residual) { return SolveResult::failure(residual.error()); }
        auto unequal = detail::query_nonzero_value(detail::node(residual.value()), facts, Domain::Real, context);
        if (!unequal) { return SolveResult::failure(unequal.error()); }
        if (unequal.value() == Tribool::False) { return SolutionSet{UniversalSolutions{}}; }
        if (unequal.value() == Tribool::True) { return SolutionSet{EmptySolutions{}}; }
        return conditional_preimage(original, var);
}

Result<std::string> capture_free_integer_parameter(const ExprPtr& original,
    const std::string& var, ComputationContext& context) {
    auto occupied = detail::all_variable_names(detail::node(original));
    occupied.insert(var);
    if (context.assumptions()) {
        /**
         * @brief 生成的绑定变量须避开全部假设中的名称，包括与方程无关的事实。
         * 遍历私有作用域栈副本，保持调用方假设栈不变。
         */
        AssumptionContext scopes = *context.assumptions();
        const auto include_expression = [&](const SymbolicExpr& expression) {
            auto names = detail::all_variable_names(detail::node(expression));
            occupied.insert(names.begin(), names.end());
        };
        for (const auto& conditional : scopes.get_active_conditionals()) {
            include_expression(conditional.condition);
            include_expression(conditional.conclusion);
        }
        for (;;) {
            auto names = scopes.current_properties().get_all_symbols();
            occupied.insert(names.begin(), names.end());
            for (const auto& relation : scopes.current_relations().get_relations()) {
                include_expression(relation.lhs);
                include_expression(relation.rhs);
            }
            if (scopes.depth() == 1) { break; }
            auto popped = scopes.pop();
            if (!popped) { return Result<std::string>::failure(popped.error()); }
        }
    }
    std::string parameter = "_k";
    for (std::size_t suffix = 1; occupied.count(parameter); ++suffix) {
        parameter = "_k" + std::to_string(suffix);
    }
    return parameter;
}

ExprPtr trigonometric_principal_value(FunctionNode::FuncType type,
    const ExprPtr& rhs, const ExprPtr& pi) {
    using Type = FunctionNode::FuncType;
    auto principal = inverse_trig(type == Type::Sin ? Type::ArcSin : Type::ArcCos, rhs);
    const auto positive_half = SymbolicExpr::number(Rational(1, 2));
    const auto negative_half = SymbolicExpr::number(Rational(-1, 2));
    if (detail::node(rhs)->equals(*detail::node(positive_half))) {
        principal = SymbolicExpr::multiply(SymbolicExpr::number(
            type == Type::Sin ? Rational(1, 6) : Rational(1, 3)), pi);
    } else if (detail::node(rhs)->equals(*detail::node(negative_half))) {
        principal = SymbolicExpr::multiply(SymbolicExpr::number(
            type == Type::Sin ? Rational(-1, 6) : Rational(2, 3)), pi);
    }
    return principal;
}

std::vector<ExprPtr> trigonometric_periodic_values(FunctionNode::FuncType type,
    const ExprPtr& rhs, const std::string& parameter) {
    using Type = FunctionNode::FuncType;
    auto pi = detail::make_expression_ptr(
        detail::make_node<VariableNode>("pi", true));
    auto integer_pi = SymbolicExpr::multiply(SymbolicExpr::variable(parameter), pi);
    auto twice_integer_pi = SymbolicExpr::multiply(SymbolicExpr::number(2), integer_pi);
    auto half_pi = SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)), pi);
    auto integer = exact_small_integer_node(detail::node(rhs->simplify()), -1, 1);
    std::vector<ExprPtr> values;
    if (type == Type::Tan) {
        auto principal = integer ?
            SymbolicExpr::multiply(SymbolicExpr::number(Rational(*integer, 4)), pi) :
            inverse_trig(Type::ArcTan, rhs);
        values.push_back(SymbolicExpr::add(principal, integer_pi));
    } else if (integer && *integer == 0) {
        values.push_back(type == Type::Sin ? integer_pi :
            SymbolicExpr::add(half_pi, integer_pi));
    } else if (integer && (*integer == 1 || *integer == -1)) {
        auto principal = type == Type::Sin ?
            SymbolicExpr::multiply(SymbolicExpr::number(*integer), half_pi) :
            (*integer == 1 ? SymbolicExpr::number(0) : pi);
        values.push_back(SymbolicExpr::add(principal, twice_integer_pi));
    } else {
        auto principal = trigonometric_principal_value(type, rhs, pi);
        values.push_back(SymbolicExpr::add(principal, twice_integer_pi));
        values.push_back(SymbolicExpr::add(
            type == Type::Sin ? subtract(pi, principal) :
                SymbolicExpr::multiply(SymbolicExpr::number(-1), principal),
            twice_integer_pi));
    }
    return values;
}

SolveResult invert_trigonometric(
    const InversePattern& pattern, const ExprPtr& rhs, const ExprPtr& original,
    const std::string& var, const FactsQuery& facts, ComputationContext& context) {
    auto affine = prepare_trigonometric_affine(pattern.inner, var, context);
    if (!affine) { return SolveResult::failure(affine.error()); }
    if (!affine.value()) { return conditional_preimage(original, var); }
    auto nonzero = detail::query_nonzero_value(detail::node(affine.value()->slope), facts, Domain::Real, context);
    if (!nonzero) { return SolveResult::failure(nonzero.error()); }
    auto real_coefficients = real_affine_coefficients(*affine.value(), nonzero.value(), facts, context);
    if (!real_coefficients) { return SolveResult::failure(real_coefficients.error()); }
    if (!real_coefficients.value()) { return conditional_preimage(original, var); }
    std::vector<ExprPtr> conditions;
    if (auto decided = trigonometric_target_range(
            pattern.func_type, rhs, original, var, facts, context, conditions)) {
        return std::move(*decided);
    }
    if (nonzero.value() == Tribool::False) {
        return constant_trigonometric_preimage(
            pattern.func_type, affine.value()->offset, rhs, original, var, facts, context);
    }
    if (nonzero.value() != Tribool::True) { return conditional_preimage(original, var); }
    auto fresh = capture_free_integer_parameter(original, var, context);
    if (!fresh) { return SolveResult::failure(fresh.error()); }
    const auto& parameter = fresh.value();
    auto values = trigonometric_periodic_values(pattern.func_type, rhs, parameter);
    ParametricSolutions solutions;
    for (const auto& value : values) {
        auto mapped = detail::simplify_expression(SymbolicExpr::divide(
            subtract(value, affine.value()->offset), affine.value()->slope), context);
        if (!mapped) { return SolveResult::failure(mapped.error()); }
        solutions.values.push_back(ParametricSolution{
            std::move(mapped.value()), {parameter}, conditions});
    }
    return SolutionSet{std::move(solutions)};
}

SolveResult invert_exponential(const ExprPtr& rhs, const ExprPtr& original,
    const std::string& var, const FactsQuery& facts, ComputationContext& context) {
    auto real = detail::query_real_value(detail::node(rhs), facts, context);
    if (!real) { return SolveResult::failure(real.error()); }
    auto positive = detail::query_positive_value(detail::node(rhs), facts, context);
    if (!positive) { return SolveResult::failure(positive.error()); }
    if (real.value() == Tribool::True && positive.value() == Tribool::False) {
        return SolutionSet{EmptySolutions{}};
    }
    if (positive.value() != Tribool::True) { return conditional_preimage(original, var); }
    auto value = SymbolicExpr::ln(rhs)->simplify();
    return SolutionSet{FiniteSolutions{{FiniteSolution{std::move(value), 1, {}}}}};
}

SolveResult invert_logarithm(const ExprPtr& rhs, const ExprPtr& original,
    const std::string& var, const FactsQuery& facts, ComputationContext& context) {
    auto real = detail::query_real_value(detail::node(rhs), facts, context);
    if (!real) { return SolveResult::failure(real.error()); }
    if (real.value() != Tribool::True) { return conditional_preimage(original, var); }
    auto value = SymbolicExpr::exp(rhs)->simplify();
    return SolutionSet{FiniteSolutions{{FiniteSolution{std::move(value), 1, {}}}}};
}

}

SolveResult invert_function(
    const InversePattern& pattern, const ExprPtr& original, const std::string& var,
    ComputationContext& context, const SolveOptions&) {
    auto step = context.consume_steps(1, "solve.inverse");
    if (!step) { return SolveResult::failure(step.error()); }
    std::optional<detail::AssumptionFacts> assumption_facts;
    if (context.assumptions()) { assumption_facts.emplace(*context.assumptions()); }
    const FactsQuery& facts = assumption_facts ? static_cast<const FactsQuery&>(*assumption_facts) :
        detail::no_facts();
    auto rhs = pattern.rhs;
    if (pattern.coefficient) {
        auto nonzero = detail::query_nonzero_value(detail::node(pattern.coefficient), facts, Domain::Real, context);
        if (!nonzero) { return SolveResult::failure(nonzero.error()); }
        if (nonzero.value() != Tribool::True) { return conditional_preimage(original, var); }
        rhs = SymbolicExpr::divide(rhs, pattern.coefficient)->simplify();
    }
    switch (pattern.func_type) {
    case FunctionNode::FuncType::Sin:
    case FunctionNode::FuncType::Cos:
    case FunctionNode::FuncType::Tan:
        return invert_trigonometric(pattern, rhs, original, var, facts, context);
    case FunctionNode::FuncType::Exp:
        return invert_exponential(rhs, original, var, facts, context);
    case FunctionNode::FuncType::Ln:
        return invert_logarithm(rhs, original, var, facts, context);
    default:
        return SolveResult::failure(CasErrc::Inconclusive,
            "Inverse requires a complete solution-set representation", "solve.inverse");
    }
}


}
