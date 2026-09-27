#include "solver.hpp"
#include "internal/solver_support.hpp"
#include "internal/exact_root.hpp"
#include "internal/normalization_utils.hpp"
#include "solve_strategies.hpp"
#include "internal/symbolic_ast.hpp"
#include "poly_utils.hpp"
#include "internal/expression_analysis.hpp"
#include "assumption_context.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <set>
#include <queue>
#include <unordered_map>
#include <optional>
#include <limits>

namespace LMCAS {
namespace {


/**
 * @brief 从解表达式提取 double 数值。
 * @return 纯数值表达式提取成功时写入 out_value 并返回 true。
 */
static bool try_get_numeric_value(const std::shared_ptr<SymbolicExpr>& expr, double& out_value) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return false;
    }
    auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr));
    if (!num) {
        return false;
    }

    return solver_detail::solver_number_value(*num, out_value);
}
static bool is_integer_value(double v) {
    if (!std::isfinite(v)) {
        return false;
    }
    int eq;
    lmmc_double_nearly_equal_tol(v, std::round(v), 1e-12, 1e-12, &eq);
    return eq != 0;
}

}

static bool domain_requires_real(Domain domain) {
    return domain == Domain::Real || domain == Domain::Algebraic ||
           domain == Domain::Rational || domain == Domain::Integer ||
           domain == Domain::Natural || domain == Domain::PositiveInt;
}

static bool numeric_domain_excludes(double value, Domain domain) {
    if (domain == Domain::Integer || domain == Domain::Natural || domain == Domain::PositiveInt) {
        if (!is_integer_value(value)) {
            return true;
        }
    }
    if (domain == Domain::Natural || domain == Domain::PositiveInt) {
        if (value < 0.0) {
            return true;
        }
    }
    if (domain == Domain::PositiveInt) {
        int equal;
        lmmc_double_nearly_equal(value, 0.0, &equal);
        if (equal || value <= 0.0) {
            return true;
        }
    }
    return false;
}

static bool numeric_sign_excludes(double value, const AssumptionContext& ctx,
                                  const std::string& variable) {
    bool exclude = false;
    if (ctx.has_sign(variable, Sign::NonNegative) && value < 0.0) {
        exclude = true;
    }
    if (ctx.has_sign(variable, Sign::Positive)) {
        int equal;
        lmmc_double_nearly_equal(value, 0.0, &equal);
        if (equal || value <= 0.0) {
            exclude = true;
        }
    }
    if (ctx.has_sign(variable, Sign::Negative) && value >= 0.0) {
        exclude = true;
    }
    if (ctx.has_sign(variable, Sign::NonPositive) && value > 0.0) {
        exclude = true;
    }
    return exclude;
}

static const NumberNode* exact_scalar_number(const SymbolicNode* node) {
    const auto* number = dynamic_cast<const NumberNode*>(node);
    if (number && (std::holds_alternative<BigInt>(number->value()) ||
                   std::holds_alternative<Rational>(number->value()))) {
        return number;
    }
    return nullptr;
}

static const NumberNode* exact_scalar_radicand(const SymbolicNode* node) {
    if (const auto* power = dynamic_cast<const PowerNode*>(node)) {
        const auto* exponent = exact_scalar_number(power->exponent().get());
        if (exponent && exact_number_as_rational(*exponent) == Rational(1, 2)) {
            return exact_scalar_number(power->base().get());
        }
    }
    if (const auto* function = dynamic_cast<const FunctionNode*>(node)) {
        if (function->type() == FunctionNode::FuncType::Sqrt &&
            function->arguments().size() == 1) {
            return exact_scalar_number(function->arguments()[0].get());
        }
    }
    return nullptr;
}

static bool collect_scalar_radical(const SymbolicNode* node, Rational& coefficient,
                                   const NumberNode*& radicand) {
    const auto* product = dynamic_cast<const MultiplyNode*>(node);
    if (!product) {
        radicand = exact_scalar_radicand(node);
        return radicand != nullptr;
    }
    for (const auto& factor : product->operands()) {
        if (const auto* number = exact_scalar_number(factor.get())) {
            coefficient = coefficient * exact_number_as_rational(*number);
        } else {
            if (radicand) {
                return false;
            }
            radicand = exact_scalar_radicand(factor.get());
            if (!radicand) {
                return false;
            }
        }
    }
    return radicand != nullptr;
}

static bool exact_rational_domain_excludes(const Rational& value, Domain domain) {
    if (domain == Domain::Rational) {
        return false;
    }
    if (!value.is_integer()) {
        return true;
    }
    if (domain == Domain::Natural) {
        return value < Rational(0);
    }
    if (domain == Domain::PositiveInt) {
        return value <= Rational(0);
    }
    return false;
}

static bool exact_radical_domain_excludes(const SymbolicNode* node, Domain domain) {
    if (domain != Domain::Rational && domain != Domain::Integer &&
        domain != Domain::Natural && domain != Domain::PositiveInt) {
        return false;
    }
    Rational coefficient(1);
    const NumberNode* radicand = nullptr;
    if (!collect_scalar_radical(node, coefficient, radicand)) {
        return false;
    }
    if (coefficient == Rational(0)) {
        return exact_rational_domain_excludes(coefficient, domain);
    }
    if (exact_number_as_rational(*radicand) < Rational(0)) {
        return false;
    }
    const auto root = normalization_exact_square_root(*radicand);
    if (!root) {
        return true;
    }
    return exact_rational_domain_excludes(
        coefficient * exact_number_as_rational(*root), domain);
}

static Result<bool> solution_is_excluded(const std::shared_ptr<SymbolicExpr>& solution,
                                        Domain domain, const AssumptionContext& ctx,
                                        const std::string& variable,
                                        ComputationContext& context) {
    if (domain_requires_real(domain)) {
        if (const auto* root = dynamic_cast<const RootOfNode*>(detail::node(solution).get())) {
            auto isolated = detail::isolate_exact_root(
                root->exact_id(), context, "solve_with_assumptions");
            if (!isolated) return Result<bool>::failure(isolated.error());
            if (std::holds_alternative<detail::ComplexIsolation>(isolated.value())) {
                return true;
            }
        } else if (solver_detail::contains_imaginary(solution)) {
            return true;
        }
    }
    if (exact_radical_domain_excludes(detail::node(solution).get(), domain)) {
        return true;
    }
    double value = 0.0;
    if (!try_get_numeric_value(solution, value)) {
        return false;
    }
    if (numeric_domain_excludes(value, domain)) {
        return true;
    }
    return numeric_sign_excludes(value, ctx, variable);
}

static AssumptionSolveResult solve_with_assumptions_impl(
    const std::shared_ptr<SymbolicExpr>& equation,
    const std::string& variable,
    const AssumptionContext* ctx,
    ComputationContext& context)
{
    auto solved = solve_finite_checked(
        std::const_pointer_cast<SymbolicExpr>(equation), variable, context);
    if (!solved) {
        return AssumptionSolveResult::failure(solved.error());
    }
    auto all_solutions = std::move(solved.value());

    if (!ctx) {
        return all_solutions;
    }

    Domain domain = ctx->get_domain(variable);

    std::vector<std::shared_ptr<SymbolicExpr>> filtered;
    filtered.reserve(all_solutions.size());

    bool has_sign_constraint = ctx->has_sign(variable, Sign::NonNegative) ||
                               ctx->has_sign(variable, Sign::Positive) ||
                               ctx->has_sign(variable, Sign::Negative) ||
                               ctx->has_sign(variable, Sign::NonPositive);

    if (domain == Domain::Complex && !has_sign_constraint) {
        return all_solutions;
    }

    for (const auto& sol : all_solutions) {
        if (!sol) {
            continue;
        }

        auto excluded = solution_is_excluded(sol, domain, *ctx, variable, context);
        if (!excluded) return AssumptionSolveResult::failure(excluded.error());
        if (!excluded.value()) {
            filtered.push_back(sol);
        }
    }

    return filtered;
}

AssumptionSolveResult solve_with_assumptions_checked(
    const std::shared_ptr<SymbolicExpr>& equation,
    const std::string& variable,
    const AssumptionContext* assumptions,
    ComputationContext& context)
{
    constexpr const char* operation = "solve_with_assumptions";
    if (!equation || variable.empty()) {
        return AssumptionSolveResult::failure(
            CasErrc::InvalidArgument,
            "assumption-aware solve requires an equation and variable",
            operation);
    }
    auto budget = context.consume_steps(1, operation);
    if (!budget) {
        return AssumptionSolveResult::failure(budget.error());
    }
    try {
        return solve_with_assumptions_impl(
            equation, variable, assumptions, context);
    } catch (const std::bad_alloc&) {
        return AssumptionSolveResult::failure(
            CasErrc::ResourceLimit,
            "allocation failed while solving with assumptions",
            operation);
    } catch (const std::exception& ex) {
        return AssumptionSolveResult::failure(
            CasErrc::InternalInvariant, ex.what(), operation);
    }
}

AssumptionSolveResult solve_with_assumptions_checked(
    const std::shared_ptr<SymbolicExpr>& equation,
    const std::string& variable,
    const AssumptionContext* assumptions)
{
    ComputationContext context;
    return solve_with_assumptions_checked(
        equation, variable, assumptions, context);
}


} // namespace LMCAS
