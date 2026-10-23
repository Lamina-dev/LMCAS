#include <stdexcept>
#include "internal/expression_transform.hpp"
#include "root_of_utils.hpp"
#include "solve_strategies.hpp"
#include "internal/calculus_utils_support.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS {

using namespace calculus_utils_detail;


static bool calculus_utils_contains_complex(const std::shared_ptr<const SymbolicNode>& node)
{
    if (!node) {
        return false;
    }

    if (std::dynamic_pointer_cast<const ComplexNode>(node)) {
        return true;
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        for (const auto& argument : function->arguments()) {
            if (calculus_utils_contains_complex(argument)) {
                return true;
            }
        }
        return false;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& operand : add->operands()) {
            if (calculus_utils_contains_complex(operand)) {
                return true;
            }
        }
        return false;
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& operand : multiply->operands()) {
            if (calculus_utils_contains_complex(operand)) {
                return true;
            }
        }
        return false;
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return calculus_utils_contains_complex(power->base()) ||
               calculus_utils_contains_complex(power->exponent());
    }
    return false;
}

static bool calculus_utils_is_negative_number(const std::shared_ptr<const NumberNode>& number)
{
    if (!number) {
        return false;
    }
    if (std::holds_alternative<BigInt>(number->value())) {
        return std::get<BigInt>(number->value()) < BigInt(0);
    }
    if (std::holds_alternative<Rational>(number->value())) {
        return std::get<Rational>(number->value()) < Rational(0);
    }
    return std::get<lmmc_real_t>(number->value()) < 0.0;
}

static bool is_nonreal_sqrt(const FunctionNode& function, ComputationContext& context)
{

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    if (function.type() != FunctionNode::FuncType::Sqrt ||
        function.arguments().size() != 1) {
        return false;
    }
    auto argument = detail::make_expression_ptr(detail::normalize_node(function.arguments()[0], transform_budget));
    auto number = argument
        ? std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(argument))
        : nullptr;
    return calculus_utils_is_negative_number(number);
}

static bool calculus_utils_contains_nonreal_sqrt(
    const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context);

static bool any_nonreal_sqrt(
    const std::vector<std::shared_ptr<const SymbolicNode>>& children, ComputationContext& context)
{
    for (const auto& child : children) {
        if (calculus_utils_contains_nonreal_sqrt(child, context)) {
            return true;
        }
    }
    return false;
}

static bool calculus_utils_contains_nonreal_sqrt(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context)
{
    if (!node) {
        return false;
    }

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (is_nonreal_sqrt(*function, context)) {
            return true;
        }
        return any_nonreal_sqrt(function->arguments(), context);
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        return any_nonreal_sqrt(add->operands(), context);
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return any_nonreal_sqrt(multiply->operands(), context);
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return calculus_utils_contains_nonreal_sqrt(power->base(), context) ||
               calculus_utils_contains_nonreal_sqrt(power->exponent(), context);
    }
    if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        return calculus_utils_contains_nonreal_sqrt(complex->real(), context) ||
               calculus_utils_contains_nonreal_sqrt(complex->imag(), context);
    }
    return false;
}

enum class CalculusRealCandidateStatus {
    Real,
    NonReal,
    Inconclusive
};

static CalculusRealCandidateStatus calculus_utils_extract_real_candidate(
    const std::shared_ptr<SymbolicExpr>& candidate,
    std::shared_ptr<SymbolicExpr>& real_candidate, ComputationContext& context)
{

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    real_candidate.reset();
    auto simplified_candidate = rootof_simplify(candidate);
    if (!simplified_candidate || !LMCAS::detail::node(simplified_candidate)) {
        return CalculusRealCandidateStatus::Inconclusive;
    }
    if (LMCAS::detail::contains_node_type<RootOfNode>(LMCAS::detail::node(simplified_candidate))) {
        return CalculusRealCandidateStatus::Inconclusive;
    }
    if (calculus_utils_contains_nonreal_sqrt(LMCAS::detail::node(simplified_candidate), context)) {
        return CalculusRealCandidateStatus::NonReal;
    }

    if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(
            LMCAS::detail::node(simplified_candidate))) {
        auto imag = detail::make_expression_ptr(detail::normalize_node(complex->imag(), transform_budget));
        if (!imag || !LMCAS::detail::node(imag)) {
            return CalculusRealCandidateStatus::Inconclusive;
        }
        if (imag->is_zero()) {
            auto real = detail::make_expression_ptr(detail::normalize_node(complex->real(), transform_budget));
            if (!real || !LMCAS::detail::node(real)) {
                return CalculusRealCandidateStatus::Inconclusive;
            }
            real_candidate = real;
            return CalculusRealCandidateStatus::Real;
        }
        return CalculusRealCandidateStatus::NonReal;
    }

    if (calculus_utils_contains_complex(LMCAS::detail::node(simplified_candidate))) {
        return CalculusRealCandidateStatus::Inconclusive;
    }

    real_candidate = simplified_candidate;
    return CalculusRealCandidateStatus::Real;
}
static bool contains_inverse_candidate(
    const std::vector<std::shared_ptr<SymbolicExpr>>& candidates,
    const std::shared_ptr<SymbolicExpr>& candidate)
{
    for (const auto& existing : candidates) {
        if (existing && LMCAS::detail::node(existing) &&
            LMCAS::detail::node(existing)->equals(*LMCAS::detail::node(candidate))) {
            return true;
        }
    }
    return false;
}

static SymbolicExprVectorResult verified_real_inverse_candidates(
    const std::vector<std::shared_ptr<SymbolicExpr>>& candidates,
    const std::shared_ptr<SymbolicExpr>& equation, const std::string& var,
    const std::string& operation, ComputationContext& context)
try {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
    std::vector<std::shared_ptr<SymbolicExpr>> real_candidates;
    real_candidates.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        std::shared_ptr<SymbolicExpr> real_candidate;
        auto status = calculus_utils_extract_real_candidate(candidate, real_candidate, context);
        if (status == CalculusRealCandidateStatus::Inconclusive) {
            return SymbolicExprVectorResult::failure(
                CasErrc::Inconclusive,
                "inverse derivative candidate reality could not be proven", operation);
        }
        if (status == CalculusRealCandidateStatus::NonReal) {
            continue;
        }
        auto residual = detail::make_expression_ptr(detail::substitute_node(detail::node(equation), var, detail::node(real_candidate), transform_budget));
        auto simplified = residual ? detail::make_expression_ptr(detail::normalize_node(detail::node(residual), transform_budget)) : nullptr;
        if (!simplified || !LMCAS::detail::node(simplified)) {
            return SymbolicExprVectorResult::failure(
                CasErrc::InternalInvariant,
                "inverse derivative candidate verification produced a null residual", operation);
        }
        if (!simplified->is_zero()) {
            continue;
        }
        if (!contains_inverse_candidate(real_candidates, real_candidate)) {
            real_candidates.push_back(real_candidate);
        }
    }
    return SymbolicExprVectorResult::success(std::move(real_candidates));
} catch (const CasError& error) {
    return SymbolicExprVectorResult::failure(error);
} catch (const std::bad_alloc&) {
    return SymbolicExprVectorResult::failure(CasErrc::ResourceLimit, "allocation failed", operation);
} catch (const std::length_error& error) {
    return SymbolicExprVectorResult::failure(CasErrc::ResourceLimit, error.what(), operation);
} catch (const std::exception& error) {
    return SymbolicExprVectorResult::failure(CasErrc::InternalInvariant, error.what(), operation);
}

static ExpressionResult inverse_branch_derivative(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::vector<std::shared_ptr<SymbolicExpr>>& candidates,
    const std::string& operation, ComputationContext& context)
try {

    detail::RewriteBudget transform_budget(context, (context).limits().max_recursion_depth,
                                            (context).limits().max_ast_nodes, "transform");
        if (candidates.empty()) {
            return ExpressionResult::failure(
                CasErrc::DomainError,
                "inverse point is outside the proven range",
                operation);
        }
        if (candidates.size() != 1) {
            return ExpressionResult::failure(
                CasErrc::Inconclusive,
                "inverse derivative requires a unique inverse branch",
                operation);
        }

        auto f_prime = detail::make_expression_ptr(detail::differentiate_node(detail::node(f), var, transform_budget));
        if (!f_prime || !LMCAS::detail::node(f_prime)) {
            return ExpressionResult::failure(
                CasErrc::Inconclusive,
                "inverse derivative could not construct the derivative",
                operation);
        }

        auto f_prime_at_x0 = detail::make_expression_ptr(detail::substitute_node(detail::node(f_prime), var, detail::node(candidates[0]), transform_budget));
        auto simplified_derivative = f_prime_at_x0 ? detail::make_expression_ptr(detail::normalize_node(detail::node(f_prime_at_x0), transform_budget)) : nullptr;
        if (!simplified_derivative || !LMCAS::detail::node(simplified_derivative)) {
            return ExpressionResult::failure(
                CasErrc::InternalInvariant,
                "inverse derivative substitution produced a null expression",
                operation);
        }

        if (simplified_derivative->is_zero()) {
            return ExpressionResult::failure(
                CasErrc::DomainError,
                "inverse derivative is undefined where f' is zero",
                operation);
        }

        auto result = SymbolicExpr::divide(SymbolicExpr::number(1), simplified_derivative);
        auto simplified = result ? detail::make_expression_ptr(detail::normalize_node(detail::node(result), transform_budget)) : nullptr;
        if (!simplified || !LMCAS::detail::node(simplified)) {
            return ExpressionResult::failure(
                CasErrc::InternalInvariant,
                "inverse derivative result construction failed",
                operation);
        }
        return ExpressionResult::success(simplified);
} catch (const CasError& error) {
    return ExpressionResult::failure(error);
} catch (const std::bad_alloc&) {
    return ExpressionResult::failure(CasErrc::ResourceLimit, "allocation failed", operation);
} catch (const std::length_error& error) {
    return ExpressionResult::failure(CasErrc::ResourceLimit, error.what(), operation);
} catch (const std::exception& error) {
    return ExpressionResult::failure(CasErrc::InternalInvariant, error.what(), operation);
}

static SymbolicExprVectorResult explicit_inverse_candidates(
    const FiniteSolutions& finite, const std::string& operation)
{
            std::vector<std::shared_ptr<SymbolicExpr>> values;
            values.reserve(finite.values.size());
            for (const auto& solution : finite.values) {
                if (!solution.value || !LMCAS::detail::node(solution.value)) {
                    return SymbolicExprVectorResult::failure(
                        CasErrc::InternalInvariant,
                        "checked solver returned a null inverse candidate",
                        operation);
                }
                auto candidate = rootof_simplify(solution.value);
                if (!candidate || !LMCAS::detail::node(candidate)) {
                    return SymbolicExprVectorResult::failure(
                        CasErrc::InternalInvariant,
                        "checked solver returned a malformed inverse candidate",
                        operation);
                }
                values.push_back(candidate);
            }
            return SymbolicExprVectorResult::success(std::move(values));
}

ExpressionResult inverse_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point,
    ComputationContext& context)
{
    const std::string operation = "inverse_derivative";
    auto input = calculus_utils_validate_expr_target(f, var, point, context, operation);
    if (!input) {
        return ExpressionResult::failure(input.error());
    }
    auto step = context.consume_steps(4, operation);
    if (!step) {
        return ExpressionResult::failure(step.error());
    }

    try {
        auto eq = SymbolicExpr::add(
            f, SymbolicExpr::multiply(SymbolicExpr::number(-1), point));

        std::vector<std::shared_ptr<SymbolicExpr>> candidates;
        auto inverse = inverse_function_checked(f, var, point, context);
        if (inverse) {
            candidates = inverse.value();
        } else {
            return ExpressionResult::failure(inverse.error());
        }
        auto verified = verified_real_inverse_candidates(candidates, eq, var, operation, context);
        if (!verified) {
            return ExpressionResult::failure(verified.error());
        }
        return inverse_branch_derivative(f, var, verified.value(), operation, context);
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(CasErrc::ResourceLimit,
                                           "inverse derivative allocation failed",
                                           operation);
    } catch (const std::exception& e) {
        return ExpressionResult::failure(CasErrc::InternalInvariant,
                                           e.what(), operation);
    }
}

ExpressionResult inverse_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point)
{
    ComputationContext context;
    return inverse_derivative_checked(f, var, point, context);
}

SymbolicExprVectorResult inverse_function_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& y, ComputationContext& context)
{
    const std::string operation = "inverse_function";
    auto input = calculus_utils_validate_expr_target(f, var, y, context, operation);
    if (!input) {
        return SymbolicExprVectorResult::failure(input.error());
    }
    auto step = context.consume_steps(4, operation);
    if (!step) {
        return SymbolicExprVectorResult::failure(step.error());
    }

    try {
        auto eq = SymbolicExpr::eq(f, y);

        auto solved = solve_equation(eq, var, context, SolveOptions{});
        if (!solved) {
            return SymbolicExprVectorResult::failure(solved.error());
        }

        const auto& solutions = solved.value();
        if (std::holds_alternative<EmptySolutions>(solutions)) {
            return SymbolicExprVectorResult::success({});
        }
        if (const auto* finite = std::get_if<FiniteSolutions>(&solutions)) {
            return explicit_inverse_candidates(*finite, operation);
        }

        return SymbolicExprVectorResult::failure(
            CasErrc::Inconclusive,
            "inverse equation is outside the finite exact support domain",
            operation);
    } catch (const std::bad_alloc&) {
        return SymbolicExprVectorResult::failure(CasErrc::ResourceLimit,
                                                 "inverse function allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return SymbolicExprVectorResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

SymbolicExprVectorResult inverse_function_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& y)
{
    ComputationContext context;
    return inverse_function_checked(f, var, y, context);
}

}
