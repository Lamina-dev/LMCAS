#include "solve_mixed_transcendental.hpp"
#include "internal/mixed_transcendental_support.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/transcendental_solver_support.hpp"
#include "poly_utils.hpp"
#include "solve_polynomial.hpp"
#include "solve_transcendental.hpp"
#include "internal/assumption_facts.hpp"
#include "newton_raphson.hpp"
#include "numeric_evaluation.hpp"
#include "residual_verification.hpp"

#include <algorithm>
#include <cmath>
#include <new>
#include <stdexcept>

namespace LMCAS {
namespace {

constexpr const char* mixed_operation = "solve_mixed_transcendental_checked";

std::optional<CasError> validate_mixed_input(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var, const SolveOptions& opts) {
    if (!expr || !detail::node(expr)) {
        return CasError{
            CasErrc::InvalidArgument, "待求解表达式不能为空", mixed_operation};
    }
    if (var.empty()) {
        return CasError{
            CasErrc::InvalidArgument, "求解变量不能为空", mixed_operation};
    }
    if (!std::isfinite(opts.tolerance) || opts.tolerance <= 0.0 ||
        opts.max_newton_iterations <= 0 || opts.max_roots == 0 || opts.max_roots < -1) {
        return CasError{
            CasErrc::InvalidArgument, "求解选项包含无效界限", mixed_operation};
    }
    if (opts.has_search_interval &&
        (!std::isfinite(opts.search_lo) || !std::isfinite(opts.search_hi) ||
         opts.search_lo >= opts.search_hi)) {
        return CasError{
            CasErrc::InvalidArgument, "搜索区间必须为有限递增区间", mixed_operation};
    }
    return std::nullopt;
}

MixedTranscendentalResult classify_constant_equation(
    const std::shared_ptr<SymbolicExpr>& expr,
    bool unresolved_parameters, const SolveOptions& opts, ComputationContext& context) {
    using Roots = MathResult<std::vector<std::shared_ptr<SymbolicExpr>>>;
    if (unresolved_parameters) {
        return MixedTranscendentalResult::success(Roots{
            {}, Completeness::Inconclusive,
            "equation depends on unresolved parameters but not the solve variable"});
    }
    auto evaluated = evaluate_numeric(*expr, {}, context);
    if (!evaluated &&
        (evaluated.error().code == CasErrc::ResourceLimit ||
         evaluated.error().code == CasErrc::Cancelled ||
         evaluated.error().code == CasErrc::InternalInvariant)) {
        return MixedTranscendentalResult::failure(evaluated.error());
    }
    if (!evaluated || !std::isfinite(evaluated.value().value)) {
        return MixedTranscendentalResult::success(Roots{
            {}, Completeness::Inconclusive,
            "constant equation could not be evaluated in the real domain"});
    }
    if (std::abs(evaluated.value().value) <= opts.tolerance) {
        return MixedTranscendentalResult::success(Roots{
            {}, Completeness::Inconclusive,
            "equation is identically zero; every real value satisfies it"});
    }
    return MixedTranscendentalResult::success(Roots{{}, Completeness::Complete, {}});
}

bool is_mixed_terminal_error(CasErrc code) {
    return code == CasErrc::ResourceLimit ||
        code == CasErrc::Cancelled ||
        code == CasErrc::InternalInvariant;
}

class MixedRootCollection {
public:
    MixedRootCollection(const std::shared_ptr<SymbolicExpr>& expr,
                        const std::string& var, const SearchInterval& interval,
                        const SolveOptions& opts, ComputationContext& context)
        : expr_(expr), var_(var), interval_(interval), opts_(opts), context_(context) {}

    std::optional<CasError> solve(ComputationContext& context) {
        auto factored = factor_transcendental(expr_, var_, context);
        if (!factored) return factored.error();
        auto factors = std::move(factored.value());
        if (factors.empty()) {
            factors.push_back(expr_);
            complete_ = false;
        }
        for (const auto& factor : factors) {
            auto step = context.consume_steps(1, mixed_operation);
            if (!step) {
                return step.error();
            }
            auto failure = solve_factor(factor, context);
            if (failure) {
                return failure;
            }
            if (failure_) return failure_;
        }
        return std::nullopt;
    }

    MixedTranscendentalResult finish() {
        struct RepresentedRoot {
            lmmc_real_t value;
            std::shared_ptr<SymbolicExpr> expression;
            bool symbolic;
        };
        auto numeric_values = deduplicate_roots(roots_, opts_.tolerance, -1);
        std::vector<RepresentedRoot> values;
        values.reserve(symbolic_roots_.size() + numeric_values.size());
        for (auto& root : symbolic_roots_) {
            values.push_back({root.value, std::move(root.expression), true});
        }
        for (const auto value : numeric_values) {
            values.push_back(
                {value, SymbolicExpr::number(static_cast<double>(value)), false});
        }
        std::sort(values.begin(), values.end(),
            [](const RepresentedRoot& left, const RepresentedRoot& right) {
                if (left.value != right.value) { return left.value < right.value; }
                return left.symbolic && !right.symbolic;
            });
        const auto threshold = 10.0 * opts_.tolerance;
        std::vector<RepresentedRoot> unique;
        unique.reserve(values.size());
        for (auto& root : values) {
            if (!unique.empty() &&
                std::fabs(root.value - unique.back().value) < threshold) {
                if (root.symbolic && !unique.back().symbolic) {
                    unique.back() = std::move(root);
                }
                continue;
            }
            unique.push_back(std::move(root));
        }
        if (opts_.max_roots > 0 &&
            unique.size() > static_cast<std::size_t>(opts_.max_roots)) {
            unique.resize(static_cast<std::size_t>(opts_.max_roots));
            complete_ = false;
        }
        std::vector<std::shared_ptr<SymbolicExpr>> results;
        results.reserve(unique.size());
        for (auto& root : unique) {
            results.push_back(std::move(root.expression));
        }
        return MixedTranscendentalResult::success(
            MathResult<std::vector<std::shared_ptr<SymbolicExpr>>>{
                std::move(results),
                complete_ ? Completeness::Complete : Completeness::Inconclusive,
                complete_ ? std::string{} : "有界数值隔离未能证明排除全部未检测根"});
    }

private:
    void polish_value(lmmc_real_t& value, lmmc_real_t& residual, int& iterations) {
        auto polish = opts_;
        polish.has_initial_guess = true;
        polish.initial_guess = value;
        auto refined = solve_numeric_checked(expr_, var_, context_, polish);
        if (!refined) {
            if (is_mixed_terminal_error(refined.error().code)) {
                failure_ = refined.error();
            }
            return;
        }
        for (const auto& candidate : refined.value()) {
            if (candidate.value < interval_.lo || candidate.value > interval_.hi) { continue; }
            const auto candidate_residual = std::fabs(
                detail::mixed_evaluate_at(expr_, var_, candidate.value));
            if (std::isfinite(candidate_residual) && candidate_residual < residual) {
                value = candidate.value;
                residual = candidate_residual;
                iterations += candidate.iterations;
            }
        }
    }

    bool accept_value(lmmc_real_t value, int iterations) {
        if (!std::isfinite(value)) { complete_ = false; return false; }
        if (value < interval_.lo || value > interval_.hi) return true;
        auto residual = std::fabs(detail::mixed_evaluate_at(expr_, var_, value));
        if (std::isfinite(residual) && residual > opts_.tolerance) {
            polish_value(value, residual, iterations);
        }
        if (!std::isfinite(residual) || residual > opts_.tolerance) {
            complete_ = false;
            return false;
        }
        roots_.push_back({value, residual, iterations});
        return true;
    }

    bool handle_symbolic_evaluation_error(const CasError& error) {
        if (error.code == CasErrc::DomainError) {
            return true;
        }
        if (is_mixed_terminal_error(error.code)) {
            failure_ = error;
        }
        complete_ = false;
        return false;
    }

    std::optional<bool> accept_proved_symbolic(
        std::shared_ptr<SymbolicExpr> simplified,
        lmmc_real_t value) {
        auto bound = detail::substitute_raw(expr_, var_, simplified);
        std::optional<detail::AssumptionFacts> assumptions;
        if (context_.assumptions()) {
            assumptions.emplace(*context_.assumptions());
        }
        const FactsQuery& facts = assumptions
            ? static_cast<const FactsQuery&>(*assumptions)
            : detail::no_facts();
        auto defined = detail::query_definedness(
            detail::node(bound), facts, Domain::Real, context_);
        if (!defined) {
            if (is_mixed_terminal_error(defined.error().code)) {
                failure_ = defined.error();
                return false;
            }
            return std::nullopt;
        }
        if (defined.value() == Tribool::False) {
            return true;
        }
        if (defined.value() != Tribool::True) {
            return std::nullopt;
        }

        auto residual = check_zero_residual(bound, context_);
        if (!residual) {
            if (is_mixed_terminal_error(residual.error().code)) {
                failure_ = residual.error();
                return false;
            }
            return std::nullopt;
        }
        if (!std::holds_alternative<ProvedZeroResidual>(
                residual.value())) {
            return std::nullopt;
        }
        symbolic_roots_.push_back(
            SymbolicRoot{value, std::move(simplified)});
        return true;
    }

    bool accept_symbolic(const std::shared_ptr<SymbolicExpr>& root) {
        if (!root) {
            complete_ = false;
            return false;
        }
        auto simplified = root->simplify();
        if (!simplified || !detail::node(simplified)) {
            complete_ = false;
            return false;
        }
        auto evaluated = evaluate_numeric(*simplified, {}, context_);
        if (!evaluated) {
            return handle_symbolic_evaluation_error(evaluated.error());
        }
        const auto value = evaluated.value().value;
        if (!std::isfinite(value)) {
            complete_ = false;
            return false;
        }
        if (value < interval_.lo || value > interval_.hi) {
            return true;
        }
        auto accepted =
            accept_proved_symbolic(std::move(simplified), value);
        if (accepted) {
            return *accepted;
        }
        return accept_value(value, 0);
    }

    Result<bool> accept_finite_solutions(
        const FiniteSolutions& solutions) {
        bool handled = true;
        for (const auto& solution : solutions.values) {
            if (!solution.conditions.empty()) {
                complete_ = false;
                handled = false;
            } else if (!accept_symbolic(solution.value)) {
                handled = false;
            }
            if (failure_) {
                return Result<bool>::failure(*failure_);
            }
        }
        return handled;
    }

    Result<bool> solve_polynomial_factor(const ExprPtr& factor, ComputationContext& context) {
        if (contains_transcendental_of_var(factor, var_)) return false;
        auto polynomial = symbolic_to_poly<SymbolicPolyCoeff>(factor, var_);
        if (!polynomial) {
            if (polynomial.error().code != CasErrc::UnsupportedExpression)
                return Result<bool>::failure(polynomial.error());
            complete_ = false;
            return false;
        }
        auto solved = solve_equation(factor, var_, context, opts_);
        if (!solved) {
            if (solved.error().code != CasErrc::Inconclusive)
                return Result<bool>::failure(solved.error());
            complete_ = false;
            return false;
        }
        if (std::holds_alternative<EmptySolutions>(solved.value())) return true;
        auto finite = std::get_if<FiniteSolutions>(&solved.value());
        if (!finite) {
            complete_ = false;
            return false;
        }
        return accept_finite_solutions(*finite);
    }

    std::optional<CasError> check_constant_factor(const ExprPtr& factor) {
        std::optional<detail::AssumptionFacts> assumptions;
        if (context_.assumptions()) { assumptions.emplace(*context_.assumptions()); }
        const FactsQuery& facts = assumptions ? static_cast<const FactsQuery&>(*assumptions) : detail::no_facts();
        auto nonzero = detail::query_nonzero_value(detail::node(factor), facts, Domain::Real, context_);
        if (!nonzero) { return nonzero.error(); }
        if (nonzero.value() != Tribool::True) { complete_ = false; }
        return std::nullopt;
    }

    Result<bool> solve_transcendental_factor(
        const ExprPtr& factor, ComputationContext& context) {
        auto symbolic =
            solve_transcendental(factor, var_, context, opts_);
        if (!symbolic) {
            if (symbolic.error().code != CasErrc::Inconclusive) {
                return Result<bool>::failure(symbolic.error());
            }
            complete_ = false;
            return false;
        }
        if (std::holds_alternative<EmptySolutions>(symbolic.value())) {
            return true;
        }
        if (auto finite =
                std::get_if<FiniteSolutions>(&symbolic.value())) {
            return accept_finite_solutions(*finite);
        }
        complete_ = false;
        return false;
    }

    std::optional<CasError> solve_factor(
        const std::shared_ptr<SymbolicExpr>& factor,
        ComputationContext& context) {
        if (!factor || !detail::node(factor)) {
            complete_ = false;
            return std::nullopt;
        }
        if (!expression_depends_on_variable(detail::node(factor), var_)) {
            return check_constant_factor(factor);
        }
        auto polynomial = solve_polynomial_factor(factor, context);
        if (!polynomial) {
            return polynomial.error();
        }
        if (polynomial.value()) {
            return std::nullopt;
        }
        auto symbolic = solve_transcendental_factor(factor, context);
        if (!symbolic) {
            return symbolic.error();
        }
        if (symbolic.value()) {
            return std::nullopt;
        }

        auto numeric = detail::mixed_numerical_path(
            factor, var_, interval_, opts_, context, complete_);
        if (!numeric) {
            return std::move(numeric.error());
        }
        for (const auto& root : numeric.value()) {
            accept_value(root.value, root.iterations);
        }
        return std::nullopt;
    }

    struct SymbolicRoot {
        lmmc_real_t value;
        std::shared_ptr<SymbolicExpr> expression;
    };

    const std::shared_ptr<SymbolicExpr>& expr_;
    const std::string& var_;
    const SearchInterval& interval_;
    const SolveOptions& opts_;
    ComputationContext& context_;
    std::optional<CasError> failure_;
    std::vector<NumericRoot> roots_;
    std::vector<SymbolicRoot> symbolic_roots_;
    bool complete_ = true;
};

MixedTranscendentalResult solve_mixed_equation(
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& var,
    const SolveOptions& opts, ComputationContext& context) {
    auto step = context.consume_steps(1, mixed_operation);
    if (!step) {
        return MixedTranscendentalResult::failure(step.error());
    }
    const auto variables = free_variables(detail::node(expr));
    if (variables.find(var) == variables.end()) {
        return classify_constant_equation(expr, !variables.empty(), opts, context);
    }
    auto interval = determine_search_interval(expr, var, opts);
    if (!interval) {
        return MixedTranscendentalResult::failure(
            CasErrc::InvalidArgument, "无法确定有效搜索区间", mixed_operation);
    }
    if (!std::isfinite(interval->lo) || !std::isfinite(interval->hi) ||
        interval->lo >= interval->hi) {
        return MixedTranscendentalResult::failure(
            CasErrc::InvalidArgument, "搜索区间必须为有限递增区间", mixed_operation);
    }
    MixedRootCollection collection(expr, var, *interval, opts, context);
    auto failure = collection.solve(context);
    if (failure) {
        return MixedTranscendentalResult::failure(std::move(*failure));
    }
    return collection.finish();
}

}

MixedTranscendentalResult solve_mixed_transcendental_checked(
    const std::shared_ptr<SymbolicExpr>& expr, const std::string& var,
    const SolveOptions& opts, ComputationContext& context) {
    auto invalid = validate_mixed_input(expr, var, opts);
    if (invalid) {
        return MixedTranscendentalResult::failure(std::move(*invalid));
    }
    try {
        return solve_mixed_equation(expr, var, opts, context);
    } catch (const std::bad_alloc&) {
        return MixedTranscendentalResult::failure(
            CasErrc::ResourceLimit, "混合求解分配失败", mixed_operation);
    } catch (const std::exception& error) {
        return MixedTranscendentalResult::failure(
            CasErrc::InternalInvariant, error.what(), mixed_operation);
    }
}

MixedTranscendentalResult solve_mixed_transcendental_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var, const SolveOptions& opts) {
    ComputationContext context;
    return solve_mixed_transcendental_checked(expr, var, opts, context);
}

}
