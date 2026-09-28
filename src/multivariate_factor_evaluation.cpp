#include "internal/multivariate_factor_checked_support.hpp"
#include "internal/multivariate_factor_support.hpp"
#include <algorithm>
#include <map>
#include <new>
#include <optional>
#include <stdexcept>
#include <utility>

namespace LMCAS::multivariate_checked_detail {

static bool restore_specialization_unit(
    const MultiPoly& poly,
    const std::string& main_var,
    const std::map<std::string, Rational>& eval_points,
    std::vector<Polynomial<Rational>>& univariate_factors)
{
    if (univariate_factors.empty()) { return false; }
    const Polynomial<Rational> specialized = poly.eval(eval_points).to_univariate();
    if (specialized.variable_name != main_var) { return false; }
    Polynomial<Rational> product({Rational(1)}, main_var);
    for (const auto& factor : univariate_factors) {
        if (factor.variable_name != main_var || factor.is_zero()) { return false; }
        product = product * factor;
    }
    auto [unit, remainder] = specialized.div_mod(product);
    if (!remainder.is_zero() || unit.is_zero() || unit.degree() != 0) { return false; }
    if (!(product * unit == specialized)) { return false; }
    univariate_factors[0] = univariate_factors[0] * unit;
    return true;
}

struct LinearSamples {
    const MultiPoly& primitive;
    const std::string& main_variable;
    const std::vector<std::string>& auxiliary;
    const std::map<std::string, Rational>& zero_evaluation;
    const std::vector<Rational>& base_constants;
    ComputationContext& context;
};

static Result<bool> collect_linear_samples(
    const LinearSamples& input, std::vector<std::vector<Rational>>& samples) {
    const size_t factor_count = input.base_constants.size();
    for (const auto& variable : input.auxiliary) {
        auto sample = input.zero_evaluation;
        sample[variable] = Rational(1);
        MultiPoly evaluated = input.primitive.eval(sample);
        if (evaluated.degree(input.main_variable) != input.primitive.degree(input.main_variable)) {
            return Result<bool>::success(false);
        }
        auto factored = factor_univariate_bridge_checked(evaluated.to_univariate(), input.context);
        if (!factored) { return Result<bool>::failure(factored.error()); }
        if (factored.value().value.size() != factor_count) { return Result<bool>::success(false); }
        std::vector<Rational> constants;
        for (const auto& factor : factored.value().value) {
            if (factor.degree() != 1) { return Result<bool>::success(false); }
            constants.push_back(factor.make_monic().coeffs[0]);
        }
        samples.push_back(std::move(constants));
    }
    return Result<bool>::success(true);
}

struct LinearAssignments {
    const LinearSamples& input;
    const std::vector<std::vector<Rational>>& samples;
    std::vector<std::vector<size_t>> assignments;
    std::vector<MultiPoly> found_factors;
};

static MultiPoly linear_assignment_factor(const LinearAssignments& search, size_t factor_index) {
    const auto& variables = search.input.primitive.variables();
    std::vector<MultiPoly::Term> terms;
    Monomial main_monomial(variables.size(), 0);
    const size_t main_index = static_cast<size_t>(std::find(
        variables.begin(), variables.end(), search.input.main_variable) - variables.begin());
    main_monomial[main_index] = 1;
    terms.emplace_back(main_monomial, Rational(1));
    Monomial constant_monomial(variables.size(), 0);
    const auto& constant = search.input.base_constants[factor_index];
    if (!constant.is_zero()) terms.emplace_back(constant_monomial, constant);
    for (size_t variable_index = 0; variable_index < search.input.auxiliary.size(); ++variable_index) {
        Rational coefficient = search.samples[variable_index]
            [search.assignments[variable_index][factor_index]] - constant;
        if (coefficient.is_zero()) { continue; }
        Monomial monomial(variables.size(), 0);
        const size_t index = static_cast<size_t>(std::find(
            variables.begin(), variables.end(), search.input.auxiliary[variable_index]) - variables.begin());
        monomial[index] = 1;
        terms.emplace_back(std::move(monomial), coefficient);
    }
    return MultiPoly(std::move(terms), variables);
}

static Result<bool> verify_linear_assignment(LinearAssignments& search) {
    auto step = search.input.context.consume_steps(1, "factor_multivariate");
    if (!step) { return Result<bool>::failure(step.error()); }
    std::vector<MultiPoly> candidates;
    for (size_t i = 0; i < search.input.base_constants.size(); ++i) {
        candidates.push_back(linear_assignment_factor(search, i));
    }
    MultiPoly product(Rational(1), search.input.primitive.variables());
    for (const auto& candidate : candidates) product = product * candidate;
    MultiPoly unit;
    try {
        unit = search.input.primitive.exact_div(product);
    } catch (const std::runtime_error&) {
        return Result<bool>::success(false);
    }
    if (unit.is_zero() || !unit.is_constant() ||
        product * unit != search.input.primitive) { return Result<bool>::success(false); }
    search.found_factors = std::move(candidates);
    return Result<bool>::success(true);
}

static Result<bool> enumerate_linear_assignments(LinearAssignments& search, size_t index) {
    if (index == search.assignments.size()) { return verify_linear_assignment(search); }
    auto& permutation = search.assignments[index];
    std::sort(permutation.begin(), permutation.end());
    do {
        auto found = enumerate_linear_assignments(search, index + 1);
        if (!found || found.value()) { return found; }
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return Result<bool>::success(false);
}

static std::optional<MultiFactorCheckedResult> factor_linear_assignments(
    const MultiPoly& original, const LinearSamples& input) {
    std::vector<std::vector<Rational>> samples;
    auto collected = collect_linear_samples(input, samples);
    if (!collected) { return MultiFactorCheckedResult::failure(collected.error()); }
    if (!collected.value()) { return std::nullopt; }
    LinearAssignments search{input, samples,
        std::vector<std::vector<size_t>>(input.auxiliary.size(),
            std::vector<size_t>(input.base_constants.size())), {}};
    for (auto& assignment : search.assignments) {
        for (size_t i = 0; i < assignment.size(); ++i) assignment[i] = i;
    }
    auto found = enumerate_linear_assignments(search, 0);
    if (!found) { return MultiFactorCheckedResult::failure(found.error()); }
    if (!found.value()) { return std::nullopt; }
    return assemble_checked_factorization(original, search.found_factors,
        std::vector<int>(search.found_factors.size(), 1), Completeness::Complete, {});
}

static std::optional<MultiFactorCheckedResult> factor_linear_hensel(
    const MultiPoly& poly, const MultiPoly& primitive,
    const std::string& main_variable, const std::vector<std::string>& auxiliary_variables,
    const std::map<std::string, Rational>& zero_evaluation,
    MathResult<std::vector<Polynomial<Rational>>>& base_value, ComputationContext& context) {
    constexpr const char* operation = "factor_multivariate";
    const auto& variables = primitive.variables();
    auto hensel_step = context.consume_steps(
        static_cast<size_t>(
            std::max(1, primitive.degree(auxiliary_variables[0]))),
        operation);
    if (!hensel_step) {
        return MultiFactorCheckedResult::failure(hensel_step.error());
    }
    std::vector<Polynomial<Rational>> prepared = base_value.value;
    if (!restore_specialization_unit(
            primitive, main_variable, zero_evaluation, prepared)) { return std::nullopt; }
    auto lifted = multivariate_hensel_lift(
        primitive, prepared, auxiliary_variables[0], Rational(0),
        primitive.degree(auxiliary_variables[0]));
    MultiPoly lifted_product(Rational(1), variables);
    for (const auto& factor : lifted) {
        lifted_product = lifted_product * factor;
    }
    if (lifted.size() > 1 && lifted_product == primitive) {
        const bool all_linear = std::all_of(lifted.begin(), lifted.end(),
            [](const MultiPoly& factor) { return factor.total_degree() == 1; });
        return assemble_checked_factorization(
            poly, lifted, std::vector<int>(lifted.size(), 1),
            all_linear ? Completeness::Complete : Completeness::Inconclusive,
            all_linear ? std::string{} : "提升因子的不可约性尚未证明");
    }
    return std::nullopt;
}

std::optional<MultiFactorCheckedResult> factor_evaluated_linear(
    const MultiPoly& poly, const MultiPoly& primitive,
    const std::string& main_variable, ComputationContext& context) {
    const auto& variables = primitive.variables();
    std::vector<std::string> auxiliary_variables;
    for (const auto& variable : variables) {
        if (variable != main_variable) { auxiliary_variables.push_back(variable); }
    }
    std::map<std::string, Rational> zero_evaluation;
    for (const auto& variable : auxiliary_variables) {
        zero_evaluation[variable] = Rational(0);
    }
    MultiPoly evaluated = primitive.eval(zero_evaluation);
    if (evaluated.degree(main_variable) != primitive.degree(main_variable)) { return std::nullopt; }
    Polynomial<Rational> base_poly = evaluated.to_univariate();
    auto base_result =
        factor_univariate_bridge_checked(base_poly, context);
    if (!base_result) {
        return MultiFactorCheckedResult::failure(base_result.error());
    }
    auto base_value = std::move(base_result.value());
    bool all_linear = base_value.value.size() > 1;
    for (const auto& factor : base_value.value) {
        all_linear = all_linear && factor.degree() == 1;
    }
    if (!all_linear) { return std::nullopt; }
    if (auxiliary_variables.size() == 1) {
        auto hensel = factor_linear_hensel(poly, primitive, main_variable,
            auxiliary_variables, zero_evaluation, base_value, context);
        if (hensel) { return hensel; }
    }
    std::vector<Rational> base_constants;
    for (const auto& factor : base_value.value) {
        base_constants.push_back(factor.make_monic().coeffs[0]);
    }
    LinearSamples input{primitive, main_variable, auxiliary_variables,
                        zero_evaluation, base_constants, context};
    return factor_linear_assignments(poly, input);
}

}
