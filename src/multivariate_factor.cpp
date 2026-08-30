/**
 * @file multivariate_factor.cpp
 * @brief 多元多项式因式分解器实现.
 */
#include "multivariate_factor.hpp"
#include "transcendental_factor.hpp"
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cmath>
#include <functional>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <new>
#include "internal/multivariate_factor_support.hpp"
namespace lamina {
static std::map<int, MultiPoly>
extract_coefficients(const MultiPoly& poly, const std::string& main_var)
{
    const auto& vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    std::vector<std::string> remaining_vars;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (static_cast<int>(i) != var_idx) remaining_vars.push_back(vars[i]);
    }
    std::map<int, MultiPoly> coeffs;
    if (var_idx < 0) { coeffs[0] = poly; return coeffs; }
    for (const auto& term : poly.terms()) {
        const Monomial& mono = term.first;
        int exp = (static_cast<size_t>(var_idx) < mono.size()) ? mono[var_idx] : 0;
        Monomial reduced_mono;
        reduced_mono.reserve(vars.size() - 1);
        for (size_t i = 0; i < mono.size(); ++i) {
            if (static_cast<int>(i) != var_idx) reduced_mono.push_back(mono[i]);
        }
        std::vector<MultiPoly::Term> t_vec = {{reduced_mono, term.second}};
        MultiPoly term_poly(std::move(t_vec), remaining_vars);
        if (coeffs.find(exp) == coeffs.end()) coeffs[exp] = term_poly;
        else coeffs[exp] = coeffs[exp] + term_poly;
    }
    return coeffs;
}
static bool divides_poly(const MultiPoly& dividend, const MultiPoly& divisor)
{
    if (divisor.is_zero()) return false;
    if (dividend.is_zero()) return true;
    try { dividend.exact_div(divisor); return true; }
    catch (...) { return false; }
}
static std::string choose_main_variable(const MultiPoly& a, const MultiPoly& b)
{
    const auto& vars = a.variables();
    if (vars.empty()) return "";
    std::string best_var = vars[0];
    int best_deg = a.degree(vars[0]) + b.degree(vars[0]);
    for (size_t i = 1; i < vars.size(); ++i) {
        int deg = a.degree(vars[i]) + b.degree(vars[i]);
        if (deg > best_deg) { best_deg = deg; best_var = vars[i]; }
    }
    return best_var;
}
static Rational rational_gcd(const Rational& a, const Rational& b)
{
    if (a.is_zero()) return b.abs();
    if (b.is_zero()) return a.abs();
    BigInt num_gcd = BigInt::gcd(a.get_numerator().Abs(), b.get_numerator().Abs());
    BigInt den_gcd = BigInt::gcd(a.get_denominator(), b.get_denominator());
    BigInt den_lcm = a.get_denominator() * b.get_denominator() / den_gcd;
    return Rational(num_gcd, den_lcm);
}
static MultiPoly lagrange_interpolate(const std::vector<Rational>& points,
    const std::vector<MultiPoly>& values, const std::string& main_var,
    const std::vector<std::string>& all_vars)
{
    size_t n = points.size();
    if (n == 0) return MultiPoly(Rational(0), all_vars);
    int var_idx = -1;
    for (size_t i = 0; i < all_vars.size(); ++i) {
        if (all_vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    MultiPoly result(Rational(0), all_vars);
    for (size_t i = 0; i < n; ++i) {
        if (values[i].is_zero()) continue;
        Rational denom(1);
        for (size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            denom = denom * (points[i] - points[j]);
        }
        if (denom.is_zero()) continue;
        std::vector<Rational> basis_coeffs = {Rational(1)};
        for (size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            std::vector<Rational> nc(basis_coeffs.size() + 1, Rational(0));
            for (size_t k = 0; k < basis_coeffs.size(); ++k) {
                nc[k + 1] = nc[k + 1] + basis_coeffs[k];
                nc[k] = nc[k] - points[j] * basis_coeffs[k];
            }
            basis_coeffs = std::move(nc);
        }
        for (auto& c : basis_coeffs) c = c / denom;
        for (size_t k = 0; k < basis_coeffs.size(); ++k) {
            if (basis_coeffs[k].is_zero()) continue;
            MultiPoly scaled = values[i] * basis_coeffs[k];
            if (scaled.is_zero()) continue;
            for (const auto& term : scaled.terms()) {
                Monomial full_mono(all_vars.size(), 0);
                size_t ri = 0;
                for (size_t vi = 0; vi < all_vars.size(); ++vi) {
                    if (static_cast<int>(vi) == var_idx) {
                        full_mono[vi] = static_cast<int>(k);
                    } else {
                        if (ri < term.first.size()) full_mono[vi] = term.first[ri];
                        ++ri;
                    }
                }
                std::vector<MultiPoly::Term> t = {{full_mono, term.second}};
                result = result + MultiPoly(std::move(t), all_vars);
            }
        }
    }
    return result;
}

/**
 * @internal
 * @brief 通过 Berlekamp、Hensel 提升和有界组合枚举执行受检一元分解。
 */
UnivariateFactorResult factor_univariate_bridge_checked(
    const Polynomial<Rational>& poly,
    ComputationContext& context)
{
    constexpr const char* operation = "factor_multivariate";
    auto step = context.consume_steps(1, operation);
    if (!step) return UnivariateFactorResult::failure(step.error());
    if (poly.is_zero() || poly.degree() <= 0) {
        return UnivariateFactorResult::success(
            MathResult<std::vector<Polynomial<Rational>>>{
                poly.is_zero() ? std::vector<Polynomial<Rational>>{}
                               : std::vector<Polynomial<Rational>>{poly},
                Completeness::Complete, {}});
    }
    if (poly.degree() == 1) {
        return UnivariateFactorResult::success(
            MathResult<std::vector<Polynomial<Rational>>>{
                {poly.make_monic()}, Completeness::Complete, {}});
    }

    try {
        TfSquareFreeResult sqf = tf_square_free(poly);
        Polynomial<Rational> work_poly = sqf.square_free.make_monic();
        static const int64_t trial_primes[] =
            {3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
        BerlekampResult berlekamp;
        bool found_valid_reduction = false;

        for (int64_t prime : trial_primes) {
            auto prime_step = context.consume_steps(1, operation);
            if (!prime_step) {
                return UnivariateFactorResult::failure(prime_step.error());
            }
            try {
                berlekamp = berlekamp_factor(work_poly, prime);
                if (berlekamp.prime > 0 && !berlekamp.factors.empty()) {
                    found_valid_reduction = true;
                    break;
                }
            } catch (const std::exception&) {
                // 该素数可能使首项系数消失；继续尝试下一个确定性素数。
            }
        }

        if (!found_valid_reduction) {
            return UnivariateFactorResult::success(
                MathResult<std::vector<Polynomial<Rational>>>{
                    {poly.make_monic()}, Completeness::Inconclusive,
                    "未找到保持次数的有效模素数"});
        }

        std::vector<Polynomial<Rational>> factors;
        Completeness completeness = Completeness::Complete;
        std::string reason;
        if (berlekamp.factors.size() <= 1) {
            factors.push_back(work_poly.make_monic());
        } else {
            BigInt denominator_lcm(1);
            for (const Rational& coefficient : work_poly.coeffs) {
                denominator_lcm = BigInt::lcm(
                    denominator_lcm, coefficient.get_denominator());
            }
            std::vector<BigInt> integer_coefficients;
            integer_coefficients.reserve(work_poly.coeffs.size());
            for (const Rational& coefficient : work_poly.coeffs) {
                integer_coefficients.push_back(
                    coefficient.get_numerator() *
                    (denominator_lcm / coefficient.get_denominator()));
            }
            Polynomial<BigInt> integer_poly(
                integer_coefficients, work_poly.variable_name);

            BigInt max_coefficient(0);
            for (const BigInt& coefficient : integer_coefficients) {
                if (coefficient.Abs() > max_coefficient) {
                    max_coefficient = coefficient.Abs();
                }
            }
            BigInt target = max_coefficient * BigInt(2);
            for (int i = 0; i < work_poly.degree(); ++i) {
                target = target * BigInt(2);
            }
            const int64_t prime = berlekamp.prime;
            int lift_bound = 1;
            BigInt modulus(static_cast<long long>(prime));
            while (modulus <= target && lift_bound < 100) {
                auto lift_step = context.consume_steps(1, operation);
                if (!lift_step) {
                    return UnivariateFactorResult::failure(lift_step.error());
                }
                modulus = modulus * BigInt(static_cast<long long>(prime));
                ++lift_bound;
            }

            auto lifted = hensel_lift_checked(
                integer_poly, berlekamp.factors, prime, lift_bound);
            if (!lifted) {
                if (lifted.error().code == CasErrc::InvalidArgument ||
                    lifted.error().code == CasErrc::Inconclusive) {
                    return UnivariateFactorResult::success(
                        MathResult<std::vector<Polynomial<Rational>>>{
                            {poly.make_monic()}, Completeness::Inconclusive,
                            "当前模素数的提升前置条件不成立"});
                }
                return UnivariateFactorResult::failure(lifted.error());
            }
            auto combined = zassenhaus_combine_checked(
                work_poly, lifted.value(), modulus, context);
            if (!combined) {
                if (combined.error().code == CasErrc::InvalidArgument ||
                    combined.error().code == CasErrc::Inconclusive) {
                    return UnivariateFactorResult::success(
                        MathResult<std::vector<Polynomial<Rational>>>{
                            {poly.make_monic()}, Completeness::Inconclusive,
                            "提升因子无法完成精确有理重构"});
                }
                return UnivariateFactorResult::failure(combined.error());
            }
            factors = std::move(combined.value().value);
            completeness = combined.value().completeness;
            reason = std::move(combined.value().reason);
        }

        if (sqf.had_repeated_factors &&
            sqf.repeated_factor.degree() >= 1) {
            auto repeated = factor_univariate_bridge_checked(
                sqf.repeated_factor, context);
            if (!repeated) {
                return UnivariateFactorResult::failure(repeated.error());
            }
            if (repeated.value().completeness == Completeness::Inconclusive) {
                completeness = Completeness::Inconclusive;
                if (reason.empty()) reason = repeated.value().reason;
            }
            auto repeated_factors = std::move(repeated.value().value);
            factors.insert(factors.end(),
                           std::make_move_iterator(repeated_factors.begin()),
                           std::make_move_iterator(repeated_factors.end()));
        }

        Polynomial<Rational> product(
            {Rational(1)}, poly.variable_name);
        for (const auto& factor : factors) product = product * factor;
        if (!(product == poly.make_monic())) {
            return UnivariateFactorResult::failure(
                CasErrc::InternalInvariant,
                "一元因子未能精确重构输入多项式", operation);
        }
        return UnivariateFactorResult::success(
            MathResult<std::vector<Polynomial<Rational>>>{
                std::move(factors), completeness, std::move(reason)});
    } catch (const std::bad_alloc&) {
        return UnivariateFactorResult::failure(
            CasErrc::ResourceLimit, "一元分解分配失败", operation);
    } catch (const std::exception& error) {
        return UnivariateFactorResult::success(
            MathResult<std::vector<Polynomial<Rational>>>{
                {poly.make_monic()}, Completeness::Inconclusive,
                std::string("一元分解路径未覆盖：") + error.what()});
    }
}


static void precompute_leading_coefficients(
    const MultiPoly& poly,
    const std::string& main_var,
    const std::map<std::string, Rational>& eval_points,
    std::vector<Polynomial<Rational>>& univariate_factors);

namespace {

struct FactorRecursionGuard {
    ComputationContext& context;
    ~FactorRecursionGuard() { context.leave_recursion(); }
};

static MultiPoly embed_univariate_factor(
    const Polynomial<Rational>& factor,
    const std::vector<std::string>& variables,
    const std::string& main_variable)
{
    size_t main_index = variables.size();
    for (size_t i = 0; i < variables.size(); ++i) {
        if (variables[i] == main_variable) {
            main_index = i;
            break;
        }
    }
    std::vector<MultiPoly::Term> terms;
    for (size_t degree = 0; degree < factor.coeffs.size(); ++degree) {
        if (factor.coeffs[degree].is_zero()) continue;
        Monomial monomial(variables.size(), 0);
        monomial[main_index] = static_cast<int>(degree);
        terms.emplace_back(std::move(monomial), factor.coeffs[degree]);
    }
    return MultiPoly(std::move(terms), variables);
}

static Result<MathResult<MultiFactorResult>> assemble_checked_factorization(
    const MultiPoly& original,
    const std::vector<MultiPoly>& raw_factors,
    const std::vector<int>& raw_multiplicities,
    Completeness completeness,
    std::string reason)
{
    constexpr const char* operation = "factor_multivariate";
    if (original.is_zero()) {
        if (!raw_factors.empty()) {
            return MultiFactorCheckedResult::failure(
                CasErrc::InternalInvariant, "零多项式收到非空因子列表", operation);
        }
        return MultiFactorCheckedResult::success(
            MathResult<MultiFactorResult>{
                MultiFactorResult{Rational(0), {}, {}},
                completeness, std::move(reason)});
    }

    MultiFactorResult result;
    for (size_t i = 0; i < raw_factors.size(); ++i) {
        const MultiPoly& raw = raw_factors[i];
        const int multiplicity =
            i < raw_multiplicities.size() ? raw_multiplicities[i] : 1;
        if (multiplicity <= 0) {
            return MultiFactorCheckedResult::failure(
                CasErrc::InternalInvariant, "因子重数必须为正", operation);
        }
        if (raw.is_zero() || raw.is_constant()) continue;
        MultiPoly factor = raw.make_primitive();
        if (!factor.terms().empty() &&
            factor.terms()[0].second < Rational(0)) {
            factor = factor * Rational(-1);
        }
        auto existing = std::find(result.factors.begin(),
                                  result.factors.end(), factor);
        if (existing == result.factors.end()) {
            result.factors.push_back(std::move(factor));
            result.multiplicities.push_back(multiplicity);
        } else {
            const size_t index =
                static_cast<size_t>(existing - result.factors.begin());
            result.multiplicities[index] += multiplicity;
        }
    }

    MultiPoly product(Rational(1), original.variables());
    for (size_t i = 0; i < result.factors.size(); ++i) {
        for (int power = 0; power < result.multiplicities[i]; ++power) {
            product = product * result.factors[i];
        }
    }
    MultiPoly quotient = original.exact_div(product);
    if (!quotient.is_constant() || quotient.terms().empty()) {
        return MultiFactorCheckedResult::failure(
            CasErrc::InternalInvariant,
            "因子乘积与输入只差常数的不变量失效", operation);
    }
    result.constant = quotient.terms()[0].second;
    if (product * result.constant != original) {
        return MultiFactorCheckedResult::failure(
            CasErrc::InternalInvariant, "最终因子未精确重构输入", operation);
    }
    return MultiFactorCheckedResult::success(
        MathResult<MultiFactorResult>{
            std::move(result), completeness, std::move(reason)});
}

static MultiFactorCheckedResult factor_multivariate_impl(
    const MultiPoly& poly,
    ComputationContext& context)
{
    constexpr const char* operation = "factor_multivariate";
    auto recursion = context.enter_recursion(operation);
    if (!recursion) return MultiFactorCheckedResult::failure(recursion.error());
    FactorRecursionGuard recursion_guard{context};

    if (poly.is_zero()) {
        return assemble_checked_factorization(
            poly, {}, {}, Completeness::Complete, {});
    }
    if (poly.is_constant()) {
        return assemble_checked_factorization(
            poly, {}, {}, Completeness::Complete, {});
    }
    const auto& variables = poly.variables();
    if (variables.empty()) {
        return MultiFactorCheckedResult::failure(
            CasErrc::InternalInvariant,
            "非常数多项式缺少变量表", operation);
    }

    std::string main_variable = variables.front();
    for (const auto& variable : variables) {
        if (poly.degree(variable) > poly.degree(main_variable)) {
            main_variable = variable;
        }
    }
    if (poly.total_degree() <= 1) {
        return assemble_checked_factorization(
            poly, {poly}, {1}, Completeness::Complete, {});
    }

    auto [common_monomial, quotient] =
        detail::extract_common_monomial(poly);
    bool has_common_monomial = false;
    for (int exponent : common_monomial) {
        if (exponent > 0) {
            has_common_monomial = true;
            break;
        }
    }
    if (has_common_monomial) {
        auto sub = factor_multivariate_impl(quotient, context);
        if (!sub) return sub;
        std::vector<MultiPoly> factors;
        std::vector<int> multiplicities;
        for (size_t i = 0; i < common_monomial.size(); ++i) {
            if (common_monomial[i] <= 0) continue;
            Monomial monomial(variables.size(), 0);
            monomial[i] = 1;
            factors.emplace_back(
                std::vector<MultiPoly::Term>{{monomial, Rational(1)}},
                variables);
            multiplicities.push_back(common_monomial[i]);
        }
        auto sub_value = std::move(sub.value());
        factors.insert(factors.end(),
                       sub_value.value.factors.begin(),
                       sub_value.value.factors.end());
        multiplicities.insert(multiplicities.end(),
                              sub_value.value.multiplicities.begin(),
                              sub_value.value.multiplicities.end());
        return assemble_checked_factorization(
            poly, factors, multiplicities, sub_value.completeness,
            std::move(sub_value.reason));
    }

    MultiPoly primitive = poly.make_primitive();
    if (auto difference =
            detail::detect_difference_of_squares(primitive)) {
        std::vector<MultiPoly> factors;
        std::vector<int> multiplicities;
        Completeness completeness = Completeness::Complete;
        std::string reason;
        for (MultiPoly candidate :
             {difference->first + difference->second,
              difference->first - difference->second}) {
            auto sub = factor_multivariate_impl(candidate, context);
            if (!sub) return sub;
            auto sub_value = std::move(sub.value());
            if (sub_value.completeness == Completeness::Inconclusive) {
                completeness = Completeness::Inconclusive;
                if (reason.empty()) reason = sub_value.reason;
            }
            factors.insert(factors.end(),
                           sub_value.value.factors.begin(),
                           sub_value.value.factors.end());
            multiplicities.insert(
                multiplicities.end(),
                sub_value.value.multiplicities.begin(),
                sub_value.value.multiplicities.end());
        }
        return assemble_checked_factorization(
            poly, factors, multiplicities, completeness, std::move(reason));
    }

    if (auto binomial = detail::detect_binomial_power(primitive)) {
        const auto& [left_variable, right_variable, exponent] = *binomial;
        (void)exponent;
        size_t left_index = 0;
        size_t right_index = 0;
        for (size_t i = 0; i < variables.size(); ++i) {
            if (variables[i] == left_variable) left_index = i;
            if (variables[i] == right_variable) right_index = i;
        }
        Monomial left_monomial(variables.size(), 0);
        Monomial right_monomial(variables.size(), 0);
        left_monomial[left_index] = 1;
        right_monomial[right_index] = 1;
        MultiPoly linear(
            std::vector<MultiPoly::Term>{
                {left_monomial, Rational(1)},
                {right_monomial, Rational(-1)}},
            variables);
        MultiPoly remaining = primitive.exact_div(linear);
        auto sub = factor_multivariate_impl(remaining, context);
        if (!sub) return sub;
        auto sub_value = std::move(sub.value());
        std::vector<MultiPoly> factors = {linear};
        std::vector<int> multiplicities = {1};
        factors.insert(factors.end(),
                       sub_value.value.factors.begin(),
                       sub_value.value.factors.end());
        multiplicities.insert(multiplicities.end(),
                              sub_value.value.multiplicities.begin(),
                              sub_value.value.multiplicities.end());
        return assemble_checked_factorization(
            poly, factors, multiplicities, sub_value.completeness,
            std::move(sub_value.reason));
    }

    if (primitive.is_univariate()) {
        Polynomial<Rational> univariate = primitive.to_univariate();
        auto factored =
            factor_univariate_bridge_checked(univariate, context);
        if (!factored) {
            return MultiFactorCheckedResult::failure(factored.error());
        }
        auto factor_value = std::move(factored.value());
        std::vector<MultiPoly> factors;
        factors.reserve(factor_value.value.size());
        for (const auto& factor : factor_value.value) {
            factors.push_back(embed_univariate_factor(
                factor, variables, main_variable));
        }
        return assemble_checked_factorization(
            poly, factors, std::vector<int>(factors.size(), 1),
            factor_value.completeness, std::move(factor_value.reason));
    }

    std::vector<std::string> auxiliary_variables;
    for (const auto& variable : variables) {
        if (variable != main_variable) auxiliary_variables.push_back(variable);
    }
    std::map<std::string, Rational> zero_evaluation;
    for (const auto& variable : auxiliary_variables) {
        zero_evaluation[variable] = Rational(0);
    }
    MultiPoly evaluated = primitive.eval(zero_evaluation);
    if (evaluated.degree(main_variable) == primitive.degree(main_variable)) {
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

        if (all_linear && auxiliary_variables.size() == 1) {
            auto hensel_step = context.consume_steps(
                static_cast<size_t>(
                    std::max(1, primitive.degree(auxiliary_variables[0]))),
                operation);
            if (!hensel_step) {
                return MultiFactorCheckedResult::failure(hensel_step.error());
            }
            std::vector<Polynomial<Rational>> prepared = base_value.value;
            precompute_leading_coefficients(
                primitive, main_variable, zero_evaluation, prepared);
            auto lifted = multivariate_hensel_lift(
                primitive, prepared, auxiliary_variables[0], Rational(0),
                primitive.degree(auxiliary_variables[0]));
            MultiPoly lifted_product(Rational(1), variables);
            for (const auto& factor : lifted) {
                lifted_product = lifted_product * factor;
            }
            if (lifted.size() > 1 && lifted_product == primitive) {
                return assemble_checked_factorization(
                    poly, lifted, std::vector<int>(lifted.size(), 1),
                    base_value.completeness, std::move(base_value.reason));
            }
        }

        if (all_linear) {
            const size_t factor_count = base_value.value.size();
            std::vector<Rational> base_constants;
            for (const auto& factor : base_value.value) {
                base_constants.push_back(factor.make_monic().coeffs[0]);
            }
            std::vector<std::vector<Rational>> sample_constants;
            bool samples_valid = true;
            for (const auto& sampled_variable : auxiliary_variables) {
                std::map<std::string, Rational> sample = zero_evaluation;
                sample[sampled_variable] = Rational(1);
                MultiPoly sample_poly = primitive.eval(sample);
                if (sample_poly.degree(main_variable) !=
                    primitive.degree(main_variable)) {
                    samples_valid = false;
                    break;
                }
                auto sample_result = factor_univariate_bridge_checked(
                    sample_poly.to_univariate(), context);
                if (!sample_result) {
                    return MultiFactorCheckedResult::failure(
                        sample_result.error());
                }
                if (sample_result.value().value.size() != factor_count) {
                    samples_valid = false;
                    break;
                }
                std::vector<Rational> constants;
                for (const auto& factor : sample_result.value().value) {
                    if (factor.degree() != 1) {
                        samples_valid = false;
                        break;
                    }
                    constants.push_back(factor.make_monic().coeffs[0]);
                }
                if (!samples_valid) break;
                sample_constants.push_back(std::move(constants));
            }

            if (samples_valid) {
                std::vector<std::vector<size_t>> assignments(
                    auxiliary_variables.size(),
                    std::vector<size_t>(factor_count));
                for (auto& assignment : assignments) {
                    for (size_t i = 0; i < factor_count; ++i) {
                        assignment[i] = i;
                    }
                }
                std::vector<MultiPoly> found_factors;
                std::optional<CasError> enumeration_error;
                std::function<bool(size_t)> enumerate =
                    [&](size_t variable_index) -> bool {
                    if (variable_index < assignments.size()) {
                        auto& permutation = assignments[variable_index];
                        std::sort(permutation.begin(), permutation.end());
                        do {
                            if (enumerate(variable_index + 1)) return true;
                            if (enumeration_error) return false;
                        } while (std::next_permutation(
                            permutation.begin(), permutation.end()));
                        return false;
                    }

                    auto candidate_step =
                        context.consume_steps(1, operation);
                    if (!candidate_step) {
                        enumeration_error = candidate_step.error();
                        return false;
                    }
                    std::vector<MultiPoly> candidates;
                    for (size_t factor_index = 0;
                         factor_index < factor_count; ++factor_index) {
                        std::vector<MultiPoly::Term> terms;
                        Monomial main_monomial(variables.size(), 0);
                        size_t main_index = static_cast<size_t>(
                            std::find(variables.begin(), variables.end(),
                                      main_variable) - variables.begin());
                        main_monomial[main_index] = 1;
                        terms.emplace_back(main_monomial, Rational(1));
                        Monomial constant_monomial(variables.size(), 0);
                        if (!base_constants[factor_index].is_zero()) {
                            terms.emplace_back(
                                constant_monomial,
                                base_constants[factor_index]);
                        }
                        for (size_t variable_index = 0;
                             variable_index < auxiliary_variables.size();
                             ++variable_index) {
                            Rational coefficient =
                                sample_constants[variable_index]
                                    [assignments[variable_index][factor_index]]
                                - base_constants[factor_index];
                            if (coefficient.is_zero()) continue;
                            Monomial monomial(variables.size(), 0);
                            const size_t index = static_cast<size_t>(
                                std::find(variables.begin(), variables.end(),
                                          auxiliary_variables[variable_index])
                                - variables.begin());
                            monomial[index] = 1;
                            terms.emplace_back(
                                std::move(monomial), coefficient);
                        }
                        candidates.emplace_back(std::move(terms), variables);
                    }
                    MultiPoly product(Rational(1), variables);
                    for (const auto& factor : candidates) {
                        product = product * factor;
                    }
                    if (product == primitive) {
                        found_factors = std::move(candidates);
                        return true;
                    }
                    return false;
                };

                const bool found = enumerate(0);
                if (enumeration_error) {
                    return MultiFactorCheckedResult::failure(
                        std::move(*enumeration_error));
                }
                if (found) {
                    return assemble_checked_factorization(
                        poly, found_factors,
                        std::vector<int>(found_factors.size(), 1),
                        Completeness::Complete, {});
                }
            }
        }
    }

    if (auto homogeneous =
            detail::factor_homogeneous_bivariate(primitive)) {
        return assemble_checked_factorization(
            poly, homogeneous->factors, homogeneous->multiplicities,
            Completeness::Complete, {});
    }

    return assemble_checked_factorization(
        poly, {poly}, {1}, Completeness::Inconclusive,
        "当前受检算法无法证明该多项式不可约");
}

} // namespace

MultiFactorCheckedResult factor_multivariate_checked(
    const MultiPoly& poly,
    ComputationContext& context)
{
    try {
        return factor_multivariate_impl(poly, context);
    } catch (const std::bad_alloc&) {
        return MultiFactorCheckedResult::failure(
            CasErrc::ResourceLimit, "多元分解分配失败",
            "factor_multivariate");
    } catch (const std::exception& error) {
        return assemble_checked_factorization(
            poly, {poly}, {1}, Completeness::Inconclusive,
            std::string("多元分解路径未覆盖：") + error.what());
    }
}

MultiFactorCheckedResult factor_multivariate_checked(const MultiPoly& poly)
{
    ComputationContext context;
    return factor_multivariate_checked(poly, context);
}
MultiPoly multivariate_content(const MultiPoly& poly, const std::string& main_var)
{
    if (poly.is_zero()) return MultiPoly();
    /// 若 main_var 不在变量列表中,或多项式在 main_var 上次数为 0,容度就是多项式本身
    const auto& vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    if (var_idx < 0) return poly;
    if (poly.degree(main_var) <= 0) return poly;
    auto coeffs = extract_coefficients(poly, main_var);
    if (coeffs.empty()) return MultiPoly();
    auto it = coeffs.begin();
    MultiPoly content_gcd = it->second;
    ++it;
    for (; it != coeffs.end(); ++it) {
        if (it->second.is_zero()) continue;
        content_gcd = multivariate_gcd(content_gcd, it->second);
        if (content_gcd.is_constant()) {
            Rational nc = content_gcd.numeric_content();
            if (nc == Rational(1)) return content_gcd;
        }
    }
    if (!content_gcd.is_zero() && !content_gcd.terms().empty()) {
        if (content_gcd.terms()[0].second < Rational(0))
            content_gcd = content_gcd * Rational(-1);
    }
    return content_gcd;
}
MultiPoly multivariate_primitive_part(const MultiPoly& poly, const std::string& main_var)
{
    if (poly.is_zero()) return poly;
    MultiPoly content = multivariate_content(poly, main_var);
    if (content.is_zero()) return MultiPoly();
    if (content.is_constant() && content.numeric_content() == Rational(1)) return poly;

    /// content 的变量集是去掉 main_var 后的辅助变量集,
    /// 需要将其嵌入到 poly 的完整变量集中才能正确执行 exact_div.
    const auto& full_vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < full_vars.size(); ++i) {
        if (full_vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }

    MultiPoly content_in_full_vars = content;
    if (var_idx >= 0 && content.variables().size() < full_vars.size()) {
        std::vector<MultiPoly::Term> new_terms;
        for (const auto& term : content.terms()) {
            Monomial full_mono(full_vars.size(), 0);
            size_t ri = 0;
            for (size_t i = 0; i < full_vars.size(); ++i) {
                if (static_cast<int>(i) == var_idx) {
                    full_mono[i] = 0;
                } else {
                    if (ri < term.first.size()) full_mono[i] = term.first[ri];
                    ++ri;
                }
            }
            new_terms.emplace_back(std::move(full_mono), term.second);
        }
        content_in_full_vars = MultiPoly(std::move(new_terms), full_vars);
    }

    return poly.exact_div(content_in_full_vars);
}
MultiPoly multivariate_gcd(const MultiPoly& a, const MultiPoly& b)
{
    if (a.is_zero()) { return b.is_zero() ? MultiPoly() : b.make_primitive(); }
    if (b.is_zero()) return a.make_primitive();
    const auto& vars = a.variables().empty() ? b.variables() : a.variables();
    /// 两者均为常数
    if (a.is_constant() && b.is_constant()) {
        Rational g = rational_gcd(a.numeric_content(), b.numeric_content());
        return MultiPoly(g, vars);
    }
    /// 基本情形:两者均为一元多项式
    if (a.is_univariate() && b.is_univariate()) {
        std::string a_var, b_var;
        for (const auto& v : a.variables()) { if (a.degree(v) > 0) { a_var = v; break; } }
        for (const auto& v : b.variables()) { if (b.degree(v) > 0) { b_var = v; break; } }
        if (a_var.empty() && b_var.empty()) {
            return MultiPoly(rational_gcd(a.numeric_content(), b.numeric_content()), vars);
        }
        if (a_var.empty() || b_var.empty()) {
            return MultiPoly(rational_gcd(a.numeric_content(), b.numeric_content()), vars);
        }
        if (a_var == b_var) {
            Rational nc_gcd = rational_gcd(a.numeric_content(), b.numeric_content());
            Polynomial<Rational> pa = a.to_univariate();
            Polynomial<Rational> pb = b.to_univariate();
            Polynomial<Rational> g = Polynomial<Rational>::gcd(pa, pb);
            if (g.is_zero()) return MultiPoly(nc_gcd, vars);
            /// monic GCD * numeric content GCD
            MultiPoly result = MultiPoly::from_univariate(g, a_var);
            result = result * nc_gcd;
            if (vars.size() > 1) {
                int avi = -1;
                for (size_t i = 0; i < vars.size(); ++i) {
                    if (vars[i] == a_var) { avi = static_cast<int>(i); break; }
                }
                std::vector<MultiPoly::Term> nt;
                for (const auto& term : result.terms()) {
                    Monomial fm(vars.size(), 0);
                    if (avi >= 0 && !term.first.empty()) fm[static_cast<size_t>(avi)] = term.first[0];
                    nt.emplace_back(std::move(fm), term.second);
                }
                result = MultiPoly(std::move(nt), vars);
            }
            if (!result.terms().empty() && result.terms()[0].second < Rational(0))
                result = result * Rational(-1);
            return result;
        }
        return MultiPoly(rational_gcd(a.numeric_content(), b.numeric_content()), vars);
    }
    /// 递归情形:求值-插值法
    std::string main_var = choose_main_variable(a, b);
    Rational nc_a = a.numeric_content();
    Rational nc_b = b.numeric_content();
    Rational nc_gcd = rational_gcd(nc_a, nc_b);
    MultiPoly a_prim = a.make_primitive();
    MultiPoly b_prim = b.make_primitive();
    int deg_a = a_prim.degree(main_var);
    int deg_b = b_prim.degree(main_var);
    int degree_bound = std::min(deg_a, deg_b);
    int num_points = degree_bound + 1;
    std::vector<Rational> eval_pts;
    eval_pts.push_back(Rational(0));
    for (int k = 1; static_cast<int>(eval_pts.size()) < num_points + 10; ++k) {
        eval_pts.push_back(Rational(k));
        eval_pts.push_back(Rational(-k));
    }
    std::vector<Rational> used_pts;
    std::vector<MultiPoly> gcd_vals;
    int target_deg = -1;
    for (size_t pi = 0; pi < eval_pts.size() &&
         static_cast<int>(used_pts.size()) < num_points; ++pi) {
        Rational pt = eval_pts[pi];
        MultiPoly lc_a = a_prim.leading_coeff(main_var);
        MultiPoly lc_b = b_prim.leading_coeff(main_var);
        MultiPoly lca_ev = lc_a.eval(main_var, pt);
        MultiPoly lcb_ev = lc_b.eval(main_var, pt);
        if (lca_ev.is_zero() || lcb_ev.is_zero()) continue;
        MultiPoly a_ev = a_prim.eval(main_var, pt);
        MultiPoly b_ev = b_prim.eval(main_var, pt);
        MultiPoly g_ev = multivariate_gcd(a_ev, b_ev);
        int gd = g_ev.is_zero() ? -1 : g_ev.total_degree();
        if (target_deg < 0) { target_deg = gd; }
        else if (gd < target_deg) { used_pts.clear(); gcd_vals.clear(); target_deg = gd; }
        else if (gd > target_deg) { continue; }
        if (!g_ev.is_zero()) g_ev = g_ev.make_primitive();
        used_pts.push_back(pt);
        gcd_vals.push_back(g_ev);
    }
    if (target_deg <= 0 || used_pts.empty()) return MultiPoly(nc_gcd, vars);
    if (used_pts.size() == 1 && degree_bound == 0) {
        MultiPoly r = gcd_vals[0].is_zero() ? MultiPoly(Rational(1), vars) : gcd_vals[0];
        if (nc_gcd != Rational(1)) r = r * nc_gcd;
        return r;
    }
    MultiPoly cand = lagrange_interpolate(used_pts, gcd_vals, main_var, vars);
    if (cand.is_zero()) return MultiPoly(nc_gcd, vars);
    cand = cand.make_primitive();
    if (divides_poly(a_prim, cand) && divides_poly(b_prim, cand)) {
        if (nc_gcd != Rational(1)) cand = cand * nc_gcd;
        if (!cand.terms().empty() && cand.terms()[0].second < Rational(0))
            cand = cand * Rational(-1);
        return cand;
    }
    /// 扩展搜索
    int ext_bound = degree_bound + 5;
    std::vector<Rational> ext_pts;
    std::vector<MultiPoly> ext_vals;
    std::vector<Rational> all_pts;
    all_pts.push_back(Rational(0));
    for (int k = 1; static_cast<int>(all_pts.size()) < ext_bound + 10; ++k) {
        all_pts.push_back(Rational(k));
        all_pts.push_back(Rational(-k));
    }
    target_deg = -1;
    for (size_t pi = 0; pi < all_pts.size() &&
         static_cast<int>(ext_pts.size()) < ext_bound + 1; ++pi) {
        Rational pt = all_pts[pi];
        MultiPoly lc_a = a_prim.leading_coeff(main_var);
        MultiPoly lc_b = b_prim.leading_coeff(main_var);
        MultiPoly lca_ev = lc_a.eval(main_var, pt);
        MultiPoly lcb_ev = lc_b.eval(main_var, pt);
        if (lca_ev.is_zero() || lcb_ev.is_zero()) continue;
        MultiPoly a_ev = a_prim.eval(main_var, pt);
        MultiPoly b_ev = b_prim.eval(main_var, pt);
        MultiPoly g_ev = multivariate_gcd(a_ev, b_ev);
        int gd = g_ev.is_zero() ? -1 : g_ev.total_degree();
        if (target_deg < 0) { target_deg = gd; }
        else if (gd < target_deg) { ext_pts.clear(); ext_vals.clear(); target_deg = gd; }
        else if (gd > target_deg) { continue; }
        if (!g_ev.is_zero()) g_ev = g_ev.make_primitive();
        ext_pts.push_back(pt);
        ext_vals.push_back(g_ev);
    }
    if (target_deg <= 0 || ext_pts.empty()) return MultiPoly(nc_gcd, vars);
    cand = lagrange_interpolate(ext_pts, ext_vals, main_var, vars);
    if (cand.is_zero()) return MultiPoly(nc_gcd, vars);
    cand = cand.make_primitive();
    if (divides_poly(a_prim, cand) && divides_poly(b_prim, cand)) {
        if (nc_gcd != Rational(1)) cand = cand * nc_gcd;
        if (!cand.terms().empty() && cand.terms()[0].second < Rational(0))
            cand = cand * Rational(-1);
        return cand;
    }
    return MultiPoly(nc_gcd, vars);
}
/**
 * @internal
 * @brief 计算多元多项式关于指定变量的形式偏导数
 *
 * 对每个含 main_var 的项,新系数 = 原系数 x 指数,新指数 = 原指数 - 1.
 * 指数为 0 的项被丢弃.
 *
 * @param[in] poly     输入多元多项式
 * @param[in] main_var 求导变量名
 * @return partialpoly/partialmain_var
 */
static MultiPoly formal_derivative(const MultiPoly& poly, const std::string& main_var)
{
    if (poly.is_zero()) return poly;

    const auto& vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    if (var_idx < 0) {
        /// main_var 不在变量列表中,导数为零
        return MultiPoly(Rational(0), vars);
    }

    std::vector<MultiPoly::Term> result_terms;
    for (const auto& term : poly.terms()) {
        const Monomial& mono = term.first;
        int exp = (static_cast<size_t>(var_idx) < mono.size()) ? mono[var_idx] : 0;
        if (exp == 0) continue; // 常数项关于该变量导数为零

        Rational new_coeff = term.second * Rational(exp);
        Monomial new_mono = mono;
        new_mono[var_idx] = exp - 1;
        result_terms.emplace_back(std::move(new_mono), std::move(new_coeff));
    }

    if (result_terms.empty()) return MultiPoly(Rational(0), vars);
    return MultiPoly(std::move(result_terms), vars);
}

/**
 * @brief 多元无平方因子分解(Yun 算法的多元推广)
 *
 * 通过计算 gcd(f, partialf/partialx_main) 检测重因子,将多项式分解为
 * f = f_1 * f_2^2 * f₃^3 * ... 的形式,其中每个 fᵢ 无平方.
 *
 * @param[in] poly     输入多元多项式
 * @param[in] main_var 主变量名
 * @return 无平方因子分解结果,components[i] 的重数为 i+1
 */
SquareFreeDecomp square_free_decompose(const MultiPoly& poly, const std::string& main_var)
{
    /// 零多项式或常数多项式:直接返回
    if (poly.is_zero() || poly.is_constant()) {
        return SquareFreeDecomp{{poly}};
    }

    /// 若 main_var 不在变量列表中或次数为 0,多项式无重因子
    const auto& vars = poly.variables();
    int var_idx = -1;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == main_var) { var_idx = static_cast<int>(i); break; }
    }
    if (var_idx < 0 || poly.degree(main_var) <= 0) {
        return SquareFreeDecomp{{poly}};
    }

    /// 步骤 1:计算 f' = partialf/partialmain_var
    MultiPoly f_prime = formal_derivative(poly, main_var);

    /// 步骤 2:计算 g = gcd(f, f');特征 0 下 f' 为零时,f 直接作为单个无平方因子.
    if (f_prime.is_zero()) {
        return SquareFreeDecomp{{poly}};
    }

    MultiPoly g = multivariate_gcd(poly, f_prime);

    /// 步骤 3:若 g 为常数,f 已无平方
    if (g.is_constant()) {
        return SquareFreeDecomp{{poly}};
    }

    /// 步骤 4:Yun 算法迭代提取无平方分量
    /// w = f / g, y = f' / g, z = y - w'
    MultiPoly w = poly.exact_div(g);
    MultiPoly y = f_prime.exact_div(g);
    MultiPoly w_prime = formal_derivative(w, main_var);
    MultiPoly z = y - w_prime;

    std::vector<MultiPoly> components;

    while (!z.is_zero()) {
        MultiPoly component_i = multivariate_gcd(w, z);
        w = w.exact_div(component_i);
        y = z.exact_div(component_i);
        w_prime = formal_derivative(w, main_var);
        z = y - w_prime;
        components.push_back(std::move(component_i));
    }

    /// w 是最后一个分量
    components.push_back(std::move(w));

    /// 过滤掉常数为 1 的平凡分量(保留结构但不影响语义)
    /// 注意:components[i] 的重数为 i+1
    return SquareFreeDecomp{std::move(components)};
}
/**
 * @internal
 * @brief 首项系数预计算(Wang's leading coefficient trick)
 *
 * 在 Hensel 提升前,将原多项式关于主变量的首项系数 lc(f, x_main) 分解为
 * 辅助变量的因子,并按求值后的首项系数值分配给各一元因子,
 * 使提升过程保持正确的首项系数与本原内容.
 *
 * 算法步骤:
 * 1. 计算 lc = lc(f, x_main),为辅助变量的多项式
 * 2. 若 lc 为常数,无需预计算
 * 3. 将 lc 在求值点处求值得到 lc_eval(有理数)
 * 4. 对每个一元因子 fᵢ,其首项系数 lc_i 应整除 lc_eval
 * 5. 按 lc_i 在 lc_eval 中的贡献比例,将 lc 的因子分配给各 fᵢ
 * 6. 将分配的 lc 部分乘入各 fᵢ,使提升后首项系数正确
 *
 * @param[in]     poly              原始多元多项式
 * @param[in]     main_var          主变量名
 * @param[in]     eval_points       辅助变量到求值点的映射
 * @param[in,out] univariate_factors 一元因子列表,函数会修改各因子的首项系数
 *
 * @see Wang, P.S. "An Improved Multivariate Polynomial Factoring Algorithm."
 *      Mathematics of Computation, 32(144), 1978.
 * @see D. Y. Y. Yun, "On Square-Free Decomposition Algorithms,"
 *      Proceedings of SYMSAC 1976.
 */
static void precompute_leading_coefficients(
    const MultiPoly& poly,
    const std::string& main_var,
    const std::map<std::string, Rational>& eval_points,
    std::vector<Polynomial<Rational>>& univariate_factors)
{
    int r = static_cast<int>(univariate_factors.size());
    if (r <= 1) return;

    /// 步骤 1:计算 lc(f, x_main) - 关于主变量的首项系数多项式
    MultiPoly lc_poly = poly.leading_coeff(main_var);

    /// 步骤 2:若 lc 为常数,无需预计算
    if (lc_poly.is_constant()) return;

    /// 步骤 3:将 lc 在求值点处求值得到 lc_eval(有理数值)
    MultiPoly lc_evaluated = lc_poly.eval(eval_points);
    if (lc_evaluated.is_zero()) return; // 退化情形:首项系数消失

    Rational lc_eval_val = lc_evaluated.is_zero() ? Rational(0)
                         : lc_evaluated.terms()[0].second;

    /// 步骤 4:收集各一元因子的首项系数
    /// 各因子 fᵢ 的首项系数 lc_i 满足 prod lc_i = lc_eval_val(至多差常数)
    std::vector<Rational> factor_lcs;
    factor_lcs.reserve(r);
    Rational product_of_lcs(1);
    for (int i = 0; i < r; ++i) {
        Rational lc_i = univariate_factors[i].lead_coeff();
        if (lc_i.is_zero()) lc_i = Rational(1);
        factor_lcs.push_back(lc_i);
        product_of_lcs = product_of_lcs * lc_i;
    }

    /// 若一元因子首项系数之积已等于 lc_eval_val,无需调整
    if (product_of_lcs == lc_eval_val) return;

    /// 步骤 5:计算分配比例
    /// 策略:将 lc_eval_val / product_of_lcs 的差额分配给第一个因子
    /// 这是简化版本--对于大多数情形(lc 为单项式或简单多项式),
    /// 将整个 lc 分配给首项系数最大的因子即可保证提升正确性.
    //
    /// 完整版本需要递归分解 lc_poly 并逐一匹配,但对于首次实现,
    /// 采用按比例分配的简化策略.

    if (product_of_lcs.is_zero()) return;

    /// 计算缩放因子:scale = lc_eval_val / product_of_lcs
    Rational scale = lc_eval_val / product_of_lcs;

    if (scale == Rational(1)) return;

    /// 策略:将缩放因子分配给第一个因子
    /// 这保证 prod lc(fᵢ) = lc_eval_val
    /// 对于更复杂的情形(lc_poly 有多个不可约因子),
    /// 需要匹配各因子的首项系数与 lc_poly 的因子.
    //
    /// 高级分配:尝试将 scale 分解为各因子的贡献
    /// 对每个因子,检查 lc_eval_val 是否能被其首项系数整除
    /// 若能,则该因子获得对应的 lc 份额

    /// 尝试精确分配:对每个因子 fᵢ,计算其应得的 lc 份额
    /// 方法:lc_eval_val = c_1 * c_2 * ... * cᵣ,其中 cᵢ 是分配给 fᵢ 的值
    /// 约束:cᵢ 的求值值等于 fᵢ 的首项系数乘以某个因子

    /// 简化实现:将整个缩放因子乘入第一个因子
    Polynomial<Rational>& first_factor = univariate_factors[0];
    std::vector<Rational> new_coeffs = first_factor.coeffs;
    for (auto& c : new_coeffs) {
        c = c * scale;
    }
    first_factor = Polynomial<Rational>(new_coeffs, first_factor.variable_name);
}

} // namespace lamina
