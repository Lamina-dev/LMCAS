#include "internal/multivariate_factor_checked_support.hpp"
#include "internal/multivariate_factor_support.hpp"
#include <algorithm>
#include <map>
#include <new>
#include <optional>
#include <utility>

namespace LMCAS {
namespace multivariate_checked_detail {

static Result<void> append_primitive_factor(MultiFactorResult& result,
                                           const MultiPoly& raw, int multiplicity) {
    if (multiplicity <= 0) {
        return Result<void>::failure(CasErrc::InternalInvariant,
            "因子重数必须为正", "factor_multivariate");
    }
    if (raw.is_zero() || raw.is_constant()) return Result<void>::success();
    MultiPoly factor = raw.make_primitive();
    if (!factor.terms().empty() && factor.terms()[0].second < Rational(0)) {
        factor = factor * Rational(-1);
    }
    auto existing = std::find(result.factors.begin(), result.factors.end(), factor);
    if (existing == result.factors.end()) {
        result.factors.push_back(std::move(factor));
        result.multiplicities.push_back(multiplicity);
    } else {
        const auto index = static_cast<size_t>(existing - result.factors.begin());
        result.multiplicities[index] += multiplicity;
    }
    return Result<void>::success();
}

MultiFactorCheckedResult assemble_checked_factorization(
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
        const int multiplicity = i < raw_multiplicities.size() ? raw_multiplicities[i] : 1;
        auto appended = append_primitive_factor(result, raw_factors[i], multiplicity);
        if (!appended) return MultiFactorCheckedResult::failure(appended.error());
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

}

namespace {
using namespace multivariate_checked_detail;

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

static MultiFactorCheckedResult factor_multivariate_impl(
    const MultiPoly& poly, ComputationContext& context);

static MultiFactorCheckedResult factor_common_monomial(
    const MultiPoly& poly, const Monomial& common_monomial,
    const MultiPoly& quotient, ComputationContext& context) {
    const auto& variables = poly.variables();
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

static MultiFactorCheckedResult factor_difference_of_squares(
    const MultiPoly& poly, const std::pair<MultiPoly, MultiPoly>& difference,
    ComputationContext& context) {
    std::vector<MultiPoly> factors;
    std::vector<int> multiplicities;
    Completeness completeness = Completeness::Complete;
    std::string reason;
    for (MultiPoly candidate :
         {difference.first + difference.second,
          difference.first - difference.second}) {
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

static MultiFactorCheckedResult factor_binomial_power(
    const MultiPoly& poly, const MultiPoly& primitive,
    const std::string& left_variable, const std::string& right_variable,
    ComputationContext& context) {
    const auto& variables = poly.variables();
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

static MultiFactorCheckedResult factor_embedded_univariate(
    const MultiPoly& poly, const MultiPoly& primitive,
    const std::string& main_variable, ComputationContext& context) {
    const auto& variables = poly.variables();
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
    const bool has_common_monomial = std::any_of(
        common_monomial.begin(), common_monomial.end(),
        [](int exponent) { return exponent > 0; });
    if (has_common_monomial) {
        return factor_common_monomial(poly, common_monomial, quotient, context);
    }
    MultiPoly primitive = poly.make_primitive();
    if (auto difference = detail::detect_difference_of_squares(primitive)) {
        return factor_difference_of_squares(poly, *difference, context);
    }
    if (auto binomial = detail::detect_binomial_power(primitive)) {
        const auto& [left, right, exponent] = *binomial;
        (void)exponent;
        return factor_binomial_power(poly, primitive, left, right, context);
    }
    if (primitive.is_univariate()) {
        return factor_embedded_univariate(poly, primitive, main_variable, context);
    }
    auto evaluated = factor_evaluated_linear(poly, primitive, main_variable, context);
    if (evaluated) return std::move(*evaluated);
    auto homogeneous =
        detail::factor_homogeneous_bivariate(primitive, context);
    if (!homogeneous) {
        return MultiFactorCheckedResult::failure(homogeneous.error());
    }
    if (homogeneous.value()) {
        auto& result = *homogeneous.value();
        return assemble_checked_factorization(
            poly, result.factors, result.multiplicities,
            Completeness::Complete, {});
    }

    return assemble_checked_factorization(
        poly, {poly}, {1}, Completeness::Inconclusive,
        "当前受检算法无法证明该多项式不可约");
}

}

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
        return MultiFactorCheckedResult::failure(
            CasErrc::InternalInvariant, error.what(),
            "factor_multivariate");
    }
}

MultiFactorCheckedResult factor_multivariate_checked(const MultiPoly& poly)
{
    ComputationContext context;
    return factor_multivariate_checked(poly, context);
}

}
