#include "internal/multivariate_factor_support.hpp"
#include "exact_factorization.hpp"
#include <exception>
#include <iterator>
#include <new>
#include <utility>

namespace LMCAS {
Polynomial<BigInt> factor_integer_polynomial(const Polynomial<Rational>& poly) {
    BigInt denominator_lcm(1);
    for (const auto& coefficient : poly.coeffs) {
        denominator_lcm = BigInt::lcm(
            denominator_lcm, coefficient.get_denominator());
    }
    std::vector<BigInt> coefficients;
    coefficients.reserve(poly.coeffs.size());
    for (const auto& coefficient : poly.coeffs) {
        coefficients.push_back(coefficient.get_numerator() *
            (denominator_lcm / coefficient.get_denominator()));
    }
    return Polynomial<BigInt>(coefficients, poly.variable_name);
}
namespace {

constexpr const char* bridge_operation = "factor_univariate";

UnivariateFactorResult bridge_inconclusive(
    const Polynomial<Rational>& poly, const char* reason) {
    return UnivariateFactorResult::success(
        MathResult<std::vector<Polynomial<Rational>>>{
            {poly.make_monic()}, Completeness::Inconclusive, reason});
}

Result<BerlekampResult> bridge_reduction(
    const Polynomial<Rational>& poly, ComputationContext& context) {
    static const int64_t primes[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
    for (int64_t prime : primes) {
        auto step = context.consume_steps(1, bridge_operation);
        if (!step) return Result<BerlekampResult>::failure(step.error());
        auto reduction = berlekamp_factor(poly, prime);
        if (reduction.prime > 0 && !reduction.factors.empty()) {
            return Result<BerlekampResult>::success(std::move(reduction));
        }
    }
    return Result<BerlekampResult>::success(BerlekampResult{});
}


Result<int> bridge_lift_bound(const Polynomial<BigInt>& poly, int64_t prime,
                              BigInt& modulus, ComputationContext& context) {
    BigInt max_coefficient(0);
    for (const auto& coefficient : poly.coeffs) {
        if (coefficient.abs() > max_coefficient) max_coefficient = coefficient.abs();
    }
    BigInt target = max_coefficient * BigInt(2);
    for (int i = 0; i < poly.degree(); ++i) target = target * BigInt(2);
    int bound = 1;
    modulus = BigInt(static_cast<std::int64_t>(prime));
    while (modulus <= target && bound < 100) {
        auto step = context.consume_steps(1, bridge_operation);
        if (!step) return Result<int>::failure(step.error());
        modulus = modulus * BigInt(static_cast<std::int64_t>(prime));
        ++bound;
    }
    return Result<int>::success(bound);
}

bool bridge_inconclusive_error(const CasError& error) {
    return error.code == CasErrc::InvalidArgument || error.code == CasErrc::Inconclusive;
}

UnivariateFactorResult bridge_lift_and_combine(
    const Polynomial<Rational>& original, const Polynomial<Rational>& work,
    BerlekampResult& reduction, ComputationContext& context, bool& whole_input) {
    auto integer_poly = factor_integer_polynomial(work);
    BigInt modulus;
    auto bound = bridge_lift_bound(integer_poly, reduction.prime, modulus, context);
    if (!bound) return UnivariateFactorResult::failure(bound.error());
    const ModInt unit(
        (integer_poly.coeffs.back() % BigInt(reduction.prime)).to_int(),
        reduction.prime);
    for (auto& coefficient : reduction.factors.front().coeffs) {
        coefficient = coefficient * unit;
    }
    auto lifted = hensel_lift_checked(
        integer_poly, reduction.factors, reduction.prime, bound.value());
    if (!lifted) {
        if (bridge_inconclusive_error(lifted.error())) {
            whole_input = true;
            return bridge_inconclusive(original, "当前模素数的提升前置条件不成立");
        }
        return UnivariateFactorResult::failure(lifted.error());
    }
    auto combined = zassenhaus_combine_checked(work, lifted.value(), modulus, context);
    if (!combined && bridge_inconclusive_error(combined.error())) {
        whole_input = true;
        return bridge_inconclusive(original, "提升因子无法完成精确有理重构");
    }
    return combined;
}

Result<void> bridge_append_repeated(
    MathResult<std::vector<Polynomial<Rational>>>& factors,
    const TfSquareFreeResult& square_free, ComputationContext& context) {
    if (!square_free.had_repeated_factors || square_free.repeated_factor.degree() < 1) {
        return Result<void>::success();
    }
    auto repeated = factor_univariate_bridge_checked(square_free.repeated_factor, context);
    if (!repeated) return Result<void>::failure(repeated.error());
    if (repeated.value().completeness == Completeness::Inconclusive) {
        factors.completeness = Completeness::Inconclusive;
        if (factors.reason.empty()) factors.reason = repeated.value().reason;
    }
    auto& values = repeated.value().value;
    factors.value.insert(factors.value.end(), std::make_move_iterator(values.begin()),
                         std::make_move_iterator(values.end()));
    return Result<void>::success();
}

UnivariateFactorResult bridge_factor(const Polynomial<Rational>& poly,
                                      ComputationContext& context) {
    auto square_free = tf_square_free(poly);
    auto work = square_free.square_free.make_monic();
    auto reduction = bridge_reduction(work, context);
    if (!reduction) return UnivariateFactorResult::failure(reduction.error());
    if (reduction.value().factors.empty()) {
        return bridge_inconclusive(poly, "未找到保持次数且无平方的有效模素数");
    }
    MathResult<std::vector<Polynomial<Rational>>> factors{
        {}, Completeness::Complete, {}};
    if (reduction.value().factors.size() <= 1) {
        factors.value.push_back(work.make_monic());
    } else {
        bool whole_input = false;
        auto combined = bridge_lift_and_combine(
            poly, work, reduction.value(), context, whole_input);
        if (!combined) return combined;
        if (whole_input) return combined;
        factors = std::move(combined.value());
    }
    auto appended = bridge_append_repeated(factors, square_free, context);
    if (!appended) return UnivariateFactorResult::failure(appended.error());
    Polynomial<Rational> product({Rational(1)}, poly.variable_name);
    for (const auto& factor : factors.value) product = product * factor;
    if (!(product == poly.make_monic())) {
        return UnivariateFactorResult::failure(CasErrc::InternalInvariant,
            "一元因子未能精确重构输入多项式", bridge_operation);
    }
    return UnivariateFactorResult::success(std::move(factors));
}

}

UnivariateFactorResult factor_univariate_bridge_checked(
    const Polynomial<Rational>& poly, ComputationContext& context) {
    auto step = context.consume_steps(1, bridge_operation);
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
        return bridge_factor(poly, context);
    } catch (const std::bad_alloc&) {
        return UnivariateFactorResult::failure(
            CasErrc::ResourceLimit, "一元分解分配失败", bridge_operation);
    } catch (const std::exception& error) {
        return UnivariateFactorResult::failure(
            CasErrc::InternalInvariant, error.what(), bridge_operation);
    }
}

}
