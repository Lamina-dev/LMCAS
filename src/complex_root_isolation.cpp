#include "internal/complex_root_isolation.hpp"
#include "internal/multivariate_factor_support.hpp"

#include <exception>
#include <iterator>
#include <optional>
#include <utility>

namespace LMCAS::detail {
namespace {

using complex_root::Poly;

std::optional<Rational> exact_rational_square_root(
    const Rational& value) {
    if (value < Rational(0)) return std::nullopt;
    const BigInt numerator = value.get_numerator();
    const BigInt denominator = value.get_denominator();
    const BigInt numerator_root = numerator.sqrt();
    const BigInt denominator_root = denominator.sqrt();
    if (numerator_root * numerator_root != numerator ||
        denominator_root * denominator_root != denominator) {
        return std::nullopt;
    }
    return Rational(numerator_root, denominator_root);
}

std::vector<Poly> factor_biquadratic(const Poly& polynomial) {
    if (polynomial.degree() != 4 ||
        polynomial.coeffs[1] != Rational(0) ||
        polynomial.coeffs[3] != Rational(0)) {
        return {};
    }
    const Rational& c = polynomial.coeffs[0];
    const Rational& b = polynomial.coeffs[2];
    const Rational& a = polynomial.coeffs[4];
    const Rational discriminant = b * b - Rational(4) * a * c;
    auto square_root = exact_rational_square_root(discriminant);
    if (!square_root) return {};
    const Rational first =
        (-b - *square_root) / (Rational(2) * a);
    const Rational second =
        (-b + *square_root) / (Rational(2) * a);
    if (first == second) return {};
    return {
        Poly({-first, Rational(0), Rational(1)}, polynomial.variable_name),
        Poly({-second, Rational(0), Rational(1)}, polynomial.variable_name)};
}

Result<int> compare_isolations(
    const RootIsolation& left,
    const RootIsolation& right,
    ComputationContext& context) {
    const auto* left_real = std::get_if<RealIsolation>(&left);
    const auto* right_real = std::get_if<RealIsolation>(&right);
    if (left_real && !right_real) {
        return Result<int>::success(-1);
    }
    if (!left_real && right_real) {
        return Result<int>::success(1);
    }
    if (left_real && right_real) {
        return compare_exact_real_algebraic(
            left_real->value, right_real->value, context);
    }
    const auto& left_complex = std::get<ComplexIsolation>(left);
    const auto& right_complex = std::get<ComplexIsolation>(right);
    auto real_comparison = compare_exact_real_algebraic(
        left_complex.real_projection,
        right_complex.real_projection,
        context);
    if (!real_comparison || real_comparison.value() != 0) {
        return real_comparison;
    }
    return compare_exact_real_algebraic(
        left_complex.imaginary_projection,
        right_complex.imaginary_projection,
        context);
}

struct RootIsolationCacheEntry {
    Poly polynomial;
    std::vector<RootIsolation> roots;
    std::size_t step_cost = 1;
};

thread_local std::vector<RootIsolationCacheEntry> root_isolation_cache;

bool same_polynomial(const Poly& left, const Poly& right) {
    return left.coeffs == right.coeffs;
}

Result<void> order_isolations(
    std::vector<RootIsolation>& roots, ComputationContext& context) {
    for (std::size_t index = 1; index < roots.size(); ++index) {
        RootIsolation current = std::move(roots[index]);
        std::size_t position = index;
        while (position > 0) {
            auto comparison = compare_isolations(roots[position - 1], current, context);
            if (!comparison) return Result<void>::failure(comparison.error());
            if (comparison.value() <= 0) break;
            roots[position] = std::move(roots[position - 1]);
            --position;
        }
        roots[position] = std::move(current);
    }
    return Result<void>::success();
}

Result<std::vector<Poly>> strict_factors(
    const Poly& polynomial, ComputationContext& context) {
    auto factors = factor_biquadratic(polynomial);
    if (factors.empty() && polynomial.degree() >= 6) {
        auto factor_result = LMCAS::factor_univariate_bridge_checked(polynomial, context);
        if (!factor_result) return Result<std::vector<Poly>>::failure(factor_result.error());
        factors = std::move(factor_result.value().value);
    }
    if (factors.empty()) factors = {polynomial};
    std::size_t factored_degree = 0;
    bool strict_factorization = factors.size() > 1;
    for (auto& factor : factors) {
        factor = factor.square_free_part().make_monic();
        factor.variable_name = "_root";
        if (factor.degree() <= 0 || factor.degree() >= polynomial.degree()) {
            strict_factorization = false;
            break;
        }
        factored_degree += static_cast<std::size_t>(factor.degree());
    }
    if (!strict_factorization ||
        factored_degree != static_cast<std::size_t>(polynomial.degree())) factors.clear();
    return Result<std::vector<Poly>>::success(std::move(factors));
}

Result<std::vector<RootIsolation>> isolate_factors(
    const std::vector<Poly>& factors, std::size_t degree,
    ComputationContext& context, const std::string& operation) {
    std::vector<RootIsolation> combined;
    combined.reserve(degree);
    for (const auto& factor : factors) {
        auto factor_roots = isolate_exact_roots(factor, context, operation + ".factor");
        if (!factor_roots) return Result<std::vector<RootIsolation>>::failure(factor_roots.error());
        combined.insert(combined.end(),
            std::make_move_iterator(factor_roots.value().begin()),
            std::make_move_iterator(factor_roots.value().end()));
    }
    auto ordered = order_isolations(combined, context);
    if (!ordered) return Result<std::vector<RootIsolation>>::failure(ordered.error());
    return Result<std::vector<RootIsolation>>::success(std::move(combined));
}

Result<std::vector<RootIsolation>> isolate_unfactored(
    const Poly& polynomial, ComputationContext& context,
    const std::string& operation) {
    auto real_roots = isolate_real_roots_exact(polynomial, context, operation + ".real");
    if (!real_roots) return Result<std::vector<RootIsolation>>::failure(real_roots.error());
    const std::size_t real_count = real_roots.value().size();
    const std::size_t degree = static_cast<std::size_t>(polynomial.degree());
    if (real_count > degree) {
        return Result<std::vector<RootIsolation>>::failure(
            CasErrc::InternalInvariant, "real-root count exceeded polynomial degree", operation);
    }
    std::vector<RootIsolation> result;
    result.reserve(degree);
    for (std::size_t index = 0; index < real_count; ++index) {
        result.push_back(RealIsolation{ExactRealAlgebraic{
            polynomial, real_roots.value()[index].first,
            real_roots.value()[index].second, index, 1}});
    }
    if (real_count == degree) {
        return Result<std::vector<RootIsolation>>::success(std::move(result));
    }
    const Poly certificate_polynomial = complex_root::clear_denominators(polynomial);
    auto complex_roots = complex_root::isolate_nonreal_roots(
        certificate_polynomial, degree - real_count, context, operation);
    if (!complex_roots) return Result<std::vector<RootIsolation>>::failure(complex_roots.error());
    for (auto& isolation : complex_roots.value()) result.push_back(std::move(isolation));
    if (result.size() != degree) {
        return Result<std::vector<RootIsolation>>::failure(
            CasErrc::InternalInvariant,
            "certified root isolation did not cover the polynomial degree", operation);
    }
    return Result<std::vector<RootIsolation>>::success(std::move(result));
}

Result<std::vector<RootIsolation>> remember_isolations(
    const Poly& polynomial, std::vector<RootIsolation> roots, std::size_t step_cost) {
    if (root_isolation_cache.size() >= 16) {
        root_isolation_cache.erase(root_isolation_cache.begin());
    }
    root_isolation_cache.push_back(RootIsolationCacheEntry{polynomial, roots, step_cost});
    return Result<std::vector<RootIsolation>>::success(std::move(roots));
}

} // namespace

Result<std::vector<RootIsolation>> isolate_exact_roots(
    const Polynomial<Rational>& input, ComputationContext& context,
    const std::string& operation) {
    const std::size_t starting_steps = context.steps_used();
    try {
        auto step = context.consume_steps(1, operation);
        if (!step) return Result<std::vector<RootIsolation>>::failure(step.error());
        if (input.degree() <= 0) {
            return Result<std::vector<RootIsolation>>::failure(
                CasErrc::InvalidArgument,
                "exact root isolation requires a nonconstant polynomial", operation);
        }
        Poly polynomial = input.square_free_part().make_monic();
        polynomial.variable_name = "_root";
        for (const auto& entry : root_isolation_cache) {
            if (!same_polynomial(entry.polynomial, polynomial)) continue;
            const std::size_t remaining_cost = entry.step_cost > 0 ? entry.step_cost - 1 : 0;
            auto cached_step = context.consume_steps(remaining_cost, operation + ".cache");
            if (!cached_step) return Result<std::vector<RootIsolation>>::failure(cached_step.error());
            return Result<std::vector<RootIsolation>>::success(entry.roots);
        }
        auto factors = strict_factors(polynomial, context);
        if (!factors) return Result<std::vector<RootIsolation>>::failure(factors.error());
        auto roots = factors.value().empty()
            ? isolate_unfactored(polynomial, context, operation)
            : isolate_factors(factors.value(), static_cast<std::size_t>(polynomial.degree()),
                              context, operation);
        if (!roots) return roots;
        return remember_isolations(polynomial, std::move(roots.value()),
                                   context.steps_used() - starting_steps);
    } catch (const std::bad_alloc&) {
        return Result<std::vector<RootIsolation>>::failure(
            CasErrc::ResourceLimit, "complex root isolation allocation failed", operation);
    } catch (const std::exception& error) {
        return Result<std::vector<RootIsolation>>::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

Result<RootIsolation> isolate_exact_root(
    const ExactRootId& root,
    ComputationContext& context,
    const std::string& operation) {
    auto roots = isolate_exact_roots(root.polynomial, context, operation);
    if (!roots) return Result<RootIsolation>::failure(roots.error());
    if (root.index >= roots.value().size()) {
        return Result<RootIsolation>::failure(
            CasErrc::InvalidArgument,
            "RootOf index exceeds the number of distinct roots",
            operation);
    }
    return Result<RootIsolation>::success(
        std::move(roots.value()[root.index]));
}

} // namespace LMCAS::detail
