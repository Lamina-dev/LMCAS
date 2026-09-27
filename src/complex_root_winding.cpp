#include "internal/complex_root_isolation.hpp"
#include "rational_polynomial.hpp"

#include <utility>

namespace LMCAS::detail::complex_root {
namespace {

std::pair<Poly, Poly> edge_image(
    const Poly& polynomial,
    const Rational& real_start,
    const Rational& imaginary_start,
    const Rational& real_delta,
    const Rational& imaginary_delta) {
    const Poly real_line(
        std::vector<Rational>{real_start, real_delta}, "t");
    const Poly imaginary_line(
        std::vector<Rational>{imaginary_start, imaginary_delta}, "t");
    Poly real("t");
    Poly imaginary("t");
    for (std::size_t position = polynomial.coeffs.size(); position-- > 0;) {
        Poly next_real = real * real_line - imaginary * imaginary_line;
        Poly next_imaginary = real * imaginary_line + imaginary * real_line;
        if (next_real.coeffs.empty()) next_real.coeffs.resize(1, Rational(0));
        next_real.coeffs[0] = next_real.coeffs[0] + polynomial.coeffs[position];
        next_real.trim();
        real = std::move(next_real);
        imaginary = std::move(next_imaginary);
    }
    return {std::move(real), std::move(imaginary)};
}

std::vector<std::pair<Poly, Poly>> rectangle_edges(
    const Poly& polynomial,
    const RationalRectangle& rectangle) {
    const Rational width = rectangle.real_upper - rectangle.real_lower;
    const Rational height =
        rectangle.imaginary_upper - rectangle.imaginary_lower;
    std::vector<std::pair<Poly, Poly>> edges;
    edges.reserve(4);
    edges.push_back(edge_image(
        polynomial, rectangle.real_lower, rectangle.imaginary_lower,
        width, Rational(0)));
    edges.push_back(edge_image(
        polynomial, rectangle.real_upper, rectangle.imaginary_lower,
        Rational(0), height));
    edges.push_back(edge_image(
        polynomial, rectangle.real_upper, rectangle.imaginary_upper,
        -width, Rational(0)));
    edges.push_back(edge_image(
        polynomial, rectangle.real_lower, rectangle.imaginary_upper,
        Rational(0), -height));
    return edges;
}

Result<bool> edge_is_clear(
    const Poly& real,
    const Poly& imaginary,
    ComputationContext& context,
    const std::string& operation) {
    if (real.eval(Rational(0)) == Rational(0) &&
        imaginary.eval(Rational(0)) == Rational(0)) {
        return Result<bool>::success(false);
    }
    if (real.eval(Rational(1)) == Rational(0) &&
        imaginary.eval(Rational(1)) == Rational(0)) {
        return Result<bool>::success(false);
    }
    Poly common = Poly::gcd(real, imaginary);
    if (common.degree() <= 0) return Result<bool>::success(true);
    auto roots = count_real_roots_exact(
        common, Rational(0), Rational(1), context, operation);
    if (!roots) return Result<bool>::failure(roots.error());
    return Result<bool>::success(roots.value() == 0);
}

Result<bool> pole_interval_is_clear(
    const Poly& denominator, const Poly& numerator,
    const RationalInterval& interval, ComputationContext& context,
    const std::string& operation) {
    auto denominator_roots = count_real_roots_exact(
        denominator, interval.first, interval.second, context, operation);
    if (!denominator_roots) return Result<bool>::failure(denominator_roots.error());
    auto numerator_roots = count_real_roots_exact(
        numerator, interval.first, interval.second, context, operation);
    if (!numerator_roots) return Result<bool>::failure(numerator_roots.error());
    return Result<bool>::success(denominator_roots.value() == 1 &&
        numerator_roots.value() == 0 &&
        denominator.eval(interval.first) != Rational(0));
}

Result<void> refine_pole_interval(
    const Poly& factor, const Poly& denominator, const Poly& numerator,
    RationalInterval& interval, ComputationContext& context,
    const std::string& operation) {
    while (true) {
        auto step = context.consume_steps(1, operation);
        if (!step) {
            return step;
        }
        if (interval.first == interval.second) {
            break;
        }
        if (interval.second <= Rational(0) || interval.first >= Rational(1)) {
            break;
        }
        if (interval.first > Rational(0) && interval.second < Rational(1)) {
            auto clear = pole_interval_is_clear(
                denominator, numerator, interval, context, operation);
            if (!clear) {
                return Result<void>::failure(clear.error());
            }
            if (clear.value()) {
                break;
            }
        }
        auto refined = refine_interval(factor, interval, context, operation);
        if (!refined) {
            return refined;
        }
    }
    return Result<void>::success();
}

Result<Rational> left_pole_sample(
    const Poly& denominator, const Rational& root,
    ComputationContext& context, const std::string& operation) {
    Rational left_sample = root / Rational(2);
    while (true) {
        auto step = context.consume_steps(1, operation);
        if (!step) return Result<Rational>::failure(step.error());
        auto count = count_real_roots_exact(
            denominator, left_sample, root, context, operation);
        if (!count) return Result<Rational>::failure(count.error());
        if (count.value() == 1 && denominator.eval(left_sample) != Rational(0)) break;
        left_sample = (left_sample + root) / Rational(2);
    }
    return Result<Rational>::success(std::move(left_sample));
}

Result<int> pole_index(
    const Poly& denominator, const Poly& numerator,
    const RationalInterval& interval, ComputationContext& context,
    const std::string& operation) {
    if (interval.second <= Rational(0) || interval.first >= Rational(1)) {
        return Result<int>::success(0);
    }
    Rational numerator_sample;
    Rational left_sample;
    if (interval.first == interval.second) {
        const Rational root = interval.first;
        if (!(root > Rational(0) && root < Rational(1))) return Result<int>::success(0);
        numerator_sample = root;
        auto sample = left_pole_sample(denominator, root, context, operation);
        if (!sample) return Result<int>::failure(sample.error());
        left_sample = std::move(sample.value());
    } else {
        numerator_sample = (interval.first + interval.second) / Rational(2);
        left_sample = interval.first;
    }
    const int denominator_left_sign = rational_sign(denominator.eval(left_sample));
    const int numerator_sign = rational_sign(numerator.eval(numerator_sample));
    if (denominator_left_sign == 0 || numerator_sign == 0) {
        return Result<int>::failure(
            CasErrc::InternalInvariant,
            "Cauchy-index pole sign could not be certified", operation);
    }
    return Result<int>::success(denominator_left_sign * numerator_sign);
}

Result<int> cauchy_index_on_edge(
    const Poly& denominator, const Poly& numerator,
    ComputationContext& context, const std::string& operation) {
    int index = 0;
    const auto factors = square_free_factorization(denominator);
    for (const auto& [factor, multiplicity] : factors) {
        if ((multiplicity & 1) == 0) continue;
        auto roots = isolate_real_roots_exact(factor, context, operation);
        if (!roots) return Result<int>::failure(roots.error());
        for (auto interval : roots.value()) {
            auto refined = refine_pole_interval(
                factor, denominator, numerator, interval, context, operation);
            if (!refined) return Result<int>::failure(refined.error());
            auto contribution = pole_index(
                denominator, numerator, interval, context, operation);
            if (!contribution) return contribution;
            index += contribution.value();
        }
    }
    return Result<int>::success(index);
}

bool shear_is_valid(const std::vector<std::pair<Poly, Poly>>& edges,
                    const Rational& shear) {
    for (const auto& [real, imaginary] : edges) {
        const Poly denominator = real + scale_poly(imaginary, shear);
        if (denominator.is_zero() ||
            denominator.eval(Rational(0)) == Rational(0) ||
            denominator.eval(Rational(1)) == Rational(0)) return false;
    }
    return true;
}

Result<RectangleCount> count_sheared_edges(
    const std::vector<std::pair<Poly, Poly>>& edges, const Rational& shear,
    int degree, ComputationContext& context, const std::string& operation) {
    int total_index = 0;
    for (const auto& [real, imaginary] : edges) {
        const Poly denominator = real + scale_poly(imaginary, shear);
        auto edge_index = cauchy_index_on_edge(
            denominator, imaginary, context, operation);
        if (!edge_index) return Result<RectangleCount>::failure(edge_index.error());
        total_index += edge_index.value();
    }
    if ((total_index & 1) != 0) {
        return Result<RectangleCount>::failure(
            CasErrc::InternalInvariant,
            "argument-principle Cauchy index was odd", operation);
    }
    const int count = total_index / 2;
    if (count < 0 || count > degree) {
        return Result<RectangleCount>::failure(
            CasErrc::InternalInvariant,
            "argument-principle root count was outside polynomial degree", operation);
    }
    return Result<RectangleCount>::success(
        RectangleCount{true, static_cast<std::size_t>(count)});
}

}

Result<RectangleCount> count_rectangle_roots(
    const Poly& polynomial, const RationalRectangle& rectangle,
    ComputationContext& context, const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) return Result<RectangleCount>::failure(step.error());
    const auto edges = rectangle_edges(polynomial, rectangle);
    for (const auto& [real, imaginary] : edges) {
        auto clear = edge_is_clear(real, imaginary, context, operation);
        if (!clear) return Result<RectangleCount>::failure(clear.error());
        if (!clear.value()) {
            return Result<RectangleCount>::success(RectangleCount{false, 0});
        }
    }
    BigInt shear_magnitude(0);
    int shear_direction = 1;
    while (true) {
        step = context.consume_steps(1, operation);
        if (!step) return Result<RectangleCount>::failure(step.error());
        auto bits = context.require_integer_bits(shear_magnitude.bit_length(), operation);
        if (!bits) return Result<RectangleCount>::failure(bits.error());
        const Rational shear(shear_direction > 0 ? shear_magnitude : -shear_magnitude);
        if (shear_is_valid(edges, shear)) {
            return count_sheared_edges(edges, shear, polynomial.degree(), context, operation);
        }
        if (shear_direction > 0) {
            shear_direction = -1;
            if (shear_magnitude.is_zero()) shear_magnitude = BigInt(1);
        } else {
            shear_direction = 1;
            shear_magnitude += BigInt(1);
        }
    }
}

}
