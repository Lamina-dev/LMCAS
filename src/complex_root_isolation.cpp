#include "internal/exact_root.hpp"

#include "internal/exact_sturm.hpp"
#include "internal/multivariate_factor_support.hpp"
#include "solve_polynomial.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <utility>
#include <optional>
#include <vector>

namespace lamina::detail {
namespace {

using Poly = Polynomial<Rational>;
using PolyMatrix = std::vector<std::vector<Poly>>;

struct RationalRectangle {
    Rational real_lower;
    Rational real_upper;
    Rational imaginary_lower;
    Rational imaginary_upper;
};

struct RectangleCount {
    bool boundary_clear = false;
    std::size_t count = 0;
};

struct BiPolynomial {
    std::vector<std::vector<Rational>> coefficients;

    void add(std::size_t x_degree, std::size_t y_degree,
             const Rational& coefficient) {
        if (coefficient == Rational(0)) return;
        if (coefficients.size() <= x_degree) {
            coefficients.resize(x_degree + 1);
        }
        if (coefficients[x_degree].size() <= y_degree) {
            coefficients[x_degree].resize(y_degree + 1, Rational(0));
        }
        coefficients[x_degree][y_degree] =
            coefficients[x_degree][y_degree] + coefficient;
    }
};

struct ComplexBivariate {
    BiPolynomial real;
    BiPolynomial imaginary;
};

int rational_sign(const Rational& value) {
    if (value < Rational(0)) return -1;
    if (value > Rational(0)) return 1;
    return 0;
}


Poly constant_poly(const Rational& value, const std::string& variable) {
    return Poly(std::vector<Rational>{value}, variable);
}

Poly scale_poly(Poly polynomial, const Rational& scalar) {
    for (auto& coefficient : polynomial.coeffs) {
        coefficient = coefficient * scalar;
    }
    polynomial.trim();
    return polynomial;
}

Poly clear_denominators(Poly polynomial) {
    BigInt common_denominator(1);
    for (const auto& coefficient : polynomial.coeffs) {
        common_denominator *= coefficient.get_denominator();
    }
    for (auto& coefficient : polynomial.coeffs) {
        coefficient = Rational(
            coefficient.get_numerator() *
                (common_denominator / coefficient.get_denominator()));
    }
    polynomial.trim();
    return polynomial;
}

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

Result<Poly> exact_poly_divide(
    const Poly& numerator,
    const Poly& denominator,
    ComputationContext& context,
    const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) return Result<Poly>::failure(step.error());
    if (denominator.is_zero()) {
        return Result<Poly>::failure(
            CasErrc::InternalInvariant,
            "fraction-free resultant encountered a zero divisor",
            operation);
    }
    auto division = numerator.div_mod(denominator);
    if (!division.second.is_zero()) {
        return Result<Poly>::failure(
            CasErrc::InternalInvariant,
            "fraction-free resultant division was not exact",
            operation);
    }
    return Result<Poly>::success(std::move(division.first));
}

BiPolynomial shift_x(const BiPolynomial& input) {
    BiPolynomial output;
    for (std::size_t x = 0; x < input.coefficients.size(); ++x) {
        for (std::size_t y = 0; y < input.coefficients[x].size(); ++y) {
            output.add(x + 1, y, input.coefficients[x][y]);
        }
    }
    return output;
}

BiPolynomial shift_y(const BiPolynomial& input) {
    BiPolynomial output;
    for (std::size_t x = 0; x < input.coefficients.size(); ++x) {
        for (std::size_t y = 0; y < input.coefficients[x].size(); ++y) {
            output.add(x, y + 1, input.coefficients[x][y]);
        }
    }
    return output;
}

void add_scaled(BiPolynomial& destination,
                const BiPolynomial& source,
                const Rational& scalar) {
    for (std::size_t x = 0; x < source.coefficients.size(); ++x) {
        for (std::size_t y = 0; y < source.coefficients[x].size(); ++y) {
            destination.add(x, y, source.coefficients[x][y] * scalar);
        }
    }
}

ComplexBivariate complex_bivariate_parts(const Poly& polynomial) {
    BiPolynomial power_real;
    power_real.add(0, 0, Rational(1));
    BiPolynomial power_imaginary;
    ComplexBivariate result;

    for (std::size_t degree = 0; degree < polynomial.coeffs.size(); ++degree) {
        add_scaled(result.real, power_real, polynomial.coeffs[degree]);
        add_scaled(result.imaginary, power_imaginary,
                   polynomial.coeffs[degree]);

        BiPolynomial next_real = shift_x(power_real);
        add_scaled(next_real, shift_y(power_imaginary), Rational(-1));
        BiPolynomial next_imaginary = shift_x(power_imaginary);
        add_scaled(next_imaginary, shift_y(power_real), Rational(1));
        power_real = std::move(next_real);
        power_imaginary = std::move(next_imaginary);
    }
    return result;
}

std::vector<Poly> coefficients_in_y(
    const BiPolynomial& polynomial,
    const std::string& remaining_variable) {
    std::size_t maximum_y = 0;
    for (const auto& x_coefficients : polynomial.coefficients) {
        if (!x_coefficients.empty()) {
            maximum_y = std::max(maximum_y, x_coefficients.size() - 1);
        }
    }
    std::vector<Poly> result(maximum_y + 1, Poly(remaining_variable));
    for (std::size_t x = 0; x < polynomial.coefficients.size(); ++x) {
        for (std::size_t y = 0; y < polynomial.coefficients[x].size(); ++y) {
            if (result[y].coeffs.size() <= x) {
                result[y].coeffs.resize(x + 1, Rational(0));
            }
            result[y].coeffs[x] = polynomial.coefficients[x][y];
        }
    }
    while (result.size() > 1 && result.back().is_zero()) result.pop_back();
    return result;
}

std::vector<Poly> coefficients_in_x(
    const BiPolynomial& polynomial,
    const std::string& remaining_variable) {
    std::vector<Poly> result(
        std::max<std::size_t>(1, polynomial.coefficients.size()),
        Poly(remaining_variable));
    for (std::size_t x = 0; x < polynomial.coefficients.size(); ++x) {
        if (result[x].coeffs.size() < polynomial.coefficients[x].size()) {
            result[x].coeffs.resize(
                polynomial.coefficients[x].size(), Rational(0));
        }
        for (std::size_t y = 0; y < polynomial.coefficients[x].size(); ++y) {
            result[x].coeffs[y] = polynomial.coefficients[x][y];
        }
        result[x].trim();
    }
    while (result.size() > 1 && result.back().is_zero()) result.pop_back();
    return result;
}

Result<Poly> sylvester_resultant(
    std::vector<Poly> left,
    std::vector<Poly> right,
    const std::string& variable,
    ComputationContext& context,
    const std::string& operation) {
    while (left.size() > 1 && left.back().is_zero()) left.pop_back();
    while (right.size() > 1 && right.back().is_zero()) right.pop_back();
    if (left.empty() || right.empty() ||
        (left.size() == 1 && left[0].is_zero()) ||
        (right.size() == 1 && right[0].is_zero())) {
        return Result<Poly>::success(Poly(variable));
    }

    const std::size_t left_degree = left.size() - 1;
    const std::size_t right_degree = right.size() - 1;
    const std::size_t size = left_degree + right_degree;
    if (size == 0) {
        return Result<Poly>::success(constant_poly(Rational(1), variable));
    }

    auto step = context.consume_steps(size * size + 1, operation);
    if (!step) return Result<Poly>::failure(step.error());
    PolyMatrix matrix(size, std::vector<Poly>(size, Poly(variable)));
    for (std::size_t row = 0; row < right_degree; ++row) {
        for (std::size_t column = 0; column <= left_degree; ++column) {
            matrix[row][row + column] = left[left_degree - column];
        }
    }
    for (std::size_t row = 0; row < left_degree; ++row) {
        for (std::size_t column = 0; column <= right_degree; ++column) {
            matrix[right_degree + row][row + column] =
                right[right_degree - column];
        }
    }

    Poly previous = constant_poly(Rational(1), variable);
    int determinant_sign = 1;
    for (std::size_t pivot = 0; pivot + 1 < size; ++pivot) {
        step = context.consume_steps(1, operation);
        if (!step) return Result<Poly>::failure(step.error());
        std::size_t replacement = pivot;
        while (replacement < size && matrix[replacement][pivot].is_zero()) {
            ++replacement;
        }
        if (replacement == size) return Result<Poly>::success(Poly(variable));
        if (replacement != pivot) {
            std::swap(matrix[pivot], matrix[replacement]);
            determinant_sign = -determinant_sign;
        }
        const Poly pivot_value = matrix[pivot][pivot];
        for (std::size_t row = pivot + 1; row < size; ++row) {
            for (std::size_t column = pivot + 1; column < size; ++column) {
                Poly numerator = matrix[row][column] * pivot_value -
                    matrix[row][pivot] * matrix[pivot][column];
                if (pivot != 0) {
                    auto quotient = exact_poly_divide(
                        numerator, previous, context, operation);
                    if (!quotient) return quotient;
                    matrix[row][column] = std::move(quotient.value());
                } else {
                    matrix[row][column] = std::move(numerator);
                }
            }
            matrix[row][pivot] = Poly(variable);
        }
        previous = pivot_value;
    }

    Poly determinant = matrix[size - 1][size - 1];
    if (determinant_sign < 0) {
        determinant = scale_poly(std::move(determinant), Rational(-1));
    }
    determinant.trim();
    return Result<Poly>::success(std::move(determinant));
}

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

Result<void> refine_interval(
    const Poly& polynomial,
    RationalInterval& interval,
    ComputationContext& context,
    const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (interval.first == interval.second) return Result<void>::success();
    const Rational midpoint =
        (interval.first + interval.second) / Rational(2);
    const Rational middle_value = polynomial.eval(midpoint);
    if (middle_value == Rational(0)) {
        interval = {midpoint, midpoint};
        return Result<void>::success();
    }
    const int lower_sign = rational_sign(polynomial.eval(interval.first));
    if (lower_sign == 0) {
        interval.second = interval.first;
    } else if (lower_sign != rational_sign(middle_value)) {
        interval.second = midpoint;
    } else {
        interval.first = midpoint;
    }
    return Result<void>::success();
}

Result<int> cauchy_index_on_edge(
    const Poly& denominator,
    const Poly& numerator,
    ComputationContext& context,
    const std::string& operation) {
    int index = 0;
    const auto factors = square_free_factorization(denominator);
    for (const auto& [factor, multiplicity] : factors) {
        if ((multiplicity & 1) == 0) continue;
        auto roots = isolate_real_roots_exact(factor, context, operation);
        if (!roots) return Result<int>::failure(roots.error());
        for (auto interval : roots.value()) {
            while (true) {
                auto step = context.consume_steps(1, operation);
                if (!step) return Result<int>::failure(step.error());
                if (interval.first == interval.second) break;
                if (interval.second <= Rational(0) ||
                    interval.first >= Rational(1)) {
                    break;
                }
                if (interval.first > Rational(0) &&
                    interval.second < Rational(1)) {
                    auto denominator_roots = count_real_roots_exact(
                        denominator, interval.first, interval.second,
                        context, operation);
                    if (!denominator_roots) {
                        return Result<int>::failure(denominator_roots.error());
                    }
                    auto numerator_roots = count_real_roots_exact(
                        numerator, interval.first, interval.second,
                        context, operation);
                    if (!numerator_roots) {
                        return Result<int>::failure(numerator_roots.error());
                    }
                    if (denominator_roots.value() == 1 &&
                        numerator_roots.value() == 0 &&
                        denominator.eval(interval.first) != Rational(0)) {
                        break;
                    }
                }
                auto refined = refine_interval(
                    factor, interval, context, operation);
                if (!refined) return Result<int>::failure(refined.error());
            }

            if (interval.second <= Rational(0) ||
                interval.first >= Rational(1)) {
                continue;
            }

            Rational numerator_sample;
            Rational left_sample;
            if (interval.first == interval.second) {
                const Rational root = interval.first;
                if (!(root > Rational(0) && root < Rational(1))) continue;
                numerator_sample = root;
                left_sample = root / Rational(2);
                while (true) {
                    auto step = context.consume_steps(1, operation);
                    if (!step) return Result<int>::failure(step.error());
                    auto count = count_real_roots_exact(
                        denominator, left_sample, root, context, operation);
                    if (!count) return Result<int>::failure(count.error());
                    if (count.value() == 1 &&
                        denominator.eval(left_sample) != Rational(0)) {
                        break;
                    }
                    left_sample = (left_sample + root) / Rational(2);
                }
            } else {
                numerator_sample =
                    (interval.first + interval.second) / Rational(2);
                left_sample = interval.first;
            }

            const int denominator_left_sign =
                rational_sign(denominator.eval(left_sample));
            const int numerator_sign =
                rational_sign(numerator.eval(numerator_sample));
            if (denominator_left_sign == 0 || numerator_sign == 0) {
                return Result<int>::failure(
                    CasErrc::InternalInvariant,
                    "Cauchy-index pole sign could not be certified",
                    operation);
            }
            index += denominator_left_sign * numerator_sign;
        }
    }
    return Result<int>::success(index);
}

Result<RectangleCount> count_rectangle_roots(
    const Poly& polynomial,
    const RationalRectangle& rectangle,
    ComputationContext& context,
    const std::string& operation) {
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

    long long shear_magnitude = 0;
    int shear_direction = 1;
    while (true) {
        step = context.consume_steps(1, operation);
        if (!step) return Result<RectangleCount>::failure(step.error());
        const long long shear_value = shear_magnitude * shear_direction;
        const Rational shear(shear_value);
        bool valid = true;
        for (const auto& [real, imaginary] : edges) {
            const Poly denominator = real + scale_poly(imaginary, shear);
            if (denominator.is_zero() ||
                denominator.eval(Rational(0)) == Rational(0) ||
                denominator.eval(Rational(1)) == Rational(0)) {
                valid = false;
                break;
            }
        }
        if (valid) {
            int total_index = 0;
            for (const auto& [real, imaginary] : edges) {
                const Poly denominator = real + scale_poly(imaginary, shear);
                auto edge_index = cauchy_index_on_edge(
                    denominator, imaginary, context, operation);
                if (!edge_index) {
                    return Result<RectangleCount>::failure(edge_index.error());
                }
                total_index += edge_index.value();
            }
            if ((total_index & 1) != 0) {
                return Result<RectangleCount>::failure(
                    CasErrc::InternalInvariant,
                    "argument-principle Cauchy index was odd",
                    operation);
            }
            const int count = total_index / 2;
            if (count < 0 || count > polynomial.degree()) {
                return Result<RectangleCount>::failure(
                    CasErrc::InternalInvariant,
                    "argument-principle root count was outside polynomial degree",
                    operation);
            }
            return Result<RectangleCount>::success(
                RectangleCount{true, static_cast<std::size_t>(count)});
        }
        if (shear_direction > 0) {
            shear_direction = -1;
            if (shear_magnitude == 0) shear_magnitude = 1;
        } else {
            shear_direction = 1;
            ++shear_magnitude;
        }
    }
}

Result<void> separate_projection_intervals(
    const Poly& projection,
    std::vector<RationalInterval>& intervals,
    ComputationContext& context,
    const std::string& operation) {
    if (intervals.size() < 2) return Result<void>::success();
    bool separated = false;
    std::size_t round = 0;
    while (!separated) {
        separated = true;
        for (std::size_t index = 1; index < intervals.size(); ++index) {
            if (intervals[index - 1].second < intervals[index].first) {
                continue;
            }
            separated = false;
            for (const std::size_t target : {index - 1, index}) {
                if (intervals[target].first == intervals[target].second) {
                    continue;
                }
                if (round < 8) {
                    auto refined = refine_interval(
                        projection, intervals[target], context, operation);
                    if (!refined) return refined;
                } else {
                    const double width =
                        (intervals[target].second -
                         intervals[target].first).to_double();
                    ExactRealAlgebraic value{
                        projection,
                        intervals[target].first,
                        intervals[target].second,
                        target,
                        1};
                    auto refined =
                        refine_exact_real_algebraic_to_tolerance(
                            value, width / 16.0, 1e-15,
                            context, operation);
                    if (!refined) return refined;
                    intervals[target] = {value.lower, value.upper};
                }
            }
        }
        ++round;
    }
    return Result<void>::success();
}

std::vector<RationalInterval> projection_bands(
    const std::vector<RationalInterval>& intervals,
    const Rational& bound) {
    std::vector<RationalInterval> bands;
    bands.reserve(intervals.size());
    for (std::size_t index = 0; index < intervals.size(); ++index) {
        const Rational lower = index == 0
            ? -bound
            : (intervals[index - 1].second + intervals[index].first) /
                  Rational(2);
        const Rational upper = index + 1 == intervals.size()
            ? bound
            : (intervals[index].second + intervals[index + 1].first) /
                  Rational(2);
        bands.emplace_back(lower, upper);
    }
    return bands;
}


Rational strict_root_bound(const Poly& polynomial) {
    Rational maximum(0);
    const Rational leading = polynomial.lead_coeff().abs();
    for (int degree = 0; degree < polynomial.degree(); ++degree) {
        const Rational ratio = polynomial.coeffs[degree].abs() / leading;
        if (maximum < ratio) maximum = ratio;
    }
    return Rational(2) + maximum;
}

Result<int> compare_isolations(
    const RootIsolation& left,
    const RootIsolation& right,
    ComputationContext& context) {
    const auto* left_real = std::get_if<RealIsolation>(&left);
    const auto* right_real = std::get_if<RealIsolation>(&right);
    if (left_real && !right_real) return Result<int>::success(-1);
    if (!left_real && right_real) return Result<int>::success(1);
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

} // namespace

Result<std::vector<RootIsolation>> isolate_exact_roots(
    const Polynomial<Rational>& input,
    ComputationContext& context,
    const std::string& operation) {
    const std::size_t starting_steps = context.steps_used();
    try {
        auto step = context.consume_steps(1, operation);
        if (!step) return Result<std::vector<RootIsolation>>::failure(step.error());
        if (input.degree() <= 0) {
            return Result<std::vector<RootIsolation>>::failure(
                CasErrc::InvalidArgument,
                "exact root isolation requires a nonconstant polynomial",
                operation);
        }
        Poly polynomial = input.square_free_part().make_monic();
        polynomial.variable_name = "_root";
        for (const auto& entry : root_isolation_cache) {
            if (!same_polynomial(entry.polynomial, polynomial)) continue;
            const std::size_t remaining_cost =
                entry.step_cost > 0 ? entry.step_cost - 1 : 0;
            auto cached_step = context.consume_steps(
                remaining_cost, operation + ".cache");
            if (!cached_step) {
                return Result<std::vector<RootIsolation>>::failure(
                    cached_step.error());
            }
            return Result<std::vector<RootIsolation>>::success(entry.roots);
        }
        auto remember = [&](std::vector<RootIsolation> roots) {
            if (root_isolation_cache.size() >= 16) {
                root_isolation_cache.erase(root_isolation_cache.begin());
            }
            root_isolation_cache.push_back(RootIsolationCacheEntry{
                polynomial, roots, context.steps_used() - starting_steps});
            return Result<std::vector<RootIsolation>>::success(
                std::move(roots));
        };

        auto factors = factor_biquadratic(polynomial);
        if (factors.empty()) {
            factors = polynomial.degree() >= 6
                ? lamina::factor_univariate_bridge(polynomial)
                : std::vector<Poly>{polynomial};
        }
        std::size_t factored_degree = 0;
        bool strict_factorization = factors.size() > 1;
        for (auto& factor : factors) {
            factor = factor.square_free_part().make_monic();
            factor.variable_name = "_root";
            if (factor.degree() <= 0 ||
                factor.degree() >= polynomial.degree()) {
                strict_factorization = false;
                break;
            }
            factored_degree += static_cast<std::size_t>(factor.degree());
        }
        if (strict_factorization &&
            factored_degree ==
                static_cast<std::size_t>(polynomial.degree())) {
            std::vector<RootIsolation> combined;
            combined.reserve(factored_degree);
            for (const auto& factor : factors) {
                auto factor_roots = isolate_exact_roots(
                    factor, context, operation + ".factor");
                if (!factor_roots) {
                    return Result<std::vector<RootIsolation>>::failure(
                        factor_roots.error());
                }
                combined.insert(
                    combined.end(),
                    std::make_move_iterator(factor_roots.value().begin()),
                    std::make_move_iterator(factor_roots.value().end()));
            }
            for (std::size_t index = 1; index < combined.size(); ++index) {
                RootIsolation current = std::move(combined[index]);
                std::size_t position = index;
                while (position > 0) {
                    auto comparison = compare_isolations(
                        combined[position - 1], current, context);
                    if (!comparison) {
                        return Result<std::vector<RootIsolation>>::failure(
                            comparison.error());
                    }
                    if (comparison.value() <= 0) break;
                    combined[position] = std::move(combined[position - 1]);
                    --position;
                }
                combined[position] = std::move(current);
            }
            return remember(std::move(combined));
        }

        auto real_roots = isolate_real_roots_exact(
            polynomial, context, operation + ".real");
        if (!real_roots) {
            return Result<std::vector<RootIsolation>>::failure(
                real_roots.error());
        }
        const std::size_t real_count = real_roots.value().size();
        const std::size_t degree = static_cast<std::size_t>(polynomial.degree());
        if (real_count > degree) {
            return Result<std::vector<RootIsolation>>::failure(
                CasErrc::InternalInvariant,
                "real-root count exceeded polynomial degree",
                operation);
        }

        std::vector<RootIsolation> result;
        result.reserve(degree);
        for (std::size_t index = 0; index < real_count; ++index) {
            result.push_back(RealIsolation{ExactRealAlgebraic{
                polynomial, real_roots.value()[index].first,
                real_roots.value()[index].second, index, 1}});
        }
        if (real_count == degree) {
            return remember(std::move(result));
        }

        const Poly certificate_polynomial = clear_denominators(polynomial);
        const ComplexBivariate parts =
            complex_bivariate_parts(certificate_polynomial);
        auto real_resultant = sylvester_resultant(
            coefficients_in_y(parts.real, "_re"),
            coefficients_in_y(parts.imaginary, "_re"),
            "_re", context, operation + ".resultant.real");
        if (!real_resultant) {
            return Result<std::vector<RootIsolation>>::failure(
                real_resultant.error());
        }
        auto imaginary_resultant = sylvester_resultant(
            coefficients_in_x(parts.real, "_im"),
            coefficients_in_x(parts.imaginary, "_im"),
            "_im", context, operation + ".resultant.imaginary");
        if (!imaginary_resultant) {
            return Result<std::vector<RootIsolation>>::failure(
                imaginary_resultant.error());
        }
        if (real_resultant.value().degree() <= 0 ||
            imaginary_resultant.value().degree() <= 0) {
            return Result<std::vector<RootIsolation>>::failure(
                CasErrc::InternalInvariant,
                "complex coordinate resultant was identically zero",
                operation);
        }
        Poly real_projection =
            real_resultant.value().square_free_part().make_monic();
        Poly imaginary_projection =
            imaginary_resultant.value().square_free_part().make_monic();
        auto real_projection_roots = isolate_real_roots_exact(
            real_projection, context, operation + ".projection.real");
        if (!real_projection_roots) {
            return Result<std::vector<RootIsolation>>::failure(
                real_projection_roots.error());
        }
        auto imaginary_projection_roots = isolate_real_roots_exact(
            imaginary_projection, context, operation + ".projection.imaginary");
        if (!imaginary_projection_roots) {
            return Result<std::vector<RootIsolation>>::failure(
                imaginary_projection_roots.error());
        }
        auto real_separated = separate_projection_intervals(
            real_projection, real_projection_roots.value(), context,
            operation + ".projection.real");
        if (!real_separated) {
            return Result<std::vector<RootIsolation>>::failure(
                real_separated.error());
        }
        auto imaginary_separated = separate_projection_intervals(
            imaginary_projection, imaginary_projection_roots.value(),
            context, operation + ".projection.imaginary");
        if (!imaginary_separated) {
            return Result<std::vector<RootIsolation>>::failure(
                imaginary_separated.error());
        }

        bool have_positive = false;
        bool have_negative = false;
        Rational positive_lower;
        Rational negative_upper;
        for (auto& interval : imaginary_projection_roots.value()) {
            while (interval.first <= Rational(0) &&
                   interval.second >= Rational(0) &&
                   !(interval.first == Rational(0) &&
                     interval.second == Rational(0))) {
                if (imaginary_projection.eval(Rational(0)) == Rational(0)) {
                    interval = {Rational(0), Rational(0)};
                    break;
                }
                auto refined = refine_interval(
                    imaginary_projection, interval, context, operation);
                if (!refined) {
                    return Result<std::vector<RootIsolation>>::failure(
                        refined.error());
                }
            }
            if (interval.first > Rational(0) &&
                (!have_positive || interval.first < positive_lower)) {
                positive_lower = interval.first;
                have_positive = true;
            }
            if (interval.second < Rational(0) &&
                (!have_negative || negative_upper < interval.second)) {
                negative_upper = interval.second;
                have_negative = true;
            }
        }

        const std::size_t complex_count = degree - real_count;
        if (!have_positive || !have_negative || (complex_count & 1U) != 0) {
            return Result<std::vector<RootIsolation>>::failure(
                CasErrc::InternalInvariant,
                "non-real roots lacked conjugate imaginary projections",
                operation);
        }
        const Rational bound = strict_root_bound(certificate_polynomial);
        const auto real_bands = projection_bands(
            real_projection_roots.value(), bound);
        const auto imaginary_bands = projection_bands(
            imaginary_projection_roots.value(), bound);
        std::vector<ComplexIsolation> complex_roots;
        complex_roots.reserve(complex_count);
        for (std::size_t real_index = 0;
             real_index < real_bands.size(); ++real_index) {
            for (std::size_t imaginary_index = 0;
                 imaginary_index < imaginary_bands.size(); ++imaginary_index) {
                const auto& imaginary_interval =
                    imaginary_projection_roots.value()[imaginary_index];
                if (imaginary_projection.eval(Rational(0)) == Rational(0) &&
                    imaginary_interval.first <= Rational(0) &&
                    imaginary_interval.second >= Rational(0)) {
                    continue;
                }
                const RationalRectangle rectangle{
                    real_bands[real_index].first,
                    real_bands[real_index].second,
                    imaginary_bands[imaginary_index].first,
                    imaginary_bands[imaginary_index].second};
                auto count = count_rectangle_roots(
                    certificate_polynomial, rectangle, context,
                    operation + ".rectangles.projection_grid");
                if (!count) {
                    return Result<std::vector<RootIsolation>>::failure(
                        count.error());
                }
                if (!count.value().boundary_clear ||
                    count.value().count > 1) {
                    return Result<std::vector<RootIsolation>>::failure(
                        CasErrc::InternalInvariant,
                        "projection-band rectangle was not a clear one-root cell",
                        operation);
                }
                if (count.value().count == 0) continue;
                complex_roots.push_back(ComplexIsolation{
                    rectangle.real_lower,
                    rectangle.real_upper,
                    rectangle.imaginary_lower,
                    rectangle.imaginary_upper,
                    ExactRealAlgebraic{
                        real_projection,
                        real_projection_roots.value()[real_index].first,
                        real_projection_roots.value()[real_index].second,
                        real_index, 1},
                    ExactRealAlgebraic{
                        imaginary_projection,
                        imaginary_interval.first,
                        imaginary_interval.second,
                        imaginary_index, 1}});
            }
        }
        if (complex_roots.size() != complex_count) {
            return Result<std::vector<RootIsolation>>::failure(
                CasErrc::InternalInvariant,
                "projection-band rectangles did not cover every non-real root",
                operation);
        }
        std::sort(complex_roots.begin(), complex_roots.end(),
                  [](const ComplexIsolation& left,
                     const ComplexIsolation& right) {
                      if (left.real_projection.root_index !=
                          right.real_projection.root_index) {
                          return left.real_projection.root_index <
                              right.real_projection.root_index;
                      }
                      return left.imaginary_projection.root_index <
                          right.imaginary_projection.root_index;
                  });
        for (auto& isolation : complex_roots) {
            result.push_back(std::move(isolation));
        }
        if (result.size() != degree) {
            return Result<std::vector<RootIsolation>>::failure(
                CasErrc::InternalInvariant,
                "certified root isolation did not cover the polynomial degree",
                operation);
        }
        return remember(std::move(result));
    } catch (const std::bad_alloc&) {
        return Result<std::vector<RootIsolation>>::failure(
            CasErrc::ResourceLimit,
            "complex root isolation allocation failed",
            operation);
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

Result<void> refine_complex_isolation(
    const Polynomial<Rational>&,
    ComplexIsolation& isolation,
    const NumericEvaluationOptions& options,
    ComputationContext& context,
    const std::string& operation) {
    auto real = refine_exact_real_algebraic_to_tolerance(
        isolation.real_projection,
        options.absolute_tolerance,
        options.relative_tolerance,
        context,
        operation);
    if (!real) return real;
    auto imaginary = refine_exact_real_algebraic_to_tolerance(
        isolation.imaginary_projection,
        options.absolute_tolerance,
        options.relative_tolerance,
        context,
        operation);
    if (!imaginary) return imaginary;
    isolation.real_lower = isolation.real_projection.lower;
    isolation.real_upper = isolation.real_projection.upper;
    isolation.imaginary_lower = isolation.imaginary_projection.lower;
    isolation.imaginary_upper = isolation.imaginary_projection.upper;
    return Result<void>::success();
}

} // namespace lamina::detail
