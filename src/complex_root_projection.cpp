#include "internal/complex_root_isolation.hpp"

#include <algorithm>
#include <utility>

namespace LMCAS::detail::complex_root {

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

namespace {

using PolyMatrix = std::vector<std::vector<Poly>>;

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

Poly constant_poly(const Rational& value, const std::string& variable) {
    return Poly(std::vector<Rational>{value}, variable);
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

PolyMatrix sylvester_matrix(const std::vector<Poly>& left,
                            const std::vector<Poly>& right,
                            const std::string& variable) {
    const std::size_t left_degree = left.size() - 1;
    const std::size_t right_degree = right.size() - 1;
    const std::size_t size = left_degree + right_degree;
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
    return matrix;
}

Result<void> eliminate_resultant_row(
    PolyMatrix& matrix, std::size_t pivot, std::size_t row,
    const Poly& previous, ComputationContext& context,
    const std::string& operation) {
    for (std::size_t column = pivot + 1; column < matrix.size(); ++column) {
        Poly numerator = matrix[row][column] * matrix[pivot][pivot] -
            matrix[row][pivot] * matrix[pivot][column];
        if (pivot != 0) {
            auto quotient = exact_poly_divide(
                numerator, previous, context, operation);
            if (!quotient) return Result<void>::failure(quotient.error());
            matrix[row][column] = std::move(quotient.value());
        } else {
            matrix[row][column] = std::move(numerator);
        }
    }
    matrix[row][pivot] = Poly(matrix[pivot][pivot].variable_name);
    return Result<void>::success();
}

Result<Poly> resultant_determinant(
    PolyMatrix matrix, const std::string& variable,
    ComputationContext& context, const std::string& operation) {
    Poly previous = constant_poly(Rational(1), variable);
    int determinant_sign = 1;
    for (std::size_t pivot = 0; pivot + 1 < matrix.size(); ++pivot) {
        auto step = context.consume_steps(1, operation);
        if (!step) return Result<Poly>::failure(step.error());
        std::size_t replacement = pivot;
        while (replacement < matrix.size() && matrix[replacement][pivot].is_zero()) {
            ++replacement;
        }
        if (replacement == matrix.size()) return Result<Poly>::success(Poly(variable));
        if (replacement != pivot) {
            std::swap(matrix[pivot], matrix[replacement]);
            determinant_sign = -determinant_sign;
        }
        const Poly pivot_value = matrix[pivot][pivot];
        for (std::size_t row = pivot + 1; row < matrix.size(); ++row) {
            auto eliminated = eliminate_resultant_row(
                matrix, pivot, row, previous, context, operation);
            if (!eliminated) return Result<Poly>::failure(eliminated.error());
        }
        previous = pivot_value;
    }
    Poly determinant = std::move(matrix.back().back());
    if (determinant_sign < 0) {
        determinant = scale_poly(std::move(determinant), Rational(-1));
    }
    determinant.trim();
    return Result<Poly>::success(std::move(determinant));
}

bool zero_projection_polynomial(const std::vector<Poly>& coefficients) {
    if (coefficients.empty()) {
        return true;
    }
    return coefficients.size() == 1 && coefficients[0].is_zero();
}

Result<Poly> sylvester_resultant(
    std::vector<Poly> left, std::vector<Poly> right,
    const std::string& variable, ComputationContext& context,
    const std::string& operation) {
    while (left.size() > 1 && left.back().is_zero()) {
        left.pop_back();
    }
    while (right.size() > 1 && right.back().is_zero()) {
        right.pop_back();
    }
    if (zero_projection_polynomial(left) || zero_projection_polynomial(right)) {
        return Result<Poly>::success(Poly(variable));
    }
    const std::size_t size = left.size() + right.size() - 2;
    if (size == 0) {
        return Result<Poly>::success(constant_poly(Rational(1), variable));
    }
    auto step = context.consume_steps(size * size + 1, operation);
    if (!step) {
        return Result<Poly>::failure(step.error());
    }
    return resultant_determinant(
        sylvester_matrix(left, right, variable), variable, context, operation);
}

}

Result<CoordinateProjections> coordinate_projections(
    const Poly& polynomial, ComputationContext& context,
    const std::string& operation) {
    const ComplexBivariate parts = complex_bivariate_parts(polynomial);
    auto real = sylvester_resultant(
        coefficients_in_y(parts.real, "_re"),
        coefficients_in_y(parts.imaginary, "_re"),
        "_re", context, operation + ".resultant.real");
    if (!real) return Result<CoordinateProjections>::failure(real.error());
    auto imaginary = sylvester_resultant(
        coefficients_in_x(parts.real, "_im"),
        coefficients_in_x(parts.imaginary, "_im"),
        "_im", context, operation + ".resultant.imaginary");
    if (!imaginary) return Result<CoordinateProjections>::failure(imaginary.error());
    if (real.value().degree() <= 0 || imaginary.value().degree() <= 0) {
        return Result<CoordinateProjections>::failure(
            CasErrc::InternalInvariant,
            "complex coordinate resultant was identically zero", operation);
    }
    return Result<CoordinateProjections>::success(CoordinateProjections{
        real.value().square_free_part().make_monic(),
        imaginary.value().square_free_part().make_monic()});
}

}
