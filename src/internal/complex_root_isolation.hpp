#pragma once

#include "internal/exact_root.hpp"
#include "internal/exact_sturm.hpp"

namespace LMCAS::detail::complex_root {

using Poly = Polynomial<Rational>;

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

struct CoordinateProjections {
    Poly real;
    Poly imaginary;
};

int rational_sign(const Rational& value);
Poly scale_poly(Poly polynomial, const Rational& scalar);
Poly clear_denominators(Poly polynomial);
Result<CoordinateProjections> coordinate_projections(
    const Poly& polynomial, ComputationContext& context,
    const std::string& operation);
Result<void> refine_interval(
    const Poly& polynomial, RationalInterval& interval,
    ComputationContext& context, const std::string& operation);
Result<RectangleCount> count_rectangle_roots(
    const Poly& polynomial, const RationalRectangle& rectangle,
    ComputationContext& context, const std::string& operation);
Result<std::vector<ComplexIsolation>> isolate_nonreal_roots(
    const Poly& polynomial, std::size_t complex_count,
    ComputationContext& context, const std::string& operation);

}
