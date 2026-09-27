#include "internal/complex_root_isolation.hpp"
#include "internal/exact_sturm.hpp"

#include <algorithm>
#include <utility>

namespace LMCAS::detail::complex_root {
int rational_sign(const Rational& value) {
    if (value < Rational(0)) return -1;
    if (value > Rational(0)) return 1;
    return 0;
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

namespace {

struct ProjectionGrid {
    CoordinateProjections projections;
    std::vector<RationalInterval> real_roots;
    std::vector<RationalInterval> imaginary_roots;
    std::vector<RationalInterval> real_bands;
    std::vector<RationalInterval> imaginary_bands;
};

Result<void> refine_projection_interval(
    const Poly& projection, RationalInterval& interval,
    std::size_t index, std::size_t round,
    ComputationContext& context, const std::string& operation) {
    if (interval.first == interval.second) return Result<void>::success();
    if (round < 8) return refine_interval(projection, interval, context, operation);
    const double width = (interval.second - interval.first).to_double();
    ExactRealAlgebraic value{projection, interval.first, interval.second, index, 1};
    auto refined = refine_exact_real_algebraic_to_tolerance(
        value, width / 16.0, 1e-15, context, operation);
    if (!refined) return refined;
    interval = {value.lower, value.upper};
    return Result<void>::success();
}

Result<void> separate_projection_intervals(
    const Poly& projection, std::vector<RationalInterval>& intervals,
    ComputationContext& context, const std::string& operation) {
    if (intervals.size() < 2) return Result<void>::success();
    bool separated = false;
    std::size_t round = 0;
    while (!separated) {
        separated = true;
        for (std::size_t index = 1; index < intervals.size(); ++index) {
            if (intervals[index - 1].second < intervals[index].first) continue;
            separated = false;
            auto left = refine_projection_interval(
                projection, intervals[index - 1], index - 1, round, context, operation);
            if (!left) return left;
            auto right = refine_projection_interval(
                projection, intervals[index], index, round, context, operation);
            if (!right) return right;
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


Result<ProjectionGrid> prepare_projection_grid(
    const Poly& polynomial, ComputationContext& context,
    const std::string& operation) {
    auto projections = coordinate_projections(polynomial, context, operation);
    if (!projections) return Result<ProjectionGrid>::failure(projections.error());
    auto real = isolate_real_roots_exact(
        projections.value().real, context, operation + ".projection.real");
    if (!real) return Result<ProjectionGrid>::failure(real.error());
    auto imaginary = isolate_real_roots_exact(
        projections.value().imaginary, context, operation + ".projection.imaginary");
    if (!imaginary) return Result<ProjectionGrid>::failure(imaginary.error());
    auto real_separated = separate_projection_intervals(
        projections.value().real, real.value(), context, operation + ".projection.real");
    if (!real_separated) return Result<ProjectionGrid>::failure(real_separated.error());
    auto imaginary_separated = separate_projection_intervals(
        projections.value().imaginary, imaginary.value(), context,
        operation + ".projection.imaginary");
    if (!imaginary_separated) {
        return Result<ProjectionGrid>::failure(imaginary_separated.error());
    }
    return Result<ProjectionGrid>::success(ProjectionGrid{
        std::move(projections.value()), std::move(real.value()),
        std::move(imaginary.value()), {}, {}});
}

Result<void> separate_imaginary_axis(
    const Poly& projection, RationalInterval& interval,
    ComputationContext& context, const std::string& operation) {
    while (interval.first <= Rational(0) && interval.second >= Rational(0) &&
           !(interval.first == Rational(0) && interval.second == Rational(0))) {
        if (projection.eval(Rational(0)) == Rational(0)) {
            interval = {Rational(0), Rational(0)};
            break;
        }
        auto refined = refine_interval(projection, interval, context, operation);
        if (!refined) return refined;
    }
    return Result<void>::success();
}

Result<void> certify_conjugate_projections(
    ProjectionGrid& grid, std::size_t complex_count,
    ComputationContext& context, const std::string& operation) {
    bool have_positive = false;
    bool have_negative = false;
    for (auto& interval : grid.imaginary_roots) {
        auto separated = separate_imaginary_axis(
            grid.projections.imaginary, interval, context, operation);
        if (!separated) return separated;
        if (interval.first > Rational(0)) have_positive = true;
        if (interval.second < Rational(0)) have_negative = true;
    }
    if (!have_positive || !have_negative || (complex_count & 1U) != 0) {
        return Result<void>::failure(
            CasErrc::InternalInvariant,
            "non-real roots lacked conjugate imaginary projections", operation);
    }
    return Result<void>::success();
}

ComplexIsolation make_cell_isolation(
    const ProjectionGrid& grid, const RationalRectangle& rectangle,
    std::size_t real_index, std::size_t imaginary_index) {
    const auto& real_interval = grid.real_roots[real_index];
    const auto& imaginary_interval = grid.imaginary_roots[imaginary_index];
    return ComplexIsolation{
        rectangle.real_lower, rectangle.real_upper,
        rectangle.imaginary_lower, rectangle.imaginary_upper,
        ExactRealAlgebraic{grid.projections.real,
            real_interval.first, real_interval.second, real_index, 1},
        ExactRealAlgebraic{grid.projections.imaginary,
            imaginary_interval.first, imaginary_interval.second, imaginary_index, 1}};
}

Result<void> isolate_grid_row(
    const Poly& polynomial, const ProjectionGrid& grid, std::size_t real_index,
    std::vector<ComplexIsolation>& roots, ComputationContext& context,
    const std::string& operation) {
    for (std::size_t imaginary_index = 0;
         imaginary_index < grid.imaginary_bands.size(); ++imaginary_index) {
        const auto& interval = grid.imaginary_roots[imaginary_index];
        if (grid.projections.imaginary.eval(Rational(0)) == Rational(0) &&
            interval.first <= Rational(0) && interval.second >= Rational(0)) continue;
        const RationalRectangle rectangle{
            grid.real_bands[real_index].first, grid.real_bands[real_index].second,
            grid.imaginary_bands[imaginary_index].first,
            grid.imaginary_bands[imaginary_index].second};
        auto count = count_rectangle_roots(
            polynomial, rectangle, context, operation + ".rectangles.projection_grid");
        if (!count) return Result<void>::failure(count.error());
        if (!count.value().boundary_clear || count.value().count > 1) {
            return Result<void>::failure(
                CasErrc::InternalInvariant,
                "projection-band rectangle was not a clear one-root cell", operation);
        }
        if (count.value().count == 0) continue;
        roots.push_back(make_cell_isolation(grid, rectangle, real_index, imaginary_index));
    }
    return Result<void>::success();
}

}

Result<std::vector<ComplexIsolation>> isolate_nonreal_roots(
    const Poly& polynomial, std::size_t complex_count,
    ComputationContext& context, const std::string& operation) {
    auto prepared = prepare_projection_grid(polynomial, context, operation);
    if (!prepared) return Result<std::vector<ComplexIsolation>>::failure(prepared.error());
    auto& grid = prepared.value();
    auto conjugate = certify_conjugate_projections(grid, complex_count, context, operation);
    if (!conjugate) return Result<std::vector<ComplexIsolation>>::failure(conjugate.error());
    const Rational bound = strict_cauchy_root_bound(polynomial);
    grid.real_bands = projection_bands(grid.real_roots, bound);
    grid.imaginary_bands = projection_bands(grid.imaginary_roots, bound);
    std::vector<ComplexIsolation> roots;
    roots.reserve(complex_count);
    for (std::size_t real_index = 0; real_index < grid.real_bands.size(); ++real_index) {
        auto isolated = isolate_grid_row(polynomial, grid, real_index, roots, context, operation);
        if (!isolated) return Result<std::vector<ComplexIsolation>>::failure(isolated.error());
    }
    if (roots.size() != complex_count) {
        return Result<std::vector<ComplexIsolation>>::failure(
            CasErrc::InternalInvariant,
            "projection-band rectangles did not cover every non-real root", operation);
    }
    std::sort(roots.begin(), roots.end(),
        [](const ComplexIsolation& left, const ComplexIsolation& right) {
            if (left.real_projection.root_index != right.real_projection.root_index) {
                return left.real_projection.root_index < right.real_projection.root_index;
            }
            return left.imaginary_projection.root_index < right.imaginary_projection.root_index;
        });
    return Result<std::vector<ComplexIsolation>>::success(std::move(roots));
}

}

namespace LMCAS::detail {

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

}
