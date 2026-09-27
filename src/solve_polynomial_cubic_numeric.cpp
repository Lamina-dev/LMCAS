#include "internal/polynomial_solver_support.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace LMCAS::polynomial_solver_detail {

static std::optional<NumericCubicDepression> depress_numeric_cubic(
    double a, double b, double c, double d) {
    double p_val = 0.0;
    double q_val = 0.0;
    double shift_val = 0.0;
    bool numeric_depression = false;
    const double max_coefficient = std::max({
        std::abs(a), std::abs(b), std::abs(c), std::abs(d)});
    /**
     * @brief 统一按二进制幂缩放，将最大系数绝对值置于 [1,2)。
     * 缩放后处于正规范围的数值保留二进制有效位。
     * @see David Goldberg, "What Every Computer Scientist Should Know About
     * Floating-Point Arithmetic" (1991)，定理 7 的证明。
     * https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
     * @see S. Ghaderpanah and S. Klasa, "Polynomial Scaling" (1990),
     * SIAM Journal on Numerical Analysis 27(1), 117-135（缩放背景）。
     * https://doi.org/10.1137/0727007
     */
    const int coefficient_exponent = std::ilogb(max_coefficient);
    const double an = std::scalbn(a, -coefficient_exponent);
    const double bn = std::scalbn(b, -coefficient_exponent);
    const double cn = std::scalbn(c, -coefficient_exponent);
    const double dn = std::scalbn(d, -coefficient_exponent);
    if (an == 0.0) return std::nullopt;
    const double an_squared = an * an;
    const double an_cubed = an_squared * an;
    const double bn_squared = bn * bn;
    const double p_left = 3.0 * an * cn;
    double p_numerator = std::fma(-bn, bn, p_left);
    const double p_roundoff =
        16.0 * std::numeric_limits<double>::epsilon() *
        (std::abs(p_left) + std::abs(bn_squared));
    if (std::abs(p_numerator) <= p_roundoff) {
        p_numerator = 0.0;
    }

    const double q_first = 2.0 * bn_squared * bn;
    const double q_second = -9.0 * an * bn * cn;
    const double q_third = 27.0 * an_squared * dn;
    double q_numerator = (q_first + q_second) + q_third;
    const double q_roundoff =
        32.0 * std::numeric_limits<double>::epsilon() *
        (std::abs(q_first) + std::abs(q_second) +
         std::abs(q_third));
    if (std::abs(q_numerator) <= q_roundoff) {
        q_numerator = 0.0;
    }

    const double p_denominator = 3.0 * an_squared;
    const double q_denominator = 27.0 * an_cubed;
    if (p_denominator == 0.0 || q_denominator == 0.0) return std::nullopt;
    p_val = p_numerator / p_denominator;
    q_val = q_numerator / q_denominator;
    shift_val = bn / (3.0 * an);
    numeric_depression =
        std::isfinite(p_val) && std::isfinite(q_val) &&
        std::isfinite(shift_val);
    if (!numeric_depression) return std::nullopt;
    return NumericCubicDepression{p_val, q_val, shift_val};
}

std::optional<NumericCubicDepression> numeric_cubic_depression(
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    const std::shared_ptr<SymbolicExpr>& c, const std::shared_ptr<SymbolicExpr>& d) {
    const bool all_numeric = is_purely_numeric(a) && is_purely_numeric(b) &&
                             is_purely_numeric(c) && is_purely_numeric(d);
    auto av = finite_numeric_value(a);
    auto bv = finite_numeric_value(b);
    auto cv = finite_numeric_value(c);
    auto dv = finite_numeric_value(d);
    if (!all_numeric || !av || !bv || !cv || !dv || *av == 0.0) return std::nullopt;
    return depress_numeric_cubic(*av, *bv, *cv, *dv);
}

static std::vector<std::shared_ptr<SymbolicExpr>> cubic_one_real_root(
    double discriminant, double q_half, double root_scale, double shift_val) {
    const double sqrt_discriminant = std::sqrt(discriminant);
    const double u = std::cbrt(-q_half + sqrt_discriminant);
    const double v = std::cbrt(-q_half - sqrt_discriminant);
    const double real_sum = (u + v) * root_scale;
    const double x1_value = real_sum - shift_val;
    const double real_part = -real_sum / 2.0 - shift_val;
    const double imaginary_part =
        std::sqrt(3.0) * (u - v) * root_scale / 2.0;

    auto x1 = SymbolicExpr::number(x1_value);
    auto i_unit = SymbolicExpr::sqrt(num(-1));
    auto x2 = SymbolicExpr::add(
        SymbolicExpr::number(real_part),
        SymbolicExpr::multiply(
            SymbolicExpr::number(imaginary_part), i_unit))->simplify();
    auto x3 = sub(
        SymbolicExpr::number(real_part),
        SymbolicExpr::multiply(
            SymbolicExpr::number(imaginary_part), i_unit))->simplify();
    return {x1, x2, x3};
}

static std::vector<std::shared_ptr<SymbolicExpr>> cubic_three_real_roots(
    double p_normalized, double q_normalized, double root_scale, double shift_val) {
    std::vector<std::shared_ptr<SymbolicExpr>> roots;
    const double radius =
        std::sqrt(-(p_normalized * p_normalized * p_normalized) / 27.0);
    double cosine_argument = -q_normalized / (2.0 * radius);
    cosine_argument = std::clamp(cosine_argument, -1.0, 1.0);
    const double theta = std::acos(cosine_argument);
    const double amplitude = 2.0 * std::cbrt(radius) * root_scale;

    for (int k = 0; k < 3; ++k) {
        const double angle =
            (theta + 2.0 * k * LMMC_CONST_PI) / 3.0;
        roots.push_back(SymbolicExpr::number(
            amplitude * std::cos(angle) - shift_val));
    }
    return roots;
}

std::vector<std::shared_ptr<SymbolicExpr>> numeric_cubic_roots(
    const NumericCubicDepression& depression) {
    const auto& [p_val, q_val, shift_val] = depression;
    const double root_scale = std::max(
        std::sqrt(std::abs(p_val)),
        std::cbrt(std::abs(q_val)));
    if (root_scale == 0.0) {
        auto root = SymbolicExpr::number(-shift_val);
        return {root, root, root};
    }
    const double p_normalized = (p_val / root_scale) / root_scale;
    const double q_normalized =
        ((q_val / root_scale) / root_scale) / root_scale;
    const double q_half = q_normalized / 2.0;
    const double p_third = p_normalized / 3.0;
    const double discriminant_q = q_half * q_half;
    const double discriminant_p = p_third * p_third * p_third;
    const double discriminant = discriminant_q + discriminant_p;
    const double discriminant_tolerance =
        32.0 * std::numeric_limits<double>::epsilon() *
        (std::abs(discriminant_q) + std::abs(discriminant_p));
    if (std::abs(discriminant) <= discriminant_tolerance) {
        const double cbrt_q_half = std::cbrt(q_half);
        const double t1 = -2.0 * cbrt_q_half * root_scale;
        const double t2 = cbrt_q_half * root_scale;

        auto x1 = SymbolicExpr::number(t1 - shift_val);
        auto x2 = SymbolicExpr::number(t2 - shift_val);
        return {x1, x2, x2};
    }
    if (discriminant > 0.0) {
        return cubic_one_real_root(discriminant, q_half, root_scale, shift_val);
    }
    return cubic_three_real_roots(p_normalized, q_normalized, root_scale, shift_val);
}

}
