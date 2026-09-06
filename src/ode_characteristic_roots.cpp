#include "internal/ode_characteristic_roots.hpp"

#include "lmmc/numeric.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace LMCAS::ode_root_detail {
/**
 * @internal
 * @brief 使用数值方法求解特征多项式的所有根.
 *
 * 对于一、二次多项式使用解析公式；对于三至六次多项式使用
 * Durand-Kerner 迭代。迭代前通过多项式与其导数的近似 GCD 提取重根，
 * 避免把一个重根的数值扰动误分类为多个复根。
 */
Result<std::vector<CharRoot>>
find_characteristic_roots(
    const std::vector<double>& coeffs, const std::string& operation)
{
    int n = static_cast<int>(coeffs.size()) - 1;
    const double root_verification_tolerance = 1e-8;
    const double roundoff_tolerance =
        256.0 * std::numeric_limits<double>::epsilon() *
        static_cast<double>(std::max(n, 1));
    if (n <= 0) {
        return Result<std::vector<CharRoot>>::failure(
            CasErrc::NumericFailure,
            "characteristic polynomial has no roots", operation);
    }

    std::vector<CharRoot> roots;
    double leading = coeffs[0];
    if (leading == 0.0) {
        return Result<std::vector<CharRoot>>::failure(
            CasErrc::NumericFailure,
            "characteristic polynomial has no usable leading coefficient",
            operation);
    }
    std::vector<double> norm_coeffs(coeffs.size());
    for (size_t i = 0; i < coeffs.size(); ++i) {
        norm_coeffs[i] = coeffs[i] / leading;
        if (!std::isfinite(norm_coeffs[i])) {
            return Result<std::vector<CharRoot>>::failure(
                CasErrc::NumericFailure,
                "characteristic polynomial normalization is non-finite",
                operation);
        }
    }

    /// 对于低阶多项式,使用解析公式
    if (n == 1) {
        /// r + norm_coeffs[1] = 0
        double r = -norm_coeffs[1];
        roots.push_back({r, 0.0, 1, false});
        return Result<std::vector<CharRoot>>::success(std::move(roots));
    }

    if (n == 2) {
        const double b = norm_coeffs[1];
        const double c = norm_coeffs[2];
        const double scale = std::max(std::abs(b), std::sqrt(std::abs(c)));
        const double scaled_b = scale == 0.0 ? 0.0 : b / scale;
        const double scaled_c =
            scale == 0.0 ? 0.0 : (c / scale) / scale;
        const double discriminant =
            std::fma(scaled_b, scaled_b, -4.0 * scaled_c);
        const double discriminant_tolerance =
            32.0 * std::numeric_limits<double>::epsilon() *
            (scaled_b * scaled_b + 4.0 * std::abs(scaled_c));
        if (std::abs(discriminant) <= discriminant_tolerance) {
            roots.push_back({-b / 2.0, 0.0, 2, false});
        } else if (discriminant > 0.0) {
            const double sqrt_discriminant =
                scale * std::sqrt(discriminant);
            const double q =
                -0.5 * (b + std::copysign(sqrt_discriminant, b));
            roots.push_back({q, 0.0, 1, false});
            roots.push_back({c / q, 0.0, 1, false});
        } else {
            roots.push_back(
                {-b / 2.0, scale * std::sqrt(-discriminant) / 2.0,
                 1, true});
        }
        for (const auto& root : roots) {
            if (polynomial_root_backward_error(
                    norm_coeffs, {root.real_part, root.imag_part}) >
                root_verification_tolerance) {
                return Result<std::vector<CharRoot>>::failure(
                    CasErrc::NumericFailure,
                    "quadratic characteristic root failed verification",
                    operation);
            }
        }
        return Result<std::vector<CharRoot>>::success(std::move(roots));
    }
    double root_scale_estimate = 1.0;
    for (int i = 1; i <= n; ++i) {
        const double magnitude = std::abs(norm_coeffs[i]);
        if (magnitude > 0.0) {
            root_scale_estimate = std::max(
                root_scale_estimate,
                std::exp(std::log(magnitude) / static_cast<double>(i)));
        }
    }
    if (!std::isfinite(root_scale_estimate)) {
        return Result<std::vector<CharRoot>>::failure(
            CasErrc::NumericFailure,
            "characteristic-root variable scaling is non-finite",
            operation);
    }
    int root_scale_exponent = std::ilogb(root_scale_estimate);
    double root_scale = std::scalbn(1.0, root_scale_exponent);
    if (root_scale_exponent < std::numeric_limits<double>::max_exponent - 1 &&
        root_scale_estimate / root_scale > std::sqrt(2.0)) {
        root_scale *= 2.0;
    }

    // For p(r)=r^n+c1*r^(n-1)+...+cn, r=root_scale*z gives
    // z^n+(c1/root_scale)z^(n-1)+...+cn/root_scale^n.
    std::vector<double> solver_coeffs = norm_coeffs;
    for (int i = 1; i <= n; ++i) {
        for (int power = 0; power < i; ++power) {
            solver_coeffs[i] /= root_scale;
        }
    }

    auto restore_root_scale =
        [&](std::vector<CharRoot> scaled_roots)
            -> Result<std::vector<CharRoot>> {
        for (auto& root : scaled_roots) {
            root.real_part *= root_scale;
            root.imag_part *= root_scale;
            if (!std::isfinite(root.real_part) ||
                !std::isfinite(root.imag_part)) {
                return Result<std::vector<CharRoot>>::failure(
                    CasErrc::NumericFailure,
                    "characteristic root exceeds the finite numeric range",
                    operation);
            }
        }
        return Result<std::vector<CharRoot>>::success(
            std::move(scaled_roots));
    };

    RealPolynomial repeated_factor = approximate_polynomial_gcd(
        solver_coeffs, polynomial_derivative(solver_coeffs));
    if (repeated_factor.size() > 1 &&
        repeated_factor.size() < solver_coeffs.size()) {
        auto square_free_division =
            divide_polynomials(solver_coeffs, repeated_factor);
        if (polynomial_remainder_is_small(
                square_free_division.second, solver_coeffs)) {
            auto distinct_result = find_characteristic_roots(
                square_free_division.first, operation);
            if (!distinct_result) {
                return Result<std::vector<CharRoot>>::failure(
                    distinct_result.error());
            }
            auto distinct_roots = std::move(distinct_result.value());
            RealPolynomial remaining = solver_coeffs;
            int recovered_degree = 0;
            bool recovered_all = !distinct_roots.empty();
            for (auto& root : distinct_roots) {
                RealPolynomial factor;
                if (root.is_complex) {
                    const double magnitude_squared =
                        root.real_part * root.real_part +
                        root.imag_part * root.imag_part;
                    if (!std::isfinite(magnitude_squared)) {
                        recovered_all = false;
                        break;
                    }
                    factor = {
                        1.0, -2.0 * root.real_part, magnitude_squared};
                } else {
                    factor = {1.0, -root.real_part};
                }

                int multiplicity = 0;
                while (remaining.size() >= factor.size()) {
                    auto division = divide_polynomials(remaining, factor);
                    if (!polynomial_remainder_is_small(
                            division.second, remaining)) {
                        break;
                    }
                    remaining = std::move(division.first);
                    ++multiplicity;
                }
                if (multiplicity == 0) {
                    recovered_all = false;
                    break;
                }
                root.multiplicity = multiplicity;
                recovered_degree +=
                    multiplicity * (root.is_complex ? 2 : 1);
            }
            if (recovered_all && recovered_degree == n) {
                return restore_root_scale(std::move(distinct_roots));
            }
        }
    }


    std::vector<Complex> z(n);
    /// 初始猜测采用不同半径,提供非对称起点以区分各根.
    for (int i = 0; i < n; ++i) {
        double angle = 2.0 * 3.14159265358979323846 * i / n + 0.1;
        double radius = 1.0 + 0.3 * i;
        z[i] = {radius * std::cos(angle), radius * std::sin(angle)};
    }

    /// 求值多项式 p(z)
    auto eval_poly = [&](const Complex& val) -> Complex {
        Complex result = {1.0, 0.0};
        for (int i = 1; i <= n; ++i) {
            result = result * val + Complex{solver_coeffs[i], 0.0};
        }
        return result;
    };

    /// Durand-Kerner 迭代
    bool converged = false;
    for (int iter = 0; iter < 1000; ++iter) {
        double max_change = 0.0;
        for (int i = 0; i < n; ++i) {
            Complex num = eval_poly(z[i]);
            Complex denom = {1.0, 0.0};
            for (int j = 0; j < n; ++j) {
                if (j != i) denom = denom * (z[i] - z[j]);
            }
            Complex delta = num / denom;
            double change = delta.magnitude();
            if (!std::isfinite(change)) {
                return Result<std::vector<CharRoot>>::failure(
                    CasErrc::NumericFailure,
                    "characteristic-root iteration became non-finite",
                    operation);
            }
            z[i] = z[i] - delta;
            if (change > max_change) max_change = change;
        }
        if (max_change < 1e-12) {
            converged = true;
            break;
        }
    }
    if (!converged) {
        return Result<std::vector<CharRoot>>::failure(
            CasErrc::NumericFailure,
            "characteristic-root iteration did not converge", operation);
    }
    for (const auto& candidate : z) {
        if (polynomial_root_backward_error(solver_coeffs, candidate) >
            root_verification_tolerance) {
            return Result<std::vector<CharRoot>>::failure(
                CasErrc::NumericFailure,
                "characteristic root failed backward-error verification",
                operation);
        }
    }

    /// 使用原多项式验证实轴投影；固定虚部容差会吞掉尺度较小的复根。
    /// 重数只由前面的多项式 GCD/因子恢复确定，不能按根间距离猜测。
    std::vector<bool> used(n, false);
    for (int i = 0; i < n; ++i) {
        if (used[i]) continue;

        const Complex real_projection{z[i].re, 0.0};
        const double candidate_error =
            polynomial_root_backward_error(solver_coeffs, z[i]);
        const double projection_error =
            polynomial_root_backward_error(solver_coeffs, real_projection);
        const double projection_tolerance =
            std::max(roundoff_tolerance, 16.0 * candidate_error);
        if (projection_error <= projection_tolerance) {
            used[i] = true;
            roots.push_back({z[i].re, 0.0, 1, false});
            continue;
        }

        int conjugate_index = -1;
        double best_distance = std::numeric_limits<double>::infinity();
        for (int j = i + 1; j < n; ++j) {
            if (used[j] || std::signbit(z[i].im) == std::signbit(z[j].im)) {
                continue;
            }
            const double scale =
                1.0 + std::max(z[i].magnitude(), z[j].magnitude());
            const double distance =
                std::hypot(z[i].re - z[j].re, z[i].im + z[j].im) /
                scale;
            if (distance < best_distance) {
                best_distance = distance;
                conjugate_index = j;
            }
        }
        if (conjugate_index < 0) {
            return Result<std::vector<CharRoot>>::failure(
                CasErrc::NumericFailure,
                "characteristic-root classification found an unpaired complex root",
                operation);
        }

        const Complex paired_root{
            0.5 * (z[i].re + z[conjugate_index].re),
            0.5 * std::abs(z[i].im - z[conjugate_index].im)};
        if (polynomial_root_backward_error(
                solver_coeffs, paired_root) > root_verification_tolerance) {
            return Result<std::vector<CharRoot>>::failure(
                CasErrc::NumericFailure,
                "conjugate characteristic-root pair failed verification",
                operation);
        }
        used[i] = true;
        used[conjugate_index] = true;
        roots.push_back(
            {paired_root.re, paired_root.im, 1, true});
    }

    int classified_degree = 0;
    for (const auto& root : roots) {
        classified_degree += root.multiplicity * (root.is_complex ? 2 : 1);
    }
    if (classified_degree != n) {
        return Result<std::vector<CharRoot>>::failure(
            CasErrc::NumericFailure,
            "characteristic-root classification lost roots", operation);
    }
    return restore_root_scale(std::move(roots));
}

} // namespace LMCAS::ode_root_detail
