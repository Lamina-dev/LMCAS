#include "newton_raphson.hpp"
#include "poly_utils.hpp"
#include <algorithm>
#include <cmath>

namespace lamina {

static int count_sign_changes(const std::vector<Rational>& values) {
    int changes = 0;
    int last_sign = 0;

    for (const auto& v : values) {
        if (v.is_zero()) continue;

        int current_sign = (v > Rational(0)) ? 1 : -1;
        if (last_sign != 0 && current_sign != last_sign) {
            ++changes;
        }
        last_sign = current_sign;
    }
    return changes;
}

static int sturm_sign_changes_at(
    const std::vector<Polynomial<Rational>>& sturm,
    const Rational& x)
{
    std::vector<Rational> values;
    values.reserve(sturm.size());
    for (const auto& p : sturm) {
        values.push_back(p.eval(x));
    }
    return count_sign_changes(values);
}

static int sturm_sign_changes_at_pos_inf(
    const std::vector<Polynomial<Rational>>& sturm)
{
    std::vector<Rational> signs;
    signs.reserve(sturm.size());
    for (const auto& p : sturm) {
        if (p.is_zero()) continue;
        signs.push_back(p.lead_coeff());
    }
    return count_sign_changes(signs);
}

static int sturm_sign_changes_at_neg_inf(
    const std::vector<Polynomial<Rational>>& sturm)
{
    std::vector<Rational> signs;
    signs.reserve(sturm.size());
    for (const auto& p : sturm) {
        if (p.is_zero()) continue;
        int deg = p.degree();
        Rational lc = p.lead_coeff();

        if (deg % 2 == 1) {
            lc = -lc;
        }
        signs.push_back(lc);
    }
    return count_sign_changes(signs);
}

static Rational cauchy_bound(const Polynomial<Rational>& poly) {
    if (poly.degree() <= 0) return Rational(1);

    Rational lc = poly.lead_coeff();
    Rational max_ratio(0);

    for (int i = 0; i < poly.degree(); ++i) {
        Rational coeff = poly.coeffs[i];
        if (coeff.is_zero()) continue;
        Rational ratio = coeff.abs() / lc.abs();
        if (ratio > max_ratio) {
            max_ratio = ratio;
        }
    }

    return Rational(1) + max_ratio;
}

static int roots_in_interval(
    const std::vector<Polynomial<Rational>>& sturm,
    const Rational& a,
    const Rational& b)
{
    int va = sturm_sign_changes_at(sturm, a);
    int vb = sturm_sign_changes_at(sturm, b);
    return va - vb;
}

std::vector<std::pair<Rational, Rational>> isolate_real_roots(
    const Polynomial<Rational>& poly)
{
    std::vector<std::pair<Rational, Rational>> result;

    if (poly.is_zero() || poly.degree() <= 0) {
        return result;
    }

    Polynomial<Rational> sqfree = poly.square_free_part();
    if (sqfree.is_zero() || sqfree.degree() <= 0) {
        return result;
    }

    std::vector<Polynomial<Rational>> sturm;
    sturm.push_back(sqfree);
    sturm.push_back(sqfree.differentiate());

    while (true) {
        size_t n = sturm.size();
        const auto& prev2 = sturm[n - 2];
        const auto& prev1 = sturm[n - 1];

        if (prev1.is_zero()) break;

        auto [quotient, remainder] = prev2.div_mod(prev1);

        if (remainder.is_zero()) break;

        Polynomial<Rational> neg_rem(remainder.variable_name);
        neg_rem.coeffs.resize(remainder.coeffs.size());
        for (size_t i = 0; i < remainder.coeffs.size(); ++i) {
            neg_rem.coeffs[i] = -remainder.coeffs[i];
        }
        neg_rem.trim();

        sturm.push_back(neg_rem);
    }

    Rational bound = cauchy_bound(sqfree);
    Rational lo = -bound;
    Rational hi = bound;

    int total_roots = roots_in_interval(sturm, lo, hi);
    if (total_roots <= 0) {
        return result;
    }

    struct Interval {
        Rational lo, hi;
        int root_count;
    };

    std::vector<Interval> work_queue;
    work_queue.push_back({lo, hi, total_roots});

    while (!work_queue.empty()) {
        Interval current = work_queue.back();
        work_queue.pop_back();

        if (current.root_count == 0) {
            continue;
        }

        if (current.root_count == 1) {
            result.push_back({current.lo, current.hi});
            continue;
        }

        Rational mid = (current.lo + current.hi) / Rational(2);

        int left_count = roots_in_interval(sturm, current.lo, mid);
        int right_count = roots_in_interval(sturm, mid, current.hi);

        if (left_count + right_count < current.root_count) {

            result.push_back({mid, mid});
        }

        if (right_count > 0) {
            work_queue.push_back({mid, current.hi, right_count});
        }
        if (left_count > 0) {
            work_queue.push_back({current.lo, mid, left_count});
        }
    }

    std::sort(result.begin(), result.end(),
        [](const std::pair<Rational, Rational>& a,
           const std::pair<Rational, Rational>& b) {
            return a.first < b.first;
        });

    return result;
}

std::optional<NumericRoot> bisection(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    lmmc_real_t lo,
    lmmc_real_t hi,
    const SolveOptions& opts)
{

    lmmc_real_t f_lo = f->substitute(var, SymbolicExpr::number(lo))->to_numeric();
    lmmc_real_t f_hi = f->substitute(var, SymbolicExpr::number(hi))->to_numeric();

    if (std::abs(f_lo) < opts.tolerance) {
        return NumericRoot{lo, std::abs(f_lo), 0};
    }
    if (std::abs(f_hi) < opts.tolerance) {
        return NumericRoot{hi, std::abs(f_hi), 0};
    }

    if (f_lo * f_hi > 0) {
        return std::nullopt;
    }

    int max_iter = opts.max_newton_iterations * 3;
    for (int i = 1; i <= max_iter; ++i) {
        lmmc_real_t mid = (lo + hi) * 0.5;
        lmmc_real_t f_mid = f->substitute(var, SymbolicExpr::number(mid))->to_numeric();

        if (std::abs(f_mid) < opts.tolerance) {
            return NumericRoot{mid, std::abs(f_mid), i};
        }

        if (std::abs(hi - lo) < opts.tolerance) {
            return NumericRoot{mid, std::abs(f_mid), i};
        }

        if (f_lo * f_mid < 0) {
            hi = mid;
            f_hi = f_mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }

    lmmc_real_t mid = (lo + hi) * 0.5;
    lmmc_real_t f_mid = f->substitute(var, SymbolicExpr::number(mid))->to_numeric();
    if (std::abs(f_mid) < opts.tolerance * 1000) {
        return NumericRoot{mid, std::abs(f_mid), max_iter};
    }
    return std::nullopt;
}

std::optional<NumericRoot> newton_raphson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    lmmc_real_t bracket_lo,
    lmmc_real_t bracket_hi,
    const SolveOptions& opts)
{
    lmmc_real_t x = x0;
    for (int i = 1; i <= opts.max_newton_iterations; ++i) {

        auto fx_expr = f->substitute(var, SymbolicExpr::number(x));
        auto dfx_expr = df->substitute(var, SymbolicExpr::number(x));

        lmmc_real_t fx = fx_expr->to_numeric();
        lmmc_real_t dfx = dfx_expr->to_numeric();

        if (std::abs(fx) < opts.tolerance) {
            return NumericRoot{x, std::abs(fx), i};
        }

        if (std::abs(dfx) < 1e-15) {

            return bisection(f, var, bracket_lo, bracket_hi, opts);
        }

        lmmc_real_t x_new = x - fx / dfx;

        if (i > 1 && std::abs(x_new - x) > 2.0 * std::abs(x - x0)) {
            x_new = x - 0.5 * fx / dfx;
        }

        x = x_new;
    }

    return std::nullopt;
}

std::optional<NumericRoot> newton_raphson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    const SolveOptions& opts)
{
    lmmc_real_t x = x0;
    for (int i = 1; i <= opts.max_newton_iterations; ++i) {

        auto fx_expr = f->substitute(var, SymbolicExpr::number(x));
        auto dfx_expr = df->substitute(var, SymbolicExpr::number(x));

        lmmc_real_t fx = fx_expr->to_numeric();
        lmmc_real_t dfx = dfx_expr->to_numeric();

        if (std::abs(fx) < opts.tolerance) {
            return NumericRoot{x, std::abs(fx), i};
        }

        if (std::abs(dfx) < 1e-15) {

            return std::nullopt;
        }

        lmmc_real_t x_new = x - fx / dfx;

        if (i > 1 && std::abs(x_new - x) > 2.0 * std::abs(x - x0)) {
            x_new = x - 0.5 * fx / dfx;
        }

        x = x_new;
    }

    return std::nullopt;
}

std::vector<NumericRoot> solve_numeric(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts)
{
    std::vector<NumericRoot> results;

    auto poly = symbolic_to_poly<Rational>(expr, var);

    if (!poly.is_zero() && poly.degree() >= 1) {

        Polynomial<Rational> current_poly = poly;

        auto df_expr = expr->differentiate(var);

        while (current_poly.degree() >= 1) {

            auto intervals = isolate_real_roots(current_poly);

            if (intervals.empty()) {
                break;
            }

            bool found_any = false;
            for (const auto& [lo_rat, hi_rat] : intervals) {

                if (opts.max_roots > 0 && (int)results.size() >= opts.max_roots) {
                    return results;
                }

                lmmc_real_t lo = lo_rat.to_double();
                lmmc_real_t hi = hi_rat.to_double();
                lmmc_real_t x0 = (lo + hi) * 0.5;

                auto current_expr = poly_to_symbolic(current_poly);
                auto current_df = current_expr->differentiate(var);

                auto root_opt = newton_raphson(current_expr, current_df, var,
                                              x0, lo, hi, opts);

                if (root_opt.has_value()) {
                    NumericRoot root = root_opt.value();

                    auto check = expr->substitute(var, SymbolicExpr::number(root.value));
                    lmmc_real_t residual = std::abs(check->to_numeric());
                    root.residual = residual;

                    if (residual < opts.tolerance * 1000) {
                        results.push_back(root);
                        found_any = true;

                        Rational r_rat = Rational::from_double(root.value);
                        Polynomial<Rational> linear_factor({-r_rat, Rational(1)},
                                                          current_poly.variable_name);
                        auto [quotient, remainder] = current_poly.div_mod(linear_factor);

                        if (!remainder.is_zero()) {

                            int deg = current_poly.degree();
                            std::vector<Rational> new_coeffs(deg);

                            new_coeffs[deg - 1] = current_poly.coeffs[deg];
                            for (int i = deg - 2; i >= 0; --i) {
                                new_coeffs[i] = current_poly.coeffs[i + 1] + new_coeffs[i + 1] * r_rat;
                            }
                            quotient = Polynomial<Rational>(new_coeffs, current_poly.variable_name);
                        }

                        current_poly = quotient;
                        break;
                    }
                }
            }

            if (!found_any) {
                break;
            }
        }
    } else {

        lmmc_real_t x0 = opts.has_initial_guess ? opts.initial_guess : 0.0;

        auto df_expr = expr->differentiate(var);

        auto root_opt = newton_raphson(expr, df_expr, var, x0, opts);

        if (root_opt.has_value()) {
            NumericRoot root = root_opt.value();

            if (root.residual < opts.tolerance) {
                results.push_back(root);
            }
        }
    }

    return results;
}

}
