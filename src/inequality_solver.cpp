#include "../include/inequality_solver.hpp"
#include "../include/poly_utils.hpp"
#include "../include/solve_polynomial.hpp"
#include "../include/solve_strategies.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <limits>
#include <functional>

namespace lamina {

static int determine_leading_sign(const Polynomial<SymbolicPolyCoeff>& poly) {
    if (poly.is_zero()) return 0;
    auto lc = poly.lead_coeff().val;
    if (!lc) return 1;
    auto simplified = lc->simplify();
    if (!simplified) return 1;

    try {
        double val = simplified->to_numeric();
        if (val > 0) return 1;
        if (val < 0) return -1;
    } catch (...) {}

    if (auto num = std::dynamic_pointer_cast<NumberNode>(simplified->root)) {
        if (std::holds_alternative<BigInt>(num->value)) {
            return std::get<BigInt>(num->value).IsNegative() ? -1 : 1;
        }
        if (std::holds_alternative<Rational>(num->value)) {
            return std::get<Rational>(num->value).get_numerator().IsNegative() ? -1 : 1;
        }
        if (std::holds_alternative<lmmc_real_t>(num->value)) {
            return std::get<lmmc_real_t>(num->value) < 0 ? -1 : 1;
        }
    }
    return 1;
}

static std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>> find_roots_with_multiplicity(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& variable) {

    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>> result;

    auto poly_rat = symbolic_to_poly<Rational>(expr, variable);
    if (poly_rat.is_zero() || poly_rat.degree() <= 0) {
        return result;
    }

    auto factors = square_free_factorization(poly_rat);
    if (factors.empty()) {

        auto poly_spc = symbolic_to_poly<SymbolicPolyCoeff>(expr, variable);
        if (!poly_spc.is_zero() && poly_spc.degree() >= 1) {
            auto roots = solve_by_factoring(poly_spc, variable);
            for (const auto& root : roots) {
                if (!root) continue;
                try {
                    double val = root->to_numeric();
                    if (std::isfinite(val)) {
                        result.push_back({root, 1});
                    }
                } catch (...) {}
            }
        }
        return result;
    }

    for (const auto& [factor, mult] : factors) {
        if (factor.degree() <= 0) continue;

        if (factor.degree() == 1) {

            Rational a = factor.coeffs[1];
            Rational b = factor.coeffs[0];
            if (a != Rational(0)) {
                Rational root_val = Rational(0) - b / a;
                auto root_expr = SymbolicExpr::number(root_val);
                try {
                    double val = root_expr->to_numeric();
                    if (std::isfinite(val)) {
                        result.push_back({root_expr, mult});
                    }
                } catch (...) {}
            }
            continue;
        }

        if (factor.degree() == 2) {
            Rational a = factor.coeffs[2];
            Rational b = factor.coeffs[1];
            Rational c = factor.coeffs[0];

            Rational disc = b * b - Rational(4) * a * c;
            double disc_val = disc.to_double();
            if (disc_val < -1e-10) continue;
            if (disc_val < 0) disc_val = 0;

            double a_val = a.to_double();
            double b_val = b.to_double();
            double sqrt_disc = std::sqrt(disc_val);

            double r1 = (-b_val + sqrt_disc) / (2.0 * a_val);
            double r2 = (-b_val - sqrt_disc) / (2.0 * a_val);

            if (std::isfinite(r1)) {
                result.push_back({SymbolicExpr::number(r1), mult});
            }
            if (std::isfinite(r2) && std::abs(r1 - r2) > 1e-10) {
                result.push_back({SymbolicExpr::number(r2), mult});
            }
            continue;
        }

        std::vector<SymbolicPolyCoeff> spc_coeffs;
        for (int i = 0; i <= factor.degree(); ++i) {
            spc_coeffs.push_back(SymbolicPolyCoeff(SymbolicExpr::number(factor.coeffs[i])));
        }
        Polynomial<SymbolicPolyCoeff> factor_spc(spc_coeffs, variable);

        auto factor_roots = solve_by_factoring(factor_spc, variable);
        for (const auto& root : factor_roots) {
            if (!root) continue;
            try {
                double val = root->to_numeric();
                if (std::isfinite(val)) {
                    result.push_back({root, mult});
                }
            } catch (...) {}
        }
    }

    return result;
}

static bool root_less_than(const std::shared_ptr<SymbolicExpr>& a,
                           const std::shared_ptr<SymbolicExpr>& b) {
    try {
        double va = a->to_numeric();
        double vb = b->to_numeric();
        return va < vb;
    } catch (...) {
        return false;
    }
}

static bool roots_equal(const std::shared_ptr<SymbolicExpr>& a,
                        const std::shared_ptr<SymbolicExpr>& b) {
    try {
        double va = a->to_numeric();
        double vb = b->to_numeric();
        return std::abs(va - vb) < 1e-10;
    } catch (...) {
        return false;
    }
}

std::vector<SignChartEntry> InequalitySolver::build_sign_chart(
    const std::shared_ptr<SymbolicExpr>& poly,
    const std::string& variable,
    const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
    const std::vector<int>& multiplicities) {

    std::vector<SignChartEntry> chart;

    if (roots.empty()) {

        auto p = symbolic_to_poly<SymbolicPolyCoeff>(poly, variable);
        int sign = determine_leading_sign(p);
        chart.push_back({Interval::entire_line(), sign});
        return chart;
    }

    auto p = symbolic_to_poly<SymbolicPolyCoeff>(poly, variable);
    int leading_sign = determine_leading_sign(p);

    int current_sign = leading_sign;

    size_t n = roots.size();
    std::vector<int> interval_signs(n + 1);

    interval_signs[n] = leading_sign;

    for (int i = (int)n - 1; i >= 0; --i) {
        interval_signs[i] = interval_signs[i + 1];
        if (multiplicities[i] % 2 != 0) {
            interval_signs[i] = -interval_signs[i];
        }
    }

    {
        Interval iv;
        iv.lower = Endpoint::neg_inf();
        iv.upper = Endpoint::open(roots[0]);
        chart.push_back({iv, interval_signs[0]});
    }

    for (size_t i = 0; i + 1 < n; ++i) {
        Interval iv;
        iv.lower = Endpoint::open(roots[i]);
        iv.upper = Endpoint::open(roots[i + 1]);
        chart.push_back({iv, interval_signs[i + 1]});
    }

    {
        Interval iv;
        iv.lower = Endpoint::open(roots[n - 1]);
        iv.upper = Endpoint::pos_inf();
        chart.push_back({iv, interval_signs[n]});
    }

    return chart;
}

IntervalUnion InequalitySolver::select_intervals(
    const std::vector<SignChartEntry>& chart,
    InequalityType type,
    const std::vector<std::shared_ptr<SymbolicExpr>>& roots,
    const std::vector<int>& multiplicities) {

    std::vector<Interval> result_intervals;

    bool want_positive = (type == InequalityType::GreaterThan || type == InequalityType::GreaterEqual);
    bool is_strict = (type == InequalityType::GreaterThan || type == InequalityType::LessThan);
    int target_sign = want_positive ? 1 : -1;

    for (const auto& entry : chart) {
        if (entry.sign == target_sign) {
            result_intervals.push_back(entry.interval);
        }
    }

    if (!is_strict) {

        for (size_t i = 0; i < roots.size(); ++i) {

            bool merged = false;
            for (auto& iv : result_intervals) {

                if (!iv.upper.is_pos_infinity && iv.upper.value) {
                    if (roots_equal(iv.upper.value, roots[i])) {
                        iv.upper.is_open = false;
                        merged = true;
                    }
                }

                if (!iv.lower.is_neg_infinity && iv.lower.value) {
                    if (roots_equal(iv.lower.value, roots[i])) {
                        iv.lower.is_open = false;
                        merged = true;
                    }
                }
            }

            if (!merged) {
                result_intervals.push_back(Interval::point(roots[i]));
            }
        }
    }

    return IntervalUnion(result_intervals);
}

IntervalUnion InequalitySolver::solve_inequality(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable) {

    if (!expr) return IntervalUnion::empty();

    auto poly = symbolic_to_poly<SymbolicPolyCoeff>(expr, variable);
    if (poly.is_zero()) {

        if (depends_on_var(expr->root, variable)) {
            return IntervalUnion::empty();
        }

        if (type == InequalityType::GreaterEqual || type == InequalityType::LessEqual) {
            return IntervalUnion::entire_line();
        }
        return IntervalUnion::empty();
    }

    if (poly.degree() <= 0) {

        auto lc = poly.lead_coeff().val;
        if (!lc) return IntervalUnion::empty();
        try {
            double val = lc->simplify()->to_numeric();
            bool satisfies = false;
            switch (type) {
                case InequalityType::GreaterThan: satisfies = (val > 0); break;
                case InequalityType::GreaterEqual: satisfies = (val >= 0); break;
                case InequalityType::LessThan: satisfies = (val < 0); break;
                case InequalityType::LessEqual: satisfies = (val <= 0); break;
            }
            return satisfies ? IntervalUnion::entire_line() : IntervalUnion::empty();
        } catch (...) {
            return IntervalUnion::empty();
        }
    }

    {
        auto poly_rat = symbolic_to_poly<Rational>(expr, variable);
        if (poly_rat.is_zero() && depends_on_var(expr->root, variable)) {

            return IntervalUnion::empty();
        }

    }

    auto roots_with_mult = find_roots_with_multiplicity(expr, variable);

    std::sort(roots_with_mult.begin(), roots_with_mult.end(),
        [](const auto& a, const auto& b) {
            return root_less_than(a.first, b.first);
        });

    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, int>> unique_roots;
    for (const auto& [root, mult] : roots_with_mult) {
        if (!unique_roots.empty() && roots_equal(unique_roots.back().first, root)) {

            unique_roots.back().second = std::max(unique_roots.back().second, mult);
        } else {
            unique_roots.push_back({root, mult});
        }
    }

    std::vector<std::shared_ptr<SymbolicExpr>> roots;
    std::vector<int> multiplicities;
    for (const auto& [root, mult] : unique_roots) {
        roots.push_back(root);
        multiplicities.push_back(mult);
    }

    auto chart = build_sign_chart(expr, variable, roots, multiplicities);

    return select_intervals(chart, type, roots, multiplicities);
}

IntervalUnion InequalitySolver::solve_rational_inequality(
    const std::shared_ptr<SymbolicExpr>& numerator,
    const std::shared_ptr<SymbolicExpr>& denominator,
    InequalityType type,
    const std::string& variable) {

    if (!numerator || !denominator) return IntervalUnion::empty();

    auto num_poly = symbolic_to_poly<SymbolicPolyCoeff>(numerator, variable);
    auto den_poly = symbolic_to_poly<SymbolicPolyCoeff>(denominator, variable);

    if (den_poly.is_zero()) return IntervalUnion::empty();

    if (num_poly.is_zero()) {
        if (type == InequalityType::GreaterThan || type == InequalityType::LessThan) {
            return IntervalUnion::empty();
        }

        auto den_roots_with_mult = find_roots_with_multiplicity(denominator, variable);
        if (den_roots_with_mult.empty()) {
            return IntervalUnion::entire_line();
        }

        std::vector<std::shared_ptr<SymbolicExpr>> den_roots;
        for (const auto& [root, mult] : den_roots_with_mult) {
            den_roots.push_back(root);
        }
        std::sort(den_roots.begin(), den_roots.end(), root_less_than);

        std::vector<Interval> intervals;

        {
            Interval iv;
            iv.lower = Endpoint::neg_inf();
            iv.upper = Endpoint::open(den_roots[0]);
            intervals.push_back(iv);
        }
        for (size_t i = 0; i + 1 < den_roots.size(); ++i) {
            Interval iv;
            iv.lower = Endpoint::open(den_roots[i]);
            iv.upper = Endpoint::open(den_roots[i + 1]);
            intervals.push_back(iv);
        }

        {
            Interval iv;
            iv.lower = Endpoint::open(den_roots.back());
            iv.upper = Endpoint::pos_inf();
            intervals.push_back(iv);
        }
        return IntervalUnion(intervals);
    }

    {
        auto num_poly_rat = symbolic_to_poly<Rational>(numerator, variable);
        if (num_poly_rat.is_zero() && depends_on_var(numerator->root, variable)) {
            return IntervalUnion::empty();
        }

    }

    {
        auto den_poly_rat = symbolic_to_poly<Rational>(denominator, variable);
        if (den_poly_rat.is_zero() && depends_on_var(denominator->root, variable)) {
            return IntervalUnion::empty();
        }

    }

    auto num_roots_with_mult = find_roots_with_multiplicity(numerator, variable);

    auto den_roots_with_mult = find_roots_with_multiplicity(denominator, variable);

    struct CriticalPoint {
        std::shared_ptr<SymbolicExpr> value;
        int num_multiplicity;
        int den_multiplicity;
    };

    std::vector<CriticalPoint> critical_points;

    for (const auto& [root, mult] : num_roots_with_mult) {
        critical_points.push_back({root, mult, 0});
    }
    for (const auto& [root, mult] : den_roots_with_mult) {

        bool found = false;
        for (auto& cp : critical_points) {
            if (roots_equal(cp.value, root)) {
                cp.den_multiplicity = mult;
                found = true;
                break;
            }
        }
        if (!found) {
            critical_points.push_back({root, 0, mult});
        }
    }

    std::sort(critical_points.begin(), critical_points.end(),
        [](const CriticalPoint& a, const CriticalPoint& b) {
            return root_less_than(a.value, b.value);
        });

    int num_leading_sign = determine_leading_sign(num_poly);
    int den_leading_sign = determine_leading_sign(den_poly);
    int combined_leading_sign = num_leading_sign * den_leading_sign;

    size_t n_cp = critical_points.size();
    std::vector<int> interval_signs(n_cp + 1);

    interval_signs[n_cp] = combined_leading_sign;

    for (int i = (int)n_cp - 1; i >= 0; --i) {
        int total_mult = critical_points[i].num_multiplicity + critical_points[i].den_multiplicity;
        interval_signs[i] = interval_signs[i + 1];
        if (total_mult % 2 != 0) {
            interval_signs[i] = -interval_signs[i];
        }
    }

    bool want_positive = (type == InequalityType::GreaterThan || type == InequalityType::GreaterEqual);
    bool is_strict = (type == InequalityType::GreaterThan || type == InequalityType::LessThan);
    int target_sign = want_positive ? 1 : -1;

    std::vector<Interval> result_intervals;

    if (n_cp == 0) {

        if (interval_signs[0] == target_sign) {
            result_intervals.push_back(Interval::entire_line());
        }
    } else {

        if (interval_signs[0] == target_sign) {
            Interval iv;
            iv.lower = Endpoint::neg_inf();
            iv.upper = Endpoint::open(critical_points[0].value);
            result_intervals.push_back(iv);
        }

        for (size_t i = 0; i + 1 < n_cp; ++i) {
            if (interval_signs[i + 1] == target_sign) {
                Interval iv;
                iv.lower = Endpoint::open(critical_points[i].value);
                iv.upper = Endpoint::open(critical_points[i + 1].value);
                result_intervals.push_back(iv);
            }
        }

        if (interval_signs[n_cp] == target_sign) {
            Interval iv;
            iv.lower = Endpoint::open(critical_points[n_cp - 1].value);
            iv.upper = Endpoint::pos_inf();
            result_intervals.push_back(iv);
        }
    }

    if (!is_strict) {
        for (const auto& cp : critical_points) {

            if (cp.num_multiplicity > 0 && cp.den_multiplicity == 0) {

                bool merged = false;
                for (auto& iv : result_intervals) {

                    if (!iv.upper.is_pos_infinity && iv.upper.value) {
                        if (roots_equal(iv.upper.value, cp.value)) {
                            iv.upper.is_open = false;
                            merged = true;
                        }
                    }

                    if (!iv.lower.is_neg_infinity && iv.lower.value) {
                        if (roots_equal(iv.lower.value, cp.value)) {
                            iv.lower.is_open = false;
                            merged = true;
                        }
                    }
                }

                if (!merged) {
                    result_intervals.push_back(Interval::point(cp.value));
                }
            }
        }
    }

    return IntervalUnion(result_intervals);
}

IntervalUnion InequalitySolver::solve_inequalities(
    const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                 InequalityType>>& inequalities,
    const std::string& variable) {

    if (inequalities.empty()) return IntervalUnion::entire_line();

    IntervalUnion result = solve_inequality(inequalities[0].first, inequalities[0].second, variable);

    for (size_t i = 1; i < inequalities.size(); ++i) {
        auto solution = solve_inequality(inequalities[i].first, inequalities[i].second, variable);
        result = result.intersect(solution);
        if (result.is_empty()) break;
    }

    return result;
}

IntervalUnion InequalitySolver::build_parametric_solution(
    const std::vector<std::shared_ptr<SymbolicExpr>>& symbolic_roots,
    const std::vector<int>& multiplicities,
    int leading_sign,
    InequalityType type) {

    if (symbolic_roots.empty()) {

        bool want_positive = (type == InequalityType::GreaterThan || type == InequalityType::GreaterEqual);
        int target_sign = want_positive ? 1 : -1;
        if (leading_sign == target_sign) {
            return IntervalUnion::entire_line();
        }
        return IntervalUnion::empty();
    }

    size_t n = symbolic_roots.size();
    std::vector<int> interval_signs(n + 1);

    interval_signs[n] = leading_sign;

    for (int i = (int)n - 1; i >= 0; --i) {
        interval_signs[i] = interval_signs[i + 1];
        if (multiplicities[i] % 2 != 0) {
            interval_signs[i] = -interval_signs[i];
        }
    }

    bool want_positive = (type == InequalityType::GreaterThan || type == InequalityType::GreaterEqual);
    bool is_strict = (type == InequalityType::GreaterThan || type == InequalityType::LessThan);
    int target_sign = want_positive ? 1 : -1;

    std::vector<Interval> result_intervals;

    if (interval_signs[0] == target_sign) {
        Interval iv;
        iv.lower = Endpoint::neg_inf();
        iv.upper = Endpoint::open(symbolic_roots[0]);
        result_intervals.push_back(iv);
    }

    for (size_t i = 0; i + 1 < n; ++i) {
        if (interval_signs[i + 1] == target_sign) {
            Interval iv;
            iv.lower = Endpoint::open(symbolic_roots[i]);
            iv.upper = Endpoint::open(symbolic_roots[i + 1]);
            result_intervals.push_back(iv);
        }
    }

    if (interval_signs[n] == target_sign) {
        Interval iv;
        iv.lower = Endpoint::open(symbolic_roots[n - 1]);
        iv.upper = Endpoint::pos_inf();
        result_intervals.push_back(iv);
    }

    if (!is_strict) {
        for (size_t i = 0; i < symbolic_roots.size(); ++i) {
            bool merged = false;
            for (auto& iv : result_intervals) {

                if (!iv.upper.is_pos_infinity && iv.upper.value) {

                    auto diff = SymbolicExpr::add(iv.upper.value,
                        SymbolicExpr::multiply(symbolic_roots[i], SymbolicExpr::number(-1)));
                    if (diff->simplify()->is_zero()) {
                        iv.upper.is_open = false;
                        merged = true;
                    }
                }

                if (!iv.lower.is_neg_infinity && iv.lower.value) {
                    auto diff = SymbolicExpr::add(iv.lower.value,
                        SymbolicExpr::multiply(symbolic_roots[i], SymbolicExpr::number(-1)));
                    if (diff->simplify()->is_zero()) {
                        iv.lower.is_open = false;
                        merged = true;
                    }
                }
            }

            if (!merged) {
                result_intervals.push_back(Interval::point(symbolic_roots[i]));
            }
        }
    }

    return IntervalUnion(result_intervals);
}

static bool depends_on_any_param(const std::shared_ptr<SymbolicExpr>& expr,
                                  const std::vector<std::string>& parameters) {
    if (!expr || !expr->root) return false;
    for (const auto& param : parameters) {
        if (depends_on_var(expr->root, param)) return true;
    }
    return false;
}

static std::vector<std::shared_ptr<SymbolicExpr>> solve_symbolic_poly(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& variable) {

    if (poly.is_zero() || poly.degree() < 1) return {};

    int deg = poly.degree();
    auto get_coeff = [&](int d) -> std::shared_ptr<SymbolicExpr> {
        if (d < 0 || d > deg) return SymbolicExpr::number(0);
        return poly.coeffs[d].val ? poly.coeffs[d].val : SymbolicExpr::number(0);
    };

    if (deg == 1) {

        auto a = get_coeff(1);
        auto b = get_coeff(0);
        auto neg_b = SymbolicExpr::multiply(b, SymbolicExpr::number(-1));
        auto root = SymbolicExpr::divide(neg_b, a)->simplify();
        return { root };
    }

    auto results = solve_by_factoring(poly, variable);
    return results;
}

PiecewiseIntervalResult InequalitySolver::solve_parametric_inequality(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable,
    const std::vector<std::string>& parameters) {

    PiecewiseIntervalResult result;

    if (!expr) return result;

    if (parameters.empty()) {
        auto solution = solve_inequality(expr, type, variable);
        PiecewiseIntervalResult::Case single_case;
        single_case.condition = nullptr;
        single_case.solution = solution;
        result.cases.push_back(single_case);
        return result;
    }

    auto poly = symbolic_to_poly<SymbolicPolyCoeff>(expr, variable);

    if (poly.is_zero()) {

        if (depends_on_var(expr->root, variable)) {

            return result;
        }

        PiecewiseIntervalResult::Case zero_case;
        zero_case.condition = nullptr;
        if (type == InequalityType::GreaterEqual || type == InequalityType::LessEqual) {
            zero_case.solution = IntervalUnion::entire_line();
        } else {
            zero_case.solution = IntervalUnion::empty();
        }
        result.cases.push_back(zero_case);
        return result;
    }

    int deg = poly.degree();
    auto leading_coeff = poly.coeffs[deg].val;
    if (!leading_coeff) leading_coeff = SymbolicExpr::number(0);
    leading_coeff = leading_coeff->simplify();

    bool lc_depends_on_params = depends_on_any_param(leading_coeff, parameters);

    if (!lc_depends_on_params) {

        int leading_sign = 1;
        try {
            double val = leading_coeff->to_numeric();
            leading_sign = (val > 0) ? 1 : -1;
        } catch (...) {

            if (auto num = std::dynamic_pointer_cast<NumberNode>(leading_coeff->root)) {
                if (std::holds_alternative<BigInt>(num->value)) {
                    leading_sign = std::get<BigInt>(num->value).IsNegative() ? -1 : 1;
                } else if (std::holds_alternative<Rational>(num->value)) {
                    leading_sign = std::get<Rational>(num->value).get_numerator().IsNegative() ? -1 : 1;
                } else if (std::holds_alternative<lmmc_real_t>(num->value)) {
                    leading_sign = std::get<lmmc_real_t>(num->value) < 0 ? -1 : 1;
                }
            }
        }

        auto symbolic_roots = solve_symbolic_poly(poly, variable);

        // If the polynomial has degree >= 1 but no roots could be obtained
        // symbolically (typical when coefficients other than the leading one
        // depend on parameters and no factoring is possible), we cannot infer
        // the sign chart safely. Signal failure.
        if (symbolic_roots.empty() && poly.degree() >= 1) {
            return PiecewiseIntervalResult{};
        }

        // solve_closed_form_poly for quadratics returns roots in the order
        // (-b+√D)/(2a), (-b-√D)/(2a). When a > 0 the first root is larger;
        // build_parametric_solution expects ascending order. Swap if needed.
        if (symbolic_roots.size() == 2 && leading_sign > 0) {
            std::swap(symbolic_roots[0], symbolic_roots[1]);
        }

        std::vector<int> multiplicities(symbolic_roots.size(), 1);

        auto solution = build_parametric_solution(symbolic_roots, multiplicities, leading_sign, type);

        PiecewiseIntervalResult::Case single_case;
        single_case.condition = nullptr;
        single_case.solution = solution;
        result.cases.push_back(single_case);
    } else {

        {
            PiecewiseIntervalResult::Case pos_case;
            pos_case.condition = std::make_shared<SymbolicExpr>(
                std::make_shared<RelationalNode>(
                    leading_coeff->root,
                    SymbolicExpr::number(0)->root,
                    RelationalNode::Op::GT));

            auto symbolic_roots = solve_symbolic_poly(poly, variable);
            // 与上面 deg>=1 时的保护一致：若拿不到符号根，无法安全构造区间表达式。
            if (symbolic_roots.empty() && poly.degree() >= 1) {
                return PiecewiseIntervalResult{};
            }
            std::vector<int> multiplicities(symbolic_roots.size(), 1);
            pos_case.solution = build_parametric_solution(symbolic_roots, multiplicities, 1, type);
            result.cases.push_back(pos_case);
        }

        {
            PiecewiseIntervalResult::Case neg_case;
            neg_case.condition = std::make_shared<SymbolicExpr>(
                std::make_shared<RelationalNode>(
                    leading_coeff->root,
                    SymbolicExpr::number(0)->root,
                    RelationalNode::Op::LT));

            auto symbolic_roots = solve_symbolic_poly(poly, variable);
            if (symbolic_roots.empty() && poly.degree() >= 1) {
                return PiecewiseIntervalResult{};
            }
            std::vector<int> multiplicities(symbolic_roots.size(), 1);
            neg_case.solution = build_parametric_solution(symbolic_roots, multiplicities, -1, type);
            result.cases.push_back(neg_case);
        }

        {
            PiecewiseIntervalResult::Case degen_case;
            degen_case.condition = std::make_shared<SymbolicExpr>(
                std::make_shared<RelationalNode>(
                    leading_coeff->root,
                    SymbolicExpr::number(0)->root,
                    RelationalNode::Op::EQ));

            // 当降阶后的多项式仍然按参数分情形时，degen_case 不再是单分支：
            // 此时 expanded_into_subcases 为 true，对应分支已直接 push 进 result.cases。
            bool expanded_into_subcases = false;

            if (deg >= 1) {
                std::vector<SymbolicPolyCoeff> reduced_coeffs;
                for (int i = 0; i < deg; ++i) {
                    reduced_coeffs.push_back(poly.coeffs[i]);
                }
                Polynomial<SymbolicPolyCoeff> reduced_poly(reduced_coeffs, variable);

                if (reduced_poly.is_zero()) {

                    if (type == InequalityType::GreaterEqual || type == InequalityType::LessEqual) {
                        degen_case.solution = IntervalUnion::entire_line();
                    } else {
                        degen_case.solution = IntervalUnion::empty();
                    }
                } else {

                    auto new_lc = reduced_poly.lead_coeff().val;
                    if (new_lc) new_lc = new_lc->simplify();

                    if (new_lc && depends_on_any_param(new_lc, parameters)) {

                        auto reduced_expr = SymbolicExpr::number(0);
                        auto var_expr = SymbolicExpr::variable(variable);
                        for (int i = reduced_poly.degree(); i >= 0; --i) {
                            auto coeff_val = reduced_poly.coeffs[i].val;
                            if (!coeff_val) continue;
                            if (i == 0) {
                                reduced_expr = SymbolicExpr::add(reduced_expr, coeff_val);
                            } else if (i == 1) {
                                reduced_expr = SymbolicExpr::add(reduced_expr,
                                    SymbolicExpr::multiply(coeff_val, var_expr));
                            } else {
                                reduced_expr = SymbolicExpr::add(reduced_expr,
                                    SymbolicExpr::multiply(coeff_val,
                                        SymbolicExpr::power(var_expr, SymbolicExpr::number(i))));
                            }
                        }
                        reduced_expr = reduced_expr->simplify();
                        auto sub_result = solve_parametric_inequality(reduced_expr, type, variable, parameters);

                        // 不能只保留 sub_result.cases[0]，否则降阶后仍按参数分情形的解会被静默丢弃。
                        // 把每个子分支的参数条件与父级 degen_case 的条件（leading coeff == 0）合取，
                        // 作为独立分支加入结果。
                        if (sub_result.cases.empty()) {
                            degen_case.solution = IntervalUnion::empty();
                        } else {
                            for (const auto& sub_case : sub_result.cases) {
                                PiecewiseIntervalResult::Case merged;
                                if (sub_case.condition) {
                                    merged.condition = std::make_shared<SymbolicExpr>(
                                        std::make_shared<LogicalNode>(
                                            degen_case.condition->root,
                                            sub_case.condition->root,
                                            LogicalNode::Op::And));
                                } else {
                                    merged.condition = degen_case.condition;
                                }
                                merged.solution = sub_case.solution;
                                result.cases.push_back(merged);
                            }
                            expanded_into_subcases = true;
                        }
                    } else {

                        int reduced_leading_sign = 1;
                        if (new_lc) {
                            try {
                                double val = new_lc->to_numeric();
                                reduced_leading_sign = (val > 0) ? 1 : -1;
                            } catch (...) {
                                if (auto num = std::dynamic_pointer_cast<NumberNode>(new_lc->root)) {
                                    if (std::holds_alternative<BigInt>(num->value)) {
                                        reduced_leading_sign = std::get<BigInt>(num->value).IsNegative() ? -1 : 1;
                                    } else if (std::holds_alternative<Rational>(num->value)) {
                                        reduced_leading_sign = std::get<Rational>(num->value).get_numerator().IsNegative() ? -1 : 1;
                                    } else if (std::holds_alternative<lmmc_real_t>(num->value)) {
                                        reduced_leading_sign = std::get<lmmc_real_t>(num->value) < 0 ? -1 : 1;
                                    }
                                }
                            }
                        }

                        auto symbolic_roots = solve_symbolic_poly(reduced_poly, variable);
                        std::vector<int> multiplicities(symbolic_roots.size(), 1);
                        degen_case.solution = build_parametric_solution(
                            symbolic_roots, multiplicities, reduced_leading_sign, type);
                    }
                }
            } else {
                degen_case.solution = IntervalUnion::empty();
            }

            if (!expanded_into_subcases) {
                result.cases.push_back(degen_case);
            }
        }
    }

    return result;
}

}
