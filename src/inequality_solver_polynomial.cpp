#include "inequality_solver.hpp"
#include "symbolic_ast.hpp"
#include "poly_utils.hpp"
#include "solve_polynomial.hpp"
#include "solve_strategies.hpp"
#include "newton_raphson.hpp"
#include "root_of_utils.hpp"
#include "internal/exact_algebraic.hpp"
#include "internal/inequality_solver_support.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace LMCAS {
namespace {

constexpr const char* kCheckedInequalityOperation = "solve_inequality_checked";

bool rational_is_negative(const Rational& value) {
    return value.get_numerator().IsNegative();
}

bool constant_satisfies(const Rational& value, InequalityType type) {
    const bool zero = value == Rational(0);
    const bool negative = rational_is_negative(value);
    switch (type) {
        case InequalityType::GreaterThan: return !negative && !zero;
        case InequalityType::GreaterEqual: return !negative;
        case InequalityType::LessThan: return negative;
        case InequalityType::LessEqual: return negative || zero;
    }
    return false;
}

Endpoint finite_endpoint(const std::shared_ptr<SymbolicExpr>& value, bool strict) {
    return strict ? Endpoint::open(value) : Endpoint::closed(value);
}

Result<void> verify_quadratic_boundary(
    const Polynomial<Rational>& polynomial,
    const std::shared_ptr<SymbolicExpr>& boundary,
    const Rational& rational_part,
    const Rational& radical_coefficient,
    const Rational& radicand,
    ComputationContext& context) {
    auto steps = context.consume_steps(2, kCheckedInequalityOperation);
    if (!steps) return steps;

    if (!boundary || !LMCAS::detail::node(boundary) || polynomial.coeffs.size() < 3) {
        return Result<void>::failure(
            CasErrc::InternalInvariant,
            "quadratic boundary verification could not construct an expression",
            kCheckedInequalityOperation);
    }

    const Rational& constant = polynomial.coeffs[0];
    const Rational& linear = polynomial.coeffs[1];
    const Rational& quadratic = polynomial.coeffs[2];
    const Rational residual_rational =
        quadratic * (rational_part * rational_part +
                     radical_coefficient * radical_coefficient * radicand) +
        linear * rational_part + constant;
    const Rational residual_radical =
        Rational(2) * quadratic * rational_part * radical_coefficient +
        linear * radical_coefficient;
    if (residual_rational != Rational(0) ||
        residual_radical != Rational(0)) {
        return Result<void>::failure(
            CasErrc::Inconclusive,
            "quadratic boundary could not be verified by exact substitution",
            kCheckedInequalityOperation);
    }
    return Result<void>::success();
}

} // namespace

Result<IntervalUnion> InequalitySolver::solve_exact_quadratic_inequality(
    const Polynomial<Rational>& polynomial,
    InequalityType type,
    ComputationContext& context) {
    const Rational& constant = polynomial.coeffs[0];
    const Rational& linear = polynomial.coeffs[1];
    const Rational& quadratic = polynomial.coeffs[2];
    if (quadratic == Rational(0)) {
        return Result<IntervalUnion>::failure(
            CasErrc::InternalInvariant,
            "quadratic polynomial has a zero leading coefficient",
            kCheckedInequalityOperation);
    }

    auto steps = context.consume_steps(3, kCheckedInequalityOperation);
    if (!steps) return Result<IntervalUnion>::failure(steps.error());
    const Rational discriminant = linear * linear - Rational(4) * quadratic * constant;
    const bool leading_positive = !rational_is_negative(quadratic);
    const bool wants_positive = type == InequalityType::GreaterThan ||
                                type == InequalityType::GreaterEqual;
    const bool strict = type == InequalityType::GreaterThan ||
                        type == InequalityType::LessThan;

    if (rational_is_negative(discriminant)) {
        const bool satisfied = wants_positive == leading_positive;
        return Result<IntervalUnion>::success(
            satisfied ? IntervalUnion::entire_line() : IntervalUnion::empty());
    }

    const Rational denominator = Rational(2) * quadratic;
    const Rational neg_linear = Rational(0) - linear;
    const Rational rational_part = neg_linear / denominator;
    if (discriminant == Rational(0)) {
        auto root = SymbolicExpr::number(rational_part);
        auto verified = verify_quadratic_boundary(
            polynomial, root, rational_part, Rational(0), Rational(0), context);
        if (!verified) return Result<IntervalUnion>::failure(verified.error());

        if (wants_positive == leading_positive) {
            if (!strict) return Result<IntervalUnion>::success(IntervalUnion::entire_line());
            std::vector<Interval> rays{
                Interval{Endpoint::neg_inf(), Endpoint::open(root)},
                Interval{Endpoint::open(root), Endpoint::pos_inf()}
            };
            return Result<IntervalUnion>::success(
                IntervalUnion::from_checked_normalized(std::move(rays)));
        }
        if (strict) return Result<IntervalUnion>::success(IntervalUnion::empty());
        return Result<IntervalUnion>::success(
            IntervalUnion::from_checked_normalized({Interval::point(root)}));
    }

    auto nodes = context.reserve_nodes(20, kCheckedInequalityOperation);
    if (!nodes) return Result<IntervalUnion>::failure(nodes.error());
    auto sqrt_discriminant = SymbolicExpr::sqrt(
        SymbolicExpr::number(discriminant))->simplify();
    const Rational first_radical_coefficient = Rational(-1) / denominator;
    const Rational second_radical_coefficient = Rational(1) / denominator;
    auto first_formula = SymbolicExpr::add(
        SymbolicExpr::number(rational_part),
        SymbolicExpr::multiply(SymbolicExpr::number(first_radical_coefficient),
                               sqrt_discriminant))->simplify();
    auto second_formula = SymbolicExpr::add(
        SymbolicExpr::number(rational_part),
        SymbolicExpr::multiply(SymbolicExpr::number(second_radical_coefficient),
                               sqrt_discriminant))->simplify();
    auto lower = leading_positive ? first_formula : second_formula;
    auto upper = leading_positive ? second_formula : first_formula;

    auto first_verified = verify_quadratic_boundary(
        polynomial, first_formula, rational_part,
        first_radical_coefficient, discriminant, context);
    if (!first_verified) return Result<IntervalUnion>::failure(first_verified.error());
    auto second_verified = verify_quadratic_boundary(
        polynomial, second_formula, rational_part,
        second_radical_coefficient, discriminant, context);
    if (!second_verified) return Result<IntervalUnion>::failure(second_verified.error());

    const bool outside = wants_positive == leading_positive;
    if (outside) {
        std::vector<Interval> rays{
            Interval{Endpoint::neg_inf(), finite_endpoint(lower, strict)},
            Interval{finite_endpoint(upper, strict), Endpoint::pos_inf()}
        };
        return Result<IntervalUnion>::success(
            IntervalUnion::from_checked_normalized(std::move(rays)));
    }
    return Result<IntervalUnion>::success(IntervalUnion::from_checked_normalized({
        Interval{finite_endpoint(lower, strict), finite_endpoint(upper, strict)}
    }));
}
namespace {

Result<IntervalUnion> solve_exact_affine_inequality_impl(
    const Polynomial<Rational>& polynomial,
    InequalityType type,
    ComputationContext& context) {
    auto step = context.consume_steps(1, kCheckedInequalityOperation);
    if (!step) return Result<IntervalUnion>::failure(step.error());

    if (polynomial.degree() <= 0) {
        const Rational value = polynomial.is_zero()
            ? Rational(0)
            : polynomial.coeffs[0];
        return Result<IntervalUnion>::success(
            constant_satisfies(value, type)
                ? IntervalUnion::entire_line()
                : IntervalUnion::empty());
    }

    if (polynomial.degree() != 1 || polynomial.coeffs.size() < 2) {
        return Result<IntervalUnion>::failure(
            CasErrc::Inconclusive,
            "checked inequality solving currently supports exact polynomials through degree two",
            kCheckedInequalityOperation);
    }

    const Rational& intercept = polynomial.coeffs[0];
    const Rational& slope = polynomial.coeffs[1];
    if (slope == Rational(0)) {
        return Result<IntervalUnion>::failure(
            CasErrc::InternalInvariant,
            "affine polynomial has a zero leading coefficient",
            kCheckedInequalityOperation);
    }

    const Rational boundary = (Rational(0) - intercept) / slope;
    auto boundary_expr = SymbolicExpr::number(boundary);
    auto nodes = context.reserve_nodes(1, kCheckedInequalityOperation);
    if (!nodes) return Result<IntervalUnion>::failure(nodes.error());

    const bool wants_greater = type == InequalityType::GreaterThan ||
                               type == InequalityType::GreaterEqual;
    const bool strict = type == InequalityType::GreaterThan ||
                        type == InequalityType::LessThan;
    const bool solution_above = wants_greater != rational_is_negative(slope);
    const Endpoint finite = strict
        ? Endpoint::open(boundary_expr)
        : Endpoint::closed(boundary_expr);
    const Interval interval = solution_above
        ? Interval{finite, Endpoint::pos_inf()}
        : Interval{Endpoint::neg_inf(), finite};
    return Result<IntervalUnion>::success(IntervalUnion::from_single(interval));
}


struct ExactInequalityRoot {
    std::shared_ptr<SymbolicExpr> expression;
    LMCAS::detail::ExactRealAlgebraic algebraic;
    int multiplicity = 1;
};

std::shared_ptr<SymbolicExpr> rational_polynomial_expression(
    const Polynomial<Rational>& polynomial) {
    auto variable = SymbolicExpr::variable(polynomial.variable_name);
    std::shared_ptr<SymbolicExpr> sum = SymbolicExpr::number(0);
    for (int degree = 0; degree <= polynomial.degree(); ++degree) {
        if (polynomial.coeffs[degree] == Rational(0)) continue;
        auto term = SymbolicExpr::number(polynomial.coeffs[degree]);
        if (degree > 0) {
            auto power = degree == 1
                ? variable
                : SymbolicExpr::power(variable, SymbolicExpr::number(degree));
            term = SymbolicExpr::multiply(term, power);
        }
        sum = SymbolicExpr::add(sum, term);
    }
    return sum->simplify();
}

Result<IntervalUnion> solve_exact_polynomial_inequality_impl(
    const Polynomial<Rational>& polynomial,
    InequalityType type,
    ComputationContext& context) {
    auto factors = square_free_factorization(polynomial);
    std::vector<ExactInequalityRoot> roots;
    try {
        for (const auto& [factor, multiplicity] : factors) {
            if (factor.degree() < 1) continue;
            auto isolated = isolate_real_roots_checked(factor, context);
            if (!isolated) {
                return Result<IntervalUnion>::failure(isolated.error());
            }
            const auto& intervals = isolated.value();
            auto factor_expression = rational_polynomial_expression(factor);
            for (std::size_t index = 0; index < intervals.size(); ++index) {
                auto algebraic = LMCAS::detail::make_exact_real_algebraic(
                    factor, index, static_cast<std::size_t>(multiplicity), context);
                if (!algebraic) {
                    return Result<IntervalUnion>::failure(algebraic.error());
                }
                roots.push_back(ExactInequalityRoot{
                    SymbolicExpr::root_of(
                        factor_expression, polynomial.variable_name,
                        static_cast<int>(index)),
                    std::move(algebraic.value()),
                    multiplicity});
            }
        }
        std::sort(roots.begin(), roots.end(),
            [](const ExactInequalityRoot& lhs, const ExactInequalityRoot& rhs) {
                return lhs.algebraic.lower < rhs.algebraic.lower;
            });
    } catch (const std::bad_alloc&) {
        return Result<IntervalUnion>::failure(
            CasErrc::ResourceLimit,
            "polynomial inequality root allocation failed",
            kCheckedInequalityOperation);
    } catch (const std::exception& error) {
        return Result<IntervalUnion>::failure(
            CasErrc::InternalInvariant, error.what(),
            kCheckedInequalityOperation);
    }

    const bool wants_positive =
        type == InequalityType::GreaterThan ||
        type == InequalityType::GreaterEqual;
    const bool strict =
        type == InequalityType::GreaterThan ||
        type == InequalityType::LessThan;
    auto sign_satisfies = [&](int sign_value) {
        return wants_positive ? sign_value > 0 : sign_value < 0;
    };

    int sign_value = rational_is_negative(polynomial.lead_coeff()) ? -1 : 1;
    if ((polynomial.degree() & 1) != 0) sign_value = -sign_value;

    std::vector<Interval> intervals;
    Endpoint lower = Endpoint::neg_inf();
    for (const auto& root : roots) {
        if (sign_satisfies(sign_value)) {
            intervals.push_back(Interval{
                lower, Endpoint::open(root.expression)});
        }
        if (!strict) {
            intervals.push_back(Interval::point(root.expression));
        }
        lower = Endpoint::open(root.expression);
        if ((root.multiplicity & 1) != 0) sign_value = -sign_value;
    }
    if (sign_satisfies(sign_value)) {
        intervals.push_back(Interval{lower, Endpoint::pos_inf()});
    }
    return IntervalUnion::from_intervals_checked(
        std::move(intervals), context);
}
} // namespace

namespace detail::inequality_support {

Result<IntervalUnion> solve_exact_affine_inequality(
    const Polynomial<Rational>& polynomial,
    InequalityType type,
    ComputationContext& context) {
    return solve_exact_affine_inequality_impl(polynomial, type, context);
}

Result<IntervalUnion> solve_exact_polynomial_inequality(
    const Polynomial<Rational>& polynomial,
    InequalityType type,
    ComputationContext& context) {
    return solve_exact_polynomial_inequality_impl(polynomial, type, context);
}

} // namespace detail::inequality_support

using detail::inequality_support::determine_leading_sign;
using detail::inequality_support::find_roots_with_multiplicity;
using detail::inequality_support::root_less_than;
using detail::inequality_support::roots_equal;

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
    const std::vector<int>&) {

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
} // namespace LMCAS
