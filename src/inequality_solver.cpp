#include "../include/inequality_solver.hpp"
#include "symbolic_ast.hpp"
#include "../include/numeric_evaluation.hpp"
#include "../include/poly_utils.hpp"
#include "poly_utils_internal.hpp"
#include "../include/solve_polynomial.hpp"
#include "../include/solve_strategies.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <limits>
#include <functional>
#include <optional>

namespace lamina {

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

    if (!boundary || !lamina::detail::node(boundary) || polynomial.coeffs.size() < 3) {
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

Result<IntervalUnion> solve_exact_affine_inequality(
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

} // namespace

Result<IntervalUnion> InequalitySolver::solve_inequality_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable,
    ComputationContext& context) {
    if (!expr || !lamina::detail::node(expr)) {
        return Result<IntervalUnion>::failure(
            CasErrc::InvalidArgument, "inequality expression cannot be null",
            kCheckedInequalityOperation);
    }
    if (variable.empty()) {
        return Result<IntervalUnion>::failure(
            CasErrc::InvalidArgument, "inequality variable cannot be empty",
            kCheckedInequalityOperation);
    }

    try {
        auto recognized = recognize_rational_polynomial(*expr, variable, context);
        if (!recognized) return Result<IntervalUnion>::failure(recognized.error());
        if (!recognized.value()) {
            return Result<IntervalUnion>::failure(
                CasErrc::Inconclusive,
                "expression is not an exact rational polynomial in the requested variable",
                kCheckedInequalityOperation);
        }
        if (recognized.value()->degree() == 2 &&
            recognized.value()->coeffs.size() >= 3) {
            return solve_exact_quadratic_inequality(
                *recognized.value(), type, context);
        }
        return solve_exact_affine_inequality(*recognized.value(), type, context);
    } catch (const std::bad_alloc&) {
        return Result<IntervalUnion>::failure(
            CasErrc::ResourceLimit, "inequality allocation failed",
            kCheckedInequalityOperation);
    } catch (const std::exception& error) {
        return Result<IntervalUnion>::failure(
            CasErrc::InternalInvariant, error.what(), kCheckedInequalityOperation);
    }
}

Result<IntervalUnion> InequalitySolver::solve_inequality_checked(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable) {
    ComputationContext context;
    return solve_inequality_checked(expr, type, variable, context);
}

Result<IntervalUnion> InequalitySolver::solve_inequalities_checked(
    const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                 InequalityType>>& inequalities,
    const std::string& variable,
    ComputationContext& context) {
    if (variable.empty()) {
        return Result<IntervalUnion>::failure(
            CasErrc::InvalidArgument, "inequality variable cannot be empty",
            "solve_inequalities_checked");
    }

    auto initial_step = context.consume_steps(1, "solve_inequalities_checked");
    if (!initial_step) return Result<IntervalUnion>::failure(initial_step.error());
    if (inequalities.empty()) {
        return Result<IntervalUnion>::success(IntervalUnion::entire_line());
    }

    std::optional<IntervalUnion> aggregate;
    for (const auto& inequality : inequalities) {
        auto solved = solve_inequality_checked(
            inequality.first, inequality.second, variable, context);
        if (!solved) return Result<IntervalUnion>::failure(solved.error());
        if (!aggregate) {
            aggregate = std::move(solved.value());
            continue;
        }
        auto intersection = aggregate->intersect_checked(solved.value(), context);
        if (!intersection) {
            CasError error = intersection.error();
            error.operation = "solve_inequalities_checked";
            return Result<IntervalUnion>::failure(std::move(error));
        }
        aggregate = std::move(intersection.value());
        if (aggregate->is_empty()) break;
    }
    return Result<IntervalUnion>::success(std::move(*aggregate));
}

Result<IntervalUnion> InequalitySolver::solve_inequalities_checked(
    const std::vector<std::pair<std::shared_ptr<SymbolicExpr>,
                                 InequalityType>>& inequalities,
    const std::string& variable) {
    ComputationContext context;
    return solve_inequalities_checked(inequalities, variable, context);
}

static std::optional<double> try_checked_numeric_constant(const SymbolicExpr& expr) {
    ComputationContext context;
    auto evaluated = evaluate_numeric(expr, NumericBindings{}, context);
    if (!evaluated || !evaluated.value().is_finite() ||
        !std::isfinite(evaluated.value().value)) {
        return std::nullopt;
    }
    return evaluated.value().value;
}

static int exact_numeric_sign(const std::shared_ptr<SymbolicExpr>& expr) {
    if (!expr || !lamina::detail::node(expr)) return 0;
    auto simplified = expr->simplify();
    if (!simplified || !lamina::detail::node(simplified)) return 0;
    auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(simplified));
    if (!num) return 0;
    if (std::holds_alternative<BigInt>(num->value())) {
        const auto& value = std::get<BigInt>(num->value());
        if (value.is_zero()) return 0;
        return value.IsNegative() ? -1 : 1;
    }
    if (std::holds_alternative<Rational>(num->value())) {
        const auto& value = std::get<Rational>(num->value());
        if (value.get_numerator().is_zero()) return 0;
        return value.get_numerator().IsNegative() ? -1 : 1;
    }
    const auto value = std::get<lmmc_real_t>(num->value());
    if (!std::isfinite(value) || value == 0.0) return 0;
    return value < 0 ? -1 : 1;
}

static int determine_leading_sign(const Polynomial<SymbolicPolyCoeff>& poly) {
    if (poly.is_zero()) return 0;
    auto lc = poly.lead_coeff().val;
    if (!lc) return 1;
    auto simplified = lc->simplify();
    if (!simplified) return 1;

    if (auto val = try_checked_numeric_constant(*simplified)) {
        if (*val > 0) return 1;
        if (*val < 0) return -1;
    }

    if (auto num = std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(simplified))) {
        if (std::holds_alternative<BigInt>(num->value())) {
            return std::get<BigInt>(num->value()).IsNegative() ? -1 : 1;
        }
        if (std::holds_alternative<Rational>(num->value())) {
            return std::get<Rational>(num->value()).get_numerator().IsNegative() ? -1 : 1;
        }
        if (std::holds_alternative<lmmc_real_t>(num->value())) {
            return std::get<lmmc_real_t>(num->value()) < 0 ? -1 : 1;
        }
    }
    return 1;
}

static std::shared_ptr<SymbolicExpr> snap_verified_integer_root(
    const Polynomial<Rational>& poly,
    const std::shared_ptr<SymbolicExpr>& root) {
    if (!root || poly.is_zero()) return root;
    auto numeric = try_checked_numeric_constant(*root);
    if (!numeric || !std::isfinite(*numeric)) return root;

    double rounded = std::round(*numeric);
    if (std::abs(*numeric - rounded) > 1e-8) return root;
    if (rounded < static_cast<double>(std::numeric_limits<long long>::min()) ||
        rounded > static_cast<double>(std::numeric_limits<long long>::max())) {
        return root;
    }

    Rational candidate(BigInt(static_cast<long long>(rounded)));
    if (poly.eval(candidate) == Rational(0)) {
        return SymbolicExpr::number(candidate);
    }
    return root;
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
                auto verified_root = snap_verified_integer_root(poly_rat, root);
                if (try_checked_numeric_constant(*verified_root)) {
                    result.push_back({verified_root, 1});
                }
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
                if (try_checked_numeric_constant(*root_expr)) {
                    result.push_back({root_expr, mult});
                }
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
                auto root = snap_verified_integer_root(factor, SymbolicExpr::number(r1));
                result.push_back({root, mult});
            }
            if (std::isfinite(r2) && std::abs(r1 - r2) > 1e-10) {
                auto root = snap_verified_integer_root(factor, SymbolicExpr::number(r2));
                result.push_back({root, mult});
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
            auto verified_root = snap_verified_integer_root(factor, root);
            if (try_checked_numeric_constant(*verified_root)) {
                result.push_back({verified_root, mult});
            }
        }
    }

    return result;
}

static bool root_less_than(const std::shared_ptr<SymbolicExpr>& a,
                           const std::shared_ptr<SymbolicExpr>& b) {
    if (!a || !b) return false;
    auto va = try_checked_numeric_constant(*a);
    auto vb = try_checked_numeric_constant(*b);
    if (va && vb) {
        return *va < *vb;
    }

    /// 当根含参数无法直接求值时，尝试通过 (a - b) 的符号判断大小。
    /// 对于二次公式的两个根，差值可化简为 ±sqrt(disc)/a 的形式。
    auto diff = SymbolicExpr::add(a, SymbolicExpr::multiply(b, SymbolicExpr::number(-1)))->simplify();
    if (auto vd = try_checked_numeric_constant(*diff)) {
        return *vd < 0;
    }

    /// 尝试判断差值表达式的符号结构：
    /// 如果差值形如 k * sqrt(...) / denom，判断各因子的符号。
    /// 这覆盖了二次公式根差 = sqrt(delta)/a 的情形。
    auto try_sign_of_node = [](const std::shared_ptr<const SymbolicNode>& node) -> int {
        if (!node) return 0;

        if (auto num_node = std::dynamic_pointer_cast<const NumberNode>(node)) {
            if (std::holds_alternative<BigInt>(num_node->value())) {
                auto& v = std::get<BigInt>(num_node->value());
                if (v.is_zero()) return 0;
                return v.IsNegative() ? -1 : 1;
            }
            if (std::holds_alternative<Rational>(num_node->value())) {
                auto& v = std::get<Rational>(num_node->value());
                if (v.get_numerator().is_zero()) return 0;
                return v.get_numerator().IsNegative() ? -1 : 1;
            }
            if (std::holds_alternative<lmmc_real_t>(num_node->value())) {
                auto v = std::get<lmmc_real_t>(num_node->value());
                if (v == 0.0) return 0;
                return v < 0 ? -1 : 1;
            }
        }

        /// sqrt(...) 非负（假设参数使判别式非负）
        if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(node)) {
            if (fn->type() == FunctionNode::FuncType::Sqrt) return 1;
        }

        /// x^(1/2) 或 x^0.5 也是平方根，非负
        if (auto pw = std::dynamic_pointer_cast<const PowerNode>(node)) {
            auto exp_expr = lamina::detail::make_expression_ptr(pw->exponent());
            if (auto ev = try_checked_numeric_constant(*exp_expr)) {
                if (*ev > 0 && *ev < 1.0) {
                    /// base^(正分数) >= 0（假设 base 为判别式等非负量）
                    return 1;
                }
            }

            /// 对于整数指数，判断底数符号
            auto base_expr = lamina::detail::make_expression_ptr(pw->base());
            auto bv_checked = try_checked_numeric_constant(*base_expr);
            auto ev_checked = try_checked_numeric_constant(*exp_expr);
            if (bv_checked && ev_checked) {
                double bv = *bv_checked;
                double ev = *ev_checked;
                int ei = static_cast<int>(ev);
                if (std::abs(ev - ei) < 1e-10) {
                    if (bv > 0) return 1;
                    if (bv < 0) return (ei % 2 == 0) ? 1 : -1;
                }
            }
        }

        return 0;
    };

    /// 对乘积节点，各因子符号之积
    auto try_sign_of_expr = [&try_sign_of_node](const std::shared_ptr<SymbolicExpr>& expr) -> int {
        if (!expr || !lamina::detail::node(expr)) return 0;

        /// 直接节点
        int s = try_sign_of_node(lamina::detail::node(expr));
        if (s != 0) return s;

        /// 乘积：各因子符号之积
        if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
            int sign = 1;
            for (const auto& op : mul->operands()) {
                int os = try_sign_of_node(op);
                if (os == 0) return 0;
                sign *= os;
            }
            return sign;
        }

        return 0;
    };

    int diff_sign = try_sign_of_expr(diff);
    if (diff_sign < 0) return true;
    if (diff_sign > 0) return false;

    return false;
}

static bool roots_equal(const std::shared_ptr<SymbolicExpr>& a,
                        const std::shared_ptr<SymbolicExpr>& b) {
    if (!a || !b) return false;
    auto va = try_checked_numeric_constant(*a);
    auto vb = try_checked_numeric_constant(*b);
    return va && vb && std::abs(*va - *vb) < 1e-10;
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

IntervalUnion InequalitySolver::solve_inequality(
    const std::shared_ptr<SymbolicExpr>& expr,
    InequalityType type,
    const std::string& variable) {

    if (!expr) return IntervalUnion::empty();

    auto poly = symbolic_to_poly<SymbolicPolyCoeff>(expr, variable);
    if (poly.is_zero()) {

        if (depends_on_var(lamina::detail::node(expr), variable)) {
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
        int sign = exact_numeric_sign(lc);
        if (sign == 0 && !lc->is_zero()) {
            auto simplified = lc->simplify();
            if (!simplified || !simplified->is_zero()) {
                return IntervalUnion::empty();
            }
        }

        bool satisfies = false;
        switch (type) {
            case InequalityType::GreaterThan: satisfies = (sign > 0); break;
            case InequalityType::GreaterEqual: satisfies = (sign >= 0); break;
            case InequalityType::LessThan: satisfies = (sign < 0); break;
            case InequalityType::LessEqual: satisfies = (sign <= 0); break;
            default:
                return IntervalUnion::empty();
        }
        return satisfies ? IntervalUnion::entire_line() : IntervalUnion::empty();
    }

    {
        auto poly_rat = symbolic_to_poly<Rational>(expr, variable);
        if (poly_rat.is_zero() && depends_on_var(lamina::detail::node(expr), variable)) {

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
        if (num_poly_rat.is_zero() && depends_on_var(lamina::detail::node(numerator), variable)) {
            return IntervalUnion::empty();
        }

    }

    {
        auto den_poly_rat = symbolic_to_poly<Rational>(denominator, variable);
        if (den_poly_rat.is_zero() && depends_on_var(lamina::detail::node(denominator), variable)) {
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
    if (!expr || !lamina::detail::node(expr)) return false;
    for (const auto& param : parameters) {
        if (depends_on_var(lamina::detail::node(expr), param)) return true;
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

        if (depends_on_var(lamina::detail::node(expr), variable)) {

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

        int leading_sign = exact_numeric_sign(leading_coeff);
        if (leading_sign == 0 && !leading_coeff->is_zero()) {
            if (auto val = try_checked_numeric_constant(*leading_coeff)) {
                leading_sign = (*val > 0) ? 1 : ((*val < 0) ? -1 : 0);
            }
        }

        auto symbolic_roots = solve_symbolic_poly(poly, variable);

        // Sort roots in ascending order
        /// 当根含参数时，root_less_than 通过差值符号判断排序。
        /// 对于二次公式根 r1=(-b+sqrt(d))/(2a), r2=(-b-sqrt(d))/(2a)，
        /// 当 a>0 时 r1>r2，需要交换为 [r2, r1]。
        if (symbolic_roots.size() == 2) {
            /// 尝试判断 root[0] > root[1]，若是则交换
            bool swapped = false;
            auto d = SymbolicExpr::add(symbolic_roots[0],
                SymbolicExpr::multiply(symbolic_roots[1], SymbolicExpr::number(-1)))->simplify();
            if (!d) {
                std::sort(symbolic_roots.begin(), symbolic_roots.end(), root_less_than);
                return;
            }
            /// 如果差值可以求值为正数，说明 root[0] > root[1]，需要交换
            if (auto dv = try_checked_numeric_constant(*d)) {
                if (*dv > 0) { std::swap(symbolic_roots[0], symbolic_roots[1]); swapped = true; }
            } else {
                /// 尝试结构化符号判断
                /// diff 为 (disc)^0.5 形式（PowerNode with exp=0.5）或含 sqrt 的乘积
                auto check_positive = [](const std::shared_ptr<SymbolicExpr>& e) -> bool {
                    if (!e || !lamina::detail::node(e)) return false;
                    /// PowerNode with exponent in (0,1) -> non-negative
                    if (auto pw = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(e))) {
                        auto exp_e = lamina::detail::make_expression_ptr(pw->exponent());
                        if (auto ev = try_checked_numeric_constant(*exp_e)) {
                            if (*ev > 0 && *ev < 1.0) return true;
                        }
                    }
                    /// FunctionNode::Sqrt -> non-negative
                    if (auto fn = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(e))) {
                        if (fn->type() == FunctionNode::FuncType::Sqrt) return true;
                    }
                    /// MultiplyNode: all factors positive
                    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(e))) {
                        int sign = 1;
                        for (const auto& op : mul->operands()) {
                            if (auto n = std::dynamic_pointer_cast<const NumberNode>(op)) {
                                if (std::holds_alternative<BigInt>(n->value())) {
                                    if (std::get<BigInt>(n->value()).IsNegative()) sign *= -1;
                                    else if (std::get<BigInt>(n->value()).is_zero()) return false;
                                } else if (std::holds_alternative<Rational>(n->value())) {
                                    if (std::get<Rational>(n->value()).get_numerator().IsNegative()) sign *= -1;
                                    else if (std::get<Rational>(n->value()).get_numerator().is_zero()) return false;
                                } else if (std::holds_alternative<lmmc_real_t>(n->value())) {
                                    if (std::get<lmmc_real_t>(n->value()) < 0) sign *= -1;
                                    else if (std::get<lmmc_real_t>(n->value()) == 0) return false;
                                }
                            } else if (auto pw2 = std::dynamic_pointer_cast<const PowerNode>(op)) {
                                auto exp_e2 = lamina::detail::make_expression_ptr(pw2->exponent());
                                auto ev2 = try_checked_numeric_constant(*exp_e2);
                                if (ev2 && *ev2 > 0 && *ev2 < 1.0) { /* positive, sign unchanged */ }
                                else return false; // can't determine
                            } else if (auto fn2 = std::dynamic_pointer_cast<const FunctionNode>(op)) {
                                if (fn2->type() == FunctionNode::FuncType::Sqrt) { /* positive */ }
                                else return false;
                            } else {
                                return false; // unknown factor
                            }
                        }
                        return sign > 0;
                    }
                    return false;
                };
                auto check_negative = [](const std::shared_ptr<SymbolicExpr>& e) -> bool {
                    if (!e || !lamina::detail::node(e)) return false;
                    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(e))) {
                        int sign = 1;
                        for (const auto& op : mul->operands()) {
                            if (auto n = std::dynamic_pointer_cast<const NumberNode>(op)) {
                                if (std::holds_alternative<BigInt>(n->value())) {
                                    if (std::get<BigInt>(n->value()).IsNegative()) sign *= -1;
                                    else if (std::get<BigInt>(n->value()).is_zero()) return false;
                                } else if (std::holds_alternative<Rational>(n->value())) {
                                    if (std::get<Rational>(n->value()).get_numerator().IsNegative()) sign *= -1;
                                    else if (std::get<Rational>(n->value()).get_numerator().is_zero()) return false;
                                } else if (std::holds_alternative<lmmc_real_t>(n->value())) {
                                    if (std::get<lmmc_real_t>(n->value()) < 0) sign *= -1;
                                    else if (std::get<lmmc_real_t>(n->value()) == 0) return false;
                                }
                            } else if (auto pw2 = std::dynamic_pointer_cast<const PowerNode>(op)) {
                                auto exp_e2 = lamina::detail::make_expression_ptr(pw2->exponent());
                                auto ev2 = try_checked_numeric_constant(*exp_e2);
                                if (ev2 && *ev2 > 0 && *ev2 < 1.0) { /* positive */ }
                                else return false;
                            } else if (auto fn2 = std::dynamic_pointer_cast<const FunctionNode>(op)) {
                                if (fn2->type() == FunctionNode::FuncType::Sqrt) { /* positive */ }
                                else return false;
                            } else {
                                return false;
                            }
                        }
                        return sign < 0;
                    }
                    return false;
                };
                if (check_positive(d)) {
                    std::swap(symbolic_roots[0], symbolic_roots[1]);
                    swapped = true;
                } else if (check_negative(d)) {
                    /// already in correct order
                    swapped = false;
                }
            }
            if (!swapped) {
                /// Fallback: 对于 a>0 的二次多项式，solve_quadratic_internal 返回
                /// [大根, 小根]，需要交换。对于 a<0 则已经是 [小根, 大根]。
                /// 这里利用 leading_sign 直接判断。
                if (leading_sign > 0 && deg == 2) {
                    std::swap(symbolic_roots[0], symbolic_roots[1]);
                }
            }
        } else {
            std::sort(symbolic_roots.begin(), symbolic_roots.end(), root_less_than);
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
            pos_case.condition = lamina::detail::make_expression_ptr(
                lamina::detail::make_node<RelationalNode>(
                    lamina::detail::node(leading_coeff),
                    lamina::detail::node(SymbolicExpr::number(0)),
                    RelationalNode::Op::GT));

            auto symbolic_roots = solve_symbolic_poly(poly, variable);
            std::vector<int> multiplicities(symbolic_roots.size(), 1);

            // Sort roots in ascending order and apply same permutation to multiplicities
            std::vector<size_t> indices(symbolic_roots.size());
            for (size_t i = 0; i < indices.size(); ++i) indices[i] = i;
            std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
                return root_less_than(symbolic_roots[a], symbolic_roots[b]);
            });
            std::vector<std::shared_ptr<SymbolicExpr>> sorted_roots;
            std::vector<int> sorted_mults;
            for (size_t idx : indices) {
                sorted_roots.push_back(symbolic_roots[idx]);
                sorted_mults.push_back(multiplicities[idx]);
            }

            pos_case.solution = build_parametric_solution(sorted_roots, sorted_mults, 1, type);
            result.cases.push_back(pos_case);
        }

        {
            PiecewiseIntervalResult::Case neg_case;
            neg_case.condition = lamina::detail::make_expression_ptr(
                lamina::detail::make_node<RelationalNode>(
                    lamina::detail::node(leading_coeff),
                    lamina::detail::node(SymbolicExpr::number(0)),
                    RelationalNode::Op::LT));

            auto symbolic_roots = solve_symbolic_poly(poly, variable);
            std::vector<int> multiplicities(symbolic_roots.size(), 1);

            // Sort roots in ascending order and apply same permutation to multiplicities
            std::vector<size_t> indices(symbolic_roots.size());
            for (size_t i = 0; i < indices.size(); ++i) indices[i] = i;
            std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
                return root_less_than(symbolic_roots[a], symbolic_roots[b]);
            });
            std::vector<std::shared_ptr<SymbolicExpr>> sorted_roots;
            std::vector<int> sorted_mults;
            for (size_t idx : indices) {
                sorted_roots.push_back(symbolic_roots[idx]);
                sorted_mults.push_back(multiplicities[idx]);
            }

            neg_case.solution = build_parametric_solution(sorted_roots, sorted_mults, -1, type);
            result.cases.push_back(neg_case);
        }

        {
            PiecewiseIntervalResult::Case degen_case;
            degen_case.condition = lamina::detail::make_expression_ptr(
                lamina::detail::make_node<RelationalNode>(
                    lamina::detail::node(leading_coeff),
                    lamina::detail::node(SymbolicExpr::number(0)),
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
                                    merged.condition = lamina::detail::make_expression_ptr(
                                        lamina::detail::make_node<LogicalNode>(
                                            lamina::detail::node(degen_case.condition),
                                            lamina::detail::node(sub_case.condition),
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

                        int reduced_leading_sign = 0;
                        if (new_lc) {
                            reduced_leading_sign = exact_numeric_sign(new_lc);
                            if (reduced_leading_sign == 0 && !new_lc->is_zero()) {
                                if (auto val = try_checked_numeric_constant(*new_lc)) {
                                    reduced_leading_sign = (*val > 0) ? 1 : ((*val < 0) ? -1 : 0);
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
