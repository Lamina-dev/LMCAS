#include "limit_result.hpp"
#include "poly_utils.hpp"
#include "root_of_utils.hpp"
#include "solve_strategies.hpp"
#include "internal/calculus_utils_support.hpp"
#include "internal/symbolic_ast.hpp"

namespace LMCAS {

using namespace calculus_utils_detail;

static std::shared_ptr<const SymbolicNode> denominator_factor(
    const std::shared_ptr<const SymbolicNode>& node)
{
    auto power = std::dynamic_pointer_cast<const PowerNode>(node);
    if (!power) {
        return nullptr;
    }
    auto number = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
    if (!number) {
        return nullptr;
    }
    double exponent = 0;
    if (std::holds_alternative<lmmc_real_t>(number->value())) {
        exponent = std::get<lmmc_real_t>(number->value());
    } else if (std::holds_alternative<BigInt>(number->value())) {
        exponent = std::get<BigInt>(number->value()).to_double();
    } else if (std::holds_alternative<Rational>(number->value())) {
        exponent = std::get<Rational>(number->value()).to_double();
    }
    if (!(exponent < 0)) {
        return nullptr;
    }
    if (exponent == -1.0) {
        return power->base();
    }
    auto positive = LMCAS::detail::make_node<NumberNode>(BigInt(static_cast<int>(-exponent)));
    return LMCAS::detail::make_node<PowerNode>(power->base(), positive);
}

static std::shared_ptr<SymbolicExpr> calculus_utils_extract_denominator(
    const std::shared_ptr<SymbolicExpr>& expr)
{
    if (!expr || !LMCAS::detail::node(expr)) {
        return nullptr;
    }
    auto factor = denominator_factor(LMCAS::detail::node(expr));
    if (factor) {
        return LMCAS::detail::make_expression_ptr(factor);
    }
    auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr));
    if (!multiply) {
        return nullptr;
    }
    std::vector<std::shared_ptr<const SymbolicNode>> factors;
    for (const auto& operand : multiply->operands()) {
        auto denominator = denominator_factor(operand);
        if (denominator) factors.push_back(std::move(denominator));
    }
    if (factors.empty()) {
        return nullptr;
    }
    if (factors.size() == 1) {
        return LMCAS::detail::make_expression_ptr(factors[0]);
    }
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<MultiplyNode>(factors));
}

static std::shared_ptr<SymbolicExpr> calculus_utils_product(
    const std::vector<std::shared_ptr<const SymbolicNode>>& factors)
{
    if (factors.empty()) {
        return SymbolicExpr::number(1);
    }
    if (factors.size() == 1) {
        return LMCAS::detail::make_expression_ptr(factors[0]);
    }
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<MultiplyNode>(factors))->simplify();
}

static bool calculus_utils_is_minus_one(const std::shared_ptr<const NumberNode>& number)
{
    if (!number) {
        return false;
    }
    if (std::holds_alternative<BigInt>(number->value())) {
        return std::get<BigInt>(number->value()) == BigInt(-1);
    }
    if (std::holds_alternative<Rational>(number->value())) {
        return std::get<Rational>(number->value()) == Rational(-1);
    }
    return std::get<lmmc_real_t>(number->value()) == -1.0;
}

static bool calculus_utils_split_rational(
    const std::shared_ptr<SymbolicExpr>& expr,
    std::shared_ptr<SymbolicExpr>& numerator,
    std::shared_ptr<SymbolicExpr>& denominator)
{
    if (!expr || !LMCAS::detail::node(expr)) {
        return false;
    }

    std::vector<std::shared_ptr<const SymbolicNode>> num_factors;
    std::vector<std::shared_ptr<const SymbolicNode>> den_factors;

    auto collect_factor = [&](const std::shared_ptr<const SymbolicNode>& factor) {
        if (auto power = std::dynamic_pointer_cast<const PowerNode>(factor)) {
            auto exponent = std::dynamic_pointer_cast<const NumberNode>(power->exponent());
            if (calculus_utils_is_minus_one(exponent)) {
                den_factors.push_back(power->base());
                return;
            }
        }
        num_factors.push_back(factor);
    };

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(expr))) {
        for (const auto& factor : multiply->operands()) collect_factor(factor);
    } else {
        collect_factor(LMCAS::detail::node(expr));
    }

    if (den_factors.empty()) {
        return false;
    }
    numerator = calculus_utils_product(num_factors);
    denominator = calculus_utils_product(den_factors);
    return numerator && denominator;
}

static bool calculus_utils_oblique_from_rational(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    std::shared_ptr<SymbolicExpr>& slope,
    std::shared_ptr<SymbolicExpr>& intercept)
{
    std::shared_ptr<SymbolicExpr> numerator;
    std::shared_ptr<SymbolicExpr> denominator;
    if (!calculus_utils_split_rational(f, numerator, denominator)) {
        return false;
    }

    try {
        auto num_poly = symbolic_to_poly<Rational>(numerator->expand()->simplify(), var);
        auto den_poly = symbolic_to_poly<Rational>(denominator->expand()->simplify(), var);
        if (!num_poly || !den_poly || den_poly.value().is_zero()) {
            return false;
        }
        auto [quotient, remainder] = num_poly.value().div_mod(den_poly.value());
        (void)remainder;
        if (quotient.degree() != 1 || quotient.coeffs.size() < 2) {
            return false;
        }
        intercept = SymbolicExpr::number(quotient.coeffs[0])->simplify();
        slope = SymbolicExpr::number(quotient.coeffs[1])->simplify();
        return !slope->is_zero();
    } catch (const std::invalid_argument&) {
        return false;
    } catch (const std::out_of_range&) {
        return false;
    }
}
static SymbolicExprVectorResult asymptote_candidates(
    const std::shared_ptr<SymbolicExpr>& denominator,
    const std::string& var, ComputationContext& context,
    const std::string& operation)
{
    auto solved = solve_equation(denominator, var, context, SolveOptions{});
    if (!solved) {
        return SymbolicExprVectorResult::failure(solved.error());
    }
    const auto& zero_set = solved.value();
    const auto* finite = std::get_if<FiniteSolutions>(&zero_set);
    if (!std::holds_alternative<EmptySolutions>(zero_set) && !finite) {
        return SymbolicExprVectorResult::failure(
            CasErrc::Inconclusive,
            "asymptote denominator zeros are outside the finite exact support domain",
            operation);
    }
    std::vector<std::shared_ptr<SymbolicExpr>> zeros;
    if (!finite) {
        return SymbolicExprVectorResult::success(std::move(zeros));
    }
    zeros.reserve(finite->values.size());
    for (const auto& solution : finite->values) {
        if (!solution.value || !LMCAS::detail::node(solution.value)) {
            return SymbolicExprVectorResult::failure(
                CasErrc::InternalInvariant,
                "checked solver returned a null vertical-asymptote candidate", operation);
        }
        auto candidate = rootof_simplify(solution.value);
        if (!candidate || !LMCAS::detail::node(candidate)) {
            return SymbolicExprVectorResult::failure(
                CasErrc::InternalInvariant,
                "checked solver returned a malformed vertical-asymptote candidate", operation);
        }
        if (LMCAS::detail::contains_node_type<RootOfNode>(LMCAS::detail::node(candidate))) {
            return SymbolicExprVectorResult::failure(
                CasErrc::Inconclusive,
                "vertical-asymptote candidate could not be reduced to an exact explicit point",
                operation);
        }
        zeros.push_back(candidate);
    }
    return SymbolicExprVectorResult::success(std::move(zeros));
}

static Result<bool> is_vertical_asymptote(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point, ComputationContext& context)
{
    auto right = limit_checked(f, var, point, LimitDirection::FromAbove, context);
    if (!right && right.error().code != CasErrc::Inconclusive) {
        return Result<bool>::failure(right.error());
    }
    auto left = limit_checked(f, var, point, LimitDirection::FromBelow, context);
    if (!left && left.error().code != CasErrc::Inconclusive) {
        return Result<bool>::failure(left.error());
    }
    auto infinite = [](const LimitResult& value) {
        return value && (std::holds_alternative<PositiveInfinityLimit>(value.value().value) ||
                         std::holds_alternative<NegativeInfinityLimit>(value.value().value));
    };
    if (infinite(right) || infinite(left)) { return true; }
    if (!right) { return Result<bool>::failure(right.error()); }
    if (!left) { return Result<bool>::failure(left.error()); }
    return false;
}

static Result<void> collect_vertical_asymptotes(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    ComputationContext& context, AsymptoteResult& result,
    const std::string& operation)
{
    auto denominator = calculus_utils_extract_denominator(f);
    if (!denominator) {
        return Result<void>::success();
    }
    auto zeros = asymptote_candidates(denominator, var, context, operation);
    if (!zeros) {
        return Result<void>::failure(zeros.error());
    }
    for (const auto& point : zeros.value()) {
        if (!point || calculus_utils_is_infinity(point)) {
            continue;
        }
        auto vertical = is_vertical_asymptote(f, var, point, context);
        if (!vertical) {
            return Result<void>::failure(vertical.error());
        }
        if (vertical.value()) result.vertical.push_back(point);
    }
    return Result<void>::success();
}

static bool has_oblique_asymptote(
    const AsymptoteResult& result, const std::shared_ptr<SymbolicExpr>& slope,
    const std::shared_ptr<SymbolicExpr>& intercept)
{
    for (const auto& entry : result.oblique) {
        if (calculus_utils_expr_equal(entry.first, slope) &&
            calculus_utils_expr_equal(entry.second, intercept)) {
            return true;
        }
    }
    return false;
}

static Result<void> collect_oblique_asymptote(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& x_expr,
    const std::shared_ptr<SymbolicExpr>& infinity, ComputationContext& context,
    AsymptoteResult& result)
{
    auto quotient = SymbolicExpr::multiply(f,
        SymbolicExpr::power(x_expr, SymbolicExpr::number(-1)));
    auto slope_result = limit_expression_checked(
        quotient, var, infinity, LimitDirection::Both, context);
    if (!slope_result) {
        return Result<void>::failure(slope_result.error());
    }
    auto slope = std::move(slope_result.value());
    if (slope) slope = slope->simplify();
    if (!slope || calculus_utils_is_infinity(slope) || slope->is_zero()) {
        return Result<void>::success();
    }
    auto linear = SymbolicExpr::multiply(slope, x_expr);
    auto remainder = SymbolicExpr::add(f,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), linear));
    auto intercept_result = limit_expression_checked(
        remainder, var, infinity, LimitDirection::Both, context);
    if (!intercept_result) {
        return Result<void>::failure(intercept_result.error());
    }
    auto intercept = std::move(intercept_result.value());
    if (intercept) intercept = intercept->simplify();
    if (intercept && !calculus_utils_is_infinity(intercept)) {
        if (!has_oblique_asymptote(result, slope, intercept)) {
            result.oblique.emplace_back(slope, intercept);
        }
    }
    return Result<void>::success();
}

static bool simplify_finite_limit(std::shared_ptr<SymbolicExpr>& value)
{
    if (value) {
        value = value->simplify();
    }
    return value && !calculus_utils_is_infinity(value);
}

static Result<void> collect_infinite_asymptotes(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& pos_inf,
    const std::shared_ptr<SymbolicExpr>& neg_inf,
    ComputationContext& context, AsymptoteResult& result)
{
    auto positive = limit_expression_checked(f, var, pos_inf, LimitDirection::Both, context);
    if (!positive) {
        return Result<void>::failure(positive.error());
    }
    auto negative = limit_expression_checked(f, var, neg_inf, LimitDirection::Both, context);
    if (!negative) {
        return Result<void>::failure(negative.error());
    }
    auto lim_pos = std::move(positive.value());
    auto lim_neg = std::move(negative.value());
    const bool horizontal_pos = simplify_finite_limit(lim_pos);
    const bool horizontal_neg = simplify_finite_limit(lim_neg);
    if (horizontal_pos) result.horizontal.push_back(lim_pos);
    if (horizontal_neg) {
        if (!horizontal_pos || !calculus_utils_expr_equal(lim_pos, lim_neg)) {
            result.horizontal.push_back(lim_neg);
        }
    }
    auto x_expr = SymbolicExpr::variable(var);
    if (!horizontal_pos && !horizontal_neg) {
        std::shared_ptr<SymbolicExpr> slope;
        std::shared_ptr<SymbolicExpr> intercept;
        if (calculus_utils_oblique_from_rational(f, var, slope, intercept)) {
            result.oblique.emplace_back(slope, intercept);
            return Result<void>::success();
        }
    }
    if (!horizontal_pos) {
        auto status = collect_oblique_asymptote(f, var, x_expr, pos_inf, context, result);
        if (!status) {
            return status;
        }
    }
    if (!horizontal_neg) {
        return collect_oblique_asymptote(f, var, x_expr, neg_inf, context, result);
    }
    return Result<void>::success();
}

AsymptoteAnalysisResult asymptotes_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    ComputationContext& context)
{
    try {
        AsymptoteResult result;
        const std::string operation = "asymptotes";
        auto input = calculus_utils_validate_expr(f, var, context, operation);
        if (!input) {
            return AsymptoteAnalysisResult::failure(input.error());
        }
        auto pos_inf = SymbolicExpr::infinity(1);
        auto neg_inf = SymbolicExpr::infinity(-1);
        auto vertical = collect_vertical_asymptotes(f, var, context, result, operation);
        if (!vertical) {
            return AsymptoteAnalysisResult::failure(vertical.error());
        }
        auto infinite = collect_infinite_asymptotes(f, var, pos_inf, neg_inf, context, result);
        if (!infinite) {
            return AsymptoteAnalysisResult::failure(infinite.error());
        }
        return AsymptoteAnalysisResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return AsymptoteAnalysisResult::failure(
            CasErrc::ResourceLimit, "asymptote analysis allocation failed", "asymptotes");
    } catch (const std::exception& ex) {
        return AsymptoteAnalysisResult::failure(
            CasErrc::InternalInvariant, ex.what(), "asymptotes");
    }
}

AsymptoteAnalysisResult asymptotes_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    ComputationContext context;
    return asymptotes_checked(f, var, context);
}

}
