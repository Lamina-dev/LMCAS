/**
 * @file calculus_utils.cpp
 * @brief 微积分工具函数实现：连续性判定、渐近线分析、对数微分、微分、全微分、反函数导数、反函数求解。
 */

#include "calculus_utils.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "integration.hpp"
#include "numeric_evaluation.hpp"
#include "poly_utils.hpp"
#include "root_of_utils.hpp"
#include "solve_strategies.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <memory>

namespace lamina {

namespace {

Result<void> calculus_utils_validate_expr(const std::shared_ptr<SymbolicExpr>& expr,
                                          const std::string& var,
                                          ComputationContext& context,
                                          const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!expr || !lamina::detail::node(expr)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "expression cannot be null", operation);
    }
    if (var.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "variable name cannot be empty", operation);
    }
    return Result<void>::success();
}

Result<void> calculus_utils_validate_two_exprs(
    const std::shared_ptr<SymbolicExpr>& first,
    const std::shared_ptr<SymbolicExpr>& second,
    const std::string& var,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!first || !lamina::detail::node(first) || !second || !lamina::detail::node(second)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "expression cannot be null", operation);
    }
    if (var.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "variable name cannot be empty", operation);
    }
    return Result<void>::success();
}

Result<void> calculus_utils_validate_expr_target(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& target,
    ComputationContext& context,
    const std::string& operation)
{
    auto input = calculus_utils_validate_expr(expr, var, context, operation);
    if (!input) return input;
    if (!target || !lamina::detail::node(target)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "target expression cannot be null",
                                     operation);
    }
    return Result<void>::success();
}

Result<void> calculus_utils_validate_expr_bounds(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    const std::string& operation)
{
    auto input = calculus_utils_validate_expr(expr, var, context, operation);
    if (!input) return input;
    if (!a || !lamina::detail::node(a) || !b || !lamina::detail::node(b)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "bounds cannot be null", operation);
    }
    return Result<void>::success();
}

} // namespace


/**
 * @internal
 * @brief 判断表达式是否包含无穷大节点。
 */
static bool calculus_utils_is_infinity(const std::shared_ptr<SymbolicExpr>& expr)
{
    if (!expr || !lamina::detail::node(expr)) return false;

    if (auto f = std::dynamic_pointer_cast<const FunctionNode>(lamina::detail::node(expr))) {
        return f->type() == FunctionNode::FuncType::Infinity;
    }
    if (auto m = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
        for (auto& op : m->operands()) {
            if (auto f = std::dynamic_pointer_cast<const FunctionNode>(op)) {
                if (f->type() == FunctionNode::FuncType::Infinity) return true;
            }
        }
    }
    return false;
}

/**
 * @internal
 * @brief 判断两个表达式是否结构相等。
 */
static bool calculus_utils_expr_equal(const std::shared_ptr<SymbolicExpr>& a,
                                      const std::shared_ptr<SymbolicExpr>& b)
{
    if (!a || !lamina::detail::node(a) || !b || !lamina::detail::node(b)) {
        return (!a || !lamina::detail::node(a)) && (!b || !lamina::detail::node(b));
    }
    return lamina::detail::node(a)->equals(*lamina::detail::node(b));
}

static bool calculus_utils_contains_rootof(const std::shared_ptr<const SymbolicNode>& node)
{
    if (!node) return false;

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (function->type() == FunctionNode::FuncType::RootOf) return true;
        for (const auto& argument : function->arguments()) {
            if (calculus_utils_contains_rootof(argument)) return true;
        }
        return false;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& operand : add->operands()) {
            if (calculus_utils_contains_rootof(operand)) return true;
        }
        return false;
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& operand : multiply->operands()) {
            if (calculus_utils_contains_rootof(operand)) return true;
        }
        return false;
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return calculus_utils_contains_rootof(power->base()) ||
               calculus_utils_contains_rootof(power->exponent());
    }
    return false;
}

static bool calculus_utils_contains_complex(const std::shared_ptr<const SymbolicNode>& node)
{
    if (!node) return false;

    if (std::dynamic_pointer_cast<const ComplexNode>(node)) return true;
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        for (const auto& argument : function->arguments()) {
            if (calculus_utils_contains_complex(argument)) return true;
        }
        return false;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& operand : add->operands()) {
            if (calculus_utils_contains_complex(operand)) return true;
        }
        return false;
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& operand : multiply->operands()) {
            if (calculus_utils_contains_complex(operand)) return true;
        }
        return false;
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return calculus_utils_contains_complex(power->base()) ||
               calculus_utils_contains_complex(power->exponent());
    }
    return false;
}

static bool calculus_utils_is_negative_number(const std::shared_ptr<const NumberNode>& number)
{
    if (!number) return false;
    if (std::holds_alternative<BigInt>(number->value())) {
        return std::get<BigInt>(number->value()) < BigInt(0);
    }
    if (std::holds_alternative<Rational>(number->value())) {
        return std::get<Rational>(number->value()) < Rational(0);
    }
    return std::get<lmmc_real_t>(number->value()) < 0.0;
}

static bool calculus_utils_contains_nonreal_sqrt(const std::shared_ptr<const SymbolicNode>& node)
{
    if (!node) return false;

    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (function->type() == FunctionNode::FuncType::Sqrt &&
            function->arguments().size() == 1) {
            auto argument = lamina::detail::make_expression_ptr(function->arguments()[0])->simplify();
            auto number = argument
                ? std::dynamic_pointer_cast<const NumberNode>(lamina::detail::node(argument))
                : nullptr;
            if (calculus_utils_is_negative_number(number)) return true;
        }
        for (const auto& argument : function->arguments()) {
            if (calculus_utils_contains_nonreal_sqrt(argument)) return true;
        }
        return false;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& operand : add->operands()) {
            if (calculus_utils_contains_nonreal_sqrt(operand)) return true;
        }
        return false;
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& operand : multiply->operands()) {
            if (calculus_utils_contains_nonreal_sqrt(operand)) return true;
        }
        return false;
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return calculus_utils_contains_nonreal_sqrt(power->base()) ||
               calculus_utils_contains_nonreal_sqrt(power->exponent());
    }
    if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        return calculus_utils_contains_nonreal_sqrt(complex->real()) ||
               calculus_utils_contains_nonreal_sqrt(complex->imag());
    }
    return false;
}

enum class CalculusRealCandidateStatus {
    Real,
    NonReal,
    Inconclusive
};

static CalculusRealCandidateStatus calculus_utils_extract_real_candidate(
    const std::shared_ptr<SymbolicExpr>& candidate,
    std::shared_ptr<SymbolicExpr>& real_candidate)
{
    real_candidate.reset();
    auto simplified_candidate = rootof_simplify(candidate);
    if (!simplified_candidate || !lamina::detail::node(simplified_candidate)) {
        return CalculusRealCandidateStatus::Inconclusive;
    }
    if (calculus_utils_contains_rootof(lamina::detail::node(simplified_candidate))) {
        return CalculusRealCandidateStatus::Inconclusive;
    }
    if (calculus_utils_contains_nonreal_sqrt(lamina::detail::node(simplified_candidate))) {
        return CalculusRealCandidateStatus::NonReal;
    }

    if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(
            lamina::detail::node(simplified_candidate))) {
        auto imag = lamina::detail::make_expression_ptr(complex->imag())->simplify();
        if (!imag || !lamina::detail::node(imag)) {
            return CalculusRealCandidateStatus::Inconclusive;
        }
        if (imag->is_zero()) {
            auto real = lamina::detail::make_expression_ptr(complex->real())->simplify();
            if (!real || !lamina::detail::node(real)) {
                return CalculusRealCandidateStatus::Inconclusive;
            }
            real_candidate = real;
            return CalculusRealCandidateStatus::Real;
        }
        return CalculusRealCandidateStatus::NonReal;
    }

    if (calculus_utils_contains_complex(lamina::detail::node(simplified_candidate))) {
        return CalculusRealCandidateStatus::Inconclusive;
    }

    real_candidate = simplified_candidate;
    return CalculusRealCandidateStatus::Real;
}

/**
 * @internal
 * @brief 从表达式中提取分母部分。
 *
 * 对于 PowerNode 且指数为负，base 即为分母。
 * 对于 MultiplyNode，查找指数为负的因子作为分母。
 */
static std::shared_ptr<SymbolicExpr> calculus_utils_extract_denominator(
    const std::shared_ptr<SymbolicExpr>& expr)
{
    if (!expr || !lamina::detail::node(expr)) return nullptr;

    /// 检查是否为 PowerNode 且指数为负（如 x^(-1) 表示 1/x）
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(lamina::detail::node(expr))) {
        if (auto num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
            double e = 0;
            if (std::holds_alternative<lmmc_real_t>(num->value()))
                e = std::get<lmmc_real_t>(num->value());
            else if (std::holds_alternative<BigInt>(num->value()))
                e = std::get<BigInt>(num->value()).to_double();
            else if (std::holds_alternative<Rational>(num->value()))
                e = std::get<Rational>(num->value()).to_double();

            if (e < 0) {
                if (e == -1.0) {
                    return lamina::detail::make_expression_ptr(pow->base());
                }
                auto pos_exp = lamina::detail::make_node<NumberNode>(BigInt(static_cast<int>(-e)));
                auto den_node = lamina::detail::make_node<PowerNode>(pow->base(), pos_exp);
                return lamina::detail::make_expression_ptr(den_node);
            }
        }
    }

    /// 检查 MultiplyNode 中是否有负指数因子
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
        std::vector<std::shared_ptr<const SymbolicNode>> den_factors;

        for (auto& op : mul->operands()) {
            if (auto pow = std::dynamic_pointer_cast<const PowerNode>(op)) {
                if (auto num = std::dynamic_pointer_cast<const NumberNode>(pow->exponent())) {
                    double e = 0;
                    if (std::holds_alternative<lmmc_real_t>(num->value()))
                        e = std::get<lmmc_real_t>(num->value());
                    else if (std::holds_alternative<BigInt>(num->value()))
                        e = std::get<BigInt>(num->value()).to_double();
                    else if (std::holds_alternative<Rational>(num->value()))
                        e = std::get<Rational>(num->value()).to_double();

                    if (e < 0) {
                        if (e == -1.0) {
                            den_factors.push_back(pow->base());
                        } else {
                            auto pos_exp = lamina::detail::make_node<NumberNode>(BigInt(static_cast<int>(-e)));
                            den_factors.push_back(lamina::detail::make_node<PowerNode>(pow->base(), pos_exp));
                        }
                    }
                }
            }
        }

        if (!den_factors.empty()) {
            std::shared_ptr<const SymbolicNode> den_node;
            if (den_factors.size() == 1) {
                den_node = den_factors[0];
            } else {
                den_node = lamina::detail::make_node<MultiplyNode>(den_factors);
            }
            return lamina::detail::make_expression_ptr(den_node);
        }
    }

    return nullptr;
}

static std::shared_ptr<SymbolicExpr> calculus_utils_product(
    const std::vector<std::shared_ptr<const SymbolicNode>>& factors)
{
    if (factors.empty()) return SymbolicExpr::number(1);
    if (factors.size() == 1) return lamina::detail::make_expression_ptr(factors[0]);
    return lamina::detail::make_expression_ptr(lamina::detail::make_node<MultiplyNode>(factors))->simplify();
}

static bool calculus_utils_is_minus_one(const std::shared_ptr<const NumberNode>& number)
{
    if (!number) return false;
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
    if (!expr || !lamina::detail::node(expr)) return false;

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

    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(lamina::detail::node(expr))) {
        for (const auto& factor : multiply->operands()) collect_factor(factor);
    } else {
        collect_factor(lamina::detail::node(expr));
    }

    if (den_factors.empty()) return false;
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
    if (!calculus_utils_split_rational(f, numerator, denominator)) return false;

    try {
        auto num_poly = symbolic_to_poly<Rational>(numerator->expand()->simplify(), var);
        auto den_poly = symbolic_to_poly<Rational>(denominator->expand()->simplify(), var);
        if (den_poly.is_zero()) return false;
        auto [quotient, remainder] = num_poly.div_mod(den_poly);
        (void)remainder;
        if (quotient.degree() != 1 || quotient.coeffs.size() < 2) return false;
        intercept = SymbolicExpr::number(quotient.coeffs[0])->simplify();
        slope = SymbolicExpr::number(quotient.coeffs[1])->simplify();
        return !slope->is_zero();
    } catch (...) {
        return false;
    }
}


ContinuityType continuity_at(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point)
{
    if (!f || !point) return ContinuityType::Essential;

    /// 计算左极限
    auto left_lim = f->limit(var, point, "-");
    /// 计算右极限
    auto right_lim = f->limit(var, point, "+");
    /// 计算函数值（直接代入）
    auto func_val = f->substitute(var, point);

    /// 简化结果
    if (left_lim) left_lim = left_lim->simplify();
    if (right_lim) right_lim = right_lim->simplify();
    if (func_val) func_val = func_val->simplify();

    /// 检查极限是否存在（不为无穷）
    bool left_exists = left_lim && !calculus_utils_is_infinity(left_lim);
    bool right_exists = right_lim && !calculus_utils_is_infinity(right_lim);

    /// 如果任一侧极限不存在或为无穷 → 本性间断点
    if (!left_exists || !right_exists) {
        return ContinuityType::Essential;
    }

    /// 检查左极限是否等于右极限
    bool limits_equal = calculus_utils_expr_equal(left_lim, right_lim);

    if (!limits_equal) {
        /// 左极限 ≠ 右极限 → 跳跃间断点
        return ContinuityType::Jump;
    }

    /// 左极限 = 右极限，检查是否等于函数值
    bool val_exists = func_val && !calculus_utils_is_infinity(func_val);

    if (!val_exists || !calculus_utils_expr_equal(left_lim, func_val)) {
        /// 极限存在但不等于函数值（或函数无定义） → 可去间断点
        return ContinuityType::Removable;
    }

    /// 左极限 = 右极限 = 函数值 → 连续
    return ContinuityType::Continuous;
}


AsymptoteAnalysisResult asymptotes_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    ComputationContext& context)
{
    AsymptoteResult result;
    const std::string operation = "asymptotes";
    auto input = calculus_utils_validate_expr(f, var, context, operation);
    if (!input) return AsymptoteAnalysisResult::failure(input.error());

    auto pos_inf = SymbolicExpr::infinity(1);
    auto neg_inf = SymbolicExpr::infinity(-1);

    /// 提取分母，求解分母 = 0 的点
    auto denominator = calculus_utils_extract_denominator(f);
    if (denominator) {
        auto solved_zeros = solve_dispatch_checked(denominator, var, context, SolveOptions{});
        if (!solved_zeros) return AsymptoteAnalysisResult::failure(solved_zeros.error());

        const auto& zero_set = solved_zeros.value();
        if (zero_set.kind() != SolutionSet::Kind::Empty &&
            zero_set.kind() != SolutionSet::Kind::Finite) {
            return AsymptoteAnalysisResult::failure(
                CasErrc::Inconclusive,
                "asymptote denominator zeros are outside the finite exact support domain",
                operation);
        }

        std::vector<std::shared_ptr<SymbolicExpr>> zeros;
        if (zero_set.kind() == SolutionSet::Kind::Finite) {
            zeros.reserve(zero_set.finite_solutions().size());
            for (const auto& solution : zero_set.finite_solutions()) {
                if (!solution.value || !lamina::detail::node(solution.value)) {
                    return AsymptoteAnalysisResult::failure(
                        CasErrc::InternalInvariant,
                        "checked solver returned a null vertical-asymptote candidate",
                        operation);
                }
                auto candidate = rootof_simplify(solution.value);
                if (!candidate || !lamina::detail::node(candidate)) {
                    return AsymptoteAnalysisResult::failure(
                        CasErrc::InternalInvariant,
                        "checked solver returned a malformed vertical-asymptote candidate",
                        operation);
                }
                if (calculus_utils_contains_rootof(lamina::detail::node(candidate))) {
                    return AsymptoteAnalysisResult::failure(
                        CasErrc::Inconclusive,
                        "vertical-asymptote candidate could not be reduced to an exact explicit point",
                        operation);
                }
                zeros.push_back(candidate);
            }
        }

        for (auto& z : zeros) {
            if (!z || calculus_utils_is_infinity(z)) continue;

            /// 验证该点处极限为 ±∞
            auto lim_at_z = f->limit(var, z);
            if (calculus_utils_is_infinity(lim_at_z)) {
                result.vertical.push_back(z);
            } else {
                /// 尝试单侧极限
                auto lim_right = f->limit(var, z, "+");
                auto lim_left = f->limit(var, z, "-");
                if (calculus_utils_is_infinity(lim_right) || calculus_utils_is_infinity(lim_left)) {
                    result.vertical.push_back(z);
                }
            }
        }
    }

    /// 计算 x→+∞ 和 x→-∞ 的极限
    auto lim_pos = f->limit(var, pos_inf);
    auto lim_neg = f->limit(var, neg_inf);

    if (lim_pos) lim_pos = lim_pos->simplify();
    if (lim_neg) lim_neg = lim_neg->simplify();

    bool has_horiz_pos = lim_pos && !calculus_utils_is_infinity(lim_pos);
    bool has_horiz_neg = lim_neg && !calculus_utils_is_infinity(lim_neg);

    if (has_horiz_pos) {
        result.horizontal.push_back(lim_pos);
    }
    if (has_horiz_neg) {
        /// 避免重复添加相同的水平渐近线
        bool duplicate = has_horiz_pos && calculus_utils_expr_equal(lim_pos, lim_neg);
        if (!duplicate) {
            result.horizontal.push_back(lim_neg);
        }
    }

    /// 仅当对应方向无水平渐近线时才检查斜渐近线
    auto x_expr = SymbolicExpr::variable(var);

    if (!has_horiz_pos && !has_horiz_neg) {
        std::shared_ptr<SymbolicExpr> slope;
        std::shared_ptr<SymbolicExpr> intercept;
        if (calculus_utils_oblique_from_rational(f, var, slope, intercept)) {
            result.oblique.emplace_back(slope, intercept);
            return AsymptoteAnalysisResult::success(std::move(result));
        }
    }

    if (!has_horiz_pos) {
        /// 计算 slope = lim(f/x) as x→+∞
        auto f_over_x = SymbolicExpr::multiply(f, SymbolicExpr::power(x_expr, SymbolicExpr::number(-1)));
        auto slope_pos = f_over_x->limit(var, pos_inf);
        if (slope_pos) slope_pos = slope_pos->simplify();

        if (slope_pos && !calculus_utils_is_infinity(slope_pos) && !slope_pos->is_zero()) {
            /// 计算 intercept = lim(f - slope*x) as x→+∞
            auto slope_times_x = SymbolicExpr::multiply(slope_pos, x_expr);
            auto f_minus_mx = SymbolicExpr::add(f, SymbolicExpr::multiply(SymbolicExpr::number(-1), slope_times_x));
            auto intercept_pos = f_minus_mx->limit(var, pos_inf);
            if (intercept_pos) intercept_pos = intercept_pos->simplify();

            if (intercept_pos && !calculus_utils_is_infinity(intercept_pos)) {
                result.oblique.emplace_back(slope_pos, intercept_pos);
            }
        }
    }

    if (!has_horiz_neg) {
        /// 计算 slope = lim(f/x) as x→-∞
        auto f_over_x = SymbolicExpr::multiply(f, SymbolicExpr::power(x_expr, SymbolicExpr::number(-1)));
        auto slope_neg = f_over_x->limit(var, neg_inf);
        if (slope_neg) slope_neg = slope_neg->simplify();

        if (slope_neg && !calculus_utils_is_infinity(slope_neg) && !slope_neg->is_zero()) {
            /// 计算 intercept = lim(f - slope*x) as x→-∞
            auto slope_times_x = SymbolicExpr::multiply(slope_neg, x_expr);
            auto f_minus_mx = SymbolicExpr::add(f, SymbolicExpr::multiply(SymbolicExpr::number(-1), slope_times_x));
            auto intercept_neg = f_minus_mx->limit(var, neg_inf);
            if (intercept_neg) intercept_neg = intercept_neg->simplify();

            if (intercept_neg && !calculus_utils_is_infinity(intercept_neg)) {
                /// 避免重复：如果 +∞ 方向已有相同斜渐近线则跳过
                bool duplicate = false;
                for (auto& [s, i] : result.oblique) {
                    if (calculus_utils_expr_equal(s, slope_neg) &&
                        calculus_utils_expr_equal(i, intercept_neg)) {
                        duplicate = true;
                        break;
                    }
                }
                if (!duplicate) {
                    result.oblique.emplace_back(slope_neg, intercept_neg);
                }
            }
        }
    }

    return AsymptoteAnalysisResult::success(std::move(result));
}

AsymptoteAnalysisResult asymptotes_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    ComputationContext context;
    return asymptotes_checked(f, var, context);
}

AsymptoteResult asymptotes(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    auto checked = asymptotes_checked(f, var);
    if (!checked) return {};
    return std::move(checked.value());
}


std::shared_ptr<SymbolicExpr> log_differentiate(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    if (!f || !lamina::detail::node(f)) return nullptr;

    /// 计算 ln(f)
    auto ln_f = SymbolicExpr::ln(f);

    /// 对 ln(f) 求导 — 化简器会自动应用对数规则
    /// （ln(a*b) = ln(a)+ln(b), ln(a^n) = n*ln(a)）
    auto ln_f_simplified = ln_f->simplify();
    auto d_ln_f = ln_f_simplified->differentiate(var);

    /// 结果 = f * d/dx[ln(f)]
    auto result = SymbolicExpr::multiply(f, d_ln_f);
    return result->simplify();
}

std::shared_ptr<SymbolicExpr> differential(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    if (!f || !lamina::detail::node(f)) return nullptr;

    /// 微分 df = f'(var) * dx，返回系数 f'(var)
    return f->differentiate(var);
}

std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::string>> total_differential(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars)
{
    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, std::string>> result;
    if (!f || !lamina::detail::node(f)) return result;

    result.reserve(vars.size());
    for (const auto& v : vars) {
        auto partial = f->differentiate(v);
        result.emplace_back(partial, v);
    }
    return result;
}

SymbolicExprVectorResult inverse_function_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& y, ComputationContext& context);

SymbolicExprResult inverse_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point,
    ComputationContext& context)
{
    const std::string operation = "inverse_derivative";
    auto input = calculus_utils_validate_expr_target(f, var, point, context, operation);
    if (!input) return SymbolicExprResult::failure(input.error());
    auto step = context.consume_steps(4, operation);
    if (!step) return SymbolicExprResult::failure(step.error());

    try {
        auto eq = SymbolicExpr::add(
            f, SymbolicExpr::multiply(SymbolicExpr::number(-1), point));
        auto simplified_eq = eq ? eq->simplify() : nullptr;
        if (simplified_eq && lamina::detail::node(simplified_eq)) eq = simplified_eq;

        std::vector<std::shared_ptr<SymbolicExpr>> candidates;
        auto inverse = inverse_function_checked(f, var, point, context);
        if (inverse) {
            candidates = inverse.value();
        } else {
            return SymbolicExprResult::failure(inverse.error());
        }

        std::vector<std::shared_ptr<SymbolicExpr>> real_candidates;
        real_candidates.reserve(candidates.size());
        for (const auto& candidate : candidates) {
            std::shared_ptr<SymbolicExpr> real_candidate;
            auto status = calculus_utils_extract_real_candidate(candidate, real_candidate);
            if (status == CalculusRealCandidateStatus::Inconclusive) {
                return SymbolicExprResult::failure(
                    CasErrc::Inconclusive,
                    "inverse derivative candidate reality could not be proven",
                    operation);
            }
            if (status == CalculusRealCandidateStatus::NonReal) {
                continue;
            }

            auto residual = eq->substitute(var, real_candidate);
            auto simplified_residual = residual ? residual->simplify() : nullptr;
            if (!simplified_residual || !lamina::detail::node(simplified_residual)) {
                return SymbolicExprResult::failure(
                    CasErrc::InternalInvariant,
                    "inverse derivative candidate verification produced a null residual",
                    operation);
            }
            if (!simplified_residual->is_zero()) {
                continue;
            }

            bool duplicate = false;
            for (const auto& existing : real_candidates) {
                if (existing && lamina::detail::node(existing) &&
                    lamina::detail::node(existing)->equals(*lamina::detail::node(real_candidate))) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) real_candidates.push_back(real_candidate);
        }
        candidates = std::move(real_candidates);

        if (candidates.empty()) {
            return SymbolicExprResult::failure(
                CasErrc::DomainError,
                "inverse point is outside the proven range",
                operation);
        }
        if (candidates.size() != 1) {
            return SymbolicExprResult::failure(
                CasErrc::Inconclusive,
                "inverse derivative requires a unique inverse branch",
                operation);
        }

        auto f_prime = f->differentiate(var);
        if (!f_prime || !lamina::detail::node(f_prime)) {
            return SymbolicExprResult::failure(
                CasErrc::Inconclusive,
                "inverse derivative could not construct the derivative",
                operation);
        }

        auto f_prime_at_x0 = f_prime->substitute(var, candidates[0]);
        auto simplified_derivative = f_prime_at_x0 ? f_prime_at_x0->simplify() : nullptr;
        if (!simplified_derivative || !lamina::detail::node(simplified_derivative)) {
            return SymbolicExprResult::failure(
                CasErrc::InternalInvariant,
                "inverse derivative substitution produced a null expression",
                operation);
        }

        if (simplified_derivative->is_zero()) {
            return SymbolicExprResult::failure(
                CasErrc::DomainError,
                "inverse derivative is undefined where f' is zero",
                operation);
        }

        auto result = SymbolicExpr::divide(SymbolicExpr::number(1), simplified_derivative);
        auto simplified = result ? result->simplify() : nullptr;
        if (!simplified || !lamina::detail::node(simplified)) {
            return SymbolicExprResult::failure(
                CasErrc::InternalInvariant,
                "inverse derivative result construction failed",
                operation);
        }
        return SymbolicExprResult::success(simplified);
    } catch (const std::bad_alloc&) {
        return SymbolicExprResult::failure(CasErrc::ResourceLimit,
                                           "inverse derivative allocation failed",
                                           operation);
    } catch (const std::exception& e) {
        return SymbolicExprResult::failure(CasErrc::InternalInvariant,
                                           e.what(), operation);
    }
}

SymbolicExprResult inverse_derivative_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point)
{
    ComputationContext context;
    return inverse_derivative_checked(f, var, point, context);
}

SymbolicExprVectorResult inverse_function_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& y, ComputationContext& context)
{
    const std::string operation = "inverse_function";
    auto input = calculus_utils_validate_expr_target(f, var, y, context, operation);
    if (!input) return SymbolicExprVectorResult::failure(input.error());
    auto step = context.consume_steps(4, operation);
    if (!step) return SymbolicExprVectorResult::failure(step.error());

    try {
        auto eq = SymbolicExpr::add(
            f, SymbolicExpr::multiply(SymbolicExpr::number(-1), y));
        auto simplified_eq = eq ? eq->simplify() : nullptr;
        if (simplified_eq && lamina::detail::node(simplified_eq)) eq = simplified_eq;

        auto solved = solve_dispatch_checked(eq, var, context, SolveOptions{});
        if (!solved) return SymbolicExprVectorResult::failure(solved.error());

        const auto& solutions = solved.value();
        if (solutions.kind() == SolutionSet::Kind::Empty) {
            return SymbolicExprVectorResult::success({});
        }
        if (solutions.kind() == SolutionSet::Kind::Finite) {
            std::vector<std::shared_ptr<SymbolicExpr>> values;
            values.reserve(solutions.finite_solutions().size());
            for (const auto& solution : solutions.finite_solutions()) {
                if (!solution.value || !lamina::detail::node(solution.value)) {
                    return SymbolicExprVectorResult::failure(
                        CasErrc::InternalInvariant,
                        "checked solver returned a null inverse candidate",
                        operation);
                }
                auto candidate = rootof_simplify(solution.value);
                if (!candidate || !lamina::detail::node(candidate)) {
                    return SymbolicExprVectorResult::failure(
                        CasErrc::InternalInvariant,
                        "checked solver returned a malformed inverse candidate",
                        operation);
                }
                values.push_back(candidate);
            }
            return SymbolicExprVectorResult::success(std::move(values));
        }

        return SymbolicExprVectorResult::failure(
            CasErrc::Inconclusive,
            "inverse equation is outside the finite exact support domain",
            operation);
    } catch (const std::bad_alloc&) {
        return SymbolicExprVectorResult::failure(CasErrc::ResourceLimit,
                                                 "inverse function allocation failed",
                                                 operation);
    } catch (const std::exception& e) {
        return SymbolicExprVectorResult::failure(CasErrc::InternalInvariant,
                                                 e.what(), operation);
    }
}

SymbolicExprVectorResult inverse_function_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& y)
{
    ComputationContext context;
    return inverse_function_checked(f, var, y, context);
}

std::shared_ptr<SymbolicExpr> inverse_derivative(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& point)
{
    auto checked = inverse_derivative_checked(f, var, point);
    if (checked) return checked.value();
    return nullptr;
}

std::vector<std::shared_ptr<SymbolicExpr>> inverse_function(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& y)
{
    auto checked = inverse_function_checked(f, var, y);
    if (!checked) return {};
    return std::move(checked.value());
}


/// 构造绝对值表达式 |expr|
static std::shared_ptr<SymbolicExpr> calculus_utils_make_abs(
    const std::shared_ptr<SymbolicExpr>& expr)
{
    if (!expr || !lamina::detail::node(expr)) return nullptr;
    return lamina::detail::make_expression_ptr(
        lamina::detail::make_node<FunctionNode>(
            FunctionNode::FuncType::Abs,
            std::vector<std::shared_ptr<const SymbolicNode>>{lamina::detail::node(expr)}));
}

static bool calculus_utils_contains_unevaluated_integral(
    const std::shared_ptr<const SymbolicNode>& node,
    std::size_t depth = 0)
{
    if (!node || depth > 200) return false;
    if (auto func = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        if (func->type() == FunctionNode::FuncType::Calculus_Integral) return true;
        for (const auto& arg : func->arguments()) {
            if (calculus_utils_contains_unevaluated_integral(arg, depth + 1)) {
                return true;
            }
        }
        return false;
    }
    if (auto add = std::dynamic_pointer_cast<const AddNode>(node)) {
        for (const auto& op : add->operands()) {
            if (calculus_utils_contains_unevaluated_integral(op, depth + 1)) {
                return true;
            }
        }
        return false;
    }
    if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        for (const auto& op : mul->operands()) {
            if (calculus_utils_contains_unevaluated_integral(op, depth + 1)) {
                return true;
            }
        }
        return false;
    }
    if (auto pow = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return calculus_utils_contains_unevaluated_integral(pow->base(), depth + 1) ||
               calculus_utils_contains_unevaluated_integral(pow->exponent(), depth + 1);
    }
    if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(node)) {
        return calculus_utils_contains_unevaluated_integral(complex->real(), depth + 1) ||
               calculus_utils_contains_unevaluated_integral(complex->imag(), depth + 1);
    }
    return false;
}

/// 尝试符号定积分，若结果仍含未求值积分节点则返回 nullptr
static std::shared_ptr<SymbolicExpr> calculus_utils_try_symbolic_definite(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    Integrator integrator;
    SymbolicExpr result = integrator.integrate_def(*integrand, var, *a, *b);

    if (calculus_utils_contains_unevaluated_integral(lamina::detail::node(result))) return nullptr;
    auto res = lamina::detail::make_expression_ptr(result);
    auto simplified = res->simplify();
    if (simplified &&
        calculus_utils_contains_unevaluated_integral(lamina::detail::node(simplified))) {
        return nullptr;
    }
    return simplified ? simplified : res;
}

/// 数值定积分回退（复合 Simpson 法）
static std::shared_ptr<SymbolicExpr> calculus_utils_numerical_definite(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    ComputationContext context;
    auto a_eval = evaluate_numeric(*a, NumericBindings{}, context);
    auto b_eval = evaluate_numeric(*b, NumericBindings{}, context);
    if (!a_eval || !b_eval ||
        !a_eval.value().is_finite() || !b_eval.value().is_finite() ||
        !std::isfinite(a_eval.value().value) ||
        !std::isfinite(b_eval.value().value)) {
        return nullptr;
    }
    double a_val = a_eval.value().value;
    double b_val = b_eval.value().value;

    int n = 1000;
    double h = (b_val - a_val) / n;
    if (!std::isfinite(h)) return nullptr;
    double sum = 0.0;

    for (int i = 0; i <= n; ++i) {
        double xi = a_val + i * h;
        auto xi_expr = SymbolicExpr::number(xi);
        auto fi = integrand->substitute(var, xi_expr);
        if (!fi || !lamina::detail::node(fi)) return nullptr;
        auto fi_eval = evaluate_numeric(*fi, NumericBindings{}, context);
        if (!fi_eval || !fi_eval.value().is_finite() ||
            !std::isfinite(fi_eval.value().value)) {
            return nullptr;
        }
        double fi_val = fi_eval.value().value;

        double weight = 1.0;
        if (i == 0 || i == n) {
            weight = 1.0;
        } else if (i % 2 == 1) {
            weight = 4.0;
        } else {
            weight = 2.0;
        }
        sum += weight * fi_val;
        if (!std::isfinite(sum)) return nullptr;
    }
    sum *= h / 3.0;
    if (!std::isfinite(sum)) return nullptr;

    return SymbolicExpr::number(sum);
}

/// 符号积分优先，失败时回退到数值积分
static std::shared_ptr<SymbolicExpr> calculus_utils_integrate_with_fallback(
    const std::shared_ptr<SymbolicExpr>& integrand,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b)
{
    auto symbolic_result = calculus_utils_try_symbolic_definite(integrand, var, a, b);
    if (symbolic_result) {
        return symbolic_result;
    }
    return calculus_utils_numerical_definite(integrand, var, a, b);
}


SymbolicExprResult curvature_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    ComputationContext& context)
{
    const std::string operation = "curvature";
    auto input = calculus_utils_validate_expr(f, var, context, operation);
    if (!input) return SymbolicExprResult::failure(input.error());

    /// f' = df/dvar
    auto f_prime = f->differentiate(var);
    if (!f_prime || !lamina::detail::node(f_prime)) {
        return SymbolicExprResult::failure(CasErrc::UnsupportedExpression,
                                           "first derivative could not be constructed",
                                           operation);
    }
    /// f'' = d²f/dvar²
    auto f_double_prime = f_prime->differentiate(var);
    if (!f_double_prime || !lamina::detail::node(f_double_prime)) {
        return SymbolicExprResult::failure(CasErrc::UnsupportedExpression,
                                           "second derivative could not be constructed",
                                           operation);
    }

    /// |f''|
    auto abs_f_pp = calculus_utils_make_abs(f_double_prime);
    if (!abs_f_pp || !lamina::detail::node(abs_f_pp)) {
        return SymbolicExprResult::failure(CasErrc::InternalInvariant,
                                           "absolute value node construction failed",
                                           operation);
    }

    /// 1 + f'²
    auto f_prime_sq = SymbolicExpr::power(f_prime, SymbolicExpr::number(2));
    auto one_plus_fp_sq = SymbolicExpr::add(SymbolicExpr::number(1), f_prime_sq);

    /// (1 + f'²)^(3/2)
    auto three_half = SymbolicExpr::number(Rational(3, 2));
    auto denom = SymbolicExpr::power(one_plus_fp_sq, three_half);

    /// κ = |f''| / (1 + f'²)^(3/2)
    auto result = SymbolicExpr::divide(abs_f_pp, denom);
    auto simplified = result->simplify();
    return SymbolicExprResult::success(simplified ? simplified : result);
}

SymbolicExprResult curvature_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    ComputationContext context;
    return curvature_checked(f, var, context);
}

std::shared_ptr<SymbolicExpr> curvature(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    auto checked = curvature_checked(f, var);
    if (!checked) return nullptr;
    return checked.value();
}

SymbolicExprResult curvature_parametric_checked(
    const std::shared_ptr<SymbolicExpr>& x_t,
    const std::shared_ptr<SymbolicExpr>& y_t, const std::string& t,
    ComputationContext& context)
{
    const std::string operation = "curvature_parametric";
    auto input = calculus_utils_validate_two_exprs(x_t, y_t, t, context, operation);
    if (!input) return SymbolicExprResult::failure(input.error());

    /// x' = dx/dt, x'' = d²x/dt²
    auto x_prime = x_t->differentiate(t);
    if (!x_prime || !lamina::detail::node(x_prime)) {
        return SymbolicExprResult::failure(CasErrc::UnsupportedExpression,
                                           "x derivative could not be constructed",
                                           operation);
    }
    auto x_double_prime = x_prime->differentiate(t);
    if (!x_double_prime || !lamina::detail::node(x_double_prime)) {
        return SymbolicExprResult::failure(CasErrc::UnsupportedExpression,
                                           "x second derivative could not be constructed",
                                           operation);
    }

    /// y' = dy/dt, y'' = d²y/dt²
    auto y_prime = y_t->differentiate(t);
    if (!y_prime || !lamina::detail::node(y_prime)) {
        return SymbolicExprResult::failure(CasErrc::UnsupportedExpression,
                                           "y derivative could not be constructed",
                                           operation);
    }
    auto y_double_prime = y_prime->differentiate(t);
    if (!y_double_prime || !lamina::detail::node(y_double_prime)) {
        return SymbolicExprResult::failure(CasErrc::UnsupportedExpression,
                                           "y second derivative could not be constructed",
                                           operation);
    }

    if (x_prime->is_zero() && y_prime->is_zero()) {
        return SymbolicExprResult::failure(CasErrc::DomainError,
                                           "parametric curvature is undefined for zero velocity",
                                           operation);
    }

    /// x'y'' - y'x''
    auto cross_term = SymbolicExpr::add(
        SymbolicExpr::multiply(x_prime, y_double_prime),
        SymbolicExpr::multiply(
            SymbolicExpr::number(-1),
            SymbolicExpr::multiply(y_prime, x_double_prime)));

    /// |x'y'' - y'x''|
    auto abs_cross = calculus_utils_make_abs(cross_term);
    if (!abs_cross || !lamina::detail::node(abs_cross)) {
        return SymbolicExprResult::failure(CasErrc::InternalInvariant,
                                           "absolute value node construction failed",
                                           operation);
    }

    /// x'² + y'²
    auto x_prime_sq = SymbolicExpr::power(x_prime, SymbolicExpr::number(2));
    auto y_prime_sq = SymbolicExpr::power(y_prime, SymbolicExpr::number(2));
    auto sum_sq = SymbolicExpr::add(x_prime_sq, y_prime_sq);

    /// (x'² + y'²)^(3/2)
    auto three_half = SymbolicExpr::number(Rational(3, 2));
    auto denom = SymbolicExpr::power(sum_sq, three_half);

    /// κ = |x'y'' - y'x''| / (x'² + y'²)^(3/2)
    auto result = SymbolicExpr::divide(abs_cross, denom);
    auto simplified = result->simplify();
    return SymbolicExprResult::success(simplified ? simplified : result);
}

SymbolicExprResult curvature_parametric_checked(
    const std::shared_ptr<SymbolicExpr>& x_t,
    const std::shared_ptr<SymbolicExpr>& y_t, const std::string& t)
{
    ComputationContext context;
    return curvature_parametric_checked(x_t, y_t, t, context);
}

std::shared_ptr<SymbolicExpr> curvature_parametric(
    const std::shared_ptr<SymbolicExpr>& x_t,
    const std::shared_ptr<SymbolicExpr>& y_t, const std::string& t)
{
    auto checked = curvature_parametric_checked(x_t, y_t, t);
    if (!checked) return nullptr;
    return checked.value();
}


SymbolicExprVectorResult inflection_points_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    ComputationContext& context)
{
    const std::string operation = "inflection_points";
    auto input = calculus_utils_validate_expr(f, var, context, operation);
    if (!input) return SymbolicExprVectorResult::failure(input.error());
    
    /// f' = df/dvar
    auto f_prime = f->differentiate(var);
    if (!f_prime || !lamina::detail::node(f_prime)) {
        return SymbolicExprVectorResult::failure(
            CasErrc::UnsupportedExpression,
            "first derivative could not be constructed", operation);
    }
    
    /// f'' = d²f/dvar²
    auto f_double_prime = f_prime->differentiate(var);
    if (!f_double_prime || !lamina::detail::node(f_double_prime)) {
        return SymbolicExprVectorResult::failure(
            CasErrc::UnsupportedExpression,
            "second derivative could not be constructed", operation);
    }

    auto equation = f_double_prime->simplify();
    if (!equation || !lamina::detail::node(equation)) equation = f_double_prime;

    auto solved = solve_dispatch_checked(equation, var, context, SolveOptions{});
    if (!solved) return SymbolicExprVectorResult::failure(solved.error());

    const auto& solutions = solved.value();
    if (solutions.kind() == SolutionSet::Kind::Empty) {
        return SymbolicExprVectorResult::success({});
    }
    if (solutions.kind() == SolutionSet::Kind::Finite) {
        std::vector<std::shared_ptr<SymbolicExpr>> points;
        points.reserve(solutions.finite_solutions().size());
        for (const auto& solution : solutions.finite_solutions()) {
            if (!solution.value || !lamina::detail::node(solution.value)) {
                return SymbolicExprVectorResult::failure(
                    CasErrc::InternalInvariant,
                    "checked solver returned a null inflection candidate",
                    operation);
            }
            points.push_back(solution.value);
        }
        return SymbolicExprVectorResult::success(std::move(points));
    }

    if (solutions.kind() == SolutionSet::Kind::Inconclusive) {
        auto legacy_candidates = SymbolicExpr::solve(equation, var);
        std::vector<std::shared_ptr<SymbolicExpr>> verified;
        verified.reserve(legacy_candidates.size());
        for (const auto& candidate : legacy_candidates) {
            if (!candidate || !lamina::detail::node(candidate)) continue;
            auto residual = equation->substitute(var, candidate);
            auto simplified = residual ? residual->simplify() : nullptr;
            if (simplified && lamina::detail::node(simplified) && simplified->is_zero()) {
                verified.push_back(candidate);
            }
        }
        if (!verified.empty()) {
            return SymbolicExprVectorResult::success(std::move(verified));
        }
    }

    return SymbolicExprVectorResult::failure(
        CasErrc::Inconclusive,
        "inflection equation is outside the finite exact support domain",
        operation);
}

SymbolicExprVectorResult inflection_points_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    ComputationContext context;
    return inflection_points_checked(f, var, context);
}

std::vector<std::shared_ptr<SymbolicExpr>> inflection_points(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var)
{
    auto checked = inflection_points_checked(f, var);
    if (!checked) return {};
    return std::move(checked.value());
}


SymbolicExprResult surface_area_revolution_x_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context)
{
    const std::string operation = "surface_area_revolution_x";
    auto input = calculus_utils_validate_expr_bounds(f, var, a, b, context, operation);
    if (!input) return SymbolicExprResult::failure(input.error());
    auto step = context.consume_steps(8, operation);
    if (!step) return SymbolicExprResult::failure(step.error());

    try {
        auto f_prime = f->differentiate(var);
        if (!f_prime || !lamina::detail::node(f_prime)) {
            return SymbolicExprResult::failure(
                CasErrc::Inconclusive,
                "surface derivative is outside the supported domain",
                operation);
        }

        auto f_prime_sq = SymbolicExpr::power(f_prime, SymbolicExpr::number(2));
        auto one_plus_fp_sq = SymbolicExpr::add(SymbolicExpr::number(1), f_prime_sq);
        auto arc_factor = SymbolicExpr::sqrt(one_plus_fp_sq);
        auto abs_f = calculus_utils_make_abs(f);
        if (!abs_f || !lamina::detail::node(abs_f)) {
            return SymbolicExprResult::failure(
                CasErrc::InternalInvariant,
                "absolute value node construction failed",
                operation);
        }
        auto integrand = SymbolicExpr::multiply(abs_f, arc_factor);

        auto integral = calculus_utils_try_symbolic_definite(integrand, var, a, b);
        if (!integral || !lamina::detail::node(integral)) {
            return SymbolicExprResult::failure(
                CasErrc::Inconclusive,
                "surface area integral could not be evaluated exactly",
                operation);
        }

        auto two_pi = SymbolicExpr::multiply(
            SymbolicExpr::number(2), SymbolicExpr::variable("pi"));
        auto result = SymbolicExpr::multiply(two_pi, integral);
        auto simplified = result->simplify();
        return SymbolicExprResult::success(simplified ? simplified : result);
    } catch (const std::bad_alloc&) {
        return SymbolicExprResult::failure(CasErrc::ResourceLimit,
                                           "surface area allocation failed",
                                           operation);
    } catch (const std::exception& e) {
        return SymbolicExprResult::failure(CasErrc::InternalInvariant,
                                           e.what(), operation);
    }
}

SymbolicExprResult surface_area_revolution_x_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b)
{
    ComputationContext context;
    return surface_area_revolution_x_checked(f, var, a, b, context);
}

std::shared_ptr<SymbolicExpr> surface_area_revolution_x(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b)
{
    /// f' = df/dvar
    auto f_prime = f->differentiate(var);

    /// √(1 + f'²)
    auto f_prime_sq = SymbolicExpr::power(f_prime, SymbolicExpr::number(2));
    auto one_plus_fp_sq = SymbolicExpr::add(SymbolicExpr::number(1), f_prime_sq);
    auto arc_factor = SymbolicExpr::sqrt(one_plus_fp_sq);

    /// |f(x)| · √(1 + f'²)
    auto abs_f = calculus_utils_make_abs(f);
    auto integrand = SymbolicExpr::multiply(abs_f, arc_factor);

    /// 2π
    auto two_pi = SymbolicExpr::multiply(
        SymbolicExpr::number(2), SymbolicExpr::variable("pi"));

    /// S = 2π ∫ₐᵇ |f(x)| · √(1 + f'²) dx
    auto integral = calculus_utils_integrate_with_fallback(integrand, var, a, b);
    if (!integral) {
        return nullptr;
    }
    auto result = SymbolicExpr::multiply(two_pi, integral);
    auto simplified = result->simplify();
    return simplified ? simplified : result;
}

SymbolicExprResult surface_area_revolution_y_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context)
{
    const std::string operation = "surface_area_revolution_y";
    auto input = calculus_utils_validate_expr_bounds(f, var, a, b, context, operation);
    if (!input) return SymbolicExprResult::failure(input.error());
    auto step = context.consume_steps(8, operation);
    if (!step) return SymbolicExprResult::failure(step.error());

    try {
        auto f_prime = f->differentiate(var);
        if (!f_prime || !lamina::detail::node(f_prime)) {
            return SymbolicExprResult::failure(
                CasErrc::Inconclusive,
                "surface derivative is outside the supported domain",
                operation);
        }

        auto f_prime_sq = SymbolicExpr::power(f_prime, SymbolicExpr::number(2));
        auto one_plus_fp_sq = SymbolicExpr::add(SymbolicExpr::number(1), f_prime_sq);
        auto arc_factor = SymbolicExpr::sqrt(one_plus_fp_sq);
        auto var_expr = SymbolicExpr::variable(var);
        auto abs_var = calculus_utils_make_abs(var_expr);
        if (!abs_var || !lamina::detail::node(abs_var)) {
            return SymbolicExprResult::failure(
                CasErrc::InternalInvariant,
                "absolute value node construction failed",
                operation);
        }
        auto integrand = SymbolicExpr::multiply(abs_var, arc_factor);

        auto integral = calculus_utils_try_symbolic_definite(integrand, var, a, b);
        if (!integral || !lamina::detail::node(integral)) {
            return SymbolicExprResult::failure(
                CasErrc::Inconclusive,
                "surface area integral could not be evaluated exactly",
                operation);
        }

        auto two_pi = SymbolicExpr::multiply(
            SymbolicExpr::number(2), SymbolicExpr::variable("pi"));
        auto result = SymbolicExpr::multiply(two_pi, integral);
        auto simplified = result->simplify();
        return SymbolicExprResult::success(simplified ? simplified : result);
    } catch (const std::bad_alloc&) {
        return SymbolicExprResult::failure(CasErrc::ResourceLimit,
                                           "surface area allocation failed",
                                           operation);
    } catch (const std::exception& e) {
        return SymbolicExprResult::failure(CasErrc::InternalInvariant,
                                           e.what(), operation);
    }
}

SymbolicExprResult surface_area_revolution_y_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b)
{
    ComputationContext context;
    return surface_area_revolution_y_checked(f, var, a, b, context);
}

std::shared_ptr<SymbolicExpr> surface_area_revolution_y(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a, const std::shared_ptr<SymbolicExpr>& b)
{
    /// f' = df/dvar
    auto f_prime = f->differentiate(var);

    /// √(1 + f'²)
    auto f_prime_sq = SymbolicExpr::power(f_prime, SymbolicExpr::number(2));
    auto one_plus_fp_sq = SymbolicExpr::add(SymbolicExpr::number(1), f_prime_sq);
    auto arc_factor = SymbolicExpr::sqrt(one_plus_fp_sq);

    /// |var| · √(1 + f'²)
    auto var_expr = SymbolicExpr::variable(var);
    auto abs_var = calculus_utils_make_abs(var_expr);
    auto integrand = SymbolicExpr::multiply(abs_var, arc_factor);

    /// 2π
    auto two_pi = SymbolicExpr::multiply(
        SymbolicExpr::number(2), SymbolicExpr::variable("pi"));

    /// S = 2π ∫ₐᵇ |x| · √(1 + f'²) dx
    auto integral = calculus_utils_integrate_with_fallback(integrand, var, a, b);
    if (!integral) {
        return nullptr;
    }
    auto result = SymbolicExpr::multiply(two_pi, integral);
    auto simplified = result->simplify();
    return simplified ? simplified : result;
}

} // namespace lamina
