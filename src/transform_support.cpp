#include "internal/transform_support.hpp"
#include <cmath>

namespace LMCAS::transform_detail {

bool te_depends_on(const std::shared_ptr<SymbolicExpr>& expression,
                          const std::string& variable) {
    return expression &&
           expression_depends_on_variable(
               LMCAS::detail::node(expression), variable);
}

Result<void> te_validate_expr_vars(const std::shared_ptr<SymbolicExpr>& expr,
                                          const std::string& input_var,
                                          const std::string& output_var,
                                          ComputationContext& context,
                                          const std::string& operation) {
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (!expr || !LMCAS::detail::node(expr)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "transform expression cannot be null",
                                     operation);
    }
    if (input_var.empty() || output_var.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "transform variable names cannot be empty",
                                     operation);
    }
    if (input_var == output_var) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "input and output variables must be distinct",
                                     operation);
    }
    return Result<void>::success();
}
bool te_contains_transform(
    const std::shared_ptr<const SymbolicNode>& node) {
    return LMCAS::detail::contains_node_type<TransformNode>(node);
}

TransformEngineResult te_wrap_transform_result(
    std::shared_ptr<SymbolicExpr> expr,
    const std::string& operation) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return TransformEngineResult::failure(
            CasErrc::InternalInvariant,
            "transform construction produced a null expression",
            operation);
    }

    if (te_contains_transform(LMCAS::detail::node(expr))) {
        return TransformEngineResult::failure(
            CasErrc::Inconclusive,
            "transform is outside the current evaluated support domain",
            operation);
    }

    return TransformEngineResult::success(
        Verified<EvaluatedTransform>{
            EvaluatedTransform{std::move(expr), {}, {}},
            ByConstructionProof{}});
}

std::shared_ptr<SymbolicExpr> te_gt_condition(
    const std::shared_ptr<SymbolicExpr>& lhs,
    const std::shared_ptr<SymbolicExpr>& rhs) {
    if (!lhs || !LMCAS::detail::node(lhs) || !rhs || !LMCAS::detail::node(rhs)) return nullptr;
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<RelationalNode>(LMCAS::detail::node(lhs), LMCAS::detail::node(rhs), RelationalNode::Op::GT));
}

bool te_is_zero(const std::shared_ptr<SymbolicExpr>& e) {
    return e && LMCAS::detail::node(e) && e->is_zero();
}
std::optional<double> te_number_value(
    const std::shared_ptr<const SymbolicNode>& node) {
    const auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return std::nullopt;
    if (std::holds_alternative<BigInt>(number->value())) {
        return std::get<BigInt>(number->value()).to_double();
    }
    if (std::holds_alternative<Rational>(number->value())) {
        return std::get<Rational>(number->value()).to_double();
    }
    return static_cast<double>(std::get<lmmc_real_t>(number->value()));
}
BigInt te_factorial(const BigInt& n, ComputationContext& context) {
    auto result = BigInt::factorial_checked(n, context);
    if (!result) throw result.error();
    return std::move(result.value());
}
std::pair<std::shared_ptr<SymbolicExpr>,std::shared_ptr<SymbolicExpr>>
te_split_coeff(const std::shared_ptr<SymbolicExpr>& e, const std::string& v) {
    if (!e || !LMCAS::detail::node(e) || !te_depends_on(e, v)) return {e ? e : SymbolicExpr::number(1), SymbolicExpr::number(1)};
    auto mul = std::dynamic_pointer_cast<const MultiplyNode>(LMCAS::detail::node(e));
    if (!mul) return {SymbolicExpr::number(1), e};
    std::vector<std::shared_ptr<const SymbolicNode>> cp, vp;
    for (auto& op : mul->operands()) {
        if (expression_depends_on_variable(op, v)) vp.push_back(op); else cp.push_back(op);
    }
    if (cp.empty()) return {SymbolicExpr::number(1), e};
    if (vp.empty()) return {e, SymbolicExpr::number(1)};
    auto c = (cp.size()==1) ? LMCAS::detail::make_expression_ptr(cp[0])
        : LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<MultiplyNode>(cp));
    auto b = (vp.size()==1) ? LMCAS::detail::make_expression_ptr(vp[0])
        : LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<MultiplyNode>(vp));
    return {c, b};
}

std::shared_ptr<SymbolicExpr> te_linear_coefficient(
    const std::shared_ptr<const SymbolicNode>& argument, const std::string& variable) {
    auto symbol = std::dynamic_pointer_cast<const VariableNode>(argument);
    if (symbol && !symbol->is_constant() && symbol->name() == variable) {
        return detail::make_expression_ptr(SymbolicFactory::create_number(BigInt(1)));
    }
    auto product = std::dynamic_pointer_cast<const MultiplyNode>(argument);
    if (!product || product->operands().size() != 2) {
        return nullptr;
    }
    for (std::size_t i = 0; i < 2; ++i) {
        symbol = std::dynamic_pointer_cast<const VariableNode>(product->operands()[i]);
        const auto& other = product->operands()[1 - i];
        if (!symbol || symbol->is_constant() || symbol->name() != variable) continue;
        if (expression_depends_on_variable(other, variable)) continue;
        return detail::make_expression_ptr(other);
    }
    return nullptr;
}

namespace {

bool te_is_variable_square(
    const std::shared_ptr<const SymbolicNode>& expression,
    const std::string& variable) {
    const auto power =
        std::dynamic_pointer_cast<const PowerNode>(expression);
    if (!power) return false;
    const auto base =
        std::dynamic_pointer_cast<const VariableNode>(power->base());
    return base && !base->is_constant() && base->name() == variable &&
        exact_small_integer_node(power->exponent(), 2, 2).has_value();
}

std::optional<double> te_atomic_gaussian_rate(
    const std::shared_ptr<const SymbolicNode>& argument,
    const std::string& variable) {
    auto parts = te_split_coeff(
        detail::make_expression_ptr(argument), variable);
    const auto coefficient = te_number_value(detail::node(parts.first));
    if (!coefficient || !std::isfinite(*coefficient) ||
        !te_is_variable_square(detail::node(parts.second), variable) ||
        !(*coefficient < 0.0)) {
        return std::nullopt;
    }
    const double rate = -*coefficient;
    return std::isfinite(rate) && rate > 0.0
        ? std::optional<double>(rate)
        : std::nullopt;
}

bool te_accumulate_gaussian_factor(
    const std::shared_ptr<const SymbolicNode>& factor,
    const std::string& variable, double& coefficient,
    bool& has_variable_square) {
    if (const auto number = te_number_value(factor)) {
        coefficient *= *number;
        return true;
    }
    if (!te_is_variable_square(factor, variable)) return false;
    has_variable_square = true;
    return true;
}

std::optional<double> te_gaussian_product_rate(
    const MultiplyNode& product, const std::string& variable) {
    double coefficient = 1.0;
    bool has_variable_square = false;
    for (const auto& factor : product.operands()) {
        if (!te_accumulate_gaussian_factor(
                factor, variable, coefficient, has_variable_square)) {
            return std::nullopt;
        }
    }
    return has_variable_square && coefficient < 0.0
        ? std::optional<double>(-coefficient)
        : std::nullopt;
}
}

std::optional<double> te_match_pure_gaussian(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable) {
    if (!expression || !detail::node(expression)) return std::nullopt;
    const auto function = std::dynamic_pointer_cast<const FunctionNode>(
        detail::node(expression));
    if (!function || function->type() != FunctionNode::FuncType::Exp ||
        function->arguments().empty()) {
        return std::nullopt;
    }
    const auto product = std::dynamic_pointer_cast<const MultiplyNode>(
        function->arguments()[0]);
    if (!product) return std::nullopt;

    return te_gaussian_product_rate(*product, variable);
}

std::optional<GaussianMatch> te_match_scaled_gaussian(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable) {
    auto parts = te_split_coeff(expression, variable);
    const auto coefficient = te_number_value(detail::node(parts.first));
    if (!coefficient || !std::isfinite(*coefficient)) {
        return std::nullopt;
    }
    const auto function = std::dynamic_pointer_cast<const FunctionNode>(
        detail::node(parts.second));
    if (!function || function->type() != FunctionNode::FuncType::Exp ||
        function->arguments().size() != 1) {
        return std::nullopt;
    }
    const auto rate =
        te_atomic_gaussian_rate(function->arguments()[0], variable);
    if (!rate) return std::nullopt;
    return GaussianMatch{std::move(parts.first), *rate};
}

bool te_is_trig(const std::shared_ptr<SymbolicExpr>& expression,
                const std::string& variable, bool& is_sin,
                std::shared_ptr<SymbolicExpr>& frequency) {
    if (!expression || !detail::node(expression)) {
        return false;
    }
    auto function = std::dynamic_pointer_cast<const FunctionNode>(detail::node(expression));
    if (!function || function->arguments().empty()) {
        return false;
    }
    if (function->type() != FunctionNode::FuncType::Sin &&
        function->type() != FunctionNode::FuncType::Cos) {
        return false;
    }
    is_sin = function->type() == FunctionNode::FuncType::Sin;
    frequency = te_linear_coefficient(function->arguments()[0], variable);
    return frequency != nullptr;
}

}
