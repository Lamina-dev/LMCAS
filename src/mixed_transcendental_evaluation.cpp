#include "internal/mixed_transcendental_support.hpp"

#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>

namespace LMCAS::detail {
namespace {

lmmc_real_t unsupported_value() {
    return std::numeric_limits<lmmc_real_t>::quiet_NaN();
}

lmmc_real_t evaluate_number(const NumberNode& number) {
    if (std::holds_alternative<lmmc_real_t>(number.value())) {
        return std::get<lmmc_real_t>(number.value());
    }
    if (std::holds_alternative<BigInt>(number.value())) {
        return static_cast<lmmc_real_t>(std::get<BigInt>(number.value()).to_double());
    }
    if (std::holds_alternative<Rational>(number.value())) {
        return static_cast<lmmc_real_t>(std::get<Rational>(number.value()).to_double());
    }
    return 0.0;
}

lmmc_real_t evaluate_inverse_trig(FunctionNode::FuncType type, lmmc_real_t arg) {
    using Type = FunctionNode::FuncType;
    switch (type) {
    case Type::ArcSin:
        if (arg < -1.0 || arg > 1.0) { return unsupported_value(); }
        return std::asin(arg);
    case Type::ArcCos:
        if (arg < -1.0 || arg > 1.0) { return unsupported_value(); }
        return std::acos(arg);
    case Type::ArcTan:
        return std::atan(arg);
    case Type::Sec: {
        const auto cosine = std::cos(arg);
        if (std::fabs(cosine) < 1e-15) { return unsupported_value(); }
        return 1.0 / cosine;
    }
    case Type::Csc:
    case Type::Cot: {
        const auto sine = std::sin(arg);
        if (std::fabs(sine) < 1e-15) { return unsupported_value(); }
        return type == Type::Csc ? 1.0 / sine : std::cos(arg) / sine;
    }
    default:
        return unsupported_value();
    }
}

lmmc_real_t evaluate_unary(FunctionNode::FuncType type, lmmc_real_t arg) {
    using Type = FunctionNode::FuncType;
    switch (type) {
    case Type::Sin: {
        return std::sin(arg);
    }
    case Type::Cos: {
        return std::cos(arg);
    }
    case Type::Tan: {
        return std::tan(arg);
    }
    case Type::Exp: {
        return std::exp(arg);
    }
    case Type::Ln: {
        return std::log(arg);
    }
    case Type::Sqrt: {
        return std::sqrt(arg);
    }
    case Type::Abs: {
        return std::abs(arg);
    }
    case Type::Sinh: {
        return std::sinh(arg);
    }
    case Type::Cosh: {
        return std::cosh(arg);
    }
    case Type::Tanh: {
        return std::tanh(arg);
    }
    default: return evaluate_inverse_trig(type, arg);
    }
}

lmmc_real_t evaluate_function(const FunctionNode& function) {
    if (function.arguments().size() == 1) {
        const auto argument = mixed_recursive_eval(function.arguments()[0]);
        if (std::isnan(argument)) { return argument; }
        return evaluate_unary(function.type(), argument);
    }
    if (function.arguments().size() == 2 &&
        function.type() == FunctionNode::FuncType::Atan2) {
        const auto y = mixed_recursive_eval(function.arguments()[0]);
        const auto x = mixed_recursive_eval(function.arguments()[1]);
        if (std::isnan(y) || std::isnan(x)) { return unsupported_value(); }
        return std::atan2(y, x);
    }
    return unsupported_value();
}

lmmc_real_t evaluate_sum(const AddNode& sum) {
    lmmc_real_t value = 0.0;
    for (const auto& operand : sum.operands()) {
        const auto term = mixed_recursive_eval(operand);
        if (std::isnan(term)) { return term; }
        value += term;
    }
    return value;
}

lmmc_real_t evaluate_product(const MultiplyNode& product) {
    lmmc_real_t value = 1.0;
    for (const auto& operand : product.operands()) {
        const auto factor = mixed_recursive_eval(operand);
        if (std::isnan(factor)) { return factor; }
        value *= factor;
    }
    return value;
}

lmmc_real_t evaluate_power(const PowerNode& power) {
    const auto base = mixed_recursive_eval(power.base());
    const auto exponent = mixed_recursive_eval(power.exponent());
    if (std::isnan(base) || std::isnan(exponent)) { return unsupported_value(); }
    if (base == 0.0 && exponent < 0.0) { return unsupported_value(); }
    return std::pow(base, exponent);
}

}

lmmc_real_t mixed_recursive_eval(const std::shared_ptr<const SymbolicNode>& node) {
    if (!node) { return 0.0; }
    if (auto number = std::dynamic_pointer_cast<const NumberNode>(node)) {
        return evaluate_number(*number);
    }
    if (auto sum = std::dynamic_pointer_cast<const AddNode>(node)) {
        return evaluate_sum(*sum);
    }
    if (auto product = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        return evaluate_product(*product);
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        return evaluate_power(*power);
    }
    if (auto function = std::dynamic_pointer_cast<const FunctionNode>(node)) {
        return evaluate_function(*function);
    }
    return unsupported_value();
}

lmmc_real_t mixed_evaluate_at(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var, lmmc_real_t x) {
    try {
        auto substituted = expr->substitute(var, SymbolicExpr::number(static_cast<double>(x)));
        if (!substituted || !node(substituted)) { return unsupported_value(); }
        const auto value = mixed_recursive_eval(node(substituted));
        if (!std::isfinite(value)) { return unsupported_value(); }
        return value;
    } catch (const std::bad_alloc&) {
        throw;
    } catch (const std::exception&) {
        return unsupported_value();
    }
}

lmmc_real_t mixed_evaluate_with_retry(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var, lmmc_real_t x,
    lmmc_real_t half_width, int max_retries) {
    auto value = mixed_evaluate_at(expr, var, x);
    if (!std::isnan(value)) { return value; }
    auto offset = half_width;
    for (int i = 0; i < max_retries; ++i) {
        offset *= 0.5;
        value = mixed_evaluate_at(expr, var, x + offset);
        if (!std::isnan(value)) { return value; }
        value = mixed_evaluate_at(expr, var, x - offset);
        if (!std::isnan(value)) { return value; }
    }
    return unsupported_value();
}

bool consume_mixed_step(
    ComputationContext* context, std::optional<CasError>* failure,
    bool* complete) {
    if (!context) { return true; }
    auto step = context->consume_steps(1, "solve_mixed_transcendental_checked");
    if (step) { return true; }
    if (failure) { *failure = step.error(); }
    if (complete) { *complete = false; }
    return false;
}

}
