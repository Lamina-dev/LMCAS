/**
 * @file series_engine.cpp
 */

#include "series_engine.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "integration.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/series_support.hpp"
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <cstdio>
#include <typeinfo>
#include <algorithm>
#include <variant>
#include <limits>
#include <optional>


namespace LMCAS {

static bool try_get_int(const std::shared_ptr<SymbolicExpr>& expr, long long& out) {
    if (!expr || !LMCAS::detail::node(expr)) return false;
    auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr));
    if (!num) return false;
    if (std::holds_alternative<BigInt>(num->value())) {
        auto value = std::get<BigInt>(num->value()).try_to_int64();
        if (!value) return false;
        out = static_cast<long long>(*value);
        return true;
    }
    if (std::holds_alternative<Rational>(num->value())) {
        auto& r = std::get<Rational>(num->value());
        if (r.get_denominator() == BigInt(1)) {
            auto value = r.get_numerator().try_to_int64();
            if (!value) return false;
            out = static_cast<long long>(*value);
            return true;
        }
        return false;
    }
    double d = static_cast<double>(std::get<lmmc_real_t>(num->value()));
    if (d == std::floor(d) && std::abs(d) < 1e15) { out = static_cast<long long>(d); return true; }
    return false;
}

} // namespace LMCAS




namespace LMCAS {

namespace {


Result<void> validate_power_series_order(
    int order,
    ComputationContext& context,
    const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (order <= 0) {
        return Result<void>::failure(
            CasErrc::InvalidArgument,
            "power series truncation order must be positive",
            operation);
    }
    return Result<void>::success();
}




} // namespace
namespace detail::series_support {

Result<void> validate_series_variable(const std::string& var,
                                      ComputationContext& context,
                                      const std::string& operation)
{
    auto step = context.consume_steps(1, operation);
    if (!step) return step;
    if (var.empty()) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "series variable name cannot be empty",
                                     operation);
    }
    return Result<void>::success();
}

std::optional<int> integer_value_from_node(const std::shared_ptr<const SymbolicNode>& node)
{
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return std::nullopt;
    if (std::holds_alternative<BigInt>(number->value())) {
        auto value = std::get<BigInt>(number->value()).try_to_int64();
        if (!value || *value < std::numeric_limits<int>::min() ||
            *value > std::numeric_limits<int>::max()) return std::nullopt;
        return static_cast<int>(*value);
    }
    if (std::holds_alternative<Rational>(number->value())) {
        const auto& rational = std::get<Rational>(number->value());
        if (rational.get_denominator() != BigInt(1)) return std::nullopt;
        auto value = rational.get_numerator().try_to_int64();
        if (!value || *value < std::numeric_limits<int>::min() ||
            *value > std::numeric_limits<int>::max()) return std::nullopt;
        return static_cast<int>(*value);
    }
    double value = static_cast<double>(std::get<lmmc_real_t>(number->value()));
    if (!std::isfinite(value) || value != std::floor(value) ||
        std::abs(value) > static_cast<double>(std::numeric_limits<int>::max())) {
        return std::nullopt;
    }
    return static_cast<int>(value);
}

std::optional<int> supported_laurent_integer_power(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& var)
{
    if (!node) return std::nullopt;
    if (std::dynamic_pointer_cast<const NumberNode>(node)) return 0;
    if (auto variable = std::dynamic_pointer_cast<const VariableNode>(node)) {
        return variable->name() == var ? std::optional<int>(1) : std::nullopt;
    }
    if (auto power = std::dynamic_pointer_cast<const PowerNode>(node)) {
        auto exponent = integer_value_from_node(power->exponent());
        if (!exponent) return std::nullopt;
        auto base = std::dynamic_pointer_cast<const VariableNode>(power->base());
        if (base && base->name() == var) return *exponent;
        auto nested = supported_laurent_integer_power(power->base(), var);
        if (!nested) return std::nullopt;
        const long long product = static_cast<long long>(*nested) * *exponent;
        if (product < std::numeric_limits<int>::min() ||
            product > std::numeric_limits<int>::max()) return std::nullopt;
        return static_cast<int>(product);
    }
    if (auto multiply = std::dynamic_pointer_cast<const MultiplyNode>(node)) {
        int total_power = 0;
        for (const auto& operand : multiply->operands()) {
            if (std::dynamic_pointer_cast<const NumberNode>(operand)) continue;
            auto power = supported_laurent_integer_power(operand, var);
            if (!power) return std::nullopt;
            total_power += *power;
        }
        return total_power;
    }
    return std::nullopt;


}

Result<void> validate_power_series_coefficients(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coeffs,
    const std::string& operation,
    const std::string& name)
{
    for (size_t i = 0; i < coeffs.size(); ++i) {
        if (!coeffs[i] || !LMCAS::detail::node(coeffs[i])) {
            return Result<void>::failure(
                CasErrc::InvalidArgument,
                name + " contains a null coefficient at index " + std::to_string(i),
                operation);
        }
    }
    return Result<void>::success();
}
} // namespace detail::series_support

using detail::series_support::validate_power_series_coefficients;

std::vector<std::shared_ptr<SymbolicExpr>> power_series_add(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b) {
    size_t len = std::max(a.size(), b.size());
    std::vector<std::shared_ptr<SymbolicExpr>> result(len);
    for (size_t i = 0; i < len; ++i) {
        auto ai = (i < a.size() && a[i]) ? a[i] : SymbolicExpr::number(0);
        auto bi = (i < b.size() && b[i]) ? b[i] : SymbolicExpr::number(0);
        result[i] = SymbolicExpr::add(ai, bi)->simplify();
    }
    return result;
}

PowerSeriesResult power_series_multiply_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b, int order,
    ComputationContext& context) {
    const std::string operation = "power_series_multiply";
    auto order_check = validate_power_series_order(order, context, operation);
    if (!order_check) return PowerSeriesResult::failure(order_check.error());
    auto a_check = validate_power_series_coefficients(a, operation, "left series");
    if (!a_check) return PowerSeriesResult::failure(a_check.error());
    auto b_check = validate_power_series_coefficients(b, operation, "right series");
    if (!b_check) return PowerSeriesResult::failure(b_check.error());

    size_t n = static_cast<size_t>(order);
    std::vector<std::shared_ptr<SymbolicExpr>> result(n, SymbolicExpr::number(0));
    for (size_t k = 0; k < n; ++k) {
        auto step = context.consume_steps(1, operation);
        if (!step) return PowerSeriesResult::failure(step.error());
        std::vector<std::shared_ptr<const SymbolicNode>> terms;
        for (size_t j = 0; j <= k; ++j) {
            if (j >= a.size() || (k-j) >= b.size()) continue;
            auto aj = a[j] ? a[j] : SymbolicExpr::number(0);
            auto bkj = b[k-j] ? b[k-j] : SymbolicExpr::number(0);
            terms.push_back(LMCAS::detail::node(SymbolicExpr::multiply(aj, bkj)));
        }
        if (terms.empty()) result[k] = SymbolicExpr::number(0);
        else if (terms.size() == 1) result[k] = LMCAS::detail::make_expression_ptr(terms[0])->simplify();
        else result[k] = LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<AddNode>(terms))->simplify();
    }
    return PowerSeriesResult::success(std::move(result));
}

PowerSeriesResult power_series_multiply_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& a,
    const std::vector<std::shared_ptr<SymbolicExpr>>& b, int order) {
    ComputationContext context;
    return power_series_multiply_checked(a, b, order, context);
}


PowerSeriesResult power_series_compose_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& g, int order,
    ComputationContext& context) {
    const std::string operation = "power_series_compose";
    auto order_check = validate_power_series_order(order, context, operation);
    if (!order_check) return PowerSeriesResult::failure(order_check.error());
    auto f_check = validate_power_series_coefficients(f, operation, "outer series");
    if (!f_check) return PowerSeriesResult::failure(f_check.error());
    auto g_check = validate_power_series_coefficients(g, operation, "inner series");
    if (!g_check) return PowerSeriesResult::failure(g_check.error());

    size_t n = static_cast<size_t>(order);
    if (!g.empty() && !g[0]->is_zero()) {
        return PowerSeriesResult::failure(
            CasErrc::DomainError,
            "power series composition around zero requires g(0) = 0",
            operation);
    }
    std::vector<std::shared_ptr<SymbolicExpr>> result(n, SymbolicExpr::number(0));
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> gp(n);
    gp[0].resize(n, SymbolicExpr::number(0)); gp[0][0] = SymbolicExpr::number(1);
    if (n > 1) { gp[1].resize(n, SymbolicExpr::number(0)); for (size_t i = 0; i < std::min(g.size(), n); ++i) gp[1][i] = g[i] ? g[i] : SymbolicExpr::number(0); }
    for (size_t k = 2; k < n; ++k) {
        auto step = context.consume_steps(1, operation);
        if (!step) return PowerSeriesResult::failure(step.error());
        auto multiplied = power_series_multiply_checked(gp[k-1], gp[1], order, context);
        if (!multiplied) return multiplied;
        gp[k] = std::move(multiplied.value());
    }
    for (size_t k = 0; k < std::min(f.size(), n); ++k) {
        auto step = context.consume_steps(1, operation);
        if (!step) return PowerSeriesResult::failure(step.error());
        auto fk = f[k] ? f[k] : SymbolicExpr::number(0);
        if (fk->is_zero()) continue;
        for (size_t i = 0; i < n && i < gp[k].size(); ++i)
            result[i] = SymbolicExpr::add(result[i], SymbolicExpr::multiply(fk, gp[k][i]))->simplify();
    }
    return PowerSeriesResult::success(std::move(result));
}

PowerSeriesResult power_series_compose_checked(
    const std::vector<std::shared_ptr<SymbolicExpr>>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& g, int order) {
    ComputationContext context;
    return power_series_compose_checked(f, g, order, context);
}



std::shared_ptr<SymbolicExpr> symbolic_sum(
    const std::shared_ptr<SymbolicExpr>& body, const std::string& index,
    const std::shared_ptr<SymbolicExpr>& lower, const std::shared_ptr<SymbolicExpr>& upper) {
    if (!body || !lower || !upper) return nullptr;
    long long lv = 0, uv = 0;
    if (try_get_int(lower, lv) && try_get_int(upper, uv)) {
        if (uv < lv) return SymbolicExpr::number(0);
        const auto span = static_cast<std::uint64_t>(uv) -
                          static_cast<std::uint64_t>(lv);
        if (span < 100) {
            std::vector<std::shared_ptr<const SymbolicNode>> terms;
            for (long long k = lv;; ++k) {
                auto val = SymbolicExpr::number(k);
                auto term = body->substitute(index, val);
                if (term) { term = term->simplify(); terms.push_back(LMCAS::detail::node(term)); }
                if (k == uv) break;
            }
            if (terms.empty()) return SymbolicExpr::number(0);
            if (terms.size() == 1) return LMCAS::detail::make_expression_ptr(terms[0])->simplify();
            return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<AddNode>(terms))->simplify();
        }
    }
    auto k_var = SymbolicExpr::variable(index);
    if (body->to_string() == k_var->to_string()) {
        auto u = upper; auto l = lower;
        auto u1 = SymbolicExpr::add(u, SymbolicExpr::number(1));
        auto su = SymbolicExpr::divide(SymbolicExpr::multiply(u, u1), SymbolicExpr::number(2));
        auto lm1 = SymbolicExpr::add(l, SymbolicExpr::number(-1));
        auto sl = SymbolicExpr::divide(SymbolicExpr::multiply(lm1, l), SymbolicExpr::number(2));
        return SymbolicExpr::add(su, SymbolicExpr::multiply(SymbolicExpr::number(-1), sl))->simplify();
    }
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<SummationNode>(LMCAS::detail::node(body), index, LMCAS::detail::node(lower), LMCAS::detail::node(upper)));
}

std::shared_ptr<SymbolicExpr> symbolic_product(
    const std::shared_ptr<SymbolicExpr>& body, const std::string& index,
    const std::shared_ptr<SymbolicExpr>& lower, const std::shared_ptr<SymbolicExpr>& upper) {
    if (!body || !lower || !upper) return nullptr;
    long long lv = 0, uv = 0;
    if (try_get_int(lower, lv) && try_get_int(upper, uv)) {
        if (uv < lv) return SymbolicExpr::number(1);
        const auto span = static_cast<std::uint64_t>(uv) -
                          static_cast<std::uint64_t>(lv);
        if (span < 50) {
            auto result = SymbolicExpr::number(1);
            for (long long k = lv;; ++k) {
                auto val = SymbolicExpr::number(k);
                auto term = body->substitute(index, val);
                if (term) { term = term->simplify(); result = SymbolicExpr::multiply(result, term); }
                if (k == uv) break;
            }
            return result->simplify();
        }
    }
    return LMCAS::detail::make_expression_ptr(LMCAS::detail::make_node<ProductNode>(LMCAS::detail::node(body), index, LMCAS::detail::node(lower), LMCAS::detail::node(upper)));
}



} // namespace LMCAS
