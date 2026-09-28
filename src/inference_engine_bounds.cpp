#include "internal/inference_engine_impl.hpp"
#include "internal/assumption_facts.hpp"
#include "internal/exact_constant_bounds.hpp"
#include "internal/interval_endpoint.hpp"
#include "internal/normalization_utils.hpp"
#include <array>
#include <limits>
#include <utility>

namespace LMCAS {


namespace {
constexpr const char* bounds_operation = "inference.bounds";
using RationalBounds = std::pair<Rational, Rational>;

Result<Interval> unknown_bounds() {
    return Result<Interval>::failure(CasErrc::Inconclusive,
        "a finite real enclosure is not proven", bounds_operation);
}

Interval closed_bounds(const Rational& lower, const Rational& upper) {
    return Interval{Endpoint::closed(SymbolicExpr::number(lower)),
                    Endpoint::closed(SymbolicExpr::number(upper))};
}
class BoundsPropagation {
public:
    BoundsPropagation(const AssumptionContext& assumptions, const InferenceEngine& engine)
        : assumptions_(assumptions), facts_(engine) {}

    Result<Interval> propagate(const std::shared_ptr<const SymbolicNode>& node) {
        auto entered = context_.enter_recursion(bounds_operation);
        if (!entered) { return Result<Interval>::failure(entered.error()); }
        struct Leave {
            ComputationContext& context;
            ~Leave() { context.leave_recursion(); }
        } leave{context_};
        if (!node) { return unknown_bounds(); }
        if (const auto* number = dynamic_cast<const NumberNode*>(node.get())) {
            if (const auto* value = std::get_if<lmmc_real_t>(&number->value())) {
                if (!std::isfinite(*value)) { return unknown_bounds(); }
            }
            return Result<Interval>::success(
                Interval::point(detail::make_expression_ptr(node)));
        }
        if (const auto* variable = dynamic_cast<const VariableNode*>(node.get())) {
            return variable_bounds(*variable);
        }
        if (const auto* add = dynamic_cast<const AddNode*>(node.get())) {
            return sum_bounds(*add);
        }
        if (const auto* multiply = dynamic_cast<const MultiplyNode*>(node.get())) {
            return product_bounds(*multiply);
        }
        if (const auto* power = dynamic_cast<const PowerNode*>(node.get())) {
            return power_bounds(*power);
        }
        if (const auto* function = dynamic_cast<const FunctionNode*>(node.get()))
            { return function_bounds(*function); }
        return unknown_bounds();
    }

private:
    Result<Interval> variable_bounds(const VariableNode& variable) {
        if (variable.is_constant()) { return unknown_bounds(); }
        auto interval = assumptions_.get_bounds(variable.name());
        if (!interval) { return unknown_bounds(); }
        auto lower = detail::comparable_endpoint(interval->lower, context_, bounds_operation);
        if (!lower) { return Result<Interval>::failure(lower.error()); }
        auto upper = detail::comparable_endpoint(interval->upper, context_, bounds_operation);
        if (!upper) { return Result<Interval>::failure(upper.error()); }
        auto order = detail::compare_comparable(
            lower.value(), upper.value(), context_);
        if (!order) { return Result<Interval>::failure(order.error()); }
        if (order.value() > 0 || (order.value() == 0 &&
            (interval->lower.is_open || interval->upper.is_open))) { return unknown_bounds(); }
        return Result<Interval>::success(std::move(*interval));
    }

    Result<Interval> sum_bounds(const AddNode& add) {
        RationalBounds total{Rational(0), Rational(0)};
        for (const auto& operand : add.operands()) {
            auto interval = propagate(operand);
            if (!interval) { return interval; }
            auto bounds = rational_bounds(interval.value());
            if (!bounds) { return Result<Interval>::failure(bounds.error()); }
            auto lower = arithmetic(total.first, bounds.value().first, '+');
            if (!lower) { return Result<Interval>::failure(lower.error()); }
            auto upper = arithmetic(total.second, bounds.value().second, '+');
            if (!upper) { return Result<Interval>::failure(upper.error()); }
            total = {std::move(lower.value()), std::move(upper.value())};
        }
        return Result<Interval>::success(closed_bounds(total.first, total.second));
    }

    Result<Interval> product_bounds(const MultiplyNode& multiply) {
        std::optional<Interval> numerator;
        std::optional<Interval> denominator;
        for (const auto& operand : multiply.operands()) {
            const auto* power = dynamic_cast<const PowerNode*>(operand.get());
            const auto exponent = power
                ? std::dynamic_pointer_cast<const NumberNode>(
                      power->exponent())
                : nullptr;
            const bool inverse = exponent && exponent->is_negative_one();
            auto next = propagate(inverse ? power->base() : operand);
            if (!next) { return next; }
            auto& accumulated = inverse ? denominator : numerator;
            if (!accumulated) {
                accumulated = std::move(next.value());
            } else {
                auto product = multiply_intervals(*accumulated, next.value(), false);
                if (!product) { return product; }
                accumulated = std::move(product.value());
            }
        }
        if (!numerator) { numerator = closed_bounds(Rational(1), Rational(1)); }
        if (!denominator) { return Result<Interval>::success(std::move(*numerator)); }
        return multiply_intervals(*numerator, *denominator, true);
    }

    Result<Interval> power_bounds(const PowerNode& power) {
        const auto exponent =
            std::dynamic_pointer_cast<const NumberNode>(power.exponent());
        BigInt integer_exponent;
        if (!try_get_integer_value(exponent, integer_exponent)) {
            return unknown_bounds();
        }
        const bool square = integer_exponent == BigInt(2);
        const bool inverse = exponent->is_negative_one();
        if (!square && !inverse) { return unknown_bounds(); }
        auto base = propagate(power.base());
        if (!base) { return base; }
        if (inverse) { return multiply_intervals(closed_bounds(Rational(1), Rational(1)), base.value(), true); }
        auto bounds = rational_bounds(base.value());
        if (!bounds) { return Result<Interval>::failure(bounds.error()); }
        auto a = arithmetic(bounds.value().first, bounds.value().first, '*');
        if (!a) { return Result<Interval>::failure(a.error()); }
        auto b = arithmetic(bounds.value().second, bounds.value().second, '*');
        if (!b) { return Result<Interval>::failure(b.error()); }
        if (bounds.value().first >= Rational(0))
            { return Result<Interval>::success(closed_bounds(a.value(), b.value())); }
        if (bounds.value().second <= Rational(0))
            { return Result<Interval>::success(closed_bounds(b.value(), a.value())); }
        auto comparison = arithmetic_budget(a.value(), b.value());
        if (!comparison) { return Result<Interval>::failure(comparison.error()); }
        return Result<Interval>::success(closed_bounds(Rational(0), std::max(a.value(), b.value())));
    }

    Result<void> arithmetic_budget(const Rational& a, const Rational& b) {
        auto step = context_.consume_steps(1, bounds_operation);
        if (!step) return step;
        const std::size_t a_bits = std::max(a.get_numerator().bit_length(), a.get_denominator().bit_length());
        const std::size_t b_bits = std::max(b.get_numerator().bit_length(), b.get_denominator().bit_length());
        if (a_bits >= std::numeric_limits<std::size_t>::max() - b_bits)
            return Result<void>::failure(CasErrc::ResourceLimit, "integer growth exceeds budget", bounds_operation);
        return context_.require_integer_bits(a_bits + b_bits + 1, bounds_operation);
    }

    Result<Rational> arithmetic(const Rational& a, const Rational& b, char operation) {
        auto growth = arithmetic_budget(a, b);
        if (!growth) return Result<Rational>::failure(growth.error());
        if (operation == '+') return Result<Rational>::success(a + b);
        if (operation == '/') return Result<Rational>::success(a / b);
        return Result<Rational>::success(a * b);
    }

    Result<Rational> rational_endpoint_bound(const Endpoint& endpoint, bool lower,
                                             bool separate_zero = false) {
        auto comparable = detail::comparable_endpoint(endpoint, context_, bounds_operation);
        if (!comparable) { return Result<Rational>::failure(comparable.error()); }
        auto& value = comparable.value();
        if (value.infinity) { return Result<Rational>::failure(CasErrc::Inconclusive,
            "unbounded endpoint has no finite rational enclosure", bounds_operation); }
        if (!value.algebraic) { return Result<Rational>::success(std::move(value.rational)); }
        auto& algebraic = *value.algebraic;
        if (separate_zero) {
            auto sign = detail::compare_comparable(
                value, detail::ComparableEndpoint{}, context_);
            if (!sign) { return Result<Rational>::failure(sign.error()); }
            if (sign.value() == 0) { return Result<Rational>::success(Rational(0)); }
            while (sign.value() > 0 ? algebraic.lower <= Rational(0) : algebraic.upper >= Rational(0)) {
                auto refined = detail::refine_exact_real_algebraic(algebraic, context_, bounds_operation);
                if (!refined) { return Result<Rational>::failure(refined.error()); }
            }
        }
        return Result<Rational>::success(lower ? algebraic.lower : algebraic.upper);
    }

    Result<RationalBounds> rational_bounds(const Interval& interval, bool separate_zero = false) {
        auto lower = rational_endpoint_bound(interval.lower, true, separate_zero);
        if (!lower) return Result<RationalBounds>::failure(lower.error());
        auto upper = rational_endpoint_bound(interval.upper, false, separate_zero);
        if (!upper) return Result<RationalBounds>::failure(upper.error());
        return Result<RationalBounds>::success({std::move(lower.value()), std::move(upper.value())});
    }

    Result<Interval> multiply_intervals(const Interval& left, const Interval& right, bool divide) {
        auto lhs = rational_bounds(left);
        if (!lhs) { return Result<Interval>::failure(lhs.error()); }
        auto rhs = rational_bounds(right, divide);
        if (!rhs) { return Result<Interval>::failure(rhs.error()); }
        if (divide && rhs.value().first <= Rational(0) && rhs.value().second >= Rational(0))
            { return unknown_bounds(); }
        std::array<Rational, 4> corners;
        std::size_t index = 0;
        for (const Rational* a : {&lhs.value().first, &lhs.value().second}) {
            for (const Rational* b : {&rhs.value().first, &rhs.value().second}) {
                auto corner = arithmetic(*a, *b, divide ? '/' : '*');
                if (!corner) { return Result<Interval>::failure(corner.error()); }
                corners[index++] = std::move(corner.value());
            }
        }
        for (const auto& corner : corners) {
            auto checked = arithmetic_budget(corner, corner);
            if (!checked) { return Result<Interval>::failure(checked.error()); }
        }
        auto extrema = std::minmax_element(corners.begin(), corners.end());
        return Result<Interval>::success(closed_bounds(*extrema.first, *extrema.second));
    }

    Result<Endpoint> arctangent_endpoint(const Endpoint& endpoint, bool lower) {
        if (lower ? endpoint.is_neg_infinity : endpoint.is_pos_infinity) {
            if (!pi_bounds_) {
                auto pi = detail::exact_pi_bounds(128, context_);
                if (!pi) { return Result<Endpoint>::failure(pi.error()); }
                pi_bounds_ = std::move(pi.value());
            }
            auto half = arithmetic(pi_bounds_->second, Rational(lower ? -2 : 2), '/');
            if (!half) { return Result<Endpoint>::failure(half.error()); }
            return Result<Endpoint>::success(Endpoint::open(SymbolicExpr::number(half.value())));
        }
        auto input = rational_endpoint_bound(endpoint, lower);
        if (!input) { return Result<Endpoint>::failure(input.error()); }
        auto bounds = detail::exact_atan_bounds(input.value(), 128, context_);
        if (!bounds) { return Result<Endpoint>::failure(bounds.error()); }
        auto value = SymbolicExpr::number(lower ? bounds.value().first : bounds.value().second);
        return Result<Endpoint>::success(endpoint.is_open ? Endpoint::open(value) : Endpoint::closed(value));
    }

    Result<Interval> function_bounds(const FunctionNode& function) {
        if (function.arguments().size() != 1) { return unknown_bounds(); }
        if (function.type() != FunctionNode::FuncType::Sin &&
            function.type() != FunctionNode::FuncType::Cos &&
            function.type() != FunctionNode::FuncType::ArcTan) { return unknown_bounds(); }
        const auto& argument = function.arguments()[0];
        auto defined = detail::query_definedness(argument, facts_, Domain::Real, context_);
        if (!defined) { return Result<Interval>::failure(defined.error()); }
        if (defined.value() != Tribool::True) { return unknown_bounds(); }
        auto argument_bounds = propagate(argument);
        if (!argument_bounds && argument_bounds.error().code != CasErrc::Inconclusive)
            { return argument_bounds; }
        if (!argument_bounds) {
            auto real = detail::query_real_value(argument, facts_, context_);
            if (!real) { return Result<Interval>::failure(real.error()); }
            if (real.value() != Tribool::True) { return unknown_bounds(); }
        }
        if (function.type() != FunctionNode::FuncType::ArcTan)
            { return Result<Interval>::success(closed_bounds(Rational(-1), Rational(1))); }
        return arctangent_bounds(argument_bounds);
    }

    Result<Interval> arctangent_bounds(Result<Interval>& argument_bounds) {
        Interval input = argument_bounds ? std::move(argument_bounds.value())
                                        : Interval{Endpoint::neg_inf(), Endpoint::pos_inf()};
        auto lower = arctangent_endpoint(input.lower, true);
        if (!lower) { return Result<Interval>::failure(lower.error()); }
        auto upper = arctangent_endpoint(input.upper, false);
        if (!upper) { return Result<Interval>::failure(upper.error()); }
        return Result<Interval>::success(Interval{std::move(lower.value()), std::move(upper.value())});
    }

    const AssumptionContext& assumptions_;
    detail::AssumptionFacts facts_;
    ComputationContext context_;
    std::optional<detail::ExactConstantInterval> pi_bounds_;
};
}

std::optional<Interval> InferenceEngine::propagate_bounds(const SymbolicExpr& expr) const {
    try {
        BoundsPropagation propagation(impl_->ctx, *this);
        auto result = propagation.propagate(detail::node(expr));
        if (result) return std::move(result.value());
    } catch (const CasError&) {
        return std::nullopt;
    } catch (const std::bad_alloc&) {
        return std::nullopt;
    } catch (const std::length_error&) {
        return std::nullopt;
    }
    return std::nullopt;
}

}
