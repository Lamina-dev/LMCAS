#include "quantity.hpp"

#include <exception>
#include <utility>

#include "symbolic_ast.hpp"

namespace lamina {
namespace {

constexpr const char* kOperation = "quantity";

struct QuantityView {
    std::shared_ptr<const SymbolicNode> value;
    DimensionSignature dimension;
    Rational scale = Rational(1);
    std::string display;
};

QuantityResult failure(CasErrc code, std::string message) {
    return QuantityResult::failure(code, std::move(message), kOperation);
}

QuantityView view(const std::shared_ptr<SymbolicExpr>& expression) {
    if (auto quantity = std::dynamic_pointer_cast<const QuantityNode>(detail::node(expression))) {
        return {quantity->value(), quantity->dimension(), quantity->scale_to_base(),
                quantity->display_unit()};
    }
    return {detail::node(expression), DimensionSignature{}, Rational(1), {}};
}

std::shared_ptr<const SymbolicNode> rational_node(const Rational& value) {
    return SymbolicFactory::create_number(value);
}

std::shared_ptr<const SymbolicNode> scaled_value(const QuantityView& quantity) {
    if (quantity.scale == Rational(1)) return quantity.value;
    return SymbolicFactory::create_multiply({quantity.value, rational_node(quantity.scale)});
}

QuantityResult make_quantity(std::shared_ptr<const SymbolicNode> value,
                             DimensionSignature dimension,
                             Rational scale,
                             std::string display) {
    try {
        if (dimension.is_dimensionless()) {
            return QuantityResult::success(detail::make_expression_ptr(
                scale == Rational(1)
                    ? std::move(value)
                    : SymbolicFactory::create_multiply({std::move(value), rational_node(scale)})));
        }
        return QuantityResult::success(detail::make_expression_ptr(
            detail::make_node<QuantityNode>(std::move(value), std::move(dimension),
                                            std::move(scale), std::move(display))));
    } catch (const std::bad_alloc&) {
        return failure(CasErrc::ResourceLimit, "quantity allocation failed");
    } catch (const std::exception& error) {
        return failure(CasErrc::InvalidArgument, error.what());
    }
}

Result<Rational> exact_rational(const std::shared_ptr<const SymbolicNode>& node) {
    auto number = std::dynamic_pointer_cast<const NumberNode>(node);
    if (!number) return Result<Rational>::failure(
        CasErrc::UnitInvalid, "quantity exponent must be an exact rational", kOperation);
    if (std::holds_alternative<BigInt>(number->value())) {
        return Result<Rational>::success(Rational(std::get<BigInt>(number->value())));
    }
    if (std::holds_alternative<Rational>(number->value())) {
        return Result<Rational>::success(std::get<Rational>(number->value()));
    }
    return Result<Rational>::failure(
        CasErrc::UnitInvalid, "approximate quantity exponents are not supported", kOperation);
}

Result<void> checked_binary(const std::shared_ptr<SymbolicExpr>& lhs,
                            const std::shared_ptr<SymbolicExpr>& rhs,
                            ComputationContext& context,
                            const char* name) {
    auto step = context.consume_steps(1, name);
    if (!step) return step;
    if (!lhs || !rhs || !detail::node(lhs) || !detail::node(rhs)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "quantity operands cannot be null", kOperation);
    }
    return Result<void>::success();
}

} // namespace

QuantityResult attach_unit(const std::shared_ptr<SymbolicExpr>& value,
                           const std::string& unit,
                           ComputationContext& context) {
    auto step = context.consume_steps(1, kOperation);
    if (!step) return QuantityResult::failure(step.error());
    if (!value || !detail::node(value)) {
        return failure(CasErrc::InvalidArgument, "quantity value cannot be null");
    }
    if (std::dynamic_pointer_cast<const QuantityNode>(detail::node(value))) {
        return failure(CasErrc::UnitInvalid, "a quantity cannot receive a second unit");
    }
    auto definition = context.units().resolve(unit);
    if (!definition) return QuantityResult::failure(definition.error());
    return make_quantity(detail::node(value), definition.value().dimension,
                         definition.value().scale_to_base, unit);
}

QuantityResult attach_unit(const std::shared_ptr<SymbolicExpr>& value,
                           std::string display_unit,
                           UnitDefinition definition,
                           ComputationContext& context) {
    auto step = context.consume_steps(1, kOperation);
    if (!step) return QuantityResult::failure(step.error());
    if (!value || !detail::node(value)) {
        return failure(CasErrc::InvalidArgument, "quantity value cannot be null");
    }
    if (std::dynamic_pointer_cast<const QuantityNode>(detail::node(value))) {
        return failure(CasErrc::UnitInvalid, "a quantity cannot receive a second unit");
    }
    if (display_unit.empty() || definition.scale_to_base <= Rational(0)) {
        return failure(CasErrc::UnitInvalid, "resolved unit definition is invalid");
    }
    return make_quantity(detail::node(value), std::move(definition.dimension),
                         std::move(definition.scale_to_base),
                         std::move(display_unit));
}

QuantityResult convert_unit(const std::shared_ptr<SymbolicExpr>& expression,
                            const std::string& target_unit,
                            ComputationContext& context) {
    auto step = context.consume_steps(1, kOperation);
    if (!step) return QuantityResult::failure(step.error());
    auto quantity = expression
        ? std::dynamic_pointer_cast<const QuantityNode>(detail::node(expression))
        : nullptr;
    if (!quantity) {
        return failure(CasErrc::UnitStripTypeMismatch,
                       "unit conversion requires a quantity expression");
    }
    auto target = context.units().resolve(target_unit);
    if (!target) return QuantityResult::failure(target.error());
    if (quantity->dimension() != target.value().dimension) {
        return failure(CasErrc::DimensionMismatch,
                       "unit conversion requires equal dimensions");
    }
    const Rational factor = quantity->scale_to_base() / target.value().scale_to_base;
    auto value = SymbolicFactory::create_multiply(
        {quantity->value(), rational_node(factor)});
    return make_quantity(std::move(value), target.value().dimension,
                         target.value().scale_to_base, target_unit);
}

QuantityResult convert_unit(const std::shared_ptr<SymbolicExpr>& expression,
                            std::string display_unit,
                            UnitDefinition definition,
                            ComputationContext& context) {
    auto step = context.consume_steps(1, kOperation);
    if (!step) return QuantityResult::failure(step.error());
    auto quantity = expression
        ? std::dynamic_pointer_cast<const QuantityNode>(detail::node(expression))
        : nullptr;
    if (!quantity) {
        return failure(CasErrc::UnitStripTypeMismatch,
                       "unit conversion requires a quantity expression");
    }
    if (display_unit.empty() || definition.scale_to_base <= Rational(0)) {
        return failure(CasErrc::UnitInvalid, "resolved unit definition is invalid");
    }
    if (quantity->dimension() != definition.dimension) {
        return failure(CasErrc::DimensionMismatch,
                       "unit conversion requires equal dimensions");
    }
    const Rational factor = quantity->scale_to_base() / definition.scale_to_base;
    auto value = SymbolicFactory::create_multiply(
        {quantity->value(), rational_node(factor)});
    return make_quantity(std::move(value), std::move(definition.dimension),
                         std::move(definition.scale_to_base),
                         std::move(display_unit));
}

QuantityResult strip_unit(const std::shared_ptr<SymbolicExpr>& expression,
                          UnitStripMode mode,
                          ComputationContext& context) {
    auto step = context.consume_steps(1, kOperation);
    if (!step) return QuantityResult::failure(step.error());
    if (!expression || !detail::node(expression)) {
        return failure(CasErrc::InvalidArgument,
                       "unit stripping requires an expression");
    }
    auto quantity = std::dynamic_pointer_cast<const QuantityNode>(
        detail::node(expression));
    if (!quantity) return QuantityResult::success(expression);
    auto value = mode == UnitStripMode::BaseValue
        ? SymbolicFactory::create_multiply(
              {quantity->value(), rational_node(quantity->scale_to_base())})
        : quantity->value();
    return QuantityResult::success(detail::make_expression_ptr(std::move(value)));
}

Result<DimensionSignature> dimension_of(const SymbolicExpr& expression) {
    if (!detail::node(expression)) {
        return Result<DimensionSignature>::failure(
            CasErrc::InvalidArgument, "expression cannot be empty", kOperation);
    }
    if (auto quantity = std::dynamic_pointer_cast<const QuantityNode>(detail::node(expression))) {
        return Result<DimensionSignature>::success(quantity->dimension());
    }
    return Result<DimensionSignature>::success(DimensionSignature{});
}

QuantityResult quantity_add(const std::shared_ptr<SymbolicExpr>& lhs,
                            const std::shared_ptr<SymbolicExpr>& rhs,
                            ComputationContext& context) {
    auto checked = checked_binary(lhs, rhs, context, "quantity.add");
    if (!checked) return QuantityResult::failure(checked.error());
    const auto left = view(lhs);
    const auto right = view(rhs);
    if (left.dimension != right.dimension) {
        return failure(CasErrc::DimensionMismatch,
                       "addition requires operands with equal dimensions");
    }
    const Rational factor = right.scale / left.scale;
    auto right_value = SymbolicFactory::create_multiply(
        {right.value, rational_node(factor)});
    return make_quantity(SymbolicFactory::create_add({left.value, right_value}),
                         left.dimension, left.scale, left.display);
}

QuantityResult quantity_subtract(const std::shared_ptr<SymbolicExpr>& lhs,
                                 const std::shared_ptr<SymbolicExpr>& rhs,
                                 ComputationContext& context) {
    auto checked = checked_binary(lhs, rhs, context, "quantity.subtract");
    if (!checked) return QuantityResult::failure(checked.error());
    const auto left = view(lhs);
    const auto right = view(rhs);
    if (left.dimension != right.dimension) {
        return failure(CasErrc::DimensionMismatch,
                       "subtraction requires operands with equal dimensions");
    }
    const Rational factor = Rational(-1) * right.scale / left.scale;
    auto right_value = SymbolicFactory::create_multiply(
        {right.value, rational_node(factor)});
    return make_quantity(SymbolicFactory::create_add({left.value, right_value}),
                         left.dimension, left.scale, left.display);
}

QuantityResult quantity_multiply(const std::shared_ptr<SymbolicExpr>& lhs,
                                 const std::shared_ptr<SymbolicExpr>& rhs,
                                 ComputationContext& context) {
    auto checked = checked_binary(lhs, rhs, context, "quantity.multiply");
    if (!checked) return QuantityResult::failure(checked.error());
    const auto left = view(lhs);
    const auto right = view(rhs);
    const auto dimension = left.dimension.multiplied_by(right.dimension);
    const std::string display = left.display.empty() ? right.display
        : right.display.empty() ? left.display : left.display + "*" + right.display;
    return make_quantity(SymbolicFactory::create_multiply({left.value, right.value}),
                         dimension, left.scale * right.scale, display);
}

QuantityResult quantity_divide(const std::shared_ptr<SymbolicExpr>& lhs,
                               const std::shared_ptr<SymbolicExpr>& rhs,
                               ComputationContext& context) {
    auto checked = checked_binary(lhs, rhs, context, "quantity.divide");
    if (!checked) return QuantityResult::failure(checked.error());
    const auto left = view(lhs);
    const auto right = view(rhs);
    auto inverse = SymbolicFactory::create_power(right.value,
                                                  SymbolicFactory::create_number(BigInt(-1)));
    const std::string display = right.display.empty() ? left.display
        : left.display.empty() ? "1/" + right.display : left.display + "/" + right.display;
    return make_quantity(SymbolicFactory::create_multiply({left.value, inverse}),
                         left.dimension.divided_by(right.dimension),
                         left.scale / right.scale, display);
}

QuantityResult quantity_power(const std::shared_ptr<SymbolicExpr>& base,
                              const std::shared_ptr<SymbolicExpr>& exponent,
                              ComputationContext& context) {
    auto checked = checked_binary(base, exponent, context, "quantity.power");
    if (!checked) return QuantityResult::failure(checked.error());
    const auto exponent_view = view(exponent);
    if (!exponent_view.dimension.is_dimensionless()) {
        return failure(CasErrc::DimensionMismatch,
                       "a quantity exponent must be dimensionless");
    }
    auto power = exact_rational(exponent_view.value);
    if (!power) return QuantityResult::failure(power.error());
    const auto base_view = view(base);
    auto base_value = scaled_value(base_view);
    return make_quantity(SymbolicFactory::create_power(base_value, exponent_view.value),
                         base_view.dimension.raised_to(power.value()),
                         Rational(1), base_view.dimension.raised_to(power.value()).to_string());
}

} // namespace lamina
