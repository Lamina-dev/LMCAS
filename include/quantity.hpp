#pragma once

#include <memory>
#include <string>

#include "computation_context.hpp"
#include "symbolic.hpp"
#include "unit.hpp"

namespace lamina {

using QuantityResult = Result<std::shared_ptr<SymbolicExpr>>;

enum class UnitStripMode {
    BaseValue,
    DisplayValue
};

/** Attaches a declared unit to a scalar symbolic expression. */
LAMINA_API QuantityResult attach_unit(
    const std::shared_ptr<SymbolicExpr>& value,
    const std::string& unit,
    ComputationContext& context);

/** Converts a quantity to a declared unit of the same dimension. */
LAMINA_API QuantityResult convert_unit(
    const std::shared_ptr<SymbolicExpr>& quantity,
    const std::string& target_unit,
    ComputationContext& context);

/** Removes a quantity annotation using an explicit scale policy. */
LAMINA_API QuantityResult strip_unit(
    const std::shared_ptr<SymbolicExpr>& quantity,
    UnitStripMode mode,
    ComputationContext& context);

LAMINA_API Result<DimensionSignature> dimension_of(const SymbolicExpr& expression);
LAMINA_API QuantityResult quantity_add(
    const std::shared_ptr<SymbolicExpr>& lhs,
    const std::shared_ptr<SymbolicExpr>& rhs,
    ComputationContext& context);
LAMINA_API QuantityResult quantity_subtract(
    const std::shared_ptr<SymbolicExpr>& lhs,
    const std::shared_ptr<SymbolicExpr>& rhs,
    ComputationContext& context);
LAMINA_API QuantityResult quantity_multiply(
    const std::shared_ptr<SymbolicExpr>& lhs,
    const std::shared_ptr<SymbolicExpr>& rhs,
    ComputationContext& context);
LAMINA_API QuantityResult quantity_divide(
    const std::shared_ptr<SymbolicExpr>& lhs,
    const std::shared_ptr<SymbolicExpr>& rhs,
    ComputationContext& context);
LAMINA_API QuantityResult quantity_power(
    const std::shared_ptr<SymbolicExpr>& base,
    const std::shared_ptr<SymbolicExpr>& exponent,
    ComputationContext& context);

} // namespace lamina
