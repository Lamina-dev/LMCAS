#pragma once

#include <memory>
#include <string>

#include "computation_context.hpp"
#include "symbolic.hpp"
#include "unit.hpp"

namespace LMCAS {

using QuantityResult = Result<std::shared_ptr<SymbolicExpr>>;

enum class UnitStripMode {
    BaseValue,
    DisplayValue
};

/** Attaches a declared unit to a scalar symbolic expression. */
LMCAS_API QuantityResult attach_unit(
    const std::shared_ptr<SymbolicExpr>& value,
    const std::string& unit,
    ComputationContext& context);

/** Attaches a resolved unit definition. */
LMCAS_API QuantityResult attach_unit(
    const std::shared_ptr<SymbolicExpr>& value,
    std::string display_unit,
    UnitDefinition definition,
    ComputationContext& context);

/** Converts a quantity to a declared unit of the same dimension. */
LMCAS_API QuantityResult convert_unit(
    const std::shared_ptr<SymbolicExpr>& quantity,
    const std::string& target_unit,
    ComputationContext& context);

/** Converts to a resolved unit definition. */
LMCAS_API QuantityResult convert_unit(
    const std::shared_ptr<SymbolicExpr>& quantity,
    std::string display_unit,
    UnitDefinition definition,
    ComputationContext& context);

/** Removes a quantity annotation using the requested scale. */
LMCAS_API QuantityResult strip_unit(
    const std::shared_ptr<SymbolicExpr>& quantity,
    UnitStripMode mode,
    ComputationContext& context);

LMCAS_API Result<DimensionSignature> dimension_of(const SymbolicExpr& expression);
LMCAS_API QuantityResult quantity_add(
    const std::shared_ptr<SymbolicExpr>& lhs,
    const std::shared_ptr<SymbolicExpr>& rhs,
    ComputationContext& context);
LMCAS_API QuantityResult quantity_subtract(
    const std::shared_ptr<SymbolicExpr>& lhs,
    const std::shared_ptr<SymbolicExpr>& rhs,
    ComputationContext& context);
LMCAS_API QuantityResult quantity_multiply(
    const std::shared_ptr<SymbolicExpr>& lhs,
    const std::shared_ptr<SymbolicExpr>& rhs,
    ComputationContext& context);
LMCAS_API QuantityResult quantity_divide(
    const std::shared_ptr<SymbolicExpr>& lhs,
    const std::shared_ptr<SymbolicExpr>& rhs,
    ComputationContext& context);
LMCAS_API QuantityResult quantity_power(
    const std::shared_ptr<SymbolicExpr>& base,
    const std::shared_ptr<SymbolicExpr>& exponent,
    ComputationContext& context);

} // namespace LMCAS
