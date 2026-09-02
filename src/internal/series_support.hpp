#pragma once

#include "result.hpp"
#include <memory>
#include <optional>
#include <string>
#include <vector>

class SymbolicNode;
class SymbolicExpr;

namespace lamina {
class ComputationContext;

namespace detail::series_support {

Result<void> validate_power_series_coefficients(
    const std::vector<std::shared_ptr<SymbolicExpr>>& coefficients,
    const std::string& operation,
    const std::string& name);

Result<void> validate_series_variable(const std::string& variable,
                                      ComputationContext& context,
                                      const std::string& operation);

std::optional<int> supported_laurent_integer_power(
    const std::shared_ptr<const SymbolicNode>& node,
    const std::string& variable);

} // namespace detail::series_support
} // namespace lamina
