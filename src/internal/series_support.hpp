#pragma once

#include "result.hpp"
#include "lmcas_export.hpp"
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace LMCAS {

class SymbolicNode;
class NumberNode;
class SymbolicExpr;


class ComputationContext;

namespace detail::series_support {

double series_number_value(const NumberNode& number);
bool series_is_number(const std::shared_ptr<SymbolicExpr>& expression);
double series_get_double(const std::shared_ptr<SymbolicExpr>& expression);
bool series_is_infinity(const std::shared_ptr<SymbolicExpr>& expression);

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

}
}
