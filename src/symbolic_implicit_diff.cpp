#include "symbolic_implicit_diff.hpp"

namespace LMCAS {

std::shared_ptr<SymbolicExpr> implicit_diff(
    std::shared_ptr<SymbolicExpr> expression,
    const std::string& independent_variable,
    const std::string& dependent_variable) {
    if (!expression || independent_variable.empty() || dependent_variable.empty()) {
        return nullptr;
    }

    auto numerator = SymbolicExpr::multiply(
        SymbolicExpr::number(-1),
        expression->differentiate(independent_variable));
    auto denominator = expression->differentiate(dependent_variable);
    return SymbolicExpr::divide(numerator, denominator);
}

} // namespace LMCAS
