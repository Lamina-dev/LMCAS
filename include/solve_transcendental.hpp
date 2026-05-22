#pragma once

#include "symbolic.hpp"
#include <vector>
#include <memory>
#include <string>
#include <optional>

namespace lamina {

struct SubstitutionResult {
    std::shared_ptr<SymbolicExpr> u_expr;
    std::shared_ptr<SymbolicExpr> poly_in_u;
    std::string u_var;
};

LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> solve_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var);

LAMINA_API std::optional<SubstitutionResult> detect_substitution(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var);

}
