#pragma once

#include "symbolic.hpp"
#include "polynomial.hpp"
#include "poly_utils.hpp"
#include "lmmc/config.h"
#include <vector>
#include <memory>
#include <string>
#include <optional>

namespace lamina {

LAMINA_API std::optional<lmmc_real_t> rootof_evaluate(
    const std::shared_ptr<SymbolicExpr>& rootof_expr);

LAMINA_API std::shared_ptr<SymbolicExpr> rootof_simplify(
    const std::shared_ptr<SymbolicExpr>& rootof_expr);

LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> make_rootof_solutions(
    const Polynomial<SymbolicPolyCoeff>& poly,
    const std::string& var);

}
