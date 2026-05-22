#pragma once

#include "symbolic.hpp"
#include "polynomial.hpp"
#include "rational.hpp"
#include "solve_strategies.hpp"
#include <vector>
#include <memory>
#include <string>
#include <optional>
#include <utility>

namespace lamina {

LAMINA_API std::vector<NumericRoot> solve_numeric(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts = {});

LAMINA_API std::vector<std::pair<Rational, Rational>> isolate_real_roots(
    const Polynomial<Rational>& poly);

LAMINA_API std::optional<NumericRoot> newton_raphson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    const SolveOptions& opts);

LAMINA_API std::optional<NumericRoot> newton_raphson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& df,
    const std::string& var,
    lmmc_real_t x0,
    lmmc_real_t bracket_lo,
    lmmc_real_t bracket_hi,
    const SolveOptions& opts);

LAMINA_API std::optional<NumericRoot> bisection(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    lmmc_real_t lo,
    lmmc_real_t hi,
    const SolveOptions& opts);

}
