#pragma once

#include "symbolic.hpp"
#include <vector>
#include <memory>
#include <string>

namespace lamina {

struct SolveOptions {
    bool allow_numeric = false;
    int max_newton_iterations = 100;
    lmmc_real_t tolerance = 1e-12;
    int max_roots = -1;
    bool return_rootof = true;
    lmmc_real_t initial_guess = 0.0;
    bool has_initial_guess = false;
};

struct NumericRoot {
    lmmc_real_t value;
    lmmc_real_t residual;
    int iterations;
};

enum class SolveStrategy {
    ClosedForm,
    Preprocessing,
    Transcendental,
    Numerical,
    RootOf
};

LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> solve_dispatch(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts);

}
