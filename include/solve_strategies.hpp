#pragma once

#include "symbolic.hpp"
#include <vector>
#include <memory>
#include <string>

namespace lamina {

// 求解选项
struct SolveOptions {
    bool allow_numeric = false;       // 是否允许数值解
    int max_newton_iterations = 100;  // Newton 最大迭代
    lmmc_real_t tolerance = 1e-12;    // 收敛容差
    int max_roots = -1;               // -1 = 全部根
    bool return_rootof = true;        // deg>=5 时返回 RootOf
    lmmc_real_t initial_guess = 0.0;  // 非多项式方程的初始猜测 x0
    bool has_initial_guess = false;   // 是否提供了初始猜测
};

// 数值根
struct NumericRoot {
    lmmc_real_t value;
    lmmc_real_t residual;  // |f(value)|
    int iterations;
};

// 求解策略枚举
enum class SolveStrategy {
    ClosedForm,       // 闭合公式 (deg 1-4)
    Preprocessing,    // 有理根 + 因式分解
    Transcendental,   // 反函数法
    Numerical,        // Newton-Raphson
    RootOf            // 符号占位
};

// Strategy dispatcher: applies strategies in priority order and returns the first non-empty result.
// Never throws; returns an empty vector when all strategies fail.
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> solve_dispatch(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var,
    const SolveOptions& opts);

} // namespace lamina
