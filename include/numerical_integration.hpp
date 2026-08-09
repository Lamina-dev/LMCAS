/**
 * @file numerical_integration.hpp
 * @brief 数值积分算法：辛普森法则和高斯求积（符号-数值混合桥接）。
 */
#pragma once
#include "numeric_evaluation.hpp"
#include "symbolic.hpp"
#include <memory>
#include <string>

namespace lamina {

/**
 * @brief Compatibility wrapper for checked Simpson integration.
 * @param f 被积函数
 * @param var 积分变量
 * @param a 积分下限
 * @param b 积分上限
 * @param n 区间等分数（必须为偶数）
 * @return 积分近似值的符号数值；checked 失败时返回 nullptr
 */
LAMINA_API std::shared_ptr<SymbolicExpr> quadrature_simpson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n = 100
);

/**
 * @brief Checked composite Simpson integration.
 *
 * Endpoints and samples are evaluated through evaluate_numeric(); invalid
 * arguments, odd subinterval counts, unbound variables, domain errors,
 * cancellation, and budget exhaustion are reported explicitly.
 */
LAMINA_API Result<ApproxReal> quadrature_simpson_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    int n = 100
);

LAMINA_API Result<ApproxReal> quadrature_simpson_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n = 100
);

/**
 * @brief Compatibility wrapper for checked Gauss-Legendre integration.
 * @param f 被积函数
 * @param var 积分变量
 * @param a 积分下限
 * @param b 积分上限
 * @param n 积分点数（通常为 5, 7, 或 10）
 * @return 积分近似值的符号数值；checked 失败时返回 nullptr
 */
LAMINA_API std::shared_ptr<SymbolicExpr> quadrature_gaussian(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n = 5
);

/**
 * @brief Checked Gauss-Legendre integration for supported point counts.
 *
 * Supports n = 1, 2, and 3 directly. Larger n falls back to checked composite
 * Simpson with 2n subintervals, matching the legacy wrapper's support domain.
 */
LAMINA_API Result<ApproxReal> quadrature_gaussian_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    int n = 5
);

LAMINA_API Result<ApproxReal> quadrature_gaussian_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n = 5
);

/**
 * @brief 自适应辛普森积分：递归细分直至误差估计低于容差。
 * @param f 被积函数
 * @param var 积分变量
 * @param a 下限（数值）
 * @param b 上限（数值）
 * @param tol 误差容差
 * @return 积分近似值（NumberNode）
 */
LAMINA_API std::shared_ptr<SymbolicExpr> adaptive_simpson(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    double tol = 1e-10
);

/**
 * @brief Checked adaptive Simpson integration.
 *
 * All endpoint and integrand evaluation goes through evaluate_numeric(); missing
 * variables, domain errors, cancellation, and budget exhaustion are reported as
 * CasError instead of being collapsed to nullptr or zero.
 */
LAMINA_API Result<ApproxReal> adaptive_simpson_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    double tol = 1e-10,
    int max_depth = 50
);

LAMINA_API Result<ApproxReal> adaptive_simpson_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    double tol = 1e-10,
    int max_depth = 50
);

/**
 * @brief Compatibility wrapper for the checked default numeric integrator.
 * @param f 被积函数
 * @param var 积分变量
 * @param a 下限
 * @param b 上限
 * @param n 子区间数（偶数）
 * @return 积分近似值的符号数值；checked 失败时返回 nullptr
 */
LAMINA_API std::shared_ptr<SymbolicExpr> numerical_integrate(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n = 100
);

/** @brief Checked default numerical integration entry point. */
LAMINA_API Result<ApproxReal> numerical_integrate_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    ComputationContext& context,
    int n = 100
);

LAMINA_API Result<ApproxReal> numerical_integrate_numeric(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& a,
    const std::shared_ptr<SymbolicExpr>& b,
    int n = 100
);

} // namespace lamina
