/**
 * @file numerical_integration.hpp
 * @brief 数值积分算法:辛普森法则和高斯求积(符号-数值混合桥接).
 */
#pragma once
#include "numeric_evaluation.hpp"
#include "symbolic.hpp"
#include <memory>
#include <string>

namespace lamina {

/**
 * @brief Checked composite Simpson integration.
 *
 * Endpoints and samples are evaluated through evaluate_numeric(); invalid
 * arguments, odd subinterval counts, unbound variables, domain errors,
 * cancellation, and budget exhaustion are reported explicitly. The returned
 * absolute_error is a step-doubling Richardson estimate, not a certified bound.
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
 * @brief Checked Gauss-Legendre integration for orders 1 through 20.
 *
 * Every supported order is evaluated by LMMC's Gauss-Legendre implementation.
 * The returned absolute_error is the difference from an adjacent-order rule,
 * not a certified bound.
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
 * @brief Checked adaptive Simpson integration.
 *
 * All endpoint and integrand evaluation goes through evaluate_numeric(); missing
 * variables, domain errors, cancellation, and budget exhaustion are represented
 * by CasError in the Result channel.
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

/** @brief Default checked numerical integration using composite Simpson. */
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
