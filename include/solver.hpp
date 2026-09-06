/**
 * @file solver.hpp
 * @brief 多项式方程组求解器 Solver，包含 Gröbner 基、线性系统求解。
 */
#pragma once
#include "symbolic.hpp"
#include <vector>
#include <string>
#include <map>

// Forward declaration for assumption-aware solving
namespace LMCAS {
class AssumptionContext;
}

namespace LMCAS {

/** @brief 多项式方程组求解器，提供线性系统求解、Gröbner 基计算及多项式系统求解功能。 */
using PolynomialSystemResult =
    Result<std::vector<std::map<std::string, SymbolicExpr>>>;
class LMCAS_API Solver {
public:

    /**
     * @brief 求解线性方程组。
     * @param equations 方程列表（每个方程视为等于零的表达式）
     * @param variables 待求解的变量名列表
     * @return 变量名到解表达式的映射
     */
    static std::map<std::string, SymbolicExpr> solve_linear_system(
        const std::vector<SymbolicExpr>& equations,
        const std::vector<std::string>& variables);

    /**
     * @brief 计算多项式集合的 Gröbner 基。
     * @param polynomials 输入多项式列表
     * @param variables 变量名列表（决定单项式序）
     * @return Gröbner 基多项式列表
     */
    static std::vector<SymbolicExpr> groebner_basis(
        const std::vector<SymbolicExpr>& polynomials,
        const std::vector<std::string>& variables);

    /**
     * @brief 求解多项式方程组，显式报告失败。
     */
    static PolynomialSystemResult solve_polynomial_system_checked(
        const std::vector<SymbolicExpr>& equations,
        const std::vector<std::string>& variables,
        ComputationContext& context);

    static PolynomialSystemResult solve_polynomial_system_checked(
        const std::vector<SymbolicExpr>& equations,
        const std::vector<std::string>& variables);

    /**
     * @brief 计算约化 Gröbner 基。
     * @param polynomials 输入多项式列表
     * @param variables 变量名列表
     * @return 约化 Gröbner 基多项式列表
     */
    static std::vector<SymbolicExpr> reduced_groebner_basis(
        const std::vector<SymbolicExpr>& polynomials,
        const std::vector<std::string>& variables);

    /**
     * @brief 判断多项式是否属于给定 Gröbner 基生成的理想。
     * @param polynomial 待检测的多项式
     * @param basis Gröbner 基
     * @param variables 变量名列表
     * @return 若多项式属于该理想则返回 true
     */
    static bool ideal_membership(
        const SymbolicExpr& polynomial,
        const std::vector<SymbolicExpr>& basis,
        const std::vector<std::string>& variables);

    /**
     * @brief 计算消元理想，消去前 elim_count 个变量。
     * @param basis Gröbner 基
     * @param variables 变量名列表（前 elim_count 个将被消去）
     * @param elim_count 要消去的变量数量
     * @return 消元理想的生成元列表
     */
    static std::vector<SymbolicExpr> elimination_ideal(
        const std::vector<SymbolicExpr>& basis,
        const std::vector<std::string>& variables,
        int elim_count);
};

using AssumptionSolveResult =
    Result<std::vector<std::shared_ptr<SymbolicExpr>>>;

/**
 * @brief Solve an equation with assumption-based domain filtering.
 *
 * Computes all solutions and applies the optional assumption context's domain
 * and sign constraints. The computation context is shared with nested solving.
 */
LMCAS_API AssumptionSolveResult solve_with_assumptions_checked(
    const std::shared_ptr<SymbolicExpr>& equation,
    const std::string& variable,
    const AssumptionContext* assumptions,
    ComputationContext& context);

LMCAS_API AssumptionSolveResult solve_with_assumptions_checked(
    const std::shared_ptr<SymbolicExpr>& equation,
    const std::string& variable,
    const AssumptionContext* assumptions = nullptr);

}
