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
namespace lamina {
class AssumptionContext;
}

namespace lamina {

/** @brief 多项式方程组求解器，提供线性系统求解、Gröbner 基计算及多项式系统求解功能。 */
class LAMINA_API Solver {
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
     * @brief 求解多项式方程组，返回所有解。
     * @param equations 方程列表（每个方程视为等于零的表达式）
     * @param variables 待求解的变量名列表
     * @return 解的列表，每个解为变量名到值的映射
     */
    static std::vector<std::map<std::string, SymbolicExpr>> solve_polynomial_system(
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

/**
 * @brief Solve an equation with assumption-based domain filtering.
 *
 * Calls the existing SymbolicExpr::solve() to compute all solutions, then
 * applies domain filtering based on the variable's declared domain in the
 * AssumptionContext. Solutions that violate the domain constraint are excluded.
 *
 * Domain filtering rules:
 * - Real: exclude solutions containing sqrt(-1) or non-real numeric values
 * - PositiveInt: exclude non-integer or <= 0 solutions
 * - NonNegative (Natural): exclude solutions evaluating to < 0
 *
 * @param equation The equation to solve (equal to zero, or RelationalNode)
 * @param variable The variable name to solve for
 * @param ctx Optional assumption context; if nullptr, all solutions returned unfiltered
 * @return Filtered solutions preserving original order
 */
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> solve_with_assumptions(
    const std::shared_ptr<SymbolicExpr>& equation,
    const std::string& variable,
    const AssumptionContext* ctx = nullptr);

}
