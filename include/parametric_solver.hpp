/**
 * @file parametric_solver.hpp
 * @brief 含参方程组求解器 ParametricSolver，支持分段解。
 */
#pragma once
#include "symbolic.hpp"
#include <vector>
#include <map>
#include <string>
#include <memory>

namespace lamina {

/** @brief 含参方程组的分段解，按参数条件分类 */
struct PiecewiseSolution {
    /** @brief 单个分段：参数条件及对应的解集 */
    struct Case {
        std::shared_ptr<SymbolicExpr> condition; ///< 参数满足的条件表达式
        std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> solutions; ///< 该条件下的解集
    };
    std::vector<Case> cases; ///< 所有分段
};

using ParametricSolutionList =
    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>;
using ParametricSolutionsResult = Result<ParametricSolutionList>;

/** @brief 含参方程组求解器，支持线性与多项式方程组的参数化求解 */
class LAMINA_API ParametricSolver {
public:

    /**
     * @brief 求解含参方程组，返回所有解
     * @param equations 方程列表（每个表达式 = 0）
     * @param unknowns 未知数名称列表
     * @param parameters 参数名称列表
     * @return 解的列表，每个解为变量名到表达式的映射
     */
    static std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
    solve_system(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters);

    /**
     * @brief 求解含参方程组，返回按参数条件分段的解
     * @param equations 方程列表（每个表达式 = 0）
     * @param unknowns 未知数名称列表
     * @param parameters 参数名称列表
     * @return 分段解集
     */
    static PiecewiseSolution solve_system_piecewise(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters);

    static ParametricSolutionsResult solve_polynomial_parametric_checked(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters,
        ComputationContext& context);

    static ParametricSolutionsResult solve_polynomial_parametric_checked(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters);

private:

    /** @brief 求解关于未知数为线性的含参方程组 */
    static std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
    solve_linear_parametric(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters);

    /** @brief 求解关于未知数为多项式的含参方程组 */
    static ParametricSolutionList solve_polynomial_parametric_impl(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters,
        ComputationContext& context);

    /**
     * @brief 判断方程组是否关于未知数为线性
     * @param equations 方程列表
     * @param unknowns 未知数名称列表
     * @return 线性返回 true
     */
    static bool is_linear_in_unknowns(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns);
};

}
