#pragma once
#include "symbolic.hpp"
#include <vector>
#include <map>
#include <string>
#include <memory>

namespace lamina {

// 分段解: 根据参数条件给出不同解
struct PiecewiseSolution {
    struct Case {
        std::shared_ptr<SymbolicExpr> condition;  // 参数条件 (如 a*e - b*d ≠ 0)
        std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> solutions;
    };
    std::vector<Case> cases;
};

// 参数方程组求解器
class LAMINA_API ParametricSolver {
public:
    // 三参数版本: 指定未知数和参数
    static std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
    solve_system(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters);

    // 带分段解的版本 (处理行列式依赖参数的情况)
    static PiecewiseSolution solve_system_piecewise(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters);

private:
    // 符号高斯消元 (线性系统)
    static std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
    solve_linear_parametric(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters);

    // 多项式系统 (Gröbner 基 + 回代)
    static std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>>
    solve_polynomial_parametric(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns,
        const std::vector<std::string>& parameters);

    // 判断系统是否为线性 (关于未知数)
    static bool is_linear_in_unknowns(
        const std::vector<std::shared_ptr<SymbolicExpr>>& equations,
        const std::vector<std::string>& unknowns);
};

} // namespace lamina
