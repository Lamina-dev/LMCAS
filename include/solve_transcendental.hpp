#pragma once

#include "symbolic.hpp"
#include <vector>
#include <memory>
#include <string>
#include <optional>

namespace lamina {

// 换元结果
struct SubstitutionResult {
    std::shared_ptr<SymbolicExpr> u_expr;    // u = h(x)
    std::shared_ptr<SymbolicExpr> poly_in_u; // 关于 u 的多项式
    std::string u_var;                        // 临时变量名
};

// 尝试用反函数法求解含超越函数的方程
LAMINA_API std::vector<std::shared_ptr<SymbolicExpr>> solve_transcendental(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var);

// 检测并执行换元: 找到 u = h(x) 使方程变为关于 u 的多项式
LAMINA_API std::optional<SubstitutionResult> detect_substitution(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::string& var);

} // namespace lamina
