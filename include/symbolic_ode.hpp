/**
 * @file symbolic_ode.hpp
 * @brief 常微分方程求解：可分离变量、一阶线性、二阶常系数。
 */
#pragma once
#include "computation_context.hpp"
#include "symbolic.hpp"
#include <memory>
#include <string>

namespace lamina {

// Forward declaration to avoid circular include.
class AssumptionContext;

/**
 * @brief 求解可分离变量型 ODE：dy/dx = rhs(x, y)
 * @param rhs 方程右端表达式
 * @param x 自变量名
 * @param y 因变量名
 * @param ctx Optional assumption context. When provided and the dependent
 *            variable y is known Positive, the solver prefers positive solution
 *            branches. When nullptr, behavior is identical to the unparameterized call.
 * @return 通解的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> solve_separable_ode(
    std::shared_ptr<SymbolicExpr> rhs,
    const std::string& x,
    const std::string& y,
    const AssumptionContext* ctx = nullptr
);

/**
 * @brief 求解一阶线性 ODE：dy/dx + P(x)*y = Q(x)
 * @param Px 系数函数 P(x)
 * @param Qx 非齐次项 Q(x)
 * @param x 自变量名
 * @param y 因变量名
 * @param ctx Optional assumption context. When provided and the dependent
 *            variable y is known Positive, the solver prefers positive solution
 *            branches. When a coefficient is known NonZero, the solver skips
 *            zero-coefficient degenerate case checks. When nullptr, behavior is
 *            identical to the unparameterized call.
 * @return 通解的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> solve_linear1_ode(
    std::shared_ptr<SymbolicExpr> Px,
    std::shared_ptr<SymbolicExpr> Qx,
    const std::string& x,
    const std::string& y,
    const AssumptionContext* ctx = nullptr
);

/**
 * @brief 求解二阶常系数线性 ODE：a*y'' + b*y' + c*y = f(x)
 * @param a 二阶导数系数
 * @param b 一阶导数系数
 * @param c 零阶项系数
 * @param fx 非齐次项 f(x)
 * @param x 自变量名
 * @param y 因变量名
 * @param ctx Optional assumption context. When provided and a coefficient is
 *            known NonZero, the solver skips the zero-coefficient degenerate
 *            case check for that coefficient. When the dependent variable y is
 *            known Positive, the solver prefers positive solution branches.
 *            When nullptr, behavior is identical to the unparameterized call.
 * @return 通解的符号表达式
 */
LAMINA_API std::shared_ptr<SymbolicExpr> solve_linear2_ode(
    double a, double b, double c,
    std::shared_ptr<SymbolicExpr> fx,
    const std::string& x,
    const std::string& y,
    const AssumptionContext* ctx = nullptr
);

/**
 * @brief Checked migration API for second-order constant-coefficient ODEs.
 *
 * Supported domain currently covers homogeneous equations and first-order
 * degeneracy when a == 0 and b != 0. Non-homogeneous second-order equations
 * are outside the current support domain and return `CasErrc::Inconclusive`
 * instead of throwing or returning the homogeneous solution as a false success.
 */
LAMINA_API Result<std::shared_ptr<SymbolicExpr>> solve_linear2_ode_checked(
    double a, double b, double c,
    std::shared_ptr<SymbolicExpr> fx,
    const std::string& x,
    const std::string& y,
    ComputationContext& context,
    const AssumptionContext* ctx = nullptr
);

LAMINA_API Result<std::shared_ptr<SymbolicExpr>> solve_linear2_ode_checked(
    double a, double b, double c,
    std::shared_ptr<SymbolicExpr> fx,
    const std::string& x,
    const std::string& y,
    const AssumptionContext* ctx = nullptr
);

}
