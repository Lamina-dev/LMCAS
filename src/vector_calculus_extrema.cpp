#include "internal/vector_calculus_support.hpp"
#include "integration.hpp"
#include "numeric_evaluation.hpp"
#include "solver.hpp"
#include "solve_strategies.hpp"
#include "symbolic_ast.hpp"

#include <cmath>
#include <exception>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace lamina {

using namespace vector_calculus_detail;

static bool vector_calculus_evaluate_hessian_at_point(
    const std::shared_ptr<SymbolicExpr>& H,
    const std::vector<std::string>&,
    const std::map<std::string, std::shared_ptr<SymbolicExpr>>& pt,
    size_t n,
    std::vector<double>& numeric_H)
{
    auto mat_node = std::dynamic_pointer_cast<const MatrixNode>(lamina::detail::node(H));
    if (!mat_node || mat_node->rows() != n || mat_node->cols() != n) {
        return false;
    }

    numeric_H.resize(n * n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            auto elem_node = mat_node->get(i, j);
            if (!elem_node) {
                numeric_H[i * n + j] = 0.0;
                continue;
            }
            auto elem = lamina::detail::make_expression_ptr(elem_node);
            /// 代入临界点坐标
            for (const auto& [var_name, val] : pt) {
                elem = elem->substitute(var_name, val);
                if (!elem) return false;
            }
            elem = elem->simplify();
            if (!elem) return false;

            if (!vector_calculus_checked_finite_numeric(
                    elem, numeric_H[i * n + j])) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @internal
 * @brief 使用特征值分析对临界点进行分类。
 *
 * 对 n×n 实对称矩阵，直接通过数值分析特征值符号判断正定性。
 * 对于小矩阵（n ≤ 3），使用解析公式；对于大矩阵，使用特征多项式求解。
 *
 * @param[in] numeric_H 数值化的海森矩阵（行优先，n×n）
 * @param[in] n         矩阵维度
 * @return 分类字符串: "minimum", "maximum", "saddle", "degenerate"
 */
static CriticalPointClassification vector_calculus_classify_critical_point(
    const std::vector<double>& numeric_H, size_t n)
{
    const double tol = 1e-10;

    /// 1×1 情况：直接判断
    if (n == 1) {
        double val = numeric_H[0];
        if (std::abs(val) < tol) return CriticalPointClassification::Degenerate;
        if (val > 0) return CriticalPointClassification::LocalMinimum;
        return CriticalPointClassification::LocalMaximum;
    }

    /// 2×2 情况：使用行列式和迹直接判断
    if (n == 2) {
        double a = numeric_H[0], b = numeric_H[1];
        double c = numeric_H[2], d = numeric_H[3];
        double det = a * d - b * c;
        double trace = a + d;

        if (std::abs(det) < tol) return CriticalPointClassification::Degenerate;
        if (det > 0 && trace > 0) return CriticalPointClassification::LocalMinimum;
        if (det > 0 && trace < 0) return CriticalPointClassification::LocalMaximum;
        return CriticalPointClassification::Saddle;
    }

    /// 一般情况：构建符号矩阵并求特征值
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            /// 使用 Rational 保持 Hessian 系数的精确表示。
            double val = numeric_H[i * n + j];
            int int_val = static_cast<int>(std::round(val));
            if (std::abs(val - int_val) < tol) {
                grid[i][j] = SymbolicExpr::number(int_val);
            } else {
                grid[i][j] = SymbolicExpr::number(val);
            }
        }
    }
    auto H_mat = SymbolicExpr::matrix(grid);

    /// 计算特征多项式并求解特征值
    auto cp = SymbolicExpr::charpoly(H_mat, "lambda");
    if (!cp) return CriticalPointClassification::Inconclusive;

    auto eigenvals = lamina::solve_finite_checked(cp, "lambda").value();
    if (eigenvals.empty()) return CriticalPointClassification::Inconclusive;

    bool all_positive = true;
    bool all_negative = true;
    bool has_zero = false;

    for (const auto& ev : eigenvals) {
        if (!ev) {
            return CriticalPointClassification::Inconclusive;
        }
        auto ev_simplified = ev->simplify();
        double val = 0.0;
        if (!vector_calculus_checked_finite_numeric(ev_simplified, val)) {
            return CriticalPointClassification::Inconclusive;
        }

        if (std::abs(val) < tol) {
            has_zero = true;
            all_positive = false;
            all_negative = false;
        } else if (val > 0) {
            all_negative = false;
        } else {
            all_positive = false;
        }
    }

    if (has_zero) return CriticalPointClassification::Degenerate;
    if (all_positive) return CriticalPointClassification::LocalMinimum;
    if (all_negative) return CriticalPointClassification::LocalMaximum;
    return CriticalPointClassification::Saddle;
}

ExtremaResult find_extrema_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    const std::string operation = "find_extrema";
    auto valid = vector_calculus_validate_expr_vars(f, vars, context, operation);
    if (!valid) return ExtremaResult::failure(valid.error());
    auto distinct = vector_calculus_validate_distinct_vars(vars, vars.size(), context, operation);
    if (!distinct) return ExtremaResult::failure(distinct.error());
    auto step = context.consume_steps(vars.size() * vars.size() + vars.size() * 8 + 8,
                                      operation);
    if (!step) return ExtremaResult::failure(step.error());

    try {
        std::vector<std::shared_ptr<SymbolicExpr>> gradient;
        gradient.reserve(vars.size());
        for (const auto& var : vars) {
            auto partial = vector_calculus_differentiate_strict(f, var, operation);
            if (!partial) return ExtremaResult::failure(partial.error());
            gradient.push_back(std::move(partial.value()));
        }

        auto extrema = find_extrema(f, vars);
        if (extrema.empty()) {
            return ExtremaResult::failure(
                CasErrc::Inconclusive,
                "extrema solver produced no verifiable critical points in the supported domain",
                operation);
        }

        for (const auto& cp : extrema) {
            if (!vector_calculus_point_has_vars(cp.point, vars)) {
                return ExtremaResult::failure(
                    CasErrc::Inconclusive,
                    "extrema candidate omits one or more variables",
                    operation);
            }
            for (const auto& partial : gradient) {
                if (!vector_calculus_expr_zero_after_substitution(partial, cp.point)) {
                    return ExtremaResult::failure(
                        CasErrc::Inconclusive,
                        "extrema candidate does not verify against the gradient",
                        operation);
                }
            }
            if (cp.classification == CriticalPointClassification::Inconclusive) {
                return ExtremaResult::failure(
                    CasErrc::Inconclusive,
                    "extrema candidate has no verified classification",
                    operation);
            }
        }
        return ExtremaResult::success(std::move(extrema));
    } catch (const std::bad_alloc&) {
        return ExtremaResult::failure(CasErrc::ResourceLimit,
                                      "extrema allocation failed",
                                      operation);
    } catch (const std::exception& e) {
        return ExtremaResult::failure(CasErrc::InternalInvariant,
                                      e.what(), operation);
    }
}

ExtremaResult find_extrema_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return find_extrema_checked(f, vars, context);
}

std::vector<CriticalPoint> find_extrema(
    const std::shared_ptr<SymbolicExpr>& f, const std::vector<std::string>& vars)
{
    if (!f || vars.empty()) {
        return {};
    }

    size_t n = vars.size();

    /// 计算梯度 ∇f
    std::vector<std::shared_ptr<SymbolicExpr>> grad_eqs;
    grad_eqs.reserve(n);
    for (const auto& var : vars) {
        auto partial = f->differentiate(var);
        if (partial) {
            partial = partial->simplify();
        }
        grad_eqs.push_back(partial);
    }

    /// 求解 ∇f = 0 系统（使用多项式系统求解器）
    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> solutions;

    std::vector<SymbolicExpr> poly_eqs;
    poly_eqs.reserve(n);
    for (const auto& eq : grad_eqs) {
        if (eq) poly_eqs.push_back(*eq);
    }
    if (poly_eqs.size() == n) {
        auto poly_solutions = Solver::solve_polynomial_system(poly_eqs, vars);
        for (const auto& sol : poly_solutions) {
            std::map<std::string, std::shared_ptr<SymbolicExpr>> pt;
            for (const auto& [name, val] : sol) {
                pt[name] = lamina::detail::make_expression_ptr(val);
            }
            solutions.push_back(pt);
        }
    }

    /// 如果多项式求解器失败，尝试线性求解器
    if (solutions.empty()) {
        solutions = SymbolicExpr::solve_system(grad_eqs, vars);
    }

    if (solutions.empty()) {
        return {};
    }

    /// 计算海森矩阵
    auto H = hessian(f, vars);

    /// 对每个临界点进行分类
    std::vector<CriticalPoint> result;
    result.reserve(solutions.size());

    for (const auto& sol : solutions) {
        CriticalPoint cp;
        cp.point = sol;

        /// 在临界点处求值海森矩阵
        std::vector<double> numeric_H;
        if (vector_calculus_evaluate_hessian_at_point(H, vars, sol, n, numeric_H)) {
            cp.classification = vector_calculus_classify_critical_point(numeric_H, n);
        } else {
            cp.classification = CriticalPointClassification::Inconclusive;
        }

        result.push_back(std::move(cp));
    }

    return result;
}


static LagrangeResult lagrange_multipliers_strict(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars,
    const std::string& operation)
{
    const size_t n = vars.size();
    const size_t m = constraints.size();

    VectorField grad_f;
    grad_f.reserve(n);
    for (const auto& var : vars) {
        auto partial = vector_calculus_differentiate_strict(f, var, operation);
        if (!partial) return LagrangeResult::failure(partial.error());
        grad_f.push_back(std::move(partial.value()));
    }

    std::vector<VectorField> grad_constraints;
    grad_constraints.reserve(m);
    for (const auto& constraint : constraints) {
        VectorField grad_g;
        grad_g.reserve(n);
        for (const auto& var : vars) {
            auto partial = vector_calculus_differentiate_strict(
                constraint, var, operation);
            if (!partial) return LagrangeResult::failure(partial.error());
            grad_g.push_back(std::move(partial.value()));
        }
        grad_constraints.push_back(std::move(grad_g));
    }

    std::vector<std::string> lambda_names;
    lambda_names.reserve(m);
    for (size_t k = 0; k < m; ++k) {
        lambda_names.push_back("lambda_" + std::to_string(k + 1));
    }

    std::vector<std::shared_ptr<SymbolicExpr>> equations;
    equations.reserve(n + m);
    for (size_t i = 0; i < n; ++i) {
        auto eq = grad_f[i];
        for (size_t k = 0; k < m; ++k) {
            auto lambda_var = SymbolicExpr::variable(lambda_names[k]);
            auto term = SymbolicExpr::multiply(lambda_var, grad_constraints[k][i]);
            eq = SymbolicExpr::add(
                eq, SymbolicExpr::multiply(SymbolicExpr::number(-1), term));
        }
        auto eq_checked = vector_calculus_simplify_strict(
            eq, operation, "Lagrange stationarity equation is outside the supported domain");
        if (!eq_checked) return LagrangeResult::failure(eq_checked.error());
        equations.push_back(std::move(eq_checked.value()));
    }
    for (const auto& constraint : constraints) {
        equations.push_back(constraint);
    }

    std::vector<std::string> all_vars;
    all_vars.reserve(n + m);
    all_vars.insert(all_vars.end(), vars.begin(), vars.end());
    all_vars.insert(all_vars.end(), lambda_names.begin(), lambda_names.end());

    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> full_solutions;
    std::vector<SymbolicExpr> poly_eqs;
    poly_eqs.reserve(equations.size());
    for (const auto& eq : equations) {
        if (!eq || !lamina::detail::node(eq)) {
            return LagrangeResult::failure(
                CasErrc::Inconclusive,
                "Lagrange equation construction failed in the supported domain",
                operation);
        }
        poly_eqs.push_back(*eq);
    }

    auto poly_solutions = Solver::solve_polynomial_system(poly_eqs, all_vars);
    for (const auto& sol : poly_solutions) {
        std::map<std::string, std::shared_ptr<SymbolicExpr>> pt;
        for (const auto& [name, val] : sol) {
            pt[name] = lamina::detail::make_expression_ptr(val);
        }
        full_solutions.push_back(std::move(pt));
    }

    if (full_solutions.empty()) {
        full_solutions = SymbolicExpr::solve_system(equations, all_vars);
    }
    if (full_solutions.empty()) {
        return LagrangeResult::failure(
            CasErrc::Inconclusive,
            "Lagrange solver produced no verifiable candidates in the supported domain",
            operation);
    }

    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> verified;
    std::set<std::string> var_set(vars.begin(), vars.end());
    for (const auto& full_solution : full_solutions) {
        if (!vector_calculus_point_has_vars(full_solution, all_vars)) {
            return LagrangeResult::failure(
                CasErrc::Inconclusive,
                "Lagrange candidate omits one or more variables or multipliers",
                operation);
        }
        for (const auto& eq : equations) {
            if (!vector_calculus_expr_zero_after_substitution(eq, full_solution)) {
                return LagrangeResult::failure(
                    CasErrc::Inconclusive,
                    "Lagrange candidate does not verify against the full stationarity system",
                    operation);
            }
        }

        std::map<std::string, std::shared_ptr<SymbolicExpr>> filtered;
        for (const auto& [name, val] : full_solution) {
            if (var_set.count(name)) {
                filtered[name] = val;
            }
        }
        if (!vector_calculus_point_has_vars(filtered, vars)) {
            return LagrangeResult::failure(
                CasErrc::Inconclusive,
                "Lagrange candidate omits one or more variables",
                operation);
        }
        verified.push_back(std::move(filtered));
    }

    if (verified.empty()) {
        return LagrangeResult::failure(
            CasErrc::Inconclusive,
            "Lagrange solver produced no verifiable candidates in the supported domain",
            operation);
    }
    return LagrangeResult::success(std::move(verified));
}

LagrangeResult lagrange_multipliers_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    const std::string operation = "lagrange_multipliers";
    auto valid = vector_calculus_validate_expr_vars(f, vars, context, operation);
    if (!valid) return LagrangeResult::failure(valid.error());
    auto distinct = vector_calculus_validate_distinct_vars(vars, vars.size(), context, operation);
    if (!distinct) return LagrangeResult::failure(distinct.error());
    if (constraints.empty()) {
        return LagrangeResult::failure(CasErrc::InvalidArgument,
                                       "constraint list cannot be empty",
                                       operation);
    }
    for (const auto& constraint : constraints) {
        if (!constraint || !lamina::detail::node(constraint)) {
            return LagrangeResult::failure(CasErrc::InvalidArgument,
                                           "constraint expressions cannot be null",
                                           operation);
        }
    }
    auto step = context.consume_steps(vars.size() * constraints.size() * 6 +
                                      vars.size() * 6 + constraints.size() * 4 + 8,
                                      operation);
    if (!step) return LagrangeResult::failure(step.error());

    try {
        return lagrange_multipliers_strict(f, constraints, vars, operation);
    } catch (const std::bad_alloc&) {
        return LagrangeResult::failure(CasErrc::ResourceLimit,
                                       "Lagrange multiplier allocation failed",
                                       operation);
    } catch (const std::exception& e) {
        return LagrangeResult::failure(CasErrc::InternalInvariant,
                                       e.what(), operation);
    }
}

LagrangeResult lagrange_multipliers_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars)
{
    ComputationContext context;
    return lagrange_multipliers_checked(f, constraints, vars, context);
}

std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> lagrange_multipliers(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::shared_ptr<SymbolicExpr>>& constraints,
    const std::vector<std::string>& vars)
{
    if (!f || constraints.empty() || vars.empty()) {
        return {};
    }

    size_t n = vars.size();
    size_t m = constraints.size();

    /// 计算目标函数的梯度
    VectorField grad_f;
    grad_f.reserve(n);
    for (const auto& var : vars) {
        auto partial = f->differentiate(var);
        if (partial) partial = partial->simplify();
        grad_f.push_back(partial);
    }

    /// 计算每个约束的梯度
    std::vector<VectorField> grad_constraints;
    grad_constraints.reserve(m);
    for (const auto& g : constraints) {
        VectorField grad_g;
        grad_g.reserve(n);
        for (const auto& var : vars) {
            auto partial = g->differentiate(var);
            if (partial) partial = partial->simplify();
            grad_g.push_back(partial);
        }
        grad_constraints.push_back(std::move(grad_g));
    }

    /// 构造乘数变量名
    std::vector<std::string> lambda_names;
    lambda_names.reserve(m);
    for (size_t k = 0; k < m; ++k) {
        lambda_names.push_back("lambda_" + std::to_string(k + 1));
    }

    /// 构造方程组: ∂f/∂xᵢ - ∑ λₖ · ∂gₖ/∂xᵢ = 0 (对每个 i)
    /// 加上约束方程: gₖ = 0 (对每个 k)
    std::vector<std::shared_ptr<SymbolicExpr>> equations;
    equations.reserve(n + m);

    for (size_t i = 0; i < n; ++i) {
        /// ∂f/∂xᵢ - λ₁·∂g₁/∂xᵢ - λ₂·∂g₂/∂xᵢ - ...
        auto eq = grad_f[i];
        for (size_t k = 0; k < m; ++k) {
            auto lambda_var = SymbolicExpr::variable(lambda_names[k]);
            auto term = SymbolicExpr::multiply(lambda_var, grad_constraints[k][i]);
            eq = SymbolicExpr::add(eq, SymbolicExpr::multiply(SymbolicExpr::number(-1), term));
        }
        if (eq) eq = eq->simplify();
        equations.push_back(eq);
    }

    /// 添加约束方程
    for (const auto& g : constraints) {
        equations.push_back(g);
    }

    /// 构造所有未知数列表: 原始变量 + 乘数
    std::vector<std::string> all_vars;
    all_vars.reserve(n + m);
    for (const auto& v : vars) {
        all_vars.push_back(v);
    }
    for (const auto& lam : lambda_names) {
        all_vars.push_back(lam);
    }

    /// 使用多项式系统求解器求解增广系统
    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> solutions;

    std::vector<SymbolicExpr> poly_eqs;
    poly_eqs.reserve(equations.size());
    for (const auto& eq : equations) {
        if (eq) poly_eqs.push_back(*eq);
    }

    if (poly_eqs.size() == equations.size()) {
        auto poly_solutions = Solver::solve_polynomial_system(poly_eqs, all_vars);
        for (const auto& sol : poly_solutions) {
            std::map<std::string, std::shared_ptr<SymbolicExpr>> pt;
            for (const auto& [name, val] : sol) {
                pt[name] = lamina::detail::make_expression_ptr(val);
            }
            solutions.push_back(pt);
        }
    }

    /// 如果多项式求解器失败，尝试线性求解器
    if (solutions.empty()) {
        solutions = SymbolicExpr::solve_system(equations, all_vars);
    }

    /// 从解中提取仅原始变量（去除乘数）
    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> result;
    std::set<std::string> var_set(vars.begin(), vars.end());

    for (const auto& sol : solutions) {
        std::map<std::string, std::shared_ptr<SymbolicExpr>> filtered;
        for (const auto& [name, val] : sol) {
            if (var_set.count(name)) {
                filtered[name] = val;
            }
        }
        if (!filtered.empty()) {
            result.push_back(std::move(filtered));
        }
    }

    return result;
}


} // namespace lamina
