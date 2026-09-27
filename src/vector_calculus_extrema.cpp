#include "internal/vector_calculus_support.hpp"
#include "vector_calculus_extrema.hpp"
#include "vector_calculus_matrices.hpp"
#include "internal/visitors/differentiation_visitor.hpp"
#include "integration.hpp"
#include "numeric_evaluation.hpp"
#include "solver.hpp"
#include "solve_strategies.hpp"
#include "symbolic_matrix.hpp"
#include "internal/symbolic_ast.hpp"

#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace LMCAS {

using namespace vector_calculus_detail;

static bool vector_calculus_evaluate_hessian_at_point(
    const std::shared_ptr<SymbolicExpr>& H,
    const std::vector<std::string>&,
    const std::map<std::string, std::shared_ptr<SymbolicExpr>>& pt,
    size_t n,
    std::vector<double>& numeric_H)
{
    auto mat_node = std::dynamic_pointer_cast<const MatrixNode>(LMCAS::detail::node(H));
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
            auto elem = LMCAS::detail::make_expression_ptr(elem_node);
            /// 代入临界点坐标
            for (const auto& [var_name, val] : pt) {
                elem = elem->substitute(var_name, val);
                if (!elem) {
                    return false;
                }
            }
            elem = elem->simplify();
            if (!elem) {
                return false;
            }

            if (!vector_calculus_checked_finite_numeric(
                    elem, numeric_H[i * n + j])) {
                return false;
            }
        }
    }
    return true;
}

static CriticalPointClassification classify_hessian_pair(
    const std::vector<double>& numeric_H, double tol)
{
        double a = numeric_H[0], b = numeric_H[1];
        double c = numeric_H[2], d = numeric_H[3];
        double det = a * d - b * c;
        double trace = a + d;

        if (std::abs(det) < tol) {
            return CriticalPointClassification::Degenerate;
        }
        if (det > 0 && trace > 0) {
            return CriticalPointClassification::LocalMinimum;
        }
        if (det > 0 && trace < 0) {
            return CriticalPointClassification::LocalMaximum;
        }
        return CriticalPointClassification::Saddle;
}

static std::shared_ptr<SymbolicExpr> symbolic_hessian(
    const std::vector<double>& numeric_H, size_t n, double tol)
{
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> grid(n,
        std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double val = numeric_H[i * n + j];
            const double rounded = std::round(val);
            const bool representable_integer =
                std::isfinite(val) &&
                rounded >=
                    static_cast<double>(std::numeric_limits<int>::min()) &&
                rounded <=
                    static_cast<double>(std::numeric_limits<int>::max()) &&
                std::abs(val - rounded) < tol;
            if (representable_integer) {
                grid[i][j] =
                    SymbolicExpr::number(static_cast<int>(rounded));
            } else {
                grid[i][j] = SymbolicExpr::number(val);
            }
        }
    }
    return SymbolicExpr::matrix(grid);
}

static CriticalPointClassification classify_hessian_eigenvalues(
    const std::vector<std::shared_ptr<SymbolicExpr>>& eigenvals, double tol)
{
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

    if (has_zero) {
        return CriticalPointClassification::Degenerate;
    }
    if (all_positive) {
        return CriticalPointClassification::LocalMinimum;
    }
    if (all_negative) {
        return CriticalPointClassification::LocalMaximum;
    }
    return CriticalPointClassification::Saddle;
}

// The 1e-10 threshold matches the existing real-eigenvalue classification
// tolerance used by the one- and two-dimensional paths.
static CriticalPointClassification vector_calculus_classify_critical_point(
    const std::vector<double>& numeric_H, size_t n)
{
    const double tol = 1e-10;
    if (n == 1) {
        const double value = numeric_H[0];
        if (std::abs(value) < tol) {
            return CriticalPointClassification::Degenerate;
        }
        if (value > 0) {
            return CriticalPointClassification::LocalMinimum;
        }
        return CriticalPointClassification::LocalMaximum;
    }
    if (n == 2) {
        return classify_hessian_pair(numeric_H, tol);
    }
    auto matrix = symbolic_hessian(numeric_H, n, tol);
    auto eigenvalues = matrix_eigenvalues_checked(matrix);
    if (!eigenvalues) {
        return CriticalPointClassification::Inconclusive;
    }
    return classify_hessian_eigenvalues(eigenvalues.value(), tol);
}


ExtremaResult find_extrema_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::vector<std::string>& vars,
    ComputationContext& context)
{
    const std::string operation = "find_extrema";
    auto valid = vector_calculus_validate_expr_vars(f, vars, context, operation);
    if (!valid) {
        return ExtremaResult::failure(valid.error());
    }
    auto distinct = vector_calculus_validate_distinct_vars(vars, vars.size(), context, operation);
    if (!distinct) {
        return ExtremaResult::failure(distinct.error());
    }
    auto step = context.consume_steps(vars.size() * vars.size() + vars.size() * 8 + 8,
                                      operation);
    if (!step) {
        return ExtremaResult::failure(step.error());
    }

    try {
        std::vector<std::shared_ptr<SymbolicExpr>> gradient;
        gradient.reserve(vars.size());
        for (const auto& var : vars) {
            auto partial = vector_calculus_differentiate_strict(f, var, operation);
            if (!partial) {
                return ExtremaResult::failure(partial.error());
            }
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
    } catch (const detail::UnsupportedDifferentiation& error) {
        return ExtremaResult::failure(
            CasErrc::UnsupportedExpression, error.what(), operation);
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

    /// 计算梯度 gradf
    std::vector<std::shared_ptr<SymbolicExpr>> grad_eqs;
    grad_eqs.reserve(n);
    for (const auto& var : vars) {
        auto partial = f->differentiate(var);
        if (partial) {
            partial = partial->simplify();
        }
        grad_eqs.push_back(partial);
    }

    /// 求解 gradf = 0 系统(使用多项式系统求解器)
    std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> solutions;

    std::vector<SymbolicExpr> poly_eqs;
    poly_eqs.reserve(n);
    for (const auto& eq : grad_eqs) {
        if (eq) poly_eqs.push_back(*eq);
    }
    if (poly_eqs.size() == n) {
        auto checked_solutions =
            Solver::solve_polynomial_system_checked(poly_eqs, vars);
        const auto poly_solutions = checked_solutions
            ? std::move(checked_solutions.value())
            : std::vector<std::map<std::string, SymbolicExpr>>{};
        for (const auto& sol : poly_solutions) {
            std::map<std::string, std::shared_ptr<SymbolicExpr>> pt;
            for (const auto& [name, val] : sol) {
                pt[name] = LMCAS::detail::make_expression_ptr(val);
            }
            solutions.push_back(pt);
        }
    }

    /// 如果多项式求解器失败,尝试线性求解器
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

}
