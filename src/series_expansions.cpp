#include "series_engine.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "integration.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/series_support.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace LMCAS {

static int detect_parity(const std::shared_ptr<SymbolicExpr>& f, const std::string& var) {
    if (!f) return 0;
    auto neg_x = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::variable(var));
    auto f_neg = f->substitute(var, neg_x);
    if (!f_neg) return 0;
    f_neg = f_neg->simplify();
    auto f_s = f->simplify();
    auto diff_even = SymbolicExpr::add(f_neg, SymbolicExpr::multiply(SymbolicExpr::number(-1), f_s));
    if (diff_even) { diff_even = diff_even->simplify(); if (diff_even->is_zero()) return 1; }
    auto diff_odd = SymbolicExpr::add(f_neg, f_s);
    if (diff_odd) { diff_odd = diff_odd->simplify(); if (diff_odd->is_zero()) return -1; }
    return 0;
}



using detail::series_support::supported_laurent_integer_power;
using detail::series_support::validate_series_variable;

Result<void> validate_series_expr_point(
    const std::shared_ptr<SymbolicExpr>& expr,
    const std::shared_ptr<SymbolicExpr>& center,
    const std::string& var,
    ComputationContext& context,
    const std::string& operation)
{
    auto var_check = validate_series_variable(var, context, operation);
    if (!var_check) return var_check;
    if (!expr || !LMCAS::detail::node(expr) || !center || !LMCAS::detail::node(center)) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "series expression and center cannot be null",
                                     operation);
    }
    return Result<void>::success();
}


Result<void> validate_laurent_orders(int order_neg,
                                     int order_pos,
                                     const std::string& operation)
{
    if (order_neg < 0 || order_pos < 0 || order_neg > 64 || order_pos > 64) {
        return Result<void>::failure(CasErrc::InvalidArgument,
                                     "Laurent truncation orders must be between 0 and 64",
                                     operation);
    }
    return Result<void>::success();
}
static LaurentSeriesResult laurent_series_full_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&,
    const std::shared_ptr<SymbolicExpr>&,
    int, int, ComputationContext&);
ExpressionResult fourier_series_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& period, int n_terms,
    ComputationContext& context) {
    constexpr const char* operation = "fourier_series";
    if (!f || !period || var.empty() || n_terms < 0) {
        return ExpressionResult::failure(
            CasErrc::InvalidArgument,
            "Fourier expansion requires an expression, variable, positive period, and non-negative term count",
            operation);
    }
    const std::size_t expansion_terms =
        1 + 2 * static_cast<std::size_t>(n_terms);
    if (expansion_terms > context.limits().max_expansion_terms) {
        return ExpressionResult::failure(
            CasErrc::ResourceLimit,
            "Fourier term count exceeds the expansion budget", operation);
    }

    try {
        auto access = context.consume_steps(1, operation);
        if (!access) return ExpressionResult::failure(access.error());

        auto T = period;
        auto pi = LMCAS::detail::make_expression_ptr(
            LMCAS::detail::make_node<VariableNode>("pi"));
        auto two = SymbolicExpr::number(2);
        auto L = SymbolicExpr::divide(T, two);
        auto half_lo = SymbolicExpr::multiply(SymbolicExpr::number(-1), L);
        auto half_hi = L;
        auto w = SymbolicExpr::divide(
            SymbolicExpr::multiply(two, pi), T);

        int parity = detect_parity(f, var);
        Integrator integrator;
        auto x = SymbolicExpr::variable(var);

        auto a0_integral = integrator.integrate_def_checked(
            *f, var, *half_lo, *half_hi, context);
        if (!a0_integral) {
            return ExpressionResult::failure(a0_integral.error());
        }
        auto a0 = SymbolicExpr::divide(
            LMCAS::detail::make_expression_ptr(a0_integral.value()), L)->simplify();
        auto result = SymbolicExpr::divide(a0, two);

        for (int k = 1; k <= n_terms; ++k) {
            auto kw = SymbolicExpr::multiply(SymbolicExpr::number(k), w);
            auto arg = SymbolicExpr::multiply(kw, x);

            if (parity != -1) {
                auto integrand = SymbolicExpr::multiply(
                    f, SymbolicExpr::cos(arg));
                auto integrated = integrator.integrate_def_checked(
                    *integrand, var, *half_lo, *half_hi, context);
                if (!integrated) {
                    return ExpressionResult::failure(integrated.error());
                }
                auto coefficient = SymbolicExpr::divide(
                    LMCAS::detail::make_expression_ptr(integrated.value()), L)->simplify();
                if (!(LMCAS::detail::node(coefficient) &&
                      LMCAS::detail::node(coefficient)->is_zero())) {
                    result = SymbolicExpr::add(
                        result,
                        SymbolicExpr::multiply(
                            coefficient, SymbolicExpr::cos(arg)));
                }
            }
            if (parity != 1) {
                auto integrand = SymbolicExpr::multiply(
                    f, SymbolicExpr::sin(arg));
                auto integrated = integrator.integrate_def_checked(
                    *integrand, var, *half_lo, *half_hi, context);
                if (!integrated) {
                    return ExpressionResult::failure(integrated.error());
                }
                auto coefficient = SymbolicExpr::divide(
                    LMCAS::detail::make_expression_ptr(integrated.value()), L)->simplify();
                if (!(LMCAS::detail::node(coefficient) &&
                      LMCAS::detail::node(coefficient)->is_zero())) {
                    result = SymbolicExpr::add(
                        result,
                        SymbolicExpr::multiply(
                            coefficient, SymbolicExpr::sin(arg)));
                }
            }
        }
        return ExpressionResult::success(result->simplify());
    } catch (const std::bad_alloc&) {
        return ExpressionResult::failure(
            CasErrc::ResourceLimit,
            "Fourier expansion allocation failed", operation);
    } catch (const std::exception& error) {
        return ExpressionResult::failure(
            CasErrc::InternalInvariant, error.what(), operation);
    }
}

ExpressionResult fourier_series_checked(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& period, int n_terms) {
    ComputationContext context;
    return fourier_series_checked(f, var, period, n_terms, context);
}


ExpressionResult laurent_series_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos,
    ComputationContext& context)
{
    auto full = laurent_series_full_checked(
        f, var, center, order_neg, order_pos, context);
    if (!full) return ExpressionResult::failure(full.error());
    auto simplified = f->simplify();
    auto supported_power = simplified
        ? supported_laurent_integer_power(
              LMCAS::detail::node(simplified), var)
        : std::nullopt;
    if (center->is_zero() && supported_power && *supported_power < 0) {
        return ExpressionResult::success(std::move(simplified));
    }
    return ExpressionResult::success(std::move(full.value().series));
}

ExpressionResult laurent_series_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos)
{
    ComputationContext context;
    return laurent_series_checked(f, var, center, order_neg, order_pos, context);
}

LaurentSeriesResult laurent_series_full_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos,
    ComputationContext& context)
{
    const std::string operation = "laurent_series_full";
    auto input = validate_series_expr_point(f, center, var, context, operation);
    if (!input) return LaurentSeriesResult::failure(input.error());
    auto order_check = validate_laurent_orders(order_neg, order_pos, operation);
    if (!order_check) return LaurentSeriesResult::failure(order_check.error());
    auto budget = context.consume_steps(static_cast<std::size_t>(order_neg + order_pos + 1) * 16 + 16,
                                        operation);
    if (!budget) return LaurentSeriesResult::failure(budget.error());
    auto simplified_input = f->simplify();
    if (!simplified_input || !LMCAS::detail::node(simplified_input)) {
        return LaurentSeriesResult::failure(
            CasErrc::Inconclusive,
            "checked Laurent input could not be simplified in the supported domain",
            operation);
    }
    auto supported_power = supported_laurent_integer_power(
        LMCAS::detail::node(simplified_input), var);
    try {
        auto built = laurent_series_full_impl(
            f, var, center, order_neg, order_pos, context);
        if (!built) return built;
        auto result = std::move(built.value());
        if (!result.series || !LMCAS::detail::node(result.series)) {
            return LaurentSeriesResult::failure(
                CasErrc::Inconclusive,
                "Laurent series could not be constructed in the supported domain",
                operation);
        }
        if (!result.residue || !LMCAS::detail::node(result.residue)) {
            return LaurentSeriesResult::failure(
                CasErrc::InternalInvariant,
                "Laurent series produced a null residue expression",
                operation);
        }
        if (supported_power && center->is_zero()) {
            int expected_pole_order =
                *supported_power < 0 ? -*supported_power : 0;
            if (result.pole_order != expected_pole_order) {
                return LaurentSeriesResult::failure(
                    CasErrc::Inconclusive,
                    "Laurent pole order could not be verified",
                    operation);
            }
        }
        return LaurentSeriesResult::success(std::move(result));
    } catch (const std::bad_alloc&) {
        return LaurentSeriesResult::failure(CasErrc::ResourceLimit,
                                            "allocation failed while calculating Laurent series",
                                            operation);
    } catch (const std::exception& ex) {
        return LaurentSeriesResult::failure(CasErrc::InternalInvariant,
                                            ex.what(),
                                            operation);
    }
}

LaurentSeriesResult laurent_series_full_checked(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    int order_neg,
    int order_pos)
{
    ComputationContext context;
    return laurent_series_full_checked(f, var, center, order_neg, order_pos, context);
}

static LaurentSeriesResult laurent_series_full_impl(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center, int order_neg, int order_pos,
    ComputationContext& context) {
    LaurentResult res{nullptr, SingularityType::Removable, 0, SymbolicExpr::number(0)};
    if (!f || !center) return res;

    auto x = SymbolicExpr::variable(var);
    /// (x - center)
    auto shift = SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), center));

    /// 寻找极点阶数 m:使 (x-c)^m * f 在 c 处解析(有限).
    int pole_order = 0;
    std::shared_ptr<SymbolicExpr> regular = f;
    {
        auto test = f;
        const int MAX_M = (order_neg > 0 ? order_neg : 5);
        /// 先判断 f 在 center 是否已解析
        auto f_at = f->substitute(var, center);
        if (f_at) f_at = f_at->simplify();
        bool finite = f_at && f_at->is_number();
        if (!finite) {
            for (int m = 1; m <= MAX_M; ++m) {
                auto powm = SymbolicExpr::power(shift, SymbolicExpr::number(m));
                auto g = SymbolicExpr::multiply(powm, f)->simplify();
                auto g_at = g->substitute(var, center);
                if (g_at) g_at = g_at->simplify();
                if (g_at && g_at->is_number()) {
                    pole_order = m;
                    regular = g;       // (x-c)^m f(x),在 c 处解析
                    break;
                }
            }
        }
    }

    /// 对 regular(解析部分)做 Taylor 展开
    int taylor_order = order_pos + pole_order + 1;
    if (taylor_order < 1) taylor_order = 1;
    auto taylor = regular->series(var, center, taylor_order);
    if (!taylor) return res;

    /// Laurent = taylor / (x-c)^pole_order
    std::shared_ptr<SymbolicExpr> laurent = taylor;
    if (pole_order > 0) {
        auto powm = SymbolicExpr::power(shift, SymbolicExpr::number(pole_order));
        laurent = SymbolicExpr::divide(taylor, powm)->simplify();
    }

    res.series = laurent->simplify();
    res.pole_order = pole_order;
    if (pole_order == 0) res.singularity = SingularityType::Removable;
    else res.singularity = SingularityType::Pole;

    /// 留数 = Taylor 展开中 (x-c)^(pole_order-1) 的系数
    if (pole_order >= 1) {
        /// a_{m-1} = (1/(m-1)!) lim_{x->c} d^{m-1}/dx^{m-1} [regular]
        auto deriv = regular;
        for (int i = 0; i < pole_order - 1; ++i) deriv = deriv->differentiate(var);
        /// 用极限求值,稳健处理 0/0 形式(如 z/(z(z+1)) 在 z=0).
        auto limited = limit_expression_checked(
            deriv, var, center, LimitDirection::Both, context);
        if (!limited) {
            return LaurentSeriesResult::failure(limited.error());
        }
        auto val = std::move(limited.value());
        if (!val) val = deriv->substitute(var, center);
        long long fact = 1;
        for (int i = 2; i <= pole_order - 1; ++i) fact *= i;
        res.residue = SymbolicExpr::divide(val, SymbolicExpr::number((int)fact))->simplify();
    }

    return res;
}

std::shared_ptr<SymbolicExpr> asymptotic_expand(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var, int order) {
    if (!f || order < 0) return nullptr;
    /// 通过 x = 1/t 替换,在 t=0 处做 Taylor 展开,再回代 t = 1/x,
    /// 得到按 x 递减幂次的渐近展开.
    auto x = SymbolicExpr::variable(var);
    std::string tname = var + "__asym_t";
    auto t = SymbolicExpr::variable(tname);
    auto inv_t = SymbolicExpr::divide(SymbolicExpr::number(1), t);
    auto g = f->substitute(var, inv_t);
    if (!g) return nullptr;
    g = g->simplify();
    auto taylor_t = g->series(tname, SymbolicExpr::number(0), order + 1);
    if (!taylor_t) return nullptr;
    /// 回代 t = 1/x
    auto inv_x = SymbolicExpr::divide(SymbolicExpr::number(1), x);
    auto back = taylor_t->substitute(tname, inv_x);
    if (!back) return nullptr;
    return back->simplify();
}
} // namespace LMCAS
