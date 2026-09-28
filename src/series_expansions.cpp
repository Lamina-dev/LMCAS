#include "limit_result.hpp"
#include "series_engine.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "integration.hpp"
#include "internal/expression_analysis.hpp"
#include "internal/series_support.hpp"
#include "internal/assumption_facts.hpp"
#include "polynomial_conversion.hpp"
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
    if (!f) {
        return 0;
    }
    auto neg_x = SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::variable(var));
    auto f_neg = f->substitute(var, neg_x);
    if (!f_neg) {
        return 0;
    }
    f_neg = f_neg->simplify();
    auto f_s = f->simplify();
    auto diff_even = SymbolicExpr::add(f_neg, SymbolicExpr::multiply(SymbolicExpr::number(-1), f_s));
    if (diff_even) { diff_even = diff_even->simplify(); if (diff_even->is_zero()) {
        return 1;
    } }
    auto diff_odd = SymbolicExpr::add(f_neg, f_s);
    if (diff_odd) { diff_odd = diff_odd->simplify(); if (diff_odd->is_zero()) {
        return -1;
    } }
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
    if (!var_check) {
        return var_check;
    }
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

static std::optional<int> essential_exponent(
    const std::shared_ptr<SymbolicExpr>& f,
    const std::shared_ptr<SymbolicExpr>& shift) {
    auto fn = std::dynamic_pointer_cast<const FunctionNode>(detail::node(f));
    if (!fn || fn->type() != FunctionNode::FuncType::Exp ||
        fn->arguments().size() != 1) {
        return std::nullopt;
    }
    auto argument = detail::make_expression_ptr(fn->arguments()[0])->simplify();
    auto power = std::dynamic_pointer_cast<const PowerNode>(detail::node(argument));
    if (!power) {
        return std::nullopt;
    }
    auto exponent = exact_small_integer_node(
        power->exponent(), -std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
    if (!exponent) {
        return std::nullopt;
    }
    auto base = detail::make_expression_ptr(power->base())->simplify();
    auto delta = shift->simplify();
    if (base->compare(delta) == 0 && *exponent < 0) {
        return -*exponent;
    }
    auto nested = std::dynamic_pointer_cast<const PowerNode>(detail::node(base));
    if (nested && *exponent == -1) {
        auto m = exact_small_integer_node(
            nested->exponent(), 1, std::numeric_limits<int>::max());
        auto nested_base = detail::make_expression_ptr(nested->base())->simplify();
        if (m && nested_base->compare(delta) == 0) {
            return m;
        }
    }
    return std::nullopt;
}

static LaurentSeriesResult essential_laurent(
    const std::shared_ptr<SymbolicExpr>& shift, int m, int order_neg,
    ComputationContext& context) {
    constexpr const char* operation = "laurent_series_full";
    auto series = SymbolicExpr::number(1);
    BigInt factorial(1);
    for (int k = 1; k <= order_neg / m; ++k) {
        auto step = context.consume_steps(1, operation);
        if (!step) {
            return LaurentSeriesResult::failure(step.error());
        }
        factorial *= BigInt(k);
        auto term = SymbolicExpr::divide(
            SymbolicExpr::number(1),
            SymbolicExpr::multiply(
                SymbolicExpr::number(factorial),
                SymbolicExpr::power(shift, SymbolicExpr::number(m * k))));
        series = SymbolicExpr::add(series, term);
    }
    return LaurentResult{
        series->simplify(), SingularityType::Essential, 0,
        SymbolicExpr::number(m == 1 ? 1 : 0)};
}
static LaurentSeriesResult laurent_series_full_impl(
    const std::shared_ptr<SymbolicExpr>&,
    const std::string&,
    const std::shared_ptr<SymbolicExpr>&,
    int, int, ComputationContext&);
static Result<void> add_fourier_harmonic(
    Integrator& integrator, const std::shared_ptr<SymbolicExpr>& f,
    const std::string& var, const std::shared_ptr<SymbolicExpr>& half_lo,
    const std::shared_ptr<SymbolicExpr>& half_hi,
    const std::shared_ptr<SymbolicExpr>& L,
    const std::shared_ptr<SymbolicExpr>& arg, bool cosine,
    std::shared_ptr<SymbolicExpr>& result, ComputationContext& context) {
    auto integrand = SymbolicExpr::multiply(
        f, cosine ? SymbolicExpr::cos(arg) : SymbolicExpr::sin(arg));
    auto integrated = integrator.integrate_def_checked(
        *integrand, var, *half_lo, *half_hi, context);
    if (!integrated) {
        return Result<void>::failure(integrated.error());
    }
    auto coefficient = SymbolicExpr::divide(
        detail::make_expression_ptr(integrated.value()), L)->simplify();
    if (!(detail::node(coefficient) && detail::node(coefficient)->is_zero())) {
        result = SymbolicExpr::add(
            result, SymbolicExpr::multiply(
                coefficient, cosine ? SymbolicExpr::cos(arg) : SymbolicExpr::sin(arg)));
    }
    return Result<void>::success();
}

struct FourierPeriod {
    std::shared_ptr<SymbolicExpr> two;
    std::shared_ptr<SymbolicExpr> half_length;
    std::shared_ptr<SymbolicExpr> lower;
    std::shared_ptr<SymbolicExpr> upper;
    std::shared_ptr<SymbolicExpr> frequency;
};

static Result<FourierPeriod> prepare_fourier_period(
    const std::shared_ptr<SymbolicExpr>& period, const std::string& var,
    ComputationContext& context) {
    /**
     * @brief 系数恒等式按输入 binary64 值精确运算，倒数与乘积也保持精确。
     * pi 的 binary64 近似值仍按该数值处理。
     */
    class ExactPeriod final : public detail::SymbolicRewriter {
    public:
        void visit(const NumberNode& node) override {
            const auto* value = std::get_if<lmmc_real_t>(&node.value());
            set_result(value && std::isfinite(*value)
                ? detail::make_node<NumberNode>(Rational::from_double(*value)) : current());
        }
    } exact_period;
    auto T = detail::make_expression_ptr(exact_period.rewrite(detail::node(period)));
    std::optional<detail::AssumptionFacts> assumed;
    if (context.assumptions()) { assumed.emplace(*context.assumptions()); }
    const auto& facts = assumed ? static_cast<const FactsQuery&>(*assumed) : detail::no_facts();
    auto positive_period = detail::query_positive_value(detail::node(T), facts, context);
    if (!positive_period) { return Result<FourierPeriod>::failure(positive_period.error()); }
    if (expression_depends_on_variable(detail::node(T), var) ||
        positive_period.value() != Tribool::True) {
        return Result<FourierPeriod>::failure(
            positive_period.value() == Tribool::False ? CasErrc::InvalidArgument : CasErrc::Inconclusive,
            "Fourier period must be a proved positive constant", "fourier_series");
    }
    auto pi = LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<VariableNode>("pi", true));
    auto two = SymbolicExpr::number(2);
    auto L = SymbolicExpr::divide(T, two);
    auto half_lo = SymbolicExpr::multiply(SymbolicExpr::number(-1), L);
    auto half_hi = L;
    auto w = SymbolicExpr::divide(SymbolicExpr::multiply(two, pi), T);
    return FourierPeriod{std::move(two), std::move(L), std::move(half_lo),
        std::move(half_hi), std::move(w)};
}

static ExpressionResult fourier_constant_coefficient(Integrator& integrator,
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const FourierPeriod& period, ComputationContext& context) {
    auto a0_integral = integrator.integrate_def_checked(
        *f, var, *period.lower, *period.upper, context);
    if (!a0_integral) { return ExpressionResult::failure(a0_integral.error()); }
    auto a0 = SymbolicExpr::divide(
        LMCAS::detail::make_expression_ptr(a0_integral.value()), period.half_length)->simplify();
    return SymbolicExpr::divide(a0, period.two);
}

static Result<void> fourier_harmonics(Integrator& integrator,
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const FourierPeriod& period, const std::shared_ptr<SymbolicExpr>& x,
    int n_terms, int parity, std::shared_ptr<SymbolicExpr>& result, ComputationContext& context) {
    for (int k = 1; k <= n_terms; ++k) {
        auto kw = SymbolicExpr::multiply(SymbolicExpr::number(k), period.frequency);
        auto arg = SymbolicExpr::multiply(kw, x);
        if (parity != -1) {
            auto added = add_fourier_harmonic(integrator, f, var, period.lower,
                period.upper, period.half_length, arg, true, result, context);
            if (!added) { return added; }
        }
        if (parity != 1) {
            auto added = add_fourier_harmonic(integrator, f, var, period.lower,
                period.upper, period.half_length, arg, false, result, context);
            if (!added) { return added; }
        }
    }
    return Result<void>::success();
}

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
        if (!access) {
            return ExpressionResult::failure(access.error());
        }

        auto prepared = prepare_fourier_period(period, var, context);
        if (!prepared) { return ExpressionResult::failure(prepared.error()); }

        int parity = detect_parity(f, var);
        Integrator integrator;
        auto x = SymbolicExpr::variable(var);

        auto constant = fourier_constant_coefficient(integrator, f, var, prepared.value(), context);
        if (!constant) { return constant; }
        auto result = std::move(constant.value());
        auto harmonics = fourier_harmonics(
            integrator, f, var, prepared.value(), x, n_terms, parity, result, context);
        if (!harmonics) { return ExpressionResult::failure(harmonics.error()); }
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
    if (!full) {
        return ExpressionResult::failure(full.error());
    }
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

static Result<void> validate_laurent_result(
    const LaurentResult& result, const std::optional<int>& supported_power,
    const std::shared_ptr<SymbolicExpr>& center, const std::string& operation) {
        if (!result.series || !LMCAS::detail::node(result.series)) {
            return Result<void>::failure(
                CasErrc::Inconclusive,
                "Laurent series could not be constructed in the supported domain",
                operation);
        }
        if (!result.residue || !LMCAS::detail::node(result.residue)) {
            return Result<void>::failure(
                CasErrc::InternalInvariant,
                "Laurent series produced a null residue expression",
                operation);
        }
        if (supported_power && center->is_zero()) {
            BigInt expected_pole_order =
                *supported_power < 0 ? -BigInt(*supported_power) : BigInt(0);
            if (BigInt(result.pole_order) != expected_pole_order) {
                return Result<void>::failure(
                    CasErrc::Inconclusive,
                    "Laurent pole order could not be verified",
                    operation);
            }
        }
    return Result<void>::success();
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
    if (!input) {
        return LaurentSeriesResult::failure(input.error());
    }
    auto order_check = validate_laurent_orders(order_neg, order_pos, operation);
    if (!order_check) {
        return LaurentSeriesResult::failure(order_check.error());
    }
    auto budget = context.consume_steps(static_cast<std::size_t>(order_neg + order_pos + 1) * 16 + 16,
                                        operation);
    if (!budget) {
        return LaurentSeriesResult::failure(budget.error());
    }
    try {
        auto simplified_input = f->simplify();
        if (!simplified_input || !LMCAS::detail::node(simplified_input)) {
            return LaurentSeriesResult::failure(
                CasErrc::Inconclusive,
                "checked Laurent input could not be simplified in the supported domain",
                operation);
        }
        auto supported_power = supported_laurent_integer_power(
            LMCAS::detail::node(simplified_input), var);
        auto built = laurent_series_full_impl(
            f, var, center, order_neg, order_pos, context);
        if (!built) {
            return built;
        }
        auto result = std::move(built.value());
        auto valid = validate_laurent_result(result, supported_power, center, operation);
        if (!valid) {
            return LaurentSeriesResult::failure(valid.error());
        }
        return LaurentSeriesResult::success(std::move(result));
    } catch (const CasError& error) {
        return LaurentSeriesResult::failure(error);
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

static Result<int> extract_laurent_regular_part(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center,
    const std::shared_ptr<SymbolicExpr>& shift,
    std::shared_ptr<SymbolicExpr>& regular, ComputationContext& context) {
    constexpr const char* operation = "laurent_series_full";
    auto center_number = std::dynamic_pointer_cast<const NumberNode>(detail::node(center));
    if (!center_number || std::holds_alternative<lmmc_real_t>(center_number->value())) {
        return Result<int>::failure(CasErrc::Inconclusive,
                                    "Laurent center must be an exact number", operation);
    }
    Rational point = std::holds_alternative<BigInt>(center_number->value())
        ? Rational(std::get<BigInt>(center_number->value()))
        : std::get<Rational>(center_number->value());
    RationalDecompositionStrategy rational;
    Polynomial<Rational> numerator(var), denominator(var);
    auto recognized = rational.extract_rational(*f, var, numerator, denominator);
    if (!recognized) {
        return Result<int>::failure(recognized.error());
    }
    if (!recognized.value()) {
        auto fn = std::dynamic_pointer_cast<const FunctionNode>(detail::node(f));
        if (fn && fn->arguments().size() == 1 &&
            (fn->type() == FunctionNode::FuncType::Sin ||
             fn->type() == FunctionNode::FuncType::Cos ||
             fn->type() == FunctionNode::FuncType::Exp)) {
            auto argument = detail::make_expression_ptr(fn->arguments()[0]);
            auto polynomial = recognize_rational_polynomial(*argument, var, context);
            if (!polynomial) {
                return Result<int>::failure(polynomial.error());
            }
            if (polynomial.value()) {
                return 0;
            }
        }
        return Result<int>::failure(CasErrc::Inconclusive,
                                    "Laurent singularity could not be proved", operation);
    }
    if (numerator.is_zero()) {
        return 0;
    }
    const Polynomial<Rational> factor({Rational(0) - point, Rational(1)}, var);
    int pole_order = 0;
    while (denominator.eval(point) == Rational(0)) {
        auto step = context.consume_steps(1, operation);
        if (!step) {
            return Result<int>::failure(step.error());
        }
        if (pole_order >= 64) {
            return Result<int>::failure(CasErrc::Inconclusive,
                                        "Laurent pole order exceeds supported range", operation);
        }
        denominator = denominator.div_mod(factor).first;
        ++pole_order;
    }
    if (pole_order > 0) {
        regular = SymbolicExpr::multiply(
            SymbolicExpr::power(shift, SymbolicExpr::number(pole_order)), f)->cancel()->simplify();
    }
    return pole_order;
}

static LaurentSeriesResult laurent_series_full_impl(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var,
    const std::shared_ptr<SymbolicExpr>& center, int order_neg, int order_pos,
    ComputationContext& context) {
    LaurentResult res{nullptr, SingularityType::Removable, 0, SymbolicExpr::number(0)};
    if (!f || !center) {
        return res;
    }

    auto x = SymbolicExpr::variable(var);
    auto shift = SymbolicExpr::add(x, SymbolicExpr::multiply(SymbolicExpr::number(-1), center));
    if (!expression_depends_on_variable(detail::node(center), var)) {
        if (auto m = essential_exponent(f, shift)) {
            return essential_laurent(shift, *m, order_neg, context);
        }
    }

    std::shared_ptr<SymbolicExpr> regular = f;
    auto extracted = extract_laurent_regular_part(
        f, var, center, shift, regular, context);
    if (!extracted) {
        return LaurentSeriesResult::failure(extracted.error());
    }
    const int pole_order = extracted.value();

    /// 对 regular(解析部分)做 Taylor 展开
    int taylor_order = order_pos + pole_order + 1;
    if (taylor_order < 1) {
        taylor_order = 1;
    }
    auto taylor = regular->series(var, center, taylor_order);
    if (!taylor) {
        return res;
    }

    /// Laurent = taylor / (x-c)^pole_order
    std::shared_ptr<SymbolicExpr> laurent = taylor;
    if (pole_order > 0) {
        auto powm = SymbolicExpr::power(shift, SymbolicExpr::number(pole_order));
        laurent = SymbolicExpr::divide(taylor, powm)->simplify();
    }

    res.series = laurent->simplify();
    res.pole_order = pole_order;
    if (pole_order == 0) {
        res.singularity = SingularityType::Removable;
    } else {
        res.singularity = SingularityType::Pole;
    }

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
        if (!val) {
            val = deriv->substitute(var, center);
        }
        auto fact = BigInt::factorial_checked(BigInt(pole_order - 1), context);
        if (!fact) {
            return LaurentSeriesResult::failure(fact.error());
        }
        res.residue = SymbolicExpr::divide(val, SymbolicExpr::number(fact.value()))->simplify();
    }

    return res;
}

std::shared_ptr<SymbolicExpr> asymptotic_expand(
    const std::shared_ptr<SymbolicExpr>& f, const std::string& var, int order) {
    if (!f || order < 0) {
        return nullptr;
    }
    /// 通过 x = 1/t 替换,在 t=0 处做 Taylor 展开,再回代 t = 1/x,
    /// 得到按 x 递减幂次的渐近展开.
    auto x = SymbolicExpr::variable(var);
    std::string tname = var + "__asym_t";
    auto t = SymbolicExpr::variable(tname);
    auto inv_t = SymbolicExpr::divide(SymbolicExpr::number(1), t);
    auto g = f->substitute(var, inv_t);
    if (!g) {
        return nullptr;
    }
    g = g->simplify();
    auto taylor_t = g->series(tname, SymbolicExpr::number(0), order + 1);
    if (!taylor_t) {
        return nullptr;
    }
    /// 回代 t = 1/x
    auto inv_x = SymbolicExpr::divide(SymbolicExpr::number(1), x);
    auto back = taylor_t->substitute(tname, inv_x);
    if (!back) {
        return nullptr;
    }
    return back->simplify();
}
}
