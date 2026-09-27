#include "test_common.hpp"
#include "symbolic_ode_engine.hpp"
#include "poly_utils.hpp"
#include "numeric_evaluation.hpp"
#include <limits>

using namespace LMCAS;

static ODESingularityType checked_singularity(
    const std::shared_ptr<SymbolicExpr> &p,
    const std::shared_ptr<SymbolicExpr> &q,
    const std::shared_ptr<SymbolicExpr> &point,
    const std::string &variable) {
    auto result =
        classify_singular_point_checked(p, q, point, variable);
    EXPECT_TRUE((result.has_value())) << "checked singularity classification succeeds";
    return result ? result.value() : ODESingularityType::IrregularSingular;
}

TEST(OdeEngineAdvanced, VariationSecant) {
    /**
     * @brief 用常数变易法构造 sec(x) 的特解。
     *
     * 齐次解 y_1=cos(x)、y_2=sin(x)，非齐次项 g(x)=sec(x)=1/cos(x)。
     * Wronskian 为 W=cos^2(x)+sin^2(x)=1。
     * u_1'=-sin(x)/cos(x)=-tan(x)，u_2'=1；积分得 u_1=ln|cos(x)|、u_2=x。
     * 特解 y_p=cos(x)*ln|cos(x)|+x*sin(x)。
     */
    auto x = SymbolicExpr::variable("x");
    auto y1 = SymbolicExpr::cos(x);
    auto y2 = SymbolicExpr::sin(x);
    auto g = SymbolicExpr::divide(SymbolicExpr::number(1), SymbolicExpr::cos(x));

    auto checked = solve_variation_of_parameters_checked(y1, y2, g, "x");
    ASSERT_TRUE(checked);
    const auto &sol = checked.value();
    ASSERT_TRUE(sol.general_solution);
    EXPECT_TRUE(sol.constants.empty())
        << "variation of parameters returns a particular solution";

    auto first_derivative = sol.general_solution->differentiate("x");
    ASSERT_TRUE(first_derivative);
    auto second_derivative = first_derivative->differentiate("x");
    ASSERT_TRUE(second_derivative);
    auto lhs = SymbolicExpr::add(second_derivative, sol.general_solution);
    for (double point : {-1.0, -0.5, 0.5, 1.0}) {
        auto value = evaluate_numeric(*lhs, {{"x", point}});
        ASSERT_TRUE(value);
        ASSERT_TRUE(value.value().is_finite());
        EXPECT_NEAR(value.value().value, 1.0 / std::cos(point), 1e-9)
            << "variation solution satisfies y''+y=sec(x) at x=" << point;
    }
}

TEST(OdeEngineAdvanced, VariationExponential) {
    /**
     * @brief 用常数变易法构造 e^x 的特解。
     *
     * 齐次解 y_1=e^x、y_2=e^(-x)，非齐次项 g(x)=e^x。
     * Wronskian 为 W=e^x*(-e^(-x))-e^(-x)*e^x=-2。
     * u_1'=1/2，u_2'=-e^(2x)/2；积分得 u_1=x/2、u_2=-e^(2x)/4。
     * 特解 y_p=(x/2)*e^x-e^x/4。
     */
    auto x = SymbolicExpr::variable("x");
    auto y1 = SymbolicExpr::exp(x);
    auto y2 = SymbolicExpr::exp(
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x));
    auto g = SymbolicExpr::exp(x);

    auto checked = solve_variation_of_parameters_checked(
        y1, y2, g, "x");
    ASSERT_TRUE(checked);
    auto solution = checked.value().general_solution;
    ASSERT_TRUE(solution);
    auto first_derivative = solution->differentiate("x");
    ASSERT_TRUE(first_derivative);
    auto second_derivative = first_derivative->differentiate("x");
    ASSERT_TRUE(second_derivative);
    auto lhs = SymbolicExpr::add(
        second_derivative,
        SymbolicExpr::multiply(SymbolicExpr::number(-1), solution));
    for (double point : {-1.0, 0.0, 1.0}) {
        auto value = evaluate_numeric(*lhs, {{"x", point}});
        ASSERT_TRUE(value);
        ASSERT_TRUE(value.value().is_finite());
        EXPECT_NEAR(value.value().value, std::exp(point), 1e-9)
            << "variation solution satisfies y''-y=exp(x)";
    }
}

TEST(OdeEngineAdvanced, VariationNull) {
    auto sol = solve_variation_of_parameters_checked(
        nullptr, nullptr, nullptr, "x");
    EXPECT_TRUE((!sol && sol.error().code == CasErrc::InvalidArgument)) << "null VoP inputs are InvalidArgument";
}

TEST(OdeEngineAdvanced, OrdinaryPoint) {
    auto p = SymbolicExpr::number(0);
    auto q = SymbolicExpr::number(1);
    auto x0 = SymbolicExpr::number(0);

    auto type = checked_singularity(p, q, x0, "x");
    EXPECT_TRUE((type == ODESingularityType::Ordinary)) << "y''+y=0 at x=0 is ordinary";
}

TEST(OdeEngineAdvanced, RegularSingularPoint) {
    /**
     * @brief 检验 n=0 的 Bessel 方程在原点的正则奇点。
     *
     * 方程 x^2*y''+x*y'+(x^2-n^2)*y=0 归一化为 y''+(1/x)*y'+(1-n^2/x^2)*y=0。
     * p=1/x；取 n=0 时 q=1，x*p=1、x^2*q=x^2 在原点均有限。
     */
    auto x = SymbolicExpr::variable("x");
    auto p = SymbolicExpr::divide(SymbolicExpr::number(1), x);
    auto q = SymbolicExpr::number(1);
    auto x0 = SymbolicExpr::number(0);

    auto type = checked_singularity(p, q, x0, "x");
    EXPECT_TRUE((type == ODESingularityType::RegularSingular)) << "Bessel at x=0 is regular singular";
}

TEST(OdeEngineAdvanced, SymbolicSingularity) {
    auto parameter = SymbolicExpr::variable("a");
    auto result = classify_singular_point_checked(
        parameter, SymbolicExpr::number(1),
        SymbolicExpr::number(0), "x");
    EXPECT_TRUE((!result)) << "symbolic point values are not mistaken for poles";
    EXPECT_TRUE((!result && result.error().code == CasErrc::Inconclusive)) << "symbolic singularity classification reports Inconclusive";
}

TEST(OdeEngineAdvanced, OrdinarySeries) {
    auto p = SymbolicExpr::number(0);
    auto q = SymbolicExpr::number(1);
    auto x0 = SymbolicExpr::number(0);

    auto checked = solve_frobenius_checked(p, q, x0, "x", 6);
    ASSERT_TRUE(checked);
    const auto &sol = checked.value();
    ASSERT_TRUE(sol.series_solution);
    EXPECT_EQ(sol.point_type, ODESingularityType::Ordinary);

    const std::vector<double> expected_coefficients = {
        1.0, 0.0, -0.5, 0.0, 1.0 / 24.0, 0.0, -1.0 / 720.0};
    auto derivative = sol.series_solution;
    double factorial = 1.0;
    for (std::size_t order = 0; order < expected_coefficients.size(); ++order) {
        if (order != 0) {
            derivative = derivative->differentiate("x");
            ASSERT_TRUE(derivative);
            factorial *= static_cast<double>(order);
        }
        auto value = evaluate_numeric(*derivative, {{"x", 0.0}});
        ASSERT_TRUE(value);
        ASSERT_TRUE(value.value().is_finite());
        EXPECT_NEAR(value.value().value / factorial,
                    expected_coefficients[order], 1e-12)
            << "coefficient of x^" << order;
    }

    auto first_derivative = sol.series_solution->differentiate("x");
    ASSERT_TRUE(first_derivative);
    auto second_derivative = first_derivative->differentiate("x");
    ASSERT_TRUE(second_derivative);
    auto residual = SymbolicExpr::add(second_derivative, sol.series_solution);
    for (double point : {-1.0, -0.5, 0.5, 1.0}) {
        auto value = evaluate_numeric(*residual, {{"x", point}});
        ASSERT_TRUE(value);
        ASSERT_TRUE(value.value().is_finite());
        EXPECT_NEAR(value.value().value, -std::pow(point, 6) / 720.0, 1e-12)
            << "the order-six truncation leaves the expected remainder";
    }
}

TEST(OdeEngineAdvanced, EulerSeries) {
    /**
     * @brief Euler 方程的指标根为 ±1。
     *
     * 归一化方程为 y''+(1/x)*y'-(1/x^2)*y=0，p=1/x、q=-1/x^2。
     * P_0=lim x*p=1，Q_0=lim x^2*q=-1；指标方程 r(r-1)+r-1=r^2-1=0。
     */
    auto x = SymbolicExpr::variable("x");
    auto p = SymbolicExpr::divide(SymbolicExpr::number(1), x);
    auto q = SymbolicExpr::divide(SymbolicExpr::number(-1),
                                  SymbolicExpr::power(x, SymbolicExpr::number(2)));
    auto x0 = SymbolicExpr::number(0);

    auto checked = solve_frobenius_checked(p, q, x0, "x", 4);
    ASSERT_TRUE(checked.has_value())
        << static_cast<int>(checked.error().code) << " "
        << checked.error().operation << ": " << checked.error().message;
    const auto &sol = checked.value();
    ASSERT_TRUE(sol.series_solution) << "Frobenius Euler eq has series solution";
    EXPECT_EQ(sol.point_type, ODESingularityType::RegularSingular)
        << "Euler eq at x=0 is regular singular";
    ASSERT_EQ(sol.indicial_roots.size(), 2u) << "has two indicial roots";
    const double expected_roots[] = {1.0, -1.0};
    for (size_t i = 0; i < sol.indicial_roots.size(); ++i) {
        SCOPED_TRACE(i);
        EXPECT_TRUE(std::isfinite(sol.indicial_roots[i]));
        EXPECT_TRUE(std::isfinite(expected_roots[i]));
        EXPECT_NEAR(sol.indicial_roots[i], expected_roots[i], 1e-9);
    }

    for (double point : {-2.0, -0.5, 0.5, 2.0}) {
        auto value = evaluate_numeric(*sol.series_solution, {{"x", point}});
        ASSERT_TRUE(value.has_value()) << "Euler series evaluates on both punctured sides";
        EXPECT_TRUE(value.value().is_finite());
        EXPECT_TRUE(std::isfinite(value.value().value));
        EXPECT_TRUE(std::isfinite(point));
        EXPECT_NEAR(value.value().value, point, 1e-9)
            << "Euler normalized larger-root solution is x";
    }
}

TEST(OdeEngineAdvanced, IrregularSingularity) {
    auto x = SymbolicExpr::variable("x");
    auto p = SymbolicExpr::number(0);
    auto q = SymbolicExpr::divide(SymbolicExpr::number(1),
                                  SymbolicExpr::power(x, SymbolicExpr::number(3)));
    auto x0 = SymbolicExpr::number(0);

    auto type = checked_singularity(p, q, x0, "x");
    EXPECT_TRUE((type == ODESingularityType::IrregularSingular)) << "1/x^3 coefficient gives irregular singular point";

    auto sol = solve_frobenius_checked(p, q, x0, "x", 4);
    EXPECT_TRUE((!sol &&
                 (sol.error().code == CasErrc::DomainError ||
                  sol.error().code == CasErrc::Inconclusive)))
        << "irregular singular Frobenius input is explicitly rejected";
}

TEST(OdeEngineAdvanced, VariationChecked) {
    auto x = SymbolicExpr::variable("x");
    auto y1 = SymbolicExpr::cos(x);
    auto y2 = SymbolicExpr::sin(x);
    auto g = SymbolicExpr::sin(x);

    auto sol = solve_variation_of_parameters_checked(y1, y2, g, "x");
    ASSERT_TRUE(sol)
        << "checked variation of parameters succeeds for independent solutions";
    ASSERT_TRUE(sol.value().general_solution);
    EXPECT_EQ(sol.value().method_used, ODEType::HigherOrder_ConstCoeff)
        << "checked variation reports the higher-order method family";

    auto dependent = solve_variation_of_parameters_checked(x, x, g, "x");
    EXPECT_TRUE((!dependent.has_value())) << "checked variation of parameters rejects zero Wronskian as unsupported";
    EXPECT_TRUE((dependent.error().code == CasErrc::Inconclusive)) << "checked variation of parameters reports Inconclusive for zero Wronskian";

    std::shared_ptr<SymbolicExpr> null_expr;
    auto invalid = solve_variation_of_parameters_checked(null_expr, y2, g, "x");
    EXPECT_TRUE((!invalid.has_value())) << "checked variation of parameters rejects null input";
    EXPECT_TRUE((invalid.error().code == CasErrc::InvalidArgument)) << "checked variation of parameters reports InvalidArgument for null input";
}

TEST(OdeEngineAdvanced, FrobeniusSymbolic) {
    auto x = SymbolicExpr::variable("x");
    auto parameter = SymbolicExpr::variable("a");
    auto p = SymbolicExpr::multiply(parameter, x);
    auto q = SymbolicExpr::number(1);
    auto x0 = SymbolicExpr::number(0);
    auto sol = solve_frobenius_checked(p, q, x0, "x", 4);
    EXPECT_TRUE((!sol.has_value())) << "checked Frobenius rejects coefficients that cannot be lowered numerically";
    if (!sol) {
        EXPECT_TRUE((sol.error().code == CasErrc::Inconclusive)) << "checked Frobenius reports Inconclusive for symbolic Taylor coefficients";
    }
}

TEST(OdeEngineAdvanced, FrobeniusInvalid) {
    auto x = SymbolicExpr::variable("x");
    auto p = SymbolicExpr::number(0);
    auto q = SymbolicExpr::divide(SymbolicExpr::number(1),
                                  SymbolicExpr::power(x, SymbolicExpr::number(3)));
    auto x0 = SymbolicExpr::number(0);
    auto sol = solve_frobenius_checked(p, q, x0, "x", 4);
    EXPECT_TRUE((!sol.has_value())) << "checked Frobenius rejects irregular singular point";
    EXPECT_TRUE((sol.error().code == CasErrc::Inconclusive)) << "checked Frobenius reports Inconclusive for irregular singular point";

    auto bad_order = solve_frobenius_checked(p, q, x0, "x", -1);
    EXPECT_TRUE((!bad_order.has_value())) << "checked Frobenius rejects negative truncation order";
    EXPECT_TRUE((bad_order.error().code == CasErrc::InvalidArgument)) << "checked Frobenius reports InvalidArgument for negative order";

    std::shared_ptr<SymbolicExpr> null_expr;
    auto invalid = solve_frobenius_checked(null_expr, q, x0, "x", 4);
    EXPECT_TRUE((!invalid.has_value())) << "checked Frobenius rejects null coefficient";
    EXPECT_TRUE((invalid.error().code == CasErrc::InvalidArgument)) << "checked Frobenius reports InvalidArgument for null coefficient";
}

TEST(OdeEngineAdvanced, AdvancedContext) {
    auto x = SymbolicExpr::variable("x");
    auto y1 = SymbolicExpr::cos(x);
    auto y2 = SymbolicExpr::sin(x);
    auto g = SymbolicExpr::sin(x);

    LMCAS::CancellationToken cancellation;
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    cancellation.cancel();
    auto cancelled = solve_variation_of_parameters_checked(
        y1, y2, g, "x", cancelled_context);
    EXPECT_TRUE((!cancelled.has_value())) << "checked variation of parameters observes cancellation";
    EXPECT_TRUE((cancelled.error().code == CasErrc::Cancelled)) << "checked variation of parameters reports Cancelled";

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto p = SymbolicExpr::number(0);
    auto q = SymbolicExpr::number(1);
    auto x0 = SymbolicExpr::number(0);
    auto limited = solve_frobenius_checked(p, q, x0, "x", 6, limited_context);
    EXPECT_TRUE((!limited.has_value())) << "checked Frobenius observes exhausted step budget";
    EXPECT_TRUE((limited.error().code == CasErrc::ResourceLimit)) << "checked Frobenius reports ResourceLimit";
}
