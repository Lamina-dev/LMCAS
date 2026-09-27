#include "test_common.hpp"
#include "symbolic_ode.hpp"
#include "numeric_evaluation.hpp"
#include "internal/expression_analysis.hpp"

using namespace LMCAS;

static void expect_separable_first_integral(
    const std::shared_ptr<SymbolicExpr>& solution,
    const std::shared_ptr<SymbolicExpr>& rhs,
    double (*expected)(double, double)) {
    ASSERT_NE(solution, nullptr);
    auto residual = SymbolicExpr::add(
        solution->differentiate("x"),
        SymbolicExpr::multiply(rhs, solution->differentiate("y")));
    EXPECT_TRUE(test_proved_equivalent(residual, SymbolicExpr::number(0)));
    const double points[][2] = {{1.0, 2.0}, {2.0, 3.0}, {3.0, 1.0}, {0.5, 4.0}};
    for (const auto& point : points) {
        SCOPED_TRACE(::testing::Message() << "x=" << point[0] << ", y=" << point[1]);
        const NumericBindings values{{"x", point[0]}, {"y", point[1]}};
        auto primitive = evaluate_numeric(*solution, values);
        ASSERT_TRUE(primitive.has_value());
        ASSERT_TRUE(std::isfinite(primitive.value().value));
        EXPECT_NEAR(primitive.value().value, expected(point[0], point[1]), 1e-12);
    }
}

TEST(SymbolicOde, SeparableXOverY) {
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        // rhs = x / y
        auto rhs = SymbolicExpr::divide(x, y);
        auto sol = LMCAS::solve_separable_ode(rhs, "x", "y");
        ASSERT_NE(sol, nullptr);
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "y"));
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "x"));
        expect_separable_first_integral(sol, rhs, [](double x, double y) {
            return (y * y - x * x) / 2.0;
        });
    }
}

TEST(SymbolicOde, SeparableXY) {
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        // rhs = x * y
        auto rhs = SymbolicExpr::multiply(x, y);
        auto sol = LMCAS::solve_separable_ode(rhs, "x", "y");
        ASSERT_NE(sol, nullptr);
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "y"));
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "x"));
        auto at_one = sol->substitute("x", SymbolicExpr::number(1));
        ASSERT_NE(at_one, nullptr);
        auto derivative = at_one->differentiate("y");
        ASSERT_NE(derivative, nullptr);
        auto at_two = evaluate_numeric(*derivative, {{"y", 2.0}});
        auto at_four = evaluate_numeric(*derivative, {{"y", 4.0}});
        ASSERT_TRUE(at_two.has_value());
        ASSERT_TRUE(at_four.has_value());
        EXPECT_TRUE(std::isfinite(at_two.value().value));
        EXPECT_TRUE(std::isfinite(at_four.value().value));
        EXPECT_NEAR(at_two.value().value, 0.5, 1e-12);
        EXPECT_NEAR(at_four.value().value, 0.25, 1e-12);
        expect_separable_first_integral(sol, rhs, [](double x, double y) {
            return std::log(y) - x * x / 2.0;
        });
    }
}

TEST(SymbolicOde, SeparableYOverX) {
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        // rhs = y / x
        auto rhs = SymbolicExpr::divide(y, x);
        auto sol = LMCAS::solve_separable_ode(rhs, "x", "y");
        ASSERT_NE(sol, nullptr);
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "y"));
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "x"));
        auto at_one = sol->substitute("x", SymbolicExpr::number(1));
        ASSERT_NE(at_one, nullptr);
        auto derivative = at_one->differentiate("y");
        ASSERT_NE(derivative, nullptr);
        auto at_two = evaluate_numeric(*derivative, {{"y", 2.0}});
        auto at_four = evaluate_numeric(*derivative, {{"y", 4.0}});
        ASSERT_TRUE(at_two.has_value());
        ASSERT_TRUE(at_four.has_value());
        EXPECT_TRUE(std::isfinite(at_two.value().value));
        EXPECT_TRUE(std::isfinite(at_four.value().value));
        EXPECT_NEAR(at_two.value().value, 0.5, 1e-12);
        EXPECT_NEAR(at_four.value().value, 0.25, 1e-12);
        expect_separable_first_integral(sol, rhs, [](double x, double y) {
            return std::log(y) - std::log(x);
        });
    }
}

TEST(SymbolicOde, SeparableConstantAndPureVariables) {
    auto zero = SymbolicExpr::number(0);
    expect_separable_first_integral(
        solve_separable_ode(zero, "x", "y"), zero,
        [](double, double y) { return y; });

    auto constant = SymbolicExpr::number(3);
    expect_separable_first_integral(
        solve_separable_ode(constant, "x", "y"), constant,
        [](double x, double y) { return y - 3.0 * x; });

    auto x = SymbolicExpr::variable("x");
    expect_separable_first_integral(
        solve_separable_ode(x, "x", "y"), x,
        [](double x, double y) { return y - x * x / 2.0; });

    auto y = SymbolicExpr::variable("y");
    expect_separable_first_integral(
        solve_separable_ode(y, "x", "y"), y,
        [](double x, double y) { return std::log(y) - x; });

    auto inverse_square = SymbolicExpr::power(y, SymbolicExpr::number(-2));
    expect_separable_first_integral(
        solve_separable_ode(inverse_square, "x", "y"), inverse_square,
        [](double x, double y) { return y * y * y / 3.0 - x; });
}

TEST(SymbolicOde, SeparableRejectsCoupledRightHandSides) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    EXPECT_EQ(solve_separable_ode(SymbolicExpr::add(x, y), "x", "y"), nullptr);
    EXPECT_EQ(solve_separable_ode(
        SymbolicExpr::sin(SymbolicExpr::multiply(x, y)), "x", "y"), nullptr);
}

TEST(SymbolicOde, Linear1Ode) {
    {
        // dy/dx + 2*y = 0 => P(x) = 2, Q(x) = 0
        auto Px = SymbolicExpr::number(2);
        auto Qx = SymbolicExpr::number(0);
        auto sol = LMCAS::solve_linear1_ode(Px, Qx, "x", "y");
        ASSERT_NE(sol, nullptr);
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "C"));
    }

    {
        // dy/dx + x*y = x => P(x) = x, Q(x) = x
        auto x = SymbolicExpr::variable("x");
        auto Px = x;
        auto Qx = SymbolicExpr::variable("x");
        auto sol = LMCAS::solve_linear1_ode(Px, Qx, "x", "y");
        ASSERT_NE(sol, nullptr);
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "C"));
    }

    {
        // dy/dx + (1/x)*y = x^2 => P(x) = 1/x, Q(x) = x^2
        auto x = SymbolicExpr::variable("x");
        auto one = SymbolicExpr::number(1);
        auto Px = SymbolicExpr::divide(one, x);
        auto Qx = SymbolicExpr::power(SymbolicExpr::variable("x"), SymbolicExpr::number(2));
        auto sol = LMCAS::solve_linear1_ode(Px, Qx, "x", "y");
        ASSERT_NE(sol, nullptr);
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "C"));
    }
}

TEST(SymbolicOde, Linear2Ode) {
    {
        // a=1, b=-3, c=2, f(x)=0
        // Characteristic: r^2 - 3r + 2 = 0 => r=1, r=2
        auto fx = SymbolicExpr::number(0);
        auto sol = LMCAS::solve_linear2_ode(1, -3, 2, fx, "x", "y");
        ASSERT_NE(sol, nullptr);
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "C1"));
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "C2"));
    }

    {
        // a=1, b=-2, c=1, f(x)=0
        // Characteristic: r^2 - 2r + 1 = 0 => r=1 (double)
        auto fx = SymbolicExpr::number(0);
        auto sol = LMCAS::solve_linear2_ode(1, -2, 1, fx, "x", "y");
        ASSERT_NE(sol, nullptr);
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "C1"));
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "C2"));
    }

    {
        // a=1, b=0, c=1, f(x)=0
        // Characteristic: r^2 + 1 = 0 => r=+/-i
        auto fx = SymbolicExpr::number(0);
        auto sol = LMCAS::solve_linear2_ode(1, 0, 1, fx, "x", "y");
        ASSERT_NE(sol, nullptr);
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "C1"));
        EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "C2"));
    }
}

TEST(SymbolicOde, Linear2Nonhomogeneous) {
    {
        // a=1, b=-3, c=2, f(x) = e^(3x)
        /// 旧版接口对支持域之外的非齐次输入抛出 std::logic_error.
        auto x = SymbolicExpr::variable("x");
        auto three_x = SymbolicExpr::multiply(SymbolicExpr::number(3), x);
        auto fx = SymbolicExpr::exp(three_x);
        bool threw = false;
        try {
            auto sol = LMCAS::solve_linear2_ode(1, -3, 2, fx, "x", "y");
            ASSERT_NE(sol, nullptr);
            EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "C1"));
            EXPECT_TRUE(expression_depends_on_variable(LMCAS::detail::node(*sol), "C2"));
        } catch (const std::logic_error &) {
            // Non-homogeneous case is not yet implemented - this is expected
            threw = true;
        }
        EXPECT_TRUE((threw)) << "legacy linear2 nonhomogeneous throws instead of returning a false solution";
    }
}

static void test_linear2_homogeneous_checked() {
    {
        auto fx = SymbolicExpr::number(0);
        auto result = LMCAS::solve_linear2_ode_checked(1, -3, 2, fx, "x", "y");
        ASSERT_TRUE(result.has_value());
        ASSERT_NE(result.value(), nullptr);
        EXPECT_TRUE(expression_depends_on_variable(
            LMCAS::detail::node(*result.value()), "C1"));
        EXPECT_TRUE(expression_depends_on_variable(
            LMCAS::detail::node(*result.value()), "C2"));
    }
}

static void test_linear2_repeated_zero() {
    {
        auto result = LMCAS::solve_linear2_ode_checked(
            1, 0, 0, SymbolicExpr::number(0), "x", "y");
        ASSERT_TRUE((result.has_value())) << "checked linear2 verifies the repeated zero characteristic root";
        if (result) {
            double initial_values[2] = {};
            double initial_slopes[2] = {};
            for (int i = 0; i < 2; ++i) {
                auto basis = result.value()
                                 ->substitute("C1", SymbolicExpr::number(i == 0 ? 1 : 0))
                                 ->substitute("C2", SymbolicExpr::number(i == 1 ? 1 : 0));
                auto derivative = basis->differentiate("x");
                EXPECT_TRUE((derivative->differentiate("x")->simplify()->is_zero())) << "each repeated-zero-root basis solves y'' = 0 exactly";
                auto value = LMCAS::evaluate_numeric(*basis, {{"x", 0.0}});
                auto slope = LMCAS::evaluate_numeric(*derivative, {{"x", 0.0}});
                EXPECT_TRUE((value.has_value() && slope.has_value())) << "repeated-zero-root basis has finite initial data";
                if (value && slope) {
                    initial_values[i] = value.value().value;
                    initial_slopes[i] = slope.value().value;
                }
            }
            const double wronskian =
                initial_values[0] * initial_slopes[1] -
                initial_values[1] * initial_slopes[0];
            EXPECT_TRUE((std::isfinite(wronskian) && std::abs(wronskian) > 1e-12)) << "repeated-zero-root solution spans arbitrary value and slope";
        }
    }
}

static void test_linear2_coefficient_scaling() {
    {
        const double scale = 1.0e200;
        auto result = LMCAS::solve_linear2_ode_checked(
            scale, 3.0 * scale, 2.0 * scale,
            SymbolicExpr::number(0), "x", "y");
        ASSERT_TRUE((result.has_value())) << "checked linear2 is invariant under finite coefficient scaling";
        if (result) {
            bool found_minus_one = false;
            bool found_minus_two = false;
            for (const char *selected : {"C1", "C2"}) {
                auto basis = result.value()
                                 ->substitute("C1", SymbolicExpr::number(
                                                        std::string(selected) == "C1" ? 1 : 0))
                                 ->substitute("C2", SymbolicExpr::number(
                                                        std::string(selected) == "C2" ? 1 : 0));
                auto derivative = LMCAS::evaluate_numeric(
                    *basis->differentiate("x"), {{"x", 0.0}});
                if (derivative) {
                    found_minus_one |=
                        std::abs(derivative.value().value + 1.0) < 1.0e-12;
                    found_minus_two |=
                        std::abs(derivative.value().value + 2.0) < 1.0e-12;
                }
            }
            EXPECT_TRUE((found_minus_one && found_minus_two)) << "scaled linear2 preserves both characteristic roots";
        }
    }
}

static void test_linear2_forcing_domain() {
    {
        auto x = SymbolicExpr::variable("x");
        auto three_x = SymbolicExpr::multiply(SymbolicExpr::number(3), x);
        auto fx = SymbolicExpr::exp(three_x);
        auto result = LMCAS::solve_linear2_ode_checked(1, -3, 2, fx, "x", "y");
        EXPECT_TRUE((!result.has_value())) << "checked linear2 nonhomogeneous is not a success";
        EXPECT_TRUE((result.error().code == LMCAS::CasErrc::Inconclusive)) << "checked linear2 nonhomogeneous reports Inconclusive";
    }
}

static void test_linear2_invalid_inputs() {
    {
        auto result = LMCAS::solve_linear2_ode_checked(1, -3, 2, nullptr, "x", "y");
        EXPECT_TRUE((!result.has_value())) << "checked linear2 rejects null forcing expression";
        EXPECT_TRUE((result.error().code == LMCAS::CasErrc::InvalidArgument)) << "checked linear2 null forcing reports InvalidArgument";
    }

    {
        auto fx = SymbolicExpr::number(0);
        auto result = LMCAS::solve_linear2_ode_checked(1, -3, 2, fx, "", "y");
        EXPECT_TRUE((!result.has_value())) << "checked linear2 rejects empty independent variable";
        EXPECT_TRUE((result.error().code == LMCAS::CasErrc::InvalidArgument)) << "checked linear2 empty variable reports InvalidArgument";
    }
}

static void test_linear2_context_errors() {
    {
        LMCAS::CancellationToken cancellation;
        LMCAS::ComputationContext context({}, cancellation);
        cancellation.cancel();
        auto result = LMCAS::solve_linear2_ode_checked(
            1, -3, 2, SymbolicExpr::number(0), "x", "y", context);
        EXPECT_TRUE((!result.has_value())) << "checked linear2 observes cancellation";
        EXPECT_TRUE((result.error().code == LMCAS::CasErrc::Cancelled)) << "checked linear2 cancellation reports Cancelled";
    }

    {
        LMCAS::ResourceLimits limits;
        limits.max_steps = 0;
        LMCAS::ComputationContext context(limits);
        auto result = LMCAS::solve_linear2_ode_checked(
            1, -3, 2, SymbolicExpr::number(0), "x", "y", context);
        EXPECT_TRUE((!result.has_value())) << "checked linear2 observes exhausted step budget";
        EXPECT_TRUE((result.error().code == LMCAS::CasErrc::ResourceLimit)) << "checked linear2 exhausted budget reports ResourceLimit";
    }
}

TEST(SymbolicOde, Linear2CheckedContracts) {
    test_linear2_homogeneous_checked();
    test_linear2_repeated_zero();
    test_linear2_coefficient_scaling();
    test_linear2_forcing_domain();
    test_linear2_invalid_inputs();
    test_linear2_context_errors();
}
