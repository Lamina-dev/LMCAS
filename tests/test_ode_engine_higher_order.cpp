#include "test_common.hpp"
#include "symbolic_ode_engine.hpp"
#include "poly_utils.hpp"
#include "numeric_evaluation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

using namespace LMCAS;
static ExprPtr derivative_order(ExprPtr expression, std::size_t order)
{
    for (std::size_t current = 0; current < order; ++current) {
        expression = expression->differentiate("x");
        if (!expression) {
            return nullptr;
        }
    }
    return expression;
}

static void expect_constant_coefficient_solution(
    const ODESolution& solution, const std::vector<double>& coefficients,
    const ExprPtr& forcing)
{
    ASSERT_TRUE(solution.general_solution);
    ASSERT_EQ(
        solution.constants.size(), coefficients.size() - 1);
    for (std::size_t selected = 0;
         selected <= solution.constants.size(); ++selected) {
        auto basis = solution.general_solution;
        for (std::size_t index = 0;
             index < solution.constants.size(); ++index) {
            basis = basis->substitute(
                solution.constants[index],
                SymbolicExpr::number(index == selected ? 1 : 0));
        }

        ExprPtr lhs = SymbolicExpr::number(0);
        const std::size_t order = coefficients.size() - 1;
        for (std::size_t index = 0; index < coefficients.size(); ++index) {
            auto derivative =
                derivative_order(basis, order - index);
            ASSERT_TRUE(derivative);
            lhs = SymbolicExpr::add(
                lhs, SymbolicExpr::multiply(
                         SymbolicExpr::number(coefficients[index]),
                         derivative));
        }
        auto residual = SymbolicExpr::add(
            lhs, SymbolicExpr::multiply(
                     SymbolicExpr::number(-1),
                     forcing ? forcing : SymbolicExpr::number(0)));
        EXPECT_TRUE(test_proved_equivalent(
            residual->simplify(), SymbolicExpr::number(0)));
    }
}


TEST(OdeEngineHigherOrder, ConstantCoefficients) {
    const std::vector<double> coefficients{1.0, 0.0, 1.0};
    auto result = solve_higher_order_ode_checked(
        coefficients, nullptr, "x", "y");
    ASSERT_TRUE(result.has_value())
        << "checked homogeneous constant-coefficient ODE succeeds";
    EXPECT_EQ(
        result.value().method_used, ODEType::HigherOrder_ConstCoeff);
    EXPECT_EQ(
        result.value().constants,
        (std::vector<std::string>{"C1", "C2"}));
    expect_constant_coefficient_solution(
        result.value(), coefficients, nullptr);
}

TEST(OdeEngineHigherOrder, RepeatedRoot) {
    auto repeated = solve_higher_order_ode_checked(
        {1.0, -5.0, 10.0, -10.0, 5.0, -1.0},
        nullptr, "x", "y");
    ASSERT_TRUE(repeated.has_value())
        << "checked fifth-order repeated-root ODE succeeds";
    ASSERT_TRUE(repeated.value().general_solution);
    EXPECT_EQ(
        repeated.value().constants,
        (std::vector<std::string>{"C1", "C2", "C3", "C4", "C5"}));
}

TEST(OdeEngineHigherOrder, NearIntegerRoot) {
    const double root = 1.0 + 1.0e-11;
    auto result = solve_higher_order_ode_checked(
        {1.0, -root}, nullptr, "x", "y");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(result.value().general_solution);
    auto basis = result.value().general_solution->substitute(
        "C1", SymbolicExpr::number(1));
    auto value = evaluate_numeric(*basis, {{"x", 0.0}});
    auto derivative = evaluate_numeric(
        *basis->differentiate("x"), {{"x", 0.0}});
    ASSERT_TRUE(value.has_value());
    ASSERT_TRUE(derivative.has_value());
    EXPECT_DOUBLE_EQ(value.value().value, 1.0);
    EXPECT_EQ(derivative.value().value, root);
    EXPECT_NEAR(
        derivative.value().value - root * value.value().value,
        0.0, 1.0e-14);
}

TEST(OdeEngineHigherOrder, CloseDistinctQuadraticRoots) {
    const double delta = std::ldexp(1.0, -24);
    const double constant = 1.0 - std::ldexp(1.0, -48);
    auto result = solve_higher_order_ode_checked(
        {1.0, -2.0, constant}, nullptr, "x", "y");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(result.value().general_solution);
    ASSERT_EQ(result.value().constants.size(), 2u);

    std::vector<double> exponents;
    for (const auto& selected : result.value().constants) {
        auto basis = result.value().general_solution;
        for (const auto& name : result.value().constants) {
            basis = basis->substitute(
                name, SymbolicExpr::number(name == selected ? 1 : 0));
        }
        auto value = evaluate_numeric(*basis, {{"x", 0.0}});
        auto first = evaluate_numeric(
            *basis->differentiate("x"), {{"x", 0.0}});
        auto second = evaluate_numeric(
            *basis->differentiate("x")->differentiate("x"),
            {{"x", 0.0}});
        ASSERT_TRUE(value.has_value());
        ASSERT_TRUE(first.has_value());
        ASSERT_TRUE(second.has_value());
        EXPECT_DOUBLE_EQ(value.value().value, 1.0);
        EXPECT_NEAR(
            second.value().value - 2.0 * first.value().value +
                constant * value.value().value,
            0.0, 1.0e-14);
        exponents.push_back(first.value().value);
    }
    std::sort(exponents.begin(), exponents.end());
    EXPECT_EQ(exponents[0], 1.0 - delta);
    EXPECT_EQ(exponents[1], 1.0 + delta);
}

TEST(OdeEngineHigherOrder, ExactRepeatedQuadraticRoot) {
    auto result = solve_higher_order_ode_checked(
        {1.0, -2.0, 1.0}, nullptr, "x", "y");
    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(result.value().general_solution);
    ASSERT_EQ(result.value().constants.size(), 2u);

    std::vector<double> values;
    for (const auto& selected : result.value().constants) {
        auto basis = result.value().general_solution;
        for (const auto& name : result.value().constants) {
            basis = basis->substitute(
                name, SymbolicExpr::number(name == selected ? 1 : 0));
        }
        auto value = evaluate_numeric(*basis, {{"x", 0.0}});
        auto first = evaluate_numeric(
            *basis->differentiate("x"), {{"x", 0.0}});
        auto second = evaluate_numeric(
            *basis->differentiate("x")->differentiate("x"),
            {{"x", 0.0}});
        ASSERT_TRUE(value.has_value());
        ASSERT_TRUE(first.has_value());
        ASSERT_TRUE(second.has_value());
        EXPECT_DOUBLE_EQ(first.value().value, 1.0);
        EXPECT_NEAR(
            second.value().value - 2.0 * first.value().value +
                value.value().value,
            0.0, 1.0e-14);
        values.push_back(value.value().value);
    }
    std::sort(values.begin(), values.end());
    EXPECT_DOUBLE_EQ(values[0], 0.0);
    EXPECT_DOUBLE_EQ(values[1], 1.0);
}

TEST(OdeEngineHigherOrder, EulerHomogeneous) {
    auto euler = solve_euler_ode_checked(
        {1.0, 1.0, -1.0}, nullptr, "x", "y");
    ASSERT_TRUE(euler.has_value())
        << "checked homogeneous Euler ODE succeeds";
    ASSERT_TRUE(euler.value().general_solution);
    EXPECT_EQ(euler.value().method_used, ODEType::Euler);
    EXPECT_EQ(
        euler.value().constants,
        (std::vector<std::string>{"C1", "C2"}));
}

TEST(OdeEngineHigherOrder, HigherOrderInvalidAndForcing) {
    auto bad_coeffs = solve_higher_order_ode_checked({}, nullptr, "x", "y");
    EXPECT_TRUE((!bad_coeffs.has_value())) << "checked higher-order ODE rejects empty coefficient list";
    EXPECT_TRUE((bad_coeffs.error().code == CasErrc::InvalidArgument)) << "checked higher-order ODE reports InvalidArgument for empty coefficients";

    auto bad_leading = solve_higher_order_ode_checked({0.0, 1.0}, nullptr, "x", "y");
    EXPECT_TRUE((!bad_leading.has_value())) << "checked higher-order ODE rejects zero leading coefficient";
    EXPECT_TRUE((bad_leading.error().code == CasErrc::InvalidArgument)) << "checked higher-order ODE reports InvalidArgument for zero leading coefficient";

    auto same_vars = solve_euler_ode_checked({1.0, 1.0, -1.0}, nullptr, "x", "x");
    EXPECT_TRUE((!same_vars.has_value())) << "checked Euler ODE rejects duplicate variables";
    EXPECT_TRUE((same_vars.error().code == CasErrc::InvalidArgument)) << "checked Euler ODE reports InvalidArgument for duplicate variables";

    const std::vector<double> coefficients{1.0, 0.0, 1.0};
    auto forcing = SymbolicExpr::number(1);
    auto constant_forcing = solve_higher_order_ode_checked(
        coefficients, forcing, "x", "y");
    ASSERT_TRUE(constant_forcing.has_value())
        << "checked higher-order ODE supports constant forcing";
    expect_constant_coefficient_solution(
        constant_forcing.value(), coefficients, forcing);

    auto euler_forcing = solve_euler_ode_checked(
        {1.0, 1.0, -1.0}, SymbolicExpr::number(1), "x", "y");
    EXPECT_TRUE((euler_forcing.has_value())) << "checked Euler ODE supports constant nonhomogeneous forcing";
}

TEST(OdeEngineHigherOrder, LargeRoot) {
    auto large_root = solve_higher_order_ode_checked(
        {1.0, -1.0e20}, nullptr, "x", "y");
    ASSERT_TRUE((large_root.has_value())) << "checked higher-order ODE supports finite roots outside int range";
    if (large_root && large_root.value().general_solution) {
        auto unit_solution =
            large_root.value().general_solution->substitute("C1", SymbolicExpr::number(1))->simplify();
        auto derivative_at_zero = evaluate_numeric(
            *unit_solution->differentiate("x"),
            {{"x", 0.0}});
        ASSERT_TRUE((derivative_at_zero.has_value())) << "large-root solution derivative is numerically evaluable";
        if (derivative_at_zero) {
            EXPECT_NEAR(derivative_at_zero.value().value, 1.0e20, 1.0e6) << "large finite characteristic root is preserved";
        }
    }
}

TEST(OdeEngineHigherOrder, SeparatedRoots) {
    auto separated_roots = solve_higher_order_ode_checked(
        {1.0, 1.0e20, 1.0}, nullptr, "x", "y");
    ASSERT_TRUE((separated_roots.has_value())) << "checked higher-order ODE supports widely separated quadratic roots";
    if (separated_roots && separated_roots.value().general_solution) {
        bool preserved_small_root = false;
        for (const auto &selected : separated_roots.value().constants) {
            auto basis = separated_roots.value().general_solution;
            for (const auto &constant : separated_roots.value().constants) {
                basis = basis->substitute(
                    constant,
                    SymbolicExpr::number(constant == selected ? 1 : 0));
            }
            auto derivative = evaluate_numeric(
                *basis->differentiate("x"), {{"x", 0.0}});
            if (derivative &&
                std::abs((derivative.value().value + 1.0e-20) / 1.0e-20) <
                    1.0e-12) {
                preserved_small_root = true;
            }
        }
        EXPECT_TRUE((preserved_small_root)) << "quadratic solution preserves the small characteristic root";
    }
}

TEST(OdeEngineHigherOrder, UnrepresentableRoot) {
    auto unrepresentable_root = solve_higher_order_ode_checked(
        {1.0e-14, std::numeric_limits<double>::max()},
        nullptr, "x", "y");
    EXPECT_TRUE((!unrepresentable_root.has_value())) << "checked higher-order ODE rejects an unrepresentable root";
    if (!unrepresentable_root) {
        EXPECT_TRUE((unrepresentable_root.error().code == CasErrc::NumericFailure)) << "unrepresentable characteristic roots report NumericFailure";
    }
}

TEST(OdeEngineHigherOrder, UnverifiedRoots) {
    auto unverified_roots = solve_higher_order_ode_checked(
        {1.0, 1.0e20, 0.0, 1.0e-40}, nullptr, "x", "y");
    EXPECT_TRUE((!unverified_roots.has_value())) << "checked higher-order ODE rejects roots with large backward error";
    if (!unverified_roots) {
        EXPECT_TRUE((unverified_roots.error().code == CasErrc::NumericFailure)) << "unverified characteristic roots report NumericFailure";
    }
}

TEST(OdeEngineHigherOrder, SmallComplexRoots) {
    auto small_complex_roots = solve_higher_order_ode_checked(
        {1.0, 1.0, 1.0e-20, 1.0e-20}, nullptr, "x", "y");
    ASSERT_TRUE((small_complex_roots.has_value())) << "checked higher-order ODE supports small nonzero complex roots";
    if (small_complex_roots &&
        small_complex_roots.value().general_solution) {
        bool preserved_cosine_basis = false;
        bool preserved_sine_basis = false;
        for (const auto &selected :
             small_complex_roots.value().constants) {
            auto basis =
                small_complex_roots.value().general_solution;
            for (const auto &constant :
                 small_complex_roots.value().constants) {
                basis = basis->substitute(
                    constant,
                    SymbolicExpr::number(
                        constant == selected ? 1 : 0));
            }
            auto first = evaluate_numeric(
                *basis->differentiate("x"), {{"x", 0.0}});
            auto second = evaluate_numeric(
                *basis->differentiate("x")->differentiate("x"),
                {{"x", 0.0}});
            if (first && second) {
                preserved_sine_basis =
                    preserved_sine_basis ||
                    std::abs(
                        (first.value().value - 1.0e-10) /
                        1.0e-10) < 1.0e-6;
                preserved_cosine_basis =
                    preserved_cosine_basis ||
                    std::abs(
                        (second.value().value + 1.0e-20) /
                        1.0e-20) < 1.0e-6;
            }
        }
        EXPECT_TRUE((preserved_cosine_basis)) << "small complex pair preserves its cosine basis";
        EXPECT_TRUE((preserved_sine_basis)) << "small complex pair preserves its sine basis";
    }
}

TEST(OdeEngineHigherOrder, CloseComplexRoots) {
    auto close_complex_roots = solve_higher_order_ode_checked(
        {1.0, -1.0, -0.9999999999, 1.0000000001},
        nullptr, "x", "y");
    ASSERT_TRUE((close_complex_roots.has_value())) << "checked higher-order ODE supports a close complex pair";
    if (close_complex_roots &&
        close_complex_roots.value().general_solution) {
        bool preserved_sine_basis = false;
        for (const auto &selected :
             close_complex_roots.value().constants) {
            auto basis =
                close_complex_roots.value().general_solution;
            for (const auto &constant :
                 close_complex_roots.value().constants) {
                basis = basis->substitute(
                    constant,
                    SymbolicExpr::number(
                        constant == selected ? 1 : 0));
            }
            auto first = evaluate_numeric(
                *basis->differentiate("x"), {{"x", 0.0}});
            auto second = evaluate_numeric(
                *basis->differentiate("x")->differentiate("x"),
                {{"x", 0.0}});
            if (first && second) {
                preserved_sine_basis =
                    preserved_sine_basis ||
                    (std::abs(
                         (first.value().value - 1.0e-5) /
                         1.0e-5) < 1.0e-5 &&
                     std::abs(
                         (second.value().value - 2.0e-5) /
                         2.0e-5) < 1.0e-5);
            }
        }
        EXPECT_TRUE((preserved_sine_basis)) << "close complex pair is not projected onto the real axis";
    }
}

TEST(OdeEngineHigherOrder, ScaledCoefficients) {
    auto uniformly_small_coefficients =
        solve_higher_order_ode_checked(
            {1.0e-300, 1.0e-300}, nullptr, "x", "y");
    ASSERT_TRUE((uniformly_small_coefficients.has_value())) << "uniform coefficient scaling does not change an ODE";
    if (uniformly_small_coefficients &&
        uniformly_small_coefficients.value().general_solution) {
        auto basis =
            uniformly_small_coefficients.value().general_solution->substitute("C1", SymbolicExpr::number(1))->simplify();
        auto derivative = evaluate_numeric(
            *basis->differentiate("x"), {{"x", 0.0}});
        ASSERT_TRUE((derivative.has_value())) << "uniformly scaled ODE solution is numerically evaluable";
        if (derivative) {
            EXPECT_NEAR(derivative.value().value, -1.0, 1.0e-12) << "uniformly scaled ODE preserves its characteristic root";
        }
    }
}

TEST(OdeEngineHigherOrder, LargeCharacteristicRoots) {
    const double root_scale = 1.0e100;
    auto large_characteristic_roots =
        solve_higher_order_ode_checked(
            {1.0, -6.0e100, 1.1e201, -6.0e300},
            nullptr, "x", "y");
    ASSERT_TRUE((large_characteristic_roots.has_value())) << "finite large characteristic roots are supported";
    if (large_characteristic_roots &&
        large_characteristic_roots.value().general_solution) {
        bool found_one = false;
        bool found_two = false;
        bool found_three = false;
        for (const auto &selected :
             large_characteristic_roots.value().constants) {
            auto basis =
                large_characteristic_roots.value().general_solution;
            for (const auto &constant :
                 large_characteristic_roots.value().constants) {
                basis = basis->substitute(
                    constant,
                    SymbolicExpr::number(
                        constant == selected ? 1 : 0));
            }
            auto derivative = evaluate_numeric(
                *basis->differentiate("x"), {{"x", 0.0}});
            if (derivative) {
                const double scaled_root =
                    derivative.value().value / root_scale;
                found_one |= std::abs(scaled_root - 1.0) < 1.0e-8;
                found_two |= std::abs(scaled_root - 2.0) < 1.0e-8;
                found_three |= std::abs(scaled_root - 3.0) < 1.0e-8;
            }
        }
        EXPECT_TRUE((found_one && found_two && found_three)) << "large characteristic roots preserve their finite scale";
    }
}

TEST(OdeEngineHigherOrder, HigherOrderContext) {
    LMCAS::CancellationToken cancellation;
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    cancellation.cancel();
    auto cancelled = solve_higher_order_ode_checked(
        {1.0, 0.0, 1.0}, nullptr, "x", "y", cancelled_context);
    EXPECT_TRUE((!cancelled.has_value())) << "checked higher-order ODE observes cancellation";
    EXPECT_TRUE((cancelled.error().code == CasErrc::Cancelled)) << "checked higher-order ODE reports Cancelled";

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = solve_euler_ode_checked(
        {1.0, 1.0, -1.0}, nullptr, "x", "y", limited_context);
    EXPECT_TRUE((!limited.has_value())) << "checked Euler ODE observes exhausted step budget";
    EXPECT_TRUE((limited.error().code == CasErrc::ResourceLimit)) << "checked Euler ODE reports ResourceLimit";
}
