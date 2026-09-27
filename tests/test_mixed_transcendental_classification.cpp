#include "test_mixed_transcendental_support.hpp"

TEST(MixedTranscendentalClassification, ClassificationSinXPlusX) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    EXPECT_TRUE((contains_transcendental_of_var(expr, "x"))) << "sin(x) + x should contain transcendental of var x";
}

TEST(MixedTranscendentalClassification, ClassificationExpXMinusXSquared) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
                               SymbolicExpr::power(x, SymbolicExpr::number(2))));

    EXPECT_TRUE((contains_transcendental_of_var(expr, "x"))) << "exp(x) - x^2 should contain transcendental of var x";
}

TEST(MixedTranscendentalClassification, ClassificationSinConstantPlusPolynomial) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::sin(SymbolicExpr::number(3)));

    EXPECT_FALSE((contains_transcendental_of_var(expr, "x"))) << "x^2 + sin(3) should NOT contain transcendental of var x (sin arg is constant)";
}

TEST(MixedTranscendentalClassification, ClassificationCosYPlusX) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto expr = SymbolicExpr::add(SymbolicExpr::cos(y), x);

    EXPECT_FALSE((contains_transcendental_of_var(expr, "x"))) << "cos(y) + x should NOT contain transcendental of var x (cos arg depends on y, not x)";
}

TEST(MixedTranscendentalClassification, ClassificationLnCosXPlusX) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::ln(x), SymbolicExpr::cos(x)),
        x);

    EXPECT_TRUE((contains_transcendental_of_var(expr, "x"))) << "ln(x) * cos(x) + x should contain transcendental of var x";
}

TEST(MixedTranscendentalClassification, ClassificationTanLinearArg) {
    auto x = SymbolicExpr::variable("x");
    auto tan_arg = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::number(1));
    auto expr = SymbolicExpr::add(
        SymbolicExpr::tan(tan_arg),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), x));

    EXPECT_TRUE((contains_transcendental_of_var(expr, "x"))) << "tan(2*x + 1) - x should contain transcendental of var x";
}

TEST(MixedTranscendentalClassification, ClassificationPurePolynomial) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(3)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-2), x),
            SymbolicExpr::number(1)));

    EXPECT_FALSE((contains_transcendental_of_var(expr, "x"))) << "x^3 - 2*x + 1 should NOT contain transcendental of var x (pure polynomial)";
}

TEST(MixedTranscendentalClassification, AllowNumericFalseIsInconclusive) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(x, SymbolicExpr::sin(x)),
        SymbolicExpr::number(-1));

    LMCAS::SolveOptions opts;
    opts.allow_numeric = false;

    auto result = solve_equation(expr, "x", opts);

    EXPECT_TRUE((!result && result.error().code == CasErrc::Inconclusive)) << "solve_equation with allow_numeric=false on x*sin(x)-1 should be Inconclusive";
}

TEST(MixedTranscendentalClassification, AllowNumericFalseCosEquationIsInconclusive) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(x, SymbolicExpr::cos(x)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(
                SymbolicExpr::power(x, SymbolicExpr::number(2)),
                SymbolicExpr::sin(x)),
            SymbolicExpr::number(-1)));

    LMCAS::SolveOptions opts;
    opts.allow_numeric = false;

    auto result = solve_equation(expr, "x", opts);

    EXPECT_TRUE((!result && result.error().code == CasErrc::Inconclusive)) << "solve_equation with allow_numeric=false on x*cos(x)+x^2*sin(x)-1 should be Inconclusive";
}

TEST(MixedTranscendentalClassification, RoutingPurePolynomialNotHybrid) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-4));

    LMCAS::SolveOptions opts;
    opts.allow_numeric = true;

    auto results = solve_vector_for_test(expr, "x", opts);
    ASSERT_EQ(results.size(), 2U);

    std::vector<double> roots;
    for (const auto &root : results) {
        auto value = test_numeric_eval(root);
        ASSERT_TRUE(value.has_value());
        roots.push_back(*value);
    }
    std::sort(roots.begin(), roots.end());
    EXPECT_NEAR(roots[0], -2.0, 1e-10);
    EXPECT_NEAR(roots[1], 2.0, 1e-10);
}

TEST(MixedTranscendentalClassification, RoutingTranscendentalSubstitutionNotHybrid) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::number(-2));

    for (bool allow_numeric : {false, true}) {
        SCOPED_TRACE(
            allow_numeric ? "numeric enabled" : "numeric disabled");
        LMCAS::SolveOptions opts;
        opts.allow_numeric = allow_numeric;

        auto results = solve_vector_for_test(expr, "x", opts);

        ASSERT_EQ(results.size(), 1U);
        auto value = test_numeric_eval(results.front());
        ASSERT_TRUE(value.has_value());
        EXPECT_NEAR(*value, std::log(2.0), 1e-10);
        auto residual = expr->substitute("x", results.front())->simplify();
        ASSERT_NE(residual, nullptr);
        auto residual_value = test_numeric_eval(residual);
        ASSERT_TRUE(residual_value.has_value());
        EXPECT_NEAR(*residual_value, 0.0, 1e-10);
    }
}

TEST(MixedTranscendentalClassification, RoutingPolynomialWithAllowNumericFalse) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-4));

    LMCAS::SolveOptions opts;
    opts.allow_numeric = false;

    auto results = solve_vector_for_test(expr, "x", opts);
    ASSERT_EQ(results.size(), 2U);

    std::vector<double> roots;
    for (const auto &root : results) {
        auto value = test_numeric_eval(root);
        ASSERT_TRUE(value.has_value());
        roots.push_back(*value);
    }
    std::sort(roots.begin(), roots.end());
    EXPECT_NEAR(roots[0], -2.0, 1e-10);
    EXPECT_NEAR(roots[1], 2.0, 1e-10);
}

TEST(MixedTranscendentalClassification, CubicPolynomialNumericEnabled) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(3)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-6),
                                   SymbolicExpr::power(x, SymbolicExpr::number(2))),
            SymbolicExpr::add(
                SymbolicExpr::multiply(SymbolicExpr::number(11), x),
                SymbolicExpr::number(-6))));
    LMCAS::SolveOptions opts_numeric;
    opts_numeric.allow_numeric = true;
    auto results_numeric = solve_vector_for_test(expr, "x", opts_numeric);

    ASSERT_EQ(results_numeric.size(), 3u)
        << "x^3-6x^2+11x-6=0 should produce 3 roots with allow_numeric=true";
    std::vector<double> root_vals;
    for (const auto &root : results_numeric) {
        auto value = test_numeric_eval(root);
        ASSERT_TRUE(value.has_value()) << "Root should be evaluable to a number";
        root_vals.push_back(*value);
    }
    std::sort(root_vals.begin(), root_vals.end());
    for (size_t i = 0; i < root_vals.size(); ++i) {
        const double expected = static_cast<double>(i + 1);
        SCOPED_TRACE(i);
        EXPECT_TRUE(std::isfinite(root_vals[i]));
        EXPECT_TRUE(std::isfinite(expected));
        EXPECT_NEAR(root_vals[i], expected, 1e-10);
    }
}

TEST(MixedTranscendentalClassification, CubicPolynomialNumericDisabled) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(3)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-6),
                                   SymbolicExpr::power(x, SymbolicExpr::number(2))),
            SymbolicExpr::add(
                SymbolicExpr::multiply(SymbolicExpr::number(11), x),
                SymbolicExpr::number(-6))));
    LMCAS::SolveOptions opts_symbolic;
    opts_symbolic.allow_numeric = false;
    auto results_symbolic = solve_vector_for_test(
        expr, "x", opts_symbolic);
    ASSERT_EQ(results_symbolic.size(), 3U);

    std::vector<double> root_vals;
    for (const auto &root : results_symbolic) {
        auto value = test_numeric_eval(root);
        ASSERT_TRUE(value.has_value());
        root_vals.push_back(*value);
    }
    std::sort(root_vals.begin(), root_vals.end());
    EXPECT_NEAR(root_vals[0], 1.0, 1e-10);
    EXPECT_NEAR(root_vals[1], 2.0, 1e-10);
    EXPECT_NEAR(root_vals[2], 3.0, 1e-10);
}
