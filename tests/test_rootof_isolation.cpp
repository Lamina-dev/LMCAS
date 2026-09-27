#include "test_rootof_support.hpp"

TEST(RootofIsolation, RootofNonPolynomialConstructionIsInconclusive) {
    auto x = SymbolicExpr::variable("x");
    auto unsupported = SymbolicExpr::add(
        SymbolicExpr::power(x, num(2)), SymbolicExpr::sin(x));
    auto root = LMCAS::make_rootof_checked(
        unsupported, "x", 0);
    EXPECT_TRUE((!root &&
                 root.error().code == LMCAS::CasErrc::Inconclusive))
        << "non-polynomial RootOf input is rejected";
}

TEST(RootofIsolation, RootofHigherDegreeSimplificationPreservesCanonicalIdentity) {
    auto x = SymbolicExpr::variable("x");
    auto poly_expr = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(x, num(3)),
                SymbolicExpr::multiply(
                    num(-6), SymbolicExpr::power(x, num(2)))),
            SymbolicExpr::multiply(num(11), x)),
        num(-6));
    for (int index = 0; index < 3; ++index) {
        auto root = SymbolicExpr::root_of(poly_expr, "x", index);
        auto simplified = LMCAS::rootof_simplify(root);
        EXPECT_TRUE((LMCAS::detail::node(simplified)->equals(*LMCAS::detail::node(root)))) << "cubic RootOf simplification preserves identity";
        auto value = LMCAS::rootof_evaluate_checked(root);
        EXPECT_TRUE((value &&
                     std::abs(value.value() -
                              static_cast<double>(index + 1)) < 1e-8))
            << "cubic RootOf exact order is stable";
    }
}

TEST(RootofIsolation, RootofDegreeFourIdentitiesRetainExactOrder) {
    auto x = SymbolicExpr::variable("x");
    auto poly_expr = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, num(4)),
            SymbolicExpr::multiply(
                num(-5), SymbolicExpr::power(x, num(2)))),
        num(4));
    const double expected[] = {-2.0, -1.0, 1.0, 2.0};
    for (int index = 0; index < 4; ++index) {
        auto root = SymbolicExpr::root_of(poly_expr, "x", index);
        auto simplified = LMCAS::rootof_simplify(root);
        EXPECT_TRUE((LMCAS::detail::node(simplified)->equals(*LMCAS::detail::node(root)))) << "quartic RootOf simplification preserves identity";
        auto value = LMCAS::rootof_evaluate_checked(root);
        EXPECT_TRUE((value &&
                     std::abs(value.value() - expected[index]) < 1e-8))
            << "quartic RootOf exact order is stable";
    }
}

TEST(RootofIsolation, RootofValidIndexOnNumericPolynomialEvaluatesCorrectly) {
    auto x = SymbolicExpr::variable("x");
    auto poly_expr = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::add(
                SymbolicExpr::power(x, num(3)),
                SymbolicExpr::multiply(num(-6), SymbolicExpr::power(x, num(2)))),
            SymbolicExpr::multiply(num(11), x)),
        num(-6));

    auto rootof_k0 = SymbolicExpr::root_of(poly_expr, "x", 0);
    auto result0 = LMCAS::rootof_evaluate_checked(rootof_k0);
    ASSERT_TRUE((result0.has_value())) << "checked RootOf evaluation with valid k=0 returns a value";
    if (result0.has_value()) {
        EXPECT_TRUE((std::abs(result0.value() - 1.0) < 1e-10)) << "checked RootOf evaluation of cubic k=0 is the smallest root";
    }

    auto rootof_k2 = SymbolicExpr::root_of(poly_expr, "x", 2);
    auto result2 = LMCAS::rootof_evaluate_checked(rootof_k2);
    ASSERT_TRUE((result2.has_value())) << "checked RootOf evaluation with valid k=2 returns a value";
    if (result2.has_value()) {
        EXPECT_TRUE((std::abs(result2.value() - 3.0) < 1e-10)) << "checked RootOf evaluation of cubic k=2 is the largest root";
    }
}

TEST(RootofIsolation, RootofCertifiedHigherDegreeComplexOrdering) {
    Polynomial<Rational> cubic(
        {Rational(-1), Rational(0), Rational(0), Rational(1)}, "x");
    auto cubic_expr = LMCAS::poly_to_symbolic(cubic);
    const double cubic_expected[][2] = {
        {1.0, 0.0},
        {-0.5, -std::sqrt(3.0) / 2.0},
        {-0.5, std::sqrt(3.0) / 2.0}};
    for (int index = 0; index < 3; ++index) {
        auto root = SymbolicExpr::root_of(cubic_expr, "x", index);
        auto value = LMCAS::rootof_evaluate_complex_checked(root);
        EXPECT_TRUE((value &&
                     std::abs(value.value().real.value -
                              cubic_expected[index][0]) < 1e-9 &&
                     std::abs(value.value().imag.value -
                              cubic_expected[index][1]) < 1e-9))
            << "x^3-1 RootOf index has certified real-first complex order";
        if (index != 0) {
            LMCAS::ComputationContext real_context;
            auto real = LMCAS::rootof_evaluate_checked(
                root, real_context);
            EXPECT_TRUE((!real &&
                         real.error().code == LMCAS::CasErrc::DomainError))
                << "real evaluation rejects a selected cubic non-real root";
        }
    }

    Polynomial<Rational> shifted(
        {Rational(10), Rational(-6), Rational(3),
         Rational(0), Rational(1)},
        "x");
    auto shifted_expr = LMCAS::poly_to_symbolic(shifted);
    const double shifted_expected[][2] = {
        {-1.0, -2.0}, {-1.0, 2.0}, {1.0, -1.0}, {1.0, 1.0}};
    for (int index = 0; index < 4; ++index) {
        auto value = LMCAS::rootof_evaluate_complex_checked(
            SymbolicExpr::root_of(shifted_expr, "x", index));
        EXPECT_TRUE((value &&
                     std::abs(value.value().real.value -
                              shifted_expected[index][0]) < 1e-9 &&
                     std::abs(value.value().imag.value -
                              shifted_expected[index][1]) < 1e-9))
            << "shifted quartic roots use exact real/imaginary lexicographic order";
    }
}

TEST(RootofIsolation, RootofClusteredComplexRectanglesRemainDistinct) {
    const Rational small(1, 1048576);
    Polynomial<Rational> polynomial =
        Polynomial<Rational>({small, Rational(0), Rational(1)}, "x") *
        Polynomial<Rational>({Rational(4) * small,
                              Rational(0), Rational(1)},
                             "x");
    auto expression = LMCAS::poly_to_symbolic(polynomial);
    const double expected[] = {
        -1.0 / 512.0, -1.0 / 1024.0,
        1.0 / 1024.0, 1.0 / 512.0};
    for (int index = 0; index < 4; ++index) {
        auto value = LMCAS::rootof_evaluate_complex_checked(
            SymbolicExpr::root_of(expression, "x", index));
        EXPECT_TRUE((value &&
                     std::abs(value.value().real.value) < 1e-12 &&
                     std::abs(value.value().imag.value - expected[index]) <
                         1e-10 &&
                     value.value().imag.absolute_error > 0.0))
            << "clustered imaginary roots retain distinct certified enclosures";
    }
}

static void test_quintic_identity_and_resource_limits(
    const std::shared_ptr<SymbolicExpr> &expression) {
    auto selected = SymbolicExpr::root_of(expression, "x", 1);
    LMCAS::CancellationToken cancellation;
    cancellation.cancel();
    LMCAS::ComputationContext cancelled({}, cancellation);
    auto cancelled_value =
        LMCAS::rootof_evaluate_complex_checked(selected, cancelled);
    EXPECT_TRUE((!cancelled_value &&
                 cancelled_value.error().code == LMCAS::CasErrc::Cancelled))
        << "higher-degree complex isolation propagates cancellation";

    LMCAS::ResourceLimits limits;
    limits.max_steps = 8;
    LMCAS::ComputationContext limited(limits);
    auto limited_value =
        LMCAS::rootof_evaluate_complex_checked(selected, limited);
    EXPECT_TRUE((!limited_value &&
                 limited_value.error().code == LMCAS::CasErrc::ResourceLimit))
        << "higher-degree complex isolation enforces step budgets";

    auto repeated_expression = SymbolicExpr::power(expression, num(2));
    for (int index = 0; index < 5; ++index) {
        auto canonical = SymbolicExpr::root_of(expression, "x", index);
        auto repeated = SymbolicExpr::root_of(
            repeated_expression, "x", index);
        EXPECT_TRUE((LMCAS::detail::node(canonical)->equals(
            *LMCAS::detail::node(repeated))))
            << "square-free canonicalization preserves quintic complex identity";
    }
}

static bool is_finite_certified_complex_enclosure(const LMCAS::ApproxComplex &value) {
    return std::isfinite(value.real.value) &&
           std::isfinite(value.imag.value) &&
           value.real.absolute_error >= 0.0 &&
           value.imag.absolute_error >= 0.0;
}

TEST(RootofIsolation, RootofGenericQuinticComplexIsolationAndResources) {
    Polynomial<Rational> polynomial(
        {Rational(1), Rational(1), Rational(0),
         Rational(0), Rational(0), Rational(1)},
        "x");
    auto expression = LMCAS::poly_to_symbolic(polynomial);
    std::vector<LMCAS::ApproxComplex> values;
    for (int index = 0; index < 5; ++index) {
        auto value = LMCAS::rootof_evaluate_complex_checked(
            SymbolicExpr::root_of(expression, "x", index));
        EXPECT_TRUE((value && is_finite_certified_complex_enclosure(value.value()))) << "generic quintic RootOf has a finite certified complex enclosure";
        if (value) {
            values.push_back(value.value());
        }
    }
    EXPECT_TRUE((values.size() == 5 &&
                 std::abs(values[0].imag.value) < 1e-12 &&
                 std::abs(values[1].real.value - values[2].real.value) <
                     1e-10 &&
                 std::abs(values[1].imag.value + values[2].imag.value) <
                     1e-10 &&
                 std::abs(values[3].real.value - values[4].real.value) <
                     1e-10 &&
                 std::abs(values[3].imag.value + values[4].imag.value) <
                     1e-10))
        << "generic quintic roots are one real root plus certified conjugate pairs";

    test_quintic_identity_and_resource_limits(expression);
}
