#include "test_rootof_support.hpp"
#include "residual_verification.hpp"

TEST(Rootof, RootofPolynomialResidualProof) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("_residual_root");
    auto root = SymbolicExpr::root_of(
        SymbolicExpr::add(SymbolicExpr::power(x, num(3)), num(-2)), "x", 1);
    auto same_root = SymbolicExpr::root_of(
        SymbolicExpr::add(
            SymbolicExpr::multiply(num(7), SymbolicExpr::power(y, num(3))),
            num(-14)),
        "_residual_root", 1);
    auto residual = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::power(root, num(3)),
                               SymbolicExpr::power(same_root, num(3))),
        num(-4));
    ComputationContext context;
    auto checked = check_zero_residual(residual, context);
    EXPECT_TRUE((checked &&
                 std::holds_alternative<ProvedZeroResidual>(checked.value())))
        << "the same exact identity reduces even with scaled and renamed presentations";
}

TEST(Rootof, RootofResidualDoesNotIdentifyDifferentRoots) {
    auto x = SymbolicExpr::variable("x");
    auto polynomial = SymbolicExpr::add(SymbolicExpr::power(x, num(3)), num(-2));
    auto first = SymbolicExpr::root_of(polynomial, "x", 0);
    auto second = SymbolicExpr::root_of(polynomial, "x", 1);
    auto difference = SymbolicExpr::add(
        first, SymbolicExpr::multiply(num(-1), second));
    ComputationContext different_context;
    auto different = check_zero_residual(difference, different_context);
    EXPECT_TRUE((different &&
                 !std::holds_alternative<ProvedZeroResidual>(different.value())))
        << "two different selected roots must not become the same dummy";

    auto dummy = SymbolicExpr::variable("_residual_root");
    auto collision = SymbolicExpr::add(
        SymbolicExpr::power(first, num(3)),
        SymbolicExpr::multiply(num(-1), SymbolicExpr::power(dummy, num(3))));
    ComputationContext collision_context;
    auto free_variable = check_zero_residual(collision, collision_context);
    EXPECT_TRUE((free_variable &&
                 !std::holds_alternative<ProvedZeroResidual>(free_variable.value())))
        << "a free variable with the dummy name must not be replaced by a root";

    auto nonzero = SymbolicExpr::add(SymbolicExpr::power(first, num(3)), num(-1));
    ComputationContext nonzero_context;
    auto remainder = check_zero_residual(nonzero, nonzero_context);
    EXPECT_TRUE((remainder &&
                 !std::holds_alternative<ProvedZeroResidual>(remainder.value())))
        << "a nonzero remainder is never certified as zero";
}

TEST(Rootof, RootofResidualPreservesResourceAndCancellationErrors) {
    auto x = SymbolicExpr::variable("x");
    auto root = SymbolicExpr::root_of(
        SymbolicExpr::add(SymbolicExpr::power(x, num(3)), num(-2)), "x", 0);
    auto residual = SymbolicExpr::add(SymbolicExpr::power(root, num(6)), num(-4));
    ResourceLimits limits;
    limits.max_expansion_terms = 3;
    ComputationContext limited_context(limits);
    auto limited = check_zero_residual(residual, limited_context);
    EXPECT_TRUE((!limited && limited.error().code == CasErrc::ResourceLimit)) << "polynomial residual recognition enforces its expansion budget";

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled_context({}, cancellation);
    auto cancelled = check_zero_residual(residual, cancelled_context);
    EXPECT_TRUE((!cancelled && cancelled.error().code == CasErrc::Cancelled)) << "cancelled exact residual proofs do not become unproved success";
}

static void test_rootof_complex_evaluation(
    const std::shared_ptr<SymbolicExpr> &x) {
    auto complex_poly = SymbolicExpr::add(SymbolicExpr::power(x, num(2)), num(1));
    LMCAS::ComputationContext complex_context;
    auto complex_expression = SymbolicExpr::root_of(complex_poly, "x", 0);
    auto complex_root = LMCAS::rootof_evaluate_checked(
        complex_expression, complex_context);
    EXPECT_TRUE((!complex_root &&
                 complex_root.error().code == LMCAS::CasErrc::DomainError))
        << "real RootOf wrapper rejects a non-real selected root";
    auto complex_value = LMCAS::rootof_evaluate_complex_checked(
        complex_expression);
    EXPECT_TRUE((complex_value &&
                 std::abs(complex_value.value().real.value) < 1e-9 &&
                 std::abs(complex_value.value().imag.value + 1.0) < 1e-8))
        << "complex RootOf evaluates the first ordered root as -i";
    auto conjugate_value = LMCAS::rootof_evaluate_complex_checked(
        SymbolicExpr::root_of(complex_poly, "x", 1));
    EXPECT_TRUE((conjugate_value &&
                 std::abs(conjugate_value.value().imag.value - 1.0) < 1e-8))
        << "complex RootOf evaluates the conjugate root as +i";
}

static void test_rootof_domain_and_resource_errors(
    const std::shared_ptr<SymbolicExpr> &x,
    const std::shared_ptr<SymbolicExpr> &exact_poly) {
    auto parametric_poly = SymbolicExpr::add(
        SymbolicExpr::power(x, num(2)), SymbolicExpr::variable("a"));
    LMCAS::ComputationContext parametric_context;
    auto parametric = LMCAS::make_rootof_checked(
        parametric_poly, "x", 0, parametric_context);
    EXPECT_TRUE((!parametric &&
                 parametric.error().code == LMCAS::CasErrc::Inconclusive))
        << "parametric RootOf construction is explicitly inconclusive";

    auto approximate_poly = SymbolicExpr::add(
        SymbolicExpr::power(x, num(2)), SymbolicExpr::number(-4.0));
    LMCAS::ComputationContext approximate_context;
    auto approximate = LMCAS::make_rootof_checked(
        approximate_poly, "x", 0, approximate_context);
    EXPECT_TRUE((!approximate &&
                 approximate.error().code == LMCAS::CasErrc::Inconclusive))
        << "ApproxReal coefficients are not silently exactified";

    LMCAS::CancellationToken cancellation;
    cancellation.cancel();
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    auto cancelled = LMCAS::rootof_evaluate_checked(
        SymbolicExpr::root_of(exact_poly, "x", 0), cancelled_context);
    EXPECT_TRUE((!cancelled &&
                 cancelled.error().code == LMCAS::CasErrc::Cancelled))
        << "checked RootOf observes cancellation";

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = LMCAS::rootof_evaluate_checked(
        SymbolicExpr::root_of(exact_poly, "x", 0), limited_context);
    EXPECT_TRUE((!limited &&
                 limited.error().code == LMCAS::CasErrc::ResourceLimit))
        << "checked RootOf enforces traversal budgets";
}

TEST(Rootof, RootofCheckedEvaluationPreservesFailureSemantics) {
    auto x = SymbolicExpr::variable("x");
    auto exact_poly = SymbolicExpr::add(SymbolicExpr::power(x, num(2)), num(-4));

    LMCAS::ComputationContext valid_context;
    auto valid = LMCAS::rootof_evaluate_checked(
        SymbolicExpr::root_of(exact_poly, "x", 0), valid_context);
    EXPECT_TRUE((valid && std::abs(valid.value() + 2.0) < 1e-10)) << "checked RootOf verifies an exact real root";

    const BigInt largest_finite_integer =
        (BigInt(1) << 1024) - (BigInt(1) << 971);
    Polynomial<Rational> endpoint_poly(
        {Rational(-largest_finite_integer), Rational(1)}, "x");
    auto endpoint_root = LMCAS::rootof_evaluate_checked(
        SymbolicExpr::root_of(
            LMCAS::poly_to_symbolic(endpoint_poly), "x", 0));
    const std::string endpoint_message =
        endpoint_root
            ? "checked RootOf represents the largest finite double root"
            : "largest finite RootOf failed: " +
                  endpoint_root.error().message;
    EXPECT_TRUE((endpoint_root &&
                 endpoint_root.value() == std::numeric_limits<double>::max()))
        << endpoint_message;

    LMCAS::ComputationContext invalid_index_context;
    auto invalid_index = LMCAS::make_rootof_checked(
        exact_poly, "x", 2, invalid_index_context);
    EXPECT_TRUE((!invalid_index &&
                 invalid_index.error().code == LMCAS::CasErrc::InvalidArgument))
        << "out-of-range RootOf construction is InvalidArgument";

    test_rootof_complex_evaluation(x);

    test_rootof_domain_and_resource_errors(x, exact_poly);
}

TEST(Rootof, RootofOutOfRangeIndicesAreRejected) {
    auto poly_expr = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(SymbolicExpr::variable("x"), num(3)),
            SymbolicExpr::variable("x")),
        num(1));
    for (const std::size_t index : {std::size_t(3), std::size_t(5),
                                    std::size_t(100)}) {
        auto result = LMCAS::make_rootof_checked(
            poly_expr, "x", index);
        EXPECT_TRUE((!result &&
                     result.error().code ==
                         LMCAS::CasErrc::InvalidArgument))
            << "out-of-range RootOf index is rejected";
    }
}

TEST(Rootof, RootofNegativeIndicesAreRejectedAtConstruction) {
    auto poly_expr = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(SymbolicExpr::variable("x"), num(3)),
            SymbolicExpr::multiply(num(-2), SymbolicExpr::variable("x"))),
        num(1));

    bool negative_one_rejected = false;
    try {
        (void)SymbolicExpr::root_of(poly_expr, "x", -1);
    } catch (const std::invalid_argument &) {
        negative_one_rejected = true;
    }
    EXPECT_TRUE((negative_one_rejected)) << "RootOf construction rejects index -1";

    bool negative_ten_rejected = false;
    try {
        (void)SymbolicExpr::root_of(poly_expr, "x", -10);
    } catch (const std::invalid_argument &) {
        negative_ten_rejected = true;
    }
    EXPECT_TRUE((negative_ten_rejected)) << "RootOf construction rejects index -10";
}

TEST(Rootof, RootofCanonicalIdentityIgnoresScaleRepetitionAndDummyName) {
    auto x = SymbolicExpr::variable("x");
    auto base = SymbolicExpr::add(
        SymbolicExpr::power(x, num(2)), num(-2));
    auto scaled = SymbolicExpr::add(
        SymbolicExpr::multiply(
            num(2), SymbolicExpr::power(x, num(2))),
        num(-4));
    auto y = SymbolicExpr::variable("y");
    auto renamed = SymbolicExpr::add(
        SymbolicExpr::power(y, num(2)), num(-2));
    auto repeated = SymbolicExpr::power(base, num(2));

    auto canonical = SymbolicExpr::root_of(base, "x", 0);
    auto scaled_root = SymbolicExpr::root_of(scaled, "x", 0);
    auto renamed_root = SymbolicExpr::root_of(renamed, "y", 0);
    auto repeated_root = SymbolicExpr::root_of(repeated, "x", 0);
    EXPECT_TRUE((LMCAS::detail::node(canonical)->equals(
                     *LMCAS::detail::node(scaled_root)) &&
                 LMCAS::detail::node(canonical)->equals(
                     *LMCAS::detail::node(renamed_root)) &&
                 LMCAS::detail::node(canonical)->equals(
                     *LMCAS::detail::node(repeated_root))))
        << "canonical RootOf identities are structurally equal";
    EXPECT_TRUE((LMCAS::detail::node(canonical)->hash() ==
                 LMCAS::detail::node(repeated_root)->hash()))
        << "canonical RootOf identities share a hash";
}

TEST(Rootof, RootofParametricCoefficientsAreInconclusive) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto poly_expr = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, num(3)),
            SymbolicExpr::multiply(a, x)),
        num(1));
    auto result = LMCAS::make_rootof_checked(
        poly_expr, "x", 0);
    EXPECT_TRUE((!result &&
                 result.error().code == LMCAS::CasErrc::Inconclusive))
        << "parametric RootOf construction is Inconclusive";
}

TEST(Rootof, RootofMultipleParametricCoefficientsAreInconclusive) {
    auto x = SymbolicExpr::variable("x");
    auto b = SymbolicExpr::variable("b");
    auto c = SymbolicExpr::variable("c");
    auto poly_expr = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, num(2)),
            SymbolicExpr::multiply(b, x)),
        c);
    auto result = LMCAS::make_rootof_checked(
        poly_expr, "x", 0);
    EXPECT_TRUE((!result &&
                 result.error().code == LMCAS::CasErrc::Inconclusive))
        << "multiple parametric coefficients are Inconclusive";
}

TEST(Rootof, RootofSimplifyDegree2PolynomialToClosedForm) {
    auto x = SymbolicExpr::variable("x");
    auto poly_expr = SymbolicExpr::add(
        SymbolicExpr::power(x, num(2)),
        num(-4));

    auto rootof_k0 = SymbolicExpr::root_of(poly_expr, "x", 0);
    auto simplified_k0 = LMCAS::rootof_simplify(rootof_k0);

    std::string s0 = simplified_k0->to_string();
    EXPECT_TRUE((s0.find("RootOf") == std::string::npos)) << "rootof_simplify(degree-2, k=0) returns non-RootOf expression";

    double val0 = eval_numeric(simplified_k0);
    EXPECT_TRUE((!std::isnan(val0) && std::abs(val0 - (-2.0)) < 1e-10)) << "rootof_simplify(x^2-4, x, 0) = -2";

    auto rootof_k1 = SymbolicExpr::root_of(poly_expr, "x", 1);
    auto simplified_k1 = LMCAS::rootof_simplify(rootof_k1);

    std::string s1 = simplified_k1->to_string();
    EXPECT_TRUE((s1.find("RootOf") == std::string::npos)) << "rootof_simplify(degree-2, k=1) returns non-RootOf expression";

    double val1 = eval_numeric(simplified_k1);
    EXPECT_TRUE((!std::isnan(val1) && std::abs(val1 - 2.0) < 1e-10)) << "rootof_simplify(x^2-4, x, 1) = 2";
}
