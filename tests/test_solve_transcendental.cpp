#include "test_common.hpp"
#include "test_expression_equivalence.hpp"
#include "numeric_evaluation.hpp"
#include "solve_transcendental.hpp"
#include "expr.hpp"
#include "assumption_context.hpp"
#include "polynomial_conversion.hpp"

using namespace LMCAS;

namespace {

void expect_periodic_family_values(const ParametricSolution &family,
                                   const std::vector<Rational> &offsets,
                                   const Rational &period, std::vector<bool> &matched) {
    auto pi_constant = LMCAS::pi();
    ASSERT_TRUE(pi_constant);
    bool found = false;
    for (std::size_t index = 0; index < offsets.size(); ++index) {
        if (matched[index]) {
            continue;
        }
        bool equal = true;
        for (int k : {-2, 0, 1}) {
            auto actual = family.value->substitute(
                family.integer_parameters.front(), SymbolicExpr::number(k));
            auto expected = SymbolicExpr::multiply(pi_constant.value(),
                                                   SymbolicExpr::add(SymbolicExpr::number(offsets[index]),
                                                                     SymbolicExpr::multiply(SymbolicExpr::number(period), SymbolicExpr::number(k))));
            equal = equal && test_expr_equivalent(actual, expected);
        }
        if (equal) {
            found = true;
            matched[index] = true;
            break;
        }
    }
    EXPECT_TRUE(found) << "family values and integer spacing are exact symbolic multiples of pi";
}

void expect_periodic(const SolveResult &solved, const std::vector<Rational> &offsets,
                     const Rational &period) {
    if (!solved) {
        std::cerr << "CasError [" << solved.error().operation << "]: "
                  << solved.error().message << '\n';
    }
    EXPECT_TRUE(solved && std::holds_alternative<ParametricSolutions>(solved.value())) << "periodic equations retain parametric solutions";
    if (!solved || !std::holds_alternative<ParametricSolutions>(solved.value())) {
        return;
    }
    const auto &families = std::get<ParametricSolutions>(solved.value()).values;
    EXPECT_TRUE(families.size() == offsets.size()) << "families are complete and endpoint duplicates are compressed";
    if (families.size() != offsets.size()) {
        for (const auto &family : families) {
            std::cerr << "family: " << family.value->to_string() << '\n';
            for (const auto &condition : family.conditions) {
                std::cerr << "condition: " << condition->to_string() << '\n';
            }
        }
    }
    if (families.size() != offsets.size()) {
        return;
    }
    std::vector<bool> matched(offsets.size(), false);
    for (const auto &family : families) {
        EXPECT_TRUE(family.integer_parameters.size() == 1) << "each periodic branch has one integer parameter";
        if (family.integer_parameters.size() != 1) {
            continue;
        }
        expect_periodic_family_values(family, offsets, period, matched);
    }
}

void expect_preimage(const ExprPtr &expression) {
    SCOPED_TRACE(expression->to_string());
    auto solved = solve_transcendental(expression, "x");
    ASSERT_TRUE(solved) << solved.error().message;
    EXPECT_TRUE(std::holds_alternative<ConditionalSolutions>(solved.value()))
        << "unresolved inverse retains a complete conditional preimage";
    if (solved && std::holds_alternative<ConditionalSolutions>(solved.value())) {
        const auto &set = std::get<ConditionalSolutions>(solved.value()).value;
        EXPECT_TRUE(set.variable == "x") << "preimage binds the requested variable";
        EXPECT_TRUE(set.predicate->compare(
                        SymbolicExpr::eq(expression, SymbolicExpr::number(0))) == 0)
            << "preimage retains the original equation without dividing";
    }
}

void expect_symbolic_range_conditions(const ParametricSolution &family,
                                      const ExprResult &lower, const ExprResult &upper) {
    bool has_lower = false;
    bool has_upper = false;
    for (const auto &condition : family.conditions) {
        has_lower = has_lower || (lower && condition->compare(lower.value()) == 0);
        has_upper = has_upper || (upper && condition->compare(upper.value()) == 0);
    }
    EXPECT_TRUE(has_lower && has_upper) << "every branch requires exact -1 <= a <= 1";
}

void test_direct_periodic_family(const ExprPtr &x) {

    auto expr = SymbolicExpr::sin(x);

    auto roots = LMCAS::solve_transcendental(expr, "x");
    expect_periodic(roots, {Rational(0)}, Rational(1));
}

void test_affine_periodic_families(const ExprPtr &twice_x) {
    expect_periodic(solve_transcendental(SymbolicExpr::sin(twice_x), "x"),
                    {Rational(0)}, Rational(1, 2));
    expect_periodic(solve_transcendental(SymbolicExpr::cos(twice_x), "x"),
                    {Rational(1, 4)}, Rational(1, 2));
    expect_periodic(solve_transcendental(SymbolicExpr::tan(twice_x), "x"),
                    {Rational(0)}, Rational(1, 2));
    auto half_equation = SymbolicExpr::add(SymbolicExpr::sin(twice_x),
                                           SymbolicExpr::number(Rational(-1, 2)));
    expect_periodic(solve_transcendental(half_equation, "x"),
                    {Rational(1, 12), Rational(5, 12)}, Rational(1));
    expect_periodic(solve_equation(half_equation, "x"),
                    {Rational(1, 12), Rational(5, 12)}, Rational(1));
}

void test_periodic_endpoints_and_unions(const ExprPtr &x, const ExprPtr &twice_x) {
    expect_periodic(solve_transcendental(SymbolicExpr::add(SymbolicExpr::sin(x),
                                                           SymbolicExpr::number(-1)),
                                         "x"),
                    {Rational(1, 2)}, Rational(2));
    expect_periodic(solve_transcendental(SymbolicExpr::add(SymbolicExpr::cos(x),
                                                           SymbolicExpr::number(1)),
                                         "x"),
                    {Rational(1)}, Rational(2));
    auto squared_sine = SymbolicExpr::add(
        SymbolicExpr::power(SymbolicExpr::sin(x), SymbolicExpr::number(2)), SymbolicExpr::number(-1));
    expect_periodic(solve_transcendental(squared_sine, "x"),
                    {Rational(1, 2)}, Rational(1));
    auto squared_affine_sine = SymbolicExpr::add(SymbolicExpr::power(
                                                     SymbolicExpr::sin(twice_x), SymbolicExpr::number(2)),
                                                 SymbolicExpr::number(-1));
    expect_periodic(solve_transcendental(squared_affine_sine, "x"),
                    {Rational(1, 4)}, Rational(1, 2));
    auto out_of_range = solve_transcendental(SymbolicExpr::add(
                                                 SymbolicExpr::sin(x), SymbolicExpr::number(-2)),
                                             "x");
    EXPECT_TRUE(out_of_range && std::holds_alternative<EmptySolutions>(out_of_range.value())) << "sine cannot attain a real target outside [-1,1]";
}

void test_symbolic_range_constraints(const ExprPtr &x, const ExprPtr &a) {
    auto symbolic_target = SymbolicExpr::add(SymbolicExpr::sin(x),
                                             SymbolicExpr::multiply(SymbolicExpr::number(-1), a));
    auto symbolic_roots = solve_transcendental(symbolic_target, "x");
    EXPECT_TRUE(symbolic_roots && std::holds_alternative<ParametricSolutions>(symbolic_roots.value())) << "symbolic sine targets retain both guarded families";
    if (symbolic_roots && std::holds_alternative<ParametricSolutions>(symbolic_roots.value())) {
        const auto &families = std::get<ParametricSolutions>(symbolic_roots.value()).values;
        EXPECT_TRUE(families.size() == 2) << "generic sine inversion has both branches";
        auto lower = LMCAS::ge(a, SymbolicExpr::number(-1));
        auto upper = LMCAS::le(a, SymbolicExpr::number(1));
        EXPECT_TRUE(lower && upper) << "exact target range predicates are constructible";
        for (const auto &family : families) {
            expect_symbolic_range_conditions(family, lower, upper);
        }
    }
}

void test_capture_free_integer_parameter(const ExprPtr &twice_x) {
    auto occupied = SymbolicExpr::variable("_k");
    auto captured_equation = SymbolicExpr::sin(SymbolicExpr::add(twice_x, occupied));
    auto capture_assumptions = std::make_shared<AssumptionContext>();
    EXPECT_TRUE(capture_assumptions->assume_domain("_k", Domain::Real).has_value()) << "the free affine offset is real";
    ComputationContext capture_context;
    EXPECT_TRUE(capture_context.set_assumptions(capture_assumptions).has_value()) << "offset facts attach to the computation";
    auto capture_free = solve_transcendental(captured_equation, "x", capture_context);
    ASSERT_TRUE(capture_free.has_value()) << "occupied parameter names do not prevent solving";
    ASSERT_TRUE(std::holds_alternative<ParametricSolutions>(capture_free.value()));
    const auto &families = std::get<ParametricSolutions>(capture_free.value()).values;
    ASSERT_EQ(families.size(), 1u) << "zero sine target needs one family";
    const auto &family = families.front();
    ASSERT_EQ(family.integer_parameters.size(), 1u);
    EXPECT_TRUE(family.integer_parameters.front() != "_k") << "input _k remains free";
    auto value = family.value->substitute(family.integer_parameters.front(), SymbolicExpr::number(0))
                     ->substitute("_k", SymbolicExpr::number(3));
    EXPECT_TRUE(test_expr_equivalent(value, SymbolicExpr::number(Rational(-3, 2)))) << "affine offset is not captured by the integer binder";
}

void test_unresolved_preimages(const ExprPtr &x, const ExprPtr &a) {
    expect_preimage(SymbolicExpr::sin(SymbolicExpr::power(x, SymbolicExpr::number(2))));
    expect_preimage(SymbolicExpr::sin(SymbolicExpr::multiply(a, x)));
    expect_preimage(SymbolicExpr::add(SymbolicExpr::multiply(a, SymbolicExpr::sin(x)),
                                      SymbolicExpr::number(-1)));
    expect_preimage(SymbolicExpr::sin(SymbolicExpr::add(x, a)));
    auto pi_constant = LMCAS::pi();
    auto imaginary = LMCAS::imaginary_unit();
    ASSERT_TRUE(pi_constant);
    ASSERT_TRUE(imaginary);
    auto complex_offset = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(Rational(1, 2)), pi_constant.value()),
        SymbolicExpr::multiply(imaginary.value(), x));
    auto nonreal_expression = SymbolicExpr::add(
        SymbolicExpr::sin(complex_offset), SymbolicExpr::number(-2));
    auto nonreal_preimage = solve_transcendental(nonreal_expression, "x");
    EXPECT_TRUE((!nonreal_preimage &&
                 nonreal_preimage.error().code == CasErrc::Inconclusive))
        << "real-valuedness of a sine with a complex argument needs proof";
    auto complex_nonzero = std::make_shared<AssumptionContext>();
    EXPECT_TRUE(complex_nonzero->assume_domain("a", Domain::Complex).has_value()) << "complex coefficients are permitted inputs";
    EXPECT_TRUE(complex_nonzero->assume_sign("a", Sign::NonZero).has_value()) << "nonzero alone need not mean real";
    ComputationContext complex_context;
    EXPECT_TRUE(complex_context.set_assumptions(complex_nonzero).has_value()) << "complex nonzero facts attach to the computation";
    auto complex_slope = solve_transcendental(
        SymbolicExpr::sin(SymbolicExpr::multiply(a, x)), "x", complex_context);
    EXPECT_TRUE(complex_slope && std::holds_alternative<ConditionalSolutions>(complex_slope.value())) << "a nonzero complex slope does not justify real inverse families";
}

void test_assumption_proved_slope(const ExprPtr &x, const ExprPtr &a) {
    auto nonzero_assumptions = std::make_shared<AssumptionContext>();
    EXPECT_TRUE(nonzero_assumptions->assume_sign("a", Sign::Positive).has_value()) << "positive slope assumption is accepted";
    ComputationContext nonzero_context;
    EXPECT_TRUE(nonzero_context.set_assumptions(nonzero_assumptions).has_value()) << "slope facts attach to the computation";
    auto proved_slope = solve_transcendental(SymbolicExpr::sin(SymbolicExpr::multiply(a, x)),
                                             "x", nonzero_context);
    ASSERT_TRUE(proved_slope.has_value()) << "proved nonzero affine slope permits family division";
    ASSERT_TRUE(std::holds_alternative<ParametricSolutions>(proved_slope.value()));
    const auto &families = std::get<ParametricSolutions>(proved_slope.value()).values;
    ASSERT_EQ(families.size(), 1u) << "proved slope retains the zero-target family";
    const auto &family = families.front();
    ASSERT_EQ(family.integer_parameters.size(), 1u);
    auto point = family.value->substitute("a", SymbolicExpr::number(2))->substitute(family.integer_parameters.front(), SymbolicExpr::number(1));
    auto pi_constant = LMCAS::pi();
    ASSERT_TRUE(pi_constant);
    EXPECT_TRUE(test_expr_equivalent(point, SymbolicExpr::multiply(
                                                SymbolicExpr::number(Rational(1, 2)), pi_constant.value())))
        << "nonzero slope rescales the period";
    auto zero_assumptions = std::make_shared<AssumptionContext>();
    EXPECT_TRUE(zero_assumptions->assume_sign("a", Sign::Zero).has_value()) << "zero slope assumption is accepted";
    ComputationContext zero_context;
    EXPECT_TRUE(zero_context.set_assumptions(zero_assumptions).has_value()) << "zero slope facts attach to the computation";
    auto constant_sine = solve_transcendental(SymbolicExpr::sin(SymbolicExpr::multiply(a, x)),
                                              "x", zero_context);
    if (!constant_sine) {
        std::cerr << "CasError [" << constant_sine.error().operation << "]: "
                  << constant_sine.error().message << '\n';
    }
    EXPECT_TRUE(constant_sine && std::holds_alternative<UniversalSolutions>(constant_sine.value())) << "proved zero slope solves the constant equation rather than dividing by zero";
    auto constant_nonzero_sine = solve_transcendental(SymbolicExpr::add(
                                                          SymbolicExpr::sin(SymbolicExpr::multiply(a, x)), SymbolicExpr::number(-1)),
                                                      "x", zero_context);
    EXPECT_TRUE(constant_nonzero_sine &&
                std::holds_alternative<EmptySolutions>(constant_nonzero_sine.value()))
        << "proved zero slope preserves a nonzero constant residual as empty";
}

void test_numeric_expression_rhs(const ExprPtr &x) {
    auto two_as_expr = SymbolicExpr::add(SymbolicExpr::number(1), SymbolicExpr::number(1));
    auto exp_eq = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), two_as_expr));
    auto exp_roots = LMCAS::solve_transcendental(exp_eq, "x");
    EXPECT_TRUE(exp_roots && std::holds_alternative<FiniteSolutions>(exp_roots.value())) << "exp(x)=2 has a finite solution set";
    if (exp_roots && std::holds_alternative<FiniteSolutions>(exp_roots.value())) {
        const auto &values = std::get<FiniteSolutions>(exp_roots.value()).values;
        EXPECT_TRUE(values.size() == 1) << "real exponential inversion is unique";
        if (values.size() == 1) {
            auto expected = SymbolicExpr::ln(SymbolicExpr::number(2))->simplify();
            const auto actual = values.front().value->simplify();
            EXPECT_EQ(actual ? actual->to_string() : "null",
                      expected ? expected->to_string() : "null")
                << "exp(x)=1+1 gives ln(2)";
        }
    }
}

void test_empty_transcendental_dispatch(const ExprPtr &x) {
    const auto positive_expression = SymbolicExpr::exp(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x));
    const auto empty_direct = solve_transcendental(positive_expression, "x");
    EXPECT_TRUE(empty_direct &&
                std::holds_alternative<EmptySolutions>(empty_direct.value()))
        << "real exponential inversion proves exp(2*x)=0 empty";
    const auto empty_set = solve_equation(positive_expression, "x");
    EXPECT_TRUE(empty_set && std::holds_alternative<EmptySolutions>(empty_set.value())) << "equation dispatch preserves the mathematical Empty variant";
    const auto empty_finite = solve_finite_checked(positive_expression, "x");
    EXPECT_TRUE(empty_finite && empty_finite.value().empty()) << "finite projection maps proved Empty to an empty vector";
}

void test_logarithmic_affine_inverse(const ExprPtr &twice_x) {
    auto logarithmic = solve_transcendental(SymbolicExpr::ln(twice_x), "x");
    EXPECT_TRUE(logarithmic && std::holds_alternative<FiniteSolutions>(logarithmic.value())) << "logarithmic inversion still yields its unique finite root";
    if (logarithmic && std::holds_alternative<FiniteSolutions>(logarithmic.value())) {
        const auto &values = std::get<FiniteSolutions>(logarithmic.value()).values;
        EXPECT_TRUE(values.size() == 1) << "ln(2*x)=0 is unique";
        if (values.size() == 1) {
            EXPECT_TRUE(test_expr_equivalent(values.front().value, SymbolicExpr::number(Rational(1, 2)))) << "ln(2*x)=0 gives x=1/2";
        }
    }
}

void test_direct_product_domains(const ExprPtr &x) {
    auto valid_quotient = solve_transcendental(
        SymbolicExpr::divide(SymbolicExpr::ln(x), x), "x");
    EXPECT_TRUE(valid_quotient && std::holds_alternative<FiniteSolutions>(valid_quotient.value())) << "direct product solving retains logarithmic roots in the denominator domain";
    if (valid_quotient && std::holds_alternative<FiniteSolutions>(valid_quotient.value())) {
        const auto &values = std::get<FiniteSolutions>(valid_quotient.value()).values;
        EXPECT_TRUE(values.size() == 1) << "ln(x)/x has exactly one real zero";
        if (values.size() == 1) {
            EXPECT_TRUE(test_expr_equivalent(values.front().value, SymbolicExpr::number(1))) << "the surviving logarithmic root is one";
        }
    }
    auto excluded_quotient = solve_transcendental(SymbolicExpr::divide(SymbolicExpr::ln(x),
                                                                       SymbolicExpr::add(x, SymbolicExpr::number(-1))),
                                                  "x");
    EXPECT_TRUE(excluded_quotient && std::holds_alternative<EmptySolutions>(excluded_quotient.value())) << "the zero of ln(x) is excluded by its original denominator x-1";
}

} // namespace

TEST(SolveTranscendental, DirectPeriodicFamily) {
    auto x = SymbolicExpr::variable("x");
    test_direct_periodic_family(x);
}

TEST(SolveTranscendental, AffinePeriodicFamilies) {
    auto x = SymbolicExpr::variable("x");
    auto twice_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
    test_affine_periodic_families(twice_x);
}

TEST(SolveTranscendental, EndpointsAndSubstitutionUnion) {
    auto x = SymbolicExpr::variable("x");
    auto twice_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
    test_periodic_endpoints_and_unions(x, twice_x);
}

TEST(SolveTranscendental, SymbolicRangeConstraints) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    test_symbolic_range_constraints(x, a);
}

TEST(SolveTranscendental, CaptureFreeIntegerParameter) {
    auto x = SymbolicExpr::variable("x");
    auto twice_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
    test_capture_free_integer_parameter(twice_x);
}

TEST(SolveTranscendental, UnresolvedPreimagesRemainComplete) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    test_unresolved_preimages(x, a);
}

TEST(SolveTranscendental, AssumptionProvedSlope) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    test_assumption_proved_slope(x, a);
}

TEST(SolveTranscendental, NumericExpressionRhs) {
    auto x = SymbolicExpr::variable("x");
    test_numeric_expression_rhs(x);
}

TEST(SolveTranscendental, ProvedEmptySurvivesTypedDispatch) {
    auto x = SymbolicExpr::variable("x");
    test_empty_transcendental_dispatch(x);
}

TEST(SolveTranscendental, LogarithmicAffineInverse) {
    auto x = SymbolicExpr::variable("x");
    auto twice_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
    test_logarithmic_affine_inverse(twice_x);
}

TEST(SolveTranscendental, DirectProductsRetainOriginalDomains) {
    auto x = SymbolicExpr::variable("x");
    test_direct_product_domains(x);
}

TEST(SolveTranscendental, SubstitutionResourceLimitSurvivesDispatch) {
    const auto x = SymbolicExpr::variable("x");
    const auto expression = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(1000), x)),
            SymbolicExpr::exp(x)),
        SymbolicExpr::number(1));

    const auto detected = detect_substitution(expression, "x");
    ASSERT_FALSE(detected);
    EXPECT_EQ(detected.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(detected.error().operation, "polynomial.convert");

    const auto transcendental = solve_transcendental(expression, "x");
    ASSERT_FALSE(transcendental);
    EXPECT_EQ(transcendental.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(transcendental.error().operation, detected.error().operation);

    const auto equation = solve_equation(expression, "x");
    ASSERT_FALSE(equation);
    EXPECT_EQ(equation.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(equation.error().operation, detected.error().operation);
}

TEST(SolveTranscendental, RedundantVariableCandidatePreservesConversionFailure) {
    const auto x = SymbolicExpr::variable("x");
    const auto expression = SymbolicExpr::add(
        SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(1000)), x),
        SymbolicExpr::number(1));

    const auto converted = symbolic_to_poly<SymbolicPolyCoeff>(expression, "x");
    ASSERT_FALSE(converted);
    EXPECT_EQ(converted.error().code, CasErrc::ResourceLimit);
    const auto detected = detect_substitution(expression, "x");
    ASSERT_FALSE(detected);
    EXPECT_EQ(detected.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(detected.error().operation, converted.error().operation);
}

TEST(SolveTranscendental, UnsupportedSubstitutionIsSuccessfulNoMatch) {
    const auto sine = SymbolicExpr::sin(SymbolicExpr::variable("x"));
    const auto expression = SymbolicExpr::add(
        SymbolicExpr::power(sine, SymbolicExpr::number(-1)),
        SymbolicExpr::power(sine, SymbolicExpr::number(2)));

    const auto detected = detect_substitution(expression, "x");
    ASSERT_TRUE(detected);
    EXPECT_FALSE(detected.value().has_value());
}

TEST(SolveTranscendental, ExponentialSubstitutionPreservesPolynomialCoefficients) {
    const auto x = SymbolicExpr::variable("x");
    const auto expression = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::exp(SymbolicExpr::multiply(SymbolicExpr::number(2), x)),
            SymbolicExpr::exp(x)),
        SymbolicExpr::number(1));

    const auto detected = detect_substitution(expression, "x");
    ASSERT_TRUE(detected);
    ASSERT_TRUE(detected.value().has_value());
    const auto &substitution = *detected.value();
    const auto polynomial = symbolic_to_poly<Rational>(
        substitution.poly_in_u, substitution.u_var);
    ASSERT_TRUE(polynomial);
    EXPECT_EQ(polynomial.value().coeffs,
              (std::vector<Rational>{Rational(1), Rational(1), Rational(1)}));
    EXPECT_TRUE(test_expr_equivalent(substitution.u_expr, SymbolicExpr::exp(x)));
}

TEST(SolveTranscendental, LambertWEquationHasUniqueRealRoot) {
    auto x = SymbolicExpr::variable("x");
    auto equation = SymbolicExpr::add(x, SymbolicExpr::exp(x));
    SolveOptions options;
    options.allow_numeric = true;
    options.tolerance = 1e-12;
    options.max_roots = -1;

    auto solved = solve_equation(equation, "x", options);
    ASSERT_TRUE(solved);
    auto finite = std::get_if<FiniteSolutions>(&solved.value());
    ASSERT_NE(finite, nullptr) << "x+exp(x)=0 has a finite real solution set";
    ASSERT_EQ(finite->values.size(), 1u)
        << "strict monotonicity gives one real root";
    ASSERT_TRUE(finite->values.front().value);
    if (!finite->values.front().conditions.empty()) {
        ADD_FAILURE() << "the defining LambertW identity must discharge "
                      << finite->values.front().conditions.front()->to_string();
    }

    auto root = evaluate_numeric(*finite->values.front().value);
    ASSERT_TRUE(root);
    ASSERT_TRUE(root.value().is_finite());
    EXPECT_NEAR(root.value().value, -0.5671432904097838, 1e-10);

    auto residual = evaluate_numeric(
        *equation, {{"x", root.value().value}});
    ASSERT_TRUE(residual);
    ASSERT_TRUE(residual.value().is_finite());
    EXPECT_LE(std::abs(residual.value().value), 1e-10);
}
