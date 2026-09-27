#define _USE_MATH_DEFINES
#include <cmath>
#include "test_common.hpp"
#include "symbolic_complex.hpp"
#include "numeric_evaluation.hpp"
#include "assumption_context.hpp"
#include "internal/visitors/normalization_visitor.hpp"

using namespace LMCAS;

static void expect_complex_components(const ComplexSymbolic &value,
                                      double real, double imaginary) {
    ASSERT_NE(value.real, nullptr);
    ASSERT_NE(value.imag, nullptr);
    auto real_value = evaluate_numeric(*value.real);
    auto imaginary_value = evaluate_numeric(*value.imag);
    EXPECT_TRUE((real_value && imaginary_value)) << "complex arithmetic components evaluate numerically";
    if (!real_value || !imaginary_value) {
        return;
    }
    {
        const double actual_value = (real_value.value().value);
        const double expected_value = (real);
        const double tolerance = (1e-12);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
    {
        const double actual_value = (imaginary_value.value().value);
        const double expected_value = (imaginary);
        const double tolerance = (1e-12);
        EXPECT_TRUE(std::isfinite(actual_value));
        EXPECT_NEAR(actual_value, expected_value, tolerance);
    }
}

TEST(SymbolicComplex, ComplexArithmetic) {
    // Pair 1: (1+2i) and (3+4i)
    auto a1 = LMCAS::make_complex(SymbolicExpr::number(1), SymbolicExpr::number(2));
    auto b1 = LMCAS::make_complex(SymbolicExpr::number(3), SymbolicExpr::number(4));

    // Pair 2: (2+0i) and (0+3i)
    auto a2 = LMCAS::make_complex(SymbolicExpr::number(2), SymbolicExpr::number(0));
    auto b2 = LMCAS::make_complex(SymbolicExpr::number(0), SymbolicExpr::number(3));

    expect_complex_components(complex_add_checked(a1, b1).value(), 4, 6);
    expect_complex_components(complex_add_checked(a2, b2).value(), 2, 3);
    expect_complex_components(complex_sub_checked(a1, b1).value(), -2, -2);
    expect_complex_components(complex_sub_checked(a2, b2).value(), 2, -3);
    expect_complex_components(complex_mul_checked(a1, b1).value(), -5, 10);
    expect_complex_components(complex_mul_checked(a2, b2).value(), 0, 6);
    expect_complex_components(complex_div_checked(a1, b1).value(), 11.0 / 25, 2.0 / 25);
    expect_complex_components(complex_div_checked(a2, b2).value(), 0, -2.0 / 3);
}

TEST(SymbolicComplex, ComplexConj) {
    auto z = LMCAS::make_complex(SymbolicExpr::number(3), SymbolicExpr::number(4));
    auto conj = LMCAS::complex_conj_checked(z);
    ASSERT_TRUE(conj.has_value());
    expect_complex_components(conj.value(), 3, -4);
}

TEST(SymbolicComplex, ComplexAbs) {
    auto exact = LMCAS::complex_abs_checked(
        LMCAS::make_complex(SymbolicExpr::number(3), SymbolicExpr::number(4)));
    ASSERT_TRUE((exact.has_value())) << "exact modulus construction succeeds";
    if (exact) {
        auto value = LMCAS::evaluate_numeric(*exact.value());
        EXPECT_TRUE((value && value.value().value == 5.0)) << "exact modulus evaluates to five";
    }
    for (double component : {1e200, 1e-200}) {
        auto modulus = LMCAS::complex_abs_checked(LMCAS::make_complex(
            SymbolicExpr::number(component), SymbolicExpr::number(component)));
        ASSERT_TRUE((modulus.has_value())) << "extreme modulus construction succeeds";
        if (modulus) {
            auto value = LMCAS::evaluate_numeric(*modulus.value());
            EXPECT_TRUE((value && value.value().is_finite() &&
                         std::abs(value.value().value / std::hypot(component, component) - 1) < 1e-14))
                << "symbolic complex modulus preserves finite nonzero extreme values";
        }
    }
}

TEST(SymbolicComplex, ComplexArg) {
    // arg(1+0i) = atan2(0, 1) = 0
    {
        auto z = LMCAS::make_complex(SymbolicExpr::number(1), SymbolicExpr::number(0));
        auto arg = LMCAS::complex_arg_checked(z);
        ASSERT_TRUE(arg.has_value());
        ASSERT_NE(arg.value(), nullptr);
        auto val = evaluate_numeric(*arg.value());
        ASSERT_TRUE(val.has_value());
        EXPECT_TRUE((std::abs(val.value().value - 0.0) < 1e-9))
            << "arg(1+0i) evaluates to 0";
    }

    // arg(0+1i) = atan2(1, 0) = pi/2
    {
        auto z = LMCAS::make_complex(SymbolicExpr::number(0), SymbolicExpr::number(1));
        auto arg = LMCAS::complex_arg_checked(z);
        ASSERT_TRUE(arg.has_value());
        ASSERT_NE(arg.value(), nullptr);
        auto val = evaluate_numeric(*arg.value());
        ASSERT_TRUE(val.has_value());
        EXPECT_TRUE((std::abs(val.value().value - M_PI / 2.0) < 1e-9))
            << "arg(0+1i) evaluates to pi/2";
    }
}

TEST(SymbolicComplex, ComplexPolarForms) {
    auto r = SymbolicExpr::number(2);
    auto theta = SymbolicExpr::number(1); // 1 radian

    auto exp_form = LMCAS::complex_exp_form_checked(r, theta);
    ASSERT_TRUE(exp_form.has_value());
    expect_complex_components(exp_form.value(), 2 * std::cos(1.0), 2 * std::sin(1.0));

    auto trig_form = LMCAS::complex_trig_form_checked(r, theta);
    ASSERT_TRUE(trig_form.has_value());
    expect_complex_components(trig_form.value(), 2 * std::cos(1.0), 2 * std::sin(1.0));
}

TEST(SymbolicComplex, ComplexNthRoot) {
    // Test n=2: square roots of 4
    {
        auto c = SymbolicExpr::number(4.0);
        auto roots = LMCAS::solve_complex_nth_root_checked(c, 2).value();
        EXPECT_TRUE((roots.size() == 2)) << "sqrt(4) returns exactly 2 roots";
    }

    // Test n=3: cube roots of 8
    {
        auto c = SymbolicExpr::number(8.0);
        auto roots = LMCAS::solve_complex_nth_root_checked(c, 3).value();
        EXPECT_TRUE((roots.size() == 3)) << "cbrt(8) returns exactly 3 roots";
    }

    // Test n=4: fourth roots of 16
    {
        auto c = SymbolicExpr::number(16.0);
        auto roots = LMCAS::solve_complex_nth_root_checked(c, 4).value();
        EXPECT_TRUE((roots.size() == 4)) << "4th root of 16 returns exactly 4 roots";
    }

    // Test n=5: fifth roots of 32
    {
        auto c = SymbolicExpr::number(32.0);
        auto roots = LMCAS::solve_complex_nth_root_checked(c, 5).value();
        EXPECT_TRUE((roots.size() == 5)) << "5th root of 32 returns exactly 5 roots";
    }
}

static void test_complex_invalid_components(const std::shared_ptr<SymbolicExpr> &one, const ComplexSymbolic &b) {
    auto bad = LMCAS::make_complex(nullptr, one);
    auto result = LMCAS::complex_mul_checked(bad, b);
    EXPECT_TRUE((!result.has_value() &&
                 result.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked complex_mul rejects null components";
}

static void test_complex_zero_division(const std::shared_ptr<SymbolicExpr> &zero, const ComplexSymbolic &a) {
    auto zero_complex = LMCAS::make_complex(zero, zero);
    auto result = LMCAS::complex_div_checked(a, zero_complex);
    EXPECT_TRUE((!result.has_value() &&
                 result.error().code == LMCAS::CasErrc::DomainError))
        << "checked complex_div rejects exact zero denominator";
}

static void test_nth_root_domains(const std::shared_ptr<SymbolicExpr> &approx_four) {
    {
        auto roots = LMCAS::solve_complex_nth_root_checked(approx_four, 2);
        EXPECT_TRUE((roots.has_value() && roots.value().size() == 2)) << "checked complex nth root accepts explicit approximate real input";
    }

    {
        auto exact = LMCAS::solve_complex_nth_root_checked(SymbolicExpr::number(4), 2);
        EXPECT_TRUE((!exact.has_value() &&
                     exact.error().code == LMCAS::CasErrc::Inconclusive))
            << "checked complex nth root does not implicitly float exact integers";
    }

    {
        auto bad_order = LMCAS::solve_complex_nth_root_checked(approx_four, 0);
        EXPECT_TRUE((!bad_order.has_value() &&
                     bad_order.error().code == LMCAS::CasErrc::InvalidArgument))
            << "checked complex nth root rejects non-positive degree";
    }
}

static void test_complex_context_errors(const ComplexSymbolic &a) {
    {
        LMCAS::CancellationToken token;
        token.cancel();
        LMCAS::ComputationContext cancelled_context({}, token);
        auto cancelled = LMCAS::complex_conj_checked(a, cancelled_context);
        EXPECT_TRUE((!cancelled.has_value() &&
                     cancelled.error().code == LMCAS::CasErrc::Cancelled))
            << "checked complex_conj observes cancellation";
    }

    {
        LMCAS::ResourceLimits limits;
        limits.max_steps = 1;
        LMCAS::ComputationContext limited_context(limits);
        auto limited = LMCAS::complex_abs_checked(a, limited_context);
        EXPECT_TRUE((!limited.has_value() &&
                     limited.error().code == LMCAS::CasErrc::ResourceLimit))
            << "checked complex_abs observes step budget";
    }
}

TEST(SymbolicComplex, CheckedComplexContracts) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto two = SymbolicExpr::number(2);
    auto approx_four = SymbolicExpr::number(4.0);

    auto a = LMCAS::make_complex(one, two);
    auto b = LMCAS::make_complex(two, one);

    {
        auto sum = LMCAS::complex_add_checked(a, b);
        EXPECT_TRUE((sum.has_value())) << "checked complex_add succeeds for valid inputs";
        EXPECT_TRUE((sum.value().real && sum.value().imag)) << "checked complex_add returns non-null components";
    }
    test_complex_invalid_components(one, b);
    test_complex_zero_division(zero, a);
    test_nth_root_domains(approx_four);
    test_complex_context_errors(a);
}

void expect_complex_zero(const ComplexSymbolic &value,
                         const AssumptionContext *facts = nullptr) {
    for (const auto &component : {value.real, value.imag}) {
        auto normalized = facts ? facts->simplify(*component) : component->simplify();
        EXPECT_TRUE((normalized && normalized->is_zero())) << "Cartesian polynomial residual component is zero under its assumptions";
    }
}

TEST(SymbolicComplex, ComplexQuadratic) {
    auto one = SymbolicExpr::number(1);
    auto zero = SymbolicExpr::number(0);
    auto roots = solve_complex_quadratic_checked(one, zero, one);
    EXPECT_TRUE((roots && roots.value().size() == 2)) << "z^2+1 returns two verified roots";
    if (!roots || roots.value().size() != 2) {
        return;
    }
    for (std::size_t index = 0; index < 2; ++index) {
        const auto &root = roots.value()[index];
        EXPECT_TRUE((root.real->is_zero())) << "root has exactly zero real component";
        auto expected_imag = SymbolicExpr::number(index == 0 ? -1 : 1);
        EXPECT_TRUE((SymbolicExpr::add(root.imag, expected_imag)->simplify()->is_zero())) << "formula ordering yields +i followed by -i";
        auto squared = complex_mul_checked(root, root);
        ASSERT_TRUE((squared.has_value())) << "root can be squared with complex arithmetic";
        if (!squared) {
            continue;
        }
        auto residual = complex_add_checked(squared.value(), make_complex(one, zero));
        ASSERT_TRUE((residual.has_value())) << "polynomial residual can be constructed";
        if (residual)
            expect_complex_zero(residual.value());
        auto modulus = complex_abs_checked(root);
        EXPECT_TRUE((modulus && modulus.value()->simplify()->is_one())) << "each root has exact modulus one";
    }
    auto conjugate = complex_conj_checked(roots.value()[0]);
    ASSERT_TRUE((conjugate.has_value())) << "conjugation succeeds";
    if (!conjugate) {
        return;
    }
    auto difference = complex_sub_checked(conjugate.value(), roots.value()[1]);
    ASSERT_TRUE((difference.has_value())) << "conjugate is comparable with the other root";
    if (difference)
        expect_complex_zero(difference.value());
}

static bool has_two_cartesian_roots(const LMCAS::ComplexRootsResult &roots) {
    return roots && roots.value().size() == 2;
}

TEST(SymbolicComplex, ComplexQuadraticRealRoots) {
    auto roots = solve_complex_quadratic_checked(
        SymbolicExpr::number(1), SymbolicExpr::number(-3), SymbolicExpr::number(2));
    EXPECT_TRUE((has_two_cartesian_roots(roots))) << "ordinary real quadratic succeeds";
    if (!has_two_cartesian_roots(roots)) {
        return;
    }
    EXPECT_TRUE((roots.value()[0].imag->is_zero() && roots.value()[1].imag->is_zero())) << "nonnegative discriminant gives real roots";
    EXPECT_TRUE((SymbolicExpr::add(roots.value()[0].real, SymbolicExpr::number(-2))
                     ->simplify()
                     ->is_zero()))
        << "first real root is two";
    EXPECT_TRUE((roots.value()[1].real->is_one())) << "second real root is one";
    auto repeated = solve_complex_quadratic_checked(
        SymbolicExpr::number(1), SymbolicExpr::number(-2), SymbolicExpr::number(1));
    EXPECT_TRUE((has_two_cartesian_roots(repeated))) << "Cartesian API retains repeated-root multiplicity";
    if (!has_two_cartesian_roots(repeated)) {
        return;
    }
    for (const auto &root : repeated.value()) {
        EXPECT_TRUE((root.real->is_one() && root.imag->is_zero())) << "repeated root is one";
    }
}

TEST(SymbolicComplex, ComplexQuadraticIrrationalRealRoots) {
    auto roots = solve_complex_quadratic_checked(
        SymbolicExpr::number(1), SymbolicExpr::number(-2), SymbolicExpr::number(-1));
    ASSERT_TRUE(has_two_cartesian_roots(roots));
    expect_complex_components(roots.value()[0], 1 + std::sqrt(2.0), 0);
    expect_complex_components(roots.value()[1], 1 - std::sqrt(2.0), 0);
}

TEST(SymbolicComplex, ComplexQuadraticNegativeLeadingCoefficient) {
    auto roots = solve_complex_quadratic_checked(
        SymbolicExpr::number(-2), SymbolicExpr::number(0), SymbolicExpr::number(-2));
    EXPECT_TRUE((roots && roots.value().size() == 2)) << "negative-leading quadratic succeeds";
    if (!roots || roots.value().size() != 2) {
        return;
    }
    EXPECT_TRUE((roots.value()[0].real->is_zero() && roots.value()[1].real->is_zero())) << "negative leading coefficient does not introduce a real component";
    EXPECT_TRUE((SymbolicExpr::add(roots.value()[0].imag, SymbolicExpr::number(1))
                     ->simplify()
                     ->is_zero()))
        << "first signed-formula root is -i";
    EXPECT_TRUE((roots.value()[1].imag->is_one())) << "second signed-formula root is +i";
}

TEST(SymbolicComplex, ComplexQuadraticRejectedCoefficients) {
    auto zero = SymbolicExpr::number(0);
    auto one = SymbolicExpr::number(1);
    auto missing = solve_complex_quadratic_checked(one, nullptr, zero);
    EXPECT_TRUE((!missing && missing.error().code == CasErrc::InvalidArgument)) << "null coefficient remains InvalidArgument";
    auto degenerate = solve_complex_quadratic_checked(zero, one, one);
    EXPECT_TRUE((!degenerate && degenerate.error().code == CasErrc::DomainError)) << "zero leading coefficient remains DomainError";
    auto nonreal = solve_complex_quadratic_checked(
        one, zero, SymbolicExpr::sqrt(SymbolicExpr::number(-1)));
    EXPECT_TRUE((!nonreal && nonreal.error().code == CasErrc::UnsupportedExpression)) << "general complex coefficients are explicitly unsupported";
}

TEST(SymbolicComplex, ComplexQuadraticUnknownSign) {
    auto facts = std::make_shared<AssumptionContext>();
    auto declared = facts->assume_domain_checked("b", Domain::Real);
    EXPECT_TRUE((declared.has_value())) << "real parameter declaration succeeds";
    ComputationContext context;
    auto attached = context.set_assumptions(facts);
    EXPECT_TRUE((attached.has_value())) << "assumptions attach to shared computation context";
    auto one = SymbolicExpr::number(1);
    auto b = SymbolicExpr::variable("b");
    auto unknown = solve_complex_quadratic_checked(one, b, one, context);
    EXPECT_TRUE((!unknown && unknown.error().code == CasErrc::Inconclusive)) << "unknown sign of b^2-4 returns Inconclusive without roots";
    auto leading = solve_complex_quadratic_checked(b, one, one, context);
    EXPECT_TRUE((!leading && leading.error().code == CasErrc::Inconclusive)) << "a possibly zero real leading coefficient returns Inconclusive";
}

TEST(SymbolicComplex, ComplexQuadraticKnownSymbolicSign) {
    auto facts = std::make_shared<AssumptionContext>();
    auto domain = facts->assume_domain_checked("p", Domain::Real);
    auto sign = facts->assume_sign_checked("p", Sign::Positive);
    EXPECT_TRUE((domain && sign)) << "positive real parameter declarations succeed";
    ComputationContext context;
    auto attached = context.set_assumptions(facts);
    EXPECT_TRUE((attached.has_value())) << "symbolic quadratic uses caller assumptions";
    auto p = SymbolicExpr::variable("p");
    auto roots = solve_complex_quadratic_checked(
        SymbolicExpr::number(1), SymbolicExpr::number(0), p, context);
    EXPECT_TRUE((roots && roots.value().size() == 2)) << "known negative symbolic discriminant returns both roots";
    if (!roots) {
        return;
    }
    for (const auto &root : roots.value()) {
        EXPECT_TRUE((root.real->is_zero())) << "symbolic imaginary roots have zero real part";
        auto squared = complex_mul_checked(root, root, context);
        ASSERT_TRUE((squared.has_value())) << "symbolic root supports complex multiplication";
        if (!squared) {
            continue;
        }
        auto residual = complex_add_checked(
            squared.value(), make_complex(p, SymbolicExpr::number(0)), context);
        ASSERT_TRUE((residual.has_value())) << "symbolic polynomial residual constructs";
        if (residual)
            expect_complex_zero(residual.value(), facts.get());
    }
}

TEST(SymbolicComplex, ComplexQuadraticContextFailures) {
    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext limited(limits);
    auto one = SymbolicExpr::number(1);
    auto exhausted = solve_complex_quadratic_checked(one, one, one, limited);
    EXPECT_TRUE((!exhausted && exhausted.error().code == CasErrc::ResourceLimit)) << "exhausted quadratic computation returns ResourceLimit";
    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled({}, cancellation);
    auto stopped = solve_complex_quadratic_checked(one, one, one, cancelled);
    EXPECT_TRUE((!stopped && stopped.error().code == CasErrc::Cancelled)) << "cancelled quadratic computation returns Cancelled";
}

TEST(SymbolicComplex, ComplexQuadraticNormalizationNodeLimit) {
    ResourceLimits limits;
    limits.max_ast_nodes = 0;
    ComputationContext context(limits);
    auto roots = solve_complex_quadratic_checked(
        SymbolicExpr::number(1), SymbolicExpr::number(0), SymbolicExpr::number(1), context);
    ASSERT_FALSE(roots);
    EXPECT_EQ(roots.error().code, CasErrc::ResourceLimit);
}

static void expect_locus_at(const std::shared_ptr<SymbolicExpr> &locus,
                            double z, bool on_locus) {
    SCOPED_TRACE(z);
    ASSERT_NE(locus, nullptr);
    auto relation = std::dynamic_pointer_cast<const RelationalNode>(detail::node(locus));
    ASSERT_NE(relation, nullptr);
    ASSERT_EQ(relation->op(), RelationalNode::Op::EQ);
    ASSERT_NE(relation->left(), nullptr);
    ASSERT_NE(relation->right(), nullptr);
    const NumericBindings bindings{{"z", z}};
    auto left = evaluate_numeric(*detail::make_expression_ptr(relation->left()), bindings);
    auto right = evaluate_numeric(*detail::make_expression_ptr(relation->right()), bindings);
    ASSERT_TRUE(left.has_value());
    ASSERT_TRUE(right.has_value());
    ASSERT_TRUE(left.value().is_finite());
    ASSERT_TRUE(right.value().is_finite());
    if (on_locus) {
        EXPECT_NEAR(left.value().value, right.value().value, 1e-12);
    } else {
        EXPECT_GT(std::abs(left.value().value - right.value().value), 1e-6);
    }
}

TEST(SymbolicComplex, ComplexLocus) {
    // Circle centered at (1+2i) with radius 3
    {
        auto center = LMCAS::make_complex(SymbolicExpr::number(1), SymbolicExpr::number(2));
        auto radius = SymbolicExpr::number(3);
        auto checked_locus =
            LMCAS::complex_locus_circle_checked(center, radius, "z");
        ASSERT_TRUE(checked_locus.has_value()) << "checked circle locus succeeds";
        expect_locus_at(checked_locus.value(), 1 + std::sqrt(5.0), true);
        expect_locus_at(checked_locus.value(), 1 - std::sqrt(5.0), true);
        expect_locus_at(checked_locus.value(), 1, false);
    }

    // Perpendicular bisector between (1+0i) and (3+0i)
    {
        auto a = LMCAS::make_complex(SymbolicExpr::number(1), SymbolicExpr::number(0));
        auto b = LMCAS::make_complex(SymbolicExpr::number(3), SymbolicExpr::number(0));
        auto checked_locus =
            LMCAS::complex_locus_perpendicular_bisector_checked(a, b, "z");
        ASSERT_TRUE(checked_locus.has_value()) << "checked perpendicular-bisector locus succeeds";
        expect_locus_at(checked_locus.value(), 2, true);
        expect_locus_at(checked_locus.value(), 0, false);
    }
}
