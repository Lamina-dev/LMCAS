#include "test_common.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include "value.hpp"
#include "irrational.hpp"
#include "symbolic.hpp"
#include "internal/exact_constant_bounds.hpp"
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

using namespace LMCAS;

TEST(IrrationalRegressions, IrrationalToSymbolicComplex) {
    // 2 + 3*pi + 1*sqrt(2)
    Irrational mix = Irrational::constant(2.0) + Irrational::pi(3.0) +
                     Irrational::sqrt_checked(BigInt(2)).value();
    EXPECT_TRUE((std::abs(mix.to_symbolic()->to_numeric() -
                          (2.0 + 3.0 * LMMC_PI + std::sqrt(2.0))) < 1e-12))
        << "symbolic reconstruction retains every combination term";

    Irrational c = Irrational::constant(5.0);
    EXPECT_TRUE((c.to_symbolic()->to_numeric() == 5.0)) << "constant-only reconstruction retains its value";
}

TEST(IrrationalRegressions, IrrationalExactSymbolicCoefficients) {
    const double tiny = std::ldexp(1.0, -60);
    for (double coefficient : {tiny, -tiny, std::nextafter(1.0, 2.0),
                               std::nextafter(1.0, 0.0)}) {
        EXPECT_TRUE((Irrational::constant(coefficient).to_symbolic()->to_numeric() == coefficient)) << "constant reconstruction retains the exact binary coefficient";

        for (auto value : {Irrational::sqrt_checked(BigInt(1), coefficient).value(), Irrational::sqrt_checked(BigInt(2), coefficient).value(),
                           Irrational::pi(coefficient), Irrational::e(coefficient)}) {
            const double expected = value.to_double();
            EXPECT_TRUE((value.to_symbolic()->to_numeric() == expected)) << "direct reconstruction neither drops nor rounds the coefficient";
            value.to_complex();
            EXPECT_TRUE((value.to_symbolic()->to_numeric() == expected)) << "combination reconstruction neither drops nor rounds the coefficient";
        }
    }

    const auto mix = Irrational::constant(tiny) + Irrational::sqrt_checked(BigInt(2), tiny).value();
    const double expected = tiny + tiny * std::sqrt(2.0);
    EXPECT_TRUE((std::abs(mix.to_symbolic()->to_numeric() - expected) <=
                 tiny * 4.0 * std::numeric_limits<double>::epsilon()))
        << "small constant and radical terms both survive in the same combination";
}

TEST(IrrationalRegressions, IrrationalLargeIntegralDisplay) {
    for (double coefficient : {4294967291.0, -4294967291.0,
                               std::ldexp(1.0, 64), -std::ldexp(1.0, 64),
                               std::numeric_limits<double>::max(),
                               -std::numeric_limits<double>::max()}) {
        const auto constant_text = Irrational::constant(coefficient).to_string();
        std::size_t consumed = 0;
        const double displayed_constant = std::stod(constant_text, &consumed);
        EXPECT_TRUE((displayed_constant == coefficient && consumed == constant_text.size())) << "displayed constant parses back to its signed finite value";

        for (auto value : {Irrational::sqrt_checked(BigInt(1), coefficient).value(), Irrational::sqrt_checked(BigInt(2), coefficient).value(),
                           Irrational::pi(coefficient), Irrational::e(coefficient)}) {
            EXPECT_TRUE((std::stod(value.to_string()) == coefficient)) << "direct basis display retains the signed integral coefficient";
            value.to_complex();
            EXPECT_TRUE((std::stod(value.to_string()) == coefficient)) << "combination basis display retains the signed integral coefficient";
        }

        const auto mix = Irrational::constant(coefficient) + Irrational::sqrt_checked(BigInt(2), -coefficient).value();
        std::istringstream displayed(mix.to_string());
        double constant = 0.0;
        double magnitude = 0.0;
        char operation = 0;
        displayed >> constant >> operation >> magnitude;
        const double signed_coefficient = operation == '-' ? -magnitude : magnitude;
        EXPECT_TRUE((!displayed.fail() && (operation == '+' || operation == '-') &&
                     constant == coefficient && signed_coefficient == -coefficient))
            << "mixed display retains both large values and their opposite signs";
    }
}

TEST(IrrationalRegressions, IrrationalLargeRadicals) {
    const BigInt prime("4294967291");
    const auto root = Irrational::sqrt_checked(prime).value();
    const auto square = root * root;
    EXPECT_TRUE((square.as_rational_checked().value() == Rational(prime))) << "equal large radicals multiply to the exact integer";
    EXPECT_TRUE((square.to_symbolic()->compare(SymbolicExpr::number(prime)) == 0)) << "cancelled radical product reconstructs the exact integer";

    const auto partial = Irrational::sqrt_checked(BigInt(2) * prime).value() *
                         Irrational::sqrt_checked(BigInt(3) * prime).value();
    EXPECT_TRUE((partial == Irrational::sqrt_checked(BigInt(6)).value() * Rational(prime))) << "shared factors cancel while distinct radical factors remain";
    const BigInt other_radicand = BigInt(3) * BigInt("2147483647");
    const auto other = Irrational::sqrt_checked(other_radicand).value();
    const auto product = root * other;
    const BigInt exact_product = prime * other_radicand;
    EXPECT_TRUE((product.pow(2).as_rational_checked().value() == Rational(exact_product))) << "the square of a formerly overflowing product is exact";
    EXPECT_TRUE((product.to_symbolic()->compare(
                     SymbolicExpr::sqrt(SymbolicExpr::number(exact_product))) == 0))
        << "symbolic reconstruction retains the complete greater-than-64-bit radicand";
    const auto zero = Irrational::sqrt_checked(BigInt(0)).value();
    EXPECT_TRUE(((zero * zero).is_zero())) << "zero radicals multiply without division by zero";
    EXPECT_TRUE(((Irrational::sqrt_checked(prime, 0.0).value() * other).is_zero())) << "zero coefficients remain exact through large products";
    const BigInt huge_root = BigInt(1) << 100;
    EXPECT_TRUE((Irrational::sqrt_checked(huge_root * huge_root).value().as_rational_checked().value() == Rational(huge_root))) << "perfect-square extraction is unbounded";
}

TEST(IrrationalRegressions, IrrationalCheckedLimits) {
    const auto negative = Irrational::sqrt_checked(-(BigInt(1) << 100));
    EXPECT_TRUE((!negative && negative.error().code == CasErrc::DomainError)) << "negative unbounded radicands report DomainError";
    ResourceLimits limits;
    limits.max_steps = 2;
    ComputationContext steps(limits);
    const auto exhausted = Irrational::sqrt_checked(BigInt(97), 1.0, steps);
    EXPECT_TRUE((!exhausted && exhausted.error().code == CasErrc::ResourceLimit)) << "factor extraction exhausts steps rather than returning an approximation";
    limits.max_steps = 100;
    limits.max_integer_bits = 5;
    ComputationContext bits(limits);
    const auto oversized = Irrational::sqrt_checked(BigInt(64), 1.0, bits);
    EXPECT_TRUE((!oversized && oversized.error().code == CasErrc::ResourceLimit)) << "input radicands obey the integer-bit budget";
    ComputationContext coefficient_bits(limits);
    const auto oversized_coefficient =
        Irrational::sqrt_checked(BigInt(4), 16.0, coefficient_bits);
    EXPECT_TRUE((!oversized_coefficient &&
                 oversized_coefficient.error().code == CasErrc::ResourceLimit))
        << "extracted rational coefficients obey the integer-bit budget";
    CancellationToken token;
    token.cancel();
    ComputationContext cancelled({}, token);
    const auto stopped = Irrational::sqrt_checked(BigInt(2), 1.0, cancelled);
    EXPECT_TRUE((!stopped && stopped.error().code == CasErrc::Cancelled)) << "cancellation is propagated without constructing an approximate radical";
}

TEST(IrrationalRegressions, IrrationalExactArithmetic) {
    const auto root = Irrational::sqrt_checked(BigInt(2)).value();
    const auto third = Irrational::constant(Rational(1, 3));
    const auto sum = root + third;
    EXPECT_TRUE(((sum * (root - third)).as_rational_checked().value() == Rational(17, 9))) << "cross terms cancel exactly in a product of radical sums";
    EXPECT_TRUE((((sum / (root - third)) * (root - third)) == sum)) << "division of combinations retains its exact denominator";
    const auto pi_e = Irrational::pi() * Irrational::e();
    EXPECT_TRUE(((pi_e / Irrational::e()) == Irrational::pi())) << "products and quotients preserve pi and e bases";
    EXPECT_TRUE((root.pow(-2).as_rational_checked().value() == Rational(1, 2))) << "negative powers retain rational results";
    EXPECT_THROW(([&] { (void)(sum / (root - root)); })(), std::domain_error) << "exact zero denominators are rejected";
    const BigInt huge = (BigInt(1) << 100) + BigInt(1);
    const auto wrapped = Value(huge).as_irrational_checked().value();
    EXPECT_TRUE((Value(wrapped).as_rational_checked().value() == Rational(huge))) << "Value preserves an integer beyond floating-point precision";
    const auto fraction = Value(Rational(1, 3)).as_irrational_checked().value();
    EXPECT_TRUE((Value(fraction).as_rational_checked().value() == Rational(1, 3))) << "Value round-trips an exact non-binary rational";
    EXPECT_TRUE((!Value(root).as_rational_checked())) << "a genuine radical cannot silently become a rounded rational";
}

static void expect_radical_quotient_signs(
    const Irrational &three, const Irrational &positive, const Irrational &exact_zero) {
    for (int numerator_sign : {-1, 1}) {
        for (int denominator_sign : {-1, 1}) {
            const auto quotient = (three * Rational(numerator_sign)) /
                                  (positive * Rational(denominator_sign));
            auto sign = quotient.sign_checked();
            EXPECT_TRUE((sign && sign.value() == (numerator_sign == denominator_sign
                                                      ? Sign::Positive
                                                      : Sign::Negative)))
                << "near-zero denominator sign is certified separately from the numerator";
            const auto magnitude = quotient.abs_checked();
            EXPECT_TRUE((magnitude && magnitude.value() ==
                                          (numerator_sign == denominator_sign ? quotient : -quotient)))
                << "quotient absolute value preserves exact arithmetic";
        }
    }
    auto sign = (exact_zero / positive).sign_checked();
    EXPECT_TRUE((sign && sign.value() == Sign::Zero)) << "zero over a nonzero radical is zero";
}

static void expect_near_cancellation_signs(
    const Irrational &negative, const Irrational &positive) {
    const auto negative_sign = negative.sign_checked();
    const auto positive_sign = positive.sign_checked();
    EXPECT_TRUE((negative_sign && negative_sign.value() == Sign::Negative)) << "the formal near-cancellation difference is strictly negative";
    EXPECT_TRUE((positive_sign && positive_sign.value() == Sign::Positive)) << "the adjacent lower rational gives a strictly positive difference";
    const auto absolute = negative.abs_checked();
    EXPECT_TRUE((absolute && absolute.value() == -negative && negative.abs() == -negative)) << "both absolute-value entries use the certified sign";
}

TEST(IrrationalRegressions, IrrationalCertifiedSign) {
    const auto root2 = Irrational::sqrt_checked(BigInt(2)).value();
    const auto root3 = Irrational::sqrt_checked(BigInt(3)).value();
    const auto root6 = Irrational::sqrt_checked(BigInt(6)).value();
    const Rational above(BigInt("6369051672525773"), BigInt("4503599627370496"));
    const Rational below(BigInt("6369051672525772"), BigInt("4503599627370496"));
    EXPECT_TRUE((above * above > Rational(2) && below * below < Rational(2))) << "the exact square inequalities establish independent sign oracles";
    const auto negative = root2 - Irrational::constant(above);
    const auto positive = root2 - Irrational::constant(below);
    expect_near_cancellation_signs(negative, positive);
    const auto three = root2 + root3 - root6;
    auto sign = three.sign_checked();
    EXPECT_TRUE((sign && sign.value() == Sign::Positive)) << "three distinct radicals retain their sign";
    const auto exact_zero = (root2 + root3).pow(2) -
                            Irrational::constant(Rational(5)) - root6 * Rational(2);
    sign = exact_zero.sign_checked();
    EXPECT_TRUE((sign && sign.value() == Sign::Zero)) << "square expansion cancels sqrt(6) exactly";
    expect_radical_quotient_signs(three, positive, exact_zero);
}

static void expect_transcendental_sign_values(
    const Irrational &pi, const Irrational &e,
    const Irrational &root2, const Irrational &root3) {
    for (const auto &value : {pi - e, pi / e, pi.pow(-2) * e.pow(3),
                              (root3 - root2) * pi / ((pi - e) * e)}) {
        const auto sign = value.sign_checked();
        EXPECT_TRUE((sign && sign.value() == Sign::Positive)) << "positive mixed constant expressions have certified signs";
        const auto opposite = (-value).sign_checked();
        EXPECT_TRUE((opposite && opposite.value() == Sign::Negative)) << "negation reverses the certified sign";
    }
    for (const auto &value : {pi - pi, e - e}) {
        const auto sign = value.sign_checked();
        EXPECT_TRUE((sign && sign.value() == Sign::Zero)) << "structural constant cancellation is exact zero";
        EXPECT_THROW(([&] { (void)(root2 / value); })(), std::domain_error) << "structural zero denominators remain errors";
    }
}

TEST(IrrationalRegressions, IrrationalTranscendentalSign) {
    const auto pi = Irrational::pi();
    const auto e = Irrational::e();
    const auto root2 = Irrational::sqrt_checked(BigInt(2)).value();
    const auto root3 = Irrational::sqrt_checked(BigInt(3)).value();
    expect_transcendental_sign_values(pi, e, root2, root3);
    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext steps(limits);
    auto sign = (pi - e).sign_checked(steps);
    EXPECT_TRUE((!sign && sign.error().code == CasErrc::ResourceLimit)) << "sign computation cannot bypass the step budget";
    limits.max_steps = 10000;
    limits.max_integer_bits = 8;
    ComputationContext bits(limits);
    sign = root2.sign_checked(bits);
    EXPECT_TRUE((!sign && sign.error().code == CasErrc::ResourceLimit)) << "dyadic refinement checks integer growth before shifting";
    const auto denominator = pi - e;
    const auto zero_quotient = Irrational::constant(Rational(0)) / denominator;
    const auto proportional = denominator / denominator;
    for (const auto &quotient : {zero_quotient, proportional}) {
        ComputationContext denominator_bits(limits);
        sign = quotient.sign_checked(denominator_bits);
        EXPECT_TRUE((!sign && sign.error().code == CasErrc::ResourceLimit)) << "normalization cannot erase an unresolved pi/e denominator";
    }
    sign = zero_quotient.sign_checked();
    EXPECT_TRUE((sign && sign.value() == Sign::Zero)) << "zero numerator succeeds after certifying the retained denominator";
    sign = proportional.sign_checked();
    EXPECT_TRUE((sign && sign.value() == Sign::Positive)) << "proportional terms succeed after certifying the retained denominator";
    limits = ResourceLimits{};
    limits.max_expansion_terms = 1;
    ComputationContext terms(limits);
    sign = (root2 + root3).sign_checked(terms);
    EXPECT_TRUE((!sign && sign.error().code == CasErrc::ResourceLimit)) << "multi-term sign computation respects the expansion budget";
    CancellationToken token;
    token.cancel();
    ComputationContext cancelled({}, token);
    const auto absolute = (pi - e).abs_checked(cancelled);
    EXPECT_TRUE((!absolute && absolute.error().code == CasErrc::Cancelled)) << "absolute value propagates cancellation rather than an approximate sign";
}

static bool encloses_reference_with_width(
    const std::pair<Rational, Rational> &bounds, const Rational &lower,
    const Rational &upper, const Rational &tolerance) {
    return bounds.first <= lower && bounds.second >= upper &&
           bounds.second - bounds.first <= tolerance;
}

static void expect_atan_series_envelopes(
    ComputationContext &context, const Rational &tolerance) {
    const Rational atan_half_lower("0.46364760900080611621425623146121440202853705428612");
    const Rational atan_half_upper("0.46364760900080611621425623146121440202853705428613");
    const Rational atan_two_lower("1.10714871779409050301706546017853704007004764540143");
    const Rational atan_two_upper("1.10714871779409050301706546017853704007004764540144");
    const struct {
        Rational x;
        const Rational &lower;
        const Rational &upper;
    } cases[] = {{Rational(1, 2), atan_half_lower, atan_half_upper},
                 {Rational(2), atan_two_lower, atan_two_upper}};
    for (const auto &item : cases) {
        const auto &x = item.x;
        auto bounds = detail::exact_atan_bounds(x, 128, context);
        EXPECT_TRUE((bounds && encloses_reference_with_width(
                                   bounds.value(), item.lower, item.upper, tolerance)))
            << "atan small-series and reciprocal branches enclose independent exact brackets";
        auto opposite = detail::exact_atan_bounds(-x, 128, context);
        EXPECT_TRUE((bounds && opposite && opposite.value().first == -bounds.value().second &&
                     opposite.value().second == -bounds.value().first))
            << "atan bounds respect odd symmetry";
    }
}

static void expect_atan_transformed_envelopes(
    ComputationContext &context, const Rational &tolerance,
    const Rational &pi_lower, const Rational &pi_upper) {
    {
        auto unit = detail::exact_atan_bounds(Rational(1), 128, context);
        EXPECT_TRUE((unit && unit.value().first <= pi_lower / Rational(4) &&
                     unit.value().second >= pi_upper / Rational(4) &&
                     unit.value().second - unit.value().first <= tolerance))
            << "the pi/4 transform retains exact endpoint enclosures";
    }
    {
        auto middle = detail::exact_atan_bounds(Rational(3, 4), 128, context);
        EXPECT_TRUE((middle &&
                     middle.value().first <= Rational("0.64350110879328438680280922871732263804151059111531") &&
                     middle.value().second >= Rational("0.64350110879328438680280922871732263804151059111532") &&
                     middle.value().second - middle.value().first <= tolerance))
            << "the nonzero pi/4 correction encloses an independent decimal bracket";
    }
    auto zero = detail::exact_atan_bounds(Rational(0), 0, context);
    EXPECT_TRUE((zero && zero.value().first.is_zero() && zero.value().second.is_zero())) << "atan zero has an exact zero interval";
}

static void expect_constant_envelope_limits() {
    ResourceLimits limits;
    limits.max_steps = 0;
    for (bool use_pi : {false, true}) {
        ComputationContext exhausted(limits);
        auto result = use_pi ? detail::exact_pi_bounds(128, exhausted)
                             : detail::exact_e_bounds(128, exhausted);
        EXPECT_TRUE((!result && result.error().code == CasErrc::ResourceLimit)) << "constant envelopes cannot return partial success on exhaustion";
    }
    ComputationContext oversized;
    auto too_precise = detail::exact_pi_bounds(std::numeric_limits<std::size_t>::max(), oversized);
    EXPECT_TRUE((!too_precise && too_precise.error().code == CasErrc::ResourceLimit)) << "precision arithmetic overflow is rejected before integer allocation";
    CancellationToken token;
    token.cancel();
    ComputationContext cancelled({}, token);
    auto stopped = detail::exact_atan_bounds(Rational(2), 128, cancelled);
    EXPECT_TRUE((!stopped && stopped.error().code == CasErrc::Cancelled)) << "atan transforms share the caller's cancellation state";
}

TEST(IrrationalRegressions, CertifiedConstantEnvelopes) {
    const Rational tolerance(BigInt(1), BigInt(1) << 128);
    const Rational pi_lower("3.14159265358979323846264338327950288419716939937510");
    const Rational pi_upper("3.14159265358979323846264338327950288419716939937511");
    const Rational e_lower("2.71828182845904523536028747135266249775724709369995");
    const Rational e_upper("2.71828182845904523536028747135266249775724709369996");
    ComputationContext context;
    auto pi = detail::exact_pi_bounds(128, context);
    auto e = detail::exact_e_bounds(128, context);
    EXPECT_TRUE((pi && e)) << "certified constant construction succeeds";
    if (pi && e) {
        EXPECT_TRUE((encloses_reference_with_width(pi.value(), pi_lower, pi_upper, tolerance))) << "pi contains independent decimal brackets and meets the width contract";
        EXPECT_TRUE((encloses_reference_with_width(e.value(), e_lower, e_upper, tolerance))) << "e contains independent decimal brackets and meets the width contract";
    }
    expect_atan_series_envelopes(context, tolerance);
    expect_atan_transformed_envelopes(context, tolerance, pi_lower, pi_upper);
    expect_constant_envelope_limits();
}
