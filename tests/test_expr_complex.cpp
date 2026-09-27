#include "expr.hpp"
#include <gtest/gtest.h>
#include <limits>

using namespace LMCAS;

namespace {

TEST(ExprComplex, ImaginaryUnitIdentity) {
    auto i = LMCAS::imaginary_unit();
    EXPECT_TRUE((i.has_value())) << "imaginary unit can be constructed";
    auto explicit_i = LMCAS::complex(SymbolicExpr::number(0),
                                     SymbolicExpr::number(1));
    EXPECT_TRUE((explicit_i.has_value())) << "complex(0, 1) can be constructed";
    EXPECT_TRUE((i && explicit_i &&
                 LMCAS::structurally_equal(*i.value(),
                                           *explicit_i.value())))
        << "imaginary unit is structurally complex(0, 1)";
    auto variable_i = SymbolicExpr::variable("i");
    EXPECT_TRUE((i && !LMCAS::structurally_equal(*i.value(), *variable_i))) << "ordinary variable(\"i\") is not structurally the imaginary unit";
}

TEST(ExprComplex, InvalidComplexComponents) {
    auto null_real_complex =
        LMCAS::complex(nullptr, SymbolicExpr::number(1));
    auto null_imag_complex =
        LMCAS::complex(SymbolicExpr::number(0), nullptr);
    EXPECT_TRUE((!null_real_complex &&
                 null_real_complex.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "complex(nullptr, 1) rejects an incompatible real part";
    EXPECT_TRUE((!null_imag_complex &&
                 null_imag_complex.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "complex(0, nullptr) rejects an incompatible imaginary part";
    EXPECT_TRUE((!null_real_complex &&
                 std::string(LMCAS::error_name(
                     null_real_complex.error())) == "ComplexTypeMismatch"))
        << "complex(nullptr, 1) exposes ComplexTypeMismatch";
    EXPECT_TRUE((!null_imag_complex &&
                 std::string(LMCAS::error_name(
                     null_imag_complex.error())) == "ComplexTypeMismatch"))
        << "complex(0, nullptr) exposes ComplexTypeMismatch";
}

TEST(ExprComplex, ImaginaryUnitMultiplication) {
    auto i = LMCAS::imaginary_unit();
    auto i_squared = SymbolicExpr::multiply(i.value(), i.value());
    LMCAS::ComputationContext complex_equivalence_context;
    auto complex_equivalent = LMCAS::equivalent_core(
        *i_squared, *SymbolicExpr::number(-1), complex_equivalence_context);
    EXPECT_TRUE((complex_equivalent && complex_equivalent.value())) << "I * I is equivalent to -1 in the LMCAS core profile";

    auto i_power_two = SymbolicExpr::power(i.value(), SymbolicExpr::number(2));
    LMCAS::ComputationContext complex_power_equivalence_context;
    auto complex_power_equivalent = LMCAS::equivalent_core(
        *i_power_two, *SymbolicExpr::number(-1),
        complex_power_equivalence_context);
    EXPECT_TRUE((complex_power_equivalent &&
                 complex_power_equivalent.value()))
        << "I^2 is equivalent to -1 in the LMCAS core profile";
}

TEST(ExprComplex, OrdinaryIRemainsSymbolic) {
    auto i = LMCAS::imaginary_unit();
    auto variable_i = SymbolicExpr::variable("i");
    auto legacy_i_squared =
        SymbolicExpr::multiply(variable_i, variable_i);
    LMCAS::ComputationContext legacy_i_equivalence_context;
    auto legacy_i_equivalent = LMCAS::equivalent_core(
        *legacy_i_squared, *SymbolicExpr::number(-1),
        legacy_i_equivalence_context);
    EXPECT_TRUE((legacy_i_equivalent && !legacy_i_equivalent.value())) << "ordinary variable(\"i\") does not follow the imaginary-unit multiplication rule";

    auto legacy_i_plus_one = SymbolicExpr::add(variable_i,
                                               SymbolicExpr::number(1));
    auto one_plus_canonical_i =
        SymbolicExpr::add(SymbolicExpr::number(1), i.value());
    LMCAS::ComputationContext legacy_i_add_context;
    auto legacy_i_add_equivalent = LMCAS::equivalent_core(
        *legacy_i_plus_one, *one_plus_canonical_i,
        legacy_i_add_context);
    EXPECT_TRUE((legacy_i_add_equivalent &&
                 !legacy_i_add_equivalent.value()))
        << "ordinary variable(\"i\") remains distinct inside additive equivalence";
}

TEST(ExprComplex, ComplexParts) {
    auto four_i = LMCAS::complex(SymbolicExpr::number(0), SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3), four_i.value());
    auto real_part = LMCAS::real(three_plus_four_i);
    EXPECT_TRUE((real_part && LMCAS::structurally_equal(
                                  *real_part.value(), *SymbolicExpr::number(3))))
        << "real(3 + 4I) returns 3";

    auto imag_part = LMCAS::imag(three_plus_four_i);
    EXPECT_TRUE((imag_part && LMCAS::structurally_equal(
                                  *imag_part.value(), *SymbolicExpr::number(4))))
        << "imag(3 + 4I) returns 4";

    auto conjugated = LMCAS::conj(three_plus_four_i);
    EXPECT_TRUE((conjugated.has_value())) << "conj(3 + 4I) succeeds";
    auto expected_conj = LMCAS::complex(SymbolicExpr::number(3),
                                        SymbolicExpr::number(-4));
    LMCAS::ComputationContext conj_context;
    auto conj_equiv = LMCAS::equivalent_core(
        *conjugated.value(), *expected_conj.value(), conj_context);
    EXPECT_TRUE((conj_equiv && conj_equiv.value())) << "conj(3 + 4I) returns 3 - 4I";
}

TEST(ExprComplex, ExactComplexModulus) {
    auto four_i = LMCAS::complex(SymbolicExpr::number(0), SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3), four_i.value());
    auto complex_abs = LMCAS::abs(three_plus_four_i);
    EXPECT_TRUE((complex_abs.has_value())) << "abs(3 + 4I) succeeds";
    auto abs_value = LMCAS::evalf(*complex_abs.value());
    EXPECT_TRUE((abs_value && abs_value.value().is_finite())) << "abs(3 + 4I) can be explicitly numerically evaluated";
    EXPECT_NEAR(abs_value.value().value, 5.0, 1e-12) << "abs(3 + 4I) evaluates to 5";
    EXPECT_TRUE((complex_abs && LMCAS::structurally_equal(
                                    *complex_abs.value(), *SymbolicExpr::number(5))))
        << "exact complex modulus remains in the exact number domain";
}

TEST(ExprComplex, ExtremeConstantModulus) {
    for (double component : {1e200, 1e-200}) {
        auto z = LMCAS::complex(
            SymbolicExpr::number(component), SymbolicExpr::number(component));
        auto magnitude = LMCAS::abs(z.value());
        ASSERT_TRUE((magnitude.has_value())) << "finite extreme modulus can be constructed";
        if (magnitude) {
            auto value = LMCAS::evalf(*magnitude.value());
            EXPECT_TRUE((value && value.value().is_finite() &&
                         std::abs(value.value().value / std::hypot(component, component) - 1) < 1e-14))
                << "constant complex modulus preserves its finite nonzero scale";
        }
    }
}

TEST(ExprComplex, ExtremeSymbolicModulus) {
    auto norm_variable = SymbolicExpr::variable("norm_component");
    auto bound_complex = LMCAS::complex(norm_variable, norm_variable);
    auto bound_modulus = LMCAS::abs(bound_complex.value());
    ASSERT_TRUE((bound_modulus.has_value())) << "symbolic modulus can be constructed";
    if (bound_modulus) {
        for (double component : {1e200, 1e-200}) {
            auto value = LMCAS::evalf(
                *bound_modulus.value(), {{"norm_component", component}});
            EXPECT_TRUE((value && value.value().is_finite() &&
                         std::abs(value.value().value / std::hypot(component, component) - 1) < 1e-14))
                << "bound complex modulus evaluates its components before squaring";
        }
    }
}

TEST(ExprComplex, MixedModulus) {
    auto norm_variable = SymbolicExpr::variable("norm_component");
    auto mixed_complex = LMCAS::complex(SymbolicExpr::number(1e200), norm_variable);
    auto mixed_modulus = LMCAS::abs(mixed_complex.value());
    ASSERT_TRUE((mixed_modulus.has_value())) << "mixed constant and symbolic modulus can be constructed";
    if (mixed_modulus) {
        auto value = LMCAS::evalf(*mixed_modulus.value(), {{"norm_component", 0.0}});
        EXPECT_TRUE((value && value.value().value == 1e200)) << "mixed modulus preserves an extreme constant before variable binding";
    }
}

TEST(ExprComplex, InvalidComplexPartInputs) {
    auto null_real = LMCAS::real(nullptr);
    auto null_imag = LMCAS::imag(nullptr);
    auto null_conj = LMCAS::conj(nullptr);
    auto null_abs = LMCAS::abs(nullptr);
    EXPECT_TRUE((!null_real &&
                 null_real.error().code == LMCAS::CasErrc::InvalidArgument))
        << "real(nullptr) rejects null LMCAS Expr";
    EXPECT_TRUE((!null_imag &&
                 null_imag.error().code == LMCAS::CasErrc::InvalidArgument))
        << "imag(nullptr) rejects null LMCAS Expr";
    EXPECT_TRUE((!null_conj &&
                 null_conj.error().code == LMCAS::CasErrc::InvalidArgument))
        << "conj(nullptr) rejects null LMCAS Expr";
    EXPECT_TRUE((!null_abs &&
                 null_abs.error().code == LMCAS::CasErrc::InvalidArgument))
        << "abs(nullptr) rejects null LMCAS Expr";
}

TEST(ExprComplex, RealValueParts) {
    auto real_number = SymbolicExpr::number(-5);
    auto real_number_part = LMCAS::real(real_number);
    EXPECT_TRUE((real_number_part &&
                 LMCAS::structurally_equal(*real_number_part.value(),
                                           *real_number)))
        << "real(-5) preserves a real value under R subset C";

    auto real_number_imag = LMCAS::imag(real_number);
    EXPECT_TRUE((real_number_imag &&
                 LMCAS::structurally_equal(*real_number_imag.value(),
                                           *SymbolicExpr::number(0))))
        << "imag(-5) returns zero for a real value under R subset C";

    auto real_number_conj = LMCAS::conj(real_number);
    EXPECT_TRUE((real_number_conj &&
                 LMCAS::structurally_equal(*real_number_conj.value(),
                                           *real_number)))
        << "conj(-5) preserves a real value under R subset C";

    auto real_number_abs = LMCAS::abs(real_number);
    auto real_number_abs_value =
        real_number_abs ? LMCAS::evalf(*real_number_abs.value())
                        : LMCAS::Result<LMCAS::ApproxReal>::failure(
                              LMCAS::CasErrc::InternalInvariant,
                              "abs(-5) construction failed", "test");
    EXPECT_TRUE((real_number_abs_value &&
                 real_number_abs_value.value().is_finite()))
        << "abs(-5) can be explicitly evaluated under R subset C";
    EXPECT_NEAR(real_number_abs_value.value().value, 5.0, 1e-12) << "abs(-5) evaluates to 5 under R subset C";
}

TEST(ExprComplex, UnsupportedComplexFunctionSplit) {
    auto four_i = LMCAS::complex(SymbolicExpr::number(0), SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3), four_i.value());
    auto unsupported_real = LMCAS::real(
        SymbolicExpr::sin(three_plus_four_i));
    EXPECT_TRUE((!unsupported_real &&
                 unsupported_real.error().code == LMCAS::CasErrc::Inconclusive))
        << "real(sin(3 + 4I)) reports unsupported complex function split";
}

void expect_exact_complex_components(
    const ExprPtr &expression, const Rational &expected_real, const Rational &expected_imag) {
    auto real_component = real(expression);
    auto imag_component = imag(expression);
    EXPECT_TRUE((real_component && imag_component)) << "complex components are extractable";
    if (!real_component || !imag_component) {
        return;
    }
    ComputationContext context;
    auto real_equal = equivalent_core(*real_component.value(),
                                      *SymbolicExpr::number(expected_real), context);
    auto imag_equal = equivalent_core(*imag_component.value(),
                                      *SymbolicExpr::number(expected_imag), context);
    EXPECT_TRUE((real_equal && real_equal.value())) << "real component equals exact oracle";
    EXPECT_TRUE((imag_equal && imag_equal.value())) << "imaginary component equals exact oracle";
}

TEST(ExprComplex, PrintedComplexComponents) {
    const auto number = [](int value) { return SymbolicExpr::number(value); };
    const auto sum = SymbolicExpr::add(number(3), number(4));
    struct Case {
        ExprPtr expression;
        Rational real;
        Rational imag;
    };
    auto basic = complex(number(2), sum).value();
    const Case cases[] = {
        {basic, Rational(2), Rational(7)},
        {complex(SymbolicExpr::add(number(2), number(3)), sum).value(),
         Rational(5), Rational(7)},
        {complex(number(2), SymbolicExpr::multiply(number(-1), sum)).value(),
         Rational(2), Rational(-7)},
        {complex(number(2), SymbolicExpr::number(Rational(3, 2))).value(),
         Rational(2), Rational(3, 2)},
        {complex(number(2), SymbolicExpr::multiply(sum, number(2))).value(),
         Rational(2), Rational(14)},
        {complex(number(2), SymbolicExpr::power(
                                SymbolicExpr::power(number(2), number(3)), number(2)))
             .value(),
         Rational(2), Rational(64)},
        {complex(number(2), SymbolicExpr::power(
                                number(2), SymbolicExpr::power(number(3), number(2))))
             .value(),
         Rational(2), Rational(512)},
        {SymbolicExpr::multiply(basic, number(3)), Rational(6), Rational(21)},
        {SymbolicExpr::power(basic, number(2)), Rational(-45), Rational(28)}};
    for (const auto &test : cases) {
        auto parsed = parse_expr(test.expression->to_string());
        ASSERT_TRUE((parsed.has_value())) << "complex expression text parses";
        if (!parsed) {
            continue;
        }
        for (const auto &expression : {test.expression, parsed.value()}) {
            expect_exact_complex_components(expression, test.real, test.imag);
        }
    }
    auto x = SymbolicExpr::variable("x");
    auto reciprocal = SymbolicExpr::power(SymbolicExpr::add(x, number(1)), number(-1));
    auto original = complex(number(2), reciprocal).value();
    auto parsed = parse_expr(original->to_string());
    ASSERT_TRUE((parsed.has_value())) << "complex reciprocal text parses";
    if (!parsed) {
        return;
    }
    for (const auto &expression : {original, parsed.value()}) {
        auto regular = eval_complex(*expression, {{"x", 1.0}});
        EXPECT_TRUE((regular && regular.value().real.value == 2.0 &&
                     regular.value().imag.value == 0.5))
            << "complex reciprocal preserves both components at regular points";
        auto singular = eval_complex(*expression, {{"x", -1.0}});
        EXPECT_TRUE((!singular && singular.error().code == CasErrc::DomainError)) << "printing cannot erase the imaginary component's singularity";
    }
}

} // namespace
