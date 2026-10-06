#include "expr.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

using namespace LMCAS;

namespace {

TEST(ExprComplexEvaluation, RealPromotion) {
    auto real_as_complex = LMCAS::eval_complex(*SymbolicExpr::number(5));
    EXPECT_TRUE((real_as_complex && real_as_complex.value().is_finite())) << "eval_complex accepts real expressions explicitly";
    EXPECT_NEAR(real_as_complex.value().real.value, 5.0, 0.0) << "eval_complex preserves real component";
    EXPECT_NEAR(real_as_complex.value().imag.value, 0.0, 0.0) << "eval_complex promotes real expression with zero imaginary component";
}

TEST(ExprComplexEvaluation, ImaginaryIdentifierEvaluation) {
    auto lowercase_i_complex =
        LMCAS::eval_complex(*SymbolicExpr::variable("i"));
    auto upper_i_complex =
        LMCAS::eval_complex(*SymbolicExpr::variable("I"));
    EXPECT_TRUE((!lowercase_i_complex &&
                 lowercase_i_complex.error().code == LMCAS::CasErrc::UnboundSymbol));
    EXPECT_TRUE((!upper_i_complex &&
                 upper_i_complex.error().code == LMCAS::CasErrc::UnboundSymbol));

    auto imaginary = LMCAS::imaginary_unit();
    ASSERT_TRUE(imaginary);
    auto evaluated = LMCAS::eval_complex(*imaginary.value());
    ASSERT_TRUE(evaluated);
    EXPECT_TRUE(evaluated.value().is_finite());
    EXPECT_NEAR(evaluated.value().imag.value, 1.0, 0.0)
        << "the explicit imaginary unit has imaginary component one";
}

TEST(ExprComplexEvaluation, ComplexSum) {
    auto four_i = LMCAS::complex(SymbolicExpr::number(0),
                                 SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3),
                                               four_i.value());
    auto lowered_complex = LMCAS::eval_complex(*three_plus_four_i);
    EXPECT_TRUE((lowered_complex && lowered_complex.value().is_finite())) << "eval_complex lowers 3 + 4I";
    EXPECT_NEAR(lowered_complex.value().real.value, 3.0, 0.0) << "eval_complex computes real part of 3 + 4I";
    EXPECT_NEAR(lowered_complex.value().imag.value, 4.0, 0.0) << "eval_complex computes imaginary part of 3 + 4I";
}

TEST(ExprComplexEvaluation, CertifiedComplexComponents) {
    auto atom = LMCAS::complex(SymbolicExpr::number(1.5),
                               SymbolicExpr::number(-0.5));
    ASSERT_TRUE(atom);
    auto evaluated = LMCAS::eval_complex(*atom.value());
    ASSERT_TRUE(evaluated);
    EXPECT_TRUE(evaluated.value().is_finite());
    EXPECT_DOUBLE_EQ(evaluated.value().real.value, 1.5);
    EXPECT_DOUBLE_EQ(evaluated.value().imag.value, -0.5);
    EXPECT_DOUBLE_EQ(evaluated.value().real.absolute_error, 0.0);
    EXPECT_DOUBLE_EQ(evaluated.value().imag.absolute_error, 0.0);

    auto mixed = LMCAS::complex(
        SymbolicExpr::number(Rational(1, 3)),
        SymbolicExpr::number(-0.5));
    ASSERT_TRUE(mixed);
    auto mixed_value = LMCAS::eval_complex(*mixed.value());
    ASSERT_TRUE(mixed_value);
    EXPECT_TRUE(mixed_value.value().is_finite());
    EXPECT_EQ(mixed_value.value().real.absolute_error,
              std::numeric_limits<double>::infinity());
    EXPECT_DOUBLE_EQ(mixed_value.value().imag.absolute_error, 0.0);
}

TEST(ExprComplexEvaluation, ComputedComplexErrorUnknown) {
    auto lhs = LMCAS::complex(SymbolicExpr::number(1.5),
                              SymbolicExpr::number(2.5));
    auto rhs = LMCAS::complex(SymbolicExpr::number(0.5),
                              SymbolicExpr::number(-0.5));
    ASSERT_TRUE(lhs);
    ASSERT_TRUE(rhs);

    auto sum = LMCAS::eval_complex(*SymbolicExpr::add(lhs.value(), rhs.value()));
    ASSERT_TRUE(sum);
    EXPECT_TRUE(sum.value().is_finite());
    EXPECT_DOUBLE_EQ(sum.value().real.value, 2.0);
    EXPECT_DOUBLE_EQ(sum.value().imag.value, 2.0);
    EXPECT_EQ(sum.value().real.absolute_error,
              std::numeric_limits<double>::infinity());
    EXPECT_EQ(sum.value().imag.absolute_error,
              std::numeric_limits<double>::infinity());

    auto product = LMCAS::eval_complex(
        *SymbolicExpr::multiply(lhs.value(), rhs.value()));
    ASSERT_TRUE(product);
    EXPECT_TRUE(product.value().is_finite());
    EXPECT_DOUBLE_EQ(product.value().real.value, 2.0);
    EXPECT_DOUBLE_EQ(product.value().imag.value, 0.5);
    EXPECT_EQ(product.value().real.absolute_error,
              std::numeric_limits<double>::infinity());
    EXPECT_EQ(product.value().imag.absolute_error,
              std::numeric_limits<double>::infinity());

    auto quotient = LMCAS::eval_complex(
        *SymbolicExpr::multiply(
            lhs.value(),
            SymbolicExpr::power(rhs.value(), SymbolicExpr::number(-1))));
    ASSERT_TRUE(quotient);
    EXPECT_TRUE(quotient.value().is_finite());
    EXPECT_DOUBLE_EQ(quotient.value().real.value, -1.0);
    EXPECT_DOUBLE_EQ(quotient.value().imag.value, 4.0);
    EXPECT_EQ(quotient.value().real.absolute_error,
              std::numeric_limits<double>::infinity());
    EXPECT_EQ(quotient.value().imag.absolute_error,
              std::numeric_limits<double>::infinity());

    auto zero_power = LMCAS::eval_complex(
        *SymbolicExpr::power(lhs.value(), SymbolicExpr::number(0)));
    ASSERT_TRUE(zero_power);
    EXPECT_DOUBLE_EQ(zero_power.value().real.value, 1.0);
    EXPECT_DOUBLE_EQ(zero_power.value().imag.value, 0.0);
    EXPECT_EQ(zero_power.value().real.absolute_error,
              std::numeric_limits<double>::infinity());
    EXPECT_EQ(zero_power.value().imag.absolute_error,
              std::numeric_limits<double>::infinity());
}

TEST(ExprComplexEvaluation, ComplexArithmeticOverflow) {
    auto large = LMCAS::complex(
        SymbolicExpr::number(std::numeric_limits<double>::max()),
        SymbolicExpr::number(0.0));
    ASSERT_TRUE(large);
    auto sum = LMCAS::eval_complex(
        *SymbolicExpr::add(large.value(), large.value()));
    ASSERT_FALSE(sum);
    EXPECT_EQ(sum.error().code, CasErrc::NumericFailure);
}

TEST(ExprComplexEvaluation, OrdinaryComplexArithmetic) {
    auto i = LMCAS::imaginary_unit();
    auto ordinary_multiply_complex = SymbolicExpr::add(
        SymbolicExpr::number(3),
        SymbolicExpr::multiply(SymbolicExpr::number(4), i.value()));
    auto lowered_ordinary_multiply =
        LMCAS::eval_complex(*ordinary_multiply_complex);
    EXPECT_TRUE((lowered_ordinary_multiply &&
                 lowered_ordinary_multiply.value().is_finite()))
        << "eval_complex lowers the LMCAS 3 + 4 * I ordinary multiplication form";
    EXPECT_NEAR(lowered_ordinary_multiply.value().real.value, 3.0, 0.0) << "ordinary multiplication complex form preserves real part";
    EXPECT_NEAR(lowered_ordinary_multiply.value().imag.value, 4.0, 0.0) << "ordinary multiplication complex form preserves imaginary part";

    auto i_power_two = SymbolicExpr::power(i.value(), SymbolicExpr::number(2));
    auto lowered_i_squared = LMCAS::eval_complex(*i_power_two);
    EXPECT_TRUE((lowered_i_squared && lowered_i_squared.value().is_finite())) << "eval_complex supports the LMCAS I^2 rule";
    EXPECT_NEAR(lowered_i_squared.value().real.value, -1.0, 0.0) << "eval_complex computes I^2 real part";
    EXPECT_NEAR(lowered_i_squared.value().imag.value, 0.0, 0.0) << "eval_complex computes I^2 imaginary part";
}

TEST(ExprComplexEvaluation, IntegerPowerBounds) {
    auto i = LMCAS::imaginary_unit();
    for (int exponent : {-65, -64, 64, 65}) {
        auto boundary_power = LMCAS::eval_complex(
            *SymbolicExpr::power(i.value(), SymbolicExpr::number(exponent)));
        if (exponent == -64 || exponent == 64) {
            EXPECT_TRUE((boundary_power &&
                         boundary_power.value().real.value == 1.0 &&
                         boundary_power.value().imag.value == 0.0))
                << "complex integer powers preserve the unit cycle at both supported bounds";
        } else {
            EXPECT_TRUE((!boundary_power &&
                         boundary_power.error().code ==
                             LMCAS::CasErrc::UnsupportedExpression))
                << "complex powers outside the LMCAS integer range remain unsupported";
        }
    }
}

TEST(ExprComplexEvaluation, ScaledMultiplication) {
    auto scaled_lhs = LMCAS::complex(
        SymbolicExpr::number(0x1.8p1023), SymbolicExpr::number(0x1p1022));
    auto scaled_rhs = LMCAS::complex(
        SymbolicExpr::number(1.375), SymbolicExpr::number(0.5));
    auto scaled_product = LMCAS::eval_complex(
        *SymbolicExpr::multiply(scaled_lhs.value(), scaled_rhs.value()));
    EXPECT_TRUE((scaled_product &&
                 scaled_product.value().real.value == 0x1.dp1023 &&
                 scaled_product.value().imag.value == 0x1.7p1023))
        << "complex multiplication retains finite components at large scales";
}

TEST(ExprComplexEvaluation, SubnormalMultiplication) {
    auto tiny_lhs = LMCAS::complex(
        SymbolicExpr::number(0x1p-1074), SymbolicExpr::number(0x1p-1074));
    auto tiny_rhs = LMCAS::complex(
        SymbolicExpr::number(0.5), SymbolicExpr::number(0.5));
    auto tiny_product = LMCAS::eval_complex(
        *SymbolicExpr::multiply(tiny_lhs.value(), tiny_rhs.value()));
    EXPECT_TRUE((tiny_product && tiny_product.value().real.value == 0.0 &&
                 tiny_product.value().imag.value == 0x1p-1074))
        << "complex multiplication combines subnormal contributions";
}

TEST(ExprComplexEvaluation, CancellingMultiplication) {
    auto cancelling_lhs = LMCAS::complex(
        SymbolicExpr::number(1.0 + 0x1p-52), SymbolicExpr::number(1.0));
    auto cancelling_rhs = LMCAS::complex(
        SymbolicExpr::number(1.0 - 0x1p-52), SymbolicExpr::number(1.0));
    auto cancelling_product = LMCAS::eval_complex(
        *SymbolicExpr::multiply(cancelling_lhs.value(), cancelling_rhs.value()));
    EXPECT_TRUE((cancelling_product &&
                 cancelling_product.value().real.value == -0x1p-104 &&
                 cancelling_product.value().imag.value == 2.0))
        << "complex multiplication preserves the nonzero remainder of cancelling products";
}

TEST(ExprComplexEvaluation, ExtremeReciprocal) {
    auto extreme_denominator = LMCAS::complex(
        SymbolicExpr::number(std::numeric_limits<double>::max()),
        SymbolicExpr::number(std::numeric_limits<double>::max()));
    ASSERT_TRUE((extreme_denominator.has_value())) << "eval_complex constructs an extreme finite denominator";
    if (extreme_denominator) {
        auto reciprocal_expression = SymbolicExpr::power(
            extreme_denominator.value(), SymbolicExpr::number(-1));
        auto extreme_reciprocal =
            LMCAS::eval_complex(*reciprocal_expression);
        EXPECT_TRUE((extreme_reciprocal &&
                     extreme_reciprocal.value().is_finite() &&
                     extreme_reciprocal.value().real.value ==
                         0.5 / std::numeric_limits<double>::max() &&
                     extreme_reciprocal.value().imag.value ==
                         -0.5 / std::numeric_limits<double>::max()))
            << "eval_complex avoids intermediate overflow in finite complex division";
    }
}

TEST(ExprComplexEvaluation, NegativePowerScale) {
    auto large_negative_power_base = LMCAS::complex(
        SymbolicExpr::number(1.0e160),
        SymbolicExpr::number(1.0e160));
    ASSERT_TRUE((large_negative_power_base.has_value())) << "eval_complex constructs a large finite negative-power base";
    if (large_negative_power_base) {
        auto negative_square_expression = SymbolicExpr::power(
            large_negative_power_base.value(), SymbolicExpr::number(-2));
        auto negative_square =
            LMCAS::eval_complex(*negative_square_expression);
        EXPECT_TRUE((negative_square &&
                     negative_square.value().is_finite()))
            << "negative complex powers avoid overflowing the positive power";
        if (negative_square) {
            EXPECT_NEAR(negative_square.value().real.value, 0.0, 0.0) << "the reciprocal square has a zero real component";
            EXPECT_TRUE((negative_square.value().imag.value < 0.0 &&
                         std::abs(
                             negative_square.value().imag.value /
                                 -5.0e-321 -
                             1.0) < 1.0e-3))
                << "the reciprocal square preserves its representable subnormal component";
        }
    }
}

TEST(ExprComplexEvaluation, UnboundComplexSymbols) {
    auto x = LMCAS::sym("x");
    auto linear = SymbolicExpr::add(x.value(), SymbolicExpr::number(2));
    auto complex_unbound = LMCAS::eval_complex(*linear);
    EXPECT_TRUE((!complex_unbound &&
                 complex_unbound.error().code == LMCAS::CasErrc::UnboundSymbol))
        << "eval_complex rejects unbound symbols during Expr to complex lowering";
    EXPECT_TRUE((!complex_unbound &&
                 std::string(LMCAS::error_name(complex_unbound.error())) ==
                     "ComplexEvalUnboundSymbol"))
        << "eval_complex exposes the LMCAS complex unbound-symbol diagnostic";
}

TEST(ExprComplexEvaluation, PrincipalComplexPowers) {
    auto i = LMCAS::imaginary_unit();
    ASSERT_TRUE(i);
    auto square_root = LMCAS::eval_complex(
        *SymbolicExpr::power(i.value(), SymbolicExpr::number(Rational(1, 2))));
    ASSERT_TRUE(square_root);
    EXPECT_NEAR(square_root.value().real.value, std::sqrt(0.5), 1e-14);
    EXPECT_NEAR(square_root.value().imag.value, std::sqrt(0.5), 1e-14);

    auto negative = LMCAS::complex(SymbolicExpr::number(-4), SymbolicExpr::number(0));
    ASSERT_TRUE(negative);
    auto branch_cut = LMCAS::eval_complex(
        *SymbolicExpr::power(negative.value(), SymbolicExpr::number(0.5)));
    ASSERT_TRUE(branch_cut);
    EXPECT_NEAR(branch_cut.value().real.value, 0.0, 1e-14);
    EXPECT_NEAR(branch_cut.value().imag.value, 2.0, 1e-14);

    auto negative_eight = LMCAS::complex(SymbolicExpr::number(-8),
                                         SymbolicExpr::number(0));
    ASSERT_TRUE(negative_eight);
    auto cube_root = LMCAS::eval_complex(
        *SymbolicExpr::power(negative_eight.value(),
                             SymbolicExpr::number(Rational(1, 3))));
    ASSERT_TRUE(cube_root);
    EXPECT_NEAR(cube_root.value().real.value, 1.0, 1e-14);
    EXPECT_NEAR(cube_root.value().imag.value, std::sqrt(3.0), 1e-14);

    auto imaginary_exponent = LMCAS::eval_complex(
        *SymbolicExpr::power(i.value(), i.value()));
    ASSERT_TRUE(imaginary_exponent);
    EXPECT_NEAR(imaginary_exponent.value().real.value,
                std::exp(-std::acos(-1.0) / 2.0), 1e-14);
    EXPECT_NEAR(imaginary_exponent.value().imag.value, 0.0, 1e-14);
}

TEST(ExprComplexEvaluation, ZeroFractionalPowers) {
    auto zero = LMCAS::complex(SymbolicExpr::number(0), SymbolicExpr::number(0));
    ASSERT_TRUE(zero);
    auto root = LMCAS::eval_complex(
        *SymbolicExpr::power(zero.value(), SymbolicExpr::number(0.5)));
    ASSERT_TRUE(root);
    EXPECT_DOUBLE_EQ(root.value().real.value, 0.0);
    EXPECT_DOUBLE_EQ(root.value().imag.value, 0.0);
    auto inverse_root = LMCAS::eval_complex(
        *SymbolicExpr::power(zero.value(), SymbolicExpr::number(-0.5)));
    ASSERT_FALSE(inverse_root);
    EXPECT_EQ(inverse_root.error().code, CasErrc::DomainError);
}

TEST(ExprComplexEvaluation, ExactFractionalExponentNearInteger) {
    auto i = LMCAS::imaginary_unit();
    ASSERT_TRUE(i);

    const BigInt exact_boundary(std::string("9007199254740992"));
    auto near_integer = SymbolicExpr::number(
        Rational(exact_boundary + BigInt(1), exact_boundary));
    auto near_integer_power = LMCAS::eval_complex(
        *SymbolicExpr::power(i.value(), near_integer));
    EXPECT_TRUE((!near_integer_power &&
                 near_integer_power.error().code ==
                     LMCAS::CasErrc::UnsupportedExpression))
        << "exact fractional exponents retain their rational domain";
}

TEST(ExprComplexEvaluation, ExactExponentArithmetic) {
    auto i = LMCAS::imaginary_unit();
    const BigInt exact_boundary(std::string("9007199254740992"));
    auto exact_difference = SymbolicExpr::add(
        SymbolicExpr::number(exact_boundary + BigInt(1)),
        SymbolicExpr::number(-exact_boundary));
    auto exact_difference_power = LMCAS::eval_complex(
        *SymbolicExpr::power(i.value(), exact_difference));
    EXPECT_TRUE((exact_difference_power &&
                 exact_difference_power.value().real.value == 0.0 &&
                 exact_difference_power.value().imag.value == 1.0))
        << "exact exponent arithmetic preserves I raised to an exact one";
}

TEST(ExprComplexEvaluation, NonfiniteComplexExponents) {
    auto i = LMCAS::imaginary_unit();
    auto nonfinite_power = SymbolicExpr::power(
        i.value(), SymbolicExpr::infinity());
    auto nonfinite_power_result =
        LMCAS::eval_complex(*nonfinite_power);
    EXPECT_TRUE((!nonfinite_power_result &&
                 nonfinite_power_result.error().code ==
                     LMCAS::CasErrc::NumericFailure))
        << "eval_complex rejects non-finite complex power exponents";
    EXPECT_TRUE((!nonfinite_power_result &&
                 std::string(LMCAS::error_name(
                     nonfinite_power_result.error())) ==
                     "NumericFailure"))
        << "eval_complex reports NumericFailure for non-finite exponents";
}

TEST(ExprComplexEvaluation, ZeroReciprocal) {
    auto zero_inverse = LMCAS::eval_complex(
        *SymbolicExpr::power(SymbolicExpr::number(0), SymbolicExpr::number(-1)));
    EXPECT_TRUE((!zero_inverse &&
                 zero_inverse.error().code == LMCAS::CasErrc::DomainError))
        << "eval_complex reports DomainError for complex reciprocal of zero";
}

TEST(ExprComplexEvaluation, ComplexEvaluationBudget) {
    auto four_i = LMCAS::complex(SymbolicExpr::number(0), SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3), four_i.value());
    LMCAS::ResourceLimits exhausted_complex_limits;
    exhausted_complex_limits.max_steps = 0;
    LMCAS::ComputationContext exhausted_complex_context(exhausted_complex_limits);
    auto exhausted_complex = LMCAS::eval_complex(
        *three_plus_four_i, {}, exhausted_complex_context);
    EXPECT_TRUE((!exhausted_complex &&
                 exhausted_complex.error().code ==
                     LMCAS::CasErrc::ResourceLimit))
        << "eval_complex reports ResourceLimit when the computation budget is exhausted";
    EXPECT_TRUE((!exhausted_complex &&
                 std::string(LMCAS::error_name(
                     exhausted_complex.error())) == "ResourceLimit"))
        << "eval_complex exposes the LMCAS resource-limit diagnostic";
}

TEST(ExprComplexEvaluation, RecursionBudgetRecovery) {
    auto four_i = LMCAS::complex(SymbolicExpr::number(0), SymbolicExpr::number(4));
    auto three_plus_four_i = SymbolicExpr::add(SymbolicExpr::number(3), four_i.value());
    LMCAS::ResourceLimits depth_limits;
    depth_limits.max_recursion_depth = 4;
    LMCAS::ComputationContext depth_context(depth_limits);
    auto imaginary = LMCAS::imaginary_unit();
    ASSERT_TRUE(imaginary);
    auto nested = imaginary.value();
    for (int depth = 0; depth < 8; ++depth) {
        nested = SymbolicExpr::power(nested, SymbolicExpr::number(1));
    }
    auto depth_limited = LMCAS::eval_complex(*nested, {}, depth_context);
    EXPECT_TRUE((!depth_limited &&
                 depth_limited.error().code == LMCAS::CasErrc::ResourceLimit))
        << "nested complex evaluation honors the configured depth budget";
    auto recovered = LMCAS::eval_complex(
        *three_plus_four_i, {}, depth_context);
    EXPECT_TRUE((recovered && recovered.value().real.value == 3.0 &&
                 recovered.value().imag.value == 4.0))
        << "the context supports shallow evaluation after a depth failure";
    auto repeated = LMCAS::eval_complex(
        *three_plus_four_i, {}, depth_context);
    EXPECT_TRUE((repeated && repeated.value().real.value == 3.0 &&
                 repeated.value().imag.value == 4.0))
        << "successful evaluation releases the active recursion frames";
}

} // namespace
