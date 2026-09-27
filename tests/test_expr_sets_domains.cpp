#include "expr.hpp"
#include <gtest/gtest.h>
#include <limits>

using namespace LMCAS;

namespace {

TEST(ExprSetsDomains, NumberDomainHierarchy) {
    auto domain_z = LMCAS::integers();
    auto domain_q = LMCAS::rationals();
    auto domain_r = LMCAS::reals();
    auto domain_c = LMCAS::complexes();
    auto domain_expr = LMCAS::expressions();
    bool names_match = std::string(domain_z.name()) == "Z";
    if (names_match) {
        names_match = std::string(domain_q.name()) == "Q";
    }
    if (names_match) {
        names_match = std::string(domain_r.name()) == "R";
    }
    if (names_match) {
        names_match = std::string(domain_c.name()) == "C";
    }
    if (names_match) {
        names_match = std::string(domain_expr.name()) == "Expr";
    }
    EXPECT_TRUE((names_match)) << "predefined number domain sets use LMCAS names";
    bool hierarchy_matches = domain_z.subset_of(domain_q);
    if (hierarchy_matches) {
        hierarchy_matches = domain_q.subset_of(domain_r);
    }
    if (hierarchy_matches) {
        hierarchy_matches = domain_r.subset_of(domain_c);
    }
    if (hierarchy_matches) {
        hierarchy_matches = domain_c.subset_of(domain_expr);
    }
    EXPECT_TRUE((hierarchy_matches)) << "predefined number domains follow Z subset Q subset R subset C subset Expr";
    auto facade_z_q = LMCAS::domain_subset(domain_z, domain_q);
    auto facade_c_r = LMCAS::domain_subset(domain_c, domain_r);
    auto facade_c_expr = LMCAS::domain_subset(domain_c, domain_expr);
    EXPECT_TRUE((facade_z_q && facade_z_q.value())) << "number domain subset facade accepts Z subset Q";
    EXPECT_TRUE((facade_c_r && !facade_c_r.value())) << "number domain subset facade rejects C subset R";
    EXPECT_TRUE((facade_c_expr && facade_c_expr.value())) << "number domain subset facade accepts C subset Expr";
}

TEST(ExprSetsDomains, ExactAndApproximateMembership) {
    auto domain_z = LMCAS::integers();
    auto domain_q = LMCAS::rationals();
    auto domain_r = LMCAS::reals();
    auto exact_two = LMCAS::integer(2);
    auto exact_half = LMCAS::rational(
        Rational(BigInt(1), BigInt(2)));
    auto domain_approx_half = LMCAS::approx_real(0.5);
    auto z_contains_two = LMCAS::domain_contains(domain_z,
                                                 exact_two.value());
    auto z_contains_half = LMCAS::domain_contains(domain_z,
                                                  exact_half.value());
    auto q_contains_half = LMCAS::domain_contains(domain_q,
                                                  exact_half.value());
    auto q_contains_approx = LMCAS::domain_contains(
        domain_q, domain_approx_half.value());
    auto r_contains_approx = LMCAS::domain_contains(
        domain_r, domain_approx_half.value());
    EXPECT_TRUE((z_contains_two && z_contains_two.value())) << "Z contains exact integer literals";
    EXPECT_TRUE((z_contains_half && !z_contains_half.value())) << "Z rejects non-integer rationals";
    EXPECT_TRUE((q_contains_half && q_contains_half.value())) << "Q contains exact rational literals";
    EXPECT_TRUE((q_contains_approx && !q_contains_approx.value())) << "Q does not claim approximate reals as exact rationals";
    EXPECT_TRUE((r_contains_approx && r_contains_approx.value())) << "R contains finite approximate real literals";
}

TEST(ExprSetsDomains, ArithmeticMembershipUsesTheExactNormalizedValue) {
    const auto half = SymbolicExpr::number(Rational(1, 2));
    const auto half_plus_half = SymbolicExpr::add(half, half);
    const auto integer_sum = domain_contains(integers(), half_plus_half);
    ASSERT_TRUE(integer_sum);
    EXPECT_TRUE(integer_sum.value())
        << "an exact sum belongs to Z when its normalized value is one";

    const auto root_two = SymbolicExpr::sqrt(SymbolicExpr::number(2));
    const auto root_product =
        SymbolicExpr::multiply(root_two, root_two);
    const auto rational_product =
        domain_contains(rationals(), root_product);
    ASSERT_TRUE(rational_product);
    EXPECT_TRUE(rational_product.value())
        << "an exactly reducible radical product belongs to Q";

    const auto sine_zero =
        domain_contains(integers(), SymbolicExpr::sin(SymbolicExpr::number(0)));
    ASSERT_TRUE(sine_zero);
    EXPECT_TRUE(sine_zero.value())
        << "exact function values are classified after normalization";
}

TEST(ExprSetsDomains, NegativeIntegerPowerMembershipDistinguishesUndefinedness) {
    const auto exponent = SymbolicExpr::number(-1);
    const auto zero_reciprocal = domain_contains(
        integers(),
        SymbolicExpr::power(SymbolicExpr::number(0), exponent));
    ASSERT_FALSE(zero_reciprocal);
    EXPECT_EQ(zero_reciprocal.error().code, CasErrc::DomainError)
        << "zero to a negative power is undefined";

    const auto negative_one_reciprocal = domain_contains(
        integers(),
        SymbolicExpr::power(SymbolicExpr::number(-1), exponent));
    ASSERT_TRUE(negative_one_reciprocal);
    EXPECT_TRUE(negative_one_reciprocal.value())
        << "negative one to a negative integer power remains an integer";
}

TEST(ExprSetsDomains, NonmemberOperandsDoNotDecideCompositeMembership) {
    const auto expression = SymbolicExpr::add(
        SymbolicExpr::number(Rational(1, 2)),
        SymbolicExpr::variable("domain_unknown_addend"));
    const auto membership = domain_contains(integers(), expression);
    ASSERT_FALSE(membership);
    EXPECT_EQ(membership.error().code, CasErrc::Inconclusive)
        << "a noninteger addend cannot prove that the complete sum is noninteger";
}

TEST(ExprSetsDomains, ImmutableMemberDomains) {
    auto source = rational(Rational(BigInt(1), BigInt(2))).value();
    auto set = expr_set({source});
    ASSERT_TRUE((set.has_value())) << "finite rational singleton constructs";
    if (!set)
        return;
    *source = *approx_real(0.5).value();
    const auto &member = set.value().elements()[0];
    auto exact_membership = rationals().contains(member);
    auto integer_membership = integers().contains(member);
    auto changed_membership = domain_contains(rationals(), source);
    EXPECT_TRUE((exact_membership && exact_membership.value() &&
                 integer_membership && !integer_membership.value()))
        << "read-only member retains its original exact rational domain";
    EXPECT_TRUE((changed_membership && !changed_membership.value())) << "mutable source independently changes to an approximate value";
    auto subset = expr_set_subset_domain(set.value(), rationals());
    EXPECT_TRUE((subset && subset.value())) << "domain subset consumes the same immutable member projection";
    ConstExprPtr null_member;
    auto invalid = reals().contains(null_member);
    EXPECT_TRUE((!invalid && invalid.error().code == CasErrc::InvalidArgument)) << "null immutable domain members are rejected";
}

TEST(ExprSetsDomains, ComplexDomainMembership) {
    auto domain_r = LMCAS::reals();
    auto domain_c = LMCAS::complexes();
    auto explicit_i_for_domain = LMCAS::imaginary_unit();
    auto r_contains_i = LMCAS::domain_contains(
        domain_r, explicit_i_for_domain.value());
    auto c_contains_i = LMCAS::domain_contains(
        domain_c, explicit_i_for_domain.value());
    auto r_contains_legacy_i = LMCAS::domain_contains(
        domain_r, SymbolicExpr::variable("i"));
    auto c_contains_legacy_i = LMCAS::domain_contains(
        domain_c, SymbolicExpr::variable("i"));
    auto c_contains_legacy_upper_i = LMCAS::domain_contains(
        domain_c, SymbolicExpr::variable("I"));
    EXPECT_TRUE((r_contains_i && !r_contains_i.value())) << "R rejects explicit non-real complex values";
    EXPECT_TRUE((c_contains_i && c_contains_i.value())) << "C contains explicit complex values";
    EXPECT_TRUE((!r_contains_legacy_i &&
                 r_contains_legacy_i.error().code ==
                     LMCAS::CasErrc::Inconclusive))
        << "R membership for ordinary i is undecidable without assumptions";
    EXPECT_TRUE((!c_contains_legacy_i &&
                 c_contains_legacy_i.error().code ==
                     LMCAS::CasErrc::Inconclusive &&
                 c_contains_legacy_upper_i &&
                 c_contains_legacy_upper_i.value()))
        << "ordinary i is undecidable while reserved I belongs to C";
}

TEST(ExprSetsDomains, SymbolicArithmeticMembership) {
    auto domain_c = LMCAS::complexes();
    auto domain_expr = LMCAS::expressions();
    auto legacy_four_i = SymbolicExpr::multiply(
        SymbolicExpr::number(4), SymbolicExpr::variable("i"));
    auto legacy_three_plus_four_i = SymbolicExpr::add(
        SymbolicExpr::number(3), legacy_four_i);
    auto c_contains_legacy_four_i =
        LMCAS::domain_contains(domain_c, legacy_four_i);
    auto c_contains_legacy_three_plus_four_i =
        LMCAS::domain_contains(domain_c, legacy_three_plus_four_i);
    auto expr_contains_symbol = LMCAS::domain_contains(
        domain_expr, LMCAS::sym("domain_expr_symbol").value());
    EXPECT_TRUE((!c_contains_legacy_four_i &&
                 c_contains_legacy_four_i.error().code ==
                     LMCAS::CasErrc::Inconclusive &&
                 !c_contains_legacy_three_plus_four_i &&
                 c_contains_legacy_three_plus_four_i.error().code ==
                     LMCAS::CasErrc::Inconclusive))
        << "arithmetic expressions built from ordinary i remain undecidable";
    EXPECT_TRUE((expr_contains_symbol && expr_contains_symbol.value())) << "Expr contains symbolic expressions";
}

TEST(ExprSetsDomains, NumericAndComplexSubsets) {
    auto domain_r = LMCAS::reals();
    auto domain_c = LMCAS::complexes();
    auto exact_two = LMCAS::integer(2);
    auto exact_half = LMCAS::rational(
        Rational(BigInt(1), BigInt(2)));
    auto domain_approx_half = LMCAS::approx_real(0.5);
    auto explicit_i_for_domain = LMCAS::imaginary_unit();
    auto numeric_domain_set = LMCAS::expr_set({exact_two.value(), exact_half.value(), domain_approx_half.value()});
    auto complex_domain_set = LMCAS::expr_set({exact_two.value(), explicit_i_for_domain.value()});
    auto numeric_subset_r =
        numeric_domain_set ? LMCAS::expr_set_subset_domain(
                                 numeric_domain_set.value(), domain_r)
                           : LMCAS::Result<bool>::failure(
                                 LMCAS::CasErrc::InternalInvariant,
                                 "numeric domain set construction failed",
                                 "test_expr");
    auto complex_subset_r =
        complex_domain_set ? LMCAS::expr_set_subset_domain(
                                 complex_domain_set.value(), domain_r)
                           : LMCAS::Result<bool>::failure(
                                 LMCAS::CasErrc::InternalInvariant,
                                 "complex domain set construction failed",
                                 "test_expr");
    auto complex_subset_c =
        complex_domain_set ? LMCAS::expr_set_subset_domain(
                                 complex_domain_set.value(), domain_c)
                           : LMCAS::Result<bool>::failure(
                                 LMCAS::CasErrc::InternalInvariant,
                                 "complex domain set construction failed",
                                 "test_expr");
    EXPECT_TRUE((numeric_subset_r && numeric_subset_r.value())) << "set<Expr> subset facade accepts numeric real sets";
    EXPECT_TRUE((complex_subset_r && !complex_subset_r.value())) << "set<Expr> subset facade rejects non-real complex members";
    EXPECT_TRUE((complex_subset_c && complex_subset_c.value())) << "set<Expr> subset facade accepts explicit complex members in C";
}

TEST(ExprSetsDomains, UndecidableComplexSubsets) {
    auto domain_c = LMCAS::complexes();
    auto legacy_four_i = SymbolicExpr::multiply(
        SymbolicExpr::number(4), SymbolicExpr::variable("i"));
    auto legacy_three_plus_four_i = SymbolicExpr::add(
        SymbolicExpr::number(3), legacy_four_i);
    auto legacy_complex_domain_set = LMCAS::expr_set({SymbolicExpr::variable("i")});
    auto legacy_complex_arithmetic_domain_set = LMCAS::expr_set({legacy_three_plus_four_i});
    auto legacy_complex_subset_c =
        legacy_complex_domain_set ? LMCAS::expr_set_subset_domain(
                                        legacy_complex_domain_set.value(),
                                        domain_c)
                                  : LMCAS::Result<bool>::failure(
                                        LMCAS::CasErrc::InternalInvariant,
                                        "legacy complex domain set construction failed",
                                        "test_expr");
    auto legacy_complex_arithmetic_subset_c =
        legacy_complex_arithmetic_domain_set
            ? LMCAS::expr_set_subset_domain(
                  legacy_complex_arithmetic_domain_set.value(), domain_c)
            : LMCAS::Result<bool>::failure(
                  LMCAS::CasErrc::InternalInvariant,
                  "legacy complex arithmetic domain set construction failed",
                  "test_expr");
    EXPECT_TRUE((!legacy_complex_subset_c &&
                 legacy_complex_subset_c.error().code ==
                     LMCAS::CasErrc::Inconclusive))
        << "set<Expr> subset facade keeps ordinary i membership undecidable";
    EXPECT_TRUE((!legacy_complex_arithmetic_subset_c &&
                 legacy_complex_arithmetic_subset_c.error().code ==
                     LMCAS::CasErrc::Inconclusive))
        << "set<Expr> subset facade keeps ordinary i arithmetic undecidable";
}

TEST(ExprSetsDomains, SymbolicExpressionSubsets) {
    auto domain_r = LMCAS::reals();
    auto domain_expr = LMCAS::expressions();
    auto unknown_domain_set = LMCAS::expr_set({LMCAS::sym("domain_set_unknown").value()});
    auto unknown_subset_r =
        unknown_domain_set ? LMCAS::expr_set_subset_domain(
                                 unknown_domain_set.value(), domain_r)
                           : LMCAS::Result<bool>::failure(
                                 LMCAS::CasErrc::InternalInvariant,
                                 "unknown domain set construction failed",
                                 "test_expr");
    auto unknown_subset_expr =
        unknown_domain_set ? LMCAS::expr_set_subset_domain(
                                 unknown_domain_set.value(), domain_expr)
                           : LMCAS::Result<bool>::failure(
                                 LMCAS::CasErrc::InternalInvariant,
                                 "unknown domain set construction failed",
                                 "test_expr");
    EXPECT_TRUE((!unknown_subset_r &&
                 unknown_subset_r.error().code ==
                     LMCAS::CasErrc::Inconclusive))
        << "set<Expr> subset facade propagates undecidable domain membership";
    EXPECT_TRUE((unknown_subset_expr && unknown_subset_expr.value())) << "set<Expr> subset facade accepts arbitrary Expr members in Expr";
}

TEST(ExprSetsDomains, InvalidAndUndecidableMembership) {
    auto domain_r = LMCAS::reals();
    auto unknown_domain_member = LMCAS::domain_contains(
        domain_r, LMCAS::sym("domain_unknown").value());
    auto null_domain_member = LMCAS::domain_contains(domain_r, nullptr);
    EXPECT_TRUE((!unknown_domain_member &&
                 unknown_domain_member.error().code ==
                     LMCAS::CasErrc::Inconclusive))
        << "number domain membership does not guess symbolic variables";
    EXPECT_TRUE((!null_domain_member &&
                 null_domain_member.error().code ==
                     LMCAS::CasErrc::InvalidArgument))
        << "number domain membership rejects null Expr values";
}

} // namespace
