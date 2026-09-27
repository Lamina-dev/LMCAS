#include "test_common.hpp"
#include "transcendental_factor.hpp"
#include "poly_utils.hpp"

using namespace LMCAS;

namespace {

void expect_factor_set(
    const std::vector<std::shared_ptr<SymbolicExpr>> &actual,
    const std::vector<std::shared_ptr<SymbolicExpr>> &expected) {
    ASSERT_EQ(actual.size(), expected.size());
    std::vector<bool> matched(expected.size(), false);
    for (const auto &factor : actual) {
        bool found = false;
        for (std::size_t i = 0; i < expected.size(); ++i) {
            if (!matched[i] &&
                test_proved_equivalent(factor, expected[i])) {
                matched[i] = true;
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found);
    }
    for (bool was_matched : matched) {
        EXPECT_TRUE(was_matched);
    }
}

} // namespace

TEST(TranscendentalFactorBackSubstitution, BackSubstituteU0PlusXToSinXPlusX) {
    auto x = SymbolicExpr::variable("x");

    auto sin_x = SymbolicExpr::sin(x);

    auto original = SymbolicExpr::add(sin_x, x);
    auto sub_result = detect_trans_substitutions(original, "x");
    ASSERT_EQ(sub_result.mappings.size(), 1U);

    auto back = sub_result.poly_expr->substitute(
        sub_result.mappings[0].indeterminate,
        sub_result.mappings[0].trans_expr);

    ASSERT_NE(back, nullptr);
    EXPECT_TRUE(test_proved_equivalent(back, original));
}

TEST(TranscendentalFactorBackSubstitution, BackSubstituteU0SquaredMinusXSquared) {
    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);

    auto original = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
                               SymbolicExpr::power(x, SymbolicExpr::number(2))));

    auto sub_result = detect_trans_substitutions(original, "x");
    ASSERT_EQ(sub_result.mappings.size(), 1U);

    auto back = sub_result.poly_expr->substitute(
        sub_result.mappings[0].indeterminate,
        sub_result.mappings[0].trans_expr);

    ASSERT_NE(back, nullptr);
    EXPECT_TRUE(test_proved_equivalent(back, original));
}

TEST(TranscendentalFactorBackSubstitution, BackSubstituteMultipleMappings) {
    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);

    auto original = SymbolicExpr::add(sin_x, cos_x);

    auto sub_result = detect_trans_substitutions(original, "x");
    ASSERT_EQ(sub_result.mappings.size(), 2U);

    auto back = sub_result.poly_expr;
    for (const auto &mapping : sub_result.mappings) {
        back = back->substitute(
            mapping.indeterminate, mapping.trans_expr);
    }

    ASSERT_NE(back, nullptr);
    EXPECT_TRUE(test_proved_equivalent(back, original));
}

TEST(TranscendentalFactorBackSubstitution, BackSubstitutePolynomialFactor) {
    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);

    auto original = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        SymbolicExpr::add(sin_x, SymbolicExpr::number(1)));

    auto sub_result = detect_trans_substitutions(original, "x");
    ASSERT_EQ(sub_result.mappings.size(), 1U);

    Polynomial<Rational> factor_poly(
        {Rational(1), Rational(1)}, "u0");
    auto factor_expr = poly_to_symbolic(factor_poly);
    auto back = factor_expr->substitute(
        sub_result.mappings[0].indeterminate,
        sub_result.mappings[0].trans_expr);
    auto expected = SymbolicExpr::add(
        SymbolicExpr::sin(x), SymbolicExpr::number(1));

    ASSERT_NE(back, nullptr);
    EXPECT_TRUE(test_proved_equivalent(back, expected));
}

TEST(TranscendentalFactorBackSubstitution, BackSubstituteNoMappings) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(1));

    auto sub_result = detect_trans_substitutions(expr, "x");
    ASSERT_TRUE(sub_result.mappings.empty());

    auto back = sub_result.poly_expr;
    for (const auto &m : sub_result.mappings) {
        back = back->substitute(m.indeterminate, m.trans_expr);
    }

    ASSERT_NE(back, nullptr);
    EXPECT_TRUE(test_proved_equivalent(back, expr));
}

TEST(TranscendentalFactorBackSubstitution, BackSubstituteRoundtrip) {
    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);

    auto original = SymbolicExpr::add(
        SymbolicExpr::multiply(sin_x, cos_x),
        sin_x);

    auto sub_result = detect_trans_substitutions(original, "x");

    auto back = sub_result.poly_expr;
    for (const auto &m : sub_result.mappings) {
        back = back->substitute(m.indeterminate, m.trans_expr);
    }

    ASSERT_NE(back, nullptr);
    EXPECT_TRUE(test_proved_equivalent(back, original));
}

TEST(TranscendentalFactorBackSubstitution, SimplifyFactorsExtractConstantFromProduct) {
    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);

    auto factor1 = SymbolicExpr::multiply(SymbolicExpr::number(2), sin_x);
    auto factor2 = SymbolicExpr::add(x, SymbolicExpr::number(1));

    std::vector<std::shared_ptr<SymbolicExpr>> factors = {factor1, factor2};
    auto result = tf_simplify_factors(factors);

    expect_factor_set(
        result,
        {SymbolicExpr::number(2), sin_x, factor2});
}

TEST(TranscendentalFactorBackSubstitution, SimplifyFactorsNoConstants) {
    auto x = SymbolicExpr::variable("x");

    auto factor1 = SymbolicExpr::add(SymbolicExpr::sin(x), x);
    auto factor2 = SymbolicExpr::add(SymbolicExpr::cos(x), SymbolicExpr::number(-1));

    std::vector<std::shared_ptr<SymbolicExpr>> factors = {factor1, factor2};
    auto result = tf_simplify_factors(factors);

    expect_factor_set(result, {factor1, factor2});
}

TEST(TranscendentalFactorBackSubstitution, SimplifyFactorsPureConstantFactor) {
    auto x = SymbolicExpr::variable("x");

    auto factor1 = SymbolicExpr::number(3);
    auto factor2 = SymbolicExpr::add(x, SymbolicExpr::number(1));

    std::vector<std::shared_ptr<SymbolicExpr>> factors = {factor1, factor2};
    auto result = tf_simplify_factors(factors);

    expect_factor_set(result, {factor1, factor2});
}

TEST(TranscendentalFactorBackSubstitution, SimplifyFactorsMultipleConstantsCombined) {
    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);

    auto factor1 = SymbolicExpr::multiply(SymbolicExpr::number(2), sin_x);
    auto factor2 = SymbolicExpr::multiply(SymbolicExpr::number(5), cos_x);

    std::vector<std::shared_ptr<SymbolicExpr>> factors = {factor1, factor2};
    auto result = tf_simplify_factors(factors);

    expect_factor_set(
        result,
        {SymbolicExpr::number(10), sin_x, cos_x});
}

TEST(TranscendentalFactorBackSubstitution, SimplifyFactorsConstantOneNotAdded) {
    auto x = SymbolicExpr::variable("x");

    auto factor1 = SymbolicExpr::sin(x);
    auto factor2 = SymbolicExpr::add(x, SymbolicExpr::number(1));

    std::vector<std::shared_ptr<SymbolicExpr>> factors = {factor1, factor2};
    auto result = tf_simplify_factors(factors);

    expect_factor_set(result, {factor1, factor2});
}

TEST(TranscendentalFactorBackSubstitution, SimplifyFactorsSimplifyCalled) {
    auto x = SymbolicExpr::variable("x");

    auto factor1 = SymbolicExpr::add(x, SymbolicExpr::number(0));

    std::vector<std::shared_ptr<SymbolicExpr>> factors = {factor1};
    auto result = tf_simplify_factors(factors);

    expect_factor_set(result, {x});
}
