/**
 * @file test_series_sequence_limits.cpp
 * @brief 数列上极限与下极限契约测试。
 */

#include "test_common.hpp"
#include "series_engine.hpp"

using namespace LMCAS;

using Expr = std::shared_ptr<SymbolicExpr>;
using Coeffs = std::vector<Expr>;

static Expr num(int n) { return SymbolicExpr::number(n); }
static Expr var(const std::string &name) { return SymbolicExpr::variable(name); }
static Expr checked_sequence_limit(LMCAS::ExpressionResult result) {
    EXPECT_TRUE((result.has_value())) << "checked sequence limit succeeds";
    return result ? std::move(result.value()) : nullptr;
}

static void expect_inconclusive_sequence_limit(LMCAS::ExpressionResult result) {
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::Inconclusive);
}

TEST(SeriesSequenceLimits, LimSupConvergentSequence) {
    auto n = var("n");
    auto a_n = SymbolicExpr::power(n, num(-1));

    auto result = checked_sequence_limit(LMCAS::lim_sup_checked(a_n, "n"));
    ASSERT_TRUE(result);
    EXPECT_EQ(result->compare(num(0)), 0) << "lim_sup(1/n) = 0";
}

TEST(SeriesSequenceLimits, LimInfConvergentSequence) {
    auto n = var("n");
    auto a_n = SymbolicExpr::power(n, num(-1));

    auto result = checked_sequence_limit(LMCAS::lim_inf_checked(a_n, "n"));
    ASSERT_TRUE(result);
    EXPECT_EQ(result->compare(num(0)), 0) << "lim_inf(1/n) = 0";
}

TEST(SeriesSequenceLimits, LimSupAlternatingSequence) {
    auto n = var("n");
    auto a_n = SymbolicExpr::power(num(-1), n);

    auto result = checked_sequence_limit(LMCAS::lim_sup_checked(a_n, "n"));
    ASSERT_TRUE(result);
    EXPECT_EQ(result->compare(num(1)), 0)
        << "lim_sup((-1)^n) is exactly one";
}

TEST(SeriesSequenceLimits, LimInfAlternatingSequence) {
    auto n = var("n");
    auto a_n = SymbolicExpr::power(num(-1), n);

    auto result = checked_sequence_limit(LMCAS::lim_inf_checked(a_n, "n"));
    ASSERT_TRUE(result);
    EXPECT_EQ(result->compare(num(-1)), 0)
        << "lim_inf((-1)^n) is exactly negative one";
}

TEST(SeriesSequenceLimits, ExactAlternatingAmplitude) {
    const auto amplitude = SymbolicExpr::number(Rational(BigInt("9007199254740993"), BigInt(10)));
    auto sequence = SymbolicExpr::multiply(amplitude, SymbolicExpr::power(num(-1), var("n")));
    auto upper = lim_sup_checked(sequence, "n");
    auto lower = lim_inf_checked(sequence, "n");
    ASSERT_TRUE(upper);
    ASSERT_TRUE(lower);
    EXPECT_EQ(upper.value()->compare(amplitude), 0)
        << "upper envelope retains all rational digits";
    EXPECT_EQ(lower.value()->compare(
                  SymbolicExpr::multiply(num(-1), amplitude)->simplify()), 0)
        << "lower envelope retains all rational digits";
}

TEST(SeriesSequenceLimits, NearNegativeOneIsNotAlternating) {
    auto n = var("n");
    auto near_negative_one = SymbolicExpr::number(Rational(
        BigInt("-9999999999999"), BigInt("10000000000000")));
    auto sequence = SymbolicExpr::power(near_negative_one, n);
    expect_inconclusive_sequence_limit(lim_sup_checked(sequence, "n"));
    expect_inconclusive_sequence_limit(lim_inf_checked(sequence, "n"));
}

TEST(SeriesSequenceLimits, ContinuousTrigPhasesAreInconclusive) {
    auto n = var("n");
    auto phase = SymbolicExpr::multiply(
        SymbolicExpr::multiply(num(2), SymbolicExpr::variable("pi")), n);
    for (const auto &sequence :
         {SymbolicExpr::sin(phase), SymbolicExpr::cos(phase)}) {
        expect_inconclusive_sequence_limit(lim_sup_checked(sequence, "n"));
        expect_inconclusive_sequence_limit(lim_inf_checked(sequence, "n"));
    }
}

TEST(SeriesSequenceLimits, LimSupMonotonePolynomial) {
    auto n = var("n");
    auto a_n = SymbolicExpr::power(n, num(2));

    auto result = checked_sequence_limit(LMCAS::lim_sup_checked(a_n, "n"));
    ASSERT_TRUE(result) << "lim_sup(n^2) has a value";
    EXPECT_EQ(result->compare(SymbolicExpr::infinity()), 0)
        << "lim_sup(n^2) is positive infinity";
}

TEST(SeriesSequenceLimits, LimSupRationalSequence) {
    auto n = var("n");
    auto a_n = SymbolicExpr::multiply(n, SymbolicExpr::power(
                                             SymbolicExpr::add(n, num(1)), num(-1)));

    auto result = checked_sequence_limit(LMCAS::lim_sup_checked(a_n, "n"));
    ASSERT_TRUE(result);
    EXPECT_EQ(result->compare(num(1)), 0) << "lim_sup(n/(n+1)) = 1";
}

TEST(SeriesSequenceLimits, LimInfRationalSequence) {
    auto n = var("n");
    auto a_n = SymbolicExpr::multiply(n, SymbolicExpr::power(
                                             SymbolicExpr::add(n, num(1)), num(-1)));

    auto result = checked_sequence_limit(LMCAS::lim_inf_checked(a_n, "n"));
    ASSERT_TRUE(result);
    EXPECT_EQ(result->compare(num(1)), 0) << "lim_inf(n/(n+1)) = 1";
}

TEST(SeriesSequenceLimits, LimSupLimInfEqualForMonotone) {
    auto n = var("n");
    auto a_n = SymbolicExpr::power(n, num(-1));

    auto sup = checked_sequence_limit(LMCAS::lim_sup_checked(a_n, "n"));
    auto inf = checked_sequence_limit(LMCAS::lim_inf_checked(a_n, "n"));

    ASSERT_TRUE(sup);
    ASSERT_TRUE(inf);
    EXPECT_EQ(sup->compare(inf), 0)
        << "lim_sup(1/n) = lim_inf(1/n) for a convergent sequence";
}
