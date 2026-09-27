/**
 * @file test_series_convergence.cpp
 * @brief 收敛半径、支持域与取消契约。
 */

#include "test_common.hpp"
#include "series_engine.hpp"

using namespace LMCAS;

using Expr = std::shared_ptr<SymbolicExpr>;
using Coeffs = std::vector<Expr>;

static Expr num(int n) { return SymbolicExpr::number(n); }
static Expr rat(int p, int q) { return SymbolicExpr::number(Rational(p, q)); }
static Expr var(const std::string &name) { return SymbolicExpr::variable(name); }
static Expr inf() { return SymbolicExpr::infinity(); }

TEST(SeriesConvergence, CoefficientRadiusContracts) {
    Coeffs geometric = {num(1), num(1), num(1)};
    auto radius = LMCAS::convergence_radius_checked(geometric, "x");
    ASSERT_TRUE(radius) << "checked convergence_radius succeeds";
    ASSERT_TRUE(radius.value()) << "checked convergence_radius returns an expression";
    EXPECT_EQ(radius.value()->compare(inf()), 0)
        << "a finite coefficient polynomial has infinite radius";

    auto sequence_index = var("n");
    auto exponential_term = SymbolicExpr::power(num(2), sequence_index);
    auto exponential_radius = LMCAS::convergence_radius_checked(
        exponential_term, "n");
    EXPECT_TRUE((exponential_radius &&
                 std::abs(exponential_radius.value()->to_numeric() - 0.5) <
                     1e-12))
        << "coefficient term 2^n has radius 1/2";

    auto polynomial_term = SymbolicExpr::power(sequence_index, num(3));
    auto polynomial_radius = LMCAS::convergence_radius_checked(
        polynomial_term, "n");
    EXPECT_TRUE((polynomial_radius &&
                 std::abs(polynomial_radius.value()->to_numeric() - 1.0) <
                     1e-12))
        << "coefficient term n^3 has radius 1";
}

TEST(SeriesConvergence, RadiusInputContracts) {
    Coeffs geometric = {num(1), num(1), num(1)};
    auto empty_coeffs = LMCAS::convergence_radius_checked(Coeffs{}, "x");
    EXPECT_TRUE((!empty_coeffs &&
                 empty_coeffs.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked convergence_radius rejects empty coefficients";

    Coeffs with_null = {num(1), nullptr};
    auto null_coeff = LMCAS::convergence_radius_checked(with_null, "x");
    EXPECT_TRUE((!null_coeff &&
                 null_coeff.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked convergence_radius rejects null coefficients";

    auto bad_var = LMCAS::convergence_radius_checked(geometric, "");
    EXPECT_TRUE((!bad_var &&
                 bad_var.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked convergence_radius rejects empty variable";

    Coeffs symbolic_coeffs = {num(1), var("a"), num(1)};
    auto symbolic_radius =
        LMCAS::convergence_radius_checked(symbolic_coeffs, "x");
    ASSERT_TRUE(symbolic_radius);
    ASSERT_TRUE(symbolic_radius.value());
    EXPECT_EQ(symbolic_radius.value()->compare(inf()), 0)
        << "finite polynomial with parameter coefficients has infinite radius";

    Coeffs variable_coeffs = {num(1), var("x"), num(1)};
    auto variable_radius =
        LMCAS::convergence_radius_checked(variable_coeffs, "x");
    EXPECT_TRUE((!variable_radius &&
                 variable_radius.error().code == LMCAS::CasErrc::InvalidArgument))
        << "coefficient cannot depend on the series variable";
}

TEST(SeriesConvergence, ConvergenceCheckedContracts) {
    auto n = var("n");
    auto p_term = SymbolicExpr::power(n, num(-2));
    auto convergence = LMCAS::convergence_test_checked(p_term, "n");
    ASSERT_TRUE((convergence.has_value())) << "checked convergence_test succeeds";
    if (convergence) {
        EXPECT_TRUE((convergence.value().result == LMCAS::ConvergenceResult::Convergent)) << "checked convergence_test detects p-series convergence";
    }

    auto null_term = LMCAS::convergence_test_checked(nullptr, "n");
    EXPECT_TRUE((!null_term &&
                 null_term.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked convergence_test rejects null term";

    auto unsupported_term = SymbolicExpr::sin(n);
    auto unsupported_convergence = LMCAS::convergence_test_checked(unsupported_term, "n");
    EXPECT_TRUE((!unsupported_convergence &&
                 unsupported_convergence.error().code == LMCAS::CasErrc::Inconclusive))
        << "checked convergence_test reports Inconclusive for unsupported terms";
    LMCAS::CancellationToken token;
    token.cancel();
    LMCAS::ComputationContext cancelled_context({}, token);
    auto cancelled = LMCAS::convergence_test_checked(p_term, "n", cancelled_context);
    EXPECT_TRUE((!cancelled &&
                 cancelled.error().code == LMCAS::CasErrc::Cancelled))
        << "checked convergence_test observes cancellation";
}

TEST(SeriesConvergence, ConvergenceRadiusGeometric) {
    Coeffs coeffs = {num(1), num(1), num(1), num(1), num(1)};
    auto result = LMCAS::convergence_radius_checked(coeffs, "x");
    ASSERT_TRUE(result) << "finite polynomial radius is constructed";
    ASSERT_TRUE(result.value());
    EXPECT_EQ(result.value()->compare(inf()), 0)
        << "finite polynomial has infinite convergence radius";
}

TEST(SeriesConvergence, ConvergenceRadiusExponential) {
    Coeffs coeffs = {num(1), num(1), rat(1, 2), rat(1, 6), rat(1, 24)};
    auto result = LMCAS::convergence_radius_checked(coeffs, "x");
    ASSERT_TRUE(result) << "finite factorial polynomial radius is constructed";
    ASSERT_TRUE(result.value());
    EXPECT_EQ(result.value()->compare(inf()), 0)
        << "finite factorial prefix denotes a finite polynomial";
}

TEST(SeriesConvergence, ConvergenceRadiusHalf) {
    Coeffs coeffs = {num(1), num(2), num(4), num(8), num(16)};
    auto result = LMCAS::convergence_radius_checked(coeffs, "x");
    ASSERT_TRUE(result) << "finite geometric polynomial radius is constructed";
    ASSERT_TRUE(result.value());
    EXPECT_EQ(result.value()->compare(inf()), 0)
        << "finite geometric prefix denotes a finite polynomial";
}

TEST(SeriesConvergence, ConvergenceRadiusSingleCoeff) {
    Coeffs coeffs = {num(5)};
    auto result = LMCAS::convergence_radius_checked(coeffs, "x");
    ASSERT_TRUE(result);
    ASSERT_TRUE(result.value());
    EXPECT_EQ(result.value()->compare(inf()), 0)
        << "a single coefficient defines a polynomial with infinite radius";
}

TEST(SeriesConvergence, ConvergenceTestGeometricConvergent) {
    /**
     * @brief a_n = (1/2)^n 的比值 |a_{n+1}/a_n| = 1/2 < 1。
     * 化简器可能保留 r^(n+1)/r^n 而未约化为 r。
     */
    auto n = var("n");
    auto a_n = SymbolicExpr::power(rat(1, 2), n);

    auto info = LMCAS::convergence_test_checked(a_n, "n");
    ASSERT_TRUE(info);
    EXPECT_EQ(info.value().result, LMCAS::ConvergenceResult::Convergent)
        << "the geometric series with ratio 1/2 converges";
}

TEST(SeriesConvergence, ConvergenceTestGeometricDivergent) {
    auto n = var("n");
    auto a_n = SymbolicExpr::power(num(2), n);

    auto info = LMCAS::convergence_test_checked(a_n, "n");
    ASSERT_TRUE(info);
    EXPECT_EQ(info.value().result, LMCAS::ConvergenceResult::Divergent)
        << "the geometric series with ratio 2 diverges";
}

TEST(SeriesConvergence, ConvergenceTestPSeries) {
    auto n = var("n");
    auto a_n = SymbolicExpr::power(n, num(-2));

    auto info = LMCAS::convergence_test_checked(a_n, "n");
    ASSERT_TRUE(info);
    EXPECT_EQ(info.value().result, LMCAS::ConvergenceResult::Convergent)
        << "the p-series with exponent -2 converges";
}

TEST(SeriesConvergence, ExactPSeriesExponentBoundary) {
    auto n = var("n");
    const BigInt denominator("1" + std::string(400, '0'));
    const auto below_negative_one = SymbolicExpr::number(
        Rational(BigInt(0) - denominator - BigInt(1), denominator));
    const auto above_negative_one = SymbolicExpr::number(
        Rational(BigInt(0) - denominator + BigInt(1), denominator));

    auto convergent = convergence_test_checked(
        SymbolicExpr::power(n, below_negative_one), "n");
    ASSERT_TRUE(convergent);
    EXPECT_EQ(convergent.value().result, ConvergenceResult::Convergent)
        << "n^(-1-10^-400) is convergent";

    auto divergent = convergence_test_checked(
        SymbolicExpr::power(n, above_negative_one), "n");
    ASSERT_TRUE(divergent);
    EXPECT_EQ(divergent.value().result, ConvergenceResult::Divergent)
        << "n^(-1+10^-400) is divergent";

    auto zero_exponent = convergence_test_checked(
        SymbolicExpr::power(n, num(0)), "n");
    ASSERT_TRUE(zero_exponent);
    EXPECT_EQ(zero_exponent.value().result, ConvergenceResult::Divergent)
        << "n^0 does not tend to zero";
}
