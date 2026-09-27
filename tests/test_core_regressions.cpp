// Cross-component regression contracts for previously corrected defects.
// Each section exercises an observable boundary or invariant.

#include "test_common.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include "value.hpp"
#include "interval.hpp"
#include "polynomial.hpp"
#include "poly_utils.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "inequality_solver.hpp"
#include "transcendental_factor.hpp"

#include <climits>
#include <stdexcept>
#include <string>

using namespace LMCAS;

using LMCAS::Endpoint;
using LMCAS::Interval;
using LMCAS::IntervalUnion;
using LMCAS::PiecewiseIntervalResult;
using LMCAS::Polynomial;
using LMCAS::symbolic_to_poly;

TEST(CoreRegressions, BigintToIntSaturation) {
    // Values inside int range round-trip cleanly.
    EXPECT_TRUE((BigInt(0).to_int() == 0)) << "0 -> 0";
    EXPECT_TRUE((BigInt(123456).to_int() == 123456)) << "small positive round-trip";
    EXPECT_TRUE((BigInt(-987654).to_int() == -987654)) << "small negative round-trip";

    // Boundary cases.
    EXPECT_TRUE((BigInt(INT_MAX).to_int() == INT_MAX)) << "INT_MAX exact";
    EXPECT_TRUE((BigInt(INT_MIN).to_int() == INT_MIN)) << "INT_MIN exact";

    // Single limb but magnitude beyond int range: must saturate, not wrap.
    BigInt big = BigInt(INT_MAX) + BigInt(1);
    EXPECT_TRUE((big.to_int() == INT_MAX)) << "INT_MAX+1 -> INT_MAX (saturate)";

    BigInt small = BigInt(INT_MIN) - BigInt(1);
    EXPECT_TRUE((small.to_int() == INT_MIN)) << "INT_MIN-1 -> INT_MIN (saturate)";

    // A clearly large positive value (single limb on 64-bit) saturates positively.
    BigInt huge = BigInt("1234567890123456789");
    EXPECT_TRUE((huge.to_int() == INT_MAX)) << "1.2e18 -> INT_MAX";

    BigInt huge_neg = BigInt("-1234567890123456789");
    EXPECT_TRUE((huge_neg.to_int() == INT_MIN)) << "-1.2e18 -> INT_MIN";
}

TEST(CoreRegressions, BigintStringValidation) {
    EXPECT_THROW(([] { BigInt("12a3"); })(), std::invalid_argument) << "12a3 -> throws (was silently 1203)";

    EXPECT_THROW(([] { BigInt("12.3"); })(), std::invalid_argument) << "12.3 -> throws (was silently 123)";

    EXPECT_THROW(([] { BigInt("1 2"); })(), std::invalid_argument) << "spaces -> throws";

    // Sign prefixes are still allowed.
    EXPECT_EQ((BigInt("-42").to_string()), ("-42")) << "negative literal works";
    EXPECT_EQ((BigInt("+42").to_string()), ("42")) << "+ prefix works";

    EXPECT_THROW(([] { BigInt(""); })(), std::invalid_argument) << "empty input has no decimal digits";
    EXPECT_THROW(([] { BigInt("-"); })(), std::invalid_argument) << "sign-only input has no decimal digits";
    EXPECT_THROW(([] { BigInt("+"); })(), std::invalid_argument) << "plus-only input has no decimal digits";
}

TEST(CoreRegressions, RationalStringSign) {
    EXPECT_EQ((Rational("-3").to_string()), ("-3")) << "-3 stays negative";
    EXPECT_EQ((Rational("-12").to_string()), ("-12")) << "-12 stays negative";
    // Scientific notation negative integer: "-2e4" was being constructed as +20000.
    EXPECT_EQ((Rational("-2e4").to_string()), ("-20000")) << "-2e4 stays negative";
    EXPECT_EQ((Rational("2e4").to_string()), ("20000")) << "+2e4 still works";
    EXPECT_EQ((Rational("-1e-2").to_string()), ("-1/100")) << "-1e-2 = -1/100";
}

TEST(CoreRegressions, EndpointDefaultInit) {
    Endpoint ep{};
    EXPECT_TRUE((ep.is_open == false)) << "is_open defaults to false";
    EXPECT_TRUE((ep.is_neg_infinity == false)) << "is_neg_infinity defaults to false";
    EXPECT_TRUE((ep.is_pos_infinity == false)) << "is_pos_infinity defaults to false";
    EXPECT_TRUE((ep.value == nullptr)) << "value defaults to null";
}

TEST(CoreRegressions, ValueStringAndNumeric) {
    Value a("hello");
    Value b("world");
    EXPECT_TRUE((a.is_string())) << "string ctor sets Type::String";
    EXPECT_TRUE((!a.is_null())) << "string is no longer Null";
    EXPECT_TRUE((!(a == b))) << "distinct strings compare unequal";
    EXPECT_TRUE((a == Value("hello"))) << "same strings compare equal";
    EXPECT_EQ((a.to_string()), ("hello")) << "to_string returns content";

    // is_numeric should not include Symbolic anymore.
    Value sym(SymbolicExpr::variable("x"));
    EXPECT_TRUE((sym.is_symbolic())) << "symbolic value is symbolic";
    EXPECT_TRUE((!sym.is_numeric())) << "symbolic value is NOT numeric";

    // Concrete numeric kinds remain numeric.
    EXPECT_TRUE((Value(42).is_numeric())) << "int is numeric";
    EXPECT_TRUE((Value(BigInt(42)).is_numeric())) << "BigInt is numeric";
    EXPECT_TRUE((Value(Rational(1, 2)).is_numeric())) << "Rational is numeric";

    // String is also not numeric.
    EXPECT_TRUE((!a.is_numeric())) << "string is NOT numeric";
}

TEST(CoreRegressions, ValueAsRationalCheckedNoTruncate) {
    BigInt huge = BigInt(INT_MAX) + BigInt(100);
    Value value(huge);
    const auto rational = value.as_rational_checked();
    ASSERT_TRUE(rational);
    EXPECT_TRUE((rational.value() == Rational(huge)))
        << "checked conversion keeps large BigInt exactly";

    Value small(BigInt(7));
    const auto small_rational = small.as_rational_checked();
    ASSERT_TRUE(small_rational);
    EXPECT_TRUE((small_rational.value() == Rational(7)))
        << "small BigInt round-trips through checked conversion";
}

TEST(CoreRegressions, PolynomialVariableMismatch) {
    Polynomial<Rational> px({Rational(1), Rational(2)}, "x"); // 1 + 2x
    Polynomial<Rational> py({Rational(3), Rational(4)}, "y"); // 3 + 4y
    Polynomial<Rational> zero_x("x");
    Polynomial<Rational> zero_y("y");

    EXPECT_THROW(([&] { auto r = px + py; (void)r; })(), std::invalid_argument) << "px + py with different vars throws";

    EXPECT_THROW(([&] { auto r = px - py; (void)r; })(), std::invalid_argument) << "px - py with different vars throws";

    // Adding zero polynomials of any name must still succeed.
    auto a = px + zero_y;
    EXPECT_TRUE((a == px)) << "px + 0_y == px";
    auto b = zero_x + py;
    EXPECT_TRUE((b.coeffs == py.coeffs)) << "0_x + py has py's coefficients";
}

TEST(CoreRegressions, MatrixCloneWithNullSlot) {
    MatrixNode::DenseStorage dense;
    dense.push_back(LMCAS::detail::make_node<NumberNode>(BigInt(1)));
    dense.push_back(nullptr); // intentional empty slot
    dense.push_back(LMCAS::detail::make_node<NumberNode>(BigInt(2)));
    dense.push_back(LMCAS::detail::make_node<NumberNode>(BigInt(3)));

    bool rejected = false;
    try {
        (void)LMCAS::detail::make_node<MatrixNode>(2, 2, std::move(dense));
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    EXPECT_TRUE((rejected)) << "MatrixNode rejects null dense storage slots";
}

TEST(CoreRegressions, NumbernodeHashConsistent) {
    auto n_int = LMCAS::detail::make_node<NumberNode>(BigInt(1));
    auto n_rat = LMCAS::detail::make_node<NumberNode>(Rational(1, 1));
    auto n_dbl = LMCAS::detail::make_node<NumberNode>((lmmc_real_t)1.0);

    // BigInt and Rational are the same exact domain; approximate doubles are distinct.
    EXPECT_TRUE((n_int->compare(*n_rat) == 0)) << "BigInt 1 == Rational 1/1";
    EXPECT_TRUE((n_int->compare(*n_dbl) != 0)) << "exact BigInt 1 differs structurally from approximate 1.0";
    EXPECT_TRUE((n_rat->compare(*n_dbl) != 0)) << "exact Rational 1 differs structurally from approximate 1.0";

    // Hash must agree (otherwise equals() short-circuits to false on hash mismatch).
    EXPECT_TRUE((n_int->hash() == n_rat->hash())) << "hash(BigInt 1) == hash(Rational 1/1)";
    EXPECT_TRUE((n_int->hash() != n_dbl->hash())) << "exact and approximate hashes are domain-separated";

    EXPECT_TRUE((n_int->equals(*n_rat))) << "equals(BigInt 1, Rational 1/1)";
    EXPECT_FALSE((n_int->equals(*n_dbl))) << "exact BigInt does not equal approximate double structurally";
}

TEST(CoreRegressions, SymbolicToPolyIntegerExponents) {
    auto x = SymbolicExpr::variable("x");

    // x^2 -> degree 2 polynomial (1 0 0 1) = 0 + 0x + 1 x^2
    auto sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto p_sq = symbolic_to_poly<Rational>(sq, "x");
    EXPECT_TRUE((p_sq && p_sq.value().degree() == 2)) << "x^2 has degree 2";

    // x^(1/2) must NOT be treated as a polynomial.
    auto half = SymbolicExpr::number(Rational(1, 2));
    auto sqrt_x = SymbolicExpr::power(x, half);
    auto p_half = symbolic_to_poly<Rational>(sqrt_x, "x");
    EXPECT_TRUE((!p_half && p_half.error().code == LMCAS::CasErrc::UnsupportedExpression)) << "x^(1/2) is explicitly unsupported";

    // x^(-1) likewise.
    auto neg = SymbolicExpr::number(BigInt(-1));
    auto x_inv = SymbolicExpr::power(x, neg);
    auto p_neg = symbolic_to_poly<Rational>(x_inv, "x");
    EXPECT_TRUE((!p_neg && p_neg.error().code == LMCAS::CasErrc::UnsupportedExpression)) << "x^(-1) is explicitly unsupported";

    // Float exponent that is exactly integer-valued should still work.
    auto pf2 = SymbolicExpr::power(x, SymbolicExpr::number(2.0));
    auto p_f2 = symbolic_to_poly<Rational>(pf2, "x");
    EXPECT_TRUE((p_f2 && p_f2.value().degree() == 2)) << "x^(2.0 as double) still works";

    // Float exponent that is non-integer must NOT be expanded.
    auto pf_half = SymbolicExpr::power(x, SymbolicExpr::number(0.5));
    auto p_fhalf = symbolic_to_poly<Rational>(pf_half, "x");
    EXPECT_TRUE((!p_fhalf && p_fhalf.error().code == LMCAS::CasErrc::UnsupportedExpression)) << "x^0.5 (double) is not a polynomial";
}

TEST(CoreRegressions, ExtractCoeffNoSilentTruncate) {
    auto e = SymbolicExpr::number(1.9);
    auto b = LMCAS::extract_coeff_value<BigInt>(e);
    EXPECT_TRUE((!b && b.error().code == LMCAS::CasErrc::UnsupportedExpression)) << "1.9 cannot be represented as an integer coefficient";

    // Integer-valued double still extracts normally.
    auto e2 = SymbolicExpr::number(7.0);
    auto b2 = LMCAS::extract_coeff_value<BigInt>(e2);
    EXPECT_TRUE((b2 && b2.value() == BigInt(7))) << "7.0 -> BigInt(7)";

    auto e3 = SymbolicExpr::number(Rational(3, 2));
    auto b3 = LMCAS::extract_coeff_value<BigInt>(e3);
    EXPECT_TRUE((!b3 && b3.error().code == LMCAS::CasErrc::UnsupportedExpression)) << "3/2 cannot be represented as an integer coefficient";
}

TEST(CoreRegressions, FactorConversionErrorsPropagate) {
    const auto x = SymbolicExpr::variable("x");
    const auto y = SymbolicExpr::variable("y");
    const auto oversized = SymbolicExpr::power(x, SymbolicExpr::number(1000));
    const auto mixed = SymbolicExpr::add(oversized, SymbolicExpr::sin(y));
    const auto factored = mixed->factor_checked();
    EXPECT_TRUE((!factored && factored.error().code == CasErrc::ResourceLimit)) << "multivariate fallback cannot turn excessive degree into unchanged success";
    const auto oversized_function = SymbolicExpr::power(
        SymbolicExpr::sin(x), SymbolicExpr::number(1000));
    const auto substituted = SymbolicExpr::add(
        oversized_function, SymbolicExpr::number(1));
    const auto trans = factor_transcendental(substituted, "x");
    EXPECT_TRUE((!trans && trans.error().code == CasErrc::ResourceLimit)) << "transcendental polynomial conversion preserves excessive-degree failure";
    const auto through_factor = oversized_function->factor_checked();
    EXPECT_TRUE((!through_factor && through_factor.error().code == CasErrc::ResourceLimit)) << "checked factoring propagates the transcendental fallback failure";

    const auto unsupported = SymbolicExpr::add(x, SymbolicExpr::sin(y));
    const auto unchanged = unsupported->factor_checked();
    ASSERT_TRUE(unchanged.has_value()) << "unsupported rational coefficients remain a strategy miss";
    const auto actual_expr = unchanged.value();
    const auto expected_expr = unsupported->simplify();
    ASSERT_TRUE(actual_expr);
    ASSERT_TRUE(expected_expr);
    EXPECT_EQ(actual_expr->to_string(), expected_expr->to_string())
        << "unsupported conversion preserves the expression";

    CancellationToken token;
    token.cancel();
    ComputationContext cancelled({}, token);
    const auto stopped = mixed->factor_checked(cancelled);
    EXPECT_TRUE((!stopped && stopped.error().code == CasErrc::Cancelled)) << "checked factoring preserves cancellation";
}

TEST(CoreRegressions, SingleSolutionStrict) {
    PiecewiseIntervalResult r;
    EXPECT_TRUE((r.single_solution().is_empty())) << "empty -> empty solution";

    PiecewiseIntervalResult::Case c1{nullptr, IntervalUnion::entire_line()};
    PiecewiseIntervalResult::Case c2{nullptr, IntervalUnion::empty()};
    r.cases.push_back(c1);
    EXPECT_TRUE((r.is_single())) << "1 case is single";
    EXPECT_TRUE((r.single_solution().is_entire_line())) << "single-case solution returned";

    r.cases.push_back(c2);
    EXPECT_TRUE((!r.is_single())) << "2 cases is not single";
    EXPECT_TRUE((r.single_solution().is_empty())) << "multi-case single_solution returns empty (not first case)";
}
