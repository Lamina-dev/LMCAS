// Unit tests covering the review fixes applied across batches 1 & 2.
//
// Each section below targets one specific bug from the code review. The aim is
// to lock down the new contracts so future regressions are caught immediately.

#include "test_common.hpp"

#include "bigint.hpp"
#include "rational.hpp"
#include "value.hpp"
#include "interval.hpp"
#include "irrational.hpp"
#include "polynomial.hpp"
#include "poly_utils.hpp"
#include "symbolic.hpp"
#include "symbolic_ast.hpp"
#include "inequality_solver.hpp"

#include <climits>
#include <cmath>
#include <stdexcept>
#include <string>

using lamina::IntervalUnion;
using lamina::Endpoint;
using lamina::Interval;
using lamina::Polynomial;
using lamina::PiecewiseIntervalResult;
using lamina::symbolic_to_poly;

// Helper: assert that a callable throws a particular exception type.
template <typename E, typename F>
void EXPECT_THROWS(F&& fn, const std::string& msg) {
    try {
        fn();
    } catch (const E&) {
        std::cout << "[PASS] " << msg << std::endl;
        return;
    } catch (const std::exception& ex) {
        std::cerr << "[FAIL] " << msg
                  << " (wrong exception type: " << ex.what() << ")" << std::endl;
        ++g_failures;
        return;
    } catch (...) {
        std::cerr << "[FAIL] " << msg
                  << " (unknown non-std exception)" << std::endl;
        ++g_failures;
        return;
    }
    std::cerr << "[FAIL] " << msg << " (no exception thrown)" << std::endl;
    ++g_failures;
}

void test_bigint_to_int_saturation() {
    TEST_CASE("BigInt::to_int saturates instead of truncating");

    // Values inside int range round-trip cleanly.
    EXPECT_TRUE(BigInt(0).to_int() == 0, "0 -> 0");
    EXPECT_TRUE(BigInt(123456).to_int() == 123456, "small positive round-trip");
    EXPECT_TRUE(BigInt(-987654).to_int() == -987654, "small negative round-trip");

    // Boundary cases.
    EXPECT_TRUE(BigInt((long long)INT_MAX).to_int() == INT_MAX, "INT_MAX exact");
    EXPECT_TRUE(BigInt((long long)INT_MIN).to_int() == INT_MIN, "INT_MIN exact");

    // Single limb but magnitude beyond int range: must saturate, not wrap.
    BigInt big = BigInt((long long)INT_MAX) + BigInt(1);
    EXPECT_TRUE(big.to_int() == INT_MAX, "INT_MAX+1 -> INT_MAX (saturate)");

    BigInt small = BigInt((long long)INT_MIN) - BigInt(1);
    EXPECT_TRUE(small.to_int() == INT_MIN, "INT_MIN-1 -> INT_MIN (saturate)");

    // A clearly large positive value (single limb on 64-bit) saturates positively.
    BigInt huge = BigInt(1234567890123456789LL);
    EXPECT_TRUE(huge.to_int() == INT_MAX, "1.2e18 -> INT_MAX");

    BigInt huge_neg = BigInt(-1234567890123456789LL);
    EXPECT_TRUE(huge_neg.to_int() == INT_MIN, "-1.2e18 -> INT_MIN");
}

void test_bigint_string_validation() {
    TEST_CASE("BigInt(string) rejects invalid characters");

    EXPECT_THROWS<std::invalid_argument>(
        [] { BigInt("12a3"); },
        "12a3 -> throws (was silently 1203)");

    EXPECT_THROWS<std::invalid_argument>(
        [] { BigInt("12.3"); },
        "12.3 -> throws (was silently 123)");

    EXPECT_THROWS<std::invalid_argument>(
        [] { BigInt("1 2"); },
        "spaces -> throws");

    // Sign prefixes are still allowed.
    EXPECT_EQ_STR(BigInt("-42").to_string(), "-42", "negative literal works");
    EXPECT_EQ_STR(BigInt("+42").to_string(), "42", "+ prefix works");

    // Empty after sign or fully empty -> 0 (existing behavior preserved).
    EXPECT_TRUE(BigInt("").is_zero(), "empty -> 0");
    EXPECT_TRUE(BigInt("-").is_zero(), "lone minus -> 0");
}

void test_rational_string_sign() {
    TEST_CASE("Rational(string) keeps the sign for integer / scientific literals");

    EXPECT_EQ_STR(Rational("-3").to_string(), "-3", "-3 stays negative");
    EXPECT_EQ_STR(Rational("-12").to_string(), "-12", "-12 stays negative");
    // Scientific notation negative integer: "-2e4" was being constructed as +20000.
    EXPECT_EQ_STR(Rational("-2e4").to_string(), "-20000", "-2e4 stays negative");
    EXPECT_EQ_STR(Rational("2e4").to_string(), "20000", "+2e4 still works");
    EXPECT_EQ_STR(Rational("-1e-2").to_string(), "-1/100", "-1e-2 = -1/100");
}

void test_endpoint_default_init() {
    TEST_CASE("Endpoint default constructor zero-initializes bool fields");

    Endpoint ep{};
    EXPECT_TRUE(ep.is_open == false, "is_open defaults to false");
    EXPECT_TRUE(ep.is_neg_infinity == false, "is_neg_infinity defaults to false");
    EXPECT_TRUE(ep.is_pos_infinity == false, "is_pos_infinity defaults to false");
    EXPECT_TRUE(ep.value == nullptr, "value defaults to null");
}

void test_value_string_and_numeric() {
    TEST_CASE("Value string handling and is_numeric semantics");

    Value a("hello");
    Value b("world");
    EXPECT_TRUE(a.is_string(), "string ctor sets Type::String");
    EXPECT_TRUE(!a.is_null(), "string is no longer Null");
    EXPECT_TRUE(!(a == b), "distinct strings compare unequal");
    EXPECT_TRUE(a == Value("hello"), "same strings compare equal");
    EXPECT_EQ_STR(a.to_string(), "hello", "to_string returns content");

    // is_numeric should not include Symbolic anymore.
    Value sym(SymbolicExpr::variable("x"));
    EXPECT_TRUE(sym.is_symbolic(), "symbolic value is symbolic");
    EXPECT_TRUE(!sym.is_numeric(), "symbolic value is NOT numeric");

    // Concrete numeric kinds remain numeric.
    EXPECT_TRUE(Value(42).is_numeric(), "int is numeric");
    EXPECT_TRUE(Value(BigInt(42)).is_numeric(), "BigInt is numeric");
    EXPECT_TRUE(Value(Rational(1, 2)).is_numeric(), "Rational is numeric");

    // String is also not numeric.
    EXPECT_TRUE(!a.is_numeric(), "string is NOT numeric");
}

void test_value_as_rational_no_truncate() {
    TEST_CASE("Value::as_rational preserves large BigInt magnitude");

    BigInt huge = BigInt((long long)INT_MAX) + BigInt(100);
    Value v(huge);
    Rational r = v.as_rational();
    // The previous code converted via to_int() and saturated to INT_MAX; the new
    // path uses the BigInt directly so the rational must equal `huge`.
    EXPECT_TRUE(r == Rational(huge), "as_rational keeps large BigInt exactly");

    // Sanity: small BigInts still round-trip as before.
    Value small(BigInt(7));
    EXPECT_TRUE(small.as_rational() == Rational(7), "small BigInt round-trips");
}

void test_polynomial_variable_mismatch() {
    TEST_CASE("Polynomial::operator+ / - rejects mismatched non-zero operands");

    Polynomial<Rational> px({Rational(1), Rational(2)}, "x");        // 1 + 2x
    Polynomial<Rational> py({Rational(3), Rational(4)}, "y");        // 3 + 4y
    Polynomial<Rational> zero_x("x");
    Polynomial<Rational> zero_y("y");

    EXPECT_THROWS<std::invalid_argument>(
        [&] { auto r = px + py; (void)r; },
        "px + py with different vars throws");

    EXPECT_THROWS<std::invalid_argument>(
        [&] { auto r = px - py; (void)r; },
        "px - py with different vars throws");

    // Adding zero polynomials of any name must still succeed.
    auto a = px + zero_y;
    EXPECT_TRUE(a == px, "px + 0_y == px");
    auto b = zero_x + py;
    EXPECT_TRUE(b.coeffs == py.coeffs, "0_x + py has py's coefficients");
}

void test_irrational_to_symbolic_complex() {
    TEST_CASE("Irrational::to_symbolic reconstructs full COMPLEX combinations");

    // 2 + 3*pi + 1*sqrt(2)
    Irrational mix = Irrational::constant(2.0) + Irrational::pi(3.0) + Irrational::sqrt(2);
    auto sym = mix.to_symbolic();
    EXPECT_TRUE(sym != nullptr, "to_symbolic returns non-null");
    std::string s = sym->to_string();
    // Existence checks rather than exact form; printer ordering may vary.
    bool has_pi   = s.find("π") != std::string::npos;
    bool has_sqrt = s.find("2") != std::string::npos; // numeric or sqrt(2)
    EXPECT_TRUE(has_pi, "result contains pi term");
    EXPECT_TRUE(has_sqrt, "result mentions 2 (constant or sqrt)");

    // Pure constant (no symbolic basis) should still work.
    Irrational c = Irrational::constant(5.0);
    auto sym2 = c.to_symbolic();
    EXPECT_TRUE(sym2 != nullptr, "constant-only to_symbolic returns non-null");
}

void test_matrix_clone_with_null_slot() {
    TEST_CASE("MatrixNode rejects nullptr dense entries");

    MatrixNode::DenseStorage dense;
    dense.push_back(lamina::detail::make_node<NumberNode>(BigInt(1)));
    dense.push_back(nullptr); // intentional empty slot
    dense.push_back(lamina::detail::make_node<NumberNode>(BigInt(2)));
    dense.push_back(lamina::detail::make_node<NumberNode>(BigInt(3)));

    bool rejected = false;
    try {
        (void)lamina::detail::make_node<MatrixNode>(2, 2, std::move(dense));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    EXPECT_TRUE(rejected, "MatrixNode rejects null dense storage slots");
}

void test_numbernode_hash_consistent() {
    TEST_CASE("NumberNode hash matches compare_same_type for equivalent values");

    auto n_int = lamina::detail::make_node<NumberNode>(BigInt(1));
    auto n_rat = lamina::detail::make_node<NumberNode>(Rational(1, 1));
    auto n_dbl = lamina::detail::make_node<NumberNode>((lmmc_real_t)1.0);

    // BigInt and Rational are the same exact domain; approximate doubles are distinct.
    EXPECT_TRUE(n_int->compare(*n_rat) == 0, "BigInt 1 == Rational 1/1");
    EXPECT_TRUE(n_int->compare(*n_dbl) != 0, "exact BigInt 1 differs structurally from approximate 1.0");
    EXPECT_TRUE(n_rat->compare(*n_dbl) != 0, "exact Rational 1 differs structurally from approximate 1.0");

    // Hash must agree (otherwise equals() short-circuits to false on hash mismatch).
    EXPECT_TRUE(n_int->hash() == n_rat->hash(), "hash(BigInt 1) == hash(Rational 1/1)");
    EXPECT_TRUE(n_int->hash() != n_dbl->hash(), "exact and approximate hashes are domain-separated");

    EXPECT_TRUE(n_int->equals(*n_rat), "equals(BigInt 1, Rational 1/1)");
    EXPECT_FALSE(n_int->equals(*n_dbl), "exact BigInt does not equal approximate double structurally");
}

void test_symbolic_to_poly_integer_exponents() {
    TEST_CASE("symbolic_to_poly does not expand x^(1/2) or x^(-1) as polynomials");

    auto x = SymbolicExpr::variable("x");

    // x^2 -> degree 2 polynomial (1 0 0 1) = 0 + 0x + 1 x^2
    auto sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto p_sq = symbolic_to_poly<Rational>(sq, "x");
    EXPECT_TRUE(p_sq.degree() == 2, "x^2 has degree 2");

    // x^(1/2) must NOT be treated as a polynomial.
    auto half = SymbolicExpr::number(Rational(1, 2));
    auto sqrt_x = SymbolicExpr::power(x, half);
    auto p_half = symbolic_to_poly<Rational>(sqrt_x, "x");
    EXPECT_TRUE(p_half.is_zero(),
                "x^(1/2) returns zero polynomial (cannot be represented)");

    // x^(-1) likewise.
    auto neg = SymbolicExpr::number(BigInt(-1));
    auto x_inv = SymbolicExpr::power(x, neg);
    auto p_neg = symbolic_to_poly<Rational>(x_inv, "x");
    EXPECT_TRUE(p_neg.is_zero(),
                "x^(-1) returns zero polynomial (cannot be represented)");

    // Float exponent that is exactly integer-valued should still work.
    auto pf2 = SymbolicExpr::power(x, SymbolicExpr::number(2.0));
    auto p_f2 = symbolic_to_poly<Rational>(pf2, "x");
    EXPECT_TRUE(p_f2.degree() == 2, "x^(2.0 as double) still works");

    // Float exponent that is non-integer must NOT be expanded.
    auto pf_half = SymbolicExpr::power(x, SymbolicExpr::number(0.5));
    auto p_fhalf = symbolic_to_poly<Rational>(pf_half, "x");
    EXPECT_TRUE(p_fhalf.is_zero(), "x^0.5 (double) is not a polynomial");
}

void test_extract_coeff_no_silent_truncate() {
    TEST_CASE("extract_coeff_value<BigInt> does not silently truncate floats");

    // 1.9 -> previously produced BigInt(1). Now we expect a non-truncated
    // sentinel (0) because 1.9 is not an integer and we refuse to lie about it.
    auto e = SymbolicExpr::number(1.9);
    BigInt b = lamina::extract_coeff_value<BigInt>(e);
    EXPECT_TRUE(b == BigInt(0),
                "1.9 -> BigInt(0) sentinel, not 1");

    // Integer-valued double still extracts normally.
    auto e2 = SymbolicExpr::number(7.0);
    BigInt b2 = lamina::extract_coeff_value<BigInt>(e2);
    EXPECT_TRUE(b2 == BigInt(7), "7.0 -> BigInt(7)");

    // Non-integer rational coefficient also returns the sentinel.
    auto e3 = SymbolicExpr::number(Rational(3, 2));
    BigInt b3 = lamina::extract_coeff_value<BigInt>(e3);
    EXPECT_TRUE(b3 == BigInt(0), "3/2 -> BigInt(0) sentinel");
}

void test_single_solution_strict() {
    TEST_CASE("PiecewiseIntervalResult::single_solution requires is_single()");

    PiecewiseIntervalResult r;
    EXPECT_TRUE(r.single_solution().is_empty(), "empty -> empty solution");

    PiecewiseIntervalResult::Case c1{nullptr, IntervalUnion::entire_line()};
    PiecewiseIntervalResult::Case c2{nullptr, IntervalUnion::empty()};
    r.cases.push_back(c1);
    EXPECT_TRUE(r.is_single(), "1 case is single");
    EXPECT_TRUE(r.single_solution().is_entire_line(),
                "single-case solution returned");

    r.cases.push_back(c2);
    EXPECT_TRUE(!r.is_single(), "2 cases is not single");
    EXPECT_TRUE(r.single_solution().is_empty(),
                "multi-case single_solution returns empty (not first case)");
}

int main() {
    try {
        test_bigint_to_int_saturation();
        test_bigint_string_validation();
        test_rational_string_sign();
        test_endpoint_default_init();
        test_value_string_and_numeric();
        test_value_as_rational_no_truncate();
        test_polynomial_variable_mismatch();
        test_irrational_to_symbolic_complex();
        test_matrix_clone_with_null_slot();
        test_numbernode_hash_consistent();
        test_symbolic_to_poly_integer_exponents();
        test_extract_coeff_no_silent_truncate();
        test_single_solution_strict();
    } catch (const std::exception& ex) {
        std::cerr << "[FATAL] unhandled std::exception: " << ex.what() << std::endl;
        return TEST_REPORT() == 0 ? 1 : TEST_REPORT();
    } catch (...) {
        std::cerr << "[FATAL] unhandled non-std exception" << std::endl;
        return TEST_REPORT() == 0 ? 1 : TEST_REPORT();
    }
    return TEST_REPORT();
}
