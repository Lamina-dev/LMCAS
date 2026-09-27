#include "polynomial.hpp"
#include "bigint.hpp"
#include "poly_utils.hpp"
#include "symbolic.hpp"
#include "test_common.hpp"
#include <vector>

using namespace LMCAS;

namespace {

template <class Coefficient>
void expect_exact_divisor(
    const Polynomial<Coefficient> &dividend,
    const Polynomial<Coefficient> &divisor) {
    auto division = dividend.div_mod(divisor);
    EXPECT_TRUE(division.second.is_zero());
    EXPECT_TRUE(division.first * divisor == dividend);
}

} // namespace

TEST(PolynomialGcd, GcdPrimitive) {
    Polynomial<BigInt> p1(
        {BigInt(-2), BigInt(1), BigInt(1)}, "x");
    Polynomial<BigInt> p2(
        {BigInt(-3), BigInt(2), BigInt(1)}, "x");
    auto gcd = Polynomial<BigInt>::gcd(p1, p2);
    EXPECT_TRUE(gcd == Polynomial<BigInt>(
        {BigInt(-1), BigInt(1)}, "x"));
    expect_exact_divisor(p1, gcd);
    expect_exact_divisor(p2, gcd);

    Polynomial<BigInt> pa(
        {BigInt(10), BigInt(21), BigInt(12), BigInt(1)}, "x");
    Polynomial<BigInt> pb(
        {BigInt(20), BigInt(41), BigInt(22), BigInt(1)}, "x");
    auto quadratic_gcd = Polynomial<BigInt>::gcd(pa, pb);
    EXPECT_TRUE(quadratic_gcd == Polynomial<BigInt>(
        {BigInt(1), BigInt(2), BigInt(1)}, "x"));
    expect_exact_divisor(pa, quadratic_gcd);
    expect_exact_divisor(pb, quadratic_gcd);

    Polynomial<BigInt> p3(
        {BigInt(-2), BigInt(2)}, "x");
    Polynomial<BigInt> p4(
        {BigInt(-6), BigInt(6)}, "x");
    auto primitive_gcd = Polynomial<BigInt>::gcd(p3, p4);
    EXPECT_TRUE(primitive_gcd.primitive_part() == Polynomial<BigInt>(
        {BigInt(-1), BigInt(1)}, "x"));
    expect_exact_divisor(p3, primitive_gcd);
    expect_exact_divisor(p4, primitive_gcd);

    Polynomial<Rational> rational_lhs(
        {Rational(-1), Rational(0), Rational(1)}, "x");
    Polynomial<Rational> rational_rhs(
        {Rational(1), Rational(2), Rational(1)}, "x");
    auto rational_gcd =
        Polynomial<Rational>::gcd(rational_lhs, rational_rhs);
    EXPECT_TRUE(rational_gcd == Polynomial<Rational>(
        {Rational(1), Rational(1)}, "x"));
    expect_exact_divisor(rational_lhs, rational_gcd);
    expect_exact_divisor(rational_rhs, rational_gcd);
}

static void test_non_linear_and_rational_gcd(
    const std::shared_ptr<SymbolicExpr> &x,
    const std::shared_ptr<SymbolicExpr> &y) {
    auto quadratic_common = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(SymbolicExpr::multiply(x, y),
                          SymbolicExpr::power(y, SymbolicExpr::number(2))));
    auto harder_lhs = SymbolicExpr::multiply(
                          quadratic_common, SymbolicExpr::add(x, SymbolicExpr::number(2)))
                          ->expand();
    auto harder_rhs = SymbolicExpr::multiply(
                          quadratic_common, SymbolicExpr::add(y, SymbolicExpr::number(3)))
                          ->expand();
    ComputationContext harder_context;
    auto harder_result = symbolic_polynomial_gcd(
        *harder_lhs, *harder_rhs, harder_context);
    ASSERT_TRUE(harder_result) << harder_result.error().message;
    EXPECT_TRUE(test_proved_equivalent(
        harder_result.value(), quadratic_common));

    auto half = SymbolicExpr::number(Rational(1, 2));
    auto rational_common = SymbolicExpr::add(x, half);
    auto rational_lhs = SymbolicExpr::multiply(
        rational_common, SymbolicExpr::add(x, SymbolicExpr::number(2)));
    auto rational_rhs = SymbolicExpr::multiply(
        rational_common, SymbolicExpr::add(x, SymbolicExpr::number(3)));
    ComputationContext rational_context;
    auto rational_result = symbolic_polynomial_gcd(
        *rational_lhs, *rational_rhs, rational_context);
    ASSERT_TRUE(rational_result) << rational_result.error().message;
    EXPECT_TRUE(test_proved_equivalent(
        rational_result.value(), rational_common));
}

TEST(PolynomialGcd, SymbolicPolynomialGcd) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto one = SymbolicExpr::number(1);
    auto common = SymbolicExpr::add(x, y);
    auto lhs = SymbolicExpr::multiply(
                   common, SymbolicExpr::add(x, one))
                   ->expand();
    auto rhs = SymbolicExpr::multiply(
                   common, SymbolicExpr::add(y, one))
                   ->expand();

    ComputationContext context;
    auto result = symbolic_polynomial_gcd(*lhs, *rhs, context);
    ASSERT_TRUE(result) << result.error().message;
    EXPECT_TRUE(test_proved_equivalent(result.value(), common));

    test_non_linear_and_rational_gcd(x, y);

    auto sine = SymbolicExpr::sin(x);
    ComputationContext unsupported_context;
    auto unsupported =
        symbolic_polynomial_gcd(*sine, *lhs, unsupported_context);
    ASSERT_FALSE(unsupported);
    EXPECT_EQ(
        unsupported.error().code, CasErrc::UnsupportedExpression);

    auto approximate = SymbolicExpr::add(
        x, SymbolicExpr::number(static_cast<lmmc_real_t>(0.5)));
    ComputationContext approximate_context;
    auto approximate_result = symbolic_polynomial_gcd(
        *approximate, *lhs, approximate_context);
    ASSERT_FALSE(approximate_result);
    EXPECT_EQ(
        approximate_result.error().code,
        CasErrc::UnsupportedExpression);

    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext limited_context(limits);
    auto limited =
        symbolic_polynomial_gcd(*lhs, *rhs, limited_context);
    ASSERT_FALSE(limited);
    EXPECT_EQ(limited.error().code, CasErrc::ResourceLimit);
}

namespace {

template <class Operation>
void expect_variable_mismatch(Operation operation, const char *message) {
    bool rejected = false;
    try {
        operation();
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    EXPECT_TRUE((rejected)) << message;
}

template <class Coefficient>
void test_named_polynomial_rings(const char *coefficient_name) {
    SCOPED_TRACE(std::string("Polynomial: named variable rings for ") + coefficient_name);
    using Poly = Polynomial<Coefficient>;
    const Poly x(std::vector<Coefficient>{Coefficient(0), Coefficient(1)}, "x");
    const Poly y(std::vector<Coefficient>{Coefficient(0), Coefficient(1)}, "y");
    const Poly one_x(Coefficient(1), "x");
    const Poly one_y(Coefficient(1), "y");
    for (const auto &operands : {std::pair{x, y}, std::pair{one_x, one_y}}) {
        const auto &lhs = operands.first;
        const auto &rhs = operands.second;
        EXPECT_FALSE((lhs == rhs)) << "nonzero polynomials in different named rings are unequal";
        expect_variable_mismatch([&] { (void)(lhs + rhs); }, "addition rejects different rings");
        expect_variable_mismatch([&] { (void)(lhs - rhs); }, "subtraction rejects different rings");
        expect_variable_mismatch([&] { (void)(lhs * rhs); }, "multiplication rejects different rings");
        expect_variable_mismatch([&] { (void)lhs.div_mod(rhs); }, "division rejects different rings");
        expect_variable_mismatch([&] { (void)lhs.pseudo_div_mod(rhs); },
                                 "pseudo-division rejects different rings");
        expect_variable_mismatch([&] { (void)lhs.pseudo_div_mod_rem(rhs); },
                                 "pseudo-remainder rejects different rings");
        expect_variable_mismatch([&] { (void)Poly::gcd(lhs, rhs); }, "GCD rejects different rings");
    }
    const Poly y_squared(std::vector<Coefficient>{Coefficient(0), Coefficient(0), Coefficient(1)}, "y");
    expect_variable_mismatch([&] { (void)x.div_mod(y_squared); },
                             "low-degree division checks the ring before returning a remainder");
    expect_variable_mismatch([&] { (void)x.pseudo_div_mod(y_squared); },
                             "low-degree pseudo-division checks the ring before returning a remainder");
}

template <class Coefficient>
void expect_zero_polynomial_arithmetic(const Polynomial<Coefficient> &zero,
                                       const Polynomial<Coefficient> &p,
                                       const Polynomial<Coefficient> &negative_p) {
    for (const auto &result : {zero + p, p + zero, p - zero}) {
        EXPECT_TRUE((result == p && result.variable_name == "y")) << "additive zero preserves the nonzero polynomial and its ring";
    }
    const auto difference = zero - p;
    EXPECT_TRUE((difference == negative_p && difference.variable_name == "y")) << "zero-p is the negated polynomial in p's ring";
    for (const auto &product : {zero * p, p * zero}) {
        EXPECT_TRUE((product.is_zero() && product.variable_name == "y")) << "zero products retain the nonzero operand's ring";
    }
}

template <class Coefficient>
void expect_zero_polynomial_division_and_gcd(const Polynomial<Coefficient> &zero,
                                             const Polynomial<Coefficient> &p) {
    using Poly = Polynomial<Coefficient>;
    for (const auto &division : {zero.div_mod(p), zero.pseudo_div_mod(p)}) {
        EXPECT_TRUE((division.first.is_zero() && division.second.is_zero() &&
                     division.first.variable_name == "y" &&
                     division.second.variable_name == "y"))
            << "zero dividend yields zero quotient and remainder in the divisor ring";
    }
    EXPECT_TRUE((Poly::gcd(zero, p) == p && Poly::gcd(p, zero) == p)) << "GCD with zero retains the nonzero operand";
}

template <class Coefficient>
void test_polynomial_zero_rings(const char *coefficient_name) {
    SCOPED_TRACE(std::string("Polynomial: zero ring embedding for ") + coefficient_name);
    using Poly = Polynomial<Coefficient>;
    const Poly zero;
    const Poly named_zero("z");
    const Poly p(std::vector<Coefficient>{Coefficient(2), Coefficient(4)}, "y");
    const Poly negative_p(std::vector<Coefficient>{Coefficient(-2), Coefficient(-4)}, "y");
    for (const auto *z : {&zero, &named_zero}) {
        expect_zero_polynomial_arithmetic(*z, p, negative_p);
        expect_zero_polynomial_division_and_gcd(*z, p);
    }
    EXPECT_TRUE((zero == named_zero)) << "zero polynomials compare equal across variable rings";
    for (const auto &result : {named_zero + zero, named_zero - zero, named_zero * zero,
                               Poly::gcd(named_zero, zero)}) {
        EXPECT_TRUE((result.is_zero() && result.variable_name == "z")) << "two zero operands retain the left ring";
    }
    bool division_rejected = false;
    try {
        (void)p.div_mod(zero);
    } catch (const std::runtime_error &) {
        division_rejected = true;
    }
    EXPECT_TRUE((division_rejected)) << "zero divisor remains runtime_error before ring checks";
    bool pseudo_division_rejected = false;
    try {
        (void)p.pseudo_div_mod(zero);
    } catch (const std::runtime_error &) {
        pseudo_division_rejected = true;
    }
    EXPECT_TRUE((pseudo_division_rejected)) << "zero pseudo-divisor remains runtime_error before ring checks";
}

} // namespace

TEST(PolynomialGcd, NamedPolynomialRingsRational) {
    test_named_polynomial_rings<Rational>("Rational");
}

TEST(PolynomialGcd, NamedPolynomialRingsBigint) {
    test_named_polynomial_rings<BigInt>("BigInt");
}

TEST(PolynomialGcd, PolynomialZeroRingsRational) {
    test_polynomial_zero_rings<Rational>("Rational");
}

TEST(PolynomialGcd, PolynomialZeroRingsBigint) {
    test_polynomial_zero_rings<BigInt>("BigInt");
}
