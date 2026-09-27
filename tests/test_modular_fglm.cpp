#include "modular_arithmetic.hpp"
#include "fglm.hpp"
#include "monomial_order.hpp"
#include "test_common.hpp"

using namespace LMCAS;

namespace {

FGLMPoly multiply_term(const FGLMPoly &polynomial,
                       const Monomial &monomial,
                       const Rational &coefficient,
                       const MonomialOrder &order) {
    FGLMPoly result(polynomial.num_vars);
    for (const auto &[term, value] : polynomial.terms) {
        Monomial product(polynomial.num_vars, 0);
        for (std::size_t index = 0; index < polynomial.num_vars; ++index) {
            product[index] = term[index] + monomial[index];
        }
        result.add_term(product, value * coefficient);
    }
    result.sort_terms(order);
    result.normalize();
    return result;
}

FGLMPoly s_polynomial(const FGLMPoly &left, const FGLMPoly &right,
                      const MonomialOrder &order) {
    const Monomial common =
        lcm_monomial(left.lead_monomial(), right.lead_monomial());
    Monomial left_multiplier(left.num_vars, 0);
    Monomial right_multiplier(right.num_vars, 0);
    for (std::size_t index = 0; index < left.num_vars; ++index) {
        left_multiplier[index] = common[index] - left.lead_monomial()[index];
        right_multiplier[index] = common[index] - right.lead_monomial()[index];
    }
    auto result = multiply_term(
        left, left_multiplier, Rational(1) / left.lead_coeff(), order);
    auto subtrahend = multiply_term(
        right, right_multiplier, Rational(-1) / right.lead_coeff(), order);
    for (const auto &[term, coefficient] : subtrahend.terms) {
        result.add_term(term, coefficient);
    }
    result.sort_terms(order);
    result.normalize();
    return result;
}

void expect_same_ideal(const std::vector<FGLMPoly> &source,
                       const MonomialOrder &source_order,
                       const std::vector<FGLMPoly> &target,
                       const MonomialOrder &target_order) {
    for (const auto &generator : source) {
        EXPECT_TRUE(normal_form(generator, target, target_order).is_zero())
            << "every source generator reduces to zero by the target basis";
    }
    for (const auto &generator : target) {
        EXPECT_TRUE(normal_form(generator, source, source_order).is_zero())
            << "every target generator reduces to zero by the source basis";
    }
    for (std::size_t left = 0; left < target.size(); ++left) {
        for (std::size_t right = left + 1; right < target.size(); ++right) {
            auto s = s_polynomial(target[left], target[right], target_order);
            EXPECT_TRUE(normal_form(s, target, target_order).is_zero())
                << "every target S-polynomial reduces to zero";
        }
    }
}

} // namespace

TEST(LmcasModularFglm, Modint) {
    int64_t p = 1000000007;

    ModInt a(5, p);
    ModInt b(3, p);

    EXPECT_TRUE(((a + b).value() == 8)) << ("ModInt addition");
    EXPECT_TRUE(((a - b).value() == 2)) << ("ModInt subtraction");
    EXPECT_TRUE(((a * b).value() == 15)) << ("ModInt multiplication");

    ModInt c = a / b;
    EXPECT_TRUE(((c * b).value() == 5)) << ("ModInt division inverse roundtrip");

    ModInt neg(-1, p);
    EXPECT_TRUE((neg.value() == p - 1)) << ("ModInt negative normalization");

    ModInt base(2, p);
    ModInt result = ModInt::pow(base, 10);
    EXPECT_TRUE((result.value() == 1024)) << ("ModInt power");
    EXPECT_TRUE((ModInt::pow(ModInt(2, 10), 5).value() == 2)) << ("ModInt power with even modulus");

    ModInt inv = b.inverse();
    EXPECT_TRUE(((b * inv).value() == 1)) << ("ModInt inverse");
}

TEST(LmcasModularFglm, Crt) {
    auto [x1, m1] = crt(2, 3, 3, 5);
    EXPECT_TRUE((m1 == 15)) << ("CRT modulus is product");
    EXPECT_TRUE((x1 % 3 == 2)) << ("CRT satisfies first residue");
    EXPECT_TRUE((x1 % 5 == 3)) << ("CRT satisfies second residue");

    std::vector<int64_t> residues = {2, 3, 2};
    std::vector<int64_t> primes = {3, 5, 7};
    auto [x2, m2] = multi_crt(residues, primes);
    EXPECT_TRUE((m2 == 105)) << ("multi CRT modulus is product");
    EXPECT_TRUE((x2 % 3 == 2)) << ("multi CRT satisfies first residue");
    EXPECT_TRUE((x2 % 5 == 3)) << ("multi CRT satisfies second residue");
    EXPECT_TRUE((x2 % 7 == 2)) << ("multi CRT satisfies third residue");
}

TEST(LmcasModularFglm, RationalReconstruction) {
    int64_t p = 1000000007;
    ModInt three(3, p);
    ModInt seven(7, p);
    ModInt encoded = three / seven;

    auto reconstructed =
        rational_reconstruction_checked(encoded.value(), p);
    ASSERT_TRUE(reconstructed) << reconstructed.error().message;
    auto [num, den] = std::move(reconstructed.value());
    EXPECT_EQ(num, 3);
    EXPECT_EQ(den, 7);

    ModInt neg2(-2, p);
    ModInt five(5, p);
    ModInt encoded2 = neg2 / five;
    auto reconstructed2 =
        rational_reconstruction_checked(encoded2.value(), p);
    ASSERT_TRUE(reconstructed2) << reconstructed2.error().message;
    auto [num2, den2] = std::move(reconstructed2.value());
    EXPECT_EQ(num2, -2);
    EXPECT_EQ(den2, 5);
}

TEST(LmcasModularFglm, FglmHelpers) {
    size_t n = 2;

    FGLMPoly g1(n);
    g1.add_term({2, 0}, Rational(1));
    g1.sort_terms(MonomialOrder::grevlex());

    FGLMPoly g2(n);
    g2.add_term({0, 2}, Rational(1));
    g2.sort_terms(MonomialOrder::grevlex());

    std::vector<FGLMPoly> basis = {g1, g2};

    EXPECT_TRUE((is_zero_dimensional(basis, n))) << ("basis <x^2,y^2> is zero-dimensional");
    int dim = quotient_dimension(basis, n);
    EXPECT_TRUE((dim == 4)) << ("quotient dimension of <x^2,y^2> is 4");

    FGLMPoly x3(n);
    x3.add_term({3, 0}, Rational(1));
    x3.sort_terms(MonomialOrder::grevlex());

    FGLMPoly nf = normal_form(x3, basis, MonomialOrder::grevlex());
    EXPECT_TRUE((nf.is_zero())) << ("x^3 reduces to zero modulo <x^2,y^2>");

    FGLMPoly xy(n);
    xy.add_term({1, 1}, Rational(1));
    xy.sort_terms(MonomialOrder::grevlex());

    FGLMPoly nf2 = normal_form(xy, basis, MonomialOrder::grevlex());
    EXPECT_TRUE((!nf2.is_zero())) << ("xy is not reduced to zero modulo <x^2,y^2>");
    EXPECT_TRUE((nf2.terms.size() == 1)) << ("normal form of xy has one term");
    EXPECT_TRUE((!nf2.terms.empty() && nf2.terms[0].first == Monomial({1, 1}))) << ("normal form of xy keeps monomial xy");
}

TEST(LmcasModularFglm, FglmConversionPreservesIdealAndCompletesBasis) {
    const std::size_t variables = 2;
    const auto source_order = MonomialOrder::grevlex();
    const auto target_order = MonomialOrder::lex();

    FGLMPoly x_squared_minus_y(variables);
    x_squared_minus_y.add_term({2, 0}, Rational(1));
    x_squared_minus_y.add_term({0, 1}, Rational(-1));
    x_squared_minus_y.sort_terms(source_order);

    FGLMPoly y_squared_minus_one(variables);
    y_squared_minus_one.add_term({0, 2}, Rational(1));
    y_squared_minus_one.add_term({0, 0}, Rational(-1));
    y_squared_minus_one.sort_terms(source_order);

    std::vector<FGLMPoly> source{
        x_squared_minus_y, y_squared_minus_one};
    auto target = fglm_convert(
        source, source_order, target_order, variables);

    ASSERT_FALSE(target.empty());
    expect_same_ideal(source, source_order, target, target_order);
}

TEST(LmcasModularFglm, ReductionBeyondFormerStepCapReturnsCompleteRemainder) {
    const auto order = MonomialOrder::lex();
    FGLMPoly x_minus_one(1);
    x_minus_one.add_term({1}, Rational(1));
    x_minus_one.add_term({0}, Rational(-1));
    x_minus_one.sort_terms(order);

    FGLMPoly high_power(1);
    high_power.add_term({10001}, Rational(1));
    auto remainder = normal_form(high_power, {x_minus_one}, order);

    ASSERT_EQ(remainder.terms.size(), 1U);
    EXPECT_EQ(remainder.terms.front().first, Monomial({0}));
    EXPECT_EQ(remainder.terms.front().second, Rational(1));
}

TEST(LmcasModularFglm, IncompleteConversionThrowsInsteadOfReturningPartialBasis) {
    const std::size_t variables = 2;
    FGLMPoly x_to_fiftieth(variables);
    x_to_fiftieth.add_term({50, 0}, Rational(1));
    x_to_fiftieth.sort_terms(MonomialOrder::lex());
    FGLMPoly y(variables);
    y.add_term({0, 1}, Rational(1));
    y.sort_terms(MonomialOrder::lex());

    EXPECT_THROW(
        (void)fglm_convert(
            {x_to_fiftieth, y}, MonomialOrder::lex(),
            MonomialOrder::grevlex(), variables),
        std::runtime_error);
}
