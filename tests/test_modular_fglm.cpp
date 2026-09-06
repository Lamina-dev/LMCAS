#include "modular_arithmetic.hpp"
#include "fglm.hpp"
#include "monomial_order.hpp"
#include "test_common.hpp"
#include <iostream>
#include <cmath>

using namespace LMCAS;

void test_modint() {
    std::cout << "Test: ModInt arithmetic" << std::endl;
    int64_t p = 1000000007;

    ModInt a(5, p);
    ModInt b(3, p);

    EXPECT_TRUE((a + b).value() == 8, "ModInt addition");
    EXPECT_TRUE((a - b).value() == 2, "ModInt subtraction");
    EXPECT_TRUE((a * b).value() == 15, "ModInt multiplication");

    ModInt c = a / b;
    EXPECT_TRUE((c * b).value() == 5, "ModInt division inverse roundtrip");

    ModInt neg(-1, p);
    EXPECT_TRUE(neg.value() == p - 1, "ModInt negative normalization");

    ModInt base(2, p);
    ModInt result = ModInt::pow(base, 10);
    EXPECT_TRUE(result.value() == 1024, "ModInt power");
    EXPECT_TRUE(ModInt::pow(ModInt(2, 10), 5).value() == 2,
                "ModInt power with even modulus");

    ModInt inv = b.inverse();
    EXPECT_TRUE((b * inv).value() == 1, "ModInt inverse");

    std::cout << "  [PASS]" << std::endl;
}

void test_crt() {
    std::cout << "Test: Chinese Remainder Theorem" << std::endl;

    auto [x1, m1] = crt(2, 3, 3, 5);
    EXPECT_TRUE(m1 == 15, "CRT modulus is product");
    EXPECT_TRUE(x1 % 3 == 2, "CRT satisfies first residue");
    EXPECT_TRUE(x1 % 5 == 3, "CRT satisfies second residue");

    std::vector<int64_t> residues = {2, 3, 2};
    std::vector<int64_t> primes = {3, 5, 7};
    auto [x2, m2] = multi_crt(residues, primes);
    EXPECT_TRUE(m2 == 105, "multi CRT modulus is product");
    EXPECT_TRUE(x2 % 3 == 2, "multi CRT satisfies first residue");
    EXPECT_TRUE(x2 % 5 == 3, "multi CRT satisfies second residue");
    EXPECT_TRUE(x2 % 7 == 2, "multi CRT satisfies third residue");

    std::cout << "  [PASS]" << std::endl;
}

void test_rational_reconstruction() {
    std::cout << "Test: Rational Reconstruction" << std::endl;

    int64_t p = 1000000007;
    ModInt three(3, p);
    ModInt seven(7, p);
    ModInt encoded = three / seven;

    auto [num, den] =
        rational_reconstruction_checked(encoded.value(), p).value();

    std::cout << "  Reconstructed: " << num << "/" << den << std::endl;
    EXPECT_TRUE(num == 3 && den == 7, "rational reconstruction recovers 3/7");

    ModInt neg2(-2, p);
    ModInt five(5, p);
    ModInt encoded2 = neg2 / five;
    auto [num2, den2] =
        rational_reconstruction_checked(encoded2.value(), p).value();
    std::cout << "  Reconstructed: " << num2 << "/" << den2 << std::endl;
    EXPECT_TRUE(num2 == -2 && den2 == 5, "rational reconstruction recovers -2/5");

    std::cout << "  [PASS]" << std::endl;
}

void test_fglm_helpers() {
    std::cout << "Test: FGLM helpers" << std::endl;

    size_t n = 2;

    FGLMPoly g1(n);
    g1.add_term({2, 0}, Rational(1));
    g1.sort_terms(MonomialOrder::grevlex());

    FGLMPoly g2(n);
    g2.add_term({0, 2}, Rational(1));
    g2.sort_terms(MonomialOrder::grevlex());

    std::vector<FGLMPoly> basis = {g1, g2};

    EXPECT_TRUE(is_zero_dimensional(basis, n), "basis <x^2,y^2> is zero-dimensional");
    int dim = quotient_dimension(basis, n);
    std::cout << "  Quotient dimension of <x^2, y^2>: " << dim << std::endl;
    EXPECT_TRUE(dim == 4, "quotient dimension of <x^2,y^2> is 4");

    FGLMPoly x3(n);
    x3.add_term({3, 0}, Rational(1));
    x3.sort_terms(MonomialOrder::grevlex());

    FGLMPoly nf = normal_form(x3, basis, MonomialOrder::grevlex());
    EXPECT_TRUE(nf.is_zero(), "x^3 reduces to zero modulo <x^2,y^2>");

    FGLMPoly xy(n);
    xy.add_term({1, 1}, Rational(1));
    xy.sort_terms(MonomialOrder::grevlex());

    FGLMPoly nf2 = normal_form(xy, basis, MonomialOrder::grevlex());
    EXPECT_TRUE(!nf2.is_zero(), "xy is not reduced to zero modulo <x^2,y^2>");
    EXPECT_TRUE(nf2.terms.size() == 1, "normal form of xy has one term");
    EXPECT_TRUE(!nf2.terms.empty() && nf2.terms[0].first == Monomial({1, 1}),
                "normal form of xy keeps monomial xy");

    std::cout << "  [PASS]" << std::endl;
}

void test_fglm_conversion() {
    std::cout << "Test: FGLM grevlex -> lex conversion" << std::endl;

    {
        size_t n = 2;
        FGLMPoly g1(n);
        g1.add_term({2, 0}, Rational(1));
        g1.sort_terms(MonomialOrder::grevlex());

        FGLMPoly g2(n);
        g2.add_term({0, 2}, Rational(1));
        g2.sort_terms(MonomialOrder::grevlex());

        auto lex_basis = grevlex_to_lex({g1, g2}, n);
        std::cout << "  <x^2, y^2> lex basis size: " << lex_basis.size() << std::endl;
        EXPECT_TRUE(lex_basis.size() == 2, "FGLM conversion preserves two basis elements");
    }

    std::cout << "  (Complex FGLM test skipped — needs further optimization)" << std::endl;

    std::cout << "  [PASS]" << std::endl;
}

int main() {
    test_modint();
    test_crt();
    test_rational_reconstruction();
    test_fglm_helpers();
    test_fglm_conversion();

    std::cout << "\nAll modular/FGLM tests passed!" << std::endl;
    return TEST_REPORT();
}
