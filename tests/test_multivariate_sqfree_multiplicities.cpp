#include "test_multivariate_sqfree_support.hpp"

TEST(MultivariateSqfreeMultiplicities, UnitXY2XYComponent0XYComponent1XY) {
    std::vector<std::string> vars = {"x", "y"};

    std::vector<MultiPoly::Term> xpy_terms = {
        make_term({1, 0}, Rational(1)),
        make_term({0, 1}, Rational(1))};
    MultiPoly xpy(xpy_terms, vars);

    std::vector<MultiPoly::Term> xmy_terms = {
        make_term({1, 0}, Rational(1)),
        make_term({0, 1}, Rational(-1))};
    MultiPoly xmy(xmy_terms, vars);

    MultiPoly f = xpy * xpy * xmy;
    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    EXPECT_TRUE((decomp.components.size() >= 2)) << "(x+y)^2*(x-y): at least 2 components";

    MultiPoly comp0_prim = decomp.components[0].make_primitive();
    MultiPoly xmy_prim = xmy.make_primitive();
    EXPECT_TRUE((comp0_prim == xmy_prim || comp0_prim == (-xmy).make_primitive())) << "(x+y)^2*(x-y): component[0] == (x-y) up to sign";

    MultiPoly comp1_prim = decomp.components[1].make_primitive();
    MultiPoly xpy_prim = xpy.make_primitive();
    EXPECT_TRUE((comp1_prim == xpy_prim || comp1_prim == (-xpy).make_primitive())) << "(x+y)^2*(x-y): component[1] == (x+y) up to sign";

    MultiPoly reconstructed = decomp.components[0] * poly_pow(decomp.components[1], 2);
    MultiPoly f_prim = f.make_primitive();
    MultiPoly recon_prim = reconstructed.make_primitive();
    EXPECT_TRUE((recon_prim == f_prim || recon_prim == (-f).make_primitive())) << "(x+y)^2*(x-y): reconstruction matches original";
}

TEST(MultivariateSqfreeMultiplicities, UnitAlreadySquareFreeX2Y2ReturnsSingleComponent) {
    std::vector<std::string> vars = {"x", "y"};
    std::vector<MultiPoly::Term> terms = {
        make_term({2, 0}, Rational(1)),
        make_term({0, 2}, Rational(-1))};
    MultiPoly f(terms, vars);

    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    EXPECT_TRUE((decomp.components.size() == 1)) << "x^2-y^2: exactly 1 component (already square-free)";

    MultiPoly comp_prim = decomp.components[0].make_primitive();
    MultiPoly f_prim = f.make_primitive();
    EXPECT_TRUE((comp_prim == f_prim || comp_prim == (-f).make_primitive())) << "x^2-y^2: component[0] == f up to sign";
}

TEST(MultivariateSqfreeMultiplicities, UnitUnivariateX12X1Component0X1Component1X1) {
    std::vector<std::string> vars = {"x"};

    std::vector<MultiPoly::Term> xm1_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(-1))};
    MultiPoly xm1(xm1_terms, vars);

    std::vector<MultiPoly::Term> xp1_terms = {
        make_term({1}, Rational(1)),
        make_term({0}, Rational(1))};
    MultiPoly xp1(xp1_terms, vars);

    MultiPoly f = xm1 * xm1 * xp1;
    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    EXPECT_TRUE((decomp.components.size() >= 2)) << "(x-1)^2*(x+1): at least 2 components";

    MultiPoly comp0_prim = decomp.components[0].make_primitive();
    MultiPoly xp1_prim = xp1.make_primitive();
    EXPECT_TRUE((comp0_prim == xp1_prim || comp0_prim == (-xp1).make_primitive())) << "(x-1)^2*(x+1): component[0] == (x+1) up to sign";

    MultiPoly comp1_prim = decomp.components[1].make_primitive();
    MultiPoly xm1_prim = xm1.make_primitive();
    EXPECT_TRUE((comp1_prim == xm1_prim || comp1_prim == (-xm1).make_primitive())) << "(x-1)^2*(x+1): component[1] == (x-1) up to sign";

    MultiPoly reconstructed = decomp.components[0] * poly_pow(decomp.components[1], 2);
    MultiPoly f_prim = f.make_primitive();
    MultiPoly recon_prim = reconstructed.make_primitive();
    EXPECT_TRUE((recon_prim == f_prim || recon_prim == (-f).make_primitive())) << "(x-1)^2*(x+1): reconstruction matches original";
}

TEST(MultivariateSqfreeMultiplicities, UnitPureSquareXY2HasConstantComponent0XYAsComponent1) {
    std::vector<std::string> vars = {"x", "y"};

    std::vector<MultiPoly::Term> xpy_terms = {
        make_term({1, 0}, Rational(1)),
        make_term({0, 1}, Rational(1))};
    MultiPoly xpy(xpy_terms, vars);

    MultiPoly f = xpy * xpy;
    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    EXPECT_TRUE((decomp.components.size() >= 2)) << "(x+y)^2: at least 2 components";

    EXPECT_TRUE((decomp.components[0].is_constant())) << "(x+y)^2: component[0] is constant";

    MultiPoly comp1_prim = decomp.components[1].make_primitive();
    MultiPoly xpy_prim = xpy.make_primitive();
    EXPECT_TRUE((comp1_prim == xpy_prim || comp1_prim == (-xpy).make_primitive())) << "(x+y)^2: component[1] == (x+y) up to sign";

    MultiPoly reconstructed = decomp.components[0] * poly_pow(decomp.components[1], 2);
    MultiPoly f_prim = f.make_primitive();
    MultiPoly recon_prim = reconstructed.make_primitive();
    EXPECT_TRUE((recon_prim == f_prim || recon_prim == (-f).make_primitive())) << "(x+y)^2: reconstruction matches original";
}

TEST(MultivariateSqfreeMultiplicities, UnitZeroPolynomialDecomposition) {
    MultiPoly zero;
    SquareFreeDecomp decomp = square_free_decompose(zero, "x");

    EXPECT_TRUE((decomp.components.size() == 1)) << "zero: exactly 1 component";
    EXPECT_TRUE((decomp.components[0].is_zero())) << "zero: component is zero polynomial";
}

TEST(MultivariateSqfreeMultiplicities, UnitConstantPolynomialDecomposition) {
    std::vector<std::string> vars = {"x", "y"};
    MultiPoly f(Rational(42), vars);

    SquareFreeDecomp decomp = square_free_decompose(f, "x");

    EXPECT_TRUE((decomp.components.size() == 1)) << "constant 42: exactly 1 component";
    EXPECT_TRUE((decomp.components[0].is_constant())) << "constant 42: component is constant";
}
