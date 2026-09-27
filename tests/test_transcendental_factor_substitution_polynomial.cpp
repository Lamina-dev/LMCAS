#include "test_common.hpp"
#include "transcendental_factor.hpp"
#include "poly_utils.hpp"

using namespace LMCAS;

TEST(TranscendentalFactorSubstitutionPolynomial, BuildPolySingleIndeterminate) {
    auto u0 = SymbolicExpr::variable("u0");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(u0, SymbolicExpr::number(2)),
        SymbolicExpr::add(u0, SymbolicExpr::number(1)));

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_TRUE((result.success)) << "conversion should succeed";
    EXPECT_EQ((result.main_variable), ("u0")) << "main variable should be u0";
    EXPECT_TRUE((result.poly.degree() == 2)) << "polynomial degree should be 2";
    EXPECT_TRUE((result.param_variables.empty())) << "no parameter variables";

    EXPECT_TRUE((result.poly.coeffs.size() == 3)) << "should have 3 coefficients";
    EXPECT_TRUE((result.poly.coeffs[0] == Rational(1))) << "constant term is 1";
    EXPECT_TRUE((result.poly.coeffs[1] == Rational(1))) << "linear term is 1";
    EXPECT_TRUE((result.poly.coeffs[2] == Rational(1))) << "quadratic term is 1";
}

TEST(TranscendentalFactorSubstitutionPolynomial, BuildPolySingleIndeterminateWithRationalCoeffs) {
    auto u0 = SymbolicExpr::variable("u0");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2),
                               SymbolicExpr::power(u0, SymbolicExpr::number(3))),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-3), u0),
            SymbolicExpr::number(5)));

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_TRUE((result.success)) << "conversion should succeed";
    EXPECT_EQ((result.main_variable), ("u0")) << "main variable should be u0";
    EXPECT_TRUE((result.poly.degree() == 3)) << "polynomial degree should be 3";

    EXPECT_TRUE((result.poly.coeffs[0] == Rational(5))) << "constant term is 5";
    EXPECT_TRUE((result.poly.coeffs[1] == Rational(-3))) << "linear term is -3";
    EXPECT_TRUE((result.poly.coeffs[2] == Rational(0))) << "quadratic term is 0";
    EXPECT_TRUE((result.poly.coeffs[3] == Rational(2))) << "cubic term is 2";
}

TEST(TranscendentalFactorSubstitutionPolynomial, BuildPolyConstantExpression) {
    auto expr = SymbolicExpr::number(42);

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_TRUE((result.success)) << "conversion should succeed for constant";
    EXPECT_TRUE((result.poly.degree() == 0)) << "constant polynomial has degree 0";
    EXPECT_TRUE((result.poly.coeffs[0] == Rational(42))) << "constant value is 42";
}

TEST(TranscendentalFactorSubstitutionPolynomial, BuildPolyNullExpression) {
    std::vector<std::string> indeterminates = {"u0"};
    auto result = tf_build_polynomial(nullptr, indeterminates, "x");

    EXPECT_TRUE((!result && result.error().code == CasErrc::InvalidArgument)) << "null expression reports InvalidArgument";
}

TEST(TranscendentalFactorSubstitutionPolynomial, BuildPolyOnlyX) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(3), x),
            SymbolicExpr::number(2)));

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_TRUE((result.success)) << "conversion should succeed";
    EXPECT_EQ((result.main_variable), ("x")) << "main variable should be x";
    EXPECT_TRUE((result.poly.degree() == 2)) << "polynomial degree should be 2";
}

TEST(TranscendentalFactorSubstitutionPolynomial, BuildPolyMultipleIndeterminatesRequireRationalCoefficients) {
    auto u0 = SymbolicExpr::variable("u0");
    auto u1 = SymbolicExpr::variable("u1");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(u0, SymbolicExpr::number(3)),
        SymbolicExpr::power(u1, SymbolicExpr::number(2)));

    std::vector<std::string> indeterminates = {"u0", "u1"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_FALSE((result.success)) << "no variable choice may silently discard the other indeterminate";
}

TEST(TranscendentalFactorSubstitutionPolynomial, BuildPolyFromSubstitutionResult) {
    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        SymbolicExpr::add(sin_x, SymbolicExpr::number(1)));

    auto sub_result = detect_trans_substitutions(expr, "x");
    EXPECT_TRUE((sub_result.mappings.size() == 1)) << "should have 1 mapping";

    std::vector<std::string> indeterminates;
    for (const auto &m : sub_result.mappings) {
        indeterminates.push_back(m.indeterminate);
    }

    auto checked = tf_build_polynomial(sub_result.poly_expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto poly_result = std::move(checked.value());

    EXPECT_TRUE((poly_result.success)) << "polynomial construction should succeed";
    EXPECT_EQ((poly_result.main_variable), ("u0")) << "main variable should be u0";
    EXPECT_TRUE((poly_result.poly.degree() == 2)) << "polynomial degree should be 2";
    EXPECT_TRUE((poly_result.poly.coeffs[0] == Rational(1))) << "constant term is 1";
    EXPECT_TRUE((poly_result.poly.coeffs[1] == Rational(1))) << "linear term is 1";
    EXPECT_TRUE((poly_result.poly.coeffs[2] == Rational(1))) << "quadratic term is 1";
}

TEST(TranscendentalFactorSubstitutionPolynomial, BuildPolyNonPolynomialFails) {
    auto u0 = SymbolicExpr::variable("u0");
    auto expr = SymbolicExpr::sin(u0);

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_FALSE((result.success)) << "non-polynomial expression should fail";
}

TEST(TranscendentalFactorSubstitutionPolynomial, ValidateValidPolynomial) {
    auto u0 = SymbolicExpr::variable("u0");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(u0, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(3), u0),
            SymbolicExpr::number(-2)));

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_TRUE((result.success)) << "valid polynomial should pass validation";
    EXPECT_TRUE((result.poly.degree() == 2)) << "degree should be 2";
}

TEST(TranscendentalFactorSubstitutionPolynomial, ValidateRemainingSinFails) {
    auto u0 = SymbolicExpr::variable("u0");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(u0), u0);

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_FALSE((result.success)) << "expression with remaining sin(u0) should fail validation";
}

TEST(TranscendentalFactorSubstitutionPolynomial, ValidateRemainingCosFails) {
    auto u0 = SymbolicExpr::variable("u0");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(u0, SymbolicExpr::number(2)),
        SymbolicExpr::cos(u0));

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_FALSE((result.success)) << "expression with remaining cos(u0) should fail validation";
}

TEST(TranscendentalFactorSubstitutionPolynomial, ValidateRemainingExpFails) {
    auto u0 = SymbolicExpr::variable("u0");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(u0),
        SymbolicExpr::power(u0, SymbolicExpr::number(2)));

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_FALSE((result.success)) << "expression with remaining exp(u0) should fail validation";
}

TEST(TranscendentalFactorSubstitutionPolynomial, ValidateFractionalExponentFails) {
    auto u0 = SymbolicExpr::variable("u0");
    auto half = SymbolicExpr::number(Rational(BigInt(1), BigInt(2)));
    auto expr = SymbolicExpr::power(u0, half);

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_FALSE((result.success)) << "fractional exponent should fail validation";
}

TEST(TranscendentalFactorSubstitutionPolynomial, ValidateNegativeExponentFails) {
    auto u0 = SymbolicExpr::variable("u0");
    auto expr = SymbolicExpr::power(u0, SymbolicExpr::number(-1));

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_FALSE((result.success)) << "negative exponent should fail validation";
}

TEST(TranscendentalFactorSubstitutionPolynomial, ValidateNegativeFractionalExponentFails) {
    auto u0 = SymbolicExpr::variable("u0");
    auto neg_three_half = SymbolicExpr::number(Rational(BigInt(-3), BigInt(2)));
    auto expr = SymbolicExpr::power(u0, neg_three_half);

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_FALSE((result.success)) << "negative fractional exponent should fail validation";
}

TEST(TranscendentalFactorSubstitutionPolynomial, ValidateTranscendentalInOriginalVarFails) {
    auto x = SymbolicExpr::variable("x");
    auto u0 = SymbolicExpr::variable("u0");
    auto expr = SymbolicExpr::add(u0, SymbolicExpr::sin(x));

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_FALSE((result.success)) << "expression with remaining sin(x) should fail validation";
}

TEST(TranscendentalFactorSubstitutionPolynomial, BuildPolyMultivariateRepresentationBoundary) {
    auto u0 = SymbolicExpr::variable("u0");
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(u0, SymbolicExpr::number(2)),
        x);

    std::vector<std::string> indeterminates = {"u0"};
    auto checked = tf_build_polynomial(expr, indeterminates, "x");
    ASSERT_TRUE((checked.has_value())) << "polynomial strategy returns without a hard error";
    if (!checked)
        return;
    auto result = std::move(checked.value());

    EXPECT_FALSE((result.success)) << "a rational univariate result cannot represent the symbolic coefficient";

    const auto symbolic = symbolic_to_poly<SymbolicPolyCoeff>(expr, "u0");
    ASSERT_TRUE((symbolic.has_value())) << "multivariate expression retains a symbolic representation";
    if (symbolic) {
        const auto restored = poly_to_symbolic(symbolic.value());
        const auto value = test_numeric_eval(
            restored->substitute("u0", SymbolicExpr::number(3))
                ->substitute("x", SymbolicExpr::number(7)));
        EXPECT_TRUE((value && *value == 16.0)) << "symbolic reconstruction retains both the square and parameter term";
    }

    const auto specialized = tf_build_polynomial(
        expr->substitute("x", SymbolicExpr::number(7)), indeterminates, "x");
    EXPECT_TRUE((specialized && specialized.value().success)) << "a rational specialization is supported";
    if (specialized && specialized.value().success) {
        EXPECT_TRUE((specialized.value().poly.eval(Rational(3)) == Rational(16))) << "rational specialization includes the former symbolic coefficient";
    }
}

TEST(TranscendentalFactorSubstitutionPolynomial, BuildPolyResourceLimit) {
    const auto u = SymbolicExpr::variable("u0");
    const auto expression = SymbolicExpr::add(
        SymbolicExpr::power(u, SymbolicExpr::number(1000)), SymbolicExpr::number(1));
    const auto result = tf_build_polynomial(expression, {"u0"}, "x");
    EXPECT_TRUE((!result && result.error().code == CasErrc::ResourceLimit)) << "excessive power propagates ResourceLimit rather than success=false";
}
