#include "test_common.hpp"
#include "transcendental_factor.hpp"

using namespace LMCAS;

TEST(TranscendentalFactorSubstitution, SingleSin) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.size() == 1)) << "should find 1 transcendental sub-expression";
    EXPECT_EQ((result.mappings[0].indeterminate), ("u0")) << "first indeterminate is u0";
    EXPECT_TRUE(test_expression_text((result.mappings[0].trans_expr), ("sin(x)"))) << "trans_expr is sin(x)";
}

TEST(TranscendentalFactorSubstitution, SinAndCos) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), SymbolicExpr::cos(x));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.size() == 2)) << "should find 2 transcendental sub-expressions";
    EXPECT_EQ((result.mappings[0].indeterminate), ("u0")) << "first indeterminate is u0";
    EXPECT_EQ((result.mappings[1].indeterminate), ("u1")) << "second indeterminate is u1";
}

TEST(TranscendentalFactorSubstitution, Deduplication) {
    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        sin_x);

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.size() == 1)) << "should deduplicate to 1 transcendental sub-expression";
    EXPECT_EQ((result.mappings[0].indeterminate), ("u0")) << "indeterminate is u0";
}

TEST(TranscendentalFactorSubstitution, ExpFunction) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::exp(x),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(x, SymbolicExpr::number(2))));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.size() == 1)) << "should find 1 transcendental sub-expression (exp)";
    EXPECT_EQ((result.mappings[0].indeterminate), ("u0")) << "indeterminate is u0";
    EXPECT_TRUE(test_expression_text((result.mappings[0].trans_expr), ("exp(x)"))) << "trans_expr is exp(x)";
}

TEST(TranscendentalFactorSubstitution, LnFunction) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::ln(x), x);

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.size() == 1)) << "should find 1 transcendental sub-expression (ln)";
    EXPECT_EQ((result.mappings[0].indeterminate), ("u0")) << "indeterminate is u0";
    EXPECT_TRUE(test_expression_text((result.mappings[0].trans_expr), ("ln(x)"))) << "trans_expr is ln(x)";
}

TEST(TranscendentalFactorSubstitution, TanFunction) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::tan(x), SymbolicExpr::number(1));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.size() == 1)) << "should find 1 transcendental sub-expression (tan)";
    EXPECT_EQ((result.mappings[0].indeterminate), ("u0")) << "indeterminate is u0";
    EXPECT_TRUE(test_expression_text((result.mappings[0].trans_expr), ("tan(x)"))) << "trans_expr is tan(x)";
}

TEST(TranscendentalFactorSubstitution, NoTranscendental) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::number(1)));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.empty())) << "should find no transcendental sub-expressions";
}

TEST(TranscendentalFactorSubstitution, IndependentVariable) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(y), x);

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.empty())) << "sin(y) should not be collected when var=x";
}

TEST(TranscendentalFactorSubstitution, NestedTranscendental) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(SymbolicExpr::exp(x)), x);

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.size() == 2)) << "should find 2 transcendental sub-expressions (outer and inner)";
    EXPECT_EQ((result.mappings[0].indeterminate), ("u0")) << "first indeterminate is u0";
    EXPECT_EQ((result.mappings[1].indeterminate), ("u1")) << "second indeterminate is u1";
}

TEST(TranscendentalFactorSubstitution, MultipleDistinct) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::add(SymbolicExpr::exp(x), SymbolicExpr::cos(x)));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.size() == 3)) << "should find 3 distinct transcendental sub-expressions";
}

TEST(TranscendentalFactorSubstitution, NullExpression) {
    auto result = detect_trans_substitutions(nullptr, "x");

    EXPECT_TRUE((result.mappings.empty())) << "null expression should return empty mappings";
    EXPECT_TRUE((result.poly_expr == nullptr)) << "poly_expr should be null";
}

TEST(TranscendentalFactorSubstitution, SinCosPythagoreanConstraint) {
    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::power(cos_x, SymbolicExpr::number(2)),
            SymbolicExpr::number(-1)));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.size() == 2)) << "should find 2 transcendental sub-expressions";
    EXPECT_TRUE((result.constraints.size() == 1)) << "should detect 1 Pythagorean constraint";

    if (!result.constraints.empty()) {
        std::string cstr = result.constraints[0]->to_string();
        EXPECT_TRUE((cstr.find("u0") != std::string::npos)) << "constraint contains u0";
        EXPECT_TRUE((cstr.find("u1") != std::string::npos)) << "constraint contains u1";
    }
}

TEST(TranscendentalFactorSubstitution, ExpInverseConstraint) {
    auto x = SymbolicExpr::variable("x");
    auto exp_x = SymbolicExpr::exp(x);
    auto neg_x = SymbolicExpr::multiply(SymbolicExpr::number(-1), x);
    auto exp_neg_x = SymbolicExpr::exp(neg_x);
    auto expr = SymbolicExpr::add(exp_x, exp_neg_x);

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.size() == 2)) << "should find 2 transcendental sub-expressions (exp(x) and exp(-x))";
    EXPECT_TRUE((result.constraints.size() == 1)) << "should detect 1 inverse constraint";

    if (!result.constraints.empty()) {
        std::string cstr = result.constraints[0]->to_string();
        EXPECT_TRUE((cstr.find("u0") != std::string::npos)) << "constraint contains u0";
        EXPECT_TRUE((cstr.find("u1") != std::string::npos)) << "constraint contains u1";
    }
}

TEST(TranscendentalFactorSubstitution, NoConstraintDifferentArgs) {
    auto x = SymbolicExpr::variable("x");
    auto two_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), SymbolicExpr::cos(two_x));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.size() == 2)) << "should find 2 transcendental sub-expressions";
    EXPECT_TRUE((result.constraints.empty())) << "no constraint when arguments differ";
}

TEST(TranscendentalFactorSubstitution, NoConstraintSameType) {
    auto x = SymbolicExpr::variable("x");
    auto two_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
    auto expr = SymbolicExpr::add(SymbolicExpr::exp(x), SymbolicExpr::exp(two_x));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.size() == 2)) << "should find 2 transcendental sub-expressions";
    EXPECT_TRUE((result.constraints.empty())) << "no constraint when exp arguments are not negations";
}

TEST(TranscendentalFactorSubstitution, MultipleConstraints) {
    auto x = SymbolicExpr::variable("x");
    auto neg_x = SymbolicExpr::multiply(SymbolicExpr::number(-1), x);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::add(
            SymbolicExpr::cos(x),
            SymbolicExpr::add(
                SymbolicExpr::exp(x),
                SymbolicExpr::exp(neg_x))));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.mappings.size() == 4)) << "should find 4 transcendental sub-expressions";
    EXPECT_TRUE((result.constraints.size() == 2)) << "should detect 2 constraints (Pythagorean + inverse)";
}

TEST(TranscendentalFactorSubstitution, SubstitutionSinPlusX) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.poly_expr != nullptr)) << "poly_expr should not be null";
    std::string poly_str = result.poly_expr->to_string();
    EXPECT_TRUE((poly_str.find("u0") != std::string::npos)) << "poly_expr should contain u0";
    EXPECT_TRUE((poly_str.find("sin") == std::string::npos)) << "poly_expr should not contain sin";
    EXPECT_TRUE((poly_str.find("x") != std::string::npos)) << "poly_expr should still contain x";
}

TEST(TranscendentalFactorSubstitution, SubstitutionSinSquaredMinusXSquared) {
    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(x, SymbolicExpr::number(2))));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.poly_expr != nullptr)) << "poly_expr should not be null";
    std::string poly_str = result.poly_expr->to_string();
    EXPECT_TRUE((poly_str.find("u0") != std::string::npos)) << "poly_expr should contain u0";
    EXPECT_TRUE((poly_str.find("sin") == std::string::npos)) << "poly_expr should not contain sin";
}

TEST(TranscendentalFactorSubstitution, SubstitutionNoTranscendental) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::number(1)));

    auto result = detect_trans_substitutions(expr, "x");

    EXPECT_TRUE((result.poly_expr != nullptr)) << "poly_expr should not be null for pure polynomial";
    std::string orig_str = expr->to_string();
    std::string poly_str = result.poly_expr->to_string();
    EXPECT_EQ((poly_str), (orig_str)) << "poly_expr should equal original expression";
}
