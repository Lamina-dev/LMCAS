#include "expr.hpp"
#include "symbolic.hpp"
#include "test_common.hpp"

using namespace LMCAS;

namespace {

void expect_resultant(const char* left_source,
                      const char* right_source,
                      const char* expected_source) {
    auto left = parse_expr(left_source);
    auto right = parse_expr(right_source);
    auto expected = parse_expr(expected_source);
    ASSERT_TRUE(left);
    ASSERT_TRUE(right);
    ASSERT_TRUE(expected);

    auto resultant = SymbolicExpr::poly_resultant(
        left.value(), right.value(), "x");
    ASSERT_TRUE(resultant);
    EXPECT_TRUE(test_proved_equivalent(resultant, expected.value()));
}

} // namespace

TEST(PolynomialResultant, SharedLinearFactorProducesZero) {
    expect_resultant("x^2-1", "x+1", "0");
}

TEST(PolynomialResultant, LinearFactorOutsideRootsProducesEvaluation) {
    expect_resultant("x^2-1", "x-2", "3");
}

TEST(PolynomialResultant, DisjointQuadraticsProduceNine) {
    expect_resultant("x^2-1", "x^2-4", "9");
}

TEST(PolynomialResultant, SharedQuadraticRootProducesZero) {
    expect_resultant("x^2-1", "x^2+x-2", "0");
}

TEST(PolynomialResultant, CubicAndQuadraticProduceTwo) {
    expect_resultant("x^3-1", "x^2+1", "2");
}
