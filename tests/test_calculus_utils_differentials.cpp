/**
 * @file test_calculus_utils_differentials.cpp
 * @brief 微分与对数求导约定。
 */

#include "test_common.hpp"
#include "calculus_utils.hpp"
#include <string>

using namespace LMCAS;

using SE = SymbolicExpr;

static auto num(int n) { return SE::number(n); }
static auto var(const std::string &name) { return SE::variable(name); }

TEST(CalculusUtilsDifferentials, LogPolynomialDerivative) {
    auto x = var("x");
    auto f = SE::power(x, num(2));
    auto result = LMCAS::log_differentiate(f, "x");
    EXPECT_TRUE((result != nullptr)) << "log_differentiate(x^2) non-null";

    auto at3 = result->substitute("x", num(3))->simplify();
    auto val = test_numeric_eval(at3);
    if (val) {
        EXPECT_NEAR(*val, 6.0, 1e-9) << "log_differentiate(x^2) at x=3 = 6";
    } else {
        ADD_FAILURE() << "log_differentiate(x^2) at x=3 numeric eval failed";
    }
}

TEST(CalculusUtilsDifferentials, LogVariableExponent) {
    auto x = var("x");
    auto f = SE::power(x, x);
    auto result = LMCAS::log_differentiate(f, "x");
    EXPECT_TRUE((result != nullptr)) << "log_differentiate(x^x) non-null";

    auto at1 = result->substitute("x", num(1))->simplify();
    auto val = test_numeric_eval(at1);
    if (val) {
        EXPECT_NEAR(*val, 1.0, 1e-9) << "log_differentiate(x^x) at x=1 = 1";
    } else {
        ADD_FAILURE() << "log_differentiate(x^x) at x=1 numeric eval failed";
    }
}

TEST(CalculusUtilsDifferentials, LogProductDerivative) {
    auto x = var("x");
    auto f = SE::multiply(x, SE::add(x, num(1)));
    auto result = LMCAS::log_differentiate(f, "x");
    EXPECT_TRUE((result != nullptr)) << "log_differentiate(x*(x+1)) non-null";

    auto at2 = result->substitute("x", num(2))->simplify();
    auto val = test_numeric_eval(at2);
    if (val) {
        EXPECT_NEAR(*val, 5.0, 1e-9) << "log_differentiate(x*(x+1)) at x=2 = 5";
    } else {
        ADD_FAILURE() << "log_differentiate(x*(x+1)) at x=2 numeric eval failed";
    }
}

TEST(CalculusUtilsDifferentials, PolynomialDifferential) {
    auto x = var("x");
    auto f = SE::power(x, num(3));
    auto result = LMCAS::differential(f, "x");
    EXPECT_TRUE((result != nullptr)) << "differential(x^3) non-null";

    auto at2 = result->substitute("x", num(2))->simplify();
    auto val = test_numeric_eval(at2);
    if (val) {
        EXPECT_NEAR(*val, 12.0, 1e-9) << "differential(x^3) at x=2 = 12";
    } else {
        ADD_FAILURE() << "differential(x^3) at x=2 numeric eval failed";
    }
}

TEST(CalculusUtilsDifferentials, SineDifferential) {
    auto x = var("x");
    auto result = LMCAS::differential(SE::sin(x), "x");
    ASSERT_NE(result, nullptr);
    EXPECT_TRUE(test_proved_equivalent(result, SE::cos(x)));
}

TEST(CalculusUtilsDifferentials, QuadraticTotalDifferential) {
    auto x = var("x");
    auto y = var("y");
    auto f = SE::add(SE::power(x, num(2)), SE::power(y, num(2)));
    auto result = LMCAS::total_differential(f, {"x", "y"});
    EXPECT_TRUE((result.size() == 2)) << "total_differential has 2 terms";

    EXPECT_EQ((result[0].second), ("x")) << "first term variable is x";

    EXPECT_EQ((result[1].second), ("y")) << "second term variable is y";

    auto dx_at = result[0].first->substitute("x", num(3))->substitute("y", num(4))->simplify();
    auto val_dx = test_numeric_eval(dx_at);
    if (val_dx) {
        EXPECT_NEAR(*val_dx, 6.0, 1e-9) << "df/dx at (3,4) = 6";
    }

    auto dy_at = result[1].first->substitute("x", num(3))->substitute("y", num(4))->simplify();
    auto val_dy = test_numeric_eval(dy_at);
    if (val_dy) {
        EXPECT_NEAR(*val_dy, 8.0, 1e-9) << "df/dy at (3,4) = 8";
    }
}

TEST(CalculusUtilsDifferentials, ProductTotalDifferential) {
    auto x = var("x");
    auto y = var("y");
    auto z = var("z");
    auto f = SE::multiply(SE::multiply(x, y), z);
    auto result = LMCAS::total_differential(f, {"x", "y", "z"});
    EXPECT_TRUE((result.size() == 3)) << "total_differential has 3 terms";
    EXPECT_EQ((result[0].second), ("x")) << "first var is x";
    EXPECT_EQ((result[1].second), ("y")) << "second var is y";
    EXPECT_EQ((result[2].second), ("z")) << "third var is z";
}

TEST(CalculusUtilsDifferentials, LogDerivativeNull) {
    auto result = LMCAS::log_differentiate(nullptr, "x");
    EXPECT_TRUE((result == nullptr)) << "log_differentiate(null) returns null";
}

TEST(CalculusUtilsDifferentials, DifferentialNull) {
    auto result = LMCAS::differential(nullptr, "x");
    EXPECT_TRUE((result == nullptr)) << "differential(null) returns null";
}

TEST(CalculusUtilsDifferentials, TotalDifferentialEmpty) {
    auto x = var("x");
    auto result = LMCAS::total_differential(x, {});
    EXPECT_TRUE((result.empty())) << "total_differential with empty vars returns empty";
}
