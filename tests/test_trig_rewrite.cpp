#include "symbolic.hpp"
#include "test_common.hpp"

using namespace LMCAS;

TEST(LmcasTrigRewrite, SinSqAddCosSq) {
    auto x = SymbolicExpr::variable("x");
    auto sinx = SymbolicExpr::sin(x);
    auto cosx = SymbolicExpr::cos(x);
    auto sin2 = SymbolicExpr::power(sinx, SymbolicExpr::number(2));
    auto cos2 = SymbolicExpr::power(cosx, SymbolicExpr::number(2));

    auto expr = SymbolicExpr::add(sin2, cos2);
    auto simplified = expr->simplify_trig();
    if (!simplified || !detail::node(simplified)) {
        FAIL() << "Scenario requirement failed";
    }

    auto expected = SymbolicExpr::number(1)->simplify();
    bool result = simplified->simplify()->compare(expected) == 0;
    for (double value : {0.0, 0.3, -0.7}) {
        auto numeric = test_numeric_eval(simplified->substitute("x", SymbolicExpr::number(value)));
        const double sine = std::sin(value);
        const double cosine = std::cos(value);
        ASSERT_TRUE(numeric.has_value());
        ASSERT_TRUE(std::isfinite(*numeric));
        EXPECT_NEAR(*numeric, sine * sine + cosine * cosine, 1e-12);
    }
    EXPECT_TRUE((result));
}

TEST(LmcasTrigRewrite, Sin2x) {
    auto x = SymbolicExpr::variable("x");
    auto two_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
    auto sin2x = SymbolicExpr::sin(two_x);

    auto simplified = sin2x->simplify_trig();
    if (!simplified || !detail::node(simplified)) {
        FAIL() << "Scenario requirement failed";
    }

    auto expected = SymbolicExpr::multiply(SymbolicExpr::number(2),
                                           SymbolicExpr::multiply(SymbolicExpr::sin(x), SymbolicExpr::cos(x)))
                        ->simplify();
    bool result = simplified->simplify()->compare(expected) == 0;
    for (double value : {0.0, 0.3, -0.7}) {
        auto numeric = test_numeric_eval(simplified->substitute("x", SymbolicExpr::number(value)));
        ASSERT_TRUE(numeric.has_value());
        ASSERT_TRUE(std::isfinite(*numeric));
        EXPECT_NEAR(*numeric, std::sin(2.0 * value), 1e-12);
    }
    EXPECT_TRUE((result));
}

TEST(LmcasTrigRewrite, Cos2x) {
    auto x = SymbolicExpr::variable("x");
    auto two_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
    auto cos2x = SymbolicExpr::cos(two_x);

    auto simplified = cos2x->simplify_trig();
    if (!simplified || !detail::node(simplified)) {
        FAIL() << "Scenario requirement failed";
    }

    auto two = SymbolicExpr::number(2);
    auto expected = SymbolicExpr::add(
                        SymbolicExpr::power(SymbolicExpr::cos(x), two),
                        SymbolicExpr::multiply(SymbolicExpr::number(-1),
                                               SymbolicExpr::power(SymbolicExpr::sin(x), two)))
                        ->simplify();
    bool result = simplified->simplify()->compare(expected) == 0;
    for (double value : {0.0, 0.3, -0.7}) {
        auto numeric = test_numeric_eval(simplified->substitute("x", SymbolicExpr::number(value)));
        ASSERT_TRUE(numeric.has_value());
        ASSERT_TRUE(std::isfinite(*numeric));
        EXPECT_NEAR(*numeric, std::cos(2.0 * value), 1e-12);
    }
    EXPECT_TRUE((result));
}
