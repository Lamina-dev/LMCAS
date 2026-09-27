#include "limit_result.hpp"
/**
 * @file test_exponential_limit.cpp
 * @brief Regression coverage for exponential limit normalization.
 */
#include "test_common.hpp"
#include "internal/visitors/limit_visitor.hpp"
#include "internal/visitors/differentiation_visitor.hpp"

using namespace LMCAS;

TEST(LmcasExponentialLimit, LogarithmicProductAtInfinity) {
    auto x = SymbolicExpr::variable("x");
    auto one = SymbolicExpr::number(1);
    auto neg_one = SymbolicExpr::number(-1);
    auto inf = SymbolicExpr::infinity(1);

    auto inv_x = SymbolicExpr::power(x, neg_one);
    auto base = SymbolicExpr::add(one, inv_x);
    auto ln_base = SymbolicExpr::ln(base);
    auto product = SymbolicExpr::multiply(x, ln_base);
    auto lim = LMCAS::limit_expression_checked(product, "x", inf).value();
    if (lim) {
        auto val = test_numeric_eval(lim);
        if (val) {
            {
                const double actual_value = *val;
                const double expected_value = 1.0;
                EXPECT_TRUE(std::isfinite(actual_value));
                EXPECT_NEAR(actual_value, expected_value, 1e-6);
            }
        } else {
            EXPECT_TRUE(test_expression_text(lim, "1")) << "lim x*ln(1+1/x) = 1";
        }
    } else {
        ADD_FAILURE() << "limit is null";
    }
}

TEST(LmcasExponentialLimit, LogarithmicDerivative) {
    auto x = SymbolicExpr::variable("x");
    auto one = SymbolicExpr::number(1);
    auto neg_one = SymbolicExpr::number(-1);

    auto inv_x = SymbolicExpr::power(x, neg_one);
    auto base = SymbolicExpr::add(one, inv_x);
    auto ln_base = SymbolicExpr::ln(base);
    auto deriv = ln_base->differentiate("x");
    EXPECT_TRUE((deriv != nullptr)) << "derivative exists";
}

TEST(LmcasExponentialLimit, ReciprocalDerivative) {
    auto x = SymbolicExpr::variable("x");
    auto neg_one = SymbolicExpr::number(-1);

    auto inv_x = SymbolicExpr::power(x, neg_one);
    auto deriv = inv_x->differentiate("x");
    EXPECT_TRUE((deriv != nullptr)) << "derivative exists";
}

TEST(LmcasExponentialLimit, ExponentialPowerAtInfinity) {
    auto x = SymbolicExpr::variable("x");
    auto one = SymbolicExpr::number(1);
    auto neg_one = SymbolicExpr::number(-1);
    auto inf = SymbolicExpr::infinity(1);

    auto inv_x = SymbolicExpr::power(x, neg_one);
    auto base = SymbolicExpr::add(one, inv_x);
    auto expr = SymbolicExpr::power(base, x);
    auto lim = LMCAS::limit_expression_checked(expr, "x", inf).value();
    EXPECT_TRUE((lim != nullptr)) << "limit is not null";
}
