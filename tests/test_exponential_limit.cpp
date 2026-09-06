/**
 * @file test_exponential_limit.cpp
 * @brief Regression coverage for exponential limit normalization.
 */
#include "test_common.hpp"
#include "visitors/limit_visitor.hpp"
#include "visitors/differentiation_visitor.hpp"

using namespace LMCAS;

int main() {
    auto x = SymbolicExpr::variable("x");
    auto one = SymbolicExpr::number(1);
    auto neg_one = SymbolicExpr::number(-1);
    auto inf = SymbolicExpr::infinity(1);

    // Test: lim(x->inf) x * ln(1+1/x) should be 1
    TEST_CASE("lim(x->inf) x*ln(1+1/x) = 1");
    {
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto base = SymbolicExpr::add(one, inv_x);
        auto ln_base = SymbolicExpr::ln(base);
        auto product = SymbolicExpr::multiply(x, ln_base);
        auto lim = LMCAS::limit_expression_checked(product, "x", inf).value();
        if (lim) {
            auto val = test_numeric_eval(lim);
            if (val) {
                EXPECT_NEAR(*val, 1.0, 1e-6, "lim x*ln(1+1/x) = 1");
            } else {
                EXPECT_EQ_EXPR_STR(lim, "1", "lim x*ln(1+1/x) = 1");
            }
        } else {
            EXPECT_TRUE(false, "limit is null");
        }
    }

    // Test: differentiate ln(1+1/x)
    TEST_CASE("d/dx[ln(1+1/x)] exists");
    {
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto base = SymbolicExpr::add(one, inv_x);
        auto ln_base = SymbolicExpr::ln(base);
        auto deriv = ln_base->differentiate("x");
        EXPECT_TRUE(deriv != nullptr, "derivative exists");
    }

    // Test: differentiate x^(-1)
    TEST_CASE("d/dx[x^(-1)] exists");
    {
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto deriv = inv_x->differentiate("x");
        EXPECT_TRUE(deriv != nullptr, "derivative exists");
    }

    // Test: the full (1+1/x)^x limit
    TEST_CASE("lim(x->inf) (1+1/x)^x exists");
    {
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto base = SymbolicExpr::add(one, inv_x);
        auto expr = SymbolicExpr::power(base, x);
        auto lim = LMCAS::limit_expression_checked(expr, "x", inf).value();
        EXPECT_TRUE(lim != nullptr, "limit is not null");
    }

    return TEST_REPORT();
}
