/**
 * @file test_debug_limit.cpp
 * @brief Debug test for (1+1/x)^x limit
 */
#include "test_common.hpp"
#include "visitors/limit_visitor.hpp"
#include "visitors/differentiation_visitor.hpp"

int main() {
    auto x = SymbolicExpr::variable("x");
    auto one = SymbolicExpr::number(1);
    auto neg_one = SymbolicExpr::number(-1);
    auto inf = SymbolicExpr::infinity(1);

    // Test: lim(x->inf) x * ln(1+1/x) should be 1
    TEST_CASE("Debug: lim(x->inf) x*ln(1+1/x) = 1");
    {
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto base = SymbolicExpr::add(one, inv_x);
        auto ln_base = SymbolicExpr::ln(base);
        auto product = SymbolicExpr::multiply(x, ln_base);
        std::cout << "  expr = " << product->to_string() << std::endl;
        auto lim = product->limit("x", inf);
        if (lim) {
            std::cout << "  result = " << lim->to_string() << std::endl;
            auto val = test_numeric_eval(lim);
            if (val) {
                std::cout << "  numeric = " << *val << std::endl;
                EXPECT_NEAR(*val, 1.0, 1e-6, "lim x*ln(1+1/x) = 1");
            } else {
                EXPECT_EQ_EXPR_STR(lim, "1", "lim x*ln(1+1/x) = 1");
            }
        } else {
            EXPECT_TRUE(false, "limit is null");
        }
    }

    // Test: differentiate ln(1+1/x)
    TEST_CASE("Debug: d/dx[ln(1+1/x)]");
    {
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto base = SymbolicExpr::add(one, inv_x);
        auto ln_base = SymbolicExpr::ln(base);
        auto deriv = ln_base->differentiate("x");
        if (deriv) {
            std::cout << "  d/dx[ln(1+1/x)] = " << deriv->to_string() << std::endl;
        }
        EXPECT_TRUE(deriv != nullptr, "derivative exists");
    }

    // Test: differentiate x^(-1)
    TEST_CASE("Debug: d/dx[x^(-1)]");
    {
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto deriv = inv_x->differentiate("x");
        if (deriv) {
            std::cout << "  d/dx[x^(-1)] = " << deriv->to_string() << std::endl;
        }
        EXPECT_TRUE(deriv != nullptr, "derivative exists");
    }

    // Test: the full (1+1/x)^x limit
    TEST_CASE("Debug: lim(x->inf) (1+1/x)^x = e");
    {
        auto inv_x = SymbolicExpr::power(x, neg_one);
        auto base = SymbolicExpr::add(one, inv_x);
        auto expr = SymbolicExpr::power(base, x);
        std::cout << "  expr = " << expr->to_string() << std::endl;
        auto lim = expr->limit("x", inf);
        if (lim) {
            std::cout << "  result = " << lim->to_string() << std::endl;
            auto val = test_numeric_eval(lim);
            if (val) {
                std::cout << "  numeric = " << *val << std::endl;
            }
        } else {
            std::cout << "  result = null" << std::endl;
        }
        EXPECT_TRUE(lim != nullptr, "limit is not null");
    }

    return TEST_REPORT();
}
