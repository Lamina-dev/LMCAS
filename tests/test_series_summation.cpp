/**
 * @file test_series_summation.cpp
 * @brief 验证有限符号求和与连乘的契约。
 */

#include "test_common.hpp"
#include "series_engine.hpp"

using namespace LMCAS;

using Expr = std::shared_ptr<SymbolicExpr>;
using Coeffs = std::vector<Expr>;

static Expr num(int n) { return SymbolicExpr::number(n); }
static Expr var(const std::string &name) { return SymbolicExpr::variable(name); }

TEST(SeriesSummation, SymbolicSumConstant) {
    auto result = LMCAS::symbolic_sum(num(3), "k", num(1), num(5));
    EXPECT_TRUE((result != nullptr)) << "result is not null";
    if (result) {
        EXPECT_TRUE(test_expression_text(result, num(15))) << "sum of constant 3 from 1 to 5 = 15";
    }
}

TEST(SeriesSummation, SymbolicSumLinear) {
    auto n = var("n");
    auto k = var("k");
    auto result = LMCAS::symbolic_sum(k, "k", num(1), n);
    ASSERT_TRUE(result) << "symbolic sum returns an expression";
    auto substituted = result->substitute("n", num(10));
    ASSERT_TRUE(substituted) << "the summation result accepts its upper bound";
    auto value = substituted->simplify();
    ASSERT_TRUE(value);
    EXPECT_EQ(value->compare(num(55)), 0)
        << "substituting n=10 into the closed form gives 55";
}

TEST(SeriesSummation, SymbolicSumQuadratic) {
    auto n = var("n");
    auto k = var("k");
    auto k_sq = SymbolicExpr::power(k, num(2));
    auto result = LMCAS::symbolic_sum(k_sq, "k", num(1), n);
    EXPECT_TRUE((result != nullptr)) << "result is not null";
    if (result) {
        auto val = result->substitute("n", num(5));
        if (val) {
            val = val->simplify();
            EXPECT_TRUE(test_expression_text(val, num(55))) << "sum k=1..5 k^2 = 55";
        }
    }
}

TEST(SeriesSummation, SymbolicSumCubic) {
    auto n = var("n");
    auto k = var("k");
    auto k_cubed = SymbolicExpr::power(k, num(3));
    auto result = LMCAS::symbolic_sum(k_cubed, "k", num(1), n);
    EXPECT_TRUE((result != nullptr)) << "result is not null";
    if (result) {
        auto val = result->substitute("n", num(4));
        if (val) {
            val = val->simplify();
            EXPECT_TRUE(test_expression_text(val, num(100))) << "sum k=1..4 k^3 = 100";
        }
    }
}

TEST(SeriesSummation, SymbolicSumGeometric) {
    auto n = var("n");
    auto k = var("k");
    auto two_pow_k = SymbolicExpr::power(num(2), k);
    auto result = LMCAS::symbolic_sum(two_pow_k, "k", num(0), n);
    EXPECT_TRUE((result != nullptr)) << "result is not null";
    if (result) {
        auto val = result->substitute("n", num(3));
        if (val) {
            val = val->simplify();
            EXPECT_TRUE(test_expression_text(val, num(15))) << "sum k=0..3 2^k = 15";
        }
    }
}

TEST(SeriesSummation, SymbolicSumDirectEval) {
    auto k = var("k");
    auto k_sq = SymbolicExpr::power(k, num(2));
    auto result = LMCAS::symbolic_sum(k_sq, "k", num(1), num(4));
    EXPECT_TRUE((result != nullptr)) << "result is not null";
    if (result) {
        result = result->simplify();
        EXPECT_TRUE(test_expression_text(result, num(30))) << "sum k=1..4 k^2 = 30";
    }
}

TEST(SeriesSummation, SymbolicSumEmptyRange) {
    auto k = var("k");
    auto result = LMCAS::symbolic_sum(k, "k", num(5), num(3));
    EXPECT_TRUE((result != nullptr)) << "result is not null";
    if (result) {
        EXPECT_TRUE(test_expression_text(result, num(0))) << "sum with upper < lower = 0";
    }
}

TEST(SeriesSummation, SymbolicSumUnevaluated) {
    auto k = var("k");
    auto n = var("n");
    auto body = SymbolicExpr::sin(k);
    auto result = LMCAS::symbolic_sum(body, "k", num(1), n);
    EXPECT_TRUE((result != nullptr)) << "result is not null";
    if (result) {
        auto sn = std::dynamic_pointer_cast<const SummationNode>(LMCAS::detail::node(result));
        EXPECT_TRUE((sn != nullptr)) << "result is SummationNode for sin(k)";
    }
}

TEST(SeriesSummation, SymbolicSumLargeExactBounds) {
    const BigInt lower_value("999999999999999999999999999999");
    const BigInt upper_value("1000000000000000000000000000000");
    auto k = var("k");
    auto body = SymbolicExpr::power(k, k);
    auto result = LMCAS::symbolic_sum(
        body, "k", SymbolicExpr::number(lower_value),
        SymbolicExpr::number(upper_value));

    EXPECT_TRUE((result != nullptr)) << "large-bound summation remains representable";
    if (result) {
        result = result->simplify();
        const std::string text = result ? result->to_string() : "";
        EXPECT_TRUE((text.find(lower_value.to_string()) != std::string::npos)) << "large lower bound remains exact after normalization";
        EXPECT_TRUE((text.find(upper_value.to_string()) != std::string::npos)) << "large upper bound remains exact after normalization";
    }
}

TEST(SeriesSummation, SymbolicProductFactorial) {
    auto k = var("k");
    auto result = LMCAS::symbolic_product(k, "k", num(1), num(5));
    EXPECT_TRUE((result != nullptr)) << "result is not null";
    if (result) {
        result = result->simplify();
        EXPECT_TRUE(test_expression_text(result, num(120))) << "prod k=1..5 k = 120";
    }
}

TEST(SeriesSummation, SymbolicProductDirectEval) {
    auto k = var("k");
    auto k_plus_1 = SymbolicExpr::add(k, num(1));
    auto result = LMCAS::symbolic_product(k_plus_1, "k", num(1), num(4));
    EXPECT_TRUE((result != nullptr)) << "result is not null";
    if (result) {
        result = result->simplify();
        EXPECT_TRUE(test_expression_text(result, num(120))) << "prod k=1..4 (k+1) = 120";
    }
}

TEST(SeriesSummation, SymbolicProductEmptyRange) {
    auto k = var("k");
    auto result = LMCAS::symbolic_product(k, "k", num(5), num(3));
    EXPECT_TRUE((result != nullptr)) << "result is not null";
    if (result) {
        EXPECT_TRUE(test_expression_text(result, num(1))) << "product with upper < lower = 1";
    }
}

TEST(SeriesSummation, SymbolicProductPochhammer) {
    auto k = var("k");
    auto body = SymbolicExpr::add(k, num(2));
    auto result = LMCAS::symbolic_product(body, "k", num(1), num(4));
    EXPECT_TRUE((result != nullptr)) << "result is not null";
    if (result) {
        result = result->simplify();
        EXPECT_TRUE(test_expression_text(result, num(360))) << "prod k=1..4 (k+2) = 360";
    }
}

TEST(SeriesSummation, SymbolicProductUnevaluated) {
    auto k = var("k");
    auto n = var("n");
    auto body = SymbolicExpr::sin(k);
    auto result = LMCAS::symbolic_product(body, "k", num(1), n);
    EXPECT_TRUE((result != nullptr)) << "result is not null";
    if (result) {
        auto pn = std::dynamic_pointer_cast<const ProductNode>(LMCAS::detail::node(result));
        EXPECT_TRUE((pn != nullptr)) << "result is ProductNode for sin(k)";
    }
}
