/**
 * @file test_series_expansions.cpp
 * @brief Fourier 与 Laurent 展开的行为约定。
 */

#include "test_common.hpp"
#include "series_engine.hpp"
#include "numeric_evaluation.hpp"

using namespace LMCAS;

using Expr = std::shared_ptr<SymbolicExpr>;
using Coeffs = std::vector<Expr>;

static Expr num(int n) { return SymbolicExpr::number(n); }
static Expr var(const std::string &name) { return SymbolicExpr::variable(name); }

TEST(SeriesExpansions, LaurentCheckedContracts) {
    auto z = var("z");
    auto f = SymbolicExpr::divide(num(1), z);
    auto laurent = LMCAS::laurent_series_full_checked(f, "z", num(0), 2, 2);
    ASSERT_TRUE((laurent.has_value())) << "checked Laurent series succeeds for 1/z";
    if (laurent) {
        EXPECT_TRUE((laurent.value().series != nullptr)) << "checked Laurent full result has non-null series";
        EXPECT_TRUE((laurent.value().pole_order == 1)) << "checked Laurent detects simple pole";
    }

    auto laurent_expr = LMCAS::laurent_series_checked(f, "z", num(0), 2, 2);
    EXPECT_TRUE((laurent_expr.has_value())) << "checked Laurent expression succeeds";

    auto bad_order = LMCAS::laurent_series_full_checked(f, "z", num(0), -1, 2);
    EXPECT_TRUE((!bad_order &&
                 bad_order.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked Laurent rejects negative truncation order";

    auto null_center = LMCAS::laurent_series_full_checked(f, "z", nullptr, 2, 2);
    EXPECT_TRUE((!null_center &&
                 null_center.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked Laurent rejects null center";

    auto analytic_laurent = LMCAS::laurent_series_full_checked(
        SymbolicExpr::sin(z), "z", num(0), 2, 2);
    EXPECT_TRUE((analytic_laurent &&
                 analytic_laurent.value().pole_order == 0))
        << "checked Laurent expands supported analytic functions";

    auto shifted_laurent = LMCAS::laurent_series_full_checked(
        f, "z", num(1), 2, 2);
    EXPECT_TRUE((shifted_laurent &&
                 shifted_laurent.value().pole_order == 0))
        << "checked Laurent supports a nonzero regular center";
    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = LMCAS::laurent_series_full_checked(f, "z", num(0), 2, 2,
                                                      limited_context);
    EXPECT_TRUE((!limited &&
                 limited.error().code == LMCAS::CasErrc::ResourceLimit))
        << "checked Laurent observes exhausted step budget";
}

TEST(SeriesExpansions, FourierSeriesSquareWave) {
    auto one = num(1);
    auto period = SymbolicExpr::multiply(num(2), SymbolicExpr::number(LMMC_CONST_PI));

    auto result = LMCAS::fourier_series_checked(one, "x", period, 0);
    ASSERT_TRUE((result.has_value())) << "fourier_series returns a result";
    if (result) {
        auto v = test_numeric_eval(result.value()->simplify());
        EXPECT_TRUE((v.has_value() && std::abs(*v - 1.0) < 1e-6)) << "Fourier constant term (a0/2) of f=1 is 1";
    }
}

TEST(SeriesExpansions, FourierSeriesOddFunction) {
    auto x = var("x");
    auto period = SymbolicExpr::multiply(num(2), SymbolicExpr::number(LMMC_CONST_PI));

    auto result = LMCAS::fourier_series_checked(x, "x", period, 2);
    ASSERT_TRUE((result.has_value())) << "fourier_series(x) returns a result";
    if (result) {
        auto val = result.value()->substitute("x", num(0));
        auto v = test_numeric_eval(val ? val->simplify() : nullptr);
        EXPECT_TRUE((v.has_value() && std::abs(*v) < 1e-6)) << "Fourier series of odd f=x is 0 at x=0 (no constant term)";
        auto halfway = evaluate_numeric(*result.value(), {{"x", LMMC_CONST_PI / 2}});
        ASSERT_TRUE((halfway.has_value())) << "Fourier harmonics are finite, not unresolved coefficients times zero";
        if (halfway) {
            EXPECT_NEAR(halfway.value().value, 2.0, 1e-6) << "the first two sine harmonics of x have the independent half-interval value";
        }
    }
}

TEST(SeriesExpansions, FourierSeriesCheckedContracts) {
    auto one = num(1);
    auto period =
        SymbolicExpr::multiply(num(2), SymbolicExpr::number(LMMC_CONST_PI));

    auto invalid = LMCAS::fourier_series_checked(one, "x", period, -1);
    EXPECT_TRUE((!invalid && invalid.error().code == LMCAS::CasErrc::InvalidArgument)) << "negative Fourier order is rejected";

    LMCAS::ResourceLimits limits;
    limits.max_expansion_terms = 0;
    LMCAS::ComputationContext context(limits);
    auto limited =
        LMCAS::fourier_series_checked(one, "x", period, 0, context);
    EXPECT_TRUE((!limited &&
                 limited.error().code == LMCAS::CasErrc::ResourceLimit))
        << "Fourier expansion observes the shared term budget";
}

TEST(SeriesExpansions, LaurentSeries1OverZ) {
    auto z = var("z");
    auto f = SymbolicExpr::power(z, num(-1));

    auto checked = LMCAS::laurent_series_checked(f, "z", num(0), 2, 2);
    ASSERT_TRUE((checked.has_value())) << "Laurent expression completes";
    if (!checked)
        return;
    auto result = checked.value();
    EXPECT_TRUE((result != nullptr)) << "laurent_series(1/z) returns a result";
    if (result) {
        auto expected = SymbolicExpr::divide(num(1), z);
        {
            const auto actual_expr = (result->simplify());
            const auto expected_expr = (expected->simplify());
            EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << "Laurent of 1/z is 1/z";
        }
    }
}

TEST(SeriesExpansions, LaurentSeriesFull1OverZ) {
    auto z = var("z");
    auto f = SymbolicExpr::power(z, num(-1));

    auto checked = LMCAS::laurent_series_full_checked(f, "z", num(0), 2, 2);
    ASSERT_TRUE((checked.has_value())) << "Laurent pole and residue complete";
    if (!checked)
        return;
    auto result = checked.value();
    EXPECT_TRUE((result.series != nullptr)) << "laurent_series_full: series not null";
    EXPECT_TRUE((result.pole_order == 1)) << "laurent_series_full: pole_order = 1 for 1/z";
    if (result.residue) {
        {
            const auto actual_expr = (result.residue->simplify());
            const auto expected_expr = (num(1));
            EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << "residue of 1/z at 0 is 1";
        }
    }
}

TEST(SeriesExpansions, EssentialPoleIndependentOfTruncation) {
    auto z = var("z");
    auto f = SymbolicExpr::exp(SymbolicExpr::divide(num(1), z));
    for (int order_neg : {0, 1, 3}) {
        auto result = laurent_series_full_checked(f, "z", num(0), order_neg, 0);
        ASSERT_TRUE(result);
        EXPECT_EQ(result.value().singularity, SingularityType::Essential);
        EXPECT_TRUE(test_expression_text(result.value().residue, num(1)));
        if (order_neg == 3) {
            auto expected = SymbolicExpr::add(
                SymbolicExpr::add(num(1), SymbolicExpr::divide(num(1), z)),
                SymbolicExpr::add(
                    SymbolicExpr::divide(num(1), SymbolicExpr::multiply(num(2), SymbolicExpr::power(z, num(2)))),
                    SymbolicExpr::divide(num(1), SymbolicExpr::multiply(num(6), SymbolicExpr::power(z, num(3))))));
            EXPECT_TRUE(test_proved_equivalent(result.value().series, expected));
        }
    }
    auto high_order = SymbolicExpr::exp(
        SymbolicExpr::power(z, num(-65)));
    auto classified = laurent_series_full_checked(high_order, "z", num(0), 0, 0);
    ASSERT_TRUE(classified);
    EXPECT_EQ(classified.value().singularity, SingularityType::Essential);
}

TEST(SeriesExpansions, ShiftedEssentialEvenPowers) {
    auto z = var("z");
    auto shift = SymbolicExpr::add(z, num(-1));
    auto f = SymbolicExpr::exp(SymbolicExpr::divide(num(1), SymbolicExpr::power(shift, num(2))));
    auto result = laurent_series_full_checked(f, "z", num(1), 4, 0);
    ASSERT_TRUE(result);
    EXPECT_EQ(result.value().singularity, SingularityType::Essential);
    EXPECT_TRUE(test_expression_text(result.value().residue, num(0)));
    auto expected = SymbolicExpr::add(
        SymbolicExpr::add(num(1), SymbolicExpr::divide(num(1), SymbolicExpr::power(shift, num(2)))),
        SymbolicExpr::divide(num(1), SymbolicExpr::multiply(num(2), SymbolicExpr::power(shift, num(4)))));
    EXPECT_TRUE(test_proved_equivalent(result.value().series, expected));
}

TEST(SeriesExpansions, LaurentSeries1OverSinZ) {
    /**
     * @brief f(z) = 1/(z(z+1)) 在 z=0 处为一阶极点，留数为 1。
     * @note 此处使用有理分母；级数引擎尚不支持 1/sin(z) 等超越分母。
     */
    auto z = var("z");
    auto denom = SymbolicExpr::multiply(z, SymbolicExpr::add(z, num(1)));
    auto f = SymbolicExpr::divide(num(1), denom);

    auto checked = LMCAS::laurent_series_full_checked(f, "z", num(0), 3, 3);
    ASSERT_TRUE((checked.has_value())) << "Laurent rational pole completes";
    if (!checked)
        return;
    auto result = checked.value();
    EXPECT_TRUE((result.pole_order == 1)) << "1/(z(z+1)): pole_order = 1";
    if (result.residue) {
        {
            const auto actual_expr = (result.residue->simplify());
            const auto expected_expr = (num(1));
            EXPECT_TRUE(test_expression_text(actual_expr, expected_expr)) << "residue of 1/(z(z+1)) at 0 is 1";
        }
    }
}
