#include "test_common.hpp"
#include "transcendental_factor.hpp"

using namespace LMCAS;

namespace {

std::shared_ptr<SymbolicExpr> repeated_mod_three_factor_fixture() {
    auto x = SymbolicExpr::variable("x");
    auto sine = SymbolicExpr::sin(x);
    return SymbolicExpr::add(
        SymbolicExpr::power(sine, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(5), sine),
            SymbolicExpr::number(4)));
}

} // namespace

TEST(TranscendentalFactorPipeline, MultStructureXTimesSinX) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::multiply(x, SymbolicExpr::sin(x));

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 2)) << "should return 2 factors";
    bool has_x = false, has_sin = false;
    for (const auto &f : factors) {
        std::string s = f->to_string();
        if (s == "x")
            has_x = true;
        if (s.find("sin") != std::string::npos && s.find("x") != std::string::npos)
            has_sin = true;
    }
    EXPECT_TRUE((has_x)) << "should contain factor x";
    EXPECT_TRUE((has_sin)) << "should contain factor sin(x)";
}

TEST(TranscendentalFactorPipeline, MultStructureX2TimesExpTimesCos) {
    auto x = SymbolicExpr::variable("x");
    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto exp_x = SymbolicExpr::exp(x);
    auto cos_x = SymbolicExpr::cos(x);

    auto expr = SymbolicExpr::multiply(x2, SymbolicExpr::multiply(exp_x, cos_x));

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 3)) << "should return 3 factors";
}

TEST(TranscendentalFactorPipeline, MultStructureWithConstant) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto expr = SymbolicExpr::multiply(two, SymbolicExpr::multiply(x, SymbolicExpr::sin(x)));

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 3)) << "should return 3 factors (constant + 2 non-constant)";
    bool has_const = false;
    for (const auto &f : factors) {
        if (f->is_number())
            has_const = true;
    }
    EXPECT_TRUE((has_const)) << "should contain numeric constant factor";
}

TEST(TranscendentalFactorPipeline, MultStructureSumNotProduct) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 1)) << "sum expression should not be split by multiplicative detection";
}

TEST(TranscendentalFactorPipeline, LinearIrreducibleSinPlusX) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 1)) << "sin(x) + x should be irreducible (linear in u0 and x)";
    EXPECT_TRUE(test_expression_text((factors[0]), (expr->to_string()))) << "factor should be original expression";
}

TEST(TranscendentalFactorPipeline, LinearIrreducible2sin3x1) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::sin(x)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(3), x),
            SymbolicExpr::number(1)));

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 1)) << "2*sin(x) + 3*x + 1 should be irreducible";
}

TEST(TranscendentalFactorPipeline, LinearIrreducibleSinCosX) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(
        SymbolicExpr::sin(x),
        SymbolicExpr::add(SymbolicExpr::cos(x), x));

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 1)) << "sin(x) + cos(x) + x should be irreducible (linear in u0, u1, x)";
}

TEST(TranscendentalFactorPipeline, LinearIrreducibleSinTimesXNotSum) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::multiply(SymbolicExpr::sin(x), x);

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 2)) << "sin(x) * x should be split into 2 factors by multiplicative detection";
}

TEST(TranscendentalFactorPipeline, ExpSeparationExpXTimesXPlusExpX) {
    auto x = SymbolicExpr::variable("x");
    auto exp_x = SymbolicExpr::exp(x);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(exp_x, x),
        exp_x);

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 2)) << "should produce 2 factors: exp(x) and (x+1)";

    bool has_exp = false;
    bool has_poly = false;
    for (const auto &f : factors) {
        std::string s = f->to_string();
        if (s.find("exp") != std::string::npos && s.find("+") == std::string::npos) {
            has_exp = true;
        }
        if (s.find("exp") == std::string::npos) {
            has_poly = true;
        }
    }
    EXPECT_TRUE((has_exp)) << "should contain exp(x) as a factor";
    EXPECT_TRUE((has_poly)) << "should contain polynomial remainder as a factor";
}

TEST(TranscendentalFactorPipeline, ExpSeparationNoCommonExp) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::exp(x), SymbolicExpr::sin(x));

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 1)) << "exp(x) + sin(x) should not be split by exponential separation";
}

TEST(TranscendentalFactorPipeline, ExpSeparationProductNotSum) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::multiply(SymbolicExpr::exp(x), SymbolicExpr::sin(x));

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 2)) << "exp(x) * sin(x) should be split by multiplicative detection";
}

TEST(TranscendentalFactorPipeline, ExpSeparationThreeTerms) {
    auto x = SymbolicExpr::variable("x");
    auto exp_x = SymbolicExpr::exp(x);
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(exp_x, x_sq),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::multiply(exp_x, x)),
            exp_x));

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 2)) << "should produce 2 factors: exp(x) and (x^2+2x+1)";

    bool has_exp = false;
    for (const auto &f : factors) {
        std::string s = f->to_string();
        if (s.find("exp") != std::string::npos && s.find("+") == std::string::npos && s.find("*") == std::string::npos) {
            has_exp = true;
        }
    }
    EXPECT_TRUE((has_exp)) << "should contain exp(x) as a factor";
}

TEST(TranscendentalFactorPipeline, PythagoreanSin2Cos2To1) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);
    auto sin2 = SymbolicExpr::power(sin_x, two);
    auto cos2 = SymbolicExpr::power(cos_x, two);
    auto expr = SymbolicExpr::add(sin2, cos2);

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 1)) << "should produce 1 factor";
    if (!factors.empty()) {
        std::string s = factors[0]->to_string();
        EXPECT_TRUE((factors[0]->is_number() || s == "1")) << "sin²(x) + cos²(x) should simplify to 1";
    }
}

TEST(TranscendentalFactorPipeline, PythagoreanSin2Cos2PlusX) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);
    auto sin2 = SymbolicExpr::power(sin_x, two);
    auto cos2 = SymbolicExpr::power(cos_x, two);
    auto expr = SymbolicExpr::add(sin2, SymbolicExpr::add(cos2, x));

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() >= 1)) << "should produce at least 1 factor";
    for (const auto &f : factors) {
        std::string s = f->to_string();
        EXPECT_TRUE((s.find("sin") == std::string::npos)) << "result should not contain sin after Pythagorean simplification";
        EXPECT_TRUE((s.find("cos") == std::string::npos)) << "result should not contain cos after Pythagorean simplification";
    }
}

TEST(TranscendentalFactorPipeline, PythagoreanCommonCoefficient) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto sin_x = SymbolicExpr::sin(x);
    auto cos_x = SymbolicExpr::cos(x);
    auto sin2 = SymbolicExpr::power(sin_x, two);
    auto cos2 = SymbolicExpr::power(cos_x, two);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(two, sin2),
        SymbolicExpr::multiply(two, cos2));

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 1)) << "should produce 1 factor";
    if (!factors.empty()) {
        EXPECT_TRUE((factors[0]->is_number())) << "2*sin²(x) + 2*cos²(x) should simplify to 2";
    }
}

TEST(TranscendentalFactorPipeline, PythagoreanDifferentArgsUnchanged) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto two = SymbolicExpr::number(2);
    auto sin2_x = SymbolicExpr::power(SymbolicExpr::sin(x), two);
    auto cos2_y = SymbolicExpr::power(SymbolicExpr::cos(y), two);
    auto expr = SymbolicExpr::add(sin2_x, cos2_y);

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    bool has_trig = false;
    for (const auto &f : factors) {
        std::string s = f->to_string();
        if (s.find("sin") != std::string::npos || s.find("cos") != std::string::npos) {
            has_trig = true;
        }
    }
    EXPECT_TRUE((has_trig)) << "sin²(x) + cos²(y) should remain unchanged (different args)";
}

TEST(TranscendentalFactorPipeline, PythagoreanNoMatchingCos) {
    auto x = SymbolicExpr::variable("x");
    auto two = SymbolicExpr::number(2);
    auto sin2_x = SymbolicExpr::power(SymbolicExpr::sin(x), two);
    auto expr = SymbolicExpr::add(sin2_x, x);

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    bool has_sin = false;
    for (const auto &f : factors) {
        std::string s = f->to_string();
        if (s.find("sin") != std::string::npos) {
            has_sin = true;
        }
    }
    EXPECT_TRUE((has_sin)) << "sin²(x) + x should remain unchanged (no matching cos²)";
}

TEST(TranscendentalFactorPipeline, FullPipelineSin2MinusX2) {
    auto x = SymbolicExpr::variable("x");
    auto sin_x = SymbolicExpr::sin(x);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(sin_x, SymbolicExpr::number(2)),
        SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::power(x, SymbolicExpr::number(2))));

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    /**
     * @brief 差平方满足 sin^2(x) - x^2 = (sin(x)+x)(sin(x)-x)。
     * 换元后为双变量多项式 u0^2 - x^2；构造成功时分为两个因子，
     * 否则保留原式。结果至少含一个因子，并保留 sin。
     */
    EXPECT_TRUE((factors.size() >= 1)) << "should produce at least 1 factor";

    if (factors.size() == 2) {
        std::string f0 = factors[0]->to_string();
        std::string f1 = factors[1]->to_string();
        bool both_have_sin = (f0.find("sin") != std::string::npos) &&
                             (f1.find("sin") != std::string::npos);
        EXPECT_TRUE((both_have_sin)) << "both factors should contain sin(x)";
    } else {
        bool has_sin = false;
        for (const auto &f : factors) {
            if (f->to_string().find("sin") != std::string::npos)
                has_sin = true;
        }
        EXPECT_TRUE((has_sin)) << "irreducible result should preserve sin(x)";
    }
}

TEST(TranscendentalFactorPipeline, FullPipelineExpTimesX) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::multiply(SymbolicExpr::exp(x), x);

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 2)) << "should produce 2 factors";

    if (factors.size() == 2) {
        bool has_exp = false;
        bool has_x = false;
        for (const auto &f : factors) {
            std::string s = f->to_string();
            if (s.find("exp") != std::string::npos)
                has_exp = true;
            if (s == "x")
                has_x = true;
        }
        EXPECT_TRUE((has_exp)) << "one factor should be exp(x)";
        EXPECT_TRUE((has_x)) << "one factor should be x";
    }
}

TEST(TranscendentalFactorPipeline, FullPipelineIrreducibleSinPlusX) {
    auto x = SymbolicExpr::variable("x");
    auto expr = SymbolicExpr::add(SymbolicExpr::sin(x), x);

    auto checked = factor_transcendental(expr, "x");
    ASSERT_TRUE(checked) << checked.error().message;

    auto factors = std::move(checked.value());

    EXPECT_TRUE((factors.size() == 1)) << "should return 1 factor (irreducible)";
}

TEST(TranscendentalFactorPipeline, UnsuitableFirstPrimeIsSkipped) {
    auto expression = repeated_mod_three_factor_fixture();
    ComputationContext context;

    auto checked = factor_transcendental(expression, "x", context);

    ASSERT_TRUE(checked) << checked.error().message;
    ASSERT_EQ(checked.value().size(), 2U);
    auto reconstructed = SymbolicExpr::number(1);
    for (const auto &factor : checked.value()) {
        reconstructed = SymbolicExpr::multiply(reconstructed, factor);
    }
    EXPECT_TRUE(test_proved_equivalent(
        reconstructed->expand(), expression->expand()));
}

TEST(TranscendentalFactorPipeline, ContextCancellationPropagates) {
    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext context({}, cancellation);

    auto checked = factor_transcendental(
        repeated_mod_three_factor_fixture(), "x", context);

    ASSERT_FALSE(checked);
    EXPECT_EQ(checked.error().code, CasErrc::Cancelled);
}

TEST(TranscendentalFactorPipeline, ContextResourceLimitPropagates) {
    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext context(limits);

    auto checked = factor_transcendental(
        repeated_mod_three_factor_fixture(), "x", context);

    ASSERT_FALSE(checked);
    EXPECT_EQ(checked.error().code, CasErrc::ResourceLimit);
}
