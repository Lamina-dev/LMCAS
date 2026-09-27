/**
 * @file test_calculus_utils.cpp
 * @brief 测试 calculus_utils 模块：对数微分、微分、全微分、反函数导数、反函数。
 */

#include "test_common.hpp"
#include "calculus_utils.hpp"
#include "assumption_context.hpp"
#include <optional>
#include <string>
#include <variant>

using namespace LMCAS;

using SE = SymbolicExpr;

static auto num(int n) { return SE::number(n); }
static auto var(const std::string &name) { return SE::variable(name); }
static LMCAS::ContinuityType checked_continuity(
    const std::shared_ptr<SymbolicExpr> &expression,
    const std::string &variable,
    const std::shared_ptr<SymbolicExpr> &point) {
    auto result =
        LMCAS::continuity_at_checked(expression, variable, point);
    EXPECT_TRUE((result.has_value())) << "checked continuity succeeds";
    return result ? result.value() : LMCAS::ContinuityType::Essential;
}

TEST(CalculusUtils, ContinuousPolynomial) {
    auto x = var("x");
    auto f = SE::power(x, num(2));
    auto result = checked_continuity(f, "x", num(1));
    EXPECT_TRUE((result == LMCAS::ContinuityType::Continuous)) << "x^2 is continuous at x=1";
}

TEST(CalculusUtils, EssentialDiscontinuity) {
    auto x = var("x");
    auto f = SE::divide(num(1), x);
    auto result = checked_continuity(f, "x", num(0));
    EXPECT_TRUE((result == LMCAS::ContinuityType::Essential)) << "1/x has essential discontinuity at x=0";
}

TEST(CalculusUtils, RemovableDiscontinuity) {
    auto x = var("x");
    auto numerator = SE::add(SE::power(x, num(2)), num(-1));
    auto denominator = SE::add(x, num(-1));
    auto f = SE::divide(numerator, denominator);
    auto result = checked_continuity(f, "x", num(1));
    EXPECT_TRUE((result == LMCAS::ContinuityType::Removable)) << "(x^2-1)/(x-1) has removable discontinuity at x=1";
}

TEST(CalculusUtils, JumpDiscontinuity) {
    auto x = var("x");
    auto one_over_x = SE::divide(num(1), x);
    auto exp_term = SE::exp(one_over_x);
    auto denom = SE::add(num(1), exp_term);
    auto f = SE::divide(num(1), denom);
    auto result = checked_continuity(f, "x", num(0));
    EXPECT_TRUE((result == LMCAS::ContinuityType::Jump)) << "1/(1+e^(1/x)) has jump discontinuity at x=0";
}

TEST(CalculusUtils, ParameterDependentContinuity) {
    auto a = var("a");
    auto zero = num(0);
    for (const auto relation : {RelationOp::LT, RelationOp::NEQ}) {
        auto condition = detail::make_node<RelationalNode>(
            detail::node(var("x")), detail::node(zero), relation);
        auto expression = detail::make_expression_ptr(detail::make_node<PiecewiseNode>(
            std::vector<PiecewiseNode::Branch>{{detail::node(a), condition}}, detail::node(zero)));
        for (const auto sign : {std::optional<Sign>{}, std::optional<Sign>{Sign::Zero},
                                std::optional<Sign>{Sign::NonZero}}) {
            auto assumptions = std::make_shared<AssumptionContext>();
            EXPECT_TRUE((assumptions->assume_domain_checked("a", Domain::Real).has_value())) << "a is real";
            if (sign) {
                EXPECT_TRUE((assumptions->assume_sign_checked("a", *sign).has_value())) << "parameter sign is attached";
            }
            ComputationContext context;
            EXPECT_TRUE((context.set_assumptions(assumptions).has_value())) << "continuity uses parameter assumptions";
            auto result = continuity_at_checked(expression, "x", zero, context);
            if (!sign) {
                EXPECT_TRUE((!result && result.error().code == CasErrc::Inconclusive)) << "a merely real parameter certifies neither a jump nor a removable point";
            } else {
                const auto expected = *sign == Sign::Zero ? ContinuityType::Continuous : relation == RelationOp::LT ? ContinuityType::Jump
                                                                                                                    : ContinuityType::Removable;
                EXPECT_TRUE((result && result.value() == expected)) << "proved zero/nonzero parameters retain the correct continuity classification";
            }
        }
    }
}

TEST(CalculusUtils, ContinuityNull) {
    auto result = LMCAS::continuity_at_checked(nullptr, "x", num(0));
    EXPECT_TRUE((!result.has_value() &&
                 result.error().code == LMCAS::CasErrc::InvalidArgument))
        << "null continuity input is rejected";
}

TEST(CalculusUtils, ReciprocalAsymptotes) {
    auto x = var("x");
    auto f = SE::divide(num(1), x);
    auto result = LMCAS::asymptotes_checked(f, "x").value();

    EXPECT_TRUE((!result.vertical.empty())) << "1/x has vertical asymptote(s)";
    if (!result.vertical.empty()) {
        auto val = test_numeric_eval(result.vertical[0]->simplify());
        if (val) {
            EXPECT_NEAR(*val, 0.0, 1e-9) << "vertical asymptote at x=0";
        } else {
            ADD_FAILURE() << "vertical asymptote numeric eval failed";
        }
    }

    EXPECT_TRUE((!result.horizontal.empty())) << "1/x has horizontal asymptote(s)";
    if (!result.horizontal.empty()) {
        auto h = result.horizontal[0]->simplify();
        bool is_zero_val = h->is_zero();
        if (!is_zero_val) {
            auto val = test_numeric_eval(h);
            is_zero_val = val && std::abs(*val) < 1e-9;
        }
        EXPECT_TRUE((is_zero_val)) << "horizontal asymptote at y=0";
    }

    EXPECT_TRUE((result.oblique.empty())) << "1/x has no oblique asymptotes";
}

TEST(CalculusUtils, RationalAsymptotes) {
    auto x = var("x");
    auto numerator = SE::add(SE::multiply(num(2), x), num(1));
    auto denominator = SE::add(x, num(-1));
    auto f = SE::divide(numerator, denominator);
    auto result = LMCAS::asymptotes_checked(f, "x").value();

    EXPECT_TRUE((!result.vertical.empty())) << "(2x+1)/(x-1) has vertical asymptote(s)";
    if (!result.vertical.empty()) {
        auto val = test_numeric_eval(result.vertical[0]->simplify());
        if (val) {
            EXPECT_NEAR(*val, 1.0, 1e-9) << "vertical asymptote at x=1";
        } else {
            ADD_FAILURE() << "vertical asymptote numeric eval failed";
        }
    }

    EXPECT_TRUE((!result.horizontal.empty())) << "(2x+1)/(x-1) has horizontal asymptote(s)";
    if (!result.horizontal.empty()) {
        auto val = test_numeric_eval(result.horizontal[0]->simplify());
        if (val) {
            EXPECT_NEAR(*val, 2.0, 1e-9) << "horizontal asymptote at y=2";
        } else {
            ADD_FAILURE() << "horizontal asymptote numeric eval failed";
        }
    }
}

TEST(CalculusUtils, ObliqueAsymptotes) {
    /**
     * @brief f(x) = (x^2+1)/(x-1) 的垂直渐近线为 x=1，斜渐近线为 y=x+1。
     * 斜率 lim(x→∞) f(x)/x = lim (x^2+1)/(x(x-1)) = 1；
     * 截距 lim(x→∞) [f(x)-x] = lim (x+1)/(x-1) = 1。
     */
    auto x = var("x");
    auto numerator = SE::add(SE::power(x, num(2)), num(1));
    auto denominator = SE::add(x, num(-1));
    auto f = SE::divide(numerator, denominator);
    auto result = LMCAS::asymptotes_checked(f, "x").value();

    ASSERT_FALSE(result.vertical.empty()) << "(x^2+1)/(x-1) has vertical asymptote(s)";
    const auto vertical = test_numeric_eval(result.vertical[0]->simplify());
    ASSERT_TRUE(vertical.has_value()) << "vertical asymptote numeric eval failed";
    EXPECT_TRUE(std::isfinite(*vertical));
    EXPECT_NEAR(*vertical, 1.0, 1e-9) << "vertical asymptote at x=1";

    ASSERT_FALSE(result.oblique.empty()) << "(x^2+1)/(x-1) has oblique asymptote(s)";
    const auto slope = test_numeric_eval(result.oblique[0].first->simplify());
    const auto intercept = test_numeric_eval(result.oblique[0].second->simplify());
    ASSERT_TRUE(slope.has_value()) << result.oblique[0].first->to_string();
    ASSERT_TRUE(intercept.has_value()) << result.oblique[0].second->to_string();
    EXPECT_TRUE(std::isfinite(*slope));
    EXPECT_NEAR(*slope, 1.0, 1e-9) << "oblique asymptote slope = 1";
    EXPECT_TRUE(std::isfinite(*intercept));
    EXPECT_NEAR(*intercept, 1.0, 1e-9) << "oblique asymptote intercept = 1";

    EXPECT_TRUE((result.horizontal.empty())) << "(x^2+1)/(x-1) has no horizontal asymptotes";
}

TEST(CalculusUtils, AsymptoteInvalidInputs) {
    auto null_result = LMCAS::asymptotes_checked(nullptr, "x");
    EXPECT_TRUE((!null_result &&
                 null_result.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked asymptotes rejects null expression";

    auto empty_var = LMCAS::asymptotes_checked(var("x"), "");
    EXPECT_TRUE((!empty_var &&
                 empty_var.error().code == LMCAS::CasErrc::InvalidArgument))
        << "checked asymptotes rejects empty variable";
}

TEST(CalculusUtils, AsymptoteBudget) {
    auto x = var("x");
    auto denominator = SE::add(x, num(-1));
    auto f = SE::divide(num(1), denominator);

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = LMCAS::asymptotes_checked(f, "x", limited_context);
    EXPECT_TRUE((!limited &&
                 limited.error().code == LMCAS::CasErrc::ResourceLimit))
        << "checked asymptotes propagates denominator solve budget failure";
}
