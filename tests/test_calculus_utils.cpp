/**
 * @file test_calculus_utils.cpp
 * @brief 测试 calculus_utils 模块：对数微分、微分、全微分、反函数导数、反函数。
 */

#include "test_common.hpp"
#include "calculus_utils.hpp"
#include <string>
#include <variant>

using namespace LMCAS;

using SE = SymbolicExpr;

static auto num(int n) { return SE::number(n); }
static auto var(const std::string& name) { return SE::variable(name); }
static std::shared_ptr<SymbolicExpr> bigint_num(const BigInt& n) {
    return LMCAS::detail::make_expression_ptr(
        LMCAS::detail::make_node<NumberNode>(
            std::variant<BigInt, Rational, lmmc_real_t>{
                std::in_place_type<BigInt>, n}));
}
static LMCAS::ContinuityType checked_continuity(
    const std::shared_ptr<SymbolicExpr>& expression,
    const std::string& variable,
    const std::shared_ptr<SymbolicExpr>& point) {
    auto result =
        LMCAS::continuity_at_checked(expression, variable, point);
    EXPECT_TRUE(result.has_value(), "checked continuity succeeds");
    return result ? result.value() : LMCAS::ContinuityType::Essential;
}

int main() {
    TEST_CASE("continuity_at: polynomial is continuous everywhere");
    {
        // f(x) = x^2 is continuous at x=1
        auto x = var("x");
        auto f = SE::power(x, num(2));
        auto result = checked_continuity(f, "x", num(1));
        EXPECT_TRUE(result == LMCAS::ContinuityType::Continuous,
                    "x^2 is continuous at x=1");
    }

    TEST_CASE("continuity_at: 1/x has essential discontinuity at x=0");
    {
        // f(x) = 1/x has essential discontinuity at x=0 (limit is ±∞)
        auto x = var("x");
        auto f = SE::divide(num(1), x);
        auto result = checked_continuity(f, "x", num(0));
        EXPECT_TRUE(result == LMCAS::ContinuityType::Essential,
                    "1/x has essential discontinuity at x=0");
    }

    TEST_CASE("continuity_at: removable discontinuity (x^2-1)/(x-1) at x=1");
    {
        // f(x) = (x^2 - 1)/(x - 1) = x + 1 for x ≠ 1, undefined at x=1
        // Left limit = right limit = 2, but f(1) is undefined → removable
        auto x = var("x");
        auto numerator = SE::add(SE::power(x, num(2)), num(-1)); // x^2 - 1
        auto denominator = SE::add(x, num(-1));                   // x - 1
        auto f = SE::divide(numerator, denominator);
        auto result = checked_continuity(f, "x", num(1));
        EXPECT_TRUE(result == LMCAS::ContinuityType::Removable,
                    "(x^2-1)/(x-1) has removable discontinuity at x=1");
    }

    TEST_CASE("continuity_at: jump discontinuity 1/(1+e^(1/x)) at x=0");
    {
        // f(x) = 1/(1 + e^(1/x))
        // As x→0+: e^(1/x)→∞, f→0
        // As x→0-: e^(1/x)→0, f→1
        // Left limit ≠ right limit → jump discontinuity
        auto x = var("x");
        auto one_over_x = SE::divide(num(1), x);
        auto exp_term = SE::exp(one_over_x);
        auto denom = SE::add(num(1), exp_term);
        auto f = SE::divide(num(1), denom);
        auto result = checked_continuity(f, "x", num(0));
        EXPECT_TRUE(result == LMCAS::ContinuityType::Jump,
                    "1/(1+e^(1/x)) has jump discontinuity at x=0");
    }

    TEST_CASE("continuity_at_checked: null input is invalid");
    {
        auto result = LMCAS::continuity_at_checked(nullptr, "x", num(0));
        EXPECT_TRUE(!result.has_value() &&
                        result.error().code == LMCAS::CasErrc::InvalidArgument,
                    "null continuity input is rejected");
    }

    TEST_CASE("asymptotes: 1/x has vertical at x=0, horizontal at y=0");
    {
        // f(x) = 1/x: vertical asymptote at x=0, horizontal at y=0
        auto x = var("x");
        auto f = SE::divide(num(1), x);
        auto result = LMCAS::asymptotes_checked(f, "x").value();

        // Vertical asymptote at x=0
        EXPECT_TRUE(!result.vertical.empty(),
                    "1/x has vertical asymptote(s)");
        if (!result.vertical.empty()) {
            auto val = test_numeric_eval(result.vertical[0]->simplify());
            if (val) {
                EXPECT_NEAR(*val, 0.0, 1e-9, "vertical asymptote at x=0");
            } else {
                EXPECT_TRUE(false, "vertical asymptote numeric eval failed");
            }
        }

        // Horizontal asymptote at y=0
        EXPECT_TRUE(!result.horizontal.empty(),
                    "1/x has horizontal asymptote(s)");
        if (!result.horizontal.empty()) {
            auto h = result.horizontal[0]->simplify();
            std::cout << "  horizontal asymptote expr: " << h->to_string() << std::endl;
            // The limit of 1/x as x→∞ is 0; check via is_zero() or numeric eval
            bool is_zero_val = h->is_zero();
            if (!is_zero_val) {
                auto val = test_numeric_eval(h);
                is_zero_val = val && std::abs(*val) < 1e-9;
            }
            EXPECT_TRUE(is_zero_val, "horizontal asymptote at y=0");
        }

        // No oblique asymptotes
        EXPECT_TRUE(result.oblique.empty(),
                    "1/x has no oblique asymptotes");
    }

    TEST_CASE("asymptotes: (2x+1)/(x-1) has vertical at x=1, horizontal at y=2");
    {
        // f(x) = (2x+1)/(x-1): vertical at x=1, horizontal at y=2
        auto x = var("x");
        auto numerator = SE::add(SE::multiply(num(2), x), num(1));
        auto denominator = SE::add(x, num(-1));
        auto f = SE::divide(numerator, denominator);
        auto result = LMCAS::asymptotes_checked(f, "x").value();

        // Vertical asymptote at x=1
        EXPECT_TRUE(!result.vertical.empty(),
                    "(2x+1)/(x-1) has vertical asymptote(s)");
        if (!result.vertical.empty()) {
            auto val = test_numeric_eval(result.vertical[0]->simplify());
            if (val) {
                EXPECT_NEAR(*val, 1.0, 1e-9, "vertical asymptote at x=1");
            } else {
                EXPECT_TRUE(false, "vertical asymptote numeric eval failed");
            }
        }

        // Horizontal asymptote at y=2
        EXPECT_TRUE(!result.horizontal.empty(),
                    "(2x+1)/(x-1) has horizontal asymptote(s)");
        if (!result.horizontal.empty()) {
            auto val = test_numeric_eval(result.horizontal[0]->simplify());
            if (val) {
                EXPECT_NEAR(*val, 2.0, 1e-9, "horizontal asymptote at y=2");
            } else {
                EXPECT_TRUE(false, "horizontal asymptote numeric eval failed");
            }
        }
    }

    TEST_CASE("asymptotes: (x^2+1)/(x-1) has oblique asymptote y=x+1");
    {
        // f(x) = (x^2+1)/(x-1): vertical at x=1, oblique y = x + 1
        // slope = lim(x→∞) f(x)/x = lim (x^2+1)/(x(x-1)) = 1
        // intercept = lim(x→∞) [f(x) - x] = lim [(x^2+1)/(x-1) - x]
        //           = lim [(x^2+1 - x(x-1))/(x-1)] = lim [(x+1)/(x-1)] = 1
        auto x = var("x");
        auto numerator = SE::add(SE::power(x, num(2)), num(1)); // x^2 + 1
        auto denominator = SE::add(x, num(-1));                   // x - 1
        auto f = SE::divide(numerator, denominator);
        auto result = LMCAS::asymptotes_checked(f, "x").value();

        // Vertical asymptote at x=1
        EXPECT_TRUE(!result.vertical.empty(),
                    "(x^2+1)/(x-1) has vertical asymptote(s)");
        if (!result.vertical.empty()) {
            auto val = test_numeric_eval(result.vertical[0]->simplify());
            if (val) {
                EXPECT_NEAR(*val, 1.0, 1e-9, "vertical asymptote at x=1");
            } else {
                EXPECT_TRUE(false, "vertical asymptote numeric eval failed");
            }
        }

        // Should have oblique asymptote with slope=1, intercept=1
        EXPECT_TRUE(!result.oblique.empty(),
                    "(x^2+1)/(x-1) has oblique asymptote(s)");
        if (!result.oblique.empty()) {
            auto slope_val = test_numeric_eval(result.oblique[0].first->simplify());
            auto intercept_val = test_numeric_eval(result.oblique[0].second->simplify());
            if (slope_val) {
                EXPECT_NEAR(*slope_val, 1.0, 1e-9, "oblique asymptote slope = 1");
            } else {
                std::cout << "  oblique slope: " << result.oblique[0].first->to_string() << std::endl;
                EXPECT_TRUE(false, "oblique slope numeric eval failed");
            }
            if (intercept_val) {
                EXPECT_NEAR(*intercept_val, 1.0, 1e-9, "oblique asymptote intercept = 1");
            } else {
                std::cout << "  oblique intercept: " << result.oblique[0].second->to_string() << std::endl;
                EXPECT_TRUE(false, "oblique intercept numeric eval failed");
            }
        }

        // No horizontal asymptotes (since oblique exists)
        EXPECT_TRUE(result.horizontal.empty(),
                    "(x^2+1)/(x-1) has no horizontal asymptotes");
    }


    TEST_CASE("asymptotes_checked: invalid arguments are explicit errors");
    {
        auto null_result = LMCAS::asymptotes_checked(nullptr, "x");
        EXPECT_TRUE(!null_result &&
                        null_result.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked asymptotes rejects null expression");

        auto empty_var = LMCAS::asymptotes_checked(var("x"), "");
        EXPECT_TRUE(!empty_var &&
                        empty_var.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked asymptotes rejects empty variable");

    }

    TEST_CASE("asymptotes_checked: denominator solving uses checked dispatcher");
    {
        auto x = var("x");
        auto denominator = SE::add(x, num(-1));
        auto f = SE::divide(num(1), denominator);

        LMCAS::ResourceLimits limits;
        limits.max_steps = 1;
        LMCAS::ComputationContext limited_context(limits);
        auto limited = LMCAS::asymptotes_checked(f, "x", limited_context);
        EXPECT_TRUE(!limited &&
                        limited.error().code == LMCAS::CasErrc::ResourceLimit,
                    "checked asymptotes propagates denominator solve budget failure");

    }

    TEST_CASE("log_differentiate: x^2");
    {
        // d/dx[x^2] = 2x
        auto x = var("x");
        auto f = SE::power(x, num(2));
        auto result = LMCAS::log_differentiate(f, "x");
        EXPECT_TRUE(result != nullptr, "log_differentiate(x^2) non-null");

        // 验证在 x=3 处的数值: f'(3) = 6
        auto at3 = result->substitute("x", num(3))->simplify();
        std::cout << "  log_differentiate(x^2) at x=3: " << at3->to_string() << std::endl;
        auto val = test_numeric_eval(at3);
        if (val) {
            EXPECT_NEAR(*val, 6.0, 1e-9, "log_differentiate(x^2) at x=3 = 6");
        } else {
            EXPECT_TRUE(false, "log_differentiate(x^2) at x=3 numeric eval failed");
        }
    }

    TEST_CASE("log_differentiate: x^x (variable exponent)");
    {
        // d/dx[x^x] = x^x * (ln(x) + 1)
        auto x = var("x");
        auto f = SE::power(x, x);
        auto result = LMCAS::log_differentiate(f, "x");
        EXPECT_TRUE(result != nullptr, "log_differentiate(x^x) non-null");
        std::cout << "  log_differentiate(x^x): " << result->to_string() << std::endl;

        // 验证在 x=1 处: f'(1) = 1^1 * (ln(1) + 1) = 1 * (0 + 1) = 1
        auto at1 = result->substitute("x", num(1))->simplify();
        std::cout << "  at x=1: " << at1->to_string() << std::endl;
        auto val = test_numeric_eval(at1);
        if (val) {
            EXPECT_NEAR(*val, 1.0, 1e-9, "log_differentiate(x^x) at x=1 = 1");
        } else {
            EXPECT_TRUE(false, "log_differentiate(x^x) at x=1 numeric eval failed");
        }
    }

    TEST_CASE("log_differentiate: product x * (x+1)");
    {
        // d/dx[x*(x+1)] = 2x + 1
        auto x = var("x");
        auto f = SE::multiply(x, SE::add(x, num(1)));
        auto result = LMCAS::log_differentiate(f, "x");
        EXPECT_TRUE(result != nullptr, "log_differentiate(x*(x+1)) non-null");
        std::cout << "  log_differentiate(x*(x+1)): " << result->to_string() << std::endl;

        // 验证在 x=2 处: f'(2) = 2*2 + 1 = 5
        auto at2 = result->substitute("x", num(2))->simplify();
        auto val = test_numeric_eval(at2);
        if (val) {
            EXPECT_NEAR(*val, 5.0, 1e-9, "log_differentiate(x*(x+1)) at x=2 = 5");
        } else {
            EXPECT_TRUE(false, "log_differentiate(x*(x+1)) at x=2 numeric eval failed");
        }
    }

    TEST_CASE("differential: x^3");
    {
        // d(x^3)/dx = 3x^2
        auto x = var("x");
        auto f = SE::power(x, num(3));
        auto result = LMCAS::differential(f, "x");
        EXPECT_TRUE(result != nullptr, "differential(x^3) non-null");
        std::cout << "  differential(x^3): " << result->to_string() << std::endl;

        // 验证在 x=2 处: 3*4 = 12
        auto at2 = result->substitute("x", num(2))->simplify();
        auto val = test_numeric_eval(at2);
        if (val) {
            EXPECT_NEAR(*val, 12.0, 1e-9, "differential(x^3) at x=2 = 12");
        } else {
            EXPECT_TRUE(false, "differential(x^3) at x=2 numeric eval failed");
        }
    }

    TEST_CASE("differential: sin(x)");
    {
        // d(sin(x))/dx = cos(x)
        auto x = var("x");
        auto f = SE::sin(x);
        auto result = LMCAS::differential(f, "x");
        EXPECT_TRUE(result != nullptr, "differential(sin(x)) non-null");
        std::cout << "  differential(sin(x)): " << result->to_string() << std::endl;
    }

    TEST_CASE("total_differential: x^2 + y^2");
    {
        // df = 2x dx + 2y dy
        auto x = var("x");
        auto y = var("y");
        auto f = SE::add(SE::power(x, num(2)), SE::power(y, num(2)));
        auto result = LMCAS::total_differential(f, {"x", "y"});
        EXPECT_TRUE(result.size() == 2, "total_differential has 2 terms");

        // 第一项: ∂f/∂x = 2x, 变量 "x"
        EXPECT_EQ_STR(result[0].second, "x", "first term variable is x");
        std::cout << "  df/dx: " << result[0].first->to_string() << std::endl;

        // 第二项: ∂f/∂y = 2y, 变量 "y"
        EXPECT_EQ_STR(result[1].second, "y", "second term variable is y");
        std::cout << "  df/dy: " << result[1].first->to_string() << std::endl;

        // 验证 ∂f/∂x at (3,4) = 6
        auto dx_at = result[0].first->substitute("x", num(3))->substitute("y", num(4))->simplify();
        auto val_dx = test_numeric_eval(dx_at);
        if (val_dx) {
            EXPECT_NEAR(*val_dx, 6.0, 1e-9, "df/dx at (3,4) = 6");
        }

        // 验证 ∂f/∂y at (3,4) = 8
        auto dy_at = result[1].first->substitute("x", num(3))->substitute("y", num(4))->simplify();
        auto val_dy = test_numeric_eval(dy_at);
        if (val_dy) {
            EXPECT_NEAR(*val_dy, 8.0, 1e-9, "df/dy at (3,4) = 8");
        }
    }

    TEST_CASE("total_differential: x*y*z (3 variables)");
    {
        auto x = var("x");
        auto y = var("y");
        auto z = var("z");
        auto f = SE::multiply(SE::multiply(x, y), z);
        auto result = LMCAS::total_differential(f, {"x", "y", "z"});
        EXPECT_TRUE(result.size() == 3, "total_differential has 3 terms");
        EXPECT_EQ_STR(result[0].second, "x", "first var is x");
        EXPECT_EQ_STR(result[1].second, "y", "second var is y");
        EXPECT_EQ_STR(result[2].second, "z", "third var is z");
    }

    TEST_CASE("inverse_function: f(x) = 2x + 1, solve for y=5");
    {
        // f(x) = 2x + 1, f^{-1}(5) = 2
        auto x = var("x");
        auto f = SE::add(SE::multiply(num(2), x), num(1));
        auto result = LMCAS::inverse_function_checked(f, "x", num(5)).value();
        EXPECT_TRUE(!result.empty(), "inverse_function(2x+1, 5) has solutions");
        if (!result.empty()) {
            auto sol = result[0]->simplify();
            std::cout << "  f^{-1}(5) = " << sol->to_string() << std::endl;
            auto val = test_numeric_eval(sol);
            if (val) {
                EXPECT_NEAR(*val, 2.0, 1e-9, "inverse_function(2x+1, 5) = 2");
            }
        }
    }

    TEST_CASE("inverse_function: f(x) = x^2, solve for y=4");
    {
        // f(x) = x^2, f^{-1}(4) = ±2
        auto x = var("x");
        auto f = SE::power(x, num(2));
        auto result = LMCAS::inverse_function_checked(f, "x", num(4)).value();
        EXPECT_TRUE(result.size() >= 1, "inverse_function(x^2, 4) has solutions");
        std::cout << "  Solutions for x^2 = 4:";
        for (auto& s : result) {
            std::cout << " " << s->to_string();
        }
        std::cout << std::endl;
    }

    TEST_CASE("inverse_function_checked: explicit solution-set semantics");
    {
        auto x = var("x");
        auto linear = SE::add(SE::multiply(num(2), x), num(1));

        auto finite = LMCAS::inverse_function_checked(linear, "x", num(5));
        EXPECT_TRUE(finite && finite.value().size() == 1,
                    "checked inverse function returns finite exact candidates");

        auto empty = LMCAS::inverse_function_checked(num(1), "x", num(2));
        EXPECT_TRUE(empty && empty.value().empty(),
                    "checked inverse function preserves a proven empty set");

        auto null_input = LMCAS::inverse_function_checked(nullptr, "x", num(1));
        EXPECT_TRUE(!null_input &&
                        null_input.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked inverse function rejects null expression");

        auto null_target = LMCAS::inverse_function_checked(x, "x", nullptr);
        EXPECT_TRUE(!null_target &&
                        null_target.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked inverse function rejects null target");

        auto trigonometric = LMCAS::inverse_function_checked(
            SE::sin(x), "x", num(1));
        EXPECT_TRUE(trigonometric && !trigonometric.value().empty(),
                    "checked inverse function preserves supported trigonometric solutions");

        LMCAS::ResourceLimits limits;
        limits.max_steps = 1;
        LMCAS::ComputationContext limited_context(limits);
        auto limited = LMCAS::inverse_function_checked(
            linear, "x", num(5), limited_context);
        EXPECT_TRUE(!limited &&
                        limited.error().code == LMCAS::CasErrc::ResourceLimit,
                    "checked inverse function observes exhausted step budget");
    }

    TEST_CASE("inverse_derivative: f(x) = 2x + 1 at point=5");
    {
        // f(x) = 2x+1, f'(x) = 2, f^{-1}(5) = 2
        // (f^{-1})'(5) = 1/f'(2) = 1/2
        auto x = var("x");
        auto f = SE::add(SE::multiply(num(2), x), num(1));
        auto result = LMCAS::inverse_derivative_checked(f, "x", num(5)).value();
        EXPECT_TRUE(result != nullptr, "inverse_derivative(2x+1, 5) non-null");
        if (result) {
            std::cout << "  (f^{-1})'(5) = " << result->to_string() << std::endl;
            auto val = test_numeric_eval(result);
            if (val) {
                EXPECT_NEAR(*val, 0.5, 1e-9, "inverse_derivative(2x+1, 5) = 1/2");
            }
        }
    }

    TEST_CASE("inverse_derivative: f(x) = x^2 at point=4 is multi-branch");
    {
        // f(x) = x^2 has two inverse branches at y=4. The checked-first
        // legacy wrapper must not silently pick one branch.
        auto x = var("x");
        auto f = SE::power(x, num(2));
        auto result = LMCAS::inverse_derivative_checked(f, "x", num(4));
        EXPECT_TRUE(!result && result.error().code == LMCAS::CasErrc::Inconclusive,
                    "checked inverse derivative rejects multiple branches");
    }

    TEST_CASE("inverse_derivative: f(x) = x^3 at point=8");
    {
        // f(x) = x^3, f'(x) = 3x^2, f^{-1}(8) = 2
        // (f^{-1})'(8) = 1/f'(2) = 1/12
        auto x = var("x");
        auto f = SE::power(x, num(3));
        auto result = LMCAS::inverse_derivative_checked(f, "x", num(8)).value();
        EXPECT_TRUE(result != nullptr, "inverse_derivative(x^3, 8) non-null");
        if (result) {
            std::cout << "  (f^{-1})'(8) = " << result->to_string() << std::endl;
            auto val = test_numeric_eval(result);
            if (val) {
                EXPECT_NEAR(*val, 1.0/12.0, 1e-9, "inverse_derivative(x^3, 8) = 1/12");
            }
        }
    }

    TEST_CASE("inverse_derivative_checked: explicit errors and support domain");
    {
        auto x = var("x");
        auto linear = SE::add(SE::multiply(num(2), x), num(1));

        auto checked = LMCAS::inverse_derivative_checked(linear, "x", num(5));
        EXPECT_TRUE(checked.has_value(),
                    "checked inverse derivative succeeds for a unique linear inverse");
        if (checked) {
            auto val = test_numeric_eval(checked.value());
            EXPECT_TRUE(val.has_value() && std::abs(*val - 0.5) < 1e-9,
                        "checked inverse derivative of 2x+1 is 1/2");
        }

        auto cubic_derivative = LMCAS::inverse_derivative_checked(
            SE::power(x, num(3)), "x", num(8));
        EXPECT_TRUE(cubic_derivative.has_value(),
                    "checked inverse derivative accepts the unique real cubic branch");
        if (cubic_derivative) {
            auto val = test_numeric_eval(cubic_derivative.value());
            EXPECT_TRUE(val.has_value() && std::abs(*val - (1.0 / 12.0)) < 1e-9,
                        "checked inverse derivative of x^3 at 8 is 1/12");
        }

        auto no_preimage = LMCAS::inverse_derivative_checked(num(1), "x", num(2));
        EXPECT_TRUE(!no_preimage &&
                        no_preimage.error().code == LMCAS::CasErrc::DomainError,
                    "checked inverse derivative reports DomainError for no inverse point");

        auto multibranch = LMCAS::inverse_derivative_checked(
            SE::power(x, num(2)), "x", num(4));
        EXPECT_TRUE(!multibranch &&
                        multibranch.error().code == LMCAS::CasErrc::Inconclusive,
                    "checked inverse derivative rejects non-unique inverse branches");

        auto zero_derivative = LMCAS::inverse_derivative_checked(
            SE::power(x, num(3)), "x", num(0));
        EXPECT_TRUE(!zero_derivative &&
                        zero_derivative.error().code == LMCAS::CasErrc::DomainError,
                    "checked inverse derivative reports DomainError when f' is zero");

        auto null_input = LMCAS::inverse_derivative_checked(nullptr, "x", num(1));
        EXPECT_TRUE(!null_input &&
                        null_input.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked inverse derivative rejects null expression");

        LMCAS::ResourceLimits limits;
        limits.max_steps = 1;
        LMCAS::ComputationContext limited_context(limits);
        auto limited = LMCAS::inverse_derivative_checked(
            linear, "x", num(5), limited_context);
        EXPECT_TRUE(!limited &&
                        limited.error().code == LMCAS::CasErrc::ResourceLimit,
                    "checked inverse derivative observes exhausted step budget");
    }

    TEST_CASE("log_differentiate: null input");
    {
        auto result = LMCAS::log_differentiate(nullptr, "x");
        EXPECT_TRUE(result == nullptr, "log_differentiate(null) returns null");
    }

    TEST_CASE("differential: null input");
    {
        auto result = LMCAS::differential(nullptr, "x");
        EXPECT_TRUE(result == nullptr, "differential(null) returns null");
    }

    TEST_CASE("total_differential: empty vars");
    {
        auto x = var("x");
        auto result = LMCAS::total_differential(x, {});
        EXPECT_TRUE(result.empty(), "total_differential with empty vars returns empty");
    }


    TEST_CASE("curvature_checked: invalid arguments are explicit errors");
    {
        auto null_result = LMCAS::curvature_checked(nullptr, "x");
        EXPECT_TRUE(!null_result &&
                        null_result.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked curvature rejects null expression");

        auto empty_var = LMCAS::curvature_checked(var("x"), "");
        EXPECT_TRUE(!empty_var &&
                        empty_var.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked curvature rejects empty variable");

    }

    TEST_CASE("curvature_parametric_checked: zero velocity is a domain error");
    {
        auto checked = LMCAS::curvature_parametric_checked(num(1), num(2), "t");
        EXPECT_TRUE(!checked &&
                        checked.error().code == LMCAS::CasErrc::DomainError,
                    "checked parametric curvature rejects zero velocity");

    }

    TEST_CASE("surface_area_revolution_x: huge exact numeric fallback fails safely");
    {
        auto x = var("x");
        auto unsupported = SE::exp(SE::power(x, num(2)));
        const BigInt huge("1" + std::string(400, '0'));
        auto result = LMCAS::surface_area_revolution_x_checked(
            unsupported, "x", bigint_num(BigInt(0)), bigint_num(huge));
        EXPECT_TRUE(!result,
                    "huge exact bounds do not fabricate a surface area");
    }

    TEST_CASE("surface_area_revolution_checked: exact results and explicit failures");
    {
        auto x = var("x");
        auto zero = num(0);
        auto one = num(1);

        auto exact = LMCAS::surface_area_revolution_x_checked(num(1), "x", zero, one);
        EXPECT_TRUE(exact.has_value(),
                    "checked x-axis surface area succeeds for constant radius");
        if (exact) {
            auto at_pi = exact.value()->substitute("pi", num(3))->simplify();
            auto val = test_numeric_eval(at_pi);
            EXPECT_TRUE(val.has_value() && std::abs(*val - 6.0) < 1e-9,
                        "checked constant surface area is 2*pi");
        }

        auto null_result = LMCAS::surface_area_revolution_x_checked(
            nullptr, "x", zero, one);
        EXPECT_TRUE(!null_result &&
                        null_result.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked surface area rejects null expression");

        auto null_bound = LMCAS::surface_area_revolution_y_checked(
            x, "x", zero, nullptr);
        EXPECT_TRUE(!null_bound &&
                        null_bound.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked surface area rejects null bounds");

        LMCAS::ResourceLimits limits;
        limits.max_steps = 1;
        LMCAS::ComputationContext limited_context(limits);
        auto limited = LMCAS::surface_area_revolution_x_checked(
            x, "x", zero, one, limited_context);
        EXPECT_TRUE(!limited &&
                        limited.error().code == LMCAS::CasErrc::ResourceLimit,
                    "checked surface area observes exhausted step budget");
    }

    TEST_CASE("surface_area_revolution_checked: implicit numeric fallback is inconclusive");
    {
        auto x = var("x");
        auto zero = num(0);
        auto one = num(1);
        auto unsupported = SE::exp(SE::power(x, num(2)));


        auto checked = LMCAS::surface_area_revolution_x_checked(
            unsupported, "x", zero, one);
        EXPECT_TRUE(!checked &&
                        checked.error().code == LMCAS::CasErrc::Inconclusive,
                    "checked surface area rejects implicit numeric fallback");
    }

    TEST_CASE("inflection_points_checked: invalid arguments are explicit errors");
    {
        auto null_result = LMCAS::inflection_points_checked(nullptr, "x");
        EXPECT_TRUE(!null_result &&
                        null_result.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked inflection points rejects null expression");

        auto empty_var = LMCAS::inflection_points_checked(var("x"), "");
        EXPECT_TRUE(!empty_var &&
                        empty_var.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked inflection points rejects empty variable");

    }

    TEST_CASE("inflection_points_checked: distinguishes empty, finite, and inconclusive");
    {
        auto x = var("x");

        auto no_inflection = LMCAS::inflection_points_checked(SE::power(x, num(2)), "x");
        EXPECT_TRUE(no_inflection && no_inflection.value().empty(),
                    "checked inflection points reports an empty set for f'' = constant nonzero");

        auto cubic = LMCAS::inflection_points_checked(SE::power(x, num(3)), "x");
        EXPECT_TRUE(cubic && cubic.value().size() == 1,
                    "checked inflection points reports finite exact candidates for x^3");

        auto exponential = LMCAS::inflection_points_checked(SE::exp(x), "x");
        EXPECT_TRUE(exponential && exponential.value().empty(),
                    "checked inflection points proves exp(x) has no real inflection point");

        LMCAS::ResourceLimits limits;
        limits.max_steps = 1;
        LMCAS::ComputationContext limited_context(limits);
        auto limited = LMCAS::inflection_points_checked(
            SE::power(x, num(3)), "x", limited_context);
        EXPECT_TRUE(!limited &&
                        limited.error().code == LMCAS::CasErrc::ResourceLimit,
                    "checked inflection points propagates solver budget failure");
    }

    return TEST_REPORT();
}
