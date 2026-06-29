/**
 * @file test_calculus_utils.cpp
 * @brief 测试 calculus_utils 模块：对数微分、微分、全微分、反函数导数、反函数。
 */

#include "test_common.hpp"
#include "calculus_utils.hpp"

using SE = SymbolicExpr;

static auto num(int n) { return SE::number(n); }
static auto var(const std::string& name) { return SE::variable(name); }

int main() {
    // =========================================================
    // continuity_at 测试 (Requirement 7)
    // =========================================================
    TEST_CASE("continuity_at: polynomial is continuous everywhere");
    {
        // f(x) = x^2 is continuous at x=1
        auto x = var("x");
        auto f = SE::power(x, num(2));
        auto result = lamina::continuity_at(f, "x", num(1));
        EXPECT_TRUE(result == lamina::ContinuityType::Continuous,
                    "x^2 is continuous at x=1");
    }

    TEST_CASE("continuity_at: 1/x has essential discontinuity at x=0");
    {
        // f(x) = 1/x has essential discontinuity at x=0 (limit is ±∞)
        auto x = var("x");
        auto f = SE::divide(num(1), x);
        auto result = lamina::continuity_at(f, "x", num(0));
        EXPECT_TRUE(result == lamina::ContinuityType::Essential,
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
        auto result = lamina::continuity_at(f, "x", num(1));
        EXPECT_TRUE(result == lamina::ContinuityType::Removable,
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
        auto result = lamina::continuity_at(f, "x", num(0));
        EXPECT_TRUE(result == lamina::ContinuityType::Jump,
                    "1/(1+e^(1/x)) has jump discontinuity at x=0");
    }

    TEST_CASE("continuity_at: null input returns Essential");
    {
        auto result = lamina::continuity_at(nullptr, "x", num(0));
        EXPECT_TRUE(result == lamina::ContinuityType::Essential,
                    "null input returns Essential");
    }

    // =========================================================
    // asymptotes 测试 (Requirement 15)
    // =========================================================
    TEST_CASE("asymptotes: 1/x has vertical at x=0, horizontal at y=0");
    {
        // f(x) = 1/x: vertical asymptote at x=0, horizontal at y=0
        auto x = var("x");
        auto f = SE::divide(num(1), x);
        auto result = lamina::asymptotes(f, "x");

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
        auto result = lamina::asymptotes(f, "x");

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
        auto result = lamina::asymptotes(f, "x");

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

    TEST_CASE("asymptotes: null input returns empty result");
    {
        auto result = lamina::asymptotes(nullptr, "x");
        EXPECT_TRUE(result.vertical.empty(), "null: no vertical");
        EXPECT_TRUE(result.horizontal.empty(), "null: no horizontal");
        EXPECT_TRUE(result.oblique.empty(), "null: no oblique");
    }

    // =========================================================
    // log_differentiate 测试
    // =========================================================
    TEST_CASE("log_differentiate: x^2");
    {
        // d/dx[x^2] = 2x
        auto x = var("x");
        auto f = SE::power(x, num(2));
        auto result = lamina::log_differentiate(f, "x");
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
        auto result = lamina::log_differentiate(f, "x");
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
        auto result = lamina::log_differentiate(f, "x");
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

    // =========================================================
    // differential 测试
    // =========================================================
    TEST_CASE("differential: x^3");
    {
        // d(x^3)/dx = 3x^2
        auto x = var("x");
        auto f = SE::power(x, num(3));
        auto result = lamina::differential(f, "x");
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
        auto result = lamina::differential(f, "x");
        EXPECT_TRUE(result != nullptr, "differential(sin(x)) non-null");
        std::cout << "  differential(sin(x)): " << result->to_string() << std::endl;
    }

    // =========================================================
    // total_differential 测试
    // =========================================================
    TEST_CASE("total_differential: x^2 + y^2");
    {
        // df = 2x dx + 2y dy
        auto x = var("x");
        auto y = var("y");
        auto f = SE::add(SE::power(x, num(2)), SE::power(y, num(2)));
        auto result = lamina::total_differential(f, {"x", "y"});
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
        auto result = lamina::total_differential(f, {"x", "y", "z"});
        EXPECT_TRUE(result.size() == 3, "total_differential has 3 terms");
        EXPECT_EQ_STR(result[0].second, "x", "first var is x");
        EXPECT_EQ_STR(result[1].second, "y", "second var is y");
        EXPECT_EQ_STR(result[2].second, "z", "third var is z");
    }

    // =========================================================
    // inverse_function 测试
    // =========================================================
    TEST_CASE("inverse_function: f(x) = 2x + 1, solve for y=5");
    {
        // f(x) = 2x + 1, f^{-1}(5) = 2
        auto x = var("x");
        auto f = SE::add(SE::multiply(num(2), x), num(1));
        auto result = lamina::inverse_function(f, "x", num(5));
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
        auto result = lamina::inverse_function(f, "x", num(4));
        EXPECT_TRUE(result.size() >= 1, "inverse_function(x^2, 4) has solutions");
        std::cout << "  Solutions for x^2 = 4:";
        for (auto& s : result) {
            std::cout << " " << s->to_string();
        }
        std::cout << std::endl;
    }

    // =========================================================
    // inverse_derivative 测试
    // =========================================================
    TEST_CASE("inverse_derivative: f(x) = 2x + 1 at point=5");
    {
        // f(x) = 2x+1, f'(x) = 2, f^{-1}(5) = 2
        // (f^{-1})'(5) = 1/f'(2) = 1/2
        auto x = var("x");
        auto f = SE::add(SE::multiply(num(2), x), num(1));
        auto result = lamina::inverse_derivative(f, "x", num(5));
        EXPECT_TRUE(result != nullptr, "inverse_derivative(2x+1, 5) non-null");
        if (result) {
            std::cout << "  (f^{-1})'(5) = " << result->to_string() << std::endl;
            auto val = test_numeric_eval(result);
            if (val) {
                EXPECT_NEAR(*val, 0.5, 1e-9, "inverse_derivative(2x+1, 5) = 1/2");
            }
        }
    }

    TEST_CASE("inverse_derivative: f(x) = x^2 at point=4");
    {
        // f(x) = x^2, f'(x) = 2x, f^{-1}(4) = 2 (first solution)
        // (f^{-1})'(4) = 1/f'(2) = 1/4
        auto x = var("x");
        auto f = SE::power(x, num(2));
        auto result = lamina::inverse_derivative(f, "x", num(4));
        EXPECT_TRUE(result != nullptr, "inverse_derivative(x^2, 4) non-null");
        if (result) {
            std::cout << "  (f^{-1})'(4) = " << result->to_string() << std::endl;
            auto val = test_numeric_eval(result);
            if (val) {
                EXPECT_NEAR(*val, 0.25, 1e-9, "inverse_derivative(x^2, 4) = 1/4");
            }
        }
    }

    TEST_CASE("inverse_derivative: f(x) = x^3 at point=8");
    {
        // f(x) = x^3, f'(x) = 3x^2, f^{-1}(8) = 2
        // (f^{-1})'(8) = 1/f'(2) = 1/12
        auto x = var("x");
        auto f = SE::power(x, num(3));
        auto result = lamina::inverse_derivative(f, "x", num(8));
        EXPECT_TRUE(result != nullptr, "inverse_derivative(x^3, 8) non-null");
        if (result) {
            std::cout << "  (f^{-1})'(8) = " << result->to_string() << std::endl;
            auto val = test_numeric_eval(result);
            if (val) {
                EXPECT_NEAR(*val, 1.0/12.0, 1e-9, "inverse_derivative(x^3, 8) = 1/12");
            }
        }
    }

    // =========================================================
    // Edge cases
    // =========================================================
    TEST_CASE("log_differentiate: null input");
    {
        auto result = lamina::log_differentiate(nullptr, "x");
        EXPECT_TRUE(result == nullptr, "log_differentiate(null) returns null");
    }

    TEST_CASE("differential: null input");
    {
        auto result = lamina::differential(nullptr, "x");
        EXPECT_TRUE(result == nullptr, "differential(null) returns null");
    }

    TEST_CASE("total_differential: empty vars");
    {
        auto x = var("x");
        auto result = lamina::total_differential(x, {});
        EXPECT_TRUE(result.empty(), "total_differential with empty vars returns empty");
    }

    TEST_CASE("inverse_function: null inputs");
    {
        auto result = lamina::inverse_function(nullptr, "x", num(1));
        EXPECT_TRUE(result.empty(), "inverse_function(null) returns empty");
    }

    TEST_CASE("inverse_derivative: null inputs");
    {
        auto result = lamina::inverse_derivative(nullptr, "x", num(1));
        EXPECT_TRUE(result == nullptr, "inverse_derivative(null) returns null");
    }

    return TEST_REPORT();
}
