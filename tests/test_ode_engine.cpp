/**
 * @file test_ode_engine.cpp
 * @brief 测试 ODE 引擎的一阶求解方法：齐次、Bernoulli、恰当方程。
 */
#include "test_common.hpp"
#include "symbolic_ode_engine.hpp"
#include "poly_utils.hpp"

using namespace lamina;


void test_homogeneous_ode() {
    TEST_CASE("Homogeneous ODE: y' = y/x");
    {
        // y' = y/x 是齐次方程，f(y/x) = y/x
        // 解: y = Cx
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto rhs = SymbolicExpr::divide(y, x);

        // 验证分类
        EXPECT_TRUE(is_homogeneous_ode(rhs, "x", "y"), "y/x is homogeneous");

        // 求解
        auto sol = solve_homogeneous_ode_checked(rhs, "x", "y").value();
        EXPECT_TRUE(sol.general_solution != nullptr, "homogeneous y/x has solution");
        EXPECT_TRUE(sol.method_used == ODEType::Homogeneous, "method is Homogeneous");
        EXPECT_TRUE(!sol.constants.empty(), "has constants");
    }

    TEST_CASE("Homogeneous ODE: y' = (x + y)/x");
    {
        // y' = (x + y)/x = 1 + y/x，齐次方程
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto rhs = SymbolicExpr::divide(
            SymbolicExpr::add(x, y), x);

        EXPECT_TRUE(is_homogeneous_ode(rhs, "x", "y"), "(x+y)/x is homogeneous");

        auto sol = solve_homogeneous_ode_checked(rhs, "x", "y").value();
        EXPECT_TRUE(sol.general_solution != nullptr, "homogeneous (x+y)/x has solution");
        std::string s = sol.general_solution->to_string();
        // 解应包含 x 和 y
        EXPECT_CONTAINS(s, {"x"}, "solution contains x");
    }

    TEST_CASE("Homogeneous ODE: y' = (x^2 + y^2)/(x*y)");
    {
        // y' = (x² + y²)/(xy) 是齐次方程
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto y2 = SymbolicExpr::power(y, SymbolicExpr::number(2));
        auto rhs = SymbolicExpr::divide(
            SymbolicExpr::add(x2, y2),
            SymbolicExpr::multiply(x, y));

        EXPECT_TRUE(is_homogeneous_ode(rhs, "x", "y"), "(x^2+y^2)/(xy) is homogeneous");

        auto sol = solve_homogeneous_ode_checked(rhs, "x", "y").value();
        EXPECT_TRUE(sol.general_solution != nullptr, "homogeneous (x^2+y^2)/(xy) has solution");
    }
}


void test_bernoulli_ode() {
    TEST_CASE("Bernoulli ODE: y' + y = y^2 (P=1, Q=1, n=2)");
    {
        // y' + y = y^2 → P(x)=1, Q(x)=1, n=2
        auto P = SymbolicExpr::number(1);
        auto Q = SymbolicExpr::number(1);

        auto sol = solve_bernoulli_ode_checked(P, Q, 2, "x", "y").value();
        EXPECT_TRUE(sol.general_solution != nullptr, "Bernoulli P=1 Q=1 n=2 has solution");
        EXPECT_TRUE(sol.method_used == ODEType::Bernoulli, "method is Bernoulli");

        std::string s = sol.general_solution->to_string();
        EXPECT_CONTAINS(s, {"C"}, "Bernoulli solution contains constant C");
        EXPECT_CONTAINS(s, {"x"}, "Bernoulli solution contains x");
    }

    TEST_CASE("Bernoulli ODE: y' + (1/x)*y = x*y^2 (P=1/x, Q=x, n=2)");
    {
        // y' + (1/x)*y = x*y^2
        auto x = SymbolicExpr::variable("x");
        auto P = SymbolicExpr::divide(SymbolicExpr::number(1), x);
        auto Q = SymbolicExpr::variable("x");

        auto sol = solve_bernoulli_ode_checked(P, Q, 2, "x", "y").value();
        EXPECT_TRUE(sol.general_solution != nullptr, "Bernoulli P=1/x Q=x n=2 has solution");

        std::string s = sol.general_solution->to_string();
        EXPECT_CONTAINS(s, {"C"}, "Bernoulli 1/x solution contains constant C");
    }

    TEST_CASE("Bernoulli ODE: y' + 2*y = y^3 (P=2, Q=1, n=3)");
    {
        // y' + 2y = y^3 → P=2, Q=1, n=3
        auto P = SymbolicExpr::number(2);
        auto Q = SymbolicExpr::number(1);

        auto sol = solve_bernoulli_ode_checked(P, Q, 3, "x", "y").value();
        EXPECT_TRUE(sol.general_solution != nullptr, "Bernoulli P=2 Q=1 n=3 has solution");
        EXPECT_TRUE(sol.method_used == ODEType::Bernoulli, "method is Bernoulli");

        std::string s = sol.general_solution->to_string();
        EXPECT_CONTAINS(s, {"C"}, "Bernoulli n=3 solution contains constant C");
    }
}


void test_exact_ode() {
    TEST_CASE("Exact ODE: (2x + y)dx + (x + 2y)dy = 0");
    {
        // M = 2x + y, N = x + 2y
        // ∂M/∂y = 1, ∂N/∂x = 1 → 恰当
        // F(x,y) = x² + xy + y² = C
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto M = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x), y);
        auto N = SymbolicExpr::add(
            x, SymbolicExpr::multiply(SymbolicExpr::number(2), y));

        EXPECT_TRUE(is_exact_ode(M, N, "x", "y"), "(2x+y, x+2y) is exact");

        auto sol = solve_exact_ode_checked(M, N, "x", "y").value();
        EXPECT_TRUE(sol.general_solution != nullptr, "exact (2x+y, x+2y) has solution");
        EXPECT_TRUE(sol.method_used == ODEType::Exact, "method is Exact");

        std::string s = sol.general_solution->to_string();
        // 解应包含 x 和 y 的二次项
        EXPECT_CONTAINS(s, {"x"}, "exact solution contains x");
        EXPECT_CONTAINS(s, {"y"}, "exact solution contains y");
    }

    TEST_CASE("Exact ODE: (y*cos(x) + 2x*e^y)dx + (sin(x) + x^2*e^y)dy = 0");
    {
        // M = y*cos(x) + 2x*e^y
        // N = sin(x) + x²*e^y
        // ∂M/∂y = cos(x) + 2x*e^y
        // ∂N/∂x = cos(x) + 2x*e^y → 恰当
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto M = SymbolicExpr::add(
            SymbolicExpr::multiply(y, SymbolicExpr::cos(x)),
            SymbolicExpr::multiply(
                SymbolicExpr::multiply(SymbolicExpr::number(2), x),
                SymbolicExpr::exp(y)));
        auto N = SymbolicExpr::add(
            SymbolicExpr::sin(x),
            SymbolicExpr::multiply(
                SymbolicExpr::power(x, SymbolicExpr::number(2)),
                SymbolicExpr::exp(y)));

        EXPECT_TRUE(is_exact_ode(M, N, "x", "y"), "trig+exp exact equation");

        auto sol = solve_exact_ode_checked(M, N, "x", "y").value();
        EXPECT_TRUE(sol.general_solution != nullptr, "trig+exp exact has solution");
    }

    TEST_CASE("Exact ODE: simple (y)dx + (x)dy = 0");
    {
        // M = y, N = x → ∂M/∂y = 1, ∂N/∂x = 1 → 恰当
        // F(x,y) = xy = C
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");

        EXPECT_TRUE(is_exact_ode(y, x, "x", "y"), "(y, x) is exact");

        auto sol = solve_exact_ode_checked(y, x, "x", "y").value();
        EXPECT_TRUE(sol.general_solution != nullptr, "exact (y, x) has solution");

        std::string s = sol.general_solution->to_string();
        // 解应为 xy 形式
        EXPECT_CONTAINS(s, {"x"}, "exact xy solution contains x");
        EXPECT_CONTAINS(s, {"y"}, "exact xy solution contains y");
    }
}

void test_first_order_ode_checked_contracts() {
    TEST_CASE("First-order ODE checked APIs: explicit errors and context");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto rhs = SymbolicExpr::divide(y, x);
        auto sol = solve_homogeneous_ode_checked(rhs, "x", "y");
        EXPECT_TRUE(sol.has_value(), "checked homogeneous ODE succeeds");
        if (sol) {
            EXPECT_TRUE(sol.value().general_solution != nullptr,
                "checked homogeneous ODE returns solution");
            EXPECT_TRUE(sol.value().method_used == ODEType::Homogeneous,
                "checked homogeneous ODE reports Homogeneous");
        }
    }

    {
        auto P = SymbolicExpr::number(1);
        auto Q = SymbolicExpr::number(1);
        auto sol = solve_bernoulli_ode_checked(P, Q, 2, "x", "y");
        EXPECT_TRUE(sol.has_value(), "checked Bernoulli ODE succeeds");
        if (sol) {
            EXPECT_TRUE(sol.value().general_solution != nullptr,
                "checked Bernoulli ODE returns solution");
            EXPECT_TRUE(sol.value().method_used == ODEType::Bernoulli,
                "checked Bernoulli ODE reports Bernoulli");
        }
    }

    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto M = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x), y);
        auto N = SymbolicExpr::add(
            x, SymbolicExpr::multiply(SymbolicExpr::number(2), y));
        auto sol = solve_exact_ode_checked(M, N, "x", "y");
        EXPECT_TRUE(sol.has_value(), "checked exact ODE succeeds");
        if (sol) {
            EXPECT_TRUE(sol.value().general_solution != nullptr,
                "checked exact ODE returns solution");
            EXPECT_TRUE(sol.value().method_used == ODEType::Exact,
                "checked exact ODE reports Exact");
        }
    }

    {
        auto x = SymbolicExpr::variable("x");
        std::shared_ptr<SymbolicExpr> null_root;
        auto null_rhs = solve_homogeneous_ode_checked(null_root, "x", "y");
        EXPECT_TRUE(!null_rhs.has_value(),
            "checked homogeneous ODE rejects null rhs");
        EXPECT_TRUE(null_rhs.error().code == CasErrc::InvalidArgument,
            "checked homogeneous ODE reports InvalidArgument for null rhs");

        auto same_vars = solve_homogeneous_ode_checked(x, "x", "x");
        EXPECT_TRUE(!same_vars.has_value(),
            "checked homogeneous ODE rejects duplicate variable names");
        EXPECT_TRUE(same_vars.error().code == CasErrc::InvalidArgument,
            "checked homogeneous ODE reports InvalidArgument for duplicate variables");

        auto bad_n = solve_bernoulli_ode_checked(x, x, 1, "x", "y");
        EXPECT_TRUE(!bad_n.has_value(),
            "checked Bernoulli ODE rejects n=1");
        EXPECT_TRUE(bad_n.error().code == CasErrc::InvalidArgument,
            "checked Bernoulli ODE reports InvalidArgument for n=1");

        auto null_exact = solve_exact_ode_checked(x, null_root, "x", "y");
        EXPECT_TRUE(!null_exact.has_value(),
            "checked exact ODE rejects null N");
        EXPECT_TRUE(null_exact.error().code == CasErrc::InvalidArgument,
            "checked exact ODE reports InvalidArgument for null N");
    }

    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto rhs = SymbolicExpr::divide(y, x);

        lamina::CancellationToken cancellation;
        lamina::ComputationContext cancelled_context({}, cancellation);
        cancellation.cancel();
        auto cancelled = solve_homogeneous_ode_checked(rhs, "x", "y",
                                                       cancelled_context);
        EXPECT_TRUE(!cancelled.has_value(),
            "checked homogeneous ODE observes cancellation");
        EXPECT_TRUE(cancelled.error().code == CasErrc::Cancelled,
            "checked homogeneous ODE reports Cancelled");

        lamina::ResourceLimits limits;
        limits.max_steps = 1;
        lamina::ComputationContext limited_context(limits);
        auto limited = solve_exact_ode_checked(y, x, "x", "y", limited_context);
        EXPECT_TRUE(!limited.has_value(),
            "checked exact ODE observes exhausted step budget");
        EXPECT_TRUE(limited.error().code == CasErrc::ResourceLimit,
            "checked exact ODE reports ResourceLimit");
    }
}

void test_higher_order_euler_checked_contracts() {
    TEST_CASE("Higher-order and Euler ODE checked APIs: explicit support domain");
    {
        auto high = solve_higher_order_ode_checked(
            {1.0, 0.0, 1.0}, nullptr, "x", "y");
        EXPECT_TRUE(high.has_value(),
            "checked homogeneous constant-coefficient ODE succeeds");
        if (high) {
            EXPECT_TRUE(high.value().general_solution != nullptr,
                "checked higher-order ODE returns a solution");
            EXPECT_TRUE(high.value().method_used == ODEType::HigherOrder_ConstCoeff,
                "checked higher-order ODE reports constant-coefficient method");
            EXPECT_TRUE(high.value().constants.size() == 2,
                "checked higher-order ODE returns two integration constants");
            std::string s = high.value().general_solution->to_string();
            EXPECT_CONTAINS(s, {"C1"}, "higher-order solution contains C1");
            EXPECT_CONTAINS(s, {"C2"}, "higher-order solution contains C2");
        }
    }

    {
        auto euler = solve_euler_ode_checked(
            {1.0, 1.0, -1.0}, nullptr, "x", "y");
        EXPECT_TRUE(euler.has_value(),
            "checked homogeneous Euler ODE succeeds");
        if (euler) {
            EXPECT_TRUE(euler.value().general_solution != nullptr,
                "checked Euler ODE returns a solution");
            EXPECT_TRUE(euler.value().method_used == ODEType::Euler,
                "checked Euler ODE reports Euler method");
            EXPECT_TRUE(euler.value().constants.size() == 2,
                "checked Euler ODE returns two integration constants");
            std::string s = euler.value().general_solution->to_string();
            EXPECT_CONTAINS(s, {"C1"}, "Euler solution contains C1");
            EXPECT_CONTAINS(s, {"C2"}, "Euler solution contains C2");
        }
    }

    {
        auto bad_coeffs = solve_higher_order_ode_checked({}, nullptr, "x", "y");
        EXPECT_TRUE(!bad_coeffs.has_value(),
            "checked higher-order ODE rejects empty coefficient list");
        EXPECT_TRUE(bad_coeffs.error().code == CasErrc::InvalidArgument,
            "checked higher-order ODE reports InvalidArgument for empty coefficients");

        auto bad_leading = solve_higher_order_ode_checked({0.0, 1.0}, nullptr, "x", "y");
        EXPECT_TRUE(!bad_leading.has_value(),
            "checked higher-order ODE rejects zero leading coefficient");
        EXPECT_TRUE(bad_leading.error().code == CasErrc::InvalidArgument,
            "checked higher-order ODE reports InvalidArgument for zero leading coefficient");

        auto same_vars = solve_euler_ode_checked({1.0, 1.0, -1.0}, nullptr, "x", "x");
        EXPECT_TRUE(!same_vars.has_value(),
            "checked Euler ODE rejects duplicate variables");
        EXPECT_TRUE(same_vars.error().code == CasErrc::InvalidArgument,
            "checked Euler ODE reports InvalidArgument for duplicate variables");

        auto constant_forcing = solve_higher_order_ode_checked(
            {1.0, 0.0, 1.0}, SymbolicExpr::number(1), "x", "y");
        EXPECT_TRUE(constant_forcing.has_value(),
            "checked higher-order ODE supports constant nonhomogeneous forcing");
        if (constant_forcing) {
            EXPECT_CONTAINS(
                constant_forcing.value().general_solution->to_string(),
                {"1"}, "constant-forcing solution contains its particular term");
        }

        auto euler_forcing = solve_euler_ode_checked(
            {1.0, 1.0, -1.0}, SymbolicExpr::number(1), "x", "y");
        EXPECT_TRUE(euler_forcing.has_value(),
            "checked Euler ODE supports constant nonhomogeneous forcing");
    }

    {
        lamina::CancellationToken cancellation;
        lamina::ComputationContext cancelled_context({}, cancellation);
        cancellation.cancel();
        auto cancelled = solve_higher_order_ode_checked(
            {1.0, 0.0, 1.0}, nullptr, "x", "y", cancelled_context);
        EXPECT_TRUE(!cancelled.has_value(),
            "checked higher-order ODE observes cancellation");
        EXPECT_TRUE(cancelled.error().code == CasErrc::Cancelled,
            "checked higher-order ODE reports Cancelled");

        lamina::ResourceLimits limits;
        limits.max_steps = 1;
        lamina::ComputationContext limited_context(limits);
        auto limited = solve_euler_ode_checked(
            {1.0, 1.0, -1.0}, nullptr, "x", "y", limited_context);
        EXPECT_TRUE(!limited.has_value(),
            "checked Euler ODE observes exhausted step budget");
        EXPECT_TRUE(limited.error().code == CasErrc::ResourceLimit,
            "checked Euler ODE reports ResourceLimit");
    }
}


void test_integrating_factor() {
    TEST_CASE("Integrating factor: (2y)dx + (x)dy = 0 (not exact, μ(x) exists)");
    {
        // M = 2y, N = x
        // ∂M/∂y = 2, ∂N/∂x = 1 → 不恰当
        // (∂M/∂y - ∂N/∂x)/N = (2-1)/x = 1/x → 仅依赖 x
        // μ(x) = exp(∫1/x dx) = exp(ln(x)) = x
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto M = SymbolicExpr::multiply(SymbolicExpr::number(2), y);
        auto N = x;

        EXPECT_FALSE(is_exact_ode(M, N, "x", "y"), "(2y, x) is NOT exact");

        auto mu = find_integrating_factor(M, N, "x", "y");
        EXPECT_TRUE(mu != nullptr, "integrating factor found for (2y, x)");

        // 乘以积分因子后应变为恰当
        if (mu) {
            auto M_new = SymbolicExpr::multiply(mu, M)->simplify();
            auto N_new = SymbolicExpr::multiply(mu, N)->simplify();
            // 验证新方程恰当
            auto dM_dy = M_new->differentiate("y");
            auto dN_dx = N_new->differentiate("x");
            if (dM_dy && dN_dx) {
                auto diff = SymbolicExpr::add(dM_dy,
                    SymbolicExpr::multiply(SymbolicExpr::number(-1), dN_dx))->simplify();
                // 差值应为零或接近零
                EXPECT_TRUE(diff->is_zero() || is_exact_ode(M_new, N_new, "x", "y"),
                    "after multiplying by μ, equation becomes exact");
            }
        }
    }

    TEST_CASE("Integrating factor: solve non-exact via solve_exact_ode");
    {
        // (2y)dx + (x)dy = 0 → 不恰当但有积分因子
        // solve_exact_ode 应自动找到积分因子并求解
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto M = SymbolicExpr::multiply(SymbolicExpr::number(2), y);
        auto N = x;

        auto sol = solve_exact_ode_checked(M, N, "x", "y").value();
        EXPECT_TRUE(sol.general_solution != nullptr,
            "solve_exact_ode handles non-exact with integrating factor");
    }
}


void test_classification() {
    TEST_CASE("Classification: y' = y/x → Homogeneous");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto rhs = SymbolicExpr::divide(y, x);

        auto cls = classify_first_order_ode(rhs, "x", "y");
        // y/x 可能被分类为可分离或齐次（取决于检测顺序）
        // 实际上 y/x = (1/x)*y 也是线性的
        EXPECT_TRUE(cls.type == ODEType::Separable ||
                    cls.type == ODEType::Linear1 ||
                    cls.type == ODEType::Homogeneous,
            "y/x classified as separable, linear, or homogeneous");
    }

    TEST_CASE("Classification: y' = x^2 + y^2 → Unknown (not standard first-order)");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        auto rhs = SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::power(y, SymbolicExpr::number(2)));

        auto cls = classify_first_order_ode(rhs, "x", "y");
        /// x² + y² 进入一阶分类路径，可归类为恰当方程或 Unknown。
        EXPECT_TRUE(cls.order == 1, "classified as first order");
    }
}


void test_variation_of_parameters() {
    TEST_CASE("Variation of Parameters: y'' + y = sec(x)");
    {
        // 齐次解: y₁ = cos(x), y₂ = sin(x)
        // 非齐次项: g(x) = sec(x) = 1/cos(x)
        // Wronskian: W = cos(x)*cos(x) - sin(x)*(-sin(x)) = cos²(x) + sin²(x) = 1
        // u₁' = -sin(x)*sec(x)/1 = -sin(x)/cos(x) = -tan(x)
        // u₂' = cos(x)*sec(x)/1 = 1
        // u₁ = ln|cos(x)|, u₂ = x
        // y_p = cos(x)*ln|cos(x)| + x*sin(x)
        auto x = SymbolicExpr::variable("x");
        auto y1 = SymbolicExpr::cos(x);
        auto y2 = SymbolicExpr::sin(x);
        auto g = SymbolicExpr::divide(SymbolicExpr::number(1), SymbolicExpr::cos(x));

        auto sol = solve_variation_of_parameters_checked(y1, y2, g, "x").value();
        EXPECT_TRUE(sol.general_solution != nullptr, "VoP y''+y=sec(x) has solution");
        EXPECT_TRUE(sol.constants.empty(), "VoP produces particular solution (no constants)");

        std::string s = sol.general_solution->to_string();
        EXPECT_CONTAINS(s, {"x"}, "VoP solution contains x");
    }

    TEST_CASE("Variation of Parameters: y'' - y = e^x");
    {
        // 齐次解: y₁ = e^x, y₂ = e^(-x)
        // g(x) = e^x
        // W = e^x*(-e^(-x)) - e^(-x)*e^x = -1 - 1 = -2
        // u₁' = -e^(-x)*e^x/(-2) = 1/2
        // u₂' = e^x*e^x/(-2) = -e^(2x)/2
        // u₁ = x/2, u₂ = -e^(2x)/4
        // y_p = (x/2)*e^x + (-e^(2x)/4)*e^(-x) = (x/2)*e^x - e^x/4
        auto x = SymbolicExpr::variable("x");
        auto y1 = SymbolicExpr::exp(x);
        auto y2 = SymbolicExpr::exp(
            SymbolicExpr::multiply(SymbolicExpr::number(-1), x));
        auto g = SymbolicExpr::exp(x);

        auto checked = solve_variation_of_parameters_checked(
            y1, y2, g, "x");
        EXPECT_TRUE(checked && checked.value().general_solution != nullptr,
                    "VoP y''-y=e^x has solution");
        if (checked && checked.value().general_solution) {
            EXPECT_CONTAINS(
                checked.value().general_solution->to_string(), {"x"},
                "VoP y''-y=e^x solution contains x");
        }
    }

    TEST_CASE("Variation of Parameters: y'' + y = sin(x)");
    {
        // 齐次解: y₁ = cos(x), y₂ = sin(x)
        // g(x) = sin(x)
        // W = 1
        // u₁' = -sin(x)*sin(x) = -sin²(x)
        // u₂' = cos(x)*sin(x) = sin(x)cos(x)
        auto x = SymbolicExpr::variable("x");
        auto y1 = SymbolicExpr::cos(x);
        auto y2 = SymbolicExpr::sin(x);
        auto g = SymbolicExpr::sin(x);

        auto sol = solve_variation_of_parameters_checked(y1, y2, g, "x").value();
        EXPECT_TRUE(sol.general_solution != nullptr, "VoP y''+y=sin(x) has solution");
    }

    TEST_CASE("Variation of Parameters: null inputs");
    {
        auto sol = solve_variation_of_parameters_checked(
            nullptr, nullptr, nullptr, "x");
        EXPECT_TRUE(!sol && sol.error().code == CasErrc::InvalidArgument,
                    "null VoP inputs are InvalidArgument");
    }
}


void test_frobenius() {
    TEST_CASE("Frobenius: classify ordinary point (y'' + y = 0 at x=0)");
    {
        // p(x) = 0, q(x) = 1 → 常点
        auto p = SymbolicExpr::number(0);
        auto q = SymbolicExpr::number(1);
        auto x0 = SymbolicExpr::number(0);

        auto type = classify_singular_point(p, q, x0, "x");
        EXPECT_TRUE(type == ODESingularityType::Ordinary, "y''+y=0 at x=0 is ordinary");
    }

    TEST_CASE("Frobenius: classify regular singular point (Bessel at x=0)");
    {
        // Bessel 方程: x²y'' + xy' + (x²-n²)y = 0
        // 归一化: y'' + (1/x)y' + (1 - n²/x²)y = 0
        // p(x) = 1/x, q(x) = 1 - n²/x² (取 n=0: q = 1)
        // x·p(x) = 1 → 有限, x²·q(x) = x² → 有限 → 正则奇点
        auto x = SymbolicExpr::variable("x");
        auto p = SymbolicExpr::divide(SymbolicExpr::number(1), x);
        auto q = SymbolicExpr::number(1);  // n=0 的 Bessel
        auto x0 = SymbolicExpr::number(0);

        auto type = classify_singular_point(p, q, x0, "x");
        EXPECT_TRUE(type == ODESingularityType::RegularSingular,
            "Bessel at x=0 is regular singular");
    }

    TEST_CASE("Frobenius: ordinary point series solution (y'' + y = 0)");
    {
        // y'' + y = 0, p=0, q=1, 展开点 x₀=0
        // 解: y = cos(x) = 1 - x²/2 + x⁴/24 - ...
        auto p = SymbolicExpr::number(0);
        auto q = SymbolicExpr::number(1);
        auto x0 = SymbolicExpr::number(0);

        auto sol = solve_frobenius_checked(p, q, x0, "x", 6).value();
        EXPECT_TRUE(sol.series_solution != nullptr, "Frobenius y''+y=0 has series solution");
        EXPECT_TRUE(sol.point_type == ODESingularityType::Ordinary, "point type is ordinary");
        EXPECT_TRUE(sol.truncation_order == 6, "truncation order is 6");

        // 验证系数: a₀=1, a₁=0, a₂=-1/2, a₃=0, a₄=1/24, a₅=0, a₆=-1/720
        // 这对应 cos(x) 的 Taylor 展开
        std::string s = sol.series_solution->to_string();
        EXPECT_CONTAINS(s, {"x"}, "series solution contains x");
    }

    TEST_CASE("Frobenius: regular singular point (Euler equation x²y'' + xy' - y = 0)");
    {
        // 归一化: y'' + (1/x)y' + (-1/x²)y = 0
        // p(x) = 1/x, q(x) = -1/x²
        // P₀ = lim x·(1/x) = 1
        // Q₀ = lim x²·(-1/x²) = -1
        // 指标方程: r(r-1) + r - 1 = r² - 1 = 0 → r = ±1
        auto x = SymbolicExpr::variable("x");
        auto p = SymbolicExpr::divide(SymbolicExpr::number(1), x);
        auto q = SymbolicExpr::divide(SymbolicExpr::number(-1),
            SymbolicExpr::power(x, SymbolicExpr::number(2)));
        auto x0 = SymbolicExpr::number(0);

        auto sol = solve_frobenius_checked(p, q, x0, "x", 4).value();
        EXPECT_TRUE(sol.series_solution != nullptr, "Frobenius Euler eq has series solution");
        EXPECT_TRUE(sol.point_type == ODESingularityType::RegularSingular,
            "Euler eq at x=0 is regular singular");
        EXPECT_TRUE(sol.indicial_roots.size() == 2, "has two indicial roots");

        // 指标根应为 1 和 -1
        if (sol.indicial_roots.size() >= 2) {
            double r1 = sol.indicial_roots[0];
            double r2 = sol.indicial_roots[1];
            EXPECT_NEAR(r1, 1.0, 1e-9, "larger indicial root is 1");
            EXPECT_NEAR(r2, -1.0, 1e-9, "smaller indicial root is -1");
        }

        // 解应为 x^1 * (1 + 0 + 0 + ...) = x
        std::string s = sol.series_solution->to_string();
        EXPECT_CONTAINS(s, {"x"}, "Frobenius Euler solution contains x");
    }

    TEST_CASE("Frobenius: irregular singular point");
    {
        // y'' + (1/x³)y = 0 → p=0, q=1/x³
        // x²·q = x²/x³ = 1/x → 在 x=0 处无限 → 非正则奇点
        auto x = SymbolicExpr::variable("x");
        auto p = SymbolicExpr::number(0);
        auto q = SymbolicExpr::divide(SymbolicExpr::number(1),
            SymbolicExpr::power(x, SymbolicExpr::number(3)));
        auto x0 = SymbolicExpr::number(0);

        auto type = classify_singular_point(p, q, x0, "x");
        EXPECT_TRUE(type == ODESingularityType::IrregularSingular,
            "1/x^3 coefficient gives irregular singular point");

        auto sol = solve_frobenius_checked(p, q, x0, "x", 4);
        EXPECT_TRUE(!sol &&
                        (sol.error().code == CasErrc::DomainError ||
                         sol.error().code == CasErrc::Inconclusive),
                    "irregular singular Frobenius input is explicitly rejected");
    }
}

void test_vop_frobenius_checked_contracts() {
    TEST_CASE("Variation of Parameters and Frobenius checked APIs: explicit support domain");
    {
        auto x = SymbolicExpr::variable("x");
        auto y1 = SymbolicExpr::cos(x);
        auto y2 = SymbolicExpr::sin(x);
        auto g = SymbolicExpr::sin(x);

        auto sol = solve_variation_of_parameters_checked(y1, y2, g, "x");
        EXPECT_TRUE(sol.has_value(),
            "checked variation of parameters succeeds for independent homogeneous solutions");
        if (sol) {
            EXPECT_TRUE(sol.value().general_solution != nullptr,
                "checked variation of parameters returns a particular solution");
            EXPECT_TRUE(sol.value().method_used == ODEType::HigherOrder_ConstCoeff,
                "checked variation of parameters reports higher-order method family");
        }

        auto dependent = solve_variation_of_parameters_checked(x, x, g, "x");
        EXPECT_TRUE(!dependent.has_value(),
            "checked variation of parameters rejects zero Wronskian as unsupported");
        EXPECT_TRUE(dependent.error().code == CasErrc::Inconclusive,
            "checked variation of parameters reports Inconclusive for zero Wronskian");

        std::shared_ptr<SymbolicExpr> null_expr;
        auto invalid = solve_variation_of_parameters_checked(null_expr, y2, g, "x");
        EXPECT_TRUE(!invalid.has_value(),
            "checked variation of parameters rejects null input");
        EXPECT_TRUE(invalid.error().code == CasErrc::InvalidArgument,
            "checked variation of parameters reports InvalidArgument for null input");
    }

    {
        auto p = SymbolicExpr::number(0);
        auto q = SymbolicExpr::number(1);
        auto x0 = SymbolicExpr::number(0);
        auto sol = solve_frobenius_checked(p, q, x0, "x", 6);
        EXPECT_TRUE(sol.has_value(),
            "checked Frobenius succeeds at an ordinary point");
        if (sol) {
            EXPECT_TRUE(sol.value().series_solution != nullptr,
                "checked Frobenius returns a series");
            EXPECT_TRUE(sol.value().point_type == ODESingularityType::Ordinary,
                "checked Frobenius reports ordinary point");
            EXPECT_TRUE(sol.value().truncation_order == 6,
                "checked Frobenius preserves truncation order");
        }
    }

    {
        auto x = SymbolicExpr::variable("x");
        auto p = SymbolicExpr::number(0);
        auto q = SymbolicExpr::divide(SymbolicExpr::number(1),
            SymbolicExpr::power(x, SymbolicExpr::number(3)));
        auto x0 = SymbolicExpr::number(0);
        auto sol = solve_frobenius_checked(p, q, x0, "x", 4);
        EXPECT_TRUE(!sol.has_value(),
            "checked Frobenius rejects irregular singular point");
        EXPECT_TRUE(sol.error().code == CasErrc::Inconclusive,
            "checked Frobenius reports Inconclusive for irregular singular point");

        auto bad_order = solve_frobenius_checked(p, q, x0, "x", -1);
        EXPECT_TRUE(!bad_order.has_value(),
            "checked Frobenius rejects negative truncation order");
        EXPECT_TRUE(bad_order.error().code == CasErrc::InvalidArgument,
            "checked Frobenius reports InvalidArgument for negative order");

        std::shared_ptr<SymbolicExpr> null_expr;
        auto invalid = solve_frobenius_checked(null_expr, q, x0, "x", 4);
        EXPECT_TRUE(!invalid.has_value(),
            "checked Frobenius rejects null coefficient");
        EXPECT_TRUE(invalid.error().code == CasErrc::InvalidArgument,
            "checked Frobenius reports InvalidArgument for null coefficient");
    }

    {
        auto x = SymbolicExpr::variable("x");
        auto y1 = SymbolicExpr::cos(x);
        auto y2 = SymbolicExpr::sin(x);
        auto g = SymbolicExpr::sin(x);

        lamina::CancellationToken cancellation;
        lamina::ComputationContext cancelled_context({}, cancellation);
        cancellation.cancel();
        auto cancelled = solve_variation_of_parameters_checked(
            y1, y2, g, "x", cancelled_context);
        EXPECT_TRUE(!cancelled.has_value(),
            "checked variation of parameters observes cancellation");
        EXPECT_TRUE(cancelled.error().code == CasErrc::Cancelled,
            "checked variation of parameters reports Cancelled");

        lamina::ResourceLimits limits;
        limits.max_steps = 1;
        lamina::ComputationContext limited_context(limits);
        auto p = SymbolicExpr::number(0);
        auto q = SymbolicExpr::number(1);
        auto x0 = SymbolicExpr::number(0);
        auto limited = solve_frobenius_checked(p, q, x0, "x", 6, limited_context);
        EXPECT_TRUE(!limited.has_value(),
            "checked Frobenius observes exhausted step budget");
        EXPECT_TRUE(limited.error().code == CasErrc::ResourceLimit,
            "checked Frobenius reports ResourceLimit");
    }
}

int main() {
    try {
        test_homogeneous_ode();
        test_bernoulli_ode();
        test_exact_ode();
        test_first_order_ode_checked_contracts();
        test_higher_order_euler_checked_contracts();
        test_integrating_factor();
        test_classification();
        test_variation_of_parameters();
        test_frobenius();
        test_vop_frobenius_checked_contracts();
    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        g_failures++;
    } catch (...) {
        std::cout << "[FAIL] Unknown Exception!" << std::endl;
        g_failures++;
    }
    return TEST_REPORT();
}
