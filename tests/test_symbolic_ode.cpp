#include "test_common.hpp"
#include "symbolic_ode.hpp"

void test_separable_ode() {
    TEST_CASE("Separable ODE: dy/dx = x/y");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        // rhs = x / y
        auto rhs = SymbolicExpr::divide(x, y);
        auto sol = lamina::solve_separable_ode(rhs, "x", "y");
        std::string s = sol ? sol->to_string() : "null";
        // Solution involves y^2 term (from integrating y dy)
        EXPECT_CONTAINS(s, {"y"}, "separable x/y contains y");
        EXPECT_CONTAINS(s, {"x"}, "separable x/y contains x");
    }

    TEST_CASE("Separable ODE: dy/dx = x*y");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        // rhs = x * y
        auto rhs = SymbolicExpr::multiply(x, y);
        auto sol = lamina::solve_separable_ode(rhs, "x", "y");
        std::string s = sol ? sol->to_string() : "null";
        // Solution involves ln(y) (from integrating 1/y dy)
        EXPECT_CONTAINS(s, {"ln(y)"}, "separable x*y contains ln(y)");
        EXPECT_CONTAINS(s, {"x"}, "separable x*y contains x");
    }

    TEST_CASE("Separable ODE: dy/dx = y/x");
    {
        auto x = SymbolicExpr::variable("x");
        auto y = SymbolicExpr::variable("y");
        // rhs = y / x
        auto rhs = SymbolicExpr::divide(y, x);
        auto sol = lamina::solve_separable_ode(rhs, "x", "y");
        std::string s = sol ? sol->to_string() : "null";
        // Solution involves ln(y) (from integrating 1/y dy)
        EXPECT_CONTAINS(s, {"ln(y)"}, "separable y/x contains ln(y)");
        EXPECT_CONTAINS(s, {"x"}, "separable y/x contains x");
    }
}

void test_linear1_ode() {
    TEST_CASE("Linear 1st Order ODE: dy/dx + 2*y = 0 (P=2, Q=0)");
    {
        // dy/dx + 2*y = 0 => P(x) = 2, Q(x) = 0
        auto Px = SymbolicExpr::number(2);
        auto Qx = SymbolicExpr::number(0);
        auto sol = lamina::solve_linear1_ode(Px, Qx, "x", "y");
        std::string s = sol ? sol->to_string() : "null";
        EXPECT_CONTAINS(s, {"C"}, "linear1 P=2 Q=0 contains integration constant C");
    }

    TEST_CASE("Linear 1st Order ODE: dy/dx + x*y = x (P=x, Q=x)");
    {
        // dy/dx + x*y = x => P(x) = x, Q(x) = x
        auto x = SymbolicExpr::variable("x");
        auto Px = x;
        auto Qx = SymbolicExpr::variable("x");
        auto sol = lamina::solve_linear1_ode(Px, Qx, "x", "y");
        std::string s = sol ? sol->to_string() : "null";
        EXPECT_CONTAINS(s, {"C"}, "linear1 P=x Q=x contains integration constant C");
    }

    TEST_CASE("Linear 1st Order ODE: dy/dx + (1/x)*y = x^2 (P=1/x, Q=x^2)");
    {
        // dy/dx + (1/x)*y = x^2 => P(x) = 1/x, Q(x) = x^2
        auto x = SymbolicExpr::variable("x");
        auto one = SymbolicExpr::number(1);
        auto Px = SymbolicExpr::divide(one, x);
        auto Qx = SymbolicExpr::power(SymbolicExpr::variable("x"), SymbolicExpr::number(2));
        auto sol = lamina::solve_linear1_ode(Px, Qx, "x", "y");
        std::string s = sol ? sol->to_string() : "null";
        EXPECT_CONTAINS(s, {"C"}, "linear1 P=1/x Q=x^2 contains integration constant C");
    }
}

void test_linear2_ode() {
    TEST_CASE("Linear 2nd Order ODE: y'' - 3y' + 2y = 0 (distinct real roots r=1,2)");
    {
        // a=1, b=-3, c=2, f(x)=0
        // Characteristic: r^2 - 3r + 2 = 0 => r=1, r=2
        auto fx = SymbolicExpr::number(0);
        auto sol = lamina::solve_linear2_ode(1, -3, 2, fx, "x", "y");
        std::string s = sol ? sol->to_string() : "null";
        // Expect exponential terms e^x and e^(2x) with constants
        EXPECT_CONTAINS(s, {"C"}, "linear2 distinct roots contains constant C");
    }

    TEST_CASE("Linear 2nd Order ODE: y'' - 2y' + y = 0 (repeated root r=1)");
    {
        // a=1, b=-2, c=1, f(x)=0
        // Characteristic: r^2 - 2r + 1 = 0 => r=1 (double)
        auto fx = SymbolicExpr::number(0);
        auto sol = lamina::solve_linear2_ode(1, -2, 1, fx, "x", "y");
        std::string s = sol ? sol->to_string() : "null";
        // Repeated root: solution is (C1 + C2*x)*e^x
        EXPECT_CONTAINS(s, {"C"}, "linear2 repeated root contains constant C");
    }

    TEST_CASE("Linear 2nd Order ODE: y'' + y = 0 (complex roots r=±i)");
    {
        // a=1, b=0, c=1, f(x)=0
        // Characteristic: r^2 + 1 = 0 => r=+/-i
        auto fx = SymbolicExpr::number(0);
        auto sol = lamina::solve_linear2_ode(1, 0, 1, fx, "x", "y");
        std::string s = sol ? sol->to_string() : "null";
        // Complex roots: solution involves sin and cos
        EXPECT_CONTAINS(s, {"C"}, "linear2 complex roots contains constant C");
    }
}

void test_linear2_nonhomogeneous() {
    TEST_CASE("Linear 2nd Order ODE Non-homogeneous: y'' - 3y' + 2y = e^(3x)");
    {
        // a=1, b=-3, c=2, f(x) = e^(3x)
        /// 旧版接口对支持域之外的非齐次输入抛出 std::logic_error.
        auto x = SymbolicExpr::variable("x");
        auto three_x = SymbolicExpr::multiply(SymbolicExpr::number(3), x);
        auto fx = SymbolicExpr::exp(three_x);
        bool threw = false;
        try {
            auto sol = lamina::solve_linear2_ode(1, -3, 2, fx, "x", "y");
            // If it doesn't throw, verify the solution contains expected structure
            std::string s = sol ? sol->to_string() : "null";
            EXPECT_CONTAINS(s, {"C"}, "linear2 nonhomogeneous contains constant C");
        } catch (const std::logic_error&) {
            // Non-homogeneous case is not yet implemented - this is expected
            threw = true;
        }
        EXPECT_TRUE(threw, "legacy linear2 nonhomogeneous throws instead of returning a false solution");
    }
}

void test_linear2_checked_contracts() {
    TEST_CASE("Checked Linear 2nd Order ODE: homogeneous success and explicit failures");
    {
        auto fx = SymbolicExpr::number(0);
        auto result = lamina::solve_linear2_ode_checked(1, -3, 2, fx, "x", "y");
        EXPECT_TRUE(result.has_value(), "checked linear2 homogeneous succeeds");
        if (result) {
            EXPECT_TRUE(result.value() != nullptr, "checked linear2 homogeneous returns expression");
            EXPECT_CONTAINS(result.value()->to_string(), {"C"},
                            "checked linear2 homogeneous contains integration constants");
        }
    }

    {
        auto x = SymbolicExpr::variable("x");
        auto three_x = SymbolicExpr::multiply(SymbolicExpr::number(3), x);
        auto fx = SymbolicExpr::exp(three_x);
        auto result = lamina::solve_linear2_ode_checked(1, -3, 2, fx, "x", "y");
        EXPECT_TRUE(!result.has_value(), "checked linear2 nonhomogeneous is not a success");
        EXPECT_TRUE(result.error().code == lamina::CasErrc::Inconclusive,
                    "checked linear2 nonhomogeneous reports Inconclusive");
    }

    {
        auto result = lamina::solve_linear2_ode_checked(1, -3, 2, nullptr, "x", "y");
        EXPECT_TRUE(!result.has_value(), "checked linear2 rejects null forcing expression");
        EXPECT_TRUE(result.error().code == lamina::CasErrc::InvalidArgument,
                    "checked linear2 null forcing reports InvalidArgument");
    }

    {
        auto fx = SymbolicExpr::number(0);
        auto result = lamina::solve_linear2_ode_checked(1, -3, 2, fx, "", "y");
        EXPECT_TRUE(!result.has_value(), "checked linear2 rejects empty independent variable");
        EXPECT_TRUE(result.error().code == lamina::CasErrc::InvalidArgument,
                    "checked linear2 empty variable reports InvalidArgument");
    }

    {
        lamina::CancellationToken cancellation;
        lamina::ComputationContext context({}, cancellation);
        cancellation.cancel();
        auto result = lamina::solve_linear2_ode_checked(
            1, -3, 2, SymbolicExpr::number(0), "x", "y", context);
        EXPECT_TRUE(!result.has_value(), "checked linear2 observes cancellation");
        EXPECT_TRUE(result.error().code == lamina::CasErrc::Cancelled,
                    "checked linear2 cancellation reports Cancelled");
    }

    {
        lamina::ResourceLimits limits;
        limits.max_steps = 0;
        lamina::ComputationContext context(limits);
        auto result = lamina::solve_linear2_ode_checked(
            1, -3, 2, SymbolicExpr::number(0), "x", "y", context);
        EXPECT_TRUE(!result.has_value(), "checked linear2 observes exhausted step budget");
        EXPECT_TRUE(result.error().code == lamina::CasErrc::ResourceLimit,
                    "checked linear2 exhausted budget reports ResourceLimit");
    }
}

int main() {
    try {
        test_separable_ode();
        test_linear1_ode();
        test_linear2_ode();
        test_linear2_nonhomogeneous();
        test_linear2_checked_contracts();
    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        g_failures++;
    } catch (...) {
        std::cout << "[FAIL] Unknown Exception!" << std::endl;
        g_failures++;
    }
    return TEST_REPORT();
}
