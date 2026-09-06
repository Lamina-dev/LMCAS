#include "test_common.hpp"
#include "symbolic_geometry.hpp"

using namespace LMCAS;

void test_volume_revolution_x() {
    TEST_CASE("Volume of Revolution X: f(x) = x on [0, h] (cone)");
    {
        // Cone: V = pi * integral(x^2, 0, h) = pi * h^3/3
        auto x = SymbolicExpr::variable("x");
        auto h = SymbolicExpr::variable("h");
        auto zero = SymbolicExpr::number(0);
        auto result = LMCAS::volume_of_revolution_x_checked(x, zero, h).value();
        std::string s = result ? result->to_string() : "null";
        std::cout << "  Cone volume result: " << s << std::endl;
        // The result should keep the exact pi/3 factor.
        EXPECT_CONTAINS(s, {"h"}, "cone volume contains h");
        EXPECT_CONTAINS(s, {"1/3", "pi"}, "cone volume contains exact pi/3 factor");
    }

    TEST_CASE("Volume of Revolution X: f(x) = sqrt(r^2 - x^2) on [-r, r] (sphere)");
    {
        // Sphere: V = pi * integral(r^2 - x^2, -r, r) = 4*pi*r^3/3
        auto x = SymbolicExpr::variable("x");
        auto r = SymbolicExpr::variable("r");
        auto neg_r = SymbolicExpr::multiply(SymbolicExpr::number(-1), r);
        // f(x) = sqrt(r^2 - x^2)
        auto r_sq = SymbolicExpr::power(r, SymbolicExpr::number(2));
        auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto inner = SymbolicExpr::add(r_sq, SymbolicExpr::multiply(SymbolicExpr::number(-1), x_sq));
        auto fx = SymbolicExpr::sqrt(inner);
        auto result = LMCAS::volume_of_revolution_x_checked(fx, neg_r, r).value();
        std::string s = result ? result->to_string() : "null";
        std::cout << "  Sphere volume result: " << s << std::endl;
        // The result should stay exact and must not leave the polynomial integral unevaluated.
        EXPECT_CONTAINS(s, {"r"}, "sphere volume contains r");
        EXPECT_CONTAINS(s, {"pi"}, "sphere volume contains exact pi factor");
        EXPECT_TRUE(s.find("integral(") == std::string::npos,
                    "sphere volume polynomial integral is evaluated");
    }
}

void test_arc_length_x() {
    TEST_CASE("Arc Length X: f(x) = 2x + 1 on [0, 3] (linear function)");
    {
        // For a linear function f(x) = 2x + 1, f'(x) = 2
        // Arc length = integral(sqrt(1 + 4), 0, 3) = 3*sqrt(5)
        auto x = SymbolicExpr::variable("x");
        auto fx = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::number(1)
        );
        auto zero = SymbolicExpr::number(0);
        auto three = SymbolicExpr::number(3);
        auto result = LMCAS::arc_length_x_checked(fx, zero, three).value();
        std::string s = result ? result->to_string() : "null";
        std::cout << "  Linear arc length result: " << s << std::endl;
        // For linear function, the integrand sqrt(1+4)=sqrt(5) is constant
        // Result should be 3*sqrt(5) ≈ 6.7082
        EXPECT_TRUE(result != nullptr, "linear arc length is not null");
    }

    TEST_CASE("Arc Length X: f(x) = x^2 on [0, 1]");
    {
        // f(x) = x^2, f'(x) = 2x
        // Arc length = integral(sqrt(1 + 4x^2), 0, 1), outside the
        // checked symbolic geometry support domain.
        auto x = SymbolicExpr::variable("x");
        auto fx = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto zero = SymbolicExpr::number(0);
        auto one = SymbolicExpr::number(1);

        auto checked = LMCAS::arc_length_x_checked(fx, zero, one);
        EXPECT_TRUE(!checked.has_value(),
                    "checked x^2 arc length rejects unsupported integral");
        EXPECT_TRUE(checked.error().code == LMCAS::CasErrc::Inconclusive,
                    "checked x^2 arc length reports Inconclusive");
    }
}

void test_volume_revolution_y() {
    TEST_CASE("Volume of Revolution Y: f(y) = y on [0, h] (cone about y-axis)");
    {
        // Cone about y-axis: V = pi * integral(y^2, 0, h) = pi * h^3/3
        auto y = SymbolicExpr::variable("y");
        auto h = SymbolicExpr::variable("h");
        auto zero = SymbolicExpr::number(0);
        auto result = LMCAS::volume_of_revolution_y_checked(y, zero, h).value();
        std::string s = result ? result->to_string() : "null";
        std::cout << "  Y-axis cone volume result: " << s << std::endl;
        // The result should keep the exact pi/3 factor.
        EXPECT_CONTAINS(s, {"h"}, "y-axis cone volume contains h");
        EXPECT_CONTAINS(s, {"1/3", "pi"}, "y-axis cone volume contains exact pi/3 factor");
    }
}

void test_arc_length_y() {
    TEST_CASE("Arc Length Y: f(y) = 2y + 1 on [0, 3] (linear function)");
    {
        // For a linear function f(y) = 2y + 1, f'(y) = 2
        // Arc length = integral(sqrt(1 + 4), 0, 3) = 3*sqrt(5)
        auto y = SymbolicExpr::variable("y");
        auto fy = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), y),
            SymbolicExpr::number(1)
        );
        auto zero = SymbolicExpr::number(0);
        auto three = SymbolicExpr::number(3);
        auto result = LMCAS::arc_length_y_checked(fy, zero, three).value();
        std::string s = result ? result->to_string() : "null";
        std::cout << "  Y-axis linear arc length result: " << s << std::endl;
        EXPECT_TRUE(result != nullptr, "y-axis linear arc length is not null");
    }
}

void test_symbolic_geometry_checked_contracts() {
    TEST_CASE("Symbolic Geometry checked APIs: explicit errors and cancellation");
    {
        auto x = SymbolicExpr::variable("x");
        auto zero = SymbolicExpr::number(0);
        auto one = SymbolicExpr::number(1);

        auto volume = LMCAS::volume_of_revolution_x_checked(x, zero, one);
        EXPECT_TRUE(volume.has_value(), "checked volume_of_revolution_x succeeds");
        if (volume) {
            EXPECT_TRUE(volume.value() != nullptr,
                        "checked volume_of_revolution_x returns an expression");
        }

        auto null_profile = LMCAS::volume_of_revolution_x_checked(nullptr, zero, one);
        EXPECT_TRUE(!null_profile.has_value(),
                    "checked volume_of_revolution_x rejects null profile");
        EXPECT_TRUE(null_profile.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked volume_of_revolution_x reports InvalidArgument");

        std::shared_ptr<SymbolicExpr> null_root;
        auto null_bound = LMCAS::arc_length_x_checked(x, null_root, one);
        EXPECT_TRUE(!null_bound.has_value(),
                    "checked arc_length_x rejects null bound");
        EXPECT_TRUE(null_bound.error().code == LMCAS::CasErrc::InvalidArgument,
                    "checked arc_length_x reports InvalidArgument for null bound");

        auto unsupported_profile = SymbolicExpr::eq(x, zero);
        auto unsupported_arc = LMCAS::arc_length_x_checked(
            unsupported_profile, zero, one);
        EXPECT_TRUE(!unsupported_arc.has_value(),
                    "checked arc_length_x rejects unsupported derivatives");
        EXPECT_TRUE(unsupported_arc.error().code == LMCAS::CasErrc::Inconclusive,
                    "checked arc_length_x reports Inconclusive for unsupported derivatives");

        LMCAS::CancellationToken cancellation;
        LMCAS::ComputationContext cancelled_context({}, cancellation);
        cancellation.cancel();
        auto cancelled = LMCAS::volume_of_revolution_y_checked(
            x, zero, one, cancelled_context);
        EXPECT_TRUE(!cancelled.has_value(),
                    "checked volume_of_revolution_y observes cancellation");
        EXPECT_TRUE(cancelled.error().code == LMCAS::CasErrc::Cancelled,
                    "checked volume_of_revolution_y reports Cancelled");

        LMCAS::ResourceLimits limits;
        limits.max_steps = 0;
        LMCAS::ComputationContext limited_context(limits);
        auto limited = LMCAS::arc_length_y_checked(x, zero, one, limited_context);
        EXPECT_TRUE(!limited.has_value(),
                    "checked arc_length_y observes exhausted step budget");
        EXPECT_TRUE(limited.error().code == LMCAS::CasErrc::ResourceLimit,
                    "checked arc_length_y reports ResourceLimit");
    }
}

int main() {
    try {
        test_volume_revolution_x();
        test_arc_length_x();
        test_volume_revolution_y();
        test_arc_length_y();
        test_symbolic_geometry_checked_contracts();
    } catch (const std::exception& e) {
        std::cout << "[FAIL] Exception: " << e.what() << std::endl;
        g_failures++;
    } catch (...) {
        std::cout << "[FAIL] Unknown Exception!" << std::endl;
        g_failures++;
    }
    return TEST_REPORT();
}
