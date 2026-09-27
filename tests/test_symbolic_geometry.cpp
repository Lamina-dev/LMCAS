#include "test_common.hpp"
#include "symbolic_geometry.hpp"

using namespace LMCAS;

TEST(SymbolicGeometry, VolumeRevolutionX) {
    {
        // Cone: V = pi * integral(x^2, 0, h) = pi * h^3/3
        auto x = SymbolicExpr::variable("x");
        auto h = SymbolicExpr::variable("h");
        auto zero = SymbolicExpr::number(0);
        auto result = LMCAS::volume_of_revolution_x_checked(x, zero, h);
        ASSERT_TRUE(result.has_value());
        ASSERT_NE(result.value(), nullptr);
        auto expected = SymbolicExpr::divide(
            SymbolicExpr::multiply(
                SymbolicExpr::variable("pi"),
                SymbolicExpr::power(h, SymbolicExpr::number(3))),
            SymbolicExpr::number(3));
        EXPECT_TRUE(test_proved_equivalent(result.value(), expected));
        EXPECT_FALSE(LMCAS::detail::contains_node_type<IntegralNode>(
            LMCAS::detail::node(*result.value())));
    }

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
        auto result = LMCAS::volume_of_revolution_x_checked(fx, neg_r, r);
        ASSERT_TRUE(result.has_value());
        ASSERT_NE(result.value(), nullptr);
        auto expected = SymbolicExpr::divide(
            SymbolicExpr::multiply(
                SymbolicExpr::multiply(
                    SymbolicExpr::number(4), SymbolicExpr::variable("pi")),
                SymbolicExpr::power(r, SymbolicExpr::number(3))),
            SymbolicExpr::number(3));
        EXPECT_TRUE(test_proved_equivalent(result.value(), expected));
        EXPECT_FALSE(LMCAS::detail::contains_node_type<IntegralNode>(
            LMCAS::detail::node(*result.value())));
    }
}

TEST(SymbolicGeometry, ArcLengthX) {
    auto x = SymbolicExpr::variable("x");
    auto linear = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::number(1));
    auto zero = SymbolicExpr::number(0);
    auto three = SymbolicExpr::number(3);
    auto result = LMCAS::arc_length_x_checked(linear, zero, three);
    ASSERT_TRUE(result.has_value()) << result.error().message;
    ASSERT_NE(result.value(), nullptr);
    auto expected = SymbolicExpr::multiply(
        three, SymbolicExpr::sqrt(SymbolicExpr::number(5)));
    EXPECT_TRUE(test_proved_equivalent(result.value(), expected));

    auto quadratic = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto unsupported = LMCAS::arc_length_x_checked(
        quadratic, zero, SymbolicExpr::number(1));
    ASSERT_FALSE(unsupported.has_value());
    EXPECT_EQ(unsupported.error().code, LMCAS::CasErrc::Inconclusive);
}

TEST(SymbolicGeometry, VolumeRevolutionY) {
    {
        // Cone about y-axis: V = pi * integral(y^2, 0, h) = pi * h^3/3
        auto y = SymbolicExpr::variable("y");
        auto h = SymbolicExpr::variable("h");
        auto zero = SymbolicExpr::number(0);
        auto result = LMCAS::volume_of_revolution_y_checked(y, zero, h);
        ASSERT_TRUE(result.has_value());
        ASSERT_NE(result.value(), nullptr);
        auto expected = SymbolicExpr::divide(
            SymbolicExpr::multiply(
                SymbolicExpr::variable("pi"),
                SymbolicExpr::power(h, SymbolicExpr::number(3))),
            SymbolicExpr::number(3));
        EXPECT_TRUE(test_proved_equivalent(result.value(), expected));
        EXPECT_FALSE(LMCAS::detail::contains_node_type<IntegralNode>(
            LMCAS::detail::node(*result.value())));
    }
}

TEST(SymbolicGeometry, ArcLengthY) {
    auto y = SymbolicExpr::variable("y");
    auto linear = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), y),
        SymbolicExpr::number(1));
    auto zero = SymbolicExpr::number(0);
    auto three = SymbolicExpr::number(3);
    auto result = LMCAS::arc_length_y_checked(linear, zero, three);
    ASSERT_TRUE(result.has_value()) << result.error().message;
    ASSERT_NE(result.value(), nullptr);
    auto expected = SymbolicExpr::multiply(
        three, SymbolicExpr::sqrt(SymbolicExpr::number(5)));
    EXPECT_TRUE(test_proved_equivalent(result.value(), expected));
}

TEST(SymbolicGeometry, SymbolicGeometryCheckedContracts) {
    {
        auto x = SymbolicExpr::variable("x");
        auto zero = SymbolicExpr::number(0);
        auto one = SymbolicExpr::number(1);

        auto volume = LMCAS::volume_of_revolution_x_checked(x, zero, one);
        ASSERT_TRUE((volume.has_value())) << "checked volume_of_revolution_x succeeds";
        if (volume) {
            EXPECT_TRUE((volume.value() != nullptr)) << "checked volume_of_revolution_x returns an expression";
        }

        auto null_profile = LMCAS::volume_of_revolution_x_checked(nullptr, zero, one);
        EXPECT_TRUE((!null_profile.has_value())) << "checked volume_of_revolution_x rejects null profile";
        EXPECT_TRUE((null_profile.error().code == LMCAS::CasErrc::InvalidArgument)) << "checked volume_of_revolution_x reports InvalidArgument";

        std::shared_ptr<SymbolicExpr> null_root;
        auto null_bound = LMCAS::arc_length_x_checked(x, null_root, one);
        EXPECT_TRUE((!null_bound.has_value())) << "checked arc_length_x rejects null bound";
        EXPECT_TRUE((null_bound.error().code == LMCAS::CasErrc::InvalidArgument)) << "checked arc_length_x reports InvalidArgument for null bound";

        auto unsupported_profile = SymbolicExpr::eq(x, zero);
        auto unsupported_arc = LMCAS::arc_length_x_checked(
            unsupported_profile, zero, one);
        EXPECT_TRUE((!unsupported_arc.has_value())) << "checked arc_length_x rejects unsupported derivatives";
        EXPECT_TRUE((unsupported_arc.error().code == LMCAS::CasErrc::Inconclusive)) << "checked arc_length_x reports Inconclusive for unsupported derivatives";

        LMCAS::CancellationToken cancellation;
        LMCAS::ComputationContext cancelled_context({}, cancellation);
        cancellation.cancel();
        auto cancelled = LMCAS::volume_of_revolution_y_checked(
            x, zero, one, cancelled_context);
        EXPECT_TRUE((!cancelled.has_value())) << "checked volume_of_revolution_y observes cancellation";
        EXPECT_TRUE((cancelled.error().code == LMCAS::CasErrc::Cancelled)) << "checked volume_of_revolution_y reports Cancelled";

        LMCAS::ResourceLimits limits;
        limits.max_steps = 0;
        LMCAS::ComputationContext limited_context(limits);
        auto limited = LMCAS::arc_length_y_checked(x, zero, one, limited_context);
        EXPECT_TRUE((!limited.has_value())) << "checked arc_length_y observes exhausted step budget";
        EXPECT_TRUE((limited.error().code == LMCAS::CasErrc::ResourceLimit)) << "checked arc_length_y reports ResourceLimit";
    }
}
