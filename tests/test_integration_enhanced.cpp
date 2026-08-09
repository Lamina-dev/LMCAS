/**
 * @file test_integration_enhanced.cpp
 * @brief Comprehensive integration test suite for enhanced integration features.
 *
 * Verifies that a returned antiderivative differentiates to the integrand.
 */
#include "test_common.hpp"
#include "integration.hpp"
#include "symbolic.hpp"

int main() {
    TEST_CASE("Integration enhanced - power round trip");
    auto x = SymbolicExpr::variable("x");
    auto integrand = SymbolicExpr::power(x, SymbolicExpr::number(2));
    lamina::Integrator integrator;
    auto antiderivative = integrator.integrate(*integrand, "x");
    auto derivative = antiderivative.differentiate("x");
    EXPECT_TRUE(derivative != nullptr, "antiderivative is differentiable");
    if (derivative) {
        auto delta = test_normalized_delta(derivative, integrand);
        EXPECT_TRUE(delta && delta->is_zero(), "d/dx integral(x^2) equals x^2");
    }
    return TEST_REPORT();
}
