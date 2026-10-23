#include "test_common.hpp"
#include "integration.hpp"
#include <string>

using namespace LMCAS;

TEST(LmcasIntegrationFallback, IndependentVariablePrimitive) {
    auto y = SymbolicExpr::variable("y");
    Integrator integrator;
    ComputationContext context;
    auto res = integrator.integrate_checked(*y, "x", context);
    ASSERT_TRUE(res) << res.error().message;
    std::string s = res.value().to_string();
    bool has_xy = (s.find("y*x") != std::string::npos) || (s.find("x*y") != std::string::npos);
    bool has_unevaluated = (s.find("integral(") != std::string::npos) || (s.find("Integral(") != std::string::npos);
    EXPECT_TRUE((has_xy)) << ("Integration of y w.r.t x should produce a product of x and y");
    EXPECT_TRUE((!has_unevaluated)) << ("Result should not be an unevaluated integral");
}
