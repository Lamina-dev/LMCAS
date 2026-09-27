#include "test_common.hpp"
#include "expr.hpp"
#include "symbolic_ode_engine.hpp"
#include "poly_utils.hpp"
#include "numeric_evaluation.hpp"
#include <limits>

using namespace LMCAS;
static void expect_numeric_zero(
    const ExprPtr& expression,
    const std::vector<NumericBindings>& samples)
{
    ASSERT_TRUE(expression);
    for (const auto& bindings : samples) {
        auto specialized = expression;
        for (const auto& [name, value] : bindings) {
            specialized = specialized->substitute(
                name, SymbolicExpr::number(value));
        }
        auto numeric = test_numeric_eval(specialized->simplify());
        ASSERT_TRUE(numeric.has_value()) << specialized->to_string();
        EXPECT_NEAR(numeric.value(), 0.0, 1e-9)
            << specialized->to_string();
    }
}

static void expect_implicit_solution(
    const ODESolution& solution, const ExprPtr& rhs)
{
    ASSERT_TRUE(solution.general_solution);
    auto dx = solution.general_solution->differentiate("x");
    auto dy = solution.general_solution->differentiate("y");
    ASSERT_TRUE(dx && dy);
    auto residual = SymbolicExpr::add(
        dx, SymbolicExpr::multiply(dy, rhs))->simplify();
    expect_numeric_zero(
        residual,
        {{{"x", 0.5}, {"y", 0.75}},
         {{"x", 1.25}, {"y", 0.5}},
         {{"x", 2.0}, {"y", 1.5}}});
}

static void expect_bernoulli_solution(
    const ODESolution& solution, const ExprPtr& p,
    const ExprPtr& q, int exponent)
{
    ASSERT_TRUE(solution.general_solution);
    auto derivative = solution.general_solution->differentiate("x");
    ASSERT_TRUE(derivative);
    auto lhs = SymbolicExpr::add(
        derivative,
        SymbolicExpr::multiply(p, solution.general_solution));
    auto rhs = SymbolicExpr::multiply(
        q, SymbolicExpr::power(
               solution.general_solution, SymbolicExpr::number(exponent)));
    auto residual = SymbolicExpr::add(
        lhs, SymbolicExpr::multiply(SymbolicExpr::number(-1), rhs))
        ->simplify();
    expect_numeric_zero(
        residual,
        {{{"x", 0.5}, {"C", 2.0}},
         {{"x", 1.0}, {"C", 3.0}}});
}

static void expect_exact_potential(
    const ODESolution& solution, const ExprPtr& m, const ExprPtr& n)
{
    ASSERT_TRUE(solution.general_solution);
    auto dx = solution.general_solution->differentiate("x");
    auto dy = solution.general_solution->differentiate("y");
    ASSERT_TRUE(dx && dy);
    EXPECT_TRUE(test_proved_equivalent(dx->simplify(), m->simplify()));
    EXPECT_TRUE(test_proved_equivalent(dy->simplify(), n->simplify()));
}


TEST(OdeEngineFirstOrder, HomogeneousRatio) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto rhs = SymbolicExpr::divide(y, x);

    EXPECT_TRUE((is_homogeneous_ode(rhs, "x", "y"))) << "y/x is homogeneous";

    auto result = solve_homogeneous_ode_checked(rhs, "x", "y");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value().method_used, ODEType::Homogeneous);
    EXPECT_EQ(result.value().constants, std::vector<std::string>{"C"});
    expect_implicit_solution(result.value(), rhs);
}

TEST(OdeEngineFirstOrder, HomogeneousSum) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto rhs = SymbolicExpr::divide(
        SymbolicExpr::add(x, y), x);

    EXPECT_TRUE((is_homogeneous_ode(rhs, "x", "y"))) << "(x+y)/x is homogeneous";

    auto result = solve_homogeneous_ode_checked(rhs, "x", "y");
    ASSERT_TRUE(result.has_value());
    expect_implicit_solution(result.value(), rhs);
}

TEST(OdeEngineFirstOrder, HomogeneousQuadratic) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto y2 = SymbolicExpr::power(y, SymbolicExpr::number(2));
    auto rhs = SymbolicExpr::divide(
        SymbolicExpr::add(x2, y2),
        SymbolicExpr::multiply(x, y));

    EXPECT_TRUE((is_homogeneous_ode(rhs, "x", "y"))) << "(x^2+y^2)/(xy) is homogeneous";

    auto result = solve_homogeneous_ode_checked(rhs, "x", "y");
    ASSERT_TRUE(result.has_value());
    expect_implicit_solution(result.value(), rhs);
}

TEST(OdeEngineFirstOrder, BernoulliQuadratic) {
    auto P = SymbolicExpr::number(1);
    auto Q = SymbolicExpr::number(1);

    auto result = solve_bernoulli_ode_checked(P, Q, 2, "x", "y");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value().method_used, ODEType::Bernoulli);
    EXPECT_EQ(result.value().constants, std::vector<std::string>{"C"});
    expect_bernoulli_solution(result.value(), P, Q, 2);
}

TEST(OdeEngineFirstOrder, BernoulliVariable) {
    auto x = SymbolicExpr::variable("x");
    auto P = SymbolicExpr::divide(SymbolicExpr::number(1), x);
    auto Q = SymbolicExpr::variable("x");

    auto result = solve_bernoulli_ode_checked(P, Q, 2, "x", "y");
    ASSERT_TRUE(result.has_value());
    expect_bernoulli_solution(result.value(), P, Q, 2);
}

TEST(OdeEngineFirstOrder, BernoulliCubic) {
    auto P = SymbolicExpr::number(2);
    auto Q = SymbolicExpr::number(1);

    auto result = solve_bernoulli_ode_checked(P, Q, 3, "x", "y");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value().method_used, ODEType::Bernoulli);
    expect_bernoulli_solution(result.value(), P, Q, 3);
}

TEST(OdeEngineFirstOrder, ExactQuadratic) {
    /**
     * @brief M = 2x+y、N = x+2y 满足 ∂M/∂y = ∂N/∂x = 1，为恰当方程。
     * 隐式解为 F(x,y) = x^2+xy+y^2 = C。
     */
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto M = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x), y);
    auto N = SymbolicExpr::add(
        x, SymbolicExpr::multiply(SymbolicExpr::number(2), y));

    EXPECT_TRUE((is_exact_ode(M, N, "x", "y"))) << "(2x+y, x+2y) is exact";

    auto result = solve_exact_ode_checked(M, N, "x", "y");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value().method_used, ODEType::Exact);
    EXPECT_EQ(result.value().constants, std::vector<std::string>{"C"});
    expect_exact_potential(result.value(), M, N);
}

TEST(OdeEngineFirstOrder, ExactTranscendental) {
    /**
     * @brief M = y*cos(x)+2x*e^y、N = sin(x)+x^2*e^y 构成恰当方程。
     * ∂M/∂y = ∂N/∂x = cos(x)+2x*e^y。
     */
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

    EXPECT_TRUE((is_exact_ode(M, N, "x", "y"))) << "trig+exp exact equation";

    auto result = solve_exact_ode_checked(M, N, "x", "y");
    ASSERT_TRUE(result.has_value());
    expect_exact_potential(result.value(), M, N);
}

TEST(OdeEngineFirstOrder, ExactProduct) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    EXPECT_TRUE((is_exact_ode(y, x, "x", "y"))) << "(y, x) is exact";

    auto result = solve_exact_ode_checked(y, x, "x", "y");
    ASSERT_TRUE(result.has_value());
    expect_exact_potential(result.value(), y, x);
}

TEST(OdeEngineFirstOrder, ExactIdentityProof) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto factor = SymbolicExpr::multiply(
        SymbolicExpr::multiply(
            SymbolicExpr::add(x, SymbolicExpr::number(-1)),
            SymbolicExpr::add(x, SymbolicExpr::number(-2))),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(2), x),
            SymbolicExpr::number(-1)));
    auto M = SymbolicExpr::multiply(y, factor);
    auto N = SymbolicExpr::number(0);

    EXPECT_FALSE((is_exact_ode(M, N, "x", "y"))) << "vanishing at the former sample points is not an exactness proof";
}

TEST(OdeEngineFirstOrder, HomogeneousChecked) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto rhs = SymbolicExpr::divide(y, x);
    auto sol = solve_homogeneous_ode_checked(rhs, "x", "y");
    ASSERT_TRUE((sol.has_value())) << "checked homogeneous ODE succeeds";
    if (sol) {
        EXPECT_TRUE((sol.value().general_solution != nullptr)) << "checked homogeneous ODE returns solution";
        EXPECT_TRUE((sol.value().method_used == ODEType::Homogeneous)) << "checked homogeneous ODE reports Homogeneous";
    }
}

TEST(OdeEngineFirstOrder, BernoulliChecked) {
    auto P = SymbolicExpr::number(1);
    auto Q = SymbolicExpr::number(1);
    auto sol = solve_bernoulli_ode_checked(P, Q, 2, "x", "y");
    ASSERT_TRUE((sol.has_value())) << "checked Bernoulli ODE succeeds";
    if (sol) {
        EXPECT_TRUE((sol.value().general_solution != nullptr)) << "checked Bernoulli ODE returns solution";
        EXPECT_TRUE((sol.value().method_used == ODEType::Bernoulli)) << "checked Bernoulli ODE reports Bernoulli";
    }
}

TEST(OdeEngineFirstOrder, ExactChecked) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto M = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x), y);
    auto N = SymbolicExpr::add(
        x, SymbolicExpr::multiply(SymbolicExpr::number(2), y));
    auto sol = solve_exact_ode_checked(M, N, "x", "y");
    ASSERT_TRUE((sol.has_value())) << "checked exact ODE succeeds";
    if (sol) {
        EXPECT_TRUE((sol.value().general_solution != nullptr)) << "checked exact ODE returns solution";
        EXPECT_TRUE((sol.value().method_used == ODEType::Exact)) << "checked exact ODE reports Exact";
    }
}

TEST(OdeEngineFirstOrder, FirstOrderInvalid) {
    auto x = SymbolicExpr::variable("x");
    std::shared_ptr<SymbolicExpr> null_root;
    auto null_rhs = solve_homogeneous_ode_checked(null_root, "x", "y");
    EXPECT_TRUE((!null_rhs.has_value())) << "checked homogeneous ODE rejects null rhs";
    EXPECT_TRUE((null_rhs.error().code == CasErrc::InvalidArgument)) << "checked homogeneous ODE reports InvalidArgument for null rhs";

    auto same_vars = solve_homogeneous_ode_checked(x, "x", "x");
    EXPECT_TRUE((!same_vars.has_value())) << "checked homogeneous ODE rejects duplicate variable names";
    EXPECT_TRUE((same_vars.error().code == CasErrc::InvalidArgument)) << "checked homogeneous ODE reports InvalidArgument for duplicate variables";

    auto bad_n = solve_bernoulli_ode_checked(x, x, 1, "x", "y");
    EXPECT_TRUE((!bad_n.has_value())) << "checked Bernoulli ODE rejects n=1";
    EXPECT_TRUE((bad_n.error().code == CasErrc::InvalidArgument)) << "checked Bernoulli ODE reports InvalidArgument for n=1";

    auto null_exact = solve_exact_ode_checked(x, null_root, "x", "y");
    EXPECT_TRUE((!null_exact.has_value())) << "checked exact ODE rejects null N";
    EXPECT_TRUE((null_exact.error().code == CasErrc::InvalidArgument)) << "checked exact ODE reports InvalidArgument for null N";
}

TEST(OdeEngineFirstOrder, FirstOrderContext) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto rhs = SymbolicExpr::divide(y, x);

    LMCAS::CancellationToken cancellation;
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    cancellation.cancel();
    auto cancelled = solve_homogeneous_ode_checked(rhs, "x", "y",
                                                   cancelled_context);
    EXPECT_TRUE((!cancelled.has_value())) << "checked homogeneous ODE observes cancellation";
    EXPECT_TRUE((cancelled.error().code == CasErrc::Cancelled)) << "checked homogeneous ODE reports Cancelled";

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = solve_exact_ode_checked(y, x, "x", "y", limited_context);
    EXPECT_TRUE((!limited.has_value())) << "checked exact ODE observes exhausted step budget";
    EXPECT_TRUE((limited.error().code == CasErrc::ResourceLimit)) << "checked exact ODE reports ResourceLimit";
}

TEST(OdeEngineFirstOrder, IntegratingFactorExactness) {
    /**
     * @brief M = 2y、N = x 的偏导分别为 2、1，通过积分因子转为恰当方程。
     * (∂M/∂y-∂N/∂x)/N = (2-1)/x = 1/x，仅依赖 x；
     * 积分因子 mu(x) = exp(∫1/x dx) = exp(ln(x)) = x。
     */
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto M = SymbolicExpr::multiply(SymbolicExpr::number(2), y);
    auto N = x;

    EXPECT_FALSE((is_exact_ode(M, N, "x", "y"))) << "(2y, x) is NOT exact";

    auto mu = find_integrating_factor(M, N, "x", "y");
    EXPECT_TRUE((mu != nullptr)) << "integrating factor found for (2y, x)";

    if (mu) {
        auto M_new = SymbolicExpr::multiply(mu, M)->simplify();
        auto N_new = SymbolicExpr::multiply(mu, N)->simplify();
        auto dM_dy = M_new->differentiate("y");
        auto dN_dx = N_new->differentiate("x");
        if (dM_dy && dN_dx) {
            auto diff = SymbolicExpr::add(dM_dy,
                                          SymbolicExpr::multiply(SymbolicExpr::number(-1), dN_dx))
                            ->simplify();
            EXPECT_TRUE((diff->is_zero() || is_exact_ode(M_new, N_new, "x", "y"))) << "after multiplying by μ, equation becomes exact";
        }
    }
}

TEST(OdeEngineFirstOrder, IntegratingFactorSolution) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto M = SymbolicExpr::multiply(SymbolicExpr::number(2), y);
    auto N = x;

    auto sol = solve_exact_ode_checked(M, N, "x", "y").value();
    EXPECT_TRUE((sol.general_solution != nullptr)) << "solve_exact_ode handles non-exact with integrating factor";
}
