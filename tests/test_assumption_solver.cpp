
#include "test_common.hpp"
#include "assumption_context.hpp"
#include "solver.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "bigint.hpp"
#include "rational.hpp"
#include "numeric_evaluation.hpp"
#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace LMCAS;

/// Try to extract a numeric double value from a solution expression.
static bool try_numeric(const std::shared_ptr<SymbolicExpr> &expr, double &out) {
    if (!expr || !LMCAS::detail::node(expr)) {
        return false;
    }
    auto num = std::dynamic_pointer_cast<const NumberNode>(LMCAS::detail::node(expr));
    if (!num) {
        return false;
    }

    if (std::holds_alternative<BigInt>(num->value())) {
        out = std::get<BigInt>(num->value()).to_double();
        return true;
    }
    if (std::holds_alternative<Rational>(num->value())) {
        out = std::get<Rational>(num->value()).to_double();
        return true;
    }
    if (std::holds_alternative<lmmc_real_t>(num->value())) {
        out = std::get<lmmc_real_t>(num->value());
        return true;
    }
    return false;
}

/// Check if a solution set contains a numeric value (within tolerance).
static bool solutions_contain_value(
    const std::vector<std::shared_ptr<SymbolicExpr>> &solutions,
    double target, double tol = 1e-9) {
    for (const auto &sol : solutions) {
        double v = 0.0;
        if (try_numeric(sol, v)) {
            if (std::abs(v - target) < tol) {
                return true;
            }
        }
    }
    return false;
}

/// Check if any solution contains imaginary components (sqrt of negative).

/// Build equation x^2 - c = 0 as a SymbolicExpr (x^2 + (-c))
static std::shared_ptr<SymbolicExpr> make_x_squared_minus(const std::string &var, int c) {
    auto x = SymbolicExpr::variable(var);
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto eq = SymbolicExpr::add(x_sq, SymbolicExpr::number(-c));
    return eq;
}

/// Build equation x^2 + c = 0 as a SymbolicExpr (x^2 + c)
static std::shared_ptr<SymbolicExpr> make_x_squared_plus(const std::string &var, int c) {
    auto x = SymbolicExpr::variable(var);
    auto x_sq = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto eq = SymbolicExpr::add(x_sq, SymbolicExpr::number(c));
    return eq;
}

TEST(LmcasAssumptionSolver, XSquaredMinus4RealDomain) {
    auto eq = make_x_squared_minus("x", 4);

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    auto solutions = solve_with_assumptions_checked(eq, "x", &ctx).value();

    // Both 2 and -2 are real, so both should be returned
    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_TRUE((has_2)) << "x²-4=0 Real domain: contains x=2";
    EXPECT_TRUE((has_neg2)) << "x²-4=0 Real domain: contains x=-2";
    EXPECT_TRUE((solutions.size() >= 2)) << "x²-4=0 Real domain: at least 2 solutions";
}

TEST(LmcasAssumptionSolver, XSquaredMinus4PositiveInt) {
    auto eq = make_x_squared_minus("x", 4);

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::PositiveInt).has_value());

    auto result = solve_with_assumptions_checked(eq, "x", &ctx);
    EXPECT_TRUE((result && result.value().size() == 1 &&
                 solutions_contain_value(result.value(), 2.0)))
        << "x²-4=0 PositiveInt: exactly x=2";
}

TEST(LmcasAssumptionSolver, XSquaredMinus4NonnegativeSign) {
    auto eq = make_x_squared_minus("x", 4);

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonNegative).has_value());

    auto solutions = solve_with_assumptions_checked(eq, "x", &ctx).value();

    // Only x=2 satisfies NonNegative; x=-2 should be excluded
    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_TRUE((has_2)) << "x²-4=0 NonNegative sign: contains x=2";
    EXPECT_FALSE((has_neg2)) << "x²-4=0 NonNegative sign: does NOT contain x=-2";
}

TEST(LmcasAssumptionSolver, XSquaredPlus1RealDomain) {
    auto eq = make_x_squared_plus("x", 1);

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());

    auto result = solve_with_assumptions_checked(eq, "x", &ctx);
    EXPECT_TRUE((result && result.value().empty())) << "x²+1=0 Real domain: successful empty set";
}

TEST(LmcasAssumptionSolver, XSquaredPlus1NoContext) {
    auto eq = make_x_squared_plus("x", 1);

    // No context (nullptr) - all solutions returned unfiltered
    auto solutions = solve_with_assumptions_checked(eq, "x", nullptr).value();

    /// 默认求解路径返回复数域中的虚数解.
    EXPECT_TRUE((solutions.size() >= 1)) << "x²+1=0 no context: at least 1 solution returned (imaginary)";
}

TEST(LmcasAssumptionSolver, XSquaredMinus1PositiveSign) {
    auto eq = make_x_squared_minus("x", 1);

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Positive).has_value());

    auto solutions = solve_with_assumptions_checked(eq, "x", &ctx).value();

    bool has_1 = solutions_contain_value(solutions, 1.0);
    bool has_neg1 = solutions_contain_value(solutions, -1.0);

    EXPECT_TRUE((has_1)) << "x²-1=0 Positive sign: contains x=1";
    EXPECT_FALSE((has_neg1)) << "x²-1=0 Positive sign: does NOT contain x=-1";
}

TEST(LmcasAssumptionSolver, XSquaredMinus1NegativeSign) {
    auto eq = make_x_squared_minus("x", 1);

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_sign("x", Sign::Negative).has_value());

    auto solutions = solve_with_assumptions_checked(eq, "x", &ctx).value();

    bool has_1 = solutions_contain_value(solutions, 1.0);
    bool has_neg1 = solutions_contain_value(solutions, -1.0);

    EXPECT_FALSE((has_1)) << "x²-1=0 Negative sign: does NOT contain x=1";
    EXPECT_TRUE((has_neg1)) << "x²-1=0 Negative sign: contains x=-1";
}

TEST(LmcasAssumptionSolver, NoContextAllSolutionsReturned) {
    // x^2 - 4 = 0 -> x=2, x=-2
    auto eq = make_x_squared_minus("x", 4);

    auto solutions = solve_with_assumptions_checked(eq, "x", nullptr).value();

    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_TRUE((has_2)) << "x²-4=0 no context: contains x=2";
    EXPECT_TRUE((has_neg2)) << "x²-4=0 no context: contains x=-2";
    EXPECT_TRUE((solutions.size() >= 2)) << "x²-4=0 no context: at least 2 solutions";
}

TEST(LmcasAssumptionSolver, NoContextXSquaredMinus1) {
    auto eq = make_x_squared_minus("x", 1);

    auto solutions = solve_with_assumptions_checked(eq, "x", nullptr).value();

    bool has_1 = solutions_contain_value(solutions, 1.0);
    bool has_neg1 = solutions_contain_value(solutions, -1.0);

    EXPECT_TRUE((has_1)) << "x²-1=0 no context: contains x=1";
    EXPECT_TRUE((has_neg1)) << "x²-1=0 no context: contains x=-1";
}

TEST(LmcasAssumptionSolver, AllSolutionsFilteredPositiveInt) {
    auto eq = make_x_squared_minus("x", 2);

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::PositiveInt).has_value());

    auto result = solve_with_assumptions_checked(eq, "x", &ctx);
    EXPECT_TRUE((result && result.value().empty())) << "x²−2=0 PositiveInt: successful empty set";
}

TEST(LmcasAssumptionSolver, NaturalDomainExcludesNegative) {
    auto eq = make_x_squared_minus("x", 4);

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Natural).has_value());

    auto solutions = solve_with_assumptions_checked(eq, "x", &ctx).value();

    // Natural = non-negative integers. x=2 is valid, x=-2 is not.
    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_TRUE((has_2)) << "x²-4=0 Natural domain: contains x=2";
    EXPECT_FALSE((has_neg2)) << "x²-4=0 Natural domain: does NOT contain x=-2";
}

TEST(LmcasAssumptionSolver, IntegerDomainBothReturned) {
    auto eq = make_x_squared_minus("x", 4);

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Integer).has_value());

    auto result = solve_with_assumptions_checked(eq, "x", &ctx);
    EXPECT_TRUE((result && result.value().size() == 2 &&
                 solutions_contain_value(result.value(), 2.0) &&
                 solutions_contain_value(result.value(), -2.0)))
        << "x²-4=0 Integer: exactly x=2 and x=-2";
}

TEST(LmcasAssumptionSolver, ComplexDomainNoFiltering) {
    // x^2 + 1 = 0 with Complex domain -> imaginary solutions should be kept
    auto eq = make_x_squared_plus("x", 1);

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Complex).has_value());

    auto solutions = solve_with_assumptions_checked(eq, "x", &ctx).value();

    // Complex is the default/least restrictive - no filtering
    EXPECT_TRUE((solutions.size() >= 1)) << "x²+1=0 Complex domain: solutions returned (no filtering)";
}

TEST(LmcasAssumptionSolver, NonpositiveSignFiltering) {
    auto eq = make_x_squared_minus("x", 4);

    AssumptionContext ctx;
    ASSERT_TRUE(ctx.assume_domain("x", Domain::Real).has_value());
    ASSERT_TRUE(ctx.assume_sign("x", Sign::NonPositive).has_value());

    auto solutions = solve_with_assumptions_checked(eq, "x", &ctx).value();

    bool has_2 = solutions_contain_value(solutions, 2.0);
    bool has_neg2 = solutions_contain_value(solutions, -2.0);

    EXPECT_FALSE((has_2)) << "x²-4=0 NonPositive sign: does NOT contain x=2";
    EXPECT_TRUE((has_neg2)) << "x²-4=0 NonPositive sign: contains x=-2";
}

TEST(LmcasAssumptionSolver, CubicExactRootDomains) {
    auto x = SymbolicExpr::variable("x");
    auto equation = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(3)),
        SymbolicExpr::number(-2));

    AssumptionContext real;
    EXPECT_TRUE((real.assume_domain("x", Domain::Real).has_value())) << "real domain assumption succeeds";
    auto real_roots = solve_with_assumptions_checked(equation, "x", &real);
    EXPECT_TRUE((real_roots && real_roots.value().size() == 1)) << "x^3-2 has exactly one real solution";
    if (real_roots && real_roots.value().size() == 1) {
        auto value = evaluate_numeric(*real_roots.value().front());
        EXPECT_TRUE((value && std::isfinite(value.value().value) &&
                     std::abs(value.value().value - std::cbrt(2.0)) < 1e-10))
            << "the retained real solution is cube-root two";
    }

    AssumptionContext complex;
    EXPECT_TRUE((complex.assume_domain("x", Domain::Complex).has_value())) << "complex domain assumption succeeds";
    auto complex_roots = solve_with_assumptions_checked(equation, "x", &complex);
    EXPECT_TRUE((complex_roots && complex_roots.value().size() == 3)) << "complex domain retains all three cubic solutions";

    auto unrestricted = solve_with_assumptions_checked(equation, "x", nullptr);
    EXPECT_TRUE((unrestricted && unrestricted.value().size() == 3)) << "absent assumptions retain all three cubic solutions";
}

TEST(LmcasAssumptionSolver, CheckedAssumptionSolverContract) {
    auto invalid = solve_with_assumptions_checked(nullptr, "x");
    EXPECT_TRUE((!invalid &&
                 invalid.error().code == CasErrc::InvalidArgument))
        << "null assumption-aware equation is invalid";

    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext context(limits);
    auto limited = solve_with_assumptions_checked(
        make_x_squared_minus("x", 1), "x", nullptr, context);
    EXPECT_TRUE((!limited &&
                 limited.error().code == CasErrc::ResourceLimit))
        << "assumption-aware solve preserves exhausted budget";
}
