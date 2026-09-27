#include "test_inequality_support.hpp"
#include <stdexcept>

TEST(InequalityExamples, CheckedPositiveSlope) {
    auto greater = InequalitySolver::solve_inequality_checked(
        linear(2, -3), InequalityType::GreaterThan, "x");
    ASSERT_TRUE((greater.has_value())) << "checked 2x-3 > 0 succeeds";
    if (greater) {
        EXPECT_TRUE((!greater.value().contains(1.5))) << "checked strict affine inequality excludes its exact boundary";
        EXPECT_TRUE((greater.value().contains(2.0))) << "checked positive-slope affine inequality selects the upper ray";
    }
}

TEST(InequalityExamples, CheckedNegativeSlope) {
    auto negative_slope = InequalitySolver::solve_inequality_checked(
        linear(-2, 4), InequalityType::LessEqual, "x");
    ASSERT_TRUE((negative_slope.has_value())) << "checked -2x+4 <= 0 succeeds";
    if (negative_slope) {
        EXPECT_TRUE((negative_slope.value().contains(2.0))) << "checked non-strict affine inequality includes its boundary";
        EXPECT_TRUE((negative_slope.value().contains(3.0))) << "checked negative-slope inequality reverses direction exactly";
        EXPECT_TRUE((!negative_slope.value().contains(1.0))) << "checked negative-slope inequality excludes the other ray";
    }
}

TEST(InequalityExamples, CheckedConstants) {
    auto true_constant = InequalitySolver::solve_inequality_checked(
        SymbolicExpr::number(Rational(1, 3)),
        InequalityType::GreaterThan, "x");
    EXPECT_TRUE((true_constant && true_constant.value().is_entire_line())) << "checked true constant inequality returns the entire real line";

    auto false_constant = InequalitySolver::solve_inequality_checked(
        SymbolicExpr::number(0), InequalityType::LessThan, "x");
    EXPECT_TRUE((false_constant && false_constant.value().is_empty())) << "checked false constant inequality returns the empty set";
}

TEST(InequalityExamples, CheckedCubic) {
    auto x = SymbolicExpr::variable("x");

    auto cubic = SymbolicExpr::power(x, SymbolicExpr::number(3));
    auto cubic_result = InequalitySolver::solve_inequality_checked(
        cubic, InequalityType::GreaterEqual, "x");
    EXPECT_TRUE((cubic_result &&
                 cubic_result.value().contains(0.0) &&
                 cubic_result.value().contains(2.0) &&
                 !cubic_result.value().contains(-1.0)))
        << "checked inequality solves exact cubic sign charts";
}

TEST(InequalityExamples, CheckedInvalidExpression) {
    auto null_result = InequalitySolver::solve_inequality_checked(
        nullptr, InequalityType::GreaterThan, "x");
    EXPECT_TRUE((!null_result && null_result.error().code == CasErrc::InvalidArgument)) << "checked inequality solving rejects null expressions";
}

TEST(InequalityExamples, CheckedCancellation) {
    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled_context(ResourceLimits{}, cancellation);
    auto cancelled = InequalitySolver::solve_inequality_checked(
        linear(1, 0), InequalityType::GreaterThan, "x", cancelled_context);
    EXPECT_TRUE((!cancelled && cancelled.error().code == CasErrc::Cancelled)) << "checked inequality solving observes cancellation";
}

TEST(InequalityExamples, CheckedStepBudget) {
    ResourceLimits limits;
    limits.max_steps = 0;
    ComputationContext limited_context(limits);
    auto limited = InequalitySolver::solve_inequality_checked(
        linear(1, 0), InequalityType::GreaterThan, "x", limited_context);
    EXPECT_TRUE((!limited && limited.error().code == CasErrc::ResourceLimit)) << "checked inequality solving observes the step budget";
}

TEST(InequalityExamples, CheckedApproximateCoefficients) {
    auto x = SymbolicExpr::variable("x");

    auto approximate = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(0.5), x),
        SymbolicExpr::number(1));
    auto approximate_result = InequalitySolver::solve_inequality_checked(
        approximate, InequalityType::GreaterThan, "x");
    EXPECT_TRUE((!approximate_result &&
                 approximate_result.error().code == CasErrc::Inconclusive))
        << "checked affine solving does not relabel approximate coefficients as exact";
}

TEST(InequalityExamples, CheckedSymbolicCoefficients) {
    auto x = SymbolicExpr::variable("x");

    auto parameterized = SymbolicExpr::add(x, SymbolicExpr::variable("a"));
    auto parameterized_result = InequalitySolver::solve_inequality_checked(
        parameterized, InequalityType::GreaterThan, "x");
    EXPECT_TRUE((!parameterized_result &&
                 parameterized_result.error().code == CasErrc::Inconclusive))
        << "checked affine solving reports symbolic coefficients as unsupported";
}

TEST(InequalityExamples, CheckedRationalQuadratic) {
    auto x = SymbolicExpr::variable("x");
    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));

    auto rational_roots = SymbolicExpr::add(x2, SymbolicExpr::number(-4));
    auto outside = InequalitySolver::solve_inequality_checked(
        rational_roots, InequalityType::GreaterEqual, "x");
    ASSERT_TRUE((outside.has_value())) << (outside ? "checked x^2-4 >= 0 succeeds exactly"
                                                   : "checked x^2-4 >= 0 failed: " + outside.error().message);
    if (outside) {
        EXPECT_TRUE((outside.value().contains(-2.0) && outside.value().contains(2.0))) << "checked non-strict quadratic includes both roots";
        EXPECT_TRUE((outside.value().contains(-3.0) && outside.value().contains(3.0))) << "checked positive quadratic selects both outside rays";
        EXPECT_TRUE((!outside.value().contains(0.0))) << "checked positive quadratic excludes its negative interior";
    }
}

TEST(InequalityExamples, CheckedIrrationalQuadratic) {
    auto x = SymbolicExpr::variable("x");
    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));

    auto irrational_roots = SymbolicExpr::add(x2, SymbolicExpr::number(-2));
    auto inside = InequalitySolver::solve_inequality_checked(
        irrational_roots, InequalityType::LessThan, "x");
    ASSERT_TRUE((inside.has_value())) << (inside ? "checked x^2-2 < 0 verifies algebraic boundaries"
                                                 : "checked x^2-2 < 0 failed: " + inside.error().message);
    if (inside) {
        EXPECT_TRUE((inside.value().contains(0.0))) << "checked irrational-root quadratic includes its interior";
        EXPECT_TRUE((!inside.value().contains(2.0))) << "checked irrational-root quadratic excludes its exterior";
        EXPECT_TRUE((!inside.value().contains(std::sqrt(2.0)))) << "checked strict irrational-root quadratic excludes its boundary";
    }
}

TEST(InequalityExamples, CheckedPositiveQuadratic) {
    auto x = SymbolicExpr::variable("x");
    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));

    auto no_real_roots = SymbolicExpr::add(x2, SymbolicExpr::number(1));
    auto always_positive = InequalitySolver::solve_inequality_checked(
        no_real_roots, InequalityType::GreaterThan, "x");
    EXPECT_TRUE((always_positive && always_positive.value().is_entire_line())) << "checked positive quadratic with negative discriminant is always positive";
    auto never_negative = InequalitySolver::solve_inequality_checked(
        no_real_roots, InequalityType::LessEqual, "x");
    EXPECT_TRUE((never_negative && never_negative.value().is_empty())) << "checked positive quadratic with negative discriminant is never non-positive";
}

TEST(InequalityExamples, CheckedRepeatedQuadratic) {
    auto x = SymbolicExpr::variable("x");
    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));

    auto repeated = SymbolicExpr::add(
        SymbolicExpr::add(x2, SymbolicExpr::multiply(SymbolicExpr::number(-2), x)),
        SymbolicExpr::number(1));
    auto repeated_nonpositive = InequalitySolver::solve_inequality_checked(
        repeated, InequalityType::LessEqual, "x");
    bool point_matches = false;
    if (repeated_nonpositive) {
        const auto &solution = repeated_nonpositive.value();
        point_matches = solution.intervals().size() == 1 &&
                        solution.contains(1.0) && !solution.contains(0.0);
    }
    EXPECT_TRUE((point_matches)) << "checked repeated-root quadratic returns exactly the root point";
    auto repeated_positive = InequalitySolver::solve_inequality_checked(
        repeated, InequalityType::GreaterThan, "x");
    bool punctured_line_matches = false;
    if (repeated_positive) {
        const auto &solution = repeated_positive.value();
        punctured_line_matches = !solution.contains(1.0) &&
                                 solution.contains(0.0) && solution.contains(2.0);
    }
    EXPECT_TRUE((punctured_line_matches)) << "checked strict repeated-root quadratic returns the punctured line";
}

TEST(InequalityExamples, CheckedDownwardQuadratic) {
    auto x = SymbolicExpr::variable("x");
    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto rational_roots = SymbolicExpr::add(x2, SymbolicExpr::number(-4));

    auto downward = SymbolicExpr::multiply(SymbolicExpr::number(-1), rational_roots);
    auto downward_positive = InequalitySolver::solve_inequality_checked(
        downward, InequalityType::GreaterThan, "x");
    EXPECT_TRUE((downward_positive && downward_positive.value().contains(0.0) &&
                 !downward_positive.value().contains(3.0)))
        << (downward_positive
                ? "checked downward quadratic reverses the sign regions exactly"
                : "checked downward quadratic failed: " +
                      downward_positive.error().message);
}

TEST(InequalityExamples, LinearExample) {
    auto expr = linear(2, -3);
    auto result = InequalitySolver::solve_inequality(expr, InequalityType::GreaterThan, "x");

    EXPECT_TRUE((!result.is_empty())) << "2x - 3 > 0 should not be empty";
    EXPECT_TRUE((!result.contains(0.0))) << "2x - 3 > 0: x=0 not in solution";
    EXPECT_TRUE((!result.contains(1.0))) << "2x - 3 > 0: x=1 not in solution";
    EXPECT_TRUE((!result.contains(1.5))) << "2x - 3 > 0: x=1.5 (root) not in solution (strict)";
    EXPECT_TRUE((result.contains(2.0))) << "2x - 3 > 0: x=2 in solution";
    EXPECT_TRUE((result.contains(100.0))) << "2x - 3 > 0: x=100 in solution";
    EXPECT_TRUE((!result.contains(-5.0))) << "2x - 3 > 0: x=-5 not in solution";
}

TEST(InequalityExamples, QuadraticExample) {
    auto x = SymbolicExpr::variable("x");

    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto expr = SymbolicExpr::add(x2, SymbolicExpr::number(-4));
    auto result = InequalitySolver::solve_inequality(expr, InequalityType::GreaterEqual, "x");

    EXPECT_TRUE((!result.is_empty())) << "x^2 - 4 >= 0 should not be empty";

    EXPECT_TRUE((result.contains(-10.0))) << "x^2 - 4 >= 0: x=-10 in solution";
    EXPECT_TRUE((result.contains(-3.0))) << "x^2 - 4 >= 0: x=-3 in solution";
    EXPECT_TRUE((result.contains(3.0))) << "x^2 - 4 >= 0: x=3 in solution";
    EXPECT_TRUE((result.contains(10.0))) << "x^2 - 4 >= 0: x=10 in solution";

    EXPECT_TRUE((!result.contains(0.0))) << "x^2 - 4 >= 0: x=0 not in solution";
    EXPECT_TRUE((!result.contains(1.0))) << "x^2 - 4 >= 0: x=1 not in solution";
    EXPECT_TRUE((!result.contains(-1.0))) << "x^2 - 4 >= 0: x=-1 not in solution";
    EXPECT_TRUE((!result.contains(1.9))) << "x^2 - 4 >= 0: x=1.9 not in solution";
    EXPECT_TRUE((!result.contains(-1.9))) << "x^2 - 4 >= 0: x=-1.9 not in solution";

    EXPECT_TRUE((result.contains(2.0))) << "x^2 - 4 >= 0: x=2 in solution (non-strict, root)";
    EXPECT_TRUE((result.contains(-2.0))) << "x^2 - 4 >= 0: x=-2 in solution (non-strict, root)";
}

TEST(InequalityExamples, CubicExample) {
    auto x = SymbolicExpr::variable("x");

    auto x3 = SymbolicExpr::power(x, SymbolicExpr::number(3));
    auto expr = SymbolicExpr::add(x3, negate(x));
    auto result = InequalitySolver::solve_inequality(expr, InequalityType::LessThan, "x");

    EXPECT_TRUE((!result.is_empty())) << "x^3 - x < 0 should not be empty";
    EXPECT_TRUE((result.contains(-5.0))) << "x^3 - x < 0: x=-5 in solution";
    EXPECT_TRUE((result.contains(-2.0))) << "x^3 - x < 0: x=-2 in solution";
    EXPECT_TRUE((!result.contains(-1.0))) << "x^3 - x < 0: x=-1 not in solution (strict, root)";
    EXPECT_TRUE((result.contains(std::nextafter(-1.0, -2.0)) &&
                 !result.contains(std::nextafter(-1.0, 0.0))))
        << "strict cubic boundary separates adjacent representable values";
    EXPECT_TRUE((!result.contains(1.0))) << "x^3 - x < 0: x=1 not in solution (strict, root)";
    EXPECT_TRUE((result.contains(0.5))) << "x^3 - x < 0: x=0.5 in solution";
    EXPECT_TRUE((!result.contains(2.0))) << "x^3 - x < 0: x=2 not in solution";
    EXPECT_TRUE((!result.contains(-0.5))) << "x^3 - x < 0: x=-0.5 not in solution (between -1 and 0)";

    EXPECT_TRUE((result.contains(0.001))) << "x^3 - x < 0: x=0.001 in solution (just above 0)";
}

TEST(InequalityExamples, RepeatedRootExample) {
    auto x = SymbolicExpr::variable("x");

    auto x_minus_1 = SymbolicExpr::add(x, SymbolicExpr::number(-1));
    auto x_minus_1_sq = SymbolicExpr::power(x_minus_1, SymbolicExpr::number(2));
    auto x_plus_2 = SymbolicExpr::add(x, SymbolicExpr::number(2));
    auto expr = SymbolicExpr::multiply(x_minus_1_sq, x_plus_2)->expand();

    auto result = InequalitySolver::solve_inequality(expr, InequalityType::GreaterThan, "x");

    EXPECT_TRUE((!result.is_empty())) << "(x-1)^2*(x+2) > 0 should not be empty";
    EXPECT_TRUE((!result.contains(-2.0))) << "(x-1)^2*(x+2) > 0: x=-2 not in solution (root, strict)";
    EXPECT_TRUE((!result.contains(-3.0))) << "(x-1)^2*(x+2) > 0: x=-3 not in solution";
    EXPECT_TRUE((result.contains(0.0))) << "(x-1)^2*(x+2) > 0: x=0 in solution";
    EXPECT_TRUE((!result.contains(1.0))) << "(x-1)^2*(x+2) > 0: x=1 not in solution (root, strict)";
    EXPECT_TRUE((result.contains(2.0))) << "(x-1)^2*(x+2) > 0: x=2 in solution";
    EXPECT_TRUE((result.contains(10.0))) << "(x-1)^2*(x+2) > 0: x=10 in solution";
    EXPECT_TRUE((result.contains(-1.0))) << "(x-1)^2*(x+2) > 0: x=-1 in solution";
}

TEST(InequalityExamples, RationalStrictExample) {
    auto x = SymbolicExpr::variable("x");

    auto numerator = SymbolicExpr::add(x, SymbolicExpr::number(-1));
    auto denominator = SymbolicExpr::add(x, SymbolicExpr::number(2));

    auto result = InequalitySolver::solve_rational_inequality(
        numerator, denominator, InequalityType::GreaterThan, "x");

    EXPECT_TRUE((!result.is_empty())) << "(x-1)/(x+2) > 0 should not be empty";
    EXPECT_TRUE((result.contains(-5.0))) << "(x-1)/(x+2) > 0: x=-5 in solution";
    EXPECT_TRUE((!result.contains(-2.0))) << "(x-1)/(x+2) > 0: x=-2 not in solution (den root)";
    EXPECT_TRUE((!result.contains(0.0))) << "(x-1)/(x+2) > 0: x=0 not in solution";
    EXPECT_TRUE((!result.contains(1.0))) << "(x-1)/(x+2) > 0: x=1 not in solution (strict, num root)";
    EXPECT_TRUE((result.contains(2.0))) << "(x-1)/(x+2) > 0: x=2 in solution";
    EXPECT_TRUE((result.contains(100.0))) << "(x-1)/(x+2) > 0: x=100 in solution";
}

TEST(InequalityExamples, RationalClosedExample) {
    auto x = SymbolicExpr::variable("x");

    auto numerator = SymbolicExpr::add(x, SymbolicExpr::number(-1));
    auto denominator = SymbolicExpr::add(x, SymbolicExpr::number(2));

    auto result = InequalitySolver::solve_rational_inequality(
        numerator, denominator, InequalityType::GreaterEqual, "x");

    EXPECT_TRUE((!result.is_empty())) << "(x-1)/(x+2) >= 0 should not be empty";
    EXPECT_TRUE((result.contains(-5.0))) << "(x-1)/(x+2) >= 0: x=-5 in solution";
    EXPECT_TRUE((!result.contains(-2.0))) << "(x-1)/(x+2) >= 0: x=-2 not in solution (den root excluded)";
    EXPECT_TRUE((!result.contains(0.0))) << "(x-1)/(x+2) >= 0: x=0 not in solution";
    EXPECT_TRUE((result.contains(1.0))) << "(x-1)/(x+2) >= 0: x=1 in solution (num root, non-strict)";
    EXPECT_TRUE((result.contains(2.0))) << "(x-1)/(x+2) >= 0: x=2 in solution";
}

TEST(InequalityExamples, SystemExample) {
    auto x = SymbolicExpr::variable("x");

    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto expr1 = SymbolicExpr::add(x2, SymbolicExpr::number(-4));

    auto expr2 = SymbolicExpr::add(x, SymbolicExpr::number(-5));

    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> system = {
        {expr1, InequalityType::GreaterThan},
        {expr2, InequalityType::LessThan}};

    auto result = InequalitySolver::solve_inequalities(system, "x");

    EXPECT_TRUE((!result.is_empty())) << "{x^2-4>0, x<5} should not be empty";
    EXPECT_TRUE((result.contains(-10.0))) << "{x^2-4>0, x<5}: x=-10 in solution";
    EXPECT_TRUE((result.contains(-3.0))) << "{x^2-4>0, x<5}: x=-3 in solution";
    EXPECT_TRUE((!result.contains(-2.0))) << "{x^2-4>0, x<5}: x=-2 not in solution (strict root)";
    EXPECT_TRUE((!result.contains(0.0))) << "{x^2-4>0, x<5}: x=0 not in solution";
    EXPECT_TRUE((!result.contains(2.0))) << "{x^2-4>0, x<5}: x=2 not in solution (strict root)";
    EXPECT_TRUE((result.contains(3.0))) << "{x^2-4>0, x<5}: x=3 in solution";
    EXPECT_TRUE((result.contains(4.0))) << "{x^2-4>0, x<5}: x=4 in solution";
    EXPECT_TRUE((!result.contains(5.0))) << "{x^2-4>0, x<5}: x=5 not in solution (strict)";
    EXPECT_TRUE((!result.contains(6.0))) << "{x^2-4>0, x<5}: x=6 not in solution";
}

TEST(InequalityExamples, ZeroStrictPositive) {
    auto zero_expr = SymbolicExpr::number(0);
    auto result = InequalitySolver::solve_inequality(zero_expr, InequalityType::GreaterThan, "x");

    EXPECT_TRUE((result.is_empty())) << "0 > 0 should be empty";
    EXPECT_TRUE((!result.contains(0.0))) << "0 > 0: contains nothing";
    EXPECT_TRUE((!result.contains(1.0))) << "0 > 0: contains nothing";
}

TEST(InequalityExamples, ZeroNonnegative) {
    auto zero_expr = SymbolicExpr::number(0);
    auto result = InequalitySolver::solve_inequality(zero_expr, InequalityType::GreaterEqual, "x");

    EXPECT_TRUE((result.is_entire_line())) << "0 >= 0 should be entire line";
    EXPECT_TRUE((result.contains(0.0))) << "0 >= 0: contains 0";
    EXPECT_TRUE((result.contains(-1000.0))) << "0 >= 0: contains -1000";
    EXPECT_TRUE((result.contains(1000.0))) << "0 >= 0: contains 1000";
}

TEST(InequalityExamples, ZeroStrictNegative) {
    auto zero_expr = SymbolicExpr::number(0);
    auto result = InequalitySolver::solve_inequality(zero_expr, InequalityType::LessThan, "x");

    EXPECT_TRUE((result.is_empty())) << "0 < 0 should be empty";
}

TEST(InequalityExamples, ZeroNonpositive) {
    auto zero_expr = SymbolicExpr::number(0);
    auto result = InequalitySolver::solve_inequality(zero_expr, InequalityType::LessEqual, "x");

    EXPECT_TRUE((result.is_entire_line())) << "0 <= 0 should be entire line";
}

TEST(InequalityExamples, HugeConstant) {
    std::string huge_digits = "1" + std::string(400, '0');
    auto huge_positive = SymbolicExpr::number(BigInt(huge_digits));
    auto positive_result = InequalitySolver::solve_inequality(
        huge_positive, InequalityType::GreaterThan, "x");
    EXPECT_TRUE((positive_result.is_entire_line())) << "huge exact positive constant > 0 should be entire line";

    auto huge_negative = SymbolicExpr::number(BigInt("-" + huge_digits));
    auto negative_result = InequalitySolver::solve_inequality(
        huge_negative, InequalityType::GreaterThan, "x");
    EXPECT_TRUE((negative_result.is_empty())) << "huge exact negative constant > 0 should be empty";
}

TEST(InequalityExamples, UnsupportedStrict) {
    auto x = SymbolicExpr::variable("x");

    auto sin_x = SymbolicExpr::sin(x);
    bool rejected = false;
    try {
        (void)InequalitySolver::solve_inequality(sin_x, InequalityType::GreaterThan, "x");
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    EXPECT_TRUE((rejected)) << "non-polynomial strict inequality is explicitly rejected, not an empty solution";
}

TEST(InequalityExamples, UnsupportedNonstrict) {
    auto x = SymbolicExpr::variable("x");

    auto sin_x = SymbolicExpr::sin(x);
    bool rejected = false;
    try {
        (void)InequalitySolver::solve_inequality(sin_x, InequalityType::GreaterEqual, "x");
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    EXPECT_TRUE((rejected)) << "non-polynomial non-strict inequality is explicitly rejected, not an empty solution";
}
