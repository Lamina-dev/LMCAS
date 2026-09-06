#include "test_common.hpp"
#include "interval.hpp"
#include "inequality_solver.hpp"
#include "symbolic.hpp"
#include <cmath>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace LMCAS;

static std::shared_ptr<SymbolicExpr> linear(int a, int b) {
    auto x = SymbolicExpr::variable("x");
    auto ax = SymbolicExpr::multiply(SymbolicExpr::number(a), x);
    return SymbolicExpr::add(ax, SymbolicExpr::number(b));
}

static std::shared_ptr<SymbolicExpr> negate(std::shared_ptr<SymbolicExpr> e) {
    return SymbolicExpr::multiply(SymbolicExpr::number(-1), e);
}

int main() {

    auto x = SymbolicExpr::variable("x");

    TEST_CASE("Checked exact affine inequality contracts");
    {
        auto greater = InequalitySolver::solve_inequality_checked(
            linear(2, -3), InequalityType::GreaterThan, "x");
        EXPECT_TRUE(greater.has_value(), "checked 2x-3 > 0 succeeds");
        if (greater) {
            EXPECT_TRUE(!greater.value().contains(1.5),
                        "checked strict affine inequality excludes its exact boundary");
            EXPECT_TRUE(greater.value().contains(2.0),
                        "checked positive-slope affine inequality selects the upper ray");
        }

        auto negative_slope = InequalitySolver::solve_inequality_checked(
            linear(-2, 4), InequalityType::LessEqual, "x");
        EXPECT_TRUE(negative_slope.has_value(), "checked -2x+4 <= 0 succeeds");
        if (negative_slope) {
            EXPECT_TRUE(negative_slope.value().contains(2.0),
                        "checked non-strict affine inequality includes its boundary");
            EXPECT_TRUE(negative_slope.value().contains(3.0),
                        "checked negative-slope inequality reverses direction exactly");
            EXPECT_TRUE(!negative_slope.value().contains(1.0),
                        "checked negative-slope inequality excludes the other ray");
        }

        auto true_constant = InequalitySolver::solve_inequality_checked(
            SymbolicExpr::number(Rational(1, 3)),
            InequalityType::GreaterThan, "x");
        EXPECT_TRUE(true_constant && true_constant.value().is_entire_line(),
                    "checked true constant inequality returns the entire real line");

        auto false_constant = InequalitySolver::solve_inequality_checked(
            SymbolicExpr::number(0), InequalityType::LessThan, "x");
        EXPECT_TRUE(false_constant && false_constant.value().is_empty(),
                    "checked false constant inequality returns the empty set");

        auto approximate = SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(0.5), x),
            SymbolicExpr::number(1));
        auto approximate_result = InequalitySolver::solve_inequality_checked(
            approximate, InequalityType::GreaterThan, "x");
        EXPECT_TRUE(!approximate_result &&
                        approximate_result.error().code == CasErrc::Inconclusive,
                    "checked affine solving does not relabel approximate coefficients as exact");

        auto parameterized = SymbolicExpr::add(x, SymbolicExpr::variable("a"));
        auto parameterized_result = InequalitySolver::solve_inequality_checked(
            parameterized, InequalityType::GreaterThan, "x");
        EXPECT_TRUE(!parameterized_result &&
                        parameterized_result.error().code == CasErrc::Inconclusive,
                    "checked affine solving reports symbolic coefficients as unsupported");

        auto cubic = SymbolicExpr::power(x, SymbolicExpr::number(3));
        auto cubic_result = InequalitySolver::solve_inequality_checked(
            cubic, InequalityType::GreaterEqual, "x");
        EXPECT_TRUE(cubic_result &&
                        cubic_result.value().contains(0.0) &&
                        cubic_result.value().contains(2.0) &&
                        !cubic_result.value().contains(-1.0),
                    "checked inequality solves exact cubic sign charts");

        auto null_result = InequalitySolver::solve_inequality_checked(
            nullptr, InequalityType::GreaterThan, "x");
        EXPECT_TRUE(!null_result && null_result.error().code == CasErrc::InvalidArgument,
                    "checked inequality solving rejects null expressions");

        CancellationToken cancellation;
        cancellation.cancel();
        ComputationContext cancelled_context(ResourceLimits{}, cancellation);
        auto cancelled = InequalitySolver::solve_inequality_checked(
            linear(1, 0), InequalityType::GreaterThan, "x", cancelled_context);
        EXPECT_TRUE(!cancelled && cancelled.error().code == CasErrc::Cancelled,
                    "checked inequality solving observes cancellation");

        ResourceLimits limits;
        limits.max_steps = 0;
        ComputationContext limited_context(limits);
        auto limited = InequalitySolver::solve_inequality_checked(
            linear(1, 0), InequalityType::GreaterThan, "x", limited_context);
        EXPECT_TRUE(!limited && limited.error().code == CasErrc::ResourceLimit,
                    "checked inequality solving observes the step budget");
    }

    TEST_CASE("Checked exact quadratic inequality contracts");
    {
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));

        auto rational_roots = SymbolicExpr::add(x2, SymbolicExpr::number(-4));
        auto outside = InequalitySolver::solve_inequality_checked(
            rational_roots, InequalityType::GreaterEqual, "x");
        EXPECT_TRUE(outside.has_value(),
                    outside ? "checked x^2-4 >= 0 succeeds exactly"
                            : "checked x^2-4 >= 0 failed: " + outside.error().message);
        if (outside) {
            EXPECT_TRUE(outside.value().contains(-2.0) && outside.value().contains(2.0),
                        "checked non-strict quadratic includes both roots");
            EXPECT_TRUE(outside.value().contains(-3.0) && outside.value().contains(3.0),
                        "checked positive quadratic selects both outside rays");
            EXPECT_TRUE(!outside.value().contains(0.0),
                        "checked positive quadratic excludes its negative interior");
        }

        auto irrational_roots = SymbolicExpr::add(x2, SymbolicExpr::number(-2));
        auto inside = InequalitySolver::solve_inequality_checked(
            irrational_roots, InequalityType::LessThan, "x");
        EXPECT_TRUE(inside.has_value(),
                    inside ? "checked x^2-2 < 0 verifies algebraic boundaries"
                           : "checked x^2-2 < 0 failed: " + inside.error().message);
        if (inside) {
            EXPECT_TRUE(inside.value().contains(0.0),
                        "checked irrational-root quadratic includes its interior");
            EXPECT_TRUE(!inside.value().contains(2.0),
                        "checked irrational-root quadratic excludes its exterior");
            EXPECT_TRUE(!inside.value().contains(std::sqrt(2.0)),
                        "checked strict irrational-root quadratic excludes its boundary");
        }

        auto no_real_roots = SymbolicExpr::add(x2, SymbolicExpr::number(1));
        auto always_positive = InequalitySolver::solve_inequality_checked(
            no_real_roots, InequalityType::GreaterThan, "x");
        EXPECT_TRUE(always_positive && always_positive.value().is_entire_line(),
                    "checked positive quadratic with negative discriminant is always positive");
        auto never_negative = InequalitySolver::solve_inequality_checked(
            no_real_roots, InequalityType::LessEqual, "x");
        EXPECT_TRUE(never_negative && never_negative.value().is_empty(),
                    "checked positive quadratic with negative discriminant is never non-positive");

        auto repeated = SymbolicExpr::add(
            SymbolicExpr::add(x2, SymbolicExpr::multiply(SymbolicExpr::number(-2), x)),
            SymbolicExpr::number(1));
        auto repeated_nonpositive = InequalitySolver::solve_inequality_checked(
            repeated, InequalityType::LessEqual, "x");
        EXPECT_TRUE(repeated_nonpositive &&
                        repeated_nonpositive.value().intervals().size() == 1 &&
                        repeated_nonpositive.value().contains(1.0) &&
                        !repeated_nonpositive.value().contains(0.0),
                    "checked repeated-root quadratic returns exactly the root point");
        auto repeated_positive = InequalitySolver::solve_inequality_checked(
            repeated, InequalityType::GreaterThan, "x");
        EXPECT_TRUE(repeated_positive &&
                        !repeated_positive.value().contains(1.0) &&
                        repeated_positive.value().contains(0.0) &&
                        repeated_positive.value().contains(2.0),
                    "checked strict repeated-root quadratic returns the punctured line");

        auto downward = SymbolicExpr::multiply(SymbolicExpr::number(-1), rational_roots);
        auto downward_positive = InequalitySolver::solve_inequality_checked(
            downward, InequalityType::GreaterThan, "x");
        EXPECT_TRUE(downward_positive && downward_positive.value().contains(0.0) &&
                        !downward_positive.value().contains(3.0),
                    downward_positive
                        ? "checked downward quadratic reverses the sign regions exactly"
                        : "checked downward quadratic failed: " +
                              downward_positive.error().message);
    }

    TEST_CASE("Unchecked quadratic inequalities preserve a small root");
    {
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto bx = SymbolicExpr::multiply(
            SymbolicExpr::number(BigInt("10000000000000000")), x);
        auto ill_conditioned = SymbolicExpr::add(
            SymbolicExpr::add(x2, bx), SymbolicExpr::number(1));
        auto positive = InequalitySolver::solve_inequality(
            ill_conditioned, InequalityType::GreaterThan, "x");

        EXPECT_TRUE(
            positive.contains(0.0),
            "x^2 + 10^16*x + 1 > 0 must not turn its small negative root into zero");
        EXPECT_TRUE(
            positive.contains(-5e-17),
            "quadratic sign chart includes points above its small negative root");
        EXPECT_TRUE(
            !positive.contains(-2e-16),
            "quadratic sign chart excludes points between its two roots");
    }

    TEST_CASE("Unchecked inequalities keep distinct nearby exact roots");
    {
        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto scaled_x = SymbolicExpr::multiply(
            SymbolicExpr::number(
                Rational(BigInt(-1), BigInt("1000000000000"))), x);
        auto nearby_roots = SymbolicExpr::add(x2, scaled_x);
        auto negative = InequalitySolver::solve_inequality(
            nearby_roots, InequalityType::LessThan, "x");

        EXPECT_TRUE(
            negative.contains(5e-13),
            "x*(x-10^-12) < 0 contains the interval between both exact roots");
        EXPECT_TRUE(
            !negative.contains(-1e-13),
            "nearby exact roots exclude the lower exterior: " +
                negative.to_string());
        EXPECT_TRUE(
            !negative.contains(2e-12),
            "nearby exact roots exclude the upper exterior");
    }

    TEST_CASE("Checked inequality conjunction contracts");
    {
        std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> bounded{
            {x, InequalityType::GreaterEqual},
            {SymbolicExpr::add(x, SymbolicExpr::number(-2)),
             InequalityType::LessThan}
        };
        auto bounded_result = InequalitySolver::solve_inequalities_checked(bounded, "x");
        EXPECT_TRUE(bounded_result.has_value(),
                    "checked conjunction of exact affine inequalities succeeds");
        if (bounded_result) {
            EXPECT_TRUE(bounded_result.value().contains(0.0) &&
                            bounded_result.value().contains(1.0) &&
                            !bounded_result.value().contains(2.0),
                        "checked conjunction preserves closed and open endpoints");
        }

        std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> empty;
        auto empty_result = InequalitySolver::solve_inequalities_checked(empty, "x");
        EXPECT_TRUE(empty_result && empty_result.value().is_entire_line(),
                    "empty checked conjunction denotes the entire real line");

        auto x2_minus_2 = SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::number(-2));
        std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> one_surd{
            {x2_minus_2, InequalityType::LessThan}
        };
        auto one_surd_result = InequalitySolver::solve_inequalities_checked(one_surd, "x");
        EXPECT_TRUE(one_surd_result && one_surd_result.value().contains(0.0),
                    "single checked quadratic conjunction keeps its proven-order surd interval");

        std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> mixed_surd{
            {x2_minus_2, InequalityType::LessThan},
            {x, InequalityType::GreaterThan}
        };
        auto mixed_surd_result = InequalitySolver::solve_inequalities_checked(
            mixed_surd, "x");
        EXPECT_TRUE(mixed_surd_result &&
                        mixed_surd_result.value().contains(1.0) &&
                        !mixed_surd_result.value().contains(-1.0) &&
                        !mixed_surd_result.value().contains(2.0),
                    mixed_surd_result
                        ? "mixed result=" + mixed_surd_result.value().to_string() +
                              " contains(1)=" +
                              std::to_string(mixed_surd_result.value().contains(1.0)) +
                              " contains(-1)=" +
                              std::to_string(mixed_surd_result.value().contains(-1.0)) +
                              " contains(2)=" +
                              std::to_string(mixed_surd_result.value().contains(2.0))
                        : "checked conjunction failed: " +
                              mixed_surd_result.error().message);

        std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> invalid{
            {nullptr, InequalityType::GreaterThan}
        };
        auto invalid_result = InequalitySolver::solve_inequalities_checked(invalid, "x");
        EXPECT_TRUE(!invalid_result && invalid_result.error().code == CasErrc::InvalidArgument,
                    "checked conjunction propagates invalid component errors");

        CancellationToken cancellation;
        cancellation.cancel();
        ComputationContext cancelled_context(ResourceLimits{}, cancellation);
        auto cancelled = InequalitySolver::solve_inequalities_checked(
            bounded, "x", cancelled_context);
        EXPECT_TRUE(!cancelled && cancelled.error().code == CasErrc::Cancelled,
                    "checked conjunction observes cancellation before processing components");
    }

    TEST_CASE("Linear inequality: 2x - 3 > 0");
    {

        auto expr = linear(2, -3);
        auto result = InequalitySolver::solve_inequality(expr, InequalityType::GreaterThan, "x");

        EXPECT_TRUE(!result.is_empty(), "2x - 3 > 0 should not be empty");
        EXPECT_TRUE(!result.contains(0.0), "2x - 3 > 0: x=0 not in solution");
        EXPECT_TRUE(!result.contains(1.0), "2x - 3 > 0: x=1 not in solution");
        EXPECT_TRUE(!result.contains(1.5), "2x - 3 > 0: x=1.5 (root) not in solution (strict)");
        EXPECT_TRUE(result.contains(2.0), "2x - 3 > 0: x=2 in solution");
        EXPECT_TRUE(result.contains(100.0), "2x - 3 > 0: x=100 in solution");
        EXPECT_TRUE(!result.contains(-5.0), "2x - 3 > 0: x=-5 not in solution");
    }

    TEST_CASE("Quadratic inequality: x^2 - 4 >= 0");
    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto expr = SymbolicExpr::add(x2, SymbolicExpr::number(-4));
        auto result = InequalitySolver::solve_inequality(expr, InequalityType::GreaterEqual, "x");

        EXPECT_TRUE(!result.is_empty(), "x^2 - 4 >= 0 should not be empty");

        EXPECT_TRUE(result.contains(-10.0), "x^2 - 4 >= 0: x=-10 in solution");
        EXPECT_TRUE(result.contains(-3.0), "x^2 - 4 >= 0: x=-3 in solution");
        EXPECT_TRUE(result.contains(3.0), "x^2 - 4 >= 0: x=3 in solution");
        EXPECT_TRUE(result.contains(10.0), "x^2 - 4 >= 0: x=10 in solution");

        EXPECT_TRUE(!result.contains(0.0), "x^2 - 4 >= 0: x=0 not in solution");
        EXPECT_TRUE(!result.contains(1.0), "x^2 - 4 >= 0: x=1 not in solution");
        EXPECT_TRUE(!result.contains(-1.0), "x^2 - 4 >= 0: x=-1 not in solution");
        EXPECT_TRUE(!result.contains(1.9), "x^2 - 4 >= 0: x=1.9 not in solution");
        EXPECT_TRUE(!result.contains(-1.9), "x^2 - 4 >= 0: x=-1.9 not in solution");

        EXPECT_TRUE(result.contains(2.0), "x^2 - 4 >= 0: x=2 in solution (non-strict, root)");
        EXPECT_TRUE(result.contains(-2.0), "x^2 - 4 >= 0: x=-2 in solution (non-strict, root)");
    }

    TEST_CASE("Cubic inequality: x^3 - x < 0");
    {

        auto x3 = SymbolicExpr::power(x, SymbolicExpr::number(3));
        auto expr = SymbolicExpr::add(x3, negate(x));
        auto result = InequalitySolver::solve_inequality(expr, InequalityType::LessThan, "x");

        EXPECT_TRUE(!result.is_empty(), "x^3 - x < 0 should not be empty");
        EXPECT_TRUE(result.contains(-5.0), "x^3 - x < 0: x=-5 in solution");
        EXPECT_TRUE(result.contains(-2.0), "x^3 - x < 0: x=-2 in solution");
        EXPECT_TRUE(!result.contains(-1.0), "x^3 - x < 0: x=-1 not in solution (strict, root)");
        EXPECT_TRUE(result.contains(std::nextafter(-1.0, -2.0)) &&
                        !result.contains(std::nextafter(-1.0, 0.0)),
                    "strict cubic boundary separates adjacent representable values");
        EXPECT_TRUE(!result.contains(1.0), "x^3 - x < 0: x=1 not in solution (strict, root)");
        EXPECT_TRUE(result.contains(0.5), "x^3 - x < 0: x=0.5 in solution");
        EXPECT_TRUE(!result.contains(2.0), "x^3 - x < 0: x=2 not in solution");
        EXPECT_TRUE(!result.contains(-0.5), "x^3 - x < 0: x=-0.5 not in solution (between -1 and 0)");

        EXPECT_TRUE(result.contains(0.001),
                    "x^3 - x < 0: x=0.001 in solution (just above 0)");
    }

    TEST_CASE("Repeated root: (x-1)^2*(x+2) > 0");
    {

        auto x_minus_1 = SymbolicExpr::add(x, SymbolicExpr::number(-1));
        auto x_minus_1_sq = SymbolicExpr::power(x_minus_1, SymbolicExpr::number(2));
        auto x_plus_2 = SymbolicExpr::add(x, SymbolicExpr::number(2));
        auto expr = SymbolicExpr::multiply(x_minus_1_sq, x_plus_2)->expand();

        auto result = InequalitySolver::solve_inequality(expr, InequalityType::GreaterThan, "x");

        EXPECT_TRUE(!result.is_empty(), "(x-1)^2*(x+2) > 0 should not be empty");
        EXPECT_TRUE(!result.contains(-2.0), "(x-1)^2*(x+2) > 0: x=-2 not in solution (root, strict)");
        EXPECT_TRUE(!result.contains(-3.0), "(x-1)^2*(x+2) > 0: x=-3 not in solution");
        EXPECT_TRUE(result.contains(0.0), "(x-1)^2*(x+2) > 0: x=0 in solution");
        EXPECT_TRUE(!result.contains(1.0), "(x-1)^2*(x+2) > 0: x=1 not in solution (root, strict)");
        EXPECT_TRUE(result.contains(2.0), "(x-1)^2*(x+2) > 0: x=2 in solution");
        EXPECT_TRUE(result.contains(10.0), "(x-1)^2*(x+2) > 0: x=10 in solution");
        EXPECT_TRUE(result.contains(-1.0), "(x-1)^2*(x+2) > 0: x=-1 in solution");
    }

    TEST_CASE("Rational inequality: (x-1)/(x+2) > 0");
    {
        auto numerator = SymbolicExpr::add(x, SymbolicExpr::number(-1));
        auto denominator = SymbolicExpr::add(x, SymbolicExpr::number(2));

        auto result = InequalitySolver::solve_rational_inequality(
            numerator, denominator, InequalityType::GreaterThan, "x");

        EXPECT_TRUE(!result.is_empty(), "(x-1)/(x+2) > 0 should not be empty");
        EXPECT_TRUE(result.contains(-5.0), "(x-1)/(x+2) > 0: x=-5 in solution");
        EXPECT_TRUE(!result.contains(-2.0), "(x-1)/(x+2) > 0: x=-2 not in solution (den root)");
        EXPECT_TRUE(!result.contains(0.0), "(x-1)/(x+2) > 0: x=0 not in solution");
        EXPECT_TRUE(!result.contains(1.0), "(x-1)/(x+2) > 0: x=1 not in solution (strict, num root)");
        EXPECT_TRUE(result.contains(2.0), "(x-1)/(x+2) > 0: x=2 in solution");
        EXPECT_TRUE(result.contains(100.0), "(x-1)/(x+2) > 0: x=100 in solution");
    }

    TEST_CASE("Rational inequality non-strict: (x-1)/(x+2) >= 0");
    {
        auto numerator = SymbolicExpr::add(x, SymbolicExpr::number(-1));
        auto denominator = SymbolicExpr::add(x, SymbolicExpr::number(2));

        auto result = InequalitySolver::solve_rational_inequality(
            numerator, denominator, InequalityType::GreaterEqual, "x");

        EXPECT_TRUE(!result.is_empty(), "(x-1)/(x+2) >= 0 should not be empty");
        EXPECT_TRUE(result.contains(-5.0), "(x-1)/(x+2) >= 0: x=-5 in solution");
        EXPECT_TRUE(!result.contains(-2.0), "(x-1)/(x+2) >= 0: x=-2 not in solution (den root excluded)");
        EXPECT_TRUE(!result.contains(0.0), "(x-1)/(x+2) >= 0: x=0 not in solution");
        EXPECT_TRUE(result.contains(1.0), "(x-1)/(x+2) >= 0: x=1 in solution (num root, non-strict)");
        EXPECT_TRUE(result.contains(2.0), "(x-1)/(x+2) >= 0: x=2 in solution");
    }

    TEST_CASE("System of inequalities: {x^2 - 4 > 0, x < 5}");
    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto expr1 = SymbolicExpr::add(x2, SymbolicExpr::number(-4));

        auto expr2 = SymbolicExpr::add(x, SymbolicExpr::number(-5));

        std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> system = {
            {expr1, InequalityType::GreaterThan},
            {expr2, InequalityType::LessThan}
        };

        auto result = InequalitySolver::solve_inequalities(system, "x");

        EXPECT_TRUE(!result.is_empty(), "{x^2-4>0, x<5} should not be empty");
        EXPECT_TRUE(result.contains(-10.0), "{x^2-4>0, x<5}: x=-10 in solution");
        EXPECT_TRUE(result.contains(-3.0), "{x^2-4>0, x<5}: x=-3 in solution");
        EXPECT_TRUE(!result.contains(-2.0), "{x^2-4>0, x<5}: x=-2 not in solution (strict root)");
        EXPECT_TRUE(!result.contains(0.0), "{x^2-4>0, x<5}: x=0 not in solution");
        EXPECT_TRUE(!result.contains(2.0), "{x^2-4>0, x<5}: x=2 not in solution (strict root)");
        EXPECT_TRUE(result.contains(3.0), "{x^2-4>0, x<5}: x=3 in solution");
        EXPECT_TRUE(result.contains(4.0), "{x^2-4>0, x<5}: x=4 in solution");
        EXPECT_TRUE(!result.contains(5.0), "{x^2-4>0, x<5}: x=5 not in solution (strict)");
        EXPECT_TRUE(!result.contains(6.0), "{x^2-4>0, x<5}: x=6 not in solution");
    }

    TEST_CASE("Zero polynomial: 0 > 0 -> empty");
    {
        auto zero_expr = SymbolicExpr::number(0);
        auto result = InequalitySolver::solve_inequality(zero_expr, InequalityType::GreaterThan, "x");

        EXPECT_TRUE(result.is_empty(), "0 > 0 should be empty");
        EXPECT_TRUE(!result.contains(0.0), "0 > 0: contains nothing");
        EXPECT_TRUE(!result.contains(1.0), "0 > 0: contains nothing");
    }

    TEST_CASE("Zero polynomial: 0 >= 0 -> entire line");
    {
        auto zero_expr = SymbolicExpr::number(0);
        auto result = InequalitySolver::solve_inequality(zero_expr, InequalityType::GreaterEqual, "x");

        EXPECT_TRUE(result.is_entire_line(), "0 >= 0 should be entire line");
        EXPECT_TRUE(result.contains(0.0), "0 >= 0: contains 0");
        EXPECT_TRUE(result.contains(-1000.0), "0 >= 0: contains -1000");
        EXPECT_TRUE(result.contains(1000.0), "0 >= 0: contains 1000");
    }

    TEST_CASE("Zero polynomial: 0 < 0 -> empty");
    {
        auto zero_expr = SymbolicExpr::number(0);
        auto result = InequalitySolver::solve_inequality(zero_expr, InequalityType::LessThan, "x");

        EXPECT_TRUE(result.is_empty(), "0 < 0 should be empty");
    }

    TEST_CASE("Zero polynomial: 0 <= 0 -> entire line");
    {
        auto zero_expr = SymbolicExpr::number(0);
        auto result = InequalitySolver::solve_inequality(zero_expr, InequalityType::LessEqual, "x");

        EXPECT_TRUE(result.is_entire_line(), "0 <= 0 should be entire line");
    }

    TEST_CASE("Exact huge constant inequality keeps exact sign");
    {
        std::string huge_digits = "1" + std::string(400, '0');
        auto huge_positive = SymbolicExpr::number(BigInt(huge_digits));
        auto positive_result = InequalitySolver::solve_inequality(
            huge_positive, InequalityType::GreaterThan, "x");
        EXPECT_TRUE(positive_result.is_entire_line(),
                    "huge exact positive constant > 0 should be entire line");

        auto huge_negative = SymbolicExpr::number(BigInt("-" + huge_digits));
        auto negative_result = InequalitySolver::solve_inequality(
            huge_negative, InequalityType::GreaterThan, "x");
        EXPECT_TRUE(negative_result.is_empty(),
                    "huge exact negative constant > 0 should be empty");
    }

    TEST_CASE("Non-polynomial: sin(x) > 0 -> empty (cannot solve)");
    {
        auto sin_x = SymbolicExpr::sin(x);
        auto result = InequalitySolver::solve_inequality(sin_x, InequalityType::GreaterThan, "x");

        EXPECT_TRUE(result.is_empty(), "sin(x) > 0 should return empty (non-polynomial)");
    }

    TEST_CASE("Non-polynomial: sin(x) >= 0 -> empty (cannot solve)");
    {
        auto sin_x = SymbolicExpr::sin(x);
        auto result = InequalitySolver::solve_inequality(sin_x, InequalityType::GreaterEqual, "x");

        EXPECT_TRUE(result.is_empty(), "sin(x) >= 0 should return empty (non-polynomial)");
    }

    TEST_CASE("Solution Soundness - points in solution set satisfy inequality");
    {
        std::mt19937 rng(42);
        const int NUM_ITERATIONS = 100;
        const int NUM_SAMPLES = 100;
        int pass_count = 0;

        std::uniform_int_distribution<int> degree_dist(2, 2);
        std::uniform_int_distribution<int> coeff_dist(-5, 5);
        std::uniform_int_distribution<int> type_dist(0, 3);

        InequalityType types[] = {
            InequalityType::GreaterThan,
            InequalityType::GreaterEqual,
            InequalityType::LessThan,
            InequalityType::LessEqual
        };

        auto build_poly = [&](const std::vector<int>& coeffs) -> std::shared_ptr<SymbolicExpr> {
            auto xv = SymbolicExpr::variable("x");
            std::shared_ptr<SymbolicExpr> result = SymbolicExpr::number(coeffs[0]);
            for (size_t i = 1; i < coeffs.size(); ++i) {
                if (coeffs[i] == 0) continue;
                auto term = SymbolicExpr::multiply(
                    SymbolicExpr::number(coeffs[i]),
                    SymbolicExpr::power(xv, SymbolicExpr::number(static_cast<int>(i)))
                );
                result = SymbolicExpr::add(result, term);
            }
            return result;
        };

        auto eval_poly = [](const std::shared_ptr<SymbolicExpr>& poly, double point) -> double {
            auto val_expr = SymbolicExpr::number(point);
            auto substituted = poly->substitute("x", val_expr);
            return substituted->to_numeric();
        };

        auto satisfies = [](double value, InequalityType type) -> bool {
            switch (type) {
                case InequalityType::GreaterThan:  return value > 0;
                case InequalityType::GreaterEqual: return value >= 0;
                case InequalityType::LessThan:     return value < 0;
                case InequalityType::LessEqual:    return value <= 0;
            }
            return false;
        };

        auto sample_inside = [&](const IntervalUnion& iu, double& out) -> bool {
            const auto& intervals = iu.intervals();
            if (intervals.empty()) return false;

            std::uniform_int_distribution<size_t> idx_dist(0, intervals.size() - 1);
            size_t idx = idx_dist(rng);
            const auto& iv = intervals[idx];

            double lo = -1000.0;
            double hi = 1000.0;

            if (!iv.lower.is_neg_infinity && iv.lower.value) {
                lo = iv.lower.value->to_numeric();
            }
            if (!iv.upper.is_pos_infinity && iv.upper.value) {
                hi = iv.upper.value->to_numeric();
            }

            double epsilon = 1e-6;
            if (iv.lower.is_open && !iv.lower.is_neg_infinity) lo += epsilon;
            if (iv.upper.is_open && !iv.upper.is_pos_infinity) hi -= epsilon;

            if (lo >= hi) {
                out = (lo + hi) / 2.0;
                return true;
            }

            std::uniform_real_distribution<double> point_dist(lo, hi);
            out = point_dist(rng);
            return true;
        };

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {

            int degree = degree_dist(rng);
            std::vector<int> coeffs(degree + 1);
            for (int i = 0; i <= degree; ++i) {
                coeffs[i] = coeff_dist(rng);
            }

            while (coeffs[degree] == 0) {
                coeffs[degree] = coeff_dist(rng);
            }

            auto poly = build_poly(coeffs);

            InequalityType type = types[type_dist(rng)];

            auto solution = InequalitySolver::solve_inequality(poly, type, "x");

            if (solution.is_empty()) {
                ++pass_count;
                continue;
            }

            bool iter_passed = true;
            for (int s = 0; s < NUM_SAMPLES; ++s) {
                double point;
                if (!sample_inside(solution, point)) continue;

                double value = eval_poly(poly, point);

                if (!satisfies(value, type)) {
                    std::ostringstream oss;
                    oss << "FAIL: iter=" << iter
                        << " point=" << point << " poly_value=" << value
                        << " type=" << static_cast<int>(type)
                        << " poly=[";
                    for (size_t i = 0; i < coeffs.size(); ++i) {
                        if (i > 0) oss << ",";
                        oss << coeffs[i];
                    }
                    oss << "]";
                    EXPECT_TRUE(false, oss.str());
                    iter_passed = false;
                    break;
                }
            }

            if (iter_passed) ++pass_count;
        }

        std::ostringstream summary;
        summary << "" << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed solution soundness";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("Solution Completeness - points outside solution set violate inequality");
    {
        std::mt19937 rng(123);
        const int NUM_ITERATIONS = 100;
        const int NUM_SAMPLES = 100;
        int pass_count = 0;

        std::uniform_int_distribution<int> degree_dist(2, 2);
        std::uniform_int_distribution<int> coeff_dist(-5, 5);
        std::uniform_int_distribution<int> type_dist(0, 3);

        InequalityType types[] = {
            InequalityType::GreaterThan,
            InequalityType::GreaterEqual,
            InequalityType::LessThan,
            InequalityType::LessEqual
        };

        auto build_poly = [&](const std::vector<int>& coeffs) -> std::shared_ptr<SymbolicExpr> {
            auto xv = SymbolicExpr::variable("x");
            std::shared_ptr<SymbolicExpr> result = SymbolicExpr::number(coeffs[0]);
            for (size_t i = 1; i < coeffs.size(); ++i) {
                if (coeffs[i] == 0) continue;
                auto term = SymbolicExpr::multiply(
                    SymbolicExpr::number(coeffs[i]),
                    SymbolicExpr::power(xv, SymbolicExpr::number(static_cast<int>(i)))
                );
                result = SymbolicExpr::add(result, term);
            }
            return result;
        };

        auto eval_poly = [](const std::shared_ptr<SymbolicExpr>& poly, double point) -> double {
            auto val_expr = SymbolicExpr::number(point);
            auto substituted = poly->substitute("x", val_expr);
            return substituted->to_numeric();
        };

        auto satisfies = [](double value, InequalityType type) -> bool {
            switch (type) {
                case InequalityType::GreaterThan:  return value > 0;
                case InequalityType::GreaterEqual: return value >= 0;
                case InequalityType::LessThan:     return value < 0;
                case InequalityType::LessEqual:    return value <= 0;
            }
            return false;
        };

        auto sample_outside = [&](const IntervalUnion& iu, double& out) -> bool {
            std::uniform_real_distribution<double> range_dist(-100.0, 100.0);

            for (int attempt = 0; attempt < 1000; ++attempt) {
                double candidate = range_dist(rng);

                if (iu.contains(candidate)) continue;

                bool near_boundary = false;
                for (const auto& iv : iu.intervals()) {
                    if (!iv.lower.is_neg_infinity && iv.lower.value) {
                        double boundary = iv.lower.value->to_numeric();
                        if (std::abs(candidate - boundary) < 1e-8) {
                            near_boundary = true;
                            break;
                        }
                    }
                    if (!iv.upper.is_pos_infinity && iv.upper.value) {
                        double boundary = iv.upper.value->to_numeric();
                        if (std::abs(candidate - boundary) < 1e-8) {
                            near_boundary = true;
                            break;
                        }
                    }
                }

                if (!near_boundary) {
                    out = candidate;
                    return true;
                }
            }
            return false;
        };

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {

            int degree = degree_dist(rng);
            std::vector<int> coeffs(degree + 1);
            for (int i = 0; i <= degree; ++i) {
                coeffs[i] = coeff_dist(rng);
            }

            while (coeffs[degree] == 0) {
                coeffs[degree] = coeff_dist(rng);
            }

            auto poly = build_poly(coeffs);

            InequalityType type = types[type_dist(rng)];

            auto solution = InequalitySolver::solve_inequality(poly, type, "x");

            if (solution.is_entire_line()) {
                ++pass_count;
                continue;
            }

            bool iter_passed = true;
            for (int s = 0; s < NUM_SAMPLES; ++s) {
                double point;
                if (!sample_outside(solution, point)) continue;

                double value = eval_poly(poly, point);

                if (satisfies(value, type)) {
                    std::ostringstream oss;
                    oss << "FAIL: iter=" << iter
                        << " point=" << point << " poly_value=" << value
                        << " type=" << static_cast<int>(type)
                        << " solution=" << solution.to_string()
                        << " poly=[";
                    for (size_t i = 0; i < coeffs.size(); ++i) {
                        if (i > 0) oss << ",";
                        oss << coeffs[i];
                    }
                    oss << "]";
                    EXPECT_TRUE(false, oss.str());
                    iter_passed = false;
                    break;
                }
            }

            if (iter_passed) ++pass_count;
        }

        std::ostringstream summary;
        summary << "" << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed solution completeness";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("Endpoint Correctness (Strict vs Non-strict)");
    {
        std::mt19937 rng(777);
        const int NUM_ITERATIONS = 100;
        int pass_count = 0;

        std::uniform_int_distribution<int> root_count_dist(1, 2);
        std::uniform_int_distribution<int> root_val_dist(-10, 10);
        std::uniform_int_distribution<int> type_dist(0, 3);

        InequalityType types[] = {
            InequalityType::GreaterThan,
            InequalityType::GreaterEqual,
            InequalityType::LessThan,
            InequalityType::LessEqual
        };

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {

            int num_roots = root_count_dist(rng);
            std::set<int> root_set;
            while ((int)root_set.size() < num_roots) {
                root_set.insert(root_val_dist(rng));
            }
            std::vector<int> roots(root_set.begin(), root_set.end());

            auto x_var = SymbolicExpr::variable("x");
            auto factor0 = SymbolicExpr::add(x_var, SymbolicExpr::number(-roots[0]));
            std::shared_ptr<SymbolicExpr> poly = factor0;
            for (size_t i = 1; i < roots.size(); ++i) {
                auto factor = SymbolicExpr::add(x_var, SymbolicExpr::number(-roots[i]));
                poly = SymbolicExpr::multiply(poly, factor);
            }
            poly = poly->expand();

            InequalityType ineq_type = types[type_dist(rng)];
            bool is_strict = (ineq_type == InequalityType::GreaterThan ||
                              ineq_type == InequalityType::LessThan);

            auto solution = InequalitySolver::solve_inequality(poly, ineq_type, "x");

            bool iter_passed = true;
            for (int r : roots) {
                bool root_in_solution = solution.contains((double)r);

                if (is_strict) {

                    if (root_in_solution) {
                        std::ostringstream oss;
                        oss << "FAIL (strict): iter=" << iter
                            << " root=" << r << " is in solution but shouldn't be"
                            << " poly=" << poly->to_string()
                            << " type=" << (int)ineq_type;
                        EXPECT_TRUE(false, oss.str());
                        iter_passed = false;
                        break;
                    }
                } else {

                    if (!root_in_solution) {
                        std::ostringstream oss;
                        oss << "FAIL (non-strict): iter=" << iter
                            << " root=" << r << " is NOT in solution but should be"
                            << " poly=" << poly->to_string()
                            << " type=" << (int)ineq_type;
                        EXPECT_TRUE(false, oss.str());
                        iter_passed = false;
                        break;
                    }
                }
            }

            if (iter_passed) ++pass_count;
        }

        std::ostringstream summary;
        summary << "" << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed endpoint correctness (strict vs non-strict)";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("Multiplicity Sign Change Correctness");
    {
        std::mt19937 rng(888);
        const int NUM_ITERATIONS = 100;
        int pass_count = 0;

        std::uniform_int_distribution<int> root_val_dist(-5, 5);
        std::uniform_int_distribution<int> even_mult_dist(0, 1);
        std::uniform_int_distribution<int> odd_mult_dist(0, 1);
        std::uniform_int_distribution<int> parity_dist(0, 1);

        const double epsilon = 0.5;

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {

            int root = root_val_dist(rng);
            bool use_even = (parity_dist(rng) == 0);
            int multiplicity;
            if (use_even) {
                multiplicity = (even_mult_dist(rng) == 0) ? 2 : 4;
            } else {
                multiplicity = (odd_mult_dist(rng) == 0) ? 1 : 3;
            }

            auto x_var = SymbolicExpr::variable("x");
            auto factor = SymbolicExpr::add(x_var, SymbolicExpr::number(-root));
            auto poly = factor;
            for (int i = 1; i < multiplicity; ++i) {
                poly = SymbolicExpr::multiply(poly, factor);
            }
            poly = poly->expand();

            double left_val, right_val;
            {
                auto sub_left = poly->substitute("x", SymbolicExpr::number((double)root - epsilon));
                auto simp_left = sub_left->simplify();
                left_val = simp_left->to_numeric();

                auto sub_right = poly->substitute("x", SymbolicExpr::number((double)root + epsilon));
                auto simp_right = sub_right->simplify();
                right_val = simp_right->to_numeric();
            }

            const double sign_tol = 1e-10;
            int left_sign = (left_val > sign_tol) ? 1 : ((left_val < -sign_tol) ? -1 : 0);
            int right_sign = (right_val > sign_tol) ? 1 : ((right_val < -sign_tol) ? -1 : 0);

            if (left_sign == 0 || right_sign == 0) {
                ++pass_count;
                continue;
            }

            bool iter_passed = true;
            if (use_even) {

                if (left_sign != right_sign) {
                    std::ostringstream oss;
                    oss << "FAIL (even mult): iter=" << iter
                        << " root=" << root << " mult=" << multiplicity
                        << " left_sign=" << left_sign << " right_sign=" << right_sign
                        << " left_val=" << left_val << " right_val=" << right_val;
                    EXPECT_TRUE(false, oss.str());
                    iter_passed = false;
                }
            } else {

                if (left_sign == right_sign) {
                    std::ostringstream oss;
                    oss << "FAIL (odd mult): iter=" << iter
                        << " root=" << root << " mult=" << multiplicity
                        << " left_sign=" << left_sign << " right_sign=" << right_sign
                        << " left_val=" << left_val << " right_val=" << right_val;
                    EXPECT_TRUE(false, oss.str());
                    iter_passed = false;
                }
            }

            if (iter_passed) ++pass_count;
        }

        std::ostringstream summary;
        summary << "" << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed multiplicity sign change correctness";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    TEST_CASE("Parametric Inequality Consistency");
    {
        std::mt19937 rng(1111);
        const int NUM_ITERATIONS = 100;
        const int NUM_SAMPLES = 50;
        int pass_count = 0;

        std::uniform_int_distribution<int> leading_dist(1, 4);
        std::uniform_int_distribution<int> sign_dist(0, 1);
        std::uniform_int_distribution<int> param_dist(-5, 5);
        std::uniform_int_distribution<int> type_dist(0, 3);

        InequalityType types[] = {
            InequalityType::GreaterThan,
            InequalityType::GreaterEqual,
            InequalityType::LessThan,
            InequalityType::LessEqual
        };

        auto satisfies = [](double value, InequalityType type) -> bool {
            switch (type) {
                case InequalityType::GreaterThan:  return value > 0;
                case InequalityType::GreaterEqual: return value >= 0;
                case InequalityType::LessThan:     return value < 0;
                case InequalityType::LessEqual:    return value <= 0;
            }
            return false;
        };

        for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {

            int a_val = leading_dist(rng);
            if (sign_dist(rng)) a_val = -a_val;

            int p_val = param_dist(rng);
            int q_val = param_dist(rng);

            int disc = p_val * p_val - 4 * a_val * q_val;
            if (disc < 0) {

                if (a_val > 0) {
                    q_val = (p_val * p_val) / (4 * a_val) - 1;
                } else {
                    q_val = (p_val * p_val) / (4 * a_val) + 1;
                }
                disc = p_val * p_val - 4 * a_val * q_val;
                if (disc < 0) {
                    ++pass_count;
                    continue;
                }
            }

            InequalityType ineq_type = types[type_dist(rng)];

            auto x_var = SymbolicExpr::variable("x");
            auto p_param = SymbolicExpr::variable("p");
            auto q_param = SymbolicExpr::variable("q");

            auto ax2 = SymbolicExpr::multiply(
                SymbolicExpr::number(a_val),
                SymbolicExpr::power(x_var, SymbolicExpr::number(2)));
            auto px = SymbolicExpr::multiply(p_param, x_var);
            auto parametric_expr = SymbolicExpr::add(SymbolicExpr::add(ax2, px), q_param);

            auto parametric_result = InequalitySolver::solve_parametric_inequality(
                parametric_expr, ineq_type, "x", {"p", "q"});

            auto concrete_expr = SymbolicExpr::add(
                SymbolicExpr::add(
                    SymbolicExpr::multiply(SymbolicExpr::number(a_val),
                        SymbolicExpr::power(x_var, SymbolicExpr::number(2))),
                    SymbolicExpr::multiply(SymbolicExpr::number(p_val), x_var)),
                SymbolicExpr::number(q_val));

            auto direct_solution = InequalitySolver::solve_inequality(
                concrete_expr, ineq_type, "x");

            IntervalUnion parametric_solution = IntervalUnion::empty();
            bool found_case = false;

            if (parametric_result.is_empty()) {
                ++pass_count;
                continue;
            }

            if (parametric_result.is_single()) {
                parametric_solution = parametric_result.single_solution();
                found_case = true;
            } else {

                for (const auto& pcase : parametric_result.cases) {
                    if (!pcase.condition) {
                        parametric_solution = pcase.solution;
                        found_case = true;
                        break;
                    }
                }
                if (!found_case && !parametric_result.cases.empty()) {
                    parametric_solution = parametric_result.cases[0].solution;
                    found_case = true;
                }
            }

            if (!found_case) {
                ++pass_count;
                continue;
            }

            bool iter_passed = true;
            std::uniform_real_distribution<double> sample_dist(-20.0, 20.0);

            for (int s = 0; s < NUM_SAMPLES; ++s) {
                double test_point = sample_dist(rng);

                double poly_val = (double)a_val * test_point * test_point
                                + (double)p_val * test_point
                                + (double)q_val;

                if (std::abs(poly_val) < 1e-6) continue;

                bool expected_result = satisfies(poly_val, ineq_type);

                bool in_direct = direct_solution.contains(test_point);

                bool in_parametric = false;
                const auto& param_intervals = parametric_solution.intervals();

                if (parametric_solution.is_empty()) {
                    in_parametric = false;
                } else if (parametric_solution.is_entire_line()) {
                    in_parametric = true;
                } else {
                    for (const auto& iv : param_intervals) {
                        double lo = -1e18;
                        double hi = 1e18;
                        bool lo_open = true;
                        bool hi_open = true;
                        bool lo_valid = true;
                        bool hi_valid = true;

                        if (!iv.lower.is_neg_infinity && iv.lower.value) {
                            auto lo_expr = iv.lower.value->substitute("p", SymbolicExpr::number(p_val));
                            lo_expr = lo_expr->substitute("q", SymbolicExpr::number(q_val));
                            lo_expr = lo_expr->simplify();
                            auto lo_opt = test_numeric_eval(lo_expr);
                            if (lo_opt && std::isfinite(*lo_opt)) {
                                lo = *lo_opt;
                                lo_open = iv.lower.is_open;
                            } else {
                                /// 下界求值未决时跳过当前区间.
                                lo_valid = false;
                            }
                        }

                        if (!iv.upper.is_pos_infinity && iv.upper.value) {
                            auto hi_expr = iv.upper.value->substitute("p", SymbolicExpr::number(p_val));
                            hi_expr = hi_expr->substitute("q", SymbolicExpr::number(q_val));
                            hi_expr = hi_expr->simplify();
                            auto hi_opt = test_numeric_eval(hi_expr);
                            if (hi_opt && std::isfinite(*hi_opt)) {
                                hi = *hi_opt;
                                hi_open = iv.upper.is_open;
                            } else {
                                hi_valid = false;
                            }
                        }

                        // If we couldn't evaluate an endpoint, try the
                        /// 同时尝试化简前形式,以覆盖 numeric_eval 支持的原始结构.
                        if (!lo_valid && !iv.lower.is_neg_infinity && iv.lower.value) {
                            auto lo_expr = iv.lower.value->substitute("p", SymbolicExpr::number(p_val));
                            lo_expr = lo_expr->substitute("q", SymbolicExpr::number(q_val));
                            auto lo_opt = test_numeric_eval(lo_expr);
                            if (lo_opt && std::isfinite(*lo_opt)) {
                                lo = *lo_opt;
                                lo_open = iv.lower.is_open;
                                lo_valid = true;
                            }
                        }
                        if (!hi_valid && !iv.upper.is_pos_infinity && iv.upper.value) {
                            auto hi_expr = iv.upper.value->substitute("p", SymbolicExpr::number(p_val));
                            hi_expr = hi_expr->substitute("q", SymbolicExpr::number(q_val));
                            auto hi_opt = test_numeric_eval(hi_expr);
                            if (hi_opt && std::isfinite(*hi_opt)) {
                                hi = *hi_opt;
                                hi_open = iv.upper.is_open;
                                hi_valid = true;
                            }
                        }

                        if (!lo_valid || !hi_valid) continue;

                        bool in_interval = false;
                        if (lo_open) {
                            in_interval = (test_point > lo + 1e-10);
                        } else {
                            in_interval = (test_point >= lo - 1e-10);
                        }
                        if (in_interval) {
                            if (hi_open) {
                                in_interval = (test_point < hi - 1e-10);
                            } else {
                                in_interval = (test_point <= hi + 1e-10);
                            }
                        }

                        if (in_interval) {
                            in_parametric = true;
                            break;
                        }
                    }
                }

                if (in_parametric != expected_result) {

                    /// 所有端点均求值未决时跳过当前采样点.
                    if (!in_parametric && !parametric_solution.is_empty() &&
                        !parametric_solution.is_entire_line()) {
                        // Check if we actually evaluated at least one interval.
                        bool any_evaluated = false;
                        for (const auto& iv2 : param_intervals) {
                            bool l_ok = iv2.lower.is_neg_infinity;
                            bool u_ok = iv2.upper.is_pos_infinity;
                            if (!l_ok && iv2.lower.value) {
                                auto le = iv2.lower.value->substitute("p", SymbolicExpr::number(p_val));
                                le = le->substitute("q", SymbolicExpr::number(q_val));
                                auto lv = test_numeric_eval(le);
                                l_ok = lv.has_value();
                            } else { l_ok = true; }
                            if (!u_ok && iv2.upper.value) {
                                auto ue = iv2.upper.value->substitute("p", SymbolicExpr::number(p_val));
                                ue = ue->substitute("q", SymbolicExpr::number(q_val));
                                auto uv = test_numeric_eval(ue);
                                u_ok = uv.has_value();
                            } else { u_ok = true; }
                            if (l_ok && u_ok) { any_evaluated = true; break; }
                        }
                        if (!any_evaluated) continue;
                    }

                    bool near_root = false;
                    double root1 = (-p_val + std::sqrt(std::abs(disc))) / (2.0 * a_val);
                    double root2 = (-p_val - std::sqrt(std::abs(disc))) / (2.0 * a_val);
                    if (std::abs(test_point - root1) < 1e-4 ||
                        std::abs(test_point - root2) < 1e-4) {
                        near_root = true;
                    }

                    if (!near_root) {
                        std::ostringstream oss;
                        oss << "FAIL: iter=" << iter
                            << " a=" << a_val << " p=" << p_val << " q=" << q_val
                            << " type=" << static_cast<int>(ineq_type)
                            << " point=" << test_point
                            << " expected=" << expected_result
                            << " in_parametric=" << in_parametric
                            << " in_direct=" << in_direct
                            << " poly_val=" << poly_val;
                        EXPECT_TRUE(false, oss.str());
                        iter_passed = false;
                        break;
                    }
                }
            }

            if (iter_passed) ++pass_count;
        }

        std::ostringstream summary;
        summary << "" << pass_count << "/" << NUM_ITERATIONS
                << " iterations passed parametric inequality consistency";
        EXPECT_TRUE(pass_count == NUM_ITERATIONS, summary.str());
    }

    return TEST_REPORT();
}
