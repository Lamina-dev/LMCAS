#include "test_common.hpp"
#include "interval.hpp"
#include "inequality_solver.hpp"
#include "symbolic.hpp"
#include "poly_utils.hpp"
#include <cmath>
#include <random>
#include <string>
#include <vector>

using namespace LMCAS;

TEST(ParametricInequality, ParametricLinearAxB0PiecewiseBySignOfA) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");

    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(a, x),
        b);

    auto result = InequalitySolver::solve_parametric_inequality(
        expr, InequalityType::GreaterThan, "x", {"a", "b"});

    EXPECT_TRUE((result.cases.size() >= 3)) << "ax + b > 0: should have at least 3 piecewise cases (a>0, a<0, then a=0 sub-cases by sign of b)";

    auto &pos_case = result.cases[0];
    EXPECT_TRUE((pos_case.condition != nullptr)) << "ax + b > 0, case a>0: condition should not be null";
    EXPECT_TRUE((!pos_case.solution.is_empty())) << "ax + b > 0, case a>0: solution should not be empty";
    EXPECT_TRUE((pos_case.solution.intervals().size() == 1)) << "ax + b > 0, case a>0: should have exactly 1 interval";

    if (!pos_case.solution.intervals().empty()) {
        auto &iv = pos_case.solution.intervals()[0];

        EXPECT_TRUE((!iv.lower.is_neg_infinity)) << "ax + b > 0, case a>0: lower bound should not be -inf";
        EXPECT_TRUE((iv.lower.is_open)) << "ax + b > 0, case a>0: lower bound should be open (strict)";
        EXPECT_TRUE((iv.upper.is_pos_infinity)) << "ax + b > 0, case a>0: upper bound should be +inf";

        if (iv.lower.value) {

            auto substituted = iv.lower.value->substitute("a", SymbolicExpr::number(2));
            substituted = substituted->substitute("b", SymbolicExpr::number(4));
            substituted = substituted->simplify();
            double val = 0.0;
            ASSERT_NO_THROW(val = substituted->to_numeric())
                << "ax + b > 0, case a>0: root should be evaluable";
            EXPECT_TRUE((std::abs(val - (-2.0)) < 1e-10)) << "ax + b > 0, case a>0: root at a=2,b=4 should be -2";
        }
    }

    auto &neg_case = result.cases[1];
    EXPECT_TRUE((neg_case.condition != nullptr)) << "ax + b > 0, case a<0: condition should not be null";
    EXPECT_TRUE((!neg_case.solution.is_empty())) << "ax + b > 0, case a<0: solution should not be empty";
    EXPECT_TRUE((neg_case.solution.intervals().size() == 1)) << "ax + b > 0, case a<0: should have exactly 1 interval";

    if (!neg_case.solution.intervals().empty()) {
        auto &iv = neg_case.solution.intervals()[0];
        EXPECT_TRUE((iv.lower.is_neg_infinity)) << "ax + b > 0, case a<0: lower bound should be -inf";
        EXPECT_TRUE((!iv.upper.is_pos_infinity)) << "ax + b > 0, case a<0: upper bound should not be +inf";
        EXPECT_TRUE((iv.upper.is_open)) << "ax + b > 0, case a<0: upper bound should be open (strict)";
    }

    auto &degen_case = result.cases[2];
    EXPECT_TRUE((degen_case.condition != nullptr)) << "ax + b > 0, case a=0: condition should not be null";
}

TEST(ParametricInequality, ParametricQuadraticWithConstantLeadingCoeffX2BxC0) {
    auto x = SymbolicExpr::variable("x");
    auto b = SymbolicExpr::variable("b");
    auto c = SymbolicExpr::variable("c");

    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto bx = SymbolicExpr::multiply(b, x);
    auto expr = SymbolicExpr::add(SymbolicExpr::add(x2, bx), c);

    auto result = InequalitySolver::solve_parametric_inequality(
        expr, InequalityType::GreaterThan, "x", {"b", "c"});

    EXPECT_TRUE((result.cases.size() == 1)) << "x^2 + bx + c > 0: should have 1 case (constant leading coeff)";

    auto &single_case = result.cases[0];
    EXPECT_TRUE((single_case.condition == nullptr)) << "x^2 + bx + c > 0: condition should be null (unconditional)";

    auto &intervals = single_case.solution.intervals();
    EXPECT_TRUE((intervals.size() == 2)) << "x^2 + bx + c > 0: should have 2 intervals (outside roots)";

    if (intervals.size() == 2) {

        EXPECT_TRUE((intervals[0].lower.is_neg_infinity)) << "x^2 + bx + c > 0: first interval starts at -inf";
        EXPECT_TRUE((intervals[0].upper.is_open)) << "x^2 + bx + c > 0: first interval upper is open (strict)";

        EXPECT_TRUE((intervals[1].upper.is_pos_infinity)) << "x^2 + bx + c > 0: second interval ends at +inf";
        EXPECT_TRUE((intervals[1].lower.is_open)) << "x^2 + bx + c > 0: second interval lower is open (strict)";

        if (intervals[0].upper.value) {
            auto root1 = intervals[0].upper.value->substitute("b", SymbolicExpr::number(-3));
            root1 = root1->substitute("c", SymbolicExpr::number(2));
            root1 = root1->simplify();
            double val = 0.0;
            ASSERT_NO_THROW(val = root1->to_numeric())
                << "x^2 + bx + c > 0: root should be evaluable";
            EXPECT_TRUE((std::abs(val - 1.0) < 1e-10 || std::abs(val - 2.0) < 1e-10)) << "x^2 + bx + c > 0: root1 at b=-3,c=2 should be 1 or 2";
        }
    }
}

TEST(ParametricInequality, ParametricInequalityKeepsExactHugeConstantLeadingSign) {
    auto x = SymbolicExpr::variable("x");

    std::string huge_digits = "1" + std::string(400, '0');
    auto huge_positive = SymbolicExpr::number(BigInt(huge_digits));
    auto positive_expr = SymbolicExpr::add(
        SymbolicExpr::multiply(huge_positive, x),
        SymbolicExpr::number(1));

    auto positive_result = InequalitySolver::solve_parametric_inequality(
        positive_expr, InequalityType::GreaterThan, "x", {"a"});

    EXPECT_TRUE((positive_result.cases.size() == 1)) << "huge positive leading coeff not depending on params should produce one case";
    if (!positive_result.cases.empty()) {
        auto &intervals = positive_result.cases[0].solution.intervals();
        EXPECT_TRUE((intervals.size() == 1)) << "huge positive leading coeff > 0 should produce one interval";
        if (!intervals.empty()) {
            EXPECT_TRUE((!intervals[0].lower.is_neg_infinity)) << "huge positive leading coeff: lower bound should be finite root";
            EXPECT_TRUE((intervals[0].upper.is_pos_infinity)) << "huge positive leading coeff: solution should extend to +inf";
        }
    }

    auto huge_negative = SymbolicExpr::number(BigInt("-" + huge_digits));
    auto negative_expr = SymbolicExpr::add(
        SymbolicExpr::multiply(huge_negative, x),
        SymbolicExpr::number(1));

    auto negative_result = InequalitySolver::solve_parametric_inequality(
        negative_expr, InequalityType::GreaterThan, "x", {"a"});

    EXPECT_TRUE((negative_result.cases.size() == 1)) << "huge negative leading coeff not depending on params should produce one case";
    if (!negative_result.cases.empty()) {
        auto &intervals = negative_result.cases[0].solution.intervals();
        EXPECT_TRUE((intervals.size() == 1)) << "huge negative leading coeff > 0 should produce one interval";
        if (!intervals.empty()) {
            EXPECT_TRUE((intervals[0].lower.is_neg_infinity)) << "huge negative leading coeff: solution should extend to -inf";
            EXPECT_TRUE((!intervals[0].upper.is_pos_infinity)) << "huge negative leading coeff: upper bound should be finite root";
        }
    }
}

TEST(ParametricInequality, ParametricQuadraticAx2BxC0PiecewiseBySignOfA) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");
    auto c = SymbolicExpr::variable("c");

    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto ax2 = SymbolicExpr::multiply(a, x2);
    auto bx = SymbolicExpr::multiply(b, x);
    auto expr = SymbolicExpr::add(SymbolicExpr::add(ax2, bx), c);

    auto result = InequalitySolver::solve_parametric_inequality(
        expr, InequalityType::GreaterEqual, "x", {"a", "b", "c"});

    EXPECT_TRUE((result.cases.size() >= 3)) << "ax^2 + bx + c >= 0: should have at least 3 piecewise cases (a>0, a<0, a=0 may further split)";

    auto &pos_case = result.cases[0];
    EXPECT_TRUE((pos_case.condition != nullptr)) << "ax^2 + bx + c >= 0, case a>0: should have condition";

    auto &neg_case = result.cases[1];
    EXPECT_TRUE((neg_case.condition != nullptr)) << "ax^2 + bx + c >= 0, case a<0: should have condition";

    // 第三个及以后的分支对应 a=0 的退化情形（可能进一步按 b、c 的符号细分）。
    if (result.cases.size() >= 3) {
        auto &degen_case = result.cases[2];
        EXPECT_TRUE((degen_case.condition != nullptr)) << "ax^2 + bx + c >= 0, case a=0: should have condition (degenerate)";
    }
}

TEST(ParametricInequality, ParametricConsistencyParametricVsNumericForSpecificValues) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");

    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(a, x),
        b);

    auto parametric_result = InequalitySolver::solve_parametric_inequality(
        expr, InequalityType::GreaterThan, "x", {"a", "b"});

    auto numeric_expr = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(3), x),
        SymbolicExpr::number(-6));
    auto numeric_result = InequalitySolver::solve_inequality(
        numeric_expr, InequalityType::GreaterThan, "x");

    EXPECT_TRUE((!numeric_result.is_empty())) << "3x - 6 > 0: numeric solution should not be empty";
    EXPECT_TRUE((numeric_result.contains(3.0))) << "3x - 6 > 0: x=3 should be in solution";
    EXPECT_TRUE((!numeric_result.contains(1.0))) << "3x - 6 > 0: x=1 should not be in solution";
    EXPECT_TRUE((!numeric_result.contains(2.0))) << "3x - 6 > 0: x=2 (root) should not be in solution (strict)";

    EXPECT_TRUE((!parametric_result.cases.empty())) << "Parametric ax + b > 0: should have cases";

    if (!parametric_result.cases.empty()) {
        auto &pos_case = parametric_result.cases[0];
        auto &intervals = pos_case.solution.intervals();
        if (!intervals.empty() && intervals[0].lower.value) {
            auto boundary = intervals[0].lower.value;
            boundary = boundary->substitute("a", SymbolicExpr::number(3));
            boundary = boundary->substitute("b", SymbolicExpr::number(-6));
            boundary = boundary->simplify();
            double val = 0.0;
            ASSERT_NO_THROW(val = boundary->to_numeric())
                << "Parametric consistency: boundary should be evaluable";
            EXPECT_TRUE((std::abs(val - 2.0) < 1e-10)) << "Parametric consistency: boundary at a=3,b=-6 should be 2";
        }
    }
}

TEST(ParametricInequality, ParametricWithEmptyParamsFallsBackToNonParametric) {
    auto x = SymbolicExpr::variable("x");

    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto expr = SymbolicExpr::add(x2, SymbolicExpr::number(-4));

    auto result = InequalitySolver::solve_parametric_inequality(
        expr, InequalityType::GreaterThan, "x", {});

    EXPECT_TRUE((result.cases.size() == 1)) << "No params: should have 1 case";
    EXPECT_TRUE((!result.cases[0].solution.is_empty())) << "No params: x^2 - 4 > 0 should not be empty";
    EXPECT_TRUE((result.cases[0].solution.contains(3.0))) << "No params: x=3 should be in solution";
    EXPECT_TRUE((result.cases[0].solution.contains(-3.0))) << "No params: x=-3 should be in solution";
    EXPECT_TRUE((!result.cases[0].solution.contains(0.0))) << "No params: x=0 should not be in solution";
}

TEST(ParametricInequality, DegenerateCaseAx22x10WhenA0) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");

    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto ax2 = SymbolicExpr::multiply(a, x2);
    auto two_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
    auto expr = SymbolicExpr::add(SymbolicExpr::add(ax2, two_x), SymbolicExpr::number(1));

    auto result = InequalitySolver::solve_parametric_inequality(
        expr, InequalityType::GreaterThan, "x", {"a"});

    EXPECT_TRUE((result.cases.size() == 3)) << "ax^2 + 2x + 1 > 0: should have 3 cases";

    if (result.cases.size() == 3) {
        auto &degen = result.cases[2];
        EXPECT_TRUE((!degen.solution.is_empty())) << "Degenerate (a=0): 2x + 1 > 0 should not be empty";

        auto &intervals = degen.solution.intervals();
        EXPECT_TRUE((intervals.size() == 1)) << "Degenerate (a=0): should have 1 interval";

        if (!intervals.empty()) {
            EXPECT_TRUE((intervals[0].upper.is_pos_infinity)) << "Degenerate (a=0): upper bound should be +inf";
            EXPECT_TRUE((intervals[0].lower.is_open)) << "Degenerate (a=0): lower bound should be open (strict)";

            if (intervals[0].lower.value) {
                auto val_expr = intervals[0].lower.value->simplify();
                try {
                    double val = val_expr->to_numeric();
                    EXPECT_TRUE((std::abs(val - (-0.5)) < 1e-10)) << "Degenerate (a=0): boundary should be -1/2";
                } catch (...) {
                    ADD_FAILURE() << "Degenerate (a=0): boundary should be evaluable";
                }
            }
        }
    }
}
