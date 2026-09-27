#include "test_common.hpp"
#include "expr.hpp"
#include "interval.hpp"
#include "inequality_solver.hpp"
#include "symbolic.hpp"
#include "poly_utils.hpp"
#include <cmath>
#include <random>
#include <string>
#include <vector>

using namespace LMCAS;

static void check_quadratic_parametric_boundary(
    const std::shared_ptr<SymbolicExpr> &boundary) {
    auto root = boundary->substitute("a", SymbolicExpr::number(1));
    root = root->substitute("b", SymbolicExpr::number(-5));
    root = root->substitute("c", SymbolicExpr::number(6));
    root = root->simplify();
    auto evaluated = evalf(*root);
    ASSERT_TRUE(evaluated.has_value())
        << "ax^2+bx+c>=0 boundary must be numerically evaluable";
    const double value = evaluated.value().value;
    EXPECT_TRUE(std::abs(value - 2.0) < 1e-10 ||
                std::abs(value - 3.0) < 1e-10)
        << "ax^2+bx+c>=0, a>0: root at a=1,b=-5,c=6 should be 2 or 3";
}

TEST(ParametricInequalityDegenerate, AxB0A0NumericVerification) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");

    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(a, x),
        b);

    auto result = InequalitySolver::solve_parametric_inequality(
        expr, InequalityType::GreaterThan, "x", {"a", "b"});

    EXPECT_TRUE((!result.cases.empty())) << "ax + b > 0: should have cases";
    EXPECT_TRUE((result.cases.size() >= 1)) << "ax + b > 0: should have at least 1 case";

    auto &pos_case = result.cases[0];
    auto &pos_intervals = pos_case.solution.intervals();

    if (!pos_intervals.empty() && pos_intervals[0].lower.value) {
        auto boundary = pos_intervals[0].lower.value;
        auto b1 = boundary->substitute("a", SymbolicExpr::number(5));
        b1 = b1->substitute("b", SymbolicExpr::number(-10));
        b1 = b1->simplify();
        try {
            double val = b1->to_numeric();
            EXPECT_TRUE((std::abs(val - 2.0) < 1e-10)) << "ax+b>0, a>0: a=5,b=-10 → boundary should be 2";
        } catch (...) {
            ADD_FAILURE() << "ax+b>0, a>0: boundary should be evaluable for a=5,b=-10";
        }

        auto b2 = boundary->substitute("a", SymbolicExpr::number(1));
        b2 = b2->substitute("b", SymbolicExpr::number(3));
        b2 = b2->simplify();
        try {
            double val = b2->to_numeric();
            EXPECT_TRUE((std::abs(val - (-3.0)) < 1e-10)) << "ax+b>0, a>0: a=1,b=3 → boundary should be -3";
        } catch (...) {
            ADD_FAILURE() << "ax+b>0, a>0: boundary should be evaluable for a=1,b=3";
        }
    }

    if (!pos_intervals.empty()) {
        EXPECT_TRUE((pos_intervals[0].upper.is_pos_infinity)) << "ax+b>0, a>0: upper bound should be +∞";
        EXPECT_TRUE((pos_intervals[0].lower.is_open)) << "ax+b>0, a>0: lower bound should be open (strict >)";
    }
}

TEST(ParametricInequalityDegenerate, AxB0A0NumericVerificationContract) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");

    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(a, x),
        b);

    auto result = InequalitySolver::solve_parametric_inequality(
        expr, InequalityType::GreaterThan, "x", {"a", "b"});

    EXPECT_TRUE((result.cases.size() >= 2)) << "ax + b > 0: should have at least 2 cases (a>0, a<0)";

    auto &neg_case = result.cases[1];
    auto &neg_intervals = neg_case.solution.intervals();

    EXPECT_TRUE((!neg_intervals.empty())) << "ax+b>0, a<0: solution should not be empty";

    if (!neg_intervals.empty()) {

        EXPECT_TRUE((neg_intervals[0].lower.is_neg_infinity)) << "ax+b>0, a<0: lower bound should be -∞";
        EXPECT_TRUE((neg_intervals[0].upper.is_open)) << "ax+b>0, a<0: upper bound should be open (strict >)";
        EXPECT_TRUE((!neg_intervals[0].upper.is_pos_infinity)) << "ax+b>0, a<0: upper bound should not be +∞";

        if (neg_intervals[0].upper.value) {
            auto boundary = neg_intervals[0].upper.value;
            auto b1 = boundary->substitute("a", SymbolicExpr::number(-2));
            b1 = b1->substitute("b", SymbolicExpr::number(6));
            b1 = b1->simplify();
            try {
                double val = b1->to_numeric();
                EXPECT_TRUE((std::abs(val - 3.0) < 1e-10)) << "ax+b>0, a<0: a=-2,b=6 → boundary should be 3";
            } catch (...) {
                ADD_FAILURE() << "ax+b>0, a<0: boundary should be evaluable for a=-2,b=6";
            }

            auto b2 = boundary->substitute("a", SymbolicExpr::number(-4));
            b2 = b2->substitute("b", SymbolicExpr::number(-8));
            b2 = b2->simplify();
            try {
                double val = b2->to_numeric();
                EXPECT_TRUE((std::abs(val - (-2.0)) < 1e-10)) << "ax+b>0, a<0: a=-4,b=-8 → boundary should be -2";
            } catch (...) {
                ADD_FAILURE() << "ax+b>0, a<0: boundary should be evaluable for a=-4,b=-8";
            }
        }
    }
}

TEST(ParametricInequalityDegenerate, Ax2BxC0DiscriminantDependsOnParameters) {
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

    EXPECT_TRUE((result.cases.size() >= 3)) << "ax^2+bx+c>=0: should have at least 3 piecewise cases (degenerate a=0 may split further)";

    if (result.cases.size() >= 1) {
        auto &pos_case = result.cases[0];

        auto &intervals = pos_case.solution.intervals();

        if (intervals.size() == 2) {

            if (intervals[0].upper.value) {
                check_quadratic_parametric_boundary(intervals[0].upper.value);
            }

            if (intervals[1].lower.value) {
                check_quadratic_parametric_boundary(intervals[1].lower.value);
            }

            EXPECT_TRUE((!intervals[0].upper.is_open)) << "ax^2+bx+c>=0, a>0: upper of first interval should be closed (non-strict)";
            EXPECT_TRUE((!intervals[1].lower.is_open)) << "ax^2+bx+c>=0, a>0: lower of second interval should be closed (non-strict)";
        }
    }

    if (result.cases.size() >= 2) {
        auto &neg_case = result.cases[1];

        auto &intervals = neg_case.solution.intervals();

        if (intervals.size() == 1) {

            EXPECT_TRUE((!intervals[0].lower.is_neg_infinity)) << "ax^2+bx+c>=0, a<0: lower should not be -∞";
            EXPECT_TRUE((!intervals[0].upper.is_pos_infinity)) << "ax^2+bx+c>=0, a<0: upper should not be +∞";
            EXPECT_TRUE((!intervals[0].lower.is_open)) << "ax^2+bx+c>=0, a<0: lower should be closed (non-strict)";
            EXPECT_TRUE((!intervals[0].upper.is_open)) << "ax^2+bx+c>=0, a<0: upper should be closed (non-strict)";
        }
    }
}

TEST(ParametricInequalityDegenerate, DegenerateCaseAx23x20WhenA0) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");

    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto ax2 = SymbolicExpr::multiply(a, x2);
    auto three_x = SymbolicExpr::multiply(SymbolicExpr::number(3), x);
    auto expr = SymbolicExpr::add(
        SymbolicExpr::add(ax2, three_x),
        SymbolicExpr::number(-2));

    auto result = InequalitySolver::solve_parametric_inequality(
        expr, InequalityType::GreaterThan, "x", {"a"});

    EXPECT_TRUE((result.cases.size() == 3)) << "ax^2+3x-2>0: should have 3 cases";

    if (result.cases.size() == 3) {

        auto &degen = result.cases[2];
        EXPECT_TRUE((!degen.solution.is_empty())) << "Degenerate (a=0): 3x-2>0 should not be empty";

        auto &intervals = degen.solution.intervals();
        EXPECT_TRUE((intervals.size() == 1)) << "Degenerate (a=0): should have 1 interval";

        if (!intervals.empty()) {
            EXPECT_TRUE((intervals[0].upper.is_pos_infinity)) << "Degenerate (a=0): upper should be +∞";
            EXPECT_TRUE((intervals[0].lower.is_open)) << "Degenerate (a=0): lower should be open (strict)";

            if (intervals[0].lower.value) {
                auto val_expr = intervals[0].lower.value->simplify();
                try {
                    double val = val_expr->to_numeric();
                    EXPECT_TRUE((std::abs(val - (2.0 / 3.0)) < 1e-10)) << "Degenerate (a=0): boundary should be 2/3";
                } catch (...) {
                    ADD_FAILURE() << "Degenerate (a=0): boundary should be evaluable";
                }
            }
        }
    }
}

TEST(ParametricInequalityDegenerate, DegenerateCaseAx250WhenA0ConstantRemainder) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");

    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto ax2 = SymbolicExpr::multiply(a, x2);
    auto expr = SymbolicExpr::add(ax2, SymbolicExpr::number(5));

    auto result = InequalitySolver::solve_parametric_inequality(
        expr, InequalityType::GreaterThan, "x", {"a"});

    EXPECT_TRUE((result.cases.size() == 3)) << "ax^2+5>0: should have 3 cases";

    if (result.cases.size() == 3) {

        auto &degen = result.cases[2];
        EXPECT_TRUE((degen.solution.is_entire_line())) << "Degenerate (a=0): 5>0 should be entire line";
    }
}

TEST(ParametricInequalityDegenerate, AxB0A0NonStrict) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");

    auto expr = SymbolicExpr::add(
        SymbolicExpr::multiply(a, x),
        b);

    auto result = InequalitySolver::solve_parametric_inequality(
        expr, InequalityType::LessEqual, "x", {"a", "b"});

    EXPECT_TRUE((!result.cases.empty())) << "ax+b<=0: should have cases";

    if (!result.cases.empty()) {
        auto &pos_case = result.cases[0];
        auto &intervals = pos_case.solution.intervals();

        EXPECT_TRUE((!intervals.empty())) << "ax+b<=0, a>0: should have intervals";

        if (!intervals.empty()) {
            EXPECT_TRUE((intervals[0].lower.is_neg_infinity)) << "ax+b<=0, a>0: lower should be -∞";
            EXPECT_TRUE((!intervals[0].upper.is_pos_infinity)) << "ax+b<=0, a>0: upper should not be +∞";

            EXPECT_TRUE((!intervals[0].upper.is_open)) << "ax+b<=0, a>0: upper should be closed (non-strict ≤)";

            if (intervals[0].upper.value) {
                auto boundary = intervals[0].upper.value;
                auto b1 = boundary->substitute("a", SymbolicExpr::number(2));
                b1 = b1->substitute("b", SymbolicExpr::number(-6));
                b1 = b1->simplify();
                try {
                    double val = b1->to_numeric();
                    EXPECT_TRUE((std::abs(val - 3.0) < 1e-10)) << "ax+b<=0, a>0: a=2,b=-6 → boundary should be 3";
                } catch (...) {
                    ADD_FAILURE() << "ax+b<=0, a>0: boundary should be evaluable";
                }
            }
        }
    }
}

TEST(ParametricInequalityDegenerate, DegenerateAx2X0WhenA0NonStrict) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");

    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto ax2 = SymbolicExpr::multiply(a, x2);
    auto expr = SymbolicExpr::add(ax2, x);

    auto result = InequalitySolver::solve_parametric_inequality(
        expr, InequalityType::GreaterEqual, "x", {"a"});

    EXPECT_TRUE((result.cases.size() == 3)) << "ax^2+x>=0: should have 3 cases";

    if (result.cases.size() == 3) {

        auto &degen = result.cases[2];
        EXPECT_TRUE((!degen.solution.is_empty())) << "Degenerate (a=0): x>=0 should not be empty";

        auto &intervals = degen.solution.intervals();
        EXPECT_TRUE((intervals.size() == 1)) << "Degenerate (a=0): x>=0 should have 1 interval";

        if (!intervals.empty()) {
            EXPECT_TRUE((intervals[0].upper.is_pos_infinity)) << "Degenerate (a=0): upper should be +∞";

            EXPECT_TRUE((!intervals[0].lower.is_open)) << "Degenerate (a=0): lower should be closed (non-strict ≥)";

            if (intervals[0].lower.value) {
                auto val_expr = intervals[0].lower.value->simplify();
                try {
                    double val = val_expr->to_numeric();
                    EXPECT_TRUE((std::abs(val - 0.0) < 1e-10)) << "Degenerate (a=0): boundary should be 0";
                } catch (...) {
                    ADD_FAILURE() << "Degenerate (a=0): boundary should be evaluable";
                }
            }
        }
    }
}
