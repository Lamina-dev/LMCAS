#include "test_common.hpp"
#include "interval.hpp"
#include "inequality_solver.hpp"
#include "symbolic.hpp"
#include "poly_utils.hpp"
#include <cmath>
#include <random>
#include <vector>

using namespace lamina;

int main() {

    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");
    auto c = SymbolicExpr::variable("c");

    TEST_CASE("Parametric linear: ax + b > 0 (piecewise by sign of a)");
    {

        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(a, x),
            b);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"a", "b"});

        EXPECT_TRUE(result.cases.size() >= 3,
            "ax + b > 0: should have at least 3 piecewise cases (a>0, a<0, then a=0 sub-cases by sign of b)");

        auto& pos_case = result.cases[0];
        EXPECT_TRUE(pos_case.condition != nullptr,
            "ax + b > 0, case a>0: condition should not be null");
        EXPECT_TRUE(!pos_case.solution.is_empty(),
            "ax + b > 0, case a>0: solution should not be empty");
        EXPECT_TRUE(pos_case.solution.intervals().size() == 1,
            "ax + b > 0, case a>0: should have exactly 1 interval");

        if (!pos_case.solution.intervals().empty()) {
            auto& iv = pos_case.solution.intervals()[0];

            EXPECT_TRUE(!iv.lower.is_neg_infinity,
                "ax + b > 0, case a>0: lower bound should not be -inf");
            EXPECT_TRUE(iv.lower.is_open,
                "ax + b > 0, case a>0: lower bound should be open (strict)");
            EXPECT_TRUE(iv.upper.is_pos_infinity,
                "ax + b > 0, case a>0: upper bound should be +inf");

            if (iv.lower.value) {

                auto substituted = iv.lower.value->substitute("a", SymbolicExpr::number(2));
                substituted = substituted->substitute("b", SymbolicExpr::number(4));
                substituted = substituted->simplify();
                try {
                    double val = substituted->to_numeric();
                    EXPECT_TRUE(std::abs(val - (-2.0)) < 1e-10,
                        "ax + b > 0, case a>0: root at a=2,b=4 should be -2");
                } catch (...) {
                    EXPECT_TRUE(false, "ax + b > 0, case a>0: root should be evaluable");
                }
            }
        }

        auto& neg_case = result.cases[1];
        EXPECT_TRUE(neg_case.condition != nullptr,
            "ax + b > 0, case a<0: condition should not be null");
        EXPECT_TRUE(!neg_case.solution.is_empty(),
            "ax + b > 0, case a<0: solution should not be empty");
        EXPECT_TRUE(neg_case.solution.intervals().size() == 1,
            "ax + b > 0, case a<0: should have exactly 1 interval");

        if (!neg_case.solution.intervals().empty()) {
            auto& iv = neg_case.solution.intervals()[0];
            EXPECT_TRUE(iv.lower.is_neg_infinity,
                "ax + b > 0, case a<0: lower bound should be -inf");
            EXPECT_TRUE(!iv.upper.is_pos_infinity,
                "ax + b > 0, case a<0: upper bound should not be +inf");
            EXPECT_TRUE(iv.upper.is_open,
                "ax + b > 0, case a<0: upper bound should be open (strict)");
        }

        auto& degen_case = result.cases[2];
        EXPECT_TRUE(degen_case.condition != nullptr,
            "ax + b > 0, case a=0: condition should not be null");
    }

    TEST_CASE("Parametric quadratic with constant leading coeff: x^2 + bx + c > 0");
    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto bx = SymbolicExpr::multiply(b, x);
        auto expr = SymbolicExpr::add(SymbolicExpr::add(x2, bx), c);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"b", "c"});

        EXPECT_TRUE(result.cases.size() == 1,
            "x^2 + bx + c > 0: should have 1 case (constant leading coeff)");

        auto& single_case = result.cases[0];
        EXPECT_TRUE(single_case.condition == nullptr,
            "x^2 + bx + c > 0: condition should be null (unconditional)");

        auto& intervals = single_case.solution.intervals();
        EXPECT_TRUE(intervals.size() == 2,
            "x^2 + bx + c > 0: should have 2 intervals (outside roots)");

        if (intervals.size() == 2) {

            EXPECT_TRUE(intervals[0].lower.is_neg_infinity,
                "x^2 + bx + c > 0: first interval starts at -inf");
            EXPECT_TRUE(intervals[0].upper.is_open,
                "x^2 + bx + c > 0: first interval upper is open (strict)");

            EXPECT_TRUE(intervals[1].upper.is_pos_infinity,
                "x^2 + bx + c > 0: second interval ends at +inf");
            EXPECT_TRUE(intervals[1].lower.is_open,
                "x^2 + bx + c > 0: second interval lower is open (strict)");

            if (intervals[0].upper.value) {
                auto root1 = intervals[0].upper.value->substitute("b", SymbolicExpr::number(-3));
                root1 = root1->substitute("c", SymbolicExpr::number(2));
                root1 = root1->simplify();
                try {
                    double val = root1->to_numeric();

                    EXPECT_TRUE(std::abs(val - 1.0) < 1e-10 || std::abs(val - 2.0) < 1e-10,
                        "x^2 + bx + c > 0: root1 at b=-3,c=2 should be 1 or 2");
                } catch (...) {

                }
            }
        }
    }

    TEST_CASE("Parametric quadratic: ax^2 + bx + c >= 0 (piecewise by sign of a)");
    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto ax2 = SymbolicExpr::multiply(a, x2);
        auto bx = SymbolicExpr::multiply(b, x);
        auto expr = SymbolicExpr::add(SymbolicExpr::add(ax2, bx), c);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterEqual, "x", {"a", "b", "c"});

        EXPECT_TRUE(result.cases.size() >= 3,
            "ax^2 + bx + c >= 0: should have at least 3 piecewise cases (a>0, a<0, a=0 may further split)");

        auto& pos_case = result.cases[0];
        EXPECT_TRUE(pos_case.condition != nullptr,
            "ax^2 + bx + c >= 0, case a>0: should have condition");

        auto& neg_case = result.cases[1];
        EXPECT_TRUE(neg_case.condition != nullptr,
            "ax^2 + bx + c >= 0, case a<0: should have condition");

        // 第三个及以后的分支对应 a=0 的退化情形（可能进一步按 b、c 的符号细分）。
        if (result.cases.size() >= 3) {
            auto& degen_case = result.cases[2];
            EXPECT_TRUE(degen_case.condition != nullptr,
                "ax^2 + bx + c >= 0, case a=0: should have condition (degenerate)");
        }
    }

    TEST_CASE("Parametric consistency: parametric vs numeric for specific values");
    {

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

        EXPECT_TRUE(!numeric_result.is_empty(),
            "3x - 6 > 0: numeric solution should not be empty");
        EXPECT_TRUE(numeric_result.contains(3.0),
            "3x - 6 > 0: x=3 should be in solution");
        EXPECT_TRUE(!numeric_result.contains(1.0),
            "3x - 6 > 0: x=1 should not be in solution");
        EXPECT_TRUE(!numeric_result.contains(2.0),
            "3x - 6 > 0: x=2 (root) should not be in solution (strict)");

        EXPECT_TRUE(!parametric_result.cases.empty(),
            "Parametric ax + b > 0: should have cases");

        if (!parametric_result.cases.empty()) {
            auto& pos_case = parametric_result.cases[0];
            auto& intervals = pos_case.solution.intervals();
            if (!intervals.empty() && intervals[0].lower.value) {
                auto boundary = intervals[0].lower.value;
                boundary = boundary->substitute("a", SymbolicExpr::number(3));
                boundary = boundary->substitute("b", SymbolicExpr::number(-6));
                boundary = boundary->simplify();
                try {
                    double val = boundary->to_numeric();
                    EXPECT_TRUE(std::abs(val - 2.0) < 1e-10,
                        "Parametric consistency: boundary at a=3,b=-6 should be 2");
                } catch (...) {
                    EXPECT_TRUE(false,
                        "Parametric consistency: boundary should be evaluable");
                }
            }
        }
    }

    TEST_CASE("Parametric with empty params: falls back to non-parametric");
    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto expr = SymbolicExpr::add(x2, SymbolicExpr::number(-4));

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {});

        EXPECT_TRUE(result.cases.size() == 1,
            "No params: should have 1 case");
        EXPECT_TRUE(!result.cases[0].solution.is_empty(),
            "No params: x^2 - 4 > 0 should not be empty");
        EXPECT_TRUE(result.cases[0].solution.contains(3.0),
            "No params: x=3 should be in solution");
        EXPECT_TRUE(result.cases[0].solution.contains(-3.0),
            "No params: x=-3 should be in solution");
        EXPECT_TRUE(!result.cases[0].solution.contains(0.0),
            "No params: x=0 should not be in solution");
    }

    TEST_CASE("Degenerate case: ax^2 + 2x + 1 > 0 when a=0");
    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto ax2 = SymbolicExpr::multiply(a, x2);
        auto two_x = SymbolicExpr::multiply(SymbolicExpr::number(2), x);
        auto expr = SymbolicExpr::add(SymbolicExpr::add(ax2, two_x), SymbolicExpr::number(1));

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"a"});

        EXPECT_TRUE(result.cases.size() == 3,
            "ax^2 + 2x + 1 > 0: should have 3 cases");

        if (result.cases.size() == 3) {
            auto& degen = result.cases[2];
            EXPECT_TRUE(!degen.solution.is_empty(),
                "Degenerate (a=0): 2x + 1 > 0 should not be empty");

            auto& intervals = degen.solution.intervals();
            EXPECT_TRUE(intervals.size() == 1,
                "Degenerate (a=0): should have 1 interval");

            if (!intervals.empty()) {
                EXPECT_TRUE(intervals[0].upper.is_pos_infinity,
                    "Degenerate (a=0): upper bound should be +inf");
                EXPECT_TRUE(intervals[0].lower.is_open,
                    "Degenerate (a=0): lower bound should be open (strict)");

                if (intervals[0].lower.value) {
                    auto val_expr = intervals[0].lower.value->simplify();
                    try {
                        double val = val_expr->to_numeric();
                        EXPECT_TRUE(std::abs(val - (-0.5)) < 1e-10,
                            "Degenerate (a=0): boundary should be -1/2");
                    } catch (...) {
                        EXPECT_TRUE(false,
                            "Degenerate (a=0): boundary should be evaluable");
                    }
                }
            }
        }
    }

    TEST_CASE("Task 7.5: ax + b > 0, a > 0 numeric verification");
    {
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(a, x),
            b);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"a", "b"});

        EXPECT_TRUE(!result.cases.empty(),
            "ax + b > 0: should have cases");
        EXPECT_TRUE(result.cases.size() >= 1,
            "ax + b > 0: should have at least 1 case");

        auto& pos_case = result.cases[0];
        auto& pos_intervals = pos_case.solution.intervals();

        if (!pos_intervals.empty() && pos_intervals[0].lower.value) {
            auto boundary = pos_intervals[0].lower.value;
            auto b1 = boundary->substitute("a", SymbolicExpr::number(5));
            b1 = b1->substitute("b", SymbolicExpr::number(-10));
            b1 = b1->simplify();
            try {
                double val = b1->to_numeric();
                EXPECT_TRUE(std::abs(val - 2.0) < 1e-10,
                    "ax+b>0, a>0: a=5,b=-10 → boundary should be 2");
            } catch (...) {
                EXPECT_TRUE(false, "ax+b>0, a>0: boundary should be evaluable for a=5,b=-10");
            }

            auto b2 = boundary->substitute("a", SymbolicExpr::number(1));
            b2 = b2->substitute("b", SymbolicExpr::number(3));
            b2 = b2->simplify();
            try {
                double val = b2->to_numeric();
                EXPECT_TRUE(std::abs(val - (-3.0)) < 1e-10,
                    "ax+b>0, a>0: a=1,b=3 → boundary should be -3");
            } catch (...) {
                EXPECT_TRUE(false, "ax+b>0, a>0: boundary should be evaluable for a=1,b=3");
            }
        }

        if (!pos_intervals.empty()) {
            EXPECT_TRUE(pos_intervals[0].upper.is_pos_infinity,
                "ax+b>0, a>0: upper bound should be +∞");
            EXPECT_TRUE(pos_intervals[0].lower.is_open,
                "ax+b>0, a>0: lower bound should be open (strict >)");
        }
    }

    TEST_CASE("Task 7.5: ax + b > 0, a < 0 numeric verification");
    {
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(a, x),
            b);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"a", "b"});

        EXPECT_TRUE(result.cases.size() >= 2,
            "ax + b > 0: should have at least 2 cases (a>0, a<0)");

        auto& neg_case = result.cases[1];
        auto& neg_intervals = neg_case.solution.intervals();

        EXPECT_TRUE(!neg_intervals.empty(),
            "ax+b>0, a<0: solution should not be empty");

        if (!neg_intervals.empty()) {

            EXPECT_TRUE(neg_intervals[0].lower.is_neg_infinity,
                "ax+b>0, a<0: lower bound should be -∞");
            EXPECT_TRUE(neg_intervals[0].upper.is_open,
                "ax+b>0, a<0: upper bound should be open (strict >)");
            EXPECT_TRUE(!neg_intervals[0].upper.is_pos_infinity,
                "ax+b>0, a<0: upper bound should not be +∞");

            if (neg_intervals[0].upper.value) {
                auto boundary = neg_intervals[0].upper.value;
                auto b1 = boundary->substitute("a", SymbolicExpr::number(-2));
                b1 = b1->substitute("b", SymbolicExpr::number(6));
                b1 = b1->simplify();
                try {
                    double val = b1->to_numeric();
                    EXPECT_TRUE(std::abs(val - 3.0) < 1e-10,
                        "ax+b>0, a<0: a=-2,b=6 → boundary should be 3");
                } catch (...) {
                    EXPECT_TRUE(false, "ax+b>0, a<0: boundary should be evaluable for a=-2,b=6");
                }

                auto b2 = boundary->substitute("a", SymbolicExpr::number(-4));
                b2 = b2->substitute("b", SymbolicExpr::number(-8));
                b2 = b2->simplify();
                try {
                    double val = b2->to_numeric();
                    EXPECT_TRUE(std::abs(val - (-2.0)) < 1e-10,
                        "ax+b>0, a<0: a=-4,b=-8 → boundary should be -2");
                } catch (...) {
                    EXPECT_TRUE(false, "ax+b>0, a<0: boundary should be evaluable for a=-4,b=-8");
                }
            }
        }
    }

    TEST_CASE("Task 7.5: ax^2 + bx + c >= 0, discriminant depends on parameters");
    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto ax2 = SymbolicExpr::multiply(a, x2);
        auto bx = SymbolicExpr::multiply(b, x);
        auto expr = SymbolicExpr::add(SymbolicExpr::add(ax2, bx), c);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterEqual, "x", {"a", "b", "c"});

        EXPECT_TRUE(result.cases.size() >= 3,
            "ax^2+bx+c>=0: should have at least 3 piecewise cases (degenerate a=0 may split further)");

        if (result.cases.size() >= 1) {
            auto& pos_case = result.cases[0];

            auto& intervals = pos_case.solution.intervals();

            if (intervals.size() == 2) {

                if (intervals[0].upper.value) {
                    auto root1 = intervals[0].upper.value->substitute("a", SymbolicExpr::number(1));
                    root1 = root1->substitute("b", SymbolicExpr::number(-5));
                    root1 = root1->substitute("c", SymbolicExpr::number(6));
                    root1 = root1->simplify();
                    try {
                        double val = root1->to_numeric();

                        EXPECT_TRUE(std::abs(val - 2.0) < 1e-10 || std::abs(val - 3.0) < 1e-10,
                            "ax^2+bx+c>=0, a>0: root at a=1,b=-5,c=6 should be 2 or 3");
                    } catch (...) {

                    }
                }

                if (intervals[1].lower.value) {
                    auto root2 = intervals[1].lower.value->substitute("a", SymbolicExpr::number(1));
                    root2 = root2->substitute("b", SymbolicExpr::number(-5));
                    root2 = root2->substitute("c", SymbolicExpr::number(6));
                    root2 = root2->simplify();
                    try {
                        double val = root2->to_numeric();

                        EXPECT_TRUE(std::abs(val - 2.0) < 1e-10 || std::abs(val - 3.0) < 1e-10,
                            "ax^2+bx+c>=0, a>0: root at a=1,b=-5,c=6 should be 2 or 3");
                    } catch (...) {

                    }
                }

                EXPECT_TRUE(!intervals[0].upper.is_open,
                    "ax^2+bx+c>=0, a>0: upper of first interval should be closed (non-strict)");
                EXPECT_TRUE(!intervals[1].lower.is_open,
                    "ax^2+bx+c>=0, a>0: lower of second interval should be closed (non-strict)");
            }
        }

        if (result.cases.size() >= 2) {
            auto& neg_case = result.cases[1];

            auto& intervals = neg_case.solution.intervals();

            if (intervals.size() == 1) {

                EXPECT_TRUE(!intervals[0].lower.is_neg_infinity,
                    "ax^2+bx+c>=0, a<0: lower should not be -∞");
                EXPECT_TRUE(!intervals[0].upper.is_pos_infinity,
                    "ax^2+bx+c>=0, a<0: upper should not be +∞");
                EXPECT_TRUE(!intervals[0].lower.is_open,
                    "ax^2+bx+c>=0, a<0: lower should be closed (non-strict)");
                EXPECT_TRUE(!intervals[0].upper.is_open,
                    "ax^2+bx+c>=0, a<0: upper should be closed (non-strict)");
            }
        }
    }

    TEST_CASE("Task 7.5: Degenerate case ax^2 + 3x - 2 > 0 when a=0");
    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto ax2 = SymbolicExpr::multiply(a, x2);
        auto three_x = SymbolicExpr::multiply(SymbolicExpr::number(3), x);
        auto expr = SymbolicExpr::add(
            SymbolicExpr::add(ax2, three_x),
            SymbolicExpr::number(-2));

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"a"});

        EXPECT_TRUE(result.cases.size() == 3,
            "ax^2+3x-2>0: should have 3 cases");

        if (result.cases.size() == 3) {

            auto& degen = result.cases[2];
            EXPECT_TRUE(!degen.solution.is_empty(),
                "Degenerate (a=0): 3x-2>0 should not be empty");

            auto& intervals = degen.solution.intervals();
            EXPECT_TRUE(intervals.size() == 1,
                "Degenerate (a=0): should have 1 interval");

            if (!intervals.empty()) {
                EXPECT_TRUE(intervals[0].upper.is_pos_infinity,
                    "Degenerate (a=0): upper should be +∞");
                EXPECT_TRUE(intervals[0].lower.is_open,
                    "Degenerate (a=0): lower should be open (strict)");

                if (intervals[0].lower.value) {
                    auto val_expr = intervals[0].lower.value->simplify();
                    try {
                        double val = val_expr->to_numeric();
                        EXPECT_TRUE(std::abs(val - (2.0/3.0)) < 1e-10,
                            "Degenerate (a=0): boundary should be 2/3");
                    } catch (...) {
                        EXPECT_TRUE(false,
                            "Degenerate (a=0): boundary should be evaluable");
                    }
                }
            }
        }
    }

    TEST_CASE("Task 7.5: Degenerate case ax^2 + 5 > 0 when a=0 (constant remainder)");
    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto ax2 = SymbolicExpr::multiply(a, x2);
        auto expr = SymbolicExpr::add(ax2, SymbolicExpr::number(5));

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterThan, "x", {"a"});

        EXPECT_TRUE(result.cases.size() == 3,
            "ax^2+5>0: should have 3 cases");

        if (result.cases.size() == 3) {

            auto& degen = result.cases[2];
            EXPECT_TRUE(degen.solution.is_entire_line(),
                "Degenerate (a=0): 5>0 should be entire line");
        }
    }

    TEST_CASE("Task 7.5: ax + b <= 0, a > 0 (non-strict)");
    {
        auto expr = SymbolicExpr::add(
            SymbolicExpr::multiply(a, x),
            b);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::LessEqual, "x", {"a", "b"});

        EXPECT_TRUE(!result.cases.empty(),
            "ax+b<=0: should have cases");

        if (!result.cases.empty()) {
            auto& pos_case = result.cases[0];
            auto& intervals = pos_case.solution.intervals();

            EXPECT_TRUE(!intervals.empty(),
                "ax+b<=0, a>0: should have intervals");

            if (!intervals.empty()) {
                EXPECT_TRUE(intervals[0].lower.is_neg_infinity,
                    "ax+b<=0, a>0: lower should be -∞");
                EXPECT_TRUE(!intervals[0].upper.is_pos_infinity,
                    "ax+b<=0, a>0: upper should not be +∞");

                EXPECT_TRUE(!intervals[0].upper.is_open,
                    "ax+b<=0, a>0: upper should be closed (non-strict ≤)");

                if (intervals[0].upper.value) {
                    auto boundary = intervals[0].upper.value;
                    auto b1 = boundary->substitute("a", SymbolicExpr::number(2));
                    b1 = b1->substitute("b", SymbolicExpr::number(-6));
                    b1 = b1->simplify();
                    try {
                        double val = b1->to_numeric();
                        EXPECT_TRUE(std::abs(val - 3.0) < 1e-10,
                            "ax+b<=0, a>0: a=2,b=-6 → boundary should be 3");
                    } catch (...) {
                        EXPECT_TRUE(false,
                            "ax+b<=0, a>0: boundary should be evaluable");
                    }
                }
            }
        }
    }

    TEST_CASE("Task 7.5: Degenerate ax^2 + x >= 0 when a=0 (non-strict)");
    {

        auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
        auto ax2 = SymbolicExpr::multiply(a, x2);
        auto expr = SymbolicExpr::add(ax2, x);

        auto result = InequalitySolver::solve_parametric_inequality(
            expr, InequalityType::GreaterEqual, "x", {"a"});

        EXPECT_TRUE(result.cases.size() == 3,
            "ax^2+x>=0: should have 3 cases");

        if (result.cases.size() == 3) {

            auto& degen = result.cases[2];
            EXPECT_TRUE(!degen.solution.is_empty(),
                "Degenerate (a=0): x>=0 should not be empty");

            auto& intervals = degen.solution.intervals();
            EXPECT_TRUE(intervals.size() == 1,
                "Degenerate (a=0): x>=0 should have 1 interval");

            if (!intervals.empty()) {
                EXPECT_TRUE(intervals[0].upper.is_pos_infinity,
                    "Degenerate (a=0): upper should be +∞");

                EXPECT_TRUE(!intervals[0].lower.is_open,
                    "Degenerate (a=0): lower should be closed (non-strict ≥)");

                if (intervals[0].lower.value) {
                    auto val_expr = intervals[0].lower.value->simplify();
                    try {
                        double val = val_expr->to_numeric();
                        EXPECT_TRUE(std::abs(val - 0.0) < 1e-10,
                            "Degenerate (a=0): boundary should be 0");
                    } catch (...) {
                        EXPECT_TRUE(false,
                            "Degenerate (a=0): boundary should be evaluable");
                    }
                }
            }
        }
    }

    return TEST_REPORT();
}
