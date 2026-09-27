#include <cmath>
#include <map>
#include "test_common.hpp"
#include "vector_calculus.hpp"
#include <memory>
#include <string>
#include <vector>

using namespace LMCAS;

TEST(VectorCalculusExtrema, FindExtremaQuadraticMin) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));

    auto extrema = find_extrema(f, {"x", "y"});

    EXPECT_TRUE((!extrema.empty())) << "found at least one critical point";
    if (!extrema.empty()) {
        EXPECT_TRUE((extrema[0].classification == CriticalPointClassification::LocalMinimum)) << "x^2+y^2 has minimum at origin";

        auto x_val = extrema[0].point.at("x");
        auto y_val = extrema[0].point.at("y");
        EXPECT_TRUE((x_val && x_val->is_zero())) << "critical point x = 0";
        EXPECT_TRUE((y_val && y_val->is_zero())) << "critical point y = 0";
    }
}

TEST(VectorCalculusExtrema, FindExtremaQuadraticMax) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto f = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
                               SymbolicExpr::power(x, SymbolicExpr::number(2))),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
                               SymbolicExpr::power(y, SymbolicExpr::number(2))));

    auto extrema = find_extrema(f, {"x", "y"});

    EXPECT_TRUE((!extrema.empty())) << "found at least one critical point";
    if (!extrema.empty()) {
        EXPECT_TRUE((extrema[0].classification == CriticalPointClassification::LocalMaximum)) << "-x^2-y^2 has maximum at origin";
    }
}

TEST(VectorCalculusExtrema, FindExtremaSaddle) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::multiply(SymbolicExpr::number(-1),
                               SymbolicExpr::power(y, SymbolicExpr::number(2))));

    auto extrema = find_extrema(f, {"x", "y"});

    EXPECT_TRUE((!extrema.empty())) << "found at least one critical point";
    if (!extrema.empty()) {
        EXPECT_TRUE((extrema[0].classification == CriticalPointClassification::Saddle)) << "x^2-y^2 has saddle at origin";
    }
}

TEST(VectorCalculusExtrema, FindExtremaDegenerate) {
    auto x = SymbolicExpr::variable("x");
    auto objective = SymbolicExpr::power(
        x, SymbolicExpr::number(2));

    auto extrema = find_extrema_checked(
        objective, {"x", "y"});
    ASSERT_TRUE(extrema.has_value()) << extrema.error().message;
    ASSERT_FALSE(extrema.value().empty());
    EXPECT_EQ(extrema.value().front().classification,
              CriticalPointClassification::Degenerate);
}

TEST(VectorCalculusExtrema, FindExtremaSingleVar) {
    auto x = SymbolicExpr::variable("x");

    auto f = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::multiply(SymbolicExpr::number(-2), x)),
        SymbolicExpr::number(1));

    auto extrema = find_extrema(f, {"x"});

    EXPECT_TRUE((!extrema.empty())) << "found critical point for (x-1)^2";
    if (!extrema.empty()) {
        EXPECT_TRUE((extrema[0].classification == CriticalPointClassification::LocalMinimum)) << "(x-1)^2 has minimum";

        auto x_val = extrema[0].point.at("x");
        if (x_val) {
            auto val = x_val->simplify();
            auto num = test_numeric_eval(val);
            EXPECT_TRUE((num.has_value() && std::abs(*num - 1.0) < 1e-10)) << "critical point at x=1";
        }
    }
}

TEST(VectorCalculusExtrema, ThreeDimensionalExtremeCoordinatesRemainInRange) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto x_center = SymbolicExpr::number(BigInt("3000000000"));
    auto y_center = SymbolicExpr::number(BigInt("-3000000000"));
    auto z_center = SymbolicExpr::number(BigInt("3000000001"));

    auto squared_offset = [](const std::shared_ptr<SymbolicExpr> &variable,
                             const std::shared_ptr<SymbolicExpr> &center) {
        return SymbolicExpr::power(
            SymbolicExpr::add(
                variable,
                SymbolicExpr::multiply(SymbolicExpr::number(-1), center)),
            SymbolicExpr::number(2));
    };
    auto bowl = SymbolicExpr::add(
        squared_offset(x, x_center),
        SymbolicExpr::add(
            squared_offset(y, y_center),
            squared_offset(z, z_center)));

    struct Case {
        std::shared_ptr<SymbolicExpr> expression;
        CriticalPointClassification classification;
    };
    const std::vector<Case> cases{
        {bowl, CriticalPointClassification::LocalMinimum},
        {SymbolicExpr::multiply(SymbolicExpr::number(-1), bowl),
         CriticalPointClassification::LocalMaximum}};
    const std::map<std::string, double> expected{
        {"x", 3000000000.0},
        {"y", -3000000000.0},
        {"z", 3000000001.0}};

    for (const auto &test_case : cases) {
        auto result = find_extrema_checked(
            test_case.expression, {"x", "y", "z"});
        ASSERT_TRUE(result.has_value()) << result.error().message;
        ASSERT_EQ(result.value().size(), 1u);
        const auto &critical = result.value().front();
        EXPECT_EQ(critical.classification, test_case.classification);
        for (const auto &[variable, coordinate] : expected) {
            auto found = critical.point.find(variable);
            ASSERT_NE(found, critical.point.end());
            auto value = test_numeric_eval(found->second);
            ASSERT_TRUE(value.has_value());
            ASSERT_TRUE(std::isfinite(*value));
            EXPECT_NEAR(*value, coordinate, 1e-6);
        }
    }
}

static void expect_circle_stationary_points(
    const std::vector<std::map<std::string, std::shared_ptr<SymbolicExpr>>> &solutions) {
    ASSERT_EQ(solutions.size(), 2u) << "circle has both antipodal stationary points";
    bool positive = false;
    bool negative = false;
    for (const auto &solution : solutions) {
        ASSERT_NE(solution.count("x"), 0u) << "stationary point assigns x";
        ASSERT_NE(solution.count("y"), 0u) << "stationary point assigns y";
        auto xn = test_numeric_eval(solution.at("x"));
        auto yn = test_numeric_eval(solution.at("y"));
        ASSERT_TRUE(xn.has_value()) << "stationary x evaluates as a real value";
        ASSERT_TRUE(yn.has_value()) << "stationary y evaluates as a real value";
        ASSERT_TRUE(std::isfinite(*xn));
        ASSERT_TRUE(std::isfinite(*yn));
        EXPECT_NEAR(*xn, *yn, 1e-10) << "stationarity requires x=y";
        const double radius_squared = *xn * *xn + *yn * *yn;
        EXPECT_TRUE(std::isfinite(radius_squared));
        EXPECT_NEAR(radius_squared, 1.0, 1e-10) << "stationary point lies on the circle";
        if (std::abs(*xn - std::sqrt(0.5)) <= 1e-10) {
            positive = true;
        }
        if (std::abs(*xn + std::sqrt(0.5)) <= 1e-10) {
            negative = true;
        }
    }
    EXPECT_TRUE(positive) << "maximum survives the algebraic back substitution";
    EXPECT_TRUE(negative) << "minimum survives the algebraic back substitution";
}

TEST(VectorCalculusExtrema, LagrangeBasic) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto f = SymbolicExpr::add(x, y);
    auto g = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::power(y, SymbolicExpr::number(2))),
        SymbolicExpr::number(-1));

    auto checked = lagrange_multipliers_checked(f, {g}, {"x", "y"});
    ASSERT_TRUE((checked.has_value())) << (checked ? "circle stationarity solve succeeds" : checked.error().operation + ": " + checked.error().message);
    if (!checked) {
        return;
    }
    const auto &solutions = checked.value();
    expect_circle_stationary_points(solutions);
}

TEST(VectorCalculusExtrema, LagrangeLinearConstraint) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto f = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));
    auto g = SymbolicExpr::add(
        SymbolicExpr::add(x, y),
        SymbolicExpr::number(-1));

    auto solutions = lagrange_multipliers(f, {g}, {"x", "y"});

    EXPECT_TRUE((!solutions.empty())) << "found critical point";

    if (!solutions.empty()) {
        auto x_val = solutions[0].at("x")->simplify();
        auto y_val = solutions[0].at("y")->simplify();
        auto xn = test_numeric_eval(x_val);
        auto yn = test_numeric_eval(y_val);
        if (xn.has_value() && yn.has_value()) {
            EXPECT_NEAR(*xn, 0.5, 1e-8) << "x = 1/2";
            EXPECT_NEAR(*yn, 0.5, 1e-8) << "y = 1/2";
        }
    }
}

TEST(VectorCalculusExtrema, ExtremaCheckedMinimum) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto parabola = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::multiply(SymbolicExpr::number(-2), x)),
        SymbolicExpr::number(1));

    auto extrema = find_extrema_checked(parabola, {"x"});
    ASSERT_TRUE((extrema.has_value())) << "checked extrema succeeds for single-variable quadratic";
    if (extrema) {
        EXPECT_TRUE((!extrema.value().empty())) << "checked extrema returns a critical point";
        EXPECT_TRUE((extrema.value()[0].classification ==
                     CriticalPointClassification::LocalMinimum))
            << "checked extrema classifies quadratic minimum";
    }
}

TEST(VectorCalculusExtrema, LagrangeCheckedCandidate) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto objective = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));
    auto constraint = SymbolicExpr::add(SymbolicExpr::add(x, y), SymbolicExpr::number(-1));

    auto lagrange = lagrange_multipliers_checked(objective, {constraint}, {"x", "y"});
    ASSERT_TRUE((lagrange.has_value())) << "checked Lagrange succeeds for linear constraint";
    if (lagrange) {
        EXPECT_TRUE((!lagrange.value().empty())) << "checked Lagrange returns candidates";
        auto x_val = lagrange.value()[0].at("x")->simplify();
        auto y_val = lagrange.value()[0].at("y")->simplify();
        auto xn = test_numeric_eval(x_val);
        auto yn = test_numeric_eval(y_val);
        EXPECT_TRUE((xn.has_value() && yn.has_value())) << "checked Lagrange candidate is numeric";
        if (xn.has_value() && yn.has_value()) {
            EXPECT_NEAR(*xn, 0.5, 1e-8) << "checked Lagrange x = 1/2";
            EXPECT_NEAR(*yn, 0.5, 1e-8) << "checked Lagrange y = 1/2";
        }
    }
}

TEST(VectorCalculusExtrema, ExtremaCheckedNull) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto null_extrema = find_extrema_checked(nullptr, {"x"});
    EXPECT_TRUE((!null_extrema.has_value())) << "checked extrema rejects null expression";
    EXPECT_TRUE((null_extrema.error().code == CasErrc::InvalidArgument)) << "checked extrema reports InvalidArgument for null expression";

    std::shared_ptr<SymbolicExpr> null_root;
    auto null_root_extrema = find_extrema_checked(null_root, {"x"});
    EXPECT_TRUE((!null_root_extrema.has_value())) << "checked extrema rejects null expression";
    EXPECT_TRUE((null_root_extrema.error().code == CasErrc::InvalidArgument)) << "checked extrema reports InvalidArgument for null expression";
}

TEST(VectorCalculusExtrema, ExtremaCheckedDuplicateVariables) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto objective = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));
    auto constraint = SymbolicExpr::add(SymbolicExpr::add(x, y), SymbolicExpr::number(-1));

    auto duplicate_vars = find_extrema_checked(objective, {"x", "x"});
    EXPECT_TRUE((!duplicate_vars.has_value())) << "checked extrema rejects duplicate variables";
    EXPECT_TRUE((duplicate_vars.error().code == CasErrc::InvalidArgument)) << "checked extrema reports InvalidArgument for duplicate variables";
}

TEST(VectorCalculusExtrema, ExtremaCheckedUnsupported) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto unsupported_objective = SymbolicExpr::eq(x, SymbolicExpr::number(0));
    auto unsupported_extrema = find_extrema_checked(unsupported_objective, {"x"});
    EXPECT_TRUE((!unsupported_extrema.has_value())) << "checked extrema rejects unsupported objective derivatives";
    EXPECT_TRUE((unsupported_extrema.error().code == CasErrc::Inconclusive)) << "checked extrema reports Inconclusive for unsupported derivatives";
}

TEST(VectorCalculusExtrema, LagrangeCheckedEmptyConstraints) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto objective = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));
    auto constraint = SymbolicExpr::add(SymbolicExpr::add(x, y), SymbolicExpr::number(-1));

    auto empty_constraints = lagrange_multipliers_checked(objective, {}, {"x", "y"});
    EXPECT_TRUE((!empty_constraints.has_value())) << "checked Lagrange rejects empty constraints";
    EXPECT_TRUE((empty_constraints.error().code == CasErrc::InvalidArgument)) << "checked Lagrange reports InvalidArgument for empty constraints";
}

TEST(VectorCalculusExtrema, LagrangeCheckedNullConstraint) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto objective = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));
    auto constraint = SymbolicExpr::add(SymbolicExpr::add(x, y), SymbolicExpr::number(-1));

    std::shared_ptr<SymbolicExpr> null_root;

    auto null_constraint = lagrange_multipliers_checked(objective, {null_root}, {"x", "y"});
    EXPECT_TRUE((!null_constraint.has_value())) << "checked Lagrange rejects null constraints";
    EXPECT_TRUE((null_constraint.error().code == CasErrc::InvalidArgument)) << "checked Lagrange reports InvalidArgument for null constraints";
}

TEST(VectorCalculusExtrema, LagrangeCheckedUnsupported) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto objective = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));
    auto constraint = SymbolicExpr::add(SymbolicExpr::add(x, y), SymbolicExpr::number(-1));

    auto unsupported_objective = SymbolicExpr::eq(x, SymbolicExpr::number(0));

    auto unsupported_lagrange = lagrange_multipliers_checked(
        unsupported_objective, {constraint}, {"x", "y"});
    EXPECT_TRUE((!unsupported_lagrange.has_value())) << "checked Lagrange rejects unsupported stationarity derivatives";
    EXPECT_TRUE((unsupported_lagrange.error().code == CasErrc::Inconclusive)) << "checked Lagrange reports Inconclusive for unsupported derivatives";
}

TEST(VectorCalculusExtrema, ExtremaCheckedCancellation) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto parabola = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::power(x, SymbolicExpr::number(2)),
            SymbolicExpr::multiply(SymbolicExpr::number(-2), x)),
        SymbolicExpr::number(1));

    LMCAS::CancellationToken cancellation;
    LMCAS::ComputationContext cancelled_context({}, cancellation);
    cancellation.cancel();
    auto cancelled = find_extrema_checked(parabola, {"x"}, cancelled_context);
    EXPECT_TRUE((!cancelled.has_value())) << "checked extrema observes cancellation";
    EXPECT_TRUE((cancelled.error().code == CasErrc::Cancelled)) << "checked extrema reports Cancelled";
}

TEST(VectorCalculusExtrema, LagrangeCheckedBudget) {
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");

    auto objective = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::power(y, SymbolicExpr::number(2)));
    auto constraint = SymbolicExpr::add(SymbolicExpr::add(x, y), SymbolicExpr::number(-1));

    LMCAS::ResourceLimits limits;
    limits.max_steps = 1;
    LMCAS::ComputationContext limited_context(limits);
    auto limited = lagrange_multipliers_checked(objective, {constraint}, {"x", "y"},
                                                limited_context);
    EXPECT_TRUE((!limited.has_value())) << "checked Lagrange observes exhausted step budget";
    EXPECT_TRUE((limited.error().code == CasErrc::ResourceLimit)) << "checked Lagrange reports ResourceLimit";
}
