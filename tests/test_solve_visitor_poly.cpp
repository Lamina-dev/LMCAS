#include <vector>
#include <string>
#include <memory>
#include "assumption_context.hpp"
#include "symbolic.hpp"
#include "test_common.hpp"

using namespace LMCAS;

TEST(LmcasSolveVisitorPoly, SolveNumeric) {
    auto x = SymbolicExpr::variable("x");
    auto eq1 = SymbolicExpr::add(
        SymbolicExpr::multiply(SymbolicExpr::number(2), x),
        SymbolicExpr::number(-4));

    auto solutions = LMCAS::solve_finite_checked(eq1, "x");
    ASSERT_TRUE(solutions) << solutions.error().message;
    ASSERT_EQ(solutions.value().size(), 1U);
    EXPECT_TRUE(test_proved_equivalent(
        solutions.value()[0], SymbolicExpr::number(2)));

    auto eq2 = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(-3), x),
            SymbolicExpr::number(2)));

    auto quadratic = LMCAS::solve_finite_checked(eq2, "x");
    ASSERT_TRUE(quadratic) << quadratic.error().message;
    ASSERT_EQ(quadratic.value().size(), 2U);
    bool found_one = false;
    bool found_two = false;
    for (const auto &root : quadratic.value()) {
        found_one = found_one || test_proved_equivalent(
                                     root, SymbolicExpr::number(1));
        found_two = found_two || test_proved_equivalent(
                                     root, SymbolicExpr::number(2));
    }
    EXPECT_TRUE(found_one);
    EXPECT_TRUE(found_two);
}

TEST(LmcasSolveVisitorPoly, SolveSymbolic) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");

    auto eq = SymbolicExpr::add(
        SymbolicExpr::multiply(a, x),
        b);

    auto solved = LMCAS::solve_equation(eq, "x");
    ASSERT_TRUE(solved) << solved.error().message;
    EXPECT_TRUE(std::holds_alternative<ConditionalSolutions>(
        solved.value()));
    auto projected = LMCAS::solve_finite_checked(eq, "x");
    ASSERT_FALSE(projected);
    EXPECT_EQ(projected.error().code, CasErrc::Inconclusive);
}

TEST(LmcasSolveVisitorPoly, QuadraticNonzeroDiscriminantHasTwoSimpleRoots) {
    auto x = SymbolicExpr::variable("x");
    auto equation = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-1));

    auto solved = solve_equation(equation, "x");

    ASSERT_TRUE(solved) << solved.error().message;
    const auto *finite = std::get_if<FiniteSolutions>(&solved.value());
    ASSERT_NE(finite, nullptr);
    ASSERT_EQ(finite->values.size(), 2U);
    bool found_negative = false;
    bool found_positive = false;
    for (const auto &root : finite->values) {
        EXPECT_EQ(root.multiplicity, 1U);
        EXPECT_TRUE(root.conditions.empty());
        found_negative =
            found_negative || test_proved_equivalent(
                                  root.value, SymbolicExpr::number(-1));
        found_positive =
            found_positive || test_proved_equivalent(
                                  root.value, SymbolicExpr::number(1));
    }
    EXPECT_TRUE(found_negative);
    EXPECT_TRUE(found_positive);
}

TEST(LmcasSolveVisitorPoly, ProvedDoubleRootHasMultiplicityTwo) {
    auto x = SymbolicExpr::variable("x");
    auto equation = SymbolicExpr::power(x, SymbolicExpr::number(2));

    auto solved = solve_equation(equation, "x");

    ASSERT_TRUE(solved) << solved.error().message;
    const auto *finite = std::get_if<FiniteSolutions>(&solved.value());
    ASSERT_NE(finite, nullptr);
    ASSERT_EQ(finite->values.size(), 1U);
    EXPECT_TRUE(test_proved_equivalent(
        finite->values.front().value, SymbolicExpr::number(0)));
    EXPECT_EQ(finite->values.front().multiplicity, 2U);
    EXPECT_TRUE(finite->values.front().conditions.empty());
}

TEST(LmcasSolveVisitorPoly, NonzeroSymbolicLeadingCoefficientKeepsDoubleRoot) {
    auto x = SymbolicExpr::variable("x");
    auto a = SymbolicExpr::variable("a");
    auto equation = SymbolicExpr::multiply(
        a, SymbolicExpr::power(x, SymbolicExpr::number(2)));
    auto assumptions = std::make_shared<AssumptionContext>();
    ASSERT_TRUE(assumptions->assume_sign_checked(
        "a", Sign::NonZero));
    ComputationContext context;
    ASSERT_TRUE(context.set_assumptions(assumptions));

    auto solved = solve_equation(
        equation, "x", context, SolveOptions{});

    ASSERT_TRUE(solved) << solved.error().message;
    const auto *finite = std::get_if<FiniteSolutions>(&solved.value());
    ASSERT_NE(finite, nullptr);
    ASSERT_EQ(finite->values.size(), 1U);
    EXPECT_TRUE(test_proved_equivalent(
        finite->values.front().value, SymbolicExpr::number(0)));
    EXPECT_EQ(finite->values.front().multiplicity, 2U);
}

static ExprPtr specialize_discriminant_condition(
    const ExprPtr &condition, int b, int c) {
    return condition->substitute("b", SymbolicExpr::number(b))
        ->substitute("c", SymbolicExpr::number(c));
}

static void expect_repeated_discriminant_branch(
    const FiniteSolution &root, const RelationalNode &relation,
    AssumptionContext &condition_evaluator) {
    EXPECT_EQ(relation.op(), RelationOp::EQ);
    auto specialized = specialize_discriminant_condition(root.value, 2, 1)
                           ->simplify();
    EXPECT_TRUE(test_proved_equivalent(
        specialized, SymbolicExpr::number(-1)));
    auto active_condition =
        specialize_discriminant_condition(root.conditions.front(), 2, 1);
    auto inactive_condition =
        specialize_discriminant_condition(root.conditions.front(), 0, -1);
    EXPECT_EQ(condition_evaluator.evaluate_condition(*active_condition),
              Tribool::True);
    EXPECT_EQ(condition_evaluator.evaluate_condition(*inactive_condition),
              Tribool::False);
}

static void expect_simple_discriminant_branch(
    const FiniteSolution &root, const RelationalNode &relation,
    AssumptionContext &condition_evaluator, bool &found_negative_one,
    bool &found_positive_one) {
    EXPECT_EQ(relation.op(), RelationOp::NEQ);
    auto specialized = specialize_discriminant_condition(root.value, 0, -1)
                           ->simplify();
    found_negative_one =
        found_negative_one ||
        test_proved_equivalent(specialized, SymbolicExpr::number(-1));
    found_positive_one =
        found_positive_one ||
        test_proved_equivalent(specialized, SymbolicExpr::number(1));
    auto active_condition =
        specialize_discriminant_condition(root.conditions.front(), 0, -1);
    auto inactive_condition =
        specialize_discriminant_condition(root.conditions.front(), 2, 1);
    EXPECT_EQ(condition_evaluator.evaluate_condition(*active_condition),
              Tribool::True);
    EXPECT_EQ(condition_evaluator.evaluate_condition(*inactive_condition),
              Tribool::False);
}

TEST(LmcasSolveVisitorPoly, UnknownDiscriminantCarriesBothMultiplicityBranches) {
    auto x = SymbolicExpr::variable("x");
    auto b = SymbolicExpr::variable("b");
    auto c = SymbolicExpr::variable("c");
    auto equation = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::add(SymbolicExpr::multiply(b, x), c));

    auto solved = solve_equation(equation, "x");

    ASSERT_TRUE(solved) << solved.error().message;
    const auto *finite = std::get_if<FiniteSolutions>(&solved.value());
    ASSERT_NE(finite, nullptr);
    ASSERT_EQ(finite->values.size(), 3U);
    std::size_t simple_roots = 0;
    std::size_t repeated_roots = 0;
    bool found_negative_one = false;
    bool found_positive_one = false;
    AssumptionContext condition_evaluator;
    for (const auto &root : finite->values) {
        ASSERT_EQ(root.conditions.size(), 1U);
        auto relation = std::dynamic_pointer_cast<const RelationalNode>(
            detail::node(root.conditions.front()));
        ASSERT_NE(relation, nullptr);
        if (root.multiplicity == 2U) {
            ++repeated_roots;
            expect_repeated_discriminant_branch(
                root, *relation, condition_evaluator);
        } else {
            ASSERT_EQ(root.multiplicity, 1U);
            ++simple_roots;
            expect_simple_discriminant_branch(
                root, *relation, condition_evaluator, found_negative_one,
                found_positive_one);
        }
    }
    EXPECT_EQ(simple_roots, 2U);
    EXPECT_EQ(repeated_roots, 1U);
    EXPECT_TRUE(found_negative_one);
    EXPECT_TRUE(found_positive_one);
}
