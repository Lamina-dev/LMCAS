#include "test_common.hpp"
#include "assumption_context.hpp"
#include "property_store.hpp"
#include "computation_context.hpp"
#include "symbolic.hpp"
#include "internal/symbolic_ast.hpp"
#include "expr.hpp"

#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <variant>

using namespace LMCAS;

static void expect_approximate_atom(const SymbolicExpr &expression, double expected) {
    auto value = evalf(expression);
    auto approximate = approx_real(expected);
    auto exact = SymbolicExpr::number(Rational::from_double(expected));
    EXPECT_TRUE(value && std::memcmp(&value.value().value, &expected, sizeof(expected)) == 0) << "serialized numeric atom preserves all binary64 bits including signed zero";
    EXPECT_TRUE(approximate &&
                detail::node(expression)->equals(*detail::node(*approximate.value())) &&
                !detail::node(expression)->equals(*detail::node(*exact)))
        << "serialized approximate atom does not silently become exact";
}

static void expect_approximate_assumption_values(AssumptionContext &restored, double value) {
    const auto intervals = restored.current_properties().get_continuity_decls("f");
    EXPECT_TRUE(intervals.size() == 1 && intervals.front().interval.lower.value) << "continuity retains its endpoint";
    if (intervals.size() == 1 && intervals.front().interval.lower.value) {
        expect_approximate_atom(*intervals.front().interval.lower.value, value);
    }
    const auto &relations = restored.current_relations().get_relations();
    EXPECT_TRUE(relations.size() == 1) << "explicit numeric relation survives";
    if (relations.size() == 1) {
        expect_approximate_atom(relations.front().rhs, value);
    }
    const auto conditionals = restored.get_active_conditionals();
    EXPECT_TRUE(conditionals.size() == 1) << "conditional numeric atoms survive";
    if (conditionals.size() == 1) {
        for (const auto *expression : {&conditionals.front().condition, &conditionals.front().conclusion}) {
            const auto *relation = dynamic_cast<const RelationalNode *>(detail::node(*expression).get());
            EXPECT_TRUE(relation != nullptr) << "conditional retains its relational structure";
            if (relation) {
                expect_approximate_atom(detail::expression_from_node(relation->right()), value);
            }
        }
    }
    if (value > 0) {
        auto period = restored.get_period("g", "x");
        EXPECT_TRUE(period.has_value()) << "approximate declared period survives";
        if (period) {
            expect_approximate_atom(*period, value);
        }
    }
}

TEST(LmcasAssumptionSerialization, CanonicalApproximateRecords) {
    for (double value : {0.5, 1.0, 1.2345678901234567,
                         std::numeric_limits<double>::denorm_min(), -0.0}) {
        auto approximate = approx_real(value);
        EXPECT_TRUE(approximate.has_value()) << "finite approximate fixture is constructible";
        if (!approximate) {
            continue;
        }
        const std::string token = approximate.value()->to_string();
        std::string source = "SCOPE 0\nCONTINUOUS f [" + token + ", approx(2)]\n" +
                             "RELATION r GT " + token + "\nCONDITIONAL (c GT " + token +
                             ") => (d LT " + token + ")\n";
        if (value > 0) {
            source += "PERIODIC g x " + token + "\n";
        }
        source += "END\n";
        auto parsed = AssumptionContext::deserialize_checked(source);
        EXPECT_TRUE(parsed.has_value()) << "approximate assumption records parse";
        if (!parsed) {
            continue;
        }
        auto restored = AssumptionContext::deserialize_checked(parsed.value().serialize());
        EXPECT_TRUE(restored.has_value()) << "approximate assumption records reparse";
        if (!restored) {
            continue;
        }
        expect_approximate_assumption_values(restored.value(), value);
    }
}

TEST(LmcasAssumptionSerialization, LegacyApproximateRecords) {
    for (const auto &item : {std::pair<const char *, double>{"0.5", 0.5}, {"1.0", 1.0}, {"1e-3", 1e-3}}) {
        auto legacy = AssumptionContext::deserialize_checked(
            std::string("SCOPE 0\nRELATION r GT ") + item.first + "\nEND\n");
        EXPECT_TRUE(legacy.has_value()) << "legacy decimal relation remains accepted";
        if (!legacy) {
            continue;
        }
        const auto &relations = legacy.value().current_relations().get_relations();
        EXPECT_TRUE(relations.size() == 1) << "legacy numeric relation exists";
        if (relations.size() == 1) {
            expect_approximate_atom(relations.front().rhs, item.second);
        }
    }
}

TEST(LmcasAssumptionSerialization, MalformedApproximateRecords) {
    for (const char *invalid : {"approx()", "approx(x)", "approx(1/2)", "approx(1,2)",
                                "approx(1e)", "approx(nan)", "approx(inf)", "approx(0x1p0)",
                                "approx(1e400)", "approx(1e-4000)", "approx(1)junk"}) {
        auto relation = AssumptionContext::deserialize_checked(
            std::string("SCOPE 0\nRELATION r GT ") + invalid + "\nEND\n");
        auto period = AssumptionContext::deserialize_checked(
            std::string("SCOPE 0\nPERIODIC f x ") + invalid + "\nEND\n");
        EXPECT_TRUE(!relation && relation.error().code == CasErrc::ParseError &&
                    !period && period.error().code == CasErrc::ParseError)
            << "all numeric atom consumers reject malformed approximate literals";
    }
}

static void expect_bounded_classifications(const AssumptionContext &context) {
    EXPECT_EQ(context.get_boundedness("nested"), Boundedness::Bounded);
    EXPECT_EQ(context.get_boundedness("exact"), Boundedness::Bounded);
    EXPECT_EQ(context.get_boundedness("approximate"), Boundedness::Bounded);
    EXPECT_EQ(context.get_boundedness("symbolic"), Boundedness::Bounded);
}

static void expect_nested_bounds(const AssumptionContext &context) {
    const auto bounds = context.get_bounds("nested");
    ASSERT_TRUE(bounds);
    ASSERT_TRUE(bounds->lower.value);
    ASSERT_TRUE(bounds->upper.value);
    auto lower = evalf(*bounds->lower.value);
    auto upper = evalf(*bounds->upper.value);
    ASSERT_TRUE(lower);
    ASSERT_TRUE(upper);
    EXPECT_NEAR(lower.value().value, 3.0, 1e-12);
    EXPECT_NEAR(upper.value().value, 4.0, 1e-12);
}

static void expect_exact_bounds(const AssumptionContext &context) {
    const auto bounds = context.get_bounds("exact");
    ASSERT_TRUE(bounds);
    ASSERT_TRUE(bounds->lower.value);
    ASSERT_TRUE(bounds->upper.value);
    const auto *lower = dynamic_cast<const NumberNode *>(
        detail::node(*bounds->lower.value).get());
    const auto *upper = dynamic_cast<const NumberNode *>(
        detail::node(*bounds->upper.value).get());
    ASSERT_NE(lower, nullptr);
    ASSERT_NE(upper, nullptr);
    ASSERT_TRUE(std::holds_alternative<Rational>(lower->value()));
    ASSERT_TRUE(std::holds_alternative<Rational>(upper->value()));
    EXPECT_EQ(std::get<Rational>(lower->value()), Rational(1, 3));
    EXPECT_EQ(std::get<Rational>(upper->value()), Rational(2, 3));
}

static void expect_approximate_bounds(const AssumptionContext &context) {
    const auto bounds = context.get_bounds("approximate");
    ASSERT_TRUE(bounds);
    ASSERT_TRUE(bounds->lower.value);
    ASSERT_TRUE(bounds->upper.value);
    const auto *lower = dynamic_cast<const NumberNode *>(
        detail::node(*bounds->lower.value).get());
    const auto *upper = dynamic_cast<const NumberNode *>(
        detail::node(*bounds->upper.value).get());
    ASSERT_NE(lower, nullptr);
    ASSERT_NE(upper, nullptr);
    ASSERT_TRUE(std::holds_alternative<lmmc_real_t>(lower->value()));
    ASSERT_TRUE(std::holds_alternative<lmmc_real_t>(upper->value()));
    EXPECT_DOUBLE_EQ(std::get<lmmc_real_t>(lower->value()), 0.125);
    EXPECT_DOUBLE_EQ(std::get<lmmc_real_t>(upper->value()), 0.875);
}

static void expect_symbolic_bounds(const AssumptionContext &context) {
    const auto bounds = context.get_bounds("symbolic");
    ASSERT_TRUE(bounds);
    ASSERT_TRUE(bounds->lower.value);
    ASSERT_TRUE(bounds->upper.value);
    auto expected = parse_expr("a + 1");
    ASSERT_TRUE(expected);
    EXPECT_TRUE(detail::node(*bounds->lower.value)
                    ->equals(*detail::node(*expected.value())));
    EXPECT_TRUE(detail::node(*bounds->upper.value)
                    ->equals(*detail::node(*expected.value())));
}

static void expect_rich_relation(const AssumptionContext &context) {
    auto lhs = parse_expr("x + sin(y)");
    auto rhs = parse_expr("sqrt(z + 1)");
    ASSERT_TRUE(lhs);
    ASSERT_TRUE(rhs);
    auto canonical_lhs = lhs.value()->simplify();
    auto canonical_rhs = rhs.value()->simplify();
    ASSERT_TRUE(canonical_lhs);
    ASSERT_TRUE(canonical_rhs);
    EXPECT_TRUE(context.has_relation(
        *canonical_lhs, *canonical_rhs, RelationOp::GT));
}

static void expect_rich_conditional(const AssumptionContext &context) {
    const auto conditionals = context.get_active_conditionals();
    ASSERT_EQ(conditionals.size(), 1u);
    auto expected_condition = parse_expr("x + cos(y) > 1/3");
    auto expected_conclusion = parse_expr("log(z, 10) < 2/3");
    ASSERT_TRUE(expected_condition);
    ASSERT_TRUE(expected_conclusion);
    auto canonical_condition = expected_condition.value()->simplify();
    auto canonical_conclusion = expected_conclusion.value()->simplify();
    ASSERT_TRUE(canonical_condition);
    ASSERT_TRUE(canonical_conclusion);
    EXPECT_TRUE(detail::node(conditionals.front().condition)
                    ->equals(*detail::node(canonical_condition)));
    EXPECT_TRUE(detail::node(conditionals.front().conclusion)
                    ->equals(*detail::node(canonical_conclusion)));
}

static void expect_rich_serialization_semantics(
    const AssumptionContext &context) {
    expect_bounded_classifications(context);
    expect_nested_bounds(context);
    expect_exact_bounds(context);
    expect_approximate_bounds(context);
    expect_symbolic_bounds(context);
    expect_rich_relation(context);
    expect_rich_conditional(context);
}

TEST(LmcasAssumptionSerialization, RichRecordsRoundtripSemantically) {
    const std::string source =
        "SCOPE 0\n"
        "BOUNDED nested Bounded [log(8, 2), log(10000, 10)]\n"
        "BOUNDED exact Bounded [1/3, 2/3]\n"
        "BOUNDED approximate Bounded [approx(0.125), approx(0.875)]\n"
        "BOUNDED symbolic Bounded [a + 1, a + 1]\n"
        "RELATION x + sin(y) GT sqrt(z + 1)\n"
        "CONDITIONAL (x + cos(y) GT 1/3) => (log(z, 10) LT 2/3)\n"
        "END\n";
    auto parsed = AssumptionContext::deserialize_checked(source);
    ASSERT_TRUE(parsed);
    expect_rich_serialization_semantics(parsed.value());

    auto restored =
        AssumptionContext::deserialize_checked(parsed.value().serialize());
    ASSERT_TRUE(restored);
    expect_rich_serialization_semantics(restored.value());
}

static void expect_parse_error_on_line_two(const std::string &record) {
    auto result = AssumptionContext::deserialize_checked(
        "SCOPE 0\n" + record + "\nEND\n");
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::ParseError);
    EXPECT_TRUE(result.error().message.find("Line 2") != std::string::npos ||
                result.error().message.find("line 2") != std::string::npos);
}

TEST(LmcasAssumptionSerialization, MalformedNestedRecordsReportTheirLine) {
    for (const char *record : {
             "RELATION f(x, g(y, 2) GT z",
             "CONDITIONAL (x GT f(y, 2)) => (g(z, 3) LT 1",
             "BOUNDED x Bounded [log(8, 2), log(100, 10)",
             "BOUNDED x Unknown [0, 1]",
             "BOUNDED x Unbounded [0, 1]"}) {
        expect_parse_error_on_line_two(record);
    }
}

TEST(LmcasAssumptionSerialization, CheckedParserPreservesContextFailures) {
    const std::string source =
        "SCOPE 0\nBOUNDED x Bounded [0, 1]\nEND\n";

    ResourceLimits limits;
    limits.max_input_bytes = source.size() - 1;
    ComputationContext limited(limits);
    auto exhausted =
        AssumptionContext::deserialize_checked(source, limited);
    ASSERT_FALSE(exhausted);
    EXPECT_EQ(exhausted.error().code, CasErrc::ResourceLimit);

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled({}, cancellation);
    auto stopped =
        AssumptionContext::deserialize_checked(source, cancelled);
    ASSERT_FALSE(stopped);
    EXPECT_EQ(stopped.error().code, CasErrc::Cancelled);
}

TEST(LmcasAssumptionSerialization, RelationAlternativeProofRoundtrip) {
    auto source = AssumptionContext::deserialize_checked(
        "SCOPE 0\nRELATION a GT b\nRELATION b GT c\n"
        "RELATION a GT d\nRELATION d GT c\nEND\n");
    ASSERT_TRUE(source);
    auto restored = AssumptionContext::deserialize_checked(source.value().serialize());
    ASSERT_TRUE(restored);
    const auto a = SymbolicExpr::variable("a");
    const auto c = SymbolicExpr::variable("c");
    for (auto* context : {&source.value(), &restored.value()}) {
        EXPECT_TRUE(context->has_relation(*a, *c, RelationOp::GT));
        context->push();
        ASSERT_TRUE(context->assume_domain("b", Domain::Integer));
        EXPECT_TRUE(context->has_relation(*a, *c, RelationOp::GT));
        ASSERT_TRUE(context->assume_domain("d", Domain::Integer));
        EXPECT_FALSE(context->has_relation(*a, *c, RelationOp::GT));
        ASSERT_TRUE(context->pop());
        EXPECT_TRUE(context->has_relation(*a, *c, RelationOp::GT));
    }
}

TEST(LmcasAssumptionSerialization, DerivedThenDeclaredRelationRoundtrip) {
    auto source = AssumptionContext::deserialize_checked(
        "SCOPE 0\nRELATION a GT b\nRELATION b GT c\n"
        "RELATION a GT d\nRELATION d GT c\nRELATION a GT c\nEND\n");
    ASSERT_TRUE(source);
    auto restored = AssumptionContext::deserialize_checked(source.value().serialize());
    ASSERT_TRUE(restored);
    const auto a = SymbolicExpr::variable("a");
    const auto c = SymbolicExpr::variable("c");
    for (auto* context : {&source.value(), &restored.value()}) {
        context->push();
        ASSERT_TRUE(context->assume_domain("b", Domain::Integer));
        ASSERT_TRUE(context->assume_domain("d", Domain::Integer));
        EXPECT_TRUE(context->has_relation(*a, *c, RelationOp::GT));
        ASSERT_TRUE(context->pop());
        EXPECT_TRUE(context->has_relation(*a, *c, RelationOp::GT));
    }
}

TEST(LmcasAssumptionSerialization, RepeatedDeclarationContinuesClosureAfterRoundtrip) {
    std::string declarations = "SCOPE 0\n";
    for (int i = 0; i < 65; ++i) {
        declarations += "RELATION star_" + std::to_string(i) + " GT center\n";
    }
    declarations += "RELATION center GT bottom\n";
    auto once = AssumptionContext::deserialize_checked(declarations + "END\n");
    ASSERT_TRUE(once);
    const auto star = SymbolicExpr::variable("star_64");
    const auto bottom = SymbolicExpr::variable("bottom");
    EXPECT_FALSE(once.value().has_relation(*star, *bottom, RelationOp::GT));

    auto repeated = AssumptionContext::deserialize_checked(
        declarations + "RELATION center GT bottom\nEND\n");
    ASSERT_TRUE(repeated);
    ASSERT_TRUE(repeated.value().has_relation(*star, *bottom, RelationOp::GT));
    auto restored = AssumptionContext::deserialize_checked(repeated.value().serialize());
    ASSERT_TRUE(restored);
    for (auto* context : {&repeated.value(), &restored.value()}) {
        EXPECT_TRUE(context->has_relation(*star, *bottom, RelationOp::GT));
        context->push();
        ASSERT_TRUE(context->assume_domain("center", Domain::Integer));
        EXPECT_FALSE(context->has_relation(*star, *bottom, RelationOp::GT));
        ASSERT_TRUE(context->pop());
        EXPECT_TRUE(context->has_relation(*star, *bottom, RelationOp::GT));
    }
}
